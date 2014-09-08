/*
   Copyright (c) 2013-2014, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI).

   RT-LPS-MRI is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   RT-LPS-MRI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with RT-LPS-MRI.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once
#ifndef RTLPSMRI_LPS_HPP
#define RTLPSMRI_LPS_HPP

namespace mri {

namespace lps {

inline void
UpdateZ
( double clipRadius,
  const DistMatrix<Complex<double>,VC,STAR>& S, 
        DistMatrix<Complex<double>,VC,STAR>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("lps::UpdateZ"))
    const int numTimesteps = S.Width();
    const int localHeight = S.LocalHeight();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        for( int j=0; j<numTimesteps-1; ++j )
        {
            const Complex<double> curr = S.GetLocal(iLoc,j);
            const Complex<double> next = S.GetLocal(iLoc,j+1);
            Z.UpdateLocal( iLoc, j, 0.25*(next-curr) ); 
            // Clip magnitude to be less than or equal to clipRadius
            const Complex<double> xi = Z.GetLocal(iLoc,j);
            const double xiAbs = Abs(xi);
            if( xiAbs > clipRadius )
                Z.SetLocal( iLoc, j, xi*(clipRadius/xiAbs) );    
        }
    }
}

inline void
SubtractAdjDz
( const DistMatrix<Complex<double>,VC,STAR>& Z,
        DistMatrix<Complex<double>,VC,STAR>& S )
{
    DEBUG_ONLY(CallStackEntry cse("lps::SubtractAdjDz"))
    const int numTimesteps = S.Width();
    const int localHeight = S.LocalHeight();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        // Add the first column of Z to the first column of S
        S.UpdateLocal( iLoc, 0, Z.GetLocal(iLoc,0) );
        // Add row-wise diff of Z to the middle
        for( int j=1; j<numTimesteps-1; ++j )
        {
            const Complex<double> curr = Z.GetLocal(iLoc,j-1);
            const Complex<double> next = Z.GetLocal(iLoc,j);
            S.UpdateLocal( iLoc, j, next-curr ); 
        }
        // Subtract the last column of Z from the last column of S
        S.UpdateLocal( iLoc, numTimesteps-1, -Z.GetLocal(iLoc,numTimesteps-2) );
    }
}

} // namespace lps

inline int
LPS
( const DistMatrix<Complex<double>,STAR,VR>& D,
        DistMatrix<Complex<double>,VC,STAR>& L,
        DistMatrix<Complex<double>,VC,STAR>& S,
  bool tv=true,
  double lambdaL=0.025, double lambdaSRelMaxM=0.5,
  double relTol=0.0025, int maxIts=100,
  bool tryTSQR=false, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("LPS"))
    using El::Timer;
    typedef double Real;
    typedef Complex<Real> F;

    const int numTimesteps = NumTimesteps();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();

    const bool amRoot = D.Grid().Rank() == 0;

    Timer initial("initial");
    initial.Start();

    // M := E' D
    DistMatrix<F,VC,STAR> M( D.Grid() );
    AdjointAcquisition( D, M, progress );

    // Set lambdaS relative to || M ||_max
    const double maxM = MaxNorm( M );
    const double lambdaS = lambdaSRelMaxM*maxM;

    if( progress )
    {
        const double frobM = FrobeniusNorm( M );
        if( amRoot )
            std::cout << "|| M= E'D ||_F = " << frobM << "\n"
                      << "lambdaL=" << lambdaL << ", lambdaS=" << lambdaS
                      << std::endl;
    }

    // Align L and S to M, and set S := 0
    L.SetGrid( M.Grid() );
    S.SetGrid( M.Grid() );
    L.AlignWith( M );
    S.AlignWith( M );
    Zeros( S, N0*N1, numTimesteps );

    // If using TV clipping, we need to accumulate data in the matrix Z 
    DistMatrix<F,VC,STAR> Z( M.Grid() );
    if( tv )
    {
        Z.AlignWith( M );
        Zeros( Z, N0*N1, numTimesteps-1 );
    }


    if( progress && amRoot )
        std::cout << "initialization time: " << initial.Stop() << std::endl;

    int numIts=0;
    Timer svt("svt"), thresh("thresh"), forward("forward"), adjoint("adjoint");
    DistMatrix<F,VC,STAR> M0( M.Grid() );
    DistMatrix<F,STAR,VR> R( M.Grid() );
    while( true )
    {
        ++numIts;

        // M0 := M
        M0 = M;

        // L := SVT(M-S,lambdaL)
        svt.Start();
        L = M;
        Axpy( F(-1), S, L );
        int rank;
        if( tryTSQR )
            rank = El::svt::TSQR( L, lambdaL, true );
        else 
            rank = El::svt::Cross( L, lambdaL, true );
        const double svtTime = svt.Stop();

        // S := TransformedST(M-L)
        thresh.Start();
        int numNonzeros;
        S = M;
        Axpy( F(-1), L, S );
        if( tv )
        {
            lps::UpdateZ( MaxNorm(M)*lambdaS/2, S, Z );
            lps::SubtractAdjDz( Z, S );
            if( progress )
                numNonzeros = ZeroNorm( S );
        }
        else
        {
            TemporalFFT( S );
            El::SoftThreshold( S, lambdaS );
            if( progress )
                numNonzeros = ZeroNorm( S );
            TemporalAdjointFFT( S );
        }
        const double threshTime = thresh.Stop();

        // M := L + S - E'(E(L+S)-D)
        forward.Start();
        M = L;
        Axpy( F(1), S, M );
        Acquisition( M, R, progress ); 
        Axpy( F(-1), D, R );
        const double forwardTime = forward.Stop();
        adjoint.Start();
        AdjointAcquisition( R, M, progress );
        Scale( F(-1), M );
        Axpy( F(1), L, M );
        Axpy( F(1), S, M );
        const double adjointTime = adjoint.Stop();

        const Real frobM0 = FrobeniusNorm( M0 );        
        Axpy( F(-1), M, M0 );
        const Real frobUpdate = FrobeniusNorm( M0 );
        if( progress )
        {
            const double frobL = FrobeniusNorm( L );
            const double frobS = FrobeniusNorm( S );
            double frobZ;
            if( tv )
                frobZ = FrobeniusNorm( Z );
            if( amRoot )
            {
                std::cout << "After " << numIts << " its: \n"
                          << "  rank(L)      = " << rank << "\n"
                          << "  nnz(TS)      = " << numNonzeros << "\n"
                          << "  || L    ||_F = " << frobL << "\n"
                          << "  || S    ||_F = " << frobS << "\n";
                if( tv )
                    std::cout 
                          << "  || Z    ||_F = " << frobZ << "\n";
                std::cout << "  || M0   ||_F = " << frobM0 << "\n"
                          << "  || M-M0 ||_F = " << frobUpdate << "\n"
                          << "  || M-M0 ||_F / || M0 ||_F = " 
                          << frobUpdate/frobM0 << "\n"
                          << "  SVT time:     " << svtTime << " seconds\n"
                          << "  Thresh time:  " << threshTime << " seconds\n"
                          << "  Forward time: " << forwardTime << " seconds\n"
                          << "  Adjoint time: " << adjointTime << " seconds\n"
                          << std::endl;
            }
        }
        if( numIts == maxIts || frobUpdate < relTol*frobM0 )
            break;
    }
    return numIts;
}

} // namespace mri

#endif // ifndef RTLPSMRI_LPS_HPP
