/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
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
#ifndef RELEASE
    CallStackEntry cse("lps::UpdateZ");
#endif
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
#ifndef RELEASE
    CallStackEntry cse("lps::SubtractAdjDz");
#endif
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
  bool tryTSQR=false )
{
#ifndef RELEASE
    CallStackEntry cse("LPS");
#endif
    typedef double Real;
    typedef Complex<Real> F;

    const int numTimesteps = NumTimesteps();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();

    // M := E' D
    DistMatrix<F,VC,STAR> M( D.Grid() );
    AdjointAcquisition( D, M );

    // Set lambdaS relative to || M ||_max
    const double maxM = MaxNorm( M );
    const double lambdaS = lambdaSRelMaxM*maxM;

#ifndef RELEASE
    DistMatrix<F,STAR,VR> R( M.Grid() );
    Acquisition( M, R );
    const double frobM = FrobeniusNorm( M );
    const double frobR = FrobeniusNorm( R );
    if( D.Grid().Rank() == 0 )
        std::cout << "|| M= E'D ||_F = " << frobM << "\n"
                  << "|| R=EE'D ||_F = " << frobR << "\n"
                  << "lambdaL=" << lambdaL << ", lambdaS=" << lambdaS
                  << std::endl;
#endif

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

    int numIts=0;
    DistMatrix<F,VC,STAR> M0( M.Grid() );
    while( true )
    {
        ++numIts;

        // M0 := M
        M0 = M;

        // L := SVT(M-S,lambdaL)
        L = M;
        Axpy( F(-1), S, L );
        int rank;
        if( tryTSQR )
            rank = elem::svt::TSQR( L, lambdaL, true );
        else 
            rank = elem::svt::TallCross( L, lambdaL, true );

        // S := TransformedST(M-L)
        int numNonzeros;
        S = M;
        Axpy( F(-1), L, S );
        if( tv )
        {
            lps::UpdateZ( lambdaS/2, S, Z );
            lps::SubtractAdjDz( Z, S );
            numNonzeros = ZeroNorm( S );
        }
        else
        {
            TemporalFFT( S );
            elem::SoftThreshold( S, lambdaS );
#ifndef RELEASE
            numNonzeros = ZeroNorm( S );
#endif
            TemporalAdjointFFT( S );
        }

        // M := L + S - E'(E(L+S)-D)
        M = L;
        Axpy( F(1), S, M );
        Acquisition( M, R );        
        Axpy( F(-1), D, R );
        AdjointAcquisition( R, M );
        Scale( F(-1), M );
        Axpy( F(1), L, M );
        Axpy( F(1), S, M );

        const Real frobM0 = FrobeniusNorm( M0 );        
        Axpy( F(-1), M, M0 );
        const Real frobUpdate = FrobeniusNorm( M0 );
#ifndef RELEASE
        const double frobL = FrobeniusNorm( L );
        const double frobS = FrobeniusNorm( S );
        if( D.Grid().Rank() == 0 )
            std::cout << "After " << numIts << " its: \n"
                      << "  rank(L)      = " << rank << "\n"
                      << "  nnz(TS)      = " << numNonzeros << "\n"
                      << "  || L    ||_F = " << frobL << "\n"
                      << "  || S    ||_F = " << frobS << "\n"
                      << "  || M0   ||_F = " << frobM0 << "\n"
                      << "  || M-M0 ||_F = " << frobUpdate << "\n"
                      << "  || M-M0 ||_F / || M ||_F = " << frobUpdate/frobM0
                      << std::endl;
#endif
        if( numIts == maxIts || frobUpdate < relTol*frobM0 )
            break;
    }
    if( numIts == maxIts )
        RuntimeError("L+S decomposition did not converge in time");

    return numIts;
}

} // namespace mri

#endif // ifndef RTLPSMRI_LPS_HPP
