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

inline int
LPS
( const DistMatrix<Complex<double>,STAR,VR>& D,
        DistMatrix<Complex<double>,VC,STAR>& L,
        DistMatrix<Complex<double>,VC,STAR>& S,
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
    const double frobM = FrobeniusNorm( M );
    if( D.Grid().Rank() == 0 )
        std::cout << "|| M=E'D ||_F = " << frobM << "\n"
                  << "lambdaL=" << lambdaL << ", lambdaS=" << lambdaS
                  << std::endl;
#endif

    // Align L and S to M, and set S := 0
    L.SetGrid( M.Grid() );
    S.SetGrid( M.Grid() );
    L.AlignWith( M );
    S.AlignWith( M );
    Zeros( S, N0*N1, numTimesteps );

    int numIts=0;
    DistMatrix<F,VC,STAR> M0( M.Grid() );
    DistMatrix<F,STAR,VR> R( M.Grid() );
    while( true )
    {
        ++numIts;

        // L := SVT(M-S,lambdaL)
        L = M;
        Axpy( F(-1), S, L );
        int rank;
        if( tryTSQR )
            rank = elem::svt::TSQR( L, lambdaL, true );
        else 
            rank = elem::svt::TallCross( L, lambdaL, true );

        // S := TransformedST(M-L)
        // TODO: Add ability to use TV instead of a temporal FFT
        S = M;
        Axpy( F(-1), L, S );
        TemporalFFT( S );
        elem::SoftThreshold( S, lambdaS );
#ifndef RELEASE
        const int numNonzeros = ZeroNorm( S );
#endif
        TemporalAdjointFFT( S );

        // M0 := M
        M0 = M;

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
                      << "  || M-M0 ||_F = " << frobUpdate << std::endl;
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
