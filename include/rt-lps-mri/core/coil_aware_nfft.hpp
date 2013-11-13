/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP
#define RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP

namespace mri {

// TODO: Decide how to multithread embarrassingly parallel local transforms

inline void
CoilAwareNFFT2D
( const DistMatrix<Complex<double>,STAR,VR>& FHat, 
        DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("CoilAwareNFFT2D");
#endif
    const int width = FHat.Width();
    const int M = NumNonUniformPoints();
#ifndef RELEASE
    const int nc = NumCoils();
    const int nt = NumTimesteps();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
    if( nc*nt != width )
        LogicError("Invalid width");
    if( N0 % 2 != 0 || N1 % 2 != 0 )
        LogicError("NFFT requires band limits to be even integers\n");
    if( FHat.Height() != N0*N1 )
        LogicError("Invalid FHat height");
    if( FHat.LocalWidth() != NumLocalPaths() )
        LogicError("Invalid alignment");
#endif
    F.AlignWith( FHat );
    Zeros( F, M, width );
    const int locWidth = F.LocalWidth();
    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        nfft_plan& p = LocalCoilPlan( jLoc );
        p.f_hat = (fftw_complex*)
            const_cast<Complex<double>*>(FHat.LockedBuffer(0,jLoc));
        p.f = (fftw_complex*)F.Buffer(0,jLoc);
        nfft_trafo_2d( &p );
    }
}

inline void
CoilAwareAdjointNFFT2D
(       DistMatrix<Complex<double>,STAR,VR>& FHat,
  const DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("CoilAwareAdjointNFFT2D");
#endif
    const int width = F.Width();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
#ifndef RELEASE
    const int nc = NumCoils();
    const int nt = NumTimesteps();
    const int M = NumNonUniformPoints();
    if( width != nc*nt )
        LogicError("Invalid width");
    if( N0 % 2 != 0 || N1 % 2 != 0 )
        LogicError("NFFT requires band limits to be even integers\n");
    if( F.Height() != M )
        LogicError("Invalid F height");
    if( F.LocalWidth() != NumLocalPaths() )
        LogicError("Invalid alignment");
#endif
    FHat.AlignWith( F );
    Zeros( FHat, N0*N1, width );
    const int locWidth = FHat.LocalWidth();
    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        nfft_plan& p = LocalCoilPlan( jLoc );
        p.f = (fftw_complex*)
            const_cast<Complex<double>*>(F.LockedBuffer(0,jLoc));
        p.f_hat = (fftw_complex*)FHat.Buffer(0,jLoc);
        nfft_adjoint( &p );
    }
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP
