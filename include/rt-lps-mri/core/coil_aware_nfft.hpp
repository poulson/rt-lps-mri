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
( const DistMatrix<int,            STAR,VR>& coils,
  const DistMatrix<Complex<double>,STAR,VR>& FHat, 
        DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("CoilAwareNFFT2D");
#endif
    const int width = coils.Width();
    const int M = NumNonUniformPoints();
#ifndef RELEASE
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
    if( N0 % 2 != 0 || N1 % 2 != 0 )
        LogicError("NFFT requires band limits to be even integers\n");
    if( FHat.Height() != N0*N1 )
        LogicError("Invalid FHat height");
    if( FHat.Width() != width )
        LogicError("FHat and coils must have the same width");
    if( FHat.RowAlign() != coils.RowAlign() )
        LogicError("FHat and coils are not aligned");
#endif
    F.AlignWith( FHat );
    Zeros( F, M, width );
    const int locWidth = F.LocalWidth();
    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        nfft_plan& p = CoilPlan( coils.GetLocal(0,jLoc) );
        p.f_hat = (fftw_complex*)
            const_cast<Complex<double>*>(FHat.LockedBuffer(0,jLoc));
        p.f = (fftw_complex*)F.Buffer(0,jLoc);
        nfft_trafo_2d( &p );
    }
}

inline void
CoilAwareAdjointNFFT2D
( const DistMatrix<int,            STAR,VR>& coils,      
        DistMatrix<Complex<double>,STAR,VR>& FHat,
  const DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("CoilAwareAdjointNFFT2D");
#endif
    const int width = coils.Width();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
#ifndef RELEASE
    const int M = NumNonUniformPoints();
    if( N0 % 2 != 0 || N1 % 2 != 0 )
        LogicError("NFFT requires band limits to be even integers\n");
    if( F.Height() != M )
        LogicError("Invalid F height");
    if( F.Width() != width )
        LogicError("F and coils must have the same width");
    if( F.RowAlign() != coils.RowAlign() )
        LogicError("F and coils are not aligned");
#endif
    FHat.AlignWith( F );
    Zeros( FHat, N0*N1, width );
    const int locWidth = FHat.LocalWidth();
    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        nfft_plan& p = CoilPlan( coils.GetLocal(0,jLoc) );
        p.f = (fftw_complex*)
            const_cast<Complex<double>*>(F.LockedBuffer(0,jLoc));
        p.f_hat = (fftw_complex*)FHat.Buffer(0,jLoc);
        nfft_adjoint( &p );
    }
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP
