/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_NFFT_HPP
#define RTLPSMRI_CORE_NFFT_HPP

namespace mri {

// TODO: Decide whether or not it is worthwhile to reuse plans. At worst, we
//       would need one plan per coil.

// TODO: Decide how to multithread embarrassingly parallel local transforms

inline void
NFFT2D
( int N0, int N1, int numNonUniform, int n0, int n1, int m, 
  const DistMatrix<Complex<double>,STAR,VR>& FHat, 
  const DistMatrix<double,         STAR,VR>& paths,
        DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("NFFT2D");
#endif
    const int dim = 2;
    const int width = paths.Width();
#ifndef RELEASE
    if( N0 % 2 != 0 || N1 % 2 != 0 )
        LogicError("NFFT requires band limits to be even integers\n");
    if( FHat.Height() != N0*N1 )
        LogicError("Invalid FHat height");
    if( paths.Height() != dim*numNonUniform )
        LogicError("Invalid paths height");
    if( FHat.Width() != paths.Width() )
        LogicError("FHat and paths must have the same width");
    if( FHat.RowAlign() != paths.RowAlign() )
        LogicError("FHat and paths are not aligned");
#endif
    F.AlignWith( FHat );
    Zeros( F, numNonUniform, width );
    const int locWidth = F.LocalWidth();

    nfft_plan plan; 
    int NN[dim] = { N0, N1 };
    int nn[dim] = { n0, n1 };

    unsigned nfftFlags = PRE_PHI_HUT| PRE_FULL_PSI| FFTW_INIT| FFT_OUT_OF_PLACE;
    unsigned fftwFlags = FFTW_MEASURE| FFTW_DESTROY_INPUT;

    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
#ifndef RELEASE
        // TODO: Ensure this column of paths is sorted
#endif
        plan.x = const_cast<double*>(paths.LockedBuffer(0,jLoc));
        nfft_init_guru
        ( &plan, dim, NN, numNonUniform, nn, m, nfftFlags, fftwFlags );
        if( plan.nfft_flags & PRE_ONE_PSI )
            nfft_precompute_one_psi( &plan ); 
        plan.f_hat = (fftw_complex*)
            const_cast<Complex<double>*>(FHat.LockedBuffer(0,jLoc));
        plan.f = (fftw_complex*)F.Buffer(0,jLoc);
        nfft_trafo_2d( &plan );
        nfft_finalize( &plan );
    }
    const double scale = 1./Sqrt(1.*N0*N1);
    Scale( scale, F );
}

inline void
AdjointNFFT2D
( int N0, int N1, int numNonUniform, int n0, int n1, int m,
  const DistMatrix<Complex<double>,STAR,VR>& F,
  const DistMatrix<double,         STAR,VR>& paths,
        DistMatrix<Complex<double>,STAR,VR>& FHat )
{
#ifndef RELEASE
    CallStackEntry cse("AdjointNFFT2D");
#endif
    const int dim = 2;
    const int width = paths.Width();
#ifndef RELEASE
    if( N0 % 2 != 0 || N1 % 2 != 0 )
        LogicError("NFFT requires band limits to be even integers\n");
    if( F.Height() != numNonUniform )
        LogicError("Invalid F height");
    if( paths.Height() != dim*numNonUniform )
        LogicError("Invalid paths height");
    if( F.Width() != paths.Width() )
        LogicError("F and paths must have the same width");
    if( F.RowAlign() != paths.RowAlign() )
        LogicError("F and paths are not aligned");
#endif
    FHat.AlignWith( F );
    Zeros( FHat, N0*N1, width );
    const int locWidth = FHat.LocalWidth();

    nfft_plan plan;
    int NN[dim] = { N0, N1 };
    int nn[dim] = { n0, n1 };

    unsigned nfftFlags = PRE_PHI_HUT| PRE_FULL_PSI| FFTW_INIT| FFT_OUT_OF_PLACE;
    unsigned fftwFlags = FFTW_MEASURE| FFTW_DESTROY_INPUT;

    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
#ifndef RELEASE
        // TODO: Ensure this column of paths is sorted
#endif
        plan.x = const_cast<double*>(paths.LockedBuffer(0,jLoc));
        nfft_init_guru
        ( &plan, dim, NN, numNonUniform, nn, m, nfftFlags, fftwFlags );
        if( plan.nfft_flags & PRE_ONE_PSI )
            nfft_precompute_one_psi( &plan ); 
        plan.f = (fftw_complex*)
            const_cast<Complex<double>*>(F.LockedBuffer(0,jLoc));
        plan.f_hat = (fftw_complex*)FHat.Buffer(0,jLoc);
        nfft_adjoint( &plan );
        nfft_finalize( &plan );
    }
    const double scale = 1./Sqrt(1.*N0*N1);
    Scale( scale, FHat );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_NFFT_HPP
