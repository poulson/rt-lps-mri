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
( int N1, int N2, int M, int n1, int n2, int m, 
  const DistMatrix<Complex<double>,STAR,VR>& FHat, 
  const DistMatrix<double,         STAR,VR>& X,
        DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("NFFT2D");
#endif
    const Int d = 2;
    const Int width = FHat.Width();
#ifndef RELEASE
    if( FHat.Height() != N1*N2 )
        LogicError("Invalid FHat height");
    if( X.Height() != d*M )
        LogicError("Invalid X height");
    if( FHat.Width() != X.Width() )
        LogicError("FHat and X must have the same width");
    if( FHat.RowAlign() != X.RowAlign() )
        LogicError("FHat and X are not aligned");
#endif
    F.AlignWith( FHat );
    Zeros( F, M, width );
    const Int locWidth = F.LocalWidth();

    nfft_plan p; 
    int NN[2] = { N1, N2 };
    int nn[2] = { n1, n2 };

    unsigned nfftFlags = PRE_PHI_HUT| PRE_FULL_PSI| FFTW_INIT| FFT_OUT_OF_PLACE;
    unsigned fftwFlags = FFTW_MEASURE| FFTW_DESTROY_INPUT;

    for( Int jLoc=0; jLoc<locWidth; ++jLoc )
    {
#ifndef RELEASE
        // TODO: Ensure this column of X is sorted
#endif
        p.x = const_cast<double*>(X.LockedBuffer(0,jLoc));
        p.f_hat = (fftw_complex*)
            const_cast<Complex<double>*>(FHat.LockedBuffer(0,jLoc));
        p.f = (fftw_complex*)F.Buffer(0,jLoc);
        nfft_init_guru( &p, d, NN, M, nn, m, nfftFlags, fftwFlags );
        if( p.nfft_flags & PRE_ONE_PSI )
            nfft_precompute_one_psi( &p );  // TODO: See if this can be hoisted
        nfft_trafo_2d( &p );
        nfft_finalize( &p );
    }
}

inline void
AdjointNFFT2D
( int N1, int N2, int M, int n1, int n2, int m,
        DistMatrix<Complex<double>,STAR,VR>& FHat,
  const DistMatrix<double,         STAR,VR>& X,
  const DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("NFFT2D");
#endif
    const Int d = 2;
    const Int width = FHat.Width();
#ifndef RELEASE
    if( F.Height() != M )
        LogicError("Invalid F height");
    if( X.Height() != d*M )
        LogicError("Invalid X height");
    if( F.Width() != X.Width() )
        LogicError("F and X must have the same width");
    if( F.RowAlign() != X.RowAlign() )
        LogicError("F and X are not aligned");
#endif
    FHat.AlignWith( F );
    Zeros( FHat, N1*N2, width );
    const Int locWidth = FHat.LocalWidth();

    nfft_plan p;
    int NN[2] = { N1, N2 };
    int nn[2] = { n1, n2 };

    unsigned nfftFlags = PRE_PHI_HUT| PRE_FULL_PSI| FFTW_INIT| FFT_OUT_OF_PLACE;
    unsigned fftwFlags = FFTW_MEASURE| FFTW_DESTROY_INPUT;

    for( Int jLoc=0; jLoc<locWidth; ++jLoc )
    {
#ifndef RELEASE
        // TODO: Ensure this column of X is sorted
#endif
        p.x = const_cast<double*>(X.LockedBuffer(0,jLoc));
        p.f = (fftw_complex*)
            const_cast<Complex<double>*>(F.LockedBuffer(0,jLoc));
        p.f_hat = (fftw_complex*)FHat.Buffer(0,jLoc);
        nfft_init_guru( &p, d, NN, M, nn, m, nfftFlags, fftwFlags );
        if( p.nfft_flags & PRE_ONE_PSI )
            nfft_precompute_one_psi( &p );  // TODO: See if this can be hoisted
        nfft_adjoint_2d( &p );
        nfft_finalize( &p );
    }
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_NFFT_HPP
