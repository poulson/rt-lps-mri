/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_TEMPORALFFT_HPP
#define RTLPSMRI_CORE_TEMPORALFFT_HPP

namespace mri {

// TODO: Decide how to avoid horrendous non-unit stride access...

// TODO: Decide how to multithread embarrassingly parallel local transforms

inline void
TemporalFFT( DistMatrix<Complex<double>,VC,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("TemporalFFT");
        if( A.Width() != NumTimesteps() )
            LogicError("Wrong number of timesteps");
    )
    const int numTimesteps = A.Width();
    fftw_complex* buf = 
        (fftw_complex*)fftw_malloc(numTimesteps*sizeof(fftw_complex));
    fftw_plan p = 
        fftw_plan_dft_1d( numTimesteps, buf, buf, FFTW_FORWARD, FFTW_ESTIMATE );
    
    const int numLocFFTs = A.LocalHeight();
    Complex<double>* ABuf = A.Buffer();
    const int ldim = A.LDim();
    for( int k=0; k<numLocFFTs; ++k )
    {
        // Copy in data   
        for( int i=0; i<numTimesteps; ++i )
            std::memcpy( &buf[i], &ABuf[k+i*ldim], sizeof(fftw_complex) );

        // Run the transform
        fftw_execute(p);

        // Copy out data
        for( int i=0; i<numTimesteps; ++i )
            std::memcpy( &ABuf[k+i*ldim], &buf[i], sizeof(fftw_complex) );
    }

    fftw_free( buf );

    const double scale = 1./Sqrt(double(numTimesteps));
    Scale( scale, A );
}

inline void
TemporalAdjointFFT( DistMatrix<Complex<double>,VC,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("TemporalAdjointFFT");
        if( A.Width() != NumTimesteps() )
            LogicError("Wrong number of timesteps");
    )
    const int numTimesteps = A.Width();
    fftw_complex* buf = 
        (fftw_complex*)fftw_malloc(numTimesteps*sizeof(fftw_complex));
    fftw_plan p = 
        fftw_plan_dft_1d
        ( numTimesteps, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE );
    
    const int numLocFFTs = A.LocalHeight();
    Complex<double>* ABuf = A.Buffer();
    const int ldim = A.LDim();
    for( int k=0; k<numLocFFTs; ++k )
    {
        // Copy in data   
        for( int i=0; i<numTimesteps; ++i )
            std::memcpy( &buf[i], &ABuf[k+i*ldim], sizeof(fftw_complex) );

        // Run the transform
        fftw_execute(p);

        // Copy out data
        for( int i=0; i<numTimesteps; ++i )
            std::memcpy( &ABuf[k+i*ldim], &buf[i], sizeof(fftw_complex) );
    }

    fftw_free( buf );

    const double scale = 1./Sqrt(double(numTimesteps));
    Scale( scale, A );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_TEMPORALFFT_HPP
