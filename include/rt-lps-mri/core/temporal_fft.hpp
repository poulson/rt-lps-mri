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
