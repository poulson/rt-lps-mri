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
#ifndef RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP
#define RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP

namespace mri {

// TODO: Decide how to multithread embarrassingly parallel local transforms

inline void
CoilAwareNFFT2D
( const DistMatrix<Complex<double>,STAR,VR>& FHat, 
        DistMatrix<Complex<double>,STAR,VR>& F )
{
    DEBUG_ONLY(CallStackEntry cse("CoilAwareNFFT2D"))
    const int width = FHat.Width();
    const int numNonUniform = NumNonUniformPoints();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
    const int numCoils = NumCoils();
    DEBUG_ONLY(
        const int numTimesteps = NumTimesteps();
        if( numCoils*numTimesteps != width )
            LogicError("Invalid width");
        if( N0 % 2 != 0 || N1 % 2 != 0 )
            LogicError("NFFT requires band limits to be even integers\n");
        if( FHat.Height() != N0*N1 )
            LogicError("Invalid FHat height");
    )
    F.AlignWith( FHat );
    Zeros( F, numNonUniform, width );
    const int locWidth = F.LocalWidth();
    const int rowShift = F.RowShift();
    const int rowStride = F.RowStride();
    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int t = j / numCoils;
        nfft_plan& p = CoilPlan( t );
        p.f_hat = (fftw_complex*)
            const_cast<Complex<double>*>(FHat.LockedBuffer(0,jLoc));
        p.f = (fftw_complex*)F.Buffer(0,jLoc);
        nfft_trafo_2d( &p );
    }
    const double scale = 1./Sqrt(1.*N0*N1);
    Scale( scale, F );
}

inline void
CoilAwareAdjointNFFT2D
( const DistMatrix<Complex<double>,STAR,VR>& F,
        DistMatrix<Complex<double>,STAR,VR>& FHat )
{
    DEBUG_ONLY(CallStackEntry cse("CoilAwareAdjointNFFT2D"))
    const int width = F.Width();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
    const int numCoils = NumCoils();
    DEBUG_ONLY(
        const int numTimesteps = NumTimesteps();
        const int numNonUniform = NumNonUniformPoints();
        if( width != numCoils*numTimesteps )
            LogicError("Invalid width");
        if( N0 % 2 != 0 || N1 % 2 != 0 )
            LogicError("NFFT requires band limits to be even integers\n");
        if( F.Height() != numNonUniform )
            LogicError("Invalid F height");
    )
    FHat.AlignWith( F );
    Zeros( FHat, N0*N1, width );
    const int locWidth = F.LocalWidth();
    const int rowShift = F.RowShift();
    const int rowStride = F.RowStride();
    for( int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int t = j / numCoils;
        nfft_plan& p = CoilPlan( t );
        p.f = (fftw_complex*)
            const_cast<Complex<double>*>(F.LockedBuffer(0,jLoc));
        p.f_hat = (fftw_complex*)FHat.Buffer(0,jLoc);
        nfft_adjoint( &p );
    }
    const double scale = 1./Sqrt(1.*N0*N1);
    Scale( scale, FHat );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_COIL_AWARE_NFFT_HPP
