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
#ifndef RTLPSMRI_CORE_DIRECT_COIL_AWARE_NFT_HPP
#define RTLPSMRI_CORE_DIRECT_COIL_AWARE_NFT_HPP

namespace mri {

inline void
CoilAwareNFT2D
( const DistMatrix<Complex<double>,STAR,VR>& FHat, 
        DistMatrix<Complex<double>,STAR,VR>& F )
{
    DEBUG_ONLY(CallStackEntry cse("CoilAwareNFT2D"))
    const double pi = 4*El::Atan( 1. );
    const int width = FHat.Width();
    const int numNonUniform = NumNonUniformPoints();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
    const int numCoils = NumCoils();
    DEBUG_ONLY(
        if( width != numCoils*NumTimesteps() )
            LogicError("Invalid width");
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
        nfft_plan& plan = CoilPlan( t );
        const double* x = plan.x;
        Complex<double>* f = F.Buffer(0,jLoc);
        const Complex<double>* fHat = FHat.LockedBuffer(0,jLoc);
        for( int xi=0; xi<numNonUniform; ++xi )
        {
            const double x0 = x[2*xi+0];
            const double x1 = x[2*xi+1];
            for( int ki=0; ki<N0; ++ki )
            {
                const double k0 = -N0/2+ki;
                for( int kj=0; kj<N1; ++kj )
                {
                    const Complex<double> phiHat = fHat[kj+ki*N1];
                    const double k1 = -N1/2+kj;
                    const double theta = -2*pi*(x0*k0+x1*k1);
                    const double realPart = cos(theta);
                    const double imagPart = sin(theta);
                    f[xi] += Complex<double>(realPart,imagPart)*phiHat;
                }
            }
        }
    }
    const double scale = 1./Sqrt(1.*N0*N1);
    Scale( scale, F );
}

inline void
CoilAwareAdjointNFT2D
( const DistMatrix<Complex<double>,STAR,VR>& F,
        DistMatrix<Complex<double>,STAR,VR>& FHat )
{
    DEBUG_ONLY(CallStackEntry cse("CoilAwareAdjointNFT2D"))
    const double pi = 4*El::Atan( 1. );
    const int width = F.Width();
    const int numNonUniform = NumNonUniformPoints();
    const int N0 = FirstBandwidth();
    const int N1 = SecondBandwidth();
    const int numCoils = NumCoils();
    DEBUG_ONLY(
        if( width != numCoils*NumTimesteps() )
            LogicError("Invalid width");
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
        nfft_plan& plan = CoilPlan( t );
        const double* x = plan.x;
        const Complex<double>* f = F.LockedBuffer(0,jLoc);
        Complex<double>* fHat = FHat.Buffer(0,jLoc);
        for( int xi=0; xi<numNonUniform; ++xi )
        {
            const double x0 = x[2*xi+0];
            const double x1 = x[2*xi+1];
            const Complex<double> phi = f[xi];
            for( int ki=0; ki<N0; ++ki )
            {
                const double k0 = -N0/2+ki;
                for( int kj=0; kj<N1; ++kj )
                {
                    const double k1 = -N1/2+kj;
                    const double theta = 2*pi*(x0*k0+x1*k1);
                    const double realPart = cos(theta);
                    const double imagPart = sin(theta);
                    fHat[kj+ki*N1] += Complex<double>(realPart,imagPart)*phi;
                }
            }
        }
    }
    const double scale = 1./Sqrt(1.*N0*N1);
    Scale( scale, FHat );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_DIRECT_COIL_AWARE_NFT_HPP
