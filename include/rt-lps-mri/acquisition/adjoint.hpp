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
#ifndef RTLPSMRI_ACQUISITION_ADJOINT_HPP
#define RTLPSMRI_ACQUISITION_ADJOINT_HPP

namespace mri {

// Adjoint application of the acquisition operator. This maps
// k-space x (coil,time) -> image x time

namespace acquisition {

inline void
ScaleByDensities
( const DistMatrix<Complex<double>,STAR,VR>& F,
        DistMatrix<Complex<double>,STAR,VR>& scaledF )
{
    DEBUG_ONLY(CallStackEntry cse("acquisition::ScaleByDensities"))
    const int numCoils = NumCoils();
    const int height = F.Height();
    const int localWidth = F.LocalWidth();
    const int rowShift = F.RowShift();
    const int rowStride = F.RowStride();
    scaledF = F;
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int time = j / numCoils; // TODO: use mapping jLoc -> time?
        auto fImage = scaledF.Buffer(0,jLoc);
        const auto density = DensityComp().LockedBuffer(0,time);
        for( int i=0; i<height; ++i )
            fImage[i] *= density[i];
    }
}

inline void
ContractionPrescaling( DistMatrix<Complex<double>,STAR,VR>& FHat )
{
    DEBUG_ONLY(CallStackEntry cse("acquisition::ContractionPrescaling"))
    const int numCoils = NumCoils();
    const int height = FHat.Height();
    const int localWidth = FHat.LocalWidth();
    const int rowShift = FHat.RowShift();
    const int rowStride = FHat.RowStride();
    const auto& sensitivity = Sensitivity();
    const auto senseScaleCol = SensitivityScalings().LockedBuffer();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int coil = j % numCoils; 
        auto fHat = FHat.Buffer(0,jLoc);
        const auto senseCol = sensitivity.LockedBuffer(0,coil);
        for( int i=0; i<height; ++i )
            fHat[i] *= Conj(senseCol[i])/senseScaleCol[i];    
    }
}

inline void
CoilContraction
( const DistMatrix<Complex<double>,STAR,VR>& FHat,
        DistMatrix<Complex<double>,VC,STAR>& images,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("acquisition::CoilContraction"))
    typedef Complex<double> F;
    const int height = FHat.Height();
    const int numCoils = NumCoils();
    const int numTimesteps = NumTimesteps();

    El::Timer timer("");
    timer.Start();
    DistMatrix<Complex<double>,VC,STAR> FHat_VC_STAR( images.Grid() );
    FHat_VC_STAR.AlignWith( images ); 
    FHat_VC_STAR = FHat;
    const double redistTime = timer.Stop();

    timer.Start();
    Zeros( images, height, numTimesteps );
    const int localHeight = images.LocalHeight();
    for( int t=0; t<numTimesteps; ++t )
    {
        for( int coil=0; coil<numCoils; ++coil )
        {
            const int jOld = coil + t*numCoils;
            blas::Axpy
            ( localHeight, F(1), 
              FHat_VC_STAR.LockedBuffer(0,jOld), 1, images.Buffer(0,t), 1 );
        }
    }
    const double axpyTime = timer.Stop();
    if( progress && FHat.Grid().Rank() == 0 )
        std::cout << "      Contract redist: " << redistTime << " seconds\n"
                  << "      Contract axpies: " << axpyTime << " seconds"
                  << std::endl;
}

} // namespace acquisition

inline void
AdjointAcquisition
( const DistMatrix<Complex<double>,STAR,VR>& F, 
        DistMatrix<Complex<double>,VC,STAR>& images, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("AdjointAcquisition"))
    // Pre-scale the k-space data by the coil sensitivities
    El::Timer timer("");
    timer.Start();
    DistMatrix<Complex<double>,STAR,VR> scaledF( F.Grid() );
    acquisition::ScaleByDensities( F, scaledF );
    const double scaleTime = timer.Stop();

    // Transform each k-space vector into the image domain
    timer.Start();
    DistMatrix<Complex<double>,STAR,VR> FHat( F.Grid() );
    CoilAwareAdjointNFFT2D( scaledF, FHat );
    const double adjNfftTime = timer.Stop();

    // Perform a contraction over the coils with a weighting related to 
    // their sensitivities
    timer.Start();
    acquisition::ContractionPrescaling( FHat );
    const double prescaleTime = timer.Stop();
    timer.Start();
    acquisition::CoilContraction( FHat, images, progress );
    const double contractTime = timer.Stop();

    if( progress && F.Grid().Rank() == 0 ) 
        std::cout << "    scale:    " << scaleTime << " seconds\n"
                  << "    adjNFFT:  " << adjNfftTime << " seconds\n"
                  << "    prescale: " << prescaleTime << " seconds\n"
                  << "    contract: " << contractTime << " seconds\n"
                  << std::endl;
}

} // namespace mri

#endif // ifndef RTLPSMRI_ACQUISITION_ADJOINT_HPP
