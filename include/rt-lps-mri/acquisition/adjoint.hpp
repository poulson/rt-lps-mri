/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_ACQUISITION_ADJOINT_HPP
#define RTLPSMRI_ACQUISITION_ADJOINT_HPP

namespace mri {

// Adjoint application of the acquisition operator. This maps
// k-space x (coil,time) -> image x time

namespace acquisition {

inline void
ScaleBySensitivities
( const DistMatrix<Complex<double>,STAR,VR>& F,
        DistMatrix<Complex<double>,STAR,VR>& scaledF )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::ScaleBySensitivities");
#endif
    // TODO
}

inline void
ContractionPrescaling( DistMatrix<Complex<double>,STAR,VR>& FHat )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::ContractionPrescaling");
#endif
    // TODO
}

inline void
CoilContraction
( const DistMatrix<Complex<double>,STAR,VR>& FHat,
        DistMatrix<Complex<double>,VC,STAR>& images )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::CoilContraction");
#endif
    // TODO
}

} // namespace acquisition

// NOTE: scaledF and FHat are large temporary matrices and are only arguments
//       in order to allow for avoiding allocating and freeing their memory.
//       scaledF is only needed in order to avoid overwriting F.
inline void
AdjointAcquisition
( const DistMatrix<Complex<double>,STAR,VR>& F, 
        DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& scaledF,
        DistMatrix<Complex<double>,STAR,VR>& FHat )
{
#ifndef RELEASE
    CallStackEntry cse("AdjointAcquisition");
#endif
    // Pre-scale the k-space data by the coil sensitivities
    acquisition::ScaleBySensitivities( F, scaledF );

    // Transform each k-space vector into the image domain
    CoilAwareAdjointNFFT2D( scaledF, FHat );
    const double M = FirstBandwidth()*SecondBandwidth();
    Scale( 1./Sqrt(M), FHat );

    // Perform a contraction over the coils with a weighting related to 
    // their sensitivities
    acquisition::ContractionPrescaling( FHat );
    acquisition::CoilContraction( FHat, images );
}

} // namespace mri

#endif // ifndef RTLPSMRI_ACQUISITION_ADJOINT_HPP
