/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_ACQUISITION_FORWARD_HPP
#define RTLPSMRI_ACQUISITION_FORWARD_HPP

namespace mri {

// Forward application of the acquisition operator. This maps
// image x time -> k-space x (coil,time)
//
// The image domain is organized in a row-major fashion due to our usage of
// NFFT, and so, for (i,j) in [0,N0) x [0,N1), pixel (i,j) corresponds to 
// index j + i*N1.
//
// On the other hand, in order to heuristically balance the communication
// required for the awkward redistributions involved during the transformations
// between the domains, the (coil,time) pair (c,t) in [0,nc) x [0,nt) is stored
// at index c + t*nc. This might eventually change.

namespace acquisition {

inline void
Scatter
( const DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& scatteredImages )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::Scatter");
#endif
    DistMatrix<Complex<double>,STAR,MR> images_STAR_MR( images );
    // TODO
}

inline void
CompensateDensities
( DistMatrix<Complex<double>,STAR,VR>& scatteredImages )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::CompensateDensities");
#endif
    // TODO
}

} // namespace acquisition

// NOTE: scatteredImages is a large temporary matrix and is only an argument
//       in order to allow for avoiding allocating and freeing its memory
inline void
Acquisition
( const DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& F,
        DistMatrix<Complex<double>,STAR,VR>& scatteredImages )
{
#ifndef RELEASE
    CallStackEntry cse("Acquisition");
#endif
    // Redundantly scatter image x time -> image x (coil,time)
    acquisition::Scatter( images, scatteredImages );

    // Apply the density compensations
    acquisition::CompensateDensities( scatteredImages ); 

    // Finish the transformation
    CoilAwareNFFT2D( scatteredImages, F );
    const double M = FirstBandwidth()*SecondBandwidth();
    Scale( 1./Sqrt(M), F );
}

} // namespace mri

#endif // ifndef RTLPSMRI_ACQUISITION_FORWARD_HPP
