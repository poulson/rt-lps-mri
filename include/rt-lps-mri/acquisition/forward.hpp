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

// TODO: Exploit redundancy in coil data to reduce amount of communication
inline void
Scatter
( const DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& scatteredImages,
  bool progress=false )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::Scatter");
#endif
    const int height = images.Height();
    const int localHeight = images.LocalHeight();
    const int numCoils = NumCoils();
    const int numTimesteps = NumTimesteps();

    elem::Timer timer("");
    timer.Start();
    DistMatrix<Complex<double>,VC,STAR> 
        scatteredImages_VC_STAR( images.Grid() );
    scatteredImages_VC_STAR.AlignWith( images );
    Zeros( scatteredImages_VC_STAR, height, numCoils*numTimesteps );
    for( int t=0; t<numTimesteps; ++t )
    {
        for( int coil=0; coil<numCoils; ++coil )
        {
            const int jNew = coil + t*numCoils;
            elem::MemCopy
            ( scatteredImages_VC_STAR.Buffer(0,jNew),
              images.LockedBuffer(0,t), localHeight ); 
        }
    }
    const double memCopyTime = timer.Stop();

    timer.Start();
    scatteredImages = scatteredImages_VC_STAR;
    const double redistTime = timer.Stop();

    if( progress && images.Grid().Rank() == 0 )
        std::cout << "      Scatter copies: " << memCopyTime << " seconds\n"
                  << "      Scatter redist: " << redistTime << " seconds"
                  << std::endl;
}

inline void
ScaleBySensitivities( DistMatrix<Complex<double>,STAR,VR>& scatteredImages )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::ScaleBySensitivities");
#endif
    const int numCoils = NumCoils();
    const int height = scatteredImages.Height();
    const int localWidth = scatteredImages.LocalWidth();
    const int rowShift = scatteredImages.RowShift();
    const int rowStride = scatteredImages.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int coil = j % numCoils; // TODO: use mapping jLoc -> coil?
        auto image = scatteredImages.Buffer(0,jLoc);
        const auto senseCol = Sensitivity().LockedBuffer(0,coil);
        for( int i=0; i<height; ++i )
            image[i] *= senseCol[i];
    }
}

} // namespace acquisition

inline void
Acquisition
( const DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& F,
  bool progress=false )
{
#ifndef RELEASE
    CallStackEntry cse("Acquisition");
#endif
    // Redundantly scatter image x time -> image x (coil,time)
    elem::Timer timer("");
    timer.Start();
    DistMatrix<Complex<double>,STAR,VR> scatteredImages( images.Grid() );
    acquisition::Scatter( images, scatteredImages, progress );
    const double scatterTime = timer.Stop();

    // Scale by the coil sensitivities
    timer.Start();
    acquisition::ScaleBySensitivities( scatteredImages ); 
    const double scaleTime = timer.Stop();

    // Finish the transformation
    timer.Start();
    CoilAwareNFFT2D( scatteredImages, F );
    const double nfftTime = timer.Stop();

    if( progress && F.Grid().Rank() == 0 )
        std::cout << "    scatter: " << scatterTime << " seconds\n"
                  << "    scale:   " << scaleTime << " seconds\n"
                  << "    NFFT:    " << nfftTime << " seconds\n"
                  << std::endl;
}

} // namespace mri

#endif // ifndef RTLPSMRI_ACQUISITION_FORWARD_HPP
