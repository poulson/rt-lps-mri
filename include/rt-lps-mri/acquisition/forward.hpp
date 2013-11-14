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

// TODO: Modify this routine to avoid sending and receiving redundant columns
inline void
Scatter
( const DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& scatteredImages )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::Scatter");
#endif
    DistMatrix<Complex<double>,STAR,VR> images_STAR_VR( images );

    const int height = images.Height();
    const int numCoils = NumCoils();
    const int numTimesteps = NumTimesteps();
    const int rowAlign = images_STAR_VR.RowAlign();
    const int rowShift = images_STAR_VR.RowShift();
    const int rowStride = images_STAR_VR.RowStride();

    scatteredImages.SetGrid( images_STAR_VR.Grid() );
    scatteredImages.AlignWith( images_STAR_VR );
    Zeros( scatteredImages, height, numCoils*numTimesteps );

    const int localOrigWidth = images_STAR_VR.LocalWidth();
    const int localWidth = scatteredImages.LocalWidth();

    // Determine the number of entries we send to each process
    std::vector<int> sendSizes( rowStride, 0 );
    for( int jLoc=0; jLoc<localOrigWidth; ++jLoc )
    {
        const int jOrig = rowShift + jLoc*rowStride;
        for( int j=jOrig*numCoils; j<(jOrig+1)*numCoils; ++j )
        {
            const int owner = (j+rowAlign) % rowStride;
            sendSizes[owner] += height;
        }
    }

    // Determine the number of entries we receive from each process
    std::vector<int> recvSizes( rowStride, 0 );
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int jOrig = j / numCoils;
        const int origOwner = (jOrig+rowAlign) % rowStride;
        recvSizes[origOwner] += height;
    }

    // Form the send and recv offset vectors
    int totalSend=0, totalRecv=0;
    std::vector<int> sendOffsets( rowStride ), recvOffsets( rowStride );
    for( int q=0; q<rowStride; ++q )
    {
        sendOffsets[q] = totalSend;
        recvOffsets[q] = totalRecv;
        totalSend += sendSizes[q];
        totalRecv += recvSizes[q];
    }

    // Pack the send buffer
    std::vector<Complex<double>> sendBuf( totalSend );
    std::vector<int> offsets = sendOffsets;
    for( int jLoc=0; jLoc<localOrigWidth; ++jLoc )
    {
        const int jOrig = rowShift + jLoc*rowStride;
        for( int j=jOrig*numCoils; j<(jOrig+1)*numCoils; ++j )
        {
            const int owner = (j+rowAlign) % rowStride;
            elem::MemCopy
            ( &sendBuf[offsets[owner]], 
              images_STAR_VR.LockedBuffer(0,jLoc), height );
            offsets[owner] += height;
        }
    }

    // Perform the non-uniform AllToAll communication
    std::vector<Complex<double>> recvBuf( totalRecv );
    mpi::AllToAll
    ( sendBuf.data(), sendSizes.data(), sendOffsets.data(), 
      recvBuf.data(), recvSizes.data(), recvOffsets.data(), 
      images_STAR_VR.RowComm() );
    std::vector<Complex<double>>().swap( sendBuf );
    std::vector<int>().swap( sendSizes );
    std::vector<int>().swap( sendOffsets );

    // Unpack the recv buffer
    offsets = recvOffsets;
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int jOrig = j / numCoils;
        const int origOwner = (jOrig+rowAlign) % rowStride;
        elem::MemCopy
        ( scatteredImages.Buffer(0,jLoc),
          &recvBuf[offsets[origOwner]], height );
        offsets[origOwner] += height;
    }
}

inline void
CompensateDensities
( DistMatrix<Complex<double>,STAR,VR>& scatteredImages )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::CompensateDensities");
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
        const auto densities = DensityComp().LockedBuffer(0,coil);
        for( int i=0; i<height; ++i )
            image[i] *= densities[i];
    }
}

} // namespace acquisition

inline void
Acquisition
( const DistMatrix<Complex<double>,VC,STAR>& images,
        DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("Acquisition");
#endif
    // Redundantly scatter image x time -> image x (coil,time)
    DistMatrix<Complex<double>,STAR,VR> scatteredImages( images.Grid() );
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
