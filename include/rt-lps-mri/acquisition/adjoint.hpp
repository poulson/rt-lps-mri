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
        const auto sensitivity = Sensitivity().LockedBuffer(0,time);
        for( int i=0; i<height; ++i )
            fImage[i] *= sensitivity[i];
    }
}

inline void
ContractionPrescaling( DistMatrix<Complex<double>,STAR,VR>& FHat )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::ContractionPrescaling");
#endif
    const int numCoils = NumCoils();
    const auto& sensitivity = Sensitivity();
    const auto& sensitivityScalings = SensitivityScalings();
    const int height = FHat.Height();
    const int localWidth = FHat.LocalWidth();
    const int rowShift = FHat.RowShift();
    const int rowStride = FHat.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int coil = j % numCoils; 
        auto fHat = FHat.Buffer(0,jLoc);
        const auto senseCol = sensitivity.LockedBuffer(0,coil);
        const auto senseScaleCol = sensitivityScalings.LockedBuffer();
        for( int i=0; i<height; ++i )
            fHat[i] *= Conj(senseCol[i])/senseScaleCol[i];    
    }
}

inline void
CoilContraction
( const DistMatrix<Complex<double>,STAR,VR>& FHat,
        DistMatrix<Complex<double>,VC,STAR>& images )
{
#ifndef RELEASE
    CallStackEntry cse("acquisition::CoilContraction");
#endif
    const int height = FHat.Height();
    const int numCoils = NumCoils();
    const int numTimesteps = NumTimesteps();
    const int rowAlign = FHat.RowAlign();
    const int rowShift = FHat.RowShift();
    const int rowStride = FHat.RowStride();

    DistMatrix<Complex<double>,STAR,VR> images_STAR_VR( FHat.Grid() );
    images_STAR_VR.AlignWith( FHat );
    Zeros( images_STAR_VR, height, numTimesteps );

    const int localNewWidth = images_STAR_VR.LocalWidth();
    const int localWidth = FHat.LocalWidth();

    // Determine the number of entries we send to each process
    std::vector<int> sendSizes( rowStride, 0 );
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int jNew = j / numCoils;
        const int newOwner = (jNew+rowAlign) % rowStride;
        sendSizes[newOwner] += height;
    }

    // Determine the number of entries we recv from each process
    std::vector<int> recvSizes( rowStride, 0 );
    for( int jLoc=0; jLoc<localNewWidth; ++jLoc )
    {
        const int jNew = rowShift + jLoc*rowStride;
        for( int j=jNew*numCoils; j<(jNew+1)*numCoils; ++j )
        {
            const int owner = (j+rowAlign) % rowStride;
            recvSizes[owner] += height;
        }
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
    std::vector<int> offsets = sendOffsets;
    std::vector<Complex<double>> sendBuf( totalSend );
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int jNew = j / numCoils;
        const int newOwner = (jNew+rowAlign) % rowStride;
        elem::MemCopy
        ( &sendBuf[offsets[newOwner]], 
          FHat.LockedBuffer(0,jLoc), height );
        offsets[newOwner] += height;
    }

    // Perform the non-uniform AllToAll communication
    std::vector<Complex<double>> recvBuf( totalRecv );
    mpi::AllToAll
    ( sendBuf.data(), sendSizes.data(), sendOffsets.data(),
      recvBuf.data(), recvSizes.data(), recvOffsets.data(),
      FHat.RowComm() );
    std::vector<Complex<double>>().swap( sendBuf );
    std::vector<int>().swap( sendSizes );
    std::vector<int>().swap( sendOffsets );

    // Unpack the recv buffer
    offsets = recvOffsets;
    for( int jLoc=0; jLoc<localNewWidth; ++jLoc )
    {
        const int jNew = rowShift + jLoc*rowStride;
        auto imageCol = images_STAR_VR.Buffer(0,jLoc);
        for( int j=jNew*numCoils; j<(jNew+1)*numCoils; ++j )
        {
            const int owner = (j+rowAlign) % rowStride;
            for( int i=0; i<height; ++i )
                imageCol[i] += recvBuf[offsets[owner]++];
        }
    }

    images = images_STAR_VR;
}

} // namespace acquisition

inline void
AdjointAcquisition
( const DistMatrix<Complex<double>,STAR,VR>& F, 
        DistMatrix<Complex<double>,VC,STAR>& images )
{
#ifndef RELEASE
    CallStackEntry cse("AdjointAcquisition");
#endif
    // Pre-scale the k-space data by the coil sensitivities
    DistMatrix<Complex<double>,STAR,VR> scaledF( F.Grid() );
    acquisition::ScaleBySensitivities( F, scaledF );

    // Transform each k-space vector into the image domain
    DistMatrix<Complex<double>,STAR,VR> FHat( F.Grid() );
    CoilAwareAdjointNFFT2D( scaledF, FHat );

    // Perform a contraction over the coils with a weighting related to 
    // their sensitivities
    acquisition::ContractionPrescaling( FHat );
    acquisition::CoilContraction( FHat, images );
}

} // namespace mri

#endif // ifndef RTLPSMRI_ACQUISITION_ADJOINT_HPP
