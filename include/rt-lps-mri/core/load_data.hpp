/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_LOADDATA_HPP
#define RTLPSMRI_CORE_LOADDATA_HPP

namespace mri {

// TODO: Some sort of error checking for the size of the file?
inline void
LoadData
( int numNonUniform, int numCoils, int numTimesteps,
  std::string filename, DistMatrix<Complex<double>,STAR,VR>& data )
{
#ifndef RELEASE
    CallStackEntry cse("LoadData");
#endif
    std::ifstream is;
    is.open( filename.c_str(), std::ios::in|std::ios::binary );
    if( !is.is_open() )
    {
        std::ostringstream os;
        os << "Could not open " << filename;
        RuntimeError( os.str() );
    }

    // Make sure the file is of the right size
    is.seekg( 0, std::ios::end ); 
    const int numBytes = is.tellg();
    const int numExpected = 
        numNonUniform*numCoils*numTimesteps*2*sizeof(double);
    if( numBytes != numExpected )
    {
        std::ostringstream os;
        os << "File was " << numBytes << " instead of "
           << numNonUniform << " x " << numCoils << " x " << numTimesteps
           << " x 2*" << sizeof(double) << " = " << numExpected << std::endl;
        RuntimeError( os.str() );
    }
    is.seekg( 0, std::ios::beg );
#ifndef RELEASE
    if( data.Grid().Rank() == 0 )
        std::cout << "Data file was " << numBytes << " as expected"
                  << std::endl;
#endif

    data.ResizeTo( numNonUniform, numCoils*numTimesteps, numNonUniform );
    const int rowShift = data.RowShift();
    const int rowStride = data.RowStride();
    const int localWidth = data.LocalWidth();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        const int localIndex = j*numNonUniform;
        const std::streamoff pos = localIndex*2*sizeof(double);
        is.seekg( pos );
        is.read( (char*)data.Buffer(0,jLoc), numNonUniform*2*sizeof(double) );
    }
#ifndef RELEASE
    const double frobNorm = FrobeniusNorm( data );
    if( data.Grid().Rank() == 0 )
        std::cout << "  Frobenius norm = " << frobNorm << std::endl;
#endif
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADDATA_HPP
