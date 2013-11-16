/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_LOADPATHS_HPP
#define RTLPSMRI_CORE_LOADPATHS_HPP

namespace mri {

// TODO: Some sort of error checking for the size of the file?
inline void
LoadPaths
( int numNonUniform, int numTimesteps, 
  std::string filename, DistMatrix<double,STAR,STAR>& paths )
{
#ifndef RELEASE
    CallStackEntry cse("LoadPaths");
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
    const int numExpected = 2*numNonUniform*numTimesteps*sizeof(double);
    if( numBytes != numExpected )
    {
        std::ostringstream os;
        os << "File was " << numBytes << " instead of 2*"
           << numNonUniform << " x " << numTimesteps
           << " x 2*" << sizeof(double) << " = " << numExpected << std::endl;
        RuntimeError( os.str() );
    }
    is.seekg( 0, std::ios::end );
#ifndef RELEASE
    if( paths.Grid().Rank() == 0 )
        std::cout << "Paths file was " << numBytes << " as expected"
                  << std::endl;
#endif
    
    paths.ResizeTo( 2*numNonUniform, numTimesteps, 2*numNonUniform );
    is.read
    ( (char*)paths.Buffer(), numNonUniform*numTimesteps*2*sizeof(double) );
#ifndef RELEASE
    const double frobNorm = FrobeniusNorm( paths );
    if( paths.Grid().Rank() == 0 )
        std::cout << "  Frobenius norm = " << frobNorm << std::endl;
#endif
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADPATHS_HPP
