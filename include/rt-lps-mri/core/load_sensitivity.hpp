/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_LOADSENSITIVITY_HPP
#define RTLPSMRI_CORE_LOADSENSITIVITY_HPP

namespace mri {

// TODO: Some sort of error checking for the size of the file?
inline void
LoadSensitivity
( int N0, int N1, int numCoils, 
  std::string filename, DistMatrix<Complex<double>,STAR,STAR>& sensitivity )
{
#ifndef RELEASE
    CallStackEntry cse("LoadSensitivity");
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
    const int numExpected = N0*N1*numCoils*2*sizeof(double);
    if( numBytes != numExpected )
    {
        std::ostringstream os;
        os << "File was " << numBytes << " instead of "
           << N0 << " x " << N1 << " x " << numCoils
           << " x 2*" << sizeof(double) << " = " << numExpected << std::endl;  
        RuntimeError( os.str() );
    }
    is.seekg( 0, std::ios::beg );
#ifndef RELEASE
    if( sensitivity.Grid().Rank() == 0 )
        std::cout << "Sensitivity file was " << numBytes << " as expected"
                  << std::endl;
#endif
    
    sensitivity.ResizeTo( N0*N1, numCoils, N0*N1 );
    is.read
    ( (char*)sensitivity.Buffer(), N0*N1*numCoils*2*sizeof(double) );
#ifndef RELEASE
    const double frobNorm = FrobeniusNorm( sensitivity );
    if( sensitivity.Grid().Rank() == 0 )
        std::cout << "  Frobenius norm = " << frobNorm << std::endl;
#endif
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADSENSITIVITY_HPP
