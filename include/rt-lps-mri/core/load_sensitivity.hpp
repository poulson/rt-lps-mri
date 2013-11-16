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
    
    sensitivity.ResizeTo( N0*N1, numCoils, N0*N1 );
    is.read
    ( (char*)sensitivity.Buffer(), N0*N1*numCoils*2*sizeof(double) );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADSENSITIVITY_HPP
