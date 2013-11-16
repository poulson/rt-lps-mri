/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_LOADDENSITY_HPP
#define RTLPSMRI_CORE_LOADDENSITY_HPP

namespace mri {

// TODO: Some sort of error checking for the size of the file?
inline void
LoadDensity
( int numNonUniform, int numTimesteps, 
  std::string filename, DistMatrix<double,STAR,STAR>& density )
{
#ifndef RELEASE
    CallStackEntry cse("LoadDensity");
#endif
    std::ifstream is;
    is.open( filename.c_str(), std::ios::in|std::ios::binary );
    if( !is.is_open() )
    {
        std::ostringstream os;
        os << "Could not open " << filename;
        RuntimeError( os.str() );
    }
    
    density.ResizeTo( numNonUniform, numTimesteps, numNonUniform );
    is.read
    ( (char*)density.Buffer(), numNonUniform*numTimesteps*sizeof(double) );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADDENSITY_HPP
