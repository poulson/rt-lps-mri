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

inline void
LoadPaths
( int numNonUniform, int numTimesteps, 
  std::string filename, DistMatrix<double,STAR,STAR>& paths )
{
    DEBUG_ONLY(CallStackEntry cse("LoadPaths"))
    const int m = 2*numNonUniform;
    const int n = numTimesteps;
    El::read::BinaryFlat( paths, m, n, filename );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADPATHS_HPP
