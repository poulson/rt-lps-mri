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

inline void
LoadDensity
( int numNonUniform, int numTimesteps, 
  std::string filename, DistMatrix<double,STAR,STAR>& density )
{
    DEBUG_ONLY(CallStackEntry cse("LoadDensity"))
    const int m = numNonUniform;
    const int n = numTimesteps;
    density.Resize( m, n );
    El::Read( density, filename, BINARY_FLAT );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADDENSITY_HPP
