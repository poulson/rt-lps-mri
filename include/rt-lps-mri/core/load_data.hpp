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

inline void
LoadData
( int numNonUniform, int numCoils, int numTimesteps,
  std::string filename, DistMatrix<Complex<double>,STAR,VR>& data )
{
    DEBUG_ONLY(CallStackEntry cse("LoadData"))
    const int m = numNonUniform;
    const int n = numCoils*numTimesteps;
    elem::read::BinaryFlat( data, m, n, filename );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADDATA_HPP
