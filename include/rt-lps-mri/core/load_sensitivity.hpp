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

inline void
LoadSensitivity
( int N0, int N1, int numCoils, 
  std::string filename, DistMatrix<Complex<double>,STAR,STAR>& sensitivity )
{
    DEBUG_ONLY(CallStackEntry cse("LoadSensitivity"))
    const int m = N0*N1;
    const int n = numCoils;
    El::read::BinaryFlat( sensitivity, m, n, filename );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADSENSITIVITY_HPP
