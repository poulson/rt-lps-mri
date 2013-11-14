/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_LPS_HPP
#define RTLPSMRI_LPS_HPP

namespace mri {

inline void
LPS
( const DistMatrix<Complex<double>,STAR,VR>& D,
  const DistMatrix<Complex<double>,VC,STAR>& M,
        DistMatrix<Complex<double>,VC,STAR>& L,
        DistMatrix<Complex<double>,VC,STAR>& S )
{
#ifndef RELEASE
    CallStackEntry cse("LPS");
#endif
    // TODO
}

} // namespace mri

#endif // ifndef RTLPSMRI_LPS_HPP
