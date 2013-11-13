/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_HPP
#define RTLPSMRI_HPP

// TODO: Replace with elemental-lite.hpp and supporting routines
#include "elemental.hpp"

#include "nfft3.h"

#include "rt-lps-mri/config.h"

//
// The core of the library
//

#include "rt-lps-mri/core/environment_decl.hpp"
#include "rt-lps-mri/core/environment_impl.hpp"
#include "rt-lps-mri/core/nfft.hpp"
#include "rt-lps-mri/core/direct_nft.hpp"

#endif // ifndef RTLPSMRI_HPP
