/*
   Copyright (c) 2013-2014, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI).

   RT-LPS-MRI is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   RT-LPS-MRI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with RT-LPS-MRI.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once
#ifndef RTLPSMRI_HPP
#define RTLPSMRI_HPP

#include "El.hpp"

#include "fftw3.h"
#include "nfft3.h"

#include "rt-lps-mri/config.h"

// The core of the library
#include "rt-lps-mri/core/environment_decl.hpp"
#include "rt-lps-mri/core/environment_impl.hpp"
#include "rt-lps-mri/core/nfft.hpp"
#include "rt-lps-mri/core/nft.hpp"
#include "rt-lps-mri/core/coil_aware_nfft.hpp"
#include "rt-lps-mri/core/coil_aware_nft.hpp"
#include "rt-lps-mri/core/temporal_fft.hpp"
#include "rt-lps-mri/core/load_data.hpp"
#include "rt-lps-mri/core/load_density.hpp"
#include "rt-lps-mri/core/load_paths.hpp"
#include "rt-lps-mri/core/load_sensitivity.hpp"

// Applying the acquisition operator and its adjoint
#include "rt-lps-mri/acquisition/forward.hpp"
#include "rt-lps-mri/acquisition/adjoint.hpp"

#include "rt-lps-mri/lps.hpp"
#include "rt-lps-mri/write_lps.hpp"

#endif // ifndef RTLPSMRI_HPP
