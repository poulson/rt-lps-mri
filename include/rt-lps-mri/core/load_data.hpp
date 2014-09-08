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
    data.Resize( m, n );
    El::Read( data, filename, BINARY_FLAT );
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_LOADDATA_HPP
