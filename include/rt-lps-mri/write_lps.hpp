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
#ifndef RTLPSMRI_WRITELPS_HPP
#define RTLPSMRI_WRITELPS_HPP

namespace mri {

inline void
WriteLPS
( const DistMatrix<Complex<double>,VC,STAR>& L,
  const DistMatrix<Complex<double>,VC,STAR>& S,
  Int N0, Int N1, Int plane,
  bool tv=true,
  FileFormat format=ASCII_MATLAB )
{
    DEBUG_ONLY(CallStackEntry cse("WriteLPS"))
    Matrix<Complex<double>> B( N0, N1, N0 );
    Complex<double>* BBuf = B.Buffer();

    // Split the timesteps among processes
    DistMatrix<Complex<double>,STAR,VR> A_STAR_VR( L );
    const int rowShift = A_STAR_VR.RowShift();
    const int rowStride = A_STAR_VR.RowStride();
    const int localWidth = A_STAR_VR.LocalWidth();
    for( Int tLoc=0; tLoc<localWidth; ++tLoc )
    {
        const Int t = rowShift + tLoc*rowStride;
        const Complex<double>* LBuf = A_STAR_VR.LockedBuffer(0,tLoc);
        for( Int j=0; j<N1; ++j )
            for( Int i=0; i<N0; ++i )    
                BBuf[i+j*N0] = LBuf[j+i*N1];

        std::ostringstream os;
        if( tv )
            os << "L-" << plane << "-tv-" << t;
        else
            os << "L-" << plane << "-temporal-" << t;
        Write( B, os.str(), format );
    }
    A_STAR_VR = S;
    for( Int tLoc=0; tLoc<localWidth; ++tLoc )
    {
        const Int t = rowShift + tLoc*rowStride;
        const Complex<double>* SBuf = A_STAR_VR.LockedBuffer(0,tLoc);
        for( Int j=0; j<N1; ++j )
            for( Int i=0; i<N0; ++i )    
                BBuf[i+j*N0] = SBuf[j+i*N1];

        std::ostringstream os;
        if( tv )
            os << "S-" << plane << "-tv-" << t;
        else
            os << "S-" << plane << "-temporal-" << t;
        Write( B, os.str(), format );
    }
}

} // namespace mri

#endif // ifndef RTLPSMRI_WRITELPS_HPP
