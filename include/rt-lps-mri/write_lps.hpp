/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_WRITELPS_HPP
#define RTLPSMRI_WRITELPS_HPP

namespace mri {

inline void
WriteLPS
( const DistMatrix<Complex<double>,VC,STAR>& L,
  const DistMatrix<Complex<double>,VC,STAR>& S,
  Int N0, Int N1,
  bool tv=true,
  elem::FileFormat format=elem::MATLAB_ASCII )
{
#ifndef RELEASE
    CallStackEntry cse("WriteLPS");
#endif
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
            os << "L-tv-" << t;
        else
            os << "L-temporal-" << t;
        Write( B, format, os.str() );
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
            os << "S-tv-" << t;
        else
            os << "S-temporal-" << t;
        Write( B, format, os.str() );
    }
}

} // namespace mri

#endif // ifndef RTLPSMRI_WRITELPS_HPP
