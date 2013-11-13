/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef RTLPSMRI_CORE_DIRECT_NFT_HPP
#define RTLPSMRI_CORE_DIRECT_NFT_HPP

namespace mri {

inline void
DirectNFT2D
( int N0, int N1, int M, 
  const DistMatrix<Complex<double>,STAR,VR>& FHat, 
  const DistMatrix<double,         STAR,VR>& X,
        DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("DirectNFT2D");
#endif
    const double pi = 4*elem::Atan( 1. );
    const Int d = 2;
    const Int width = X.Width();
#ifndef RELEASE
    if( FHat.Height() != N0*N1 )
        LogicError("Invalid FHat height");
    if( X.Height() != d*M )
        LogicError("Invalid X height");
    if( FHat.Width() != X.Width() )
        LogicError("FHat and X must have the same width");
    if( FHat.RowAlign() != X.RowAlign() )
        LogicError("FHat and X are not aligned");
#endif
    F.AlignWith( FHat );
    Zeros( F, M, width );
    const Int locWidth = F.LocalWidth();
    for( Int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        Complex<double>* FCol = F.Buffer(0,jLoc);
        const double* XCol = X.LockedBuffer(0,jLoc);
        const Complex<double>* FHatCol = FHat.LockedBuffer(0,jLoc);
        for( Int xi=0; xi<M; ++xi )
        {
            const double x0 = XCol[2*xi+0];
            const double x1 = XCol[2*xi+1];
            for( Int ki=0; ki<N0; ++ki )
            {
                const double k0 = -N0/2+ki;
                for( Int kj=0; kj<N1; ++kj )
                {
                    const Complex<double> fhat = FHatCol[kj+ki*N1];
                    const double k1 = -N1/2+kj;
                    const double theta = -2*pi*(x0*k0+x1*k1);
                    const double realPart = cos(theta);
                    const double imagPart = sin(theta);
                    FCol[xi] += Complex<double>(realPart,imagPart)*fhat;
                }
            }
        }
    }
}

inline void
DirectAdjointNFT2D
( int N0, int N1, int M, 
        DistMatrix<Complex<double>,STAR,VR>& FHat,
  const DistMatrix<double,         STAR,VR>& X,
  const DistMatrix<Complex<double>,STAR,VR>& F )
{
#ifndef RELEASE
    CallStackEntry cse("DirectAdjointNFT2D");
#endif
    const double pi = 4*elem::Atan( 1. );
    const Int d = 2;
    const Int width = X.Width();
#ifndef RELEASE
    if( F.Height() != M )
        LogicError("Invalid F height");
    if( X.Height() != d*M )
        LogicError("Invalid X height");
    if( F.Width() != X.Width() )
        LogicError("F and X must have the same width");
    if( F.RowAlign() != X.RowAlign() )
        LogicError("F and X are not aligned");
#endif
    FHat.AlignWith( F );
    Zeros( FHat, N0*N1, width );
    const Int locWidth = FHat.LocalWidth();

    for( Int jLoc=0; jLoc<locWidth; ++jLoc )
    {
        const double* XCol = X.LockedBuffer(0,jLoc);
        const Complex<double>* FCol = F.LockedBuffer(0,jLoc);
        Complex<double>* FHatCol = FHat.Buffer(0,jLoc);
        for( Int xi=0; xi<M; ++xi )
        {
            const double x0 = XCol[2*xi+0];
            const double x1 = XCol[2*xi+1];
            const Complex<double> f = FCol[xi];
            for( Int ki=0; ki<N0; ++ki )
            {
                const double k0 = -N0/2+ki;
                for( Int kj=0; kj<N1; ++kj )
                {
                    const double k1 = -N1/2+kj;
                    const double theta = 2*pi*(x0*k0+x1*k1);
                    const double realPart = cos(theta);
                    const double imagPart = sin(theta);
                    FHatCol[kj+ki*N1] += Complex<double>(realPart,imagPart)*f;
                }
            }
        }

    }
}

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_DIRECT_NFT_HPP
