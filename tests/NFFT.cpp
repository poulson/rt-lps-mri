/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#include "rt-lps-mri.hpp"
using namespace mri;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const int N1 = Input("--N1","bandwidth in x direction",5);
        const int N2 = Input("--N2","bandwidth in y direction",5);
        const int M  = Input("--M","number of non-uniform nodes",25);
        const int n1 = Input("--n1","FFT size in x direction",16);
        const int n2 = Input("--n2","FFT size in y direction",16);
        const int m = Input("--m","cutoff parameter",2);
        const int width = Input("--width","number of indep. vectors",10);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double,         STAR,VR> X;
        DistMatrix<Complex<double>,STAR,VR> F, FHat;

        Uniform(X,2*M,width,0.,0.5);
        Sort(X); 
        Uniform(FHat,N1*N2,width);
        Zeros(F,M,width);

        Print( X, "X" );
        Print( FHat, "FHat" );

        NFFT2D( N1, N2, M, n1, n2, m, FHat, X, F );
        Print( F, "F after forward" );

        AdjointNFFT2D( N1, N2, M, n1, n2, m, FHat, X, F );
        Print( FHat, "FHat after adjoint" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
