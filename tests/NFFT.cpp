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
        const int N0 = Input("--N0","bandwidth in x direction",6);
        const int N1 = Input("--N1","bandwidth in y direction",6);
        const int M  = Input("--M","number of non-uniform nodes",36);
        const int n0 = Input("--n0","FFT size in x direction",16);
        const int n1 = Input("--n1","FFT size in y direction",16);
        const int m = Input("--m","cutoff parameter",2);
        const int width = Input("--width","number of indep. vectors",10);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double,         STAR,VR> X;
        DistMatrix<Complex<double>,STAR,VR> F, FDirect, FHat, FHatDirect;

        Uniform(X,2*M,width,0.,0.5);
        Sort(X); 
        Uniform(FHat,N0*N1,width);
        if( print )
        {
            Print( X, "X" );
            Print( FHat, "FHat" );
        }

        NFFT2D( N0, N1, M, n0, n1, m, FHat, X, F );
        if( print )
            Print( F, "F after forward" );

        DirectNFT2D( N0, N1, M, FHat, X, FDirect );
        if( print )
            Print( FDirect, "F after direct forward" );

        AdjointNFFT2D( N0, N1, M, n0, n1, m, FHat, X, F );
        if( print )
            Print( FHat, "FHat after adjoint" );

        DirectAdjointNFT2D( N0, N1, M, FHatDirect, X, F );
        if( print )
            Print( FHatDirect, "FHat after direct adjoint" );

        const double frobFDir = FrobeniusNorm( FDirect );
        const double frobFHatDir = FrobeniusNorm( FHatDirect );
        Axpy( -1., F, FDirect );
        Axpy( -1., FHat, FHatDirect );
        const double frobE = FrobeniusNorm( FDirect );
        const double frobEHat = FrobeniusNorm( FHatDirect );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| F ||_F = " << frobFDir << "\n"
                      << "|| E ||_F = " << frobE << "\n"
                      << "|| E ||_F / || F ||_F = " << frobE/frobFDir <<"\n"
                      << "\n"
                      << "|| FHat ||_F = " << frobFHatDir << "\n"
                      << "|| EHat ||_F = " << frobEHat << "\n"
                      << "|| EHat ||_F / || FHat ||_F = " 
                      << frobEHat/frobFHatDir << "\n"
                      << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
