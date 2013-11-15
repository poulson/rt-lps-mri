/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#include "rt-lps-mri.hpp"
using namespace mri;

typedef double Real;
typedef Complex<Real> F;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const int nc = Input("--nc","number of coils",16);
        const int nt = Input("--nt","number of timesteps",10);
        const int N0 = Input("--N0","bandwidth in x direction",6);
        const int N1 = Input("--N1","bandwidth in y direction",6);
        const int nnu = Input("--nnu","number of non-uniform nodes",36);
        const int n0 = Input("--n0","FFT size in x direction",16);
        const int n1 = Input("--n1","FFT size in y direction",16);
        const int m = Input("--m","cutoff parameter",2);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Initialize the coil plans
        DistMatrix<double,STAR,VR> X;
        Uniform( X, 2*nnu, nc*nt, 0., 0.5 );
        Sort( X );
        InitializeCoilPlans( X, nc, nt, N0, N1, n0, n1, m );

        // Generate a random matrix to apply the temporal FFT's to
        DistMatrix<F,VC,STAR> A;
        Uniform( A, N0*N1, nt );
        DistMatrix<F,VC,STAR> ACopy( A );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );

        // Apply the forward FFT
        TemporalFFT( A );
        const Real frobFA = FrobeniusNorm( A );
        if( print )
            Print( A, "F(A)" );
        if( display )
            Display( A, "F(A)" );

        // Apply the backward FFT
        TemporalAdjointFFT( A );
        const Real frobFAdjFA = FrobeniusNorm( A );
        if( print )
            Print( A, "F'(F(A))" );
        if( display )
            Display( A, "F'(F(A))" );

        Axpy( F(-1), ACopy, A );
        const Real frobE = FrobeniusNorm( A );

        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n"
                      << "|| F(A) ||_F = " << frobFA << "\n"
                      << "|| F'(F(A)) ||_F = " << frobFAdjFA << "\n"
                      << "|| F'(F(A))-A ||_F = " << frobE << "\n"
                      << "|| F'(F(A))-A ||_F / || A ||_F = " << frobE/frobA 
                      << "\n" << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
