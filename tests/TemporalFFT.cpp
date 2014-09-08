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
        DistMatrix<double,STAR,STAR> paths;
        Uniform( paths, 2*nnu, nt, 0., 0.5 );
        InitializeCoilPlans( paths, nc, N0, N1, n0, n1, m );

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
