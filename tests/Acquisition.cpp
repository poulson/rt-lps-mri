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
        const int nc = Input("--nc","number of coils",16);
        const int nt = Input("--nt","number of timesteps",10);
        const int N0 = Input("--N0","bandwidth in x direction",6);
        const int N1 = Input("--N1","bandwidth in y direction",6);
        const int nnu  = Input("--nnu","number of non-uniform nodes",36);
        const int n0 = Input("--n0","FFT size in x direction",16);
        const int n1 = Input("--n1","FFT size in y direction",16);
        const int m = Input("--m","cutoff parameter",2);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Sample from [0,1] and then have the top half of the matrix sorted
        // downwards, and the bottom-half upwards. 
        DistMatrix<double,STAR,STAR> densityComp;
        Uniform( densityComp, nnu, nt, 0.5, 0.5 );
        {
            DistMatrix<double,STAR,STAR> densTop, densBot;
            PartitionDown( densityComp, densTop, densBot, nnu/2 ); 
            Sort( densTop, DESCENDING );
            Sort( densBot, ASCENDING );
        }

        // Sample from the complex unit ball
        DistMatrix<Complex<double>,STAR,STAR> sensitivity;
        Uniform( sensitivity, N0*N1, nc, Complex<double>(0.,0.), 1. );

        // Initialize each column to a 2*nnu length vector of samples from the
        // real ball of radius 0.5 centered at the origin
        DistMatrix<double,STAR,STAR> paths;
        Uniform( paths, 2*nnu, nt, 0., 0.5 );

        // Generate original data matrix
        DistMatrix<Complex<double>,STAR,VR> data;
        Uniform( data, nnu, nc*nt );

        if( print )
        {
            Print( densityComp, "density compensation" );
            Print( sensitivity, "coil sensitivities" );
            Print( paths, "paths" );
            Print( data, "data" );
        }
        if( display )
        {
            Display( densityComp, "density compensation" );
            Display( sensitivity, "coil sensitivities" );
            Display( paths, "paths" );
            Display( data, "data" );
        }

        // Initialize acquisition operator and its adjoint
        InitializeAcquisition
        ( densityComp, sensitivity, paths, nc, N0, N1, n0, n1, m );

        // Apply the adjoint of the acquisition operator
        DistMatrix<Complex<double>,VC,STAR> M;
        AdjointAcquisition( data, M );
        if( print )
            Print( M, "M := E' D" );
        if( display )
            Display( M, "M := E' D" );

        // Apply acquisition operator
        DistMatrix<Complex<double>,STAR,VR> R;
        Acquisition( M, R );
        if( print )
            Print( R, "R := E M = E E' D" );
        if( display )
            Display( R, "R := E M = E E' D" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
