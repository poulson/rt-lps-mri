/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the
   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#include "rt-lps-mri.hpp"
using namespace mri;
using std::string;

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
        const bool tv = Input("--tv","TV clipping for sparsity",true);
        const double lambdaL = Input("--lambdaL","low-rank scale",0.025);
        const double lambdaSRel = Input("--lambdaSRel","sparse rel scale",0.5);
        const double relTol = Input("--relTol","relative L+S tolerance",0.0025);
        const int maxIts = Input("--maxIts","max L+S iterations",100);
        const bool tryTSQR = Input("--tryTSQR","try Tall-Skinny QR?",false);
        const string sensName = 
            Input("--sens","sens. filename",string("sensitivity.bin"));
        const string densName = 
            Input("--dens","density filename",string("density.bin"));
        const string pathsName = 
            Input("--path","paths filename",string("paths.bin"));
        const string dataName = 
            Input("--data","data filename",string("data.bin"));
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double,STAR,STAR> densityComp;
        LoadDensity( nnu, nt, densName, densityComp );

        DistMatrix<Complex<double>,STAR,STAR> sensitivity;
        LoadSensitivity( N0, N1, nc, sensName, sensitivity );

        DistMatrix<double,STAR,STAR> paths;
        LoadPaths( nnu, nt, pathsName, paths );

        DistMatrix<Complex<double>,STAR,VR> data;
        LoadData( nnu, nc, nt, dataName, data );

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

        DistMatrix<Complex<double>,VC,STAR> L, S;
        LPS( data, L, S, tv, lambdaL, lambdaSRel, relTol, maxIts, tryTSQR );
        if( print )
        {
            Print( L, "L" );
            Print( S, "S" );
        }
        if( display )
        {
            Display( L, "L" );
            Display( S, "S" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
