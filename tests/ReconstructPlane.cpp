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
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

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
        const bool progress = Input("--progress","print L+S progress",true);
        const string sensName = 
            Input("--sens","sens. filename",string("sensitivity.bin"));
        const string densName = 
            Input("--dens","density filename",string("density.bin"));
        const string pathsName = 
            Input("--path","paths filename",string("paths.bin"));
        const string dataName = 
            Input("--data","data filename",string("data.bin"));
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",true);
#ifdef HAVE_QT5
        const int formatInt = Input("--format","format to store matrices",5);
#else
        const int formatInt = Input("--format","format to store matrices",1);
#endif
        ProcessInput();
        PrintInputReport();

        const elem::FileFormat format = 
            static_cast<elem::FileFormat>(formatInt);

        mpi::Barrier( comm );
        const double loadStart = mpi::Time();
        if( commRank == 0 )
        {
            std::cout << "Loading files...";
            std::cout.flush();
        }

        DistMatrix<double,STAR,STAR> densityComp;
        LoadDensity( nnu, nt, densName, densityComp );

        DistMatrix<Complex<double>,STAR,STAR> sensitivity;
        LoadSensitivity( N0, N1, nc, sensName, sensitivity );

        DistMatrix<double,STAR,STAR> paths;
        LoadPaths( nnu, nt, pathsName, paths );

        DistMatrix<Complex<double>,STAR,VR> data;
        LoadData( nnu, nc, nt, dataName, data );

        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "DONE. " << mpi::Time()-loadStart << " seconds" 
                      << std::endl;

        if( display )
        {
            Display( densityComp, "density compensation" );
            Display( sensitivity, "coil sensitivities" );
            Display( paths, "paths" );
            Display( data, "data" );
        }
        if( write )
        {
            Write( densityComp, format, "density" );
            Write( sensitivity, format, "sensitivity" );
            Write( paths, format, "paths" );
            Write( data, format, "data" );
        }

        // Initialize acquisition operator and its adjoint
        mpi::Barrier( comm );
        const double startInit = mpi::Time();
        if( commRank == 0 )
        {
            std::cout << "Initializing acquisition operator...";
            std::cout.flush();
        }
        InitializeAcquisition
        ( densityComp, sensitivity, paths, nc, N0, N1, n0, n1, m );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "DONE. " << mpi::Time()-startInit << " seconds"
                      << std::endl;

        mpi::Barrier( comm );
        const double startLPS = mpi::Time();
        if( commRank == 0 )
        {
            std::cout << "Starting L+S decomposition...";
            std::cout.flush();
        }
        DistMatrix<Complex<double>,VC,STAR> L, S;
        LPS
        ( data, L, S, tv, lambdaL, lambdaSRel, relTol, maxIts, 
          tryTSQR, progress );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "DONE. " << mpi::Time()-startLPS << " seconds"
                      << std::endl;

        if( write )
            WriteLPS( L, S, N0, N1, 0, tv, format );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
