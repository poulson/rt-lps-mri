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
using std::string;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

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

        if( formatInt < 1 || formatInt >= FileFormat_MAX )
            LogicError("Format integer must be in [1,",FileFormat_MAX,")");
        const auto format = static_cast<El::FileFormat>(formatInt);

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
            Write( densityComp, "density", format );
            Write( sensitivity, "sensitivity", format );
            Write( paths, "paths", format );
            Write( data, "data", format );
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
