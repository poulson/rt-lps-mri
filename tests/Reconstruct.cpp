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
    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );

    try
    {
        const int np = Input("--np","number of planes to process",38);
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
        const string dataBase = 
            Input("--data","data base filename",string("data"));
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",true);
#ifdef HAVE_QT5
        const int formatInt = Input("--format","format to store matrices",7);
#else
        const int formatInt = Input("--format","format to store matrices",1);
#endif
        ProcessInput();
        PrintInputReport();

        if( formatInt < 1 || formatInt >= FileFormat_MAX )
            LogicError("Format integer must be in [1,",FileFormat_MAX,")");
        const auto format = static_cast<El::FileFormat>(formatInt);

        // Load and possibly display and write the plane-independent data
        mpi::Barrier( comm );
        if( commRank == 0 )
        {
            std::cout << "Loading plane-independent data...";     
            std::cout.flush();
        }
        const double startLoad = mpi::Time();
        DistMatrix<double,STAR,STAR> paths, densityComp;
        DistMatrix<Complex<double>,STAR,STAR> sensitivity;
        LoadPaths( nnu, nt, pathsName, paths );
        LoadDensity( nnu, nt, densName, densityComp );
        LoadSensitivity( N0, N1, nc, sensName, sensitivity );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "DONE. " << mpi::Time()-startLoad << " seconds" 
                      << std::endl;
        if( display )
        {
            Display( densityComp, "density compensation" );
            Display( sensitivity, "coil sensitivities" );
            Display( paths, "paths" );
        }
        if( write )
        {
            Write( densityComp, "density", format );
            Write( sensitivity, "sensitivity", format );
            Write( paths, "paths", format );
        }

        // Initialize the acquisition operator and its adjoint
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

        // Exploit the available trivial plane parallelism
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "Starting trivially parallel work..." << std::endl;
        const double trivialStart = mpi::Time();
        const int numSeqRounds = np / commSize;
        Grid selfGrid( mpi::COMM_SELF );
        DistMatrix<Complex<double>,STAR,VR> data(selfGrid);
        DistMatrix<Complex<double>,VC,STAR> L(selfGrid), S(selfGrid);
        for( Int round=0; round<numSeqRounds; ++round )
        {
            const int plane = commRank + round*commSize;
            std::ostringstream os, binOs;
            os << dataBase << "-" << plane;
            LoadData( nnu, nc, nt, os.str()+".bin", data );
            if( display )
                Display( data, os.str() );
            if( write )
                Write( data, os.str(), format );

            LPS
            ( data, L, S, tv, lambdaL, lambdaSRel, relTol, maxIts, tryTSQR, 
              progress );
            if( write )
                WriteLPS( L, S, N0, N1, plane, tv, format );
        }
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "Finished trivially parallel section: " 
                      << mpi::Time()-trivialStart << " seconds" << std::endl;

        // Process the remaining planes using subteams
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "Starting parallel work..." << std::endl;
        const double parallelStart = mpi::Time(); 
        const int numSeqPlanes = numSeqRounds*commSize; 
        const int numParPlanes = np - numSeqPlanes;
        if( numParPlanes > 0 )
        {
            // Split the communicator
            int color, key;
            const int mainTeamSize = commSize / numParPlanes;
            if( commRank < mainTeamSize*(numParPlanes-1) )
            {
                color = commRank / mainTeamSize;
                key = commRank % mainTeamSize;
            } 
            else
            {
                color = numParPlanes-1;
                key = commRank - mainTeamSize*(numParPlanes-1);
            } 
            mpi::Comm subComm;
            mpi::Split( comm, color, key, subComm );
            Grid subGrid( subComm );
            data.SetGrid( subGrid );
            L.SetGrid( subGrid );
            S.SetGrid( subGrid );

            const int plane = numSeqPlanes + color;
            std::ostringstream os, binOs;
            os << dataBase << "-" << plane;
            LoadData( nnu, nc, nt, os.str()+".bin", data );
            if( display )
                Display( data, os.str() );
            if( write )
                Write( data, os.str(), format );

            mpi::Barrier( comm );
            const double lpsStart = mpi::Time();
            LPS
            ( data, L, S, tv, lambdaL, lambdaSRel, relTol, maxIts, tryTSQR, 
              progress );
            mpi::Barrier( comm );
            if( commRank == 0 )
                std::cout << "  Parallel LPS's took " << mpi::Time()-lpsStart 
                          << " seconds" << std::endl;
            if( write )
                WriteLPS( L, S, N0, N1, plane, tv, format );
        }
        if( commRank == 0 )
            std::cout << "Finished parallel section: " 
                      << mpi::Time()-parallelStart << " seconds" << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
