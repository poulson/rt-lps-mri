/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.

   This file is part of RealTime-MRI and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "rt-mri.hpp"
using namespace rtmri;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::CommRank( comm );

    if( commRank == 0 )
    {
        PrintVersion();
        PrintCCompilerInfo();
        PrintCxxCompilerInfo();
    }

    Finalize();
    return 0;
}
