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

namespace { 
bool mriInitializedElemental; 
int numMriInits = 0;
mri::Args* args = 0;
DEBUG_ONLY(std::stack<std::string> callStack)

bool initializedCoilPlans = false;
int numCoils;
int numTimesteps;
int numNonUniformPoints;
int firstBandwidth;
int secondBandwidth;
fftw_plan fftwForward, fftwBackward;
fftw_complex *g1, *g2;
El::DistMatrix<double,El::STAR,El::STAR>* coilPaths;
std::vector<nfft_plan> coilPlans;

bool initializedAcquisition = false;
El::DistMatrix<double,El::STAR,El::STAR>* densityComp;
El::DistMatrix<El::Complex<double>,El::STAR,El::STAR>* sensitivity;
El::DistMatrix<double,El::STAR,El::STAR>* sensitivityScalings;
}

namespace mri {

void PrintVersion( std::ostream& os )
{
    os << "RT-LPS-MRI version information:\n"
       << "  Git revision: " << EL_GIT_SHA1 << "\n"
       << "  Version:      " << RTLPSMRI_VERSION_MAJOR << "."
                             << RTLPSMRI_VERSION_MINOR << "\n"
       << "  Build type:   " << EL_CMAKE_BUILD_TYPE << "\n"
       << std::endl;
}

void PrintConfig( std::ostream& os )
{
    os << "RT-LPS-MRI configuration:\n";
    os << "  NFFT_INC_DIR: " << NFFT_INC_DIR << "\n";
    El::PrintConfig( os );
}

void PrintCCompilerInfo( std::ostream& os )
{
    os << "RT-LPS-MRI's C compiler info:\n"
       << "  CMAKE_C_COMPILER:    " << EL_CMAKE_C_COMPILER << "\n"
       << "  MPI_C_COMPILER:      " << EL_MPI_C_COMPILER << "\n"
       << "  MPI_C_INCLUDE_PATH:  " << EL_MPI_C_INCLUDE_PATH << "\n"
       << "  MPI_C_COMPILE_FLAGS: " << EL_MPI_C_COMPILE_FLAGS << "\n"
       << "  MPI_C_LINK_FLAGS:    " << EL_MPI_C_LINK_FLAGS << "\n"
       << "  MPI_C_LIBRARIES:     " << EL_MPI_C_LIBRARIES << "\n"
       << std::endl;
}

void PrintCxxCompilerInfo( std::ostream& os )
{
    os << "RT-LPS-MRI's C++ compiler info:\n"
       << "  CMAKE_CXX_COMPILER:    " << EL_CMAKE_CXX_COMPILER << "\n"
       << "  CXX_FLAGS:             " << EL_CXX_FLAGS << "\n"
       << "  MPI_CXX_COMPILER:      " << EL_MPI_CXX_COMPILER << "\n"
       << "  MPI_CXX_INCLUDE_PATH:  " << EL_MPI_CXX_INCLUDE_PATH << "\n"
       << "  MPI_CXX_COMPILE_FLAGS: " << EL_MPI_CXX_COMPILE_FLAGS << "\n"
       << "  MPI_CXX_LINK_FLAGS:    " << EL_MPI_CXX_LINK_FLAGS << "\n"
       << "  MPI_CXX_LIBRARIES:     " << EL_MPI_CXX_LIBRARIES << "\n"
       << std::endl;
}

bool Initialized()
{ return ::numMriInits > 0; }

void Initialize( int& argc, char**& argv )
{
    // If RT-LPS-MRI has already been initialized, this is a no-op
    if( ::numMriInits > 0 )
    {
        ++::numMriInits;
        return;
    }

    ::args = new Args( argc, argv );
    ::numMriInits = 1;
    if( !El::Initialized() )
    {
        El::Initialize( argc, argv );
        ::mriInitializedElemental = true;
    }
    else
    {
        ::mriInitializedElemental = false;
    }
    El::SetColorMap( El::GRAYSCALE );
}

void Finalize()
{
    if( ::numMriInits <= 0 )
        throw std::logic_error("Finalized RT-LPS-MRI more than initialized");
    --::numMriInits;
    if( ::mriInitializedElemental )
        El::Finalize();

    if( ::numMriInits == 0 )
    {
        if( InitializedAcquisition() )
            FinalizeAcquisition();
        else if( InitializedCoilPlans() )
            FinalizeCoilPlans();
        delete ::args;    
        ::args = 0;
    }
}

DEBUG_ONLY(
    void PushCallStack( std::string s )
    { ::callStack.push( s ); }

    void PopCallStack()
    { ::callStack.pop(); }

    void DumpCallStack()
    {
        std::ostringstream msg;
        while( !::callStack.empty() )
        {
            msg << "[" << ::callStack.size() << "]: " << ::callStack.top() 
                << "\n";
            ::callStack.pop();
        }
        std::cerr << msg.str();;
        std::cerr.flush();
    }
)

void ReportException( std::exception& e )
{
    El::ReportException( e );
    DEBUG_ONLY(DumpCallStack())
}

Args& GetArgs()
{
    if( args == 0 )
        throw std::runtime_error("No available instance of Args");
    return *::args;
}

bool InitializedCoilPlans()
{ return ::initializedCoilPlans; }

bool InitializedAcquisition()
{ return ::initializedAcquisition; }

// Each column of X corresponds to the Fourier-domain path for each timestep.
// The trajectories are the same for each coil.
void InitializeCoilPlans
( const DistMatrix<double,STAR,STAR>& X, 
  int numCoils, int N0, int N1, int n0, int n1, int m )
{
    DEBUG_ONLY(
        CallStackEntry cse("InitializeCoilPlans");
        if( InitializedCoilPlans() )
            LogicError("Already initialized coil plans");
    )
    const int dim = 2;
    const int numNonUniform = X.Height()/dim;
    const int numTimesteps = X.Width();

    ::numCoils = numCoils;
    ::numTimesteps = numTimesteps;
    ::numNonUniformPoints = numNonUniform;
    ::firstBandwidth = N0;
    ::secondBandwidth = N1;

    ::coilPaths = new DistMatrix<double,STAR,STAR>( X );

    int NN[dim] = { N0, N1 };
    int nn[dim] = { n0, n1 };
    const int nTotal = n0*n1;

    unsigned nfftFlags = PRE_PHI_HUT| PRE_FULL_PSI| NFFT_SORT_NODES;
    unsigned fftwFlags = FFTW_MEASURE| FFTW_DESTROY_INPUT;

    ::g1 = (fftw_complex*)nfft_malloc( nTotal*sizeof(fftw_complex) );
    ::g2 = (fftw_complex*)nfft_malloc( nTotal*sizeof(fftw_complex) );
    ::fftwForward = 
        fftw_plan_dft_2d( n0, n1, ::g1, ::g2, FFTW_FORWARD, fftwFlags );
    ::fftwBackward = 
        fftw_plan_dft_2d( n0, n1, ::g2, ::g1, FFTW_BACKWARD, fftwFlags );

    ::coilPlans.clear();
    ::coilPlans.resize( numTimesteps );
    for( int t=0; t<numTimesteps; ++t )
    {
        nfft_plan& plan = ::coilPlans[t];
        plan.x = ::coilPaths->Buffer(0,t);
        plan.g1 = ::g1;
        plan.g2 = ::g2;
        plan.my_fftw_plan1 = ::fftwForward;
        plan.my_fftw_plan2 = ::fftwBackward;
        nfft_init_guru
        ( &plan, dim, NN, numNonUniform, nn, m, nfftFlags, fftwFlags );
        if( plan.nfft_flags & PRE_ONE_PSI )
            nfft_precompute_one_psi( &plan );
    }

    ::initializedCoilPlans = true;
}

void InitializeAcquisition
( const DistMatrix<double,         STAR,STAR>& dens, 
  const DistMatrix<Complex<double>,STAR,STAR>& sens,
  const DistMatrix<double,         STAR,STAR>& X, 
  int numCoils, int N0, int N1, int n0, int n1, int m )
{
    DEBUG_ONLY(
        CallStackEntry cse("InitializeAcquisition");
        if( InitializedAcquisition() )
            LogicError("Already initialized acquisition operator");
        if( sens.Height() != N0*N1 || sens.Width() != numCoils )
            LogicError("Coil sensitivity matrix of the wrong size");
        if( dens.Height() != X.Height()/2 || dens.Width() != X.Width() )
            LogicError("Density composition matrix of the wrong size");
    )
    InitializeCoilPlans( X, numCoils, N0, N1, n0, n1, m );

    ::densityComp = new DistMatrix<double,STAR,STAR>( dens );

    // We have to rearrange sensitivity to be row-major instead of column-major
    // in order to be compatible with NFFT3's row-major image ordering
    ::sensitivity = new DistMatrix<Complex<double>,STAR,STAR>( sens.Grid() );
    Zeros( *::sensitivity, N0*N1, numCoils );
    for( int c=0; c<numCoils; ++c )
    {
        auto newCol = ::sensitivity->Buffer(0,c);
        const auto oldCol = sens.LockedBuffer(0,c);
        for( int j0=0; j0<N0; ++j0 )
            for( int j1=0; j1<N1; ++j1 )
                newCol[j1+j0*N1] = oldCol[j0+j1*N0];
    } 

    ::sensitivityScalings = new DistMatrix<double,STAR,STAR>( sens.Grid() );
    Zeros( *::sensitivityScalings, N0*N1, 1 );
    for( int i=0; i<N0*N1; ++i )
    {
        for( int j=0; j<numCoils; ++j )
        {
            const Complex<double> eta = ::sensitivity->GetLocal(i,j);
            ::sensitivityScalings->UpdateLocal( i, 0, RealPart(eta*Conj(eta)) );
        }
    }

    ::initializedAcquisition = true;
}

void FinalizeCoilPlans()
{
    DEBUG_ONLY(
        CallStackEntry cse("FinalizeCoilPlans");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    const int numTimesteps = NumTimesteps();
    for( int t=0; t<numTimesteps; ++t )
        nfft_finalize( &::coilPlans[t] );
    ::coilPlans.clear();
    delete ::coilPaths;

    fftw_destroy_plan( ::fftwBackward );
    fftw_destroy_plan( ::fftwForward );
    nfft_free( ::g2 );
    nfft_free( ::g1 );

    ::initializedCoilPlans = false;
}

void FinalizeAcquisition()
{
    DEBUG_ONLY(
        CallStackEntry cse("FinalizeAcquisition");
        if( !InitializedAcquisition() )
            LogicError("Have not yet initialized acquisition operator");
    )
    ::initializedAcquisition = false;
    delete ::densityComp;
    delete ::sensitivity;
    delete ::sensitivityScalings;
    FinalizeCoilPlans();
}

int NumCoils()
{
    DEBUG_ONLY(
        CallStackEntry cse("NumCoils");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return ::numCoils;
}

int NumTimesteps()
{
    DEBUG_ONLY(
        CallStackEntry cse("NumTimesteps");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return ::numTimesteps;
}

int NumNonUniformPoints()
{ 
    DEBUG_ONLY(
        CallStackEntry cse("NumNonUniformPoints");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return ::numNonUniformPoints; 
}

int FirstBandwidth()
{ 
    DEBUG_ONLY(
        CallStackEntry cse("FirstBandwidth");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return ::firstBandwidth; 
}

int SecondBandwidth()
{ 
    DEBUG_ONLY(
        CallStackEntry cse("SecondBandwidth");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return ::secondBandwidth; 
}

nfft_plan& CoilPlan( int path )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("CoilPlan");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return ::coilPlans[path]; 
}

const DistMatrix<double,STAR,STAR>& CoilPaths()
{
    DEBUG_ONLY(
        CallStackEntry cse("CoilPaths");
        if( !InitializedCoilPlans() )
            LogicError("Have not yet initialized coil plans");
    )
    return *::coilPaths;
}

const DistMatrix<double,STAR,STAR>& DensityComp()
{
    DEBUG_ONLY(
        CallStackEntry cse("DensityComp");
        if( !InitializedAcquisition() )
            LogicError("Have not yet initialized acquisition operator");
    )
    return *::densityComp;
}

const DistMatrix<El::Complex<double>,STAR,STAR>& Sensitivity()
{
    DEBUG_ONLY(
        CallStackEntry cse("Sensitivity");
        if( !InitializedAcquisition() )
            LogicError("Have not yet initialized acquisition operator");
    )
    return *::sensitivity;
}

const DistMatrix<double,STAR,STAR>& SensitivityScalings()
{
    DEBUG_ONLY(
        CallStackEntry cse("SensitivityScalings");
        if( !InitializedAcquisition() )
            LogicError("Have not yet initialized acquisition operator");
    )
    return *::sensitivityScalings;
}

} // namespace mri
