/*
   Copyright (c) 2013, Jack Poulson, Ricardo Otazo, and Emmanuel Candes
   All rights reserved.
 
   This file is part of Real-Time Low-rank Plus Sparse MRI (RT-LPS-MRI) and is 
   under the BSD 2-Clause License, which can be found in the LICENSE file in the   root directory, or at http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef RTLPSMRI_CORE_ENVIRONMENT_DECL_HPP
#define RTLPSMRI_CORE_ENVIRONMENT_DECL_HPP

namespace mri {

typedef unsigned char byte;

// Pull in some of Elemental's imported libraries
namespace blas = elem::blas;
namespace lapack = elem::lapack;
namespace mpi = elem::mpi;

// Pull in a number of useful enums from Elemental
using namespace elem::distribution_wrapper;
using namespace elem::left_or_right_wrapper;
using namespace elem::orientation_wrapper;
using namespace elem::unit_or_non_unit_wrapper;
using namespace elem::upper_or_lower_wrapper;
using namespace elem::sort_type_wrapper;

// For scalar operations
using elem::Int;
using elem::Base;
using elem::Complex;
using elem::Abs;
using elem::Conj;
using elem::RealPart;
using elem::Sqrt;
using elem::SampleUniform;

// Pull in a few classes from Elemental
using elem::Matrix;
using elem::Grid;
using elem::DistMatrix;
using elem::View;
using elem::LockedView;

// Pull in a few indexing routines
using elem::Shift;
using elem::Length;

// Matrix generation, sorting, and norm calculation
using elem::Uniform;
using elem::Sort;
using elem::FrobeniusNorm;

// Exception generation
using elem::LogicError;
using elem::RuntimeError;

void PrintVersion( std::ostream& os=std::cout );
void PrintConfig( std::ostream& os=std::cout );
void PrintCCompilerInfo( std::ostream& os=std::cout );
void PrintCxxCompilerInfo( std::ostream& os=std::cout );
 
bool Initialized();
void Initialize( int& argc, char**& argv );
void Finalize();

// For getting the MPI argument instance (for internal usage)
class Args : public elem::choice::MpiArgs
{
public:
    Args
    ( int argc, char** argv,
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr )
    : elem::choice::MpiArgs(argc,argv,comm,error)
    { }
    virtual ~Args() { }
protected:
    virtual void HandleVersion( std::ostream& os=std::cout ) const;
    virtual void HandleBuild( std::ostream& os=std::cout ) const;
};
Args& GetArgs();

// For processing command-line arguments
template<typename T>
T Input( std::string name, std::string desc );
template<typename T>
T Input( std::string name, std::string desc, T defaultVal );
void ProcessInput();
void PrintInputReport();

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();

class CallStackEntry
{
public:
    CallStackEntry( std::string s ) 
    { 
        if( !std::uncaught_exception() )
            PushCallStack(s);
    }
    ~CallStackEntry() 
    { 
        if( !std::uncaught_exception() )
            PopCallStack(); 
    }
};
#endif

void ReportException( std::exception& e );

bool InitializedCoilPlans();
bool InitializedAcquisition();
void InitializeCoilPlans
( const DistMatrix<double,STAR,VR>& X, int numCoils, int numTimesteps, 
  int N0, int N1, int n0, int n1, int m );
void InitializeAcquisition
( const DistMatrix<double,         STAR,STAR>& densityComp,
  const DistMatrix<Complex<double>,STAR,STAR>& sensitivity,
  const DistMatrix<double,STAR,VR>& X, int numCoils, int numTimesteps, 
  int N0, int N1, int n0, int n1, int m );
void FinalizeCoilPlans();
void FinalizeAcquisition();

int NumCoils();
int NumTimesteps();
int NumLocalPaths();
int NumNonUniformPoints();
int FirstBandwidth();
int SecondBandwidth();

nfft_plan& LocalCoilPlan( int localPath );

// M x numCoils*numTimesteps
const DistMatrix<double,STAR,VR>& CoilPaths();

// M x numTimesteps
const DistMatrix<double,STAR,STAR>& DensityComp();

// N0*N1 x numCoils
const DistMatrix<Complex<double>,STAR,STAR>& Sensitivity();

// N0*N1 x 1
const DistMatrix<double,STAR,STAR>& SensitivityScalings();

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_ENVIRONMENT_DECL_HPP
