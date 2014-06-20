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
namespace blas = El::blas;
namespace lapack = El::lapack;
namespace mpi = El::mpi;

// Pull in a number of useful enums from Elemental
using namespace El::DistNS;
using namespace El::LeftOrRightNS;
using namespace El::OrientationNS;
using namespace El::UnitOrNonUnitNS;
using namespace El::UpperOrLowerNS;
using namespace El::SortTypeNS;
using namespace El::FileFormatNS;

// For scalar operations
using El::Int;
using El::Base;
using El::Complex;
using El::Abs;
using El::Conj;
using El::RealPart;
using El::Sqrt;
using El::SampleUniform;

// Pull in a few classes from Elemental
using El::Matrix;
using El::Grid;
using El::DistMatrix;
using El::View;
using El::LockedView;

// Pull in a few indexing routines
using El::Shift;
using El::Length;

// Matrix generation, sorting, and norm calculation
using El::Uniform;
using El::Sort;
using El::FrobeniusNorm;
using El::MaxNorm;

// Utilities
using El::Read;
using El::Write;
using El::FileSize;
using El::LogicError;
using El::RuntimeError;

void PrintVersion( std::ostream& os=std::cout );
void PrintConfig( std::ostream& os=std::cout );
void PrintCCompilerInfo( std::ostream& os=std::cout );
void PrintCxxCompilerInfo( std::ostream& os=std::cout );
 
bool Initialized();
void Initialize( int& argc, char**& argv );
void Finalize();

// For getting the MPI argument instance (for internal usage)
class Args : public El::choice::MpiArgs
{
public:
    Args
    ( int argc, char** argv,
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr )
    : El::choice::MpiArgs(argc,argv,comm,error)
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

DEBUG_ONLY(
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
)

void ReportException( std::exception& e );

bool InitializedCoilPlans();
bool InitializedAcquisition();
void InitializeCoilPlans
( const DistMatrix<double,STAR,STAR>& paths, 
  int numCoils, int N0, int N1, int n0, int n1, int m );
void InitializeAcquisition
( const DistMatrix<double,         STAR,STAR>& densityComp,
  const DistMatrix<Complex<double>,STAR,STAR>& sensitivity,
  const DistMatrix<double,         STAR,STAR>& paths, 
  int numCoils, int N0, int N1, int n0, int n1, int m );
void FinalizeCoilPlans();
void FinalizeAcquisition();

int NumCoils();
int NumTimesteps();
int NumNonUniformPoints();
int FirstBandwidth();
int SecondBandwidth();

nfft_plan& CoilPlan( int path );

// 2*M x numTimesteps
const DistMatrix<double,STAR,STAR>& CoilPaths();

// M x numTimesteps
const DistMatrix<double,STAR,STAR>& DensityComp();

// N0*N1 x numCoils
const DistMatrix<Complex<double>,STAR,STAR>& Sensitivity();

// N0*N1 x 1
const DistMatrix<double,STAR,STAR>& SensitivityScalings();

} // namespace mri

#endif // ifndef RTLPSMRI_CORE_ENVIRONMENT_DECL_HPP
