// This is the main DLL file.

#include "stdafx.h"

#include "WinBinding.h"


using namespace ilPSP::ExternalBinding;


void Common_ilPSPInitialize() {
    Common_::ilPSPInitialize();
}


void ISparseSolver_Solve(int* SolverRef,int* N,double* x,double* rhs,int* ierr){
    ISparseSolver_::Solve( (*SolverRef), (*N), x, rhs, (*ierr));
}





