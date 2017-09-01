#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>


#include <mpi.h>

#include "monkey.h"

// size of the system to solve
#define NN 100

int main(int argc, char** argv) {
        int i, ierr, i0, iE, L, iRow, iCol;
        int N, n;
        int MtxRef;
        double d;
        int MPIrank, MPIsize, partRef;
        int RefSolver;
        double* X;
        double* rhs;
        double err;
        char solver[] = "<sparsesolver name=\"pressure-solver\"> \
                             <type>PCG</type> \
                             <library>monkey</library> \
                             <specific> \
                                  <MaxIterations>20000</MaxIterations> \
                                  <ConvergenceType>Relative</ConvergenceType> \
                                  <Tolerance>1.0e-9</Tolerance> \
                                  <DevType>CPU</DevType> \
                                  <MatrixType>ELLPACKcache</MatrixType> \
                             </specific> \
                        </sparsesolver>";


	ilPSPinit(); // If MPI is not initilized, ilPSP will do that;
                     // in this case, ilPSP will do also the MPI finalization.

	
        MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);
        MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);

        i0 = NN*MPIrank/MPIsize;
        iE = NN*(MPIrank + 1)/MPIsize;
        L=iE-i0;
        
        n = 1; N = NN;
        MsrMatrix_New(&MtxRef,&L,&N,&n,&n,&ierr); // local no. of rows, but global no. of columns!
        assert(ierr == 0);

        // Not necessary, just instructional ...
        ISparseMatrix_GetRowPart(&MtxRef,&partRef,&ierr);
        assert(ierr == 0);

        Partition_GetMPIrank(&partRef, &MPIrank, &ierr); assert(ierr == 0);
        Partition_GetMPIsize(&partRef, &MPIsize, &ierr); assert(ierr == 0);


        Common_ReleaseObject(&partRef,&ierr); // don't forget to release objects, as soon as you
                                              // don't need them anymore. Otherwise, the garbage 
                                              // collector wont be able to track them.
        // end of instructional.


        // set the simplest system, an identity matrix.
        d = 1.0;
        for( i = 0; i < L; i++) {
                iRow = i + i0;
                iCol = iRow;
                IMutableMatrix_SetEntry(&MtxRef,&iRow,&iCol,&d,&ierr);
                assert(ierr == 0);
        }

        // store matrix in a text file (e.g. for MATLAB import)
        IMutableMatrixEx_SaveToTextFileSparse(&MtxRef,"meinematrix.txt",&ierr); // error in MPI parallel mode, has to be investigated.
	assert(ierr == 0);

        // create solver:
        ISparseSolver_FromXML(&RefSolver,solver,&ierr);
        assert(ierr == 0);
    
        // pass matrix to solver:
        ISparseSolver_DefineMatrix(&RefSolver,&MtxRef, &ierr); 	assert(ierr == 0);
        Common_ReleaseObject(&MtxRef,&ierr); // MSR matrix is not required anymore, the solver 
                                             // creates his internal, more efficient format.
        assert(ierr == 0);

        // solve system
        X = malloc(sizeof(double)*L);    // alloc mem. for solution
        rhs = malloc(sizeof(double)*L);  // alloc mem. for right-hand-side

        for( i = 0; i < L; i++) {
                rhs[i] = 1;  // rhs: also simple
                X[i] = 0;    // initial guess: should be somehow sane (no NAN's, ...)
        }

        ISparseSolver_Solve(&RefSolver,&L,X,rhs,&ierr);
        assert(ierr == 0);

        // check: compute L1 error norm against exact solution.
        err = 0; 
        for( i = 0; i < L; i++) {
                d = X[i] - 1.0;
                if( d < 0) d*= -1.0;
                err += d;
        }
        printf("accumulated error: %f \n",err);

        // finalize: if ilPSP has done MPI_Initialize(...), it also does MPI_Finalize();
        //           otherwise, it leaves MPI as it is.
	ilPSPshutdown();
        return 0;
}
