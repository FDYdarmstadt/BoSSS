using System;
using ilPSP.ExternalBinding;
using MPI.Wrappers;
using System.Runtime.InteropServices;

namespace CsExtBinding {
    class MainClass {
        const int NN = 100;

        unsafe public static void Main(string[] args) {

            int i, ierr, i0, iE, L, iRow, iCol;
            int N, n;
            int MtxRef;
            double d;
            int MPIrank, MPIsize;
            int RefSolver;
            double[] X;
            double[] rhs;
            double err;
            string solver1 = "<sparsesolver name=\"pressure-solver\"> "
                             + "<type>PCG</type> "
                             + "<library>monkey</library> "
                             + "<specific> "
                             + "     <MaxIterations>20000</MaxIterations> "
                             + "     <ConvergenceType>Relative</ConvergenceType> "
                             + "     <Tolerance>1.0e-9</Tolerance> "
                             + "     <DevType>CPU</DevType> "
                             + "     <MatrixType>ELLPACKcache</MatrixType> "
                             + "</specific> "
                        + "</sparsesolver>";
            string solver2 = "<sparsesolver name=\"pressure-solver\"> "
                             + "<type>PCG</type> "
                             + "<library>hypre</library> "
                             + "<specific> "
                             + "     <MaxIterations>20000</MaxIterations> "
                             + "     <ConvergenceType>Relative</ConvergenceType> "
                             + "     <Tolerance>1.0e-9</Tolerance> "
                             + "</specific> "
                        + "</sparsesolver>";
			

			Common_.ilPSPInitialize(); // If MPI is not initilized, ilPSP will do that;
                                       // in this case, ilPSP will do also the MPI finalization.


            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MPIsize);

            i0 = NN * MPIrank / MPIsize;
            iE = NN * (MPIrank + 1) / MPIsize;
            L = iE - i0;

            n = 1; N = NN;
            MsrMatrix_.New(out MtxRef, ref L, ref N, ref n, ref n, out ierr); // local no. of rows, but global no. of columns!
            //assert(ierr == 0);

            //        // Not necessary, just instructional ...
            //	ISparseMatrix_GetRowPart(&MtxRef,&partRef,&ierr);
            //	assert(ierr == 0);
            //
            //        Partition_GetMPIrank(&partRef, &MPIrank, &ierr); assert(ierr == 0);
            //        Partition_GetMPIsize(&partRef, &MPIsize, &ierr); assert(ierr == 0);
            //
            //
            //        Common_ReleaseObject(&partRef,&ierr); // don't forget to release objects, as soon as you
            //                                              // don't need them anymore. Otherwise, the garbage 
            //                                              // collector wont be able to track them.
            //        // end of instructional.


            // set the simplest system, an identity matrix.
            d = 1.0;
            for (i = 0; i < L; i++) {
                iRow = i + i0;
                iCol = iRow;
                IMutableMatrix_.SetEntry(ref MtxRef, ref iRow, ref iCol, ref d, out ierr);
                //assert(ierr == 0);
            }

            // store matrix in a text file (e.g. for MATLAB import)
            byte* path = (byte*)Marshal.StringToCoTaskMemAnsi("meinematrix.txt");
            IMutableMatrixEx_.SaveToTextFileSparse(ref MtxRef, path, out ierr); // error in MPI parallel mode, has to be investigated.

            //assert(ierr == 0);

            // create solver:
            unsafe {
                byte* _solver = (byte*)Marshal.StringToCoTaskMemAnsi(solver2);
                ISparseSolver_.FromXML(out RefSolver, _solver, out ierr);
                Marshal.FreeCoTaskMem((IntPtr)_solver);

            }
            // pass matrix to solver:
            ISparseSolver_.DefineMatrix(ref RefSolver, ref MtxRef, out ierr);
            Common_.ReleaseObject(ref MtxRef, out ierr); // MSR matrix is not required anymore, the solver 

            // solve system
            X = new double[L];  // alloc mem. for solution
            rhs = new double[L];  // alloc mem. for right-hand-side

            for (i = 0; i < L; i++) {
                rhs[i] = 1;  // rhs: also simple
                X[i] = 0;    // initial guess: should be somehow sane (no NAN's, ...)
            }

            unsafe {
                fixed (double* pX = X, pRhs = rhs) {
                    ISparseSolver_.Solve(ref RefSolver, ref L, pX, pRhs, out ierr);
                }
            }

            // check: compute L1 error norm against exact solution.
            err = 0;
            for (i = 0; i < L; i++) {
                d = X[i] - 1.0;
                if (d < 0) d *= -1.0;
                err += d;
            }
            Console.WriteLine("accumulated error: "+err);

            // finalize: if ilPSP has done MPI_Initialize(...), it also does MPI_Finalize();
            //           otherwise, it leaves MPI as it is.
            Common_.ilPSPFinalize();


        }
    }
}

