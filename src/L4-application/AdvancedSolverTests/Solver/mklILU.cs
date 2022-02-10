using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AdvancedSolverTests.Solver {
    public class mklILU {
        private static (MultidimensionalArray L_matlab, MultidimensionalArray U_matlab, MultidimensionalArray A) GenMatrix() {

            int basis = 5;
            int size = (int)Math.Floor(Math.Pow(basis, 2));
            var L_matlab = MultidimensionalArray.Create(size, size);
            var U_matlab = MultidimensionalArray.Create(size, size);
            var A = MultidimensionalArray.Create(size, size);

            using (BatchmodeConnector matlab = new BatchmodeConnector()) {
                //note: BatchmodeCon maybe working on proc0 but savetotxt file, etc. (I/O) is full mpi parallel
                //so concider this as full mpi-parallel

                matlab.Cmd($"A = gallery('neumann', {size}) + speye({size});");
                matlab.Cmd($"options.type = 'nofill';");
                matlab.Cmd($"[L, U] = ilu(A, options);");
                matlab.Cmd($"Lfull = full(L);");
                matlab.Cmd($"Ufull = full(U);");
                matlab.Cmd($"Afull = full(A);");
                matlab.GetMatrix(A, "Afull");
                matlab.GetMatrix(L_matlab, "Lfull");
                matlab.GetMatrix(U_matlab, "Ufull");
                matlab.Execute();
            }
            Console.WriteLine(L_matlab.InfNorm());
            return (L_matlab, U_matlab, A);
        }

        private static MsrMatrix ExecuteILU(MultidimensionalArray A) {
            var ilu = new ilPSP.LinSolvers.ILU.ILUSolver();
            var A_BMsr = new MsrMatrix(A.Lengths[0],A.Lengths[1]);
            A_BMsr.AccDenseMatrix(1.0, A);
            ilu.DefineMatrix(A_BMsr);
            var LU=ilu.GetILUFactorization;
            return LU;
        }

        public static void CompareFactorization() {
            var matrices=GenMatrix();
            var ILU = ExecuteILU(matrices.A);
            Console.WriteLine($"dim mkl ILU: {ILU.RowPartitioning.LocalLength}, dim matlab L: {matrices.L_matlab.Lengths[0]}, dim matlab U: {matrices.U_matlab.Lengths[0]}");
            //matrices.L_matlab.SaveToTextFile("Lmatlab");
            //matrices.U_matlab.SaveToTextFile("Umatlab");
            //ILU.SaveToTextFileSparse("ILU");
            ILU.AccBlock(0,0,-1.0, matrices.L_matlab);
            ILU.AccBlock(0,0,-1.0, matrices.U_matlab);
            Console.WriteLine(ILU.InfNorm()-1<1E-14);
        }
    }
}
