/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.LinSolvers;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using BoSSS.Foundation;
using ilPSP.Connectors.Matlab;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BoSSS.Solution.Multigrid {
    public class SchurPrecond : ISolverSmootherTemplate, ISolverWithCallback {
        public int IterationsInNested {
            get {
                throw new NotImplementedException();
            }
        }

        public int ThisLevelIterations {
            get {
                throw new NotImplementedException();
            }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get {
                throw new NotImplementedException();
            }
            set {
                throw new NotImplementedException();
            }
        }

        MultigridOperator m_mgop;

        BlockMsrMatrix Mtx;

        MsrMatrix P;
        MsrMatrix ConvDiff, pGrad, divVel, ConvDiffpGrad, SchurMtx;
        int[] Uidx, Pidx;

        public void Init(MultigridOperator op) {
            int D = op.Mapping.GridData.SpatialDimension;
            var M = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_mgop = op;

            if (!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            Uidx = MgMap.ProblemMapping.GetSubvectorIndices(true, D.ForLoop(i => i));
            Pidx = MgMap.ProblemMapping.GetSubvectorIndices(true, D);

            int Upart = Uidx.Length;
            int Ppart = Pidx.Length;

            ConvDiff = new MsrMatrix(Upart, Upart, 1, 1);
            pGrad = new MsrMatrix(Upart, Ppart, 1, 1);
            ConvDiffpGrad = new MsrMatrix(Upart, Upart + Ppart);
            divVel = new MsrMatrix(Ppart, Upart, 1, 1);

            M.AccSubMatrixTo(1.0, ConvDiff, Uidx, default(int[]), Uidx, default(int[]));
            M.AccSubMatrixTo(1.0, pGrad, Uidx, default(int[]), Pidx, default(int[]));
            M.AccSubMatrixTo(1.0, divVel, Pidx, default(int[]), Uidx, default(int[]));

            Mtx = M;

            int L = M.RowPartitioning.LocalLength;

            int i0 = Mtx.RowPartitioning.i0;

            P = new MsrMatrix(Mtx);
            P.Clear();

            // Debugging output
            //ConvDiff.SaveToTextFileSparse("ConvDiff");
            //divVel.SaveToTextFileSparse("divVel");
            //pGrad.SaveToTextFileSparse("pGrad");

            // Building the Schur complement
            //MultidimensionalArray ConvDiffInv = MultidimensionalArray.Create(Uidx.Length, Uidx.Length);
            MultidimensionalArray Schur = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);

            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(ConvDiff, "ConvDiff");
                bmc.PutSparseMatrix(divVel, "divVel");
                bmc.PutSparseMatrix(pGrad, "pGrad");
                bmc.Cmd("ConvDiffInv = inv(full(ConvDiff))");
                bmc.Cmd("Schur = divVel*ConvDiffInv");
                bmc.Cmd("Schur = Schur*pGrad");
                //bmc.Cmd("SchurInv = inv(Schur)");
                //bmc.GetMatrix(ConvDiffInv, "ConvDiffInv");
                bmc.GetMatrix(Schur, "-Schur");
                bmc.Execute(false);
            }



            SchurMtx = Schur.ToMsrMatrix();
            //var ConvDiffInvMtx = ConvDiffInv.ToMsrMatrix();

            //ConvDiff.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Uidx);
            //pGrad.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Pidx);
            //SchurMtx.AccSubMatrixTo(1.0, P, default(int[]), Pidx, default(int[]), Pidx);
            //// x= inv(P)*b !!!!! To be done with approximate Inverse
            // P.SpMV(1, B, 0, X);
        }

        public void ResetStat() {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {

            SolveSubproblems(X, B);

            //Directly invert Preconditioning Matrix
            //using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
            //    solver.DefineMatrix(P);
            //    solver.Solve(X, B);
            //}

        }


        /// <summary>
        /// Solve Preconditioning Matrix in Subsystems with ConvDiff, pGrad and Schur
        /// </summary>
        /// <typeparam name="U"></typeparam>
        /// <typeparam name="V"></typeparam>
        /// <param name="X"></param>
        /// <param name="B"></param>
        public void SolveSubproblems<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {

            var Bu = new double[Uidx.Length];
            var Xu = Bu.CloneAs();
            Bu = B.GetSubVector(Uidx,default(int[]));
            var Bp = new double[Pidx.Length];
            var Xp = Bp.CloneAs();
            Bp = B.GetSubVector(Pidx, default(int[]));

            Xu = X.GetSubVector(Uidx, default(int[]));
            Xp = X.GetSubVector(Pidx, default(int[]));

            // Solve Schur*s=q
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
                solver.DefineMatrix(SchurMtx);
                solver.Solve(Xp, Bp);
            }

            // Solve ConvDiff*w=v-q*pGrad
            pGrad.SpMVpara(-1, Xp, 1, Bu);
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
                solver.DefineMatrix(ConvDiff);
                solver.Solve(Xu, Bu);
            }

            var temp = new double[Uidx.Length + Pidx.Length];

            for (int i = 0; i < Uidx.Length; i++) {
                temp[Uidx[i]] = Xu[i];
            }

            for (int i = 0; i < Pidx.Length; i++) {
                temp[Pidx[i]] = Xp[i];
            }

            X.SetV(temp);
        }
    }
}
