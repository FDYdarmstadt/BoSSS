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
    public class SchurPrec : ISolverSmootherTemplate, ISolverWithCallback {
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
        MsrMatrix ConvDiff, pGrad, divVel;
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

            //CoordinateMapping Umap = this.Velocity.Mapping;
            //CoordinateMapping Pmap = this.Pressure.Mapping;
            int Upart = Uidx.Length;
            int Ppart = Pidx.Length;


            ConvDiff = new MsrMatrix(Upart, Upart, 1, 1);
            pGrad = new MsrMatrix(Upart, Ppart, 1, 1);
            divVel = new MsrMatrix(Ppart, Upart, 1, 1);

            M.AccSubMatrixTo(1.0, ConvDiff, Uidx, default(int[]), Uidx, default(int[]));
            M.AccSubMatrixTo(1.0, pGrad, Uidx, default(int[]), Pidx, default(int[]));
            M.AccSubMatrixTo(1.0, divVel, Pidx, default(int[]), Uidx, default(int[]));

            Mtx = M;

            int L = M.RowPartitioning.LocalLength;

            int i0 = Mtx.RowPartitioning.i0;
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

            MsrMatrix P = new MsrMatrix(Mtx);
            P.Clear();

            // A and pressure
            ConvDiff.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Uidx);
            pGrad.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Pidx);

            // Debugging output
            //ConvDiff.SaveToTextFileSparse("ConvDiff");
            //divVel.SaveToTextFileSparse("divVel");
            //pGrad.SaveToTextFileSparse("pGrad");

            // Building the Schur complement
            MultidimensionalArray ConvDiffInv = MultidimensionalArray.Create(Uidx.Length, Uidx.Length);
            MultidimensionalArray SchurInv = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);

            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutMatrix(ConvDiff.ToFullMatrixOnProc0(), "ConvDiff");
                bmc.PutMatrix(divVel.ToFullMatrixOnProc0(), "divVel");
                bmc.PutMatrix(pGrad.ToFullMatrixOnProc0(), "pGrad");
                bmc.Cmd("ConvDiffInv = inv(ConvDiff)");
                bmc.Cmd("Schur = divVel*ConvDiffInv");
                bmc.Cmd("Schur = Schur*pGrad");
                //bmc.Cmd("SchurInv = inv(Schur)");
                //bmc.GetMatrix(ConvDiffInv, "ConvDiffInv");
                bmc.GetMatrix(SchurInv, "-Schur");

                bmc.Execute(false);
            }


            var SchurInvMtx = SchurInv.ToMsrMatrix();
            var ConvDiffInvMtx = ConvDiffInv.ToMsrMatrix();

            //ConvDiffInvMtx.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Uidx);
            //pGrad.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Pidx);
            SchurInvMtx.AccSubMatrixTo(1.0, P, default(int[]), Pidx, default(int[]), Pidx);
            //// x= inv(P)*b !!!!! To be done with approximate Inverse
            // P.SpMV(1, B, 0, X);

            //// Building the exact inverse of the Preconditioning Matrix x=inv(P)*b
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
                solver.DefineMatrix(P);
                solver.Solve(X, B);
            }
        }
    }
}
