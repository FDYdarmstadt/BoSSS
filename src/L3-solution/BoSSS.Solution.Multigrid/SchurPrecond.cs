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

namespace BoSSS.Solution.Multigrid
{
    public class SchurPrecond : ISolverSmootherTemplate, ISolverWithCallback
    {
        public int IterationsInNested {
            get {
                return 0;
                //throw new NotImplementedException();
            }
        }

        public int ThisLevelIterations {
            get {
                return 0;
                //throw new NotImplementedException();
            }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get {
                return null;
                //throw new NotImplementedException();
            }
            set {
             
                //throw new NotImplementedException();
            }
        }

        MultigridOperator m_mgop;

        BlockMsrMatrix Mtx;

        MsrMatrix P;
        MsrMatrix ConvDiff, pGrad, divVel, ConvDiffPoissonMtx, SchurMtx, PoissonMtx, PoissonMtx_T, PoissonMtx_H, SchurConvMtx, invVelMassMatrix, invVelMassMatrixSqrt, simpleSchur, velMassMatrix;
        int[] Uidx, Pidx;

        public enum SchurOptions { exact = 1, decoupledApprox = 2, SIMPLE = 3 }

        public bool ApproxScaling = false;

        public SchurOptions SchurOpt = SchurOptions.exact;

        public void Init(MultigridOperator op)
        {
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
            divVel = new MsrMatrix(Ppart, Upart, 1, 1);
            var PxP = new MsrMatrix(Ppart, Ppart, 1, 1);

            M.AccSubMatrixTo(1.0, ConvDiff, Uidx, default(int[]), Uidx, default(int[]));
            M.AccSubMatrixTo(1.0, pGrad, Uidx, default(int[]), Pidx, default(int[]));
            M.AccSubMatrixTo(1.0, divVel, Pidx, default(int[]), Uidx, default(int[]));
            M.AccSubMatrixTo(1.0, PxP, Pidx, default(int[]), Pidx, default(int[]));

            Mtx = M;

            int L = M.RowPartitioning.LocalLength;

            int i0 = Mtx.RowPartitioning.i0;

            P = new MsrMatrix(Mtx);
            P.Clear();

            // Debugging output
            //ConvDiff.SaveToTextFileSparse("ConvDiff");
            //divVel.SaveToTextFileSparse("divVel");
            //pGrad.SaveToTextFileSparse("pGrad");
            //PxP.SaveToTextFileSparse("PxP");


            velMassMatrix = new MsrMatrix(Upart, Upart, 1, 1);
            op.MassMatrix.AccSubMatrixTo(1.0, velMassMatrix, Uidx, default(int[]), Uidx, default(int[]));

            switch (SchurOpt)
            {
                case SchurOptions.exact:
                    {
                        // Building complete Schur and Approximate Schur
                        MultidimensionalArray Poisson = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
                        MultidimensionalArray SchurConvPart = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
                        MultidimensionalArray Schur = MultidimensionalArray.Create(Pidx.Length, Pidx.Length);
                        using (BatchmodeConnector bmc = new BatchmodeConnector())
                        {
                            bmc.PutSparseMatrix(ConvDiff, "ConvDiff");
                            bmc.PutSparseMatrix(velMassMatrix, "MassMatrix");
                            bmc.PutSparseMatrix(divVel, "divVel");
                            bmc.PutSparseMatrix(pGrad, "pGrad");
                            bmc.Cmd("Qdiag = diag(diag(MassMatrix))");
                            bmc.Cmd("invT= inv(Qdiag)");
                            bmc.Cmd("Poisson = full(invT)*pGrad");
                            bmc.Cmd("ConvPart = ConvDiff*Poisson");
                            bmc.Cmd("ConvPart= full(invT)*ConvPart");
                            bmc.Cmd("ConvPart= divVel*ConvPart");
                            bmc.Cmd("Poisson = divVel*Poisson");
                            bmc.Cmd("ConvDiffInv = inv(full(ConvDiff))");
                            bmc.Cmd("Schur = divVel*ConvDiffInv");
                            bmc.Cmd("Schur = Schur*pGrad");
                            bmc.GetMatrix(Poisson, "Poisson");
                            bmc.GetMatrix(SchurConvPart, "ConvPart");
                            bmc.GetMatrix(Schur, "-Schur");
                            bmc.Execute(false);
                        }
                        PoissonMtx_T = Poisson.ToMsrMatrix();
                        PoissonMtx_H = Poisson.ToMsrMatrix();
                        SchurConvMtx = SchurConvPart.ToMsrMatrix();
                        SchurMtx = Schur.ToMsrMatrix();
                        SchurMtx.Acc(PxP, 1);

                        ConvDiff.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Uidx);
                        pGrad.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Pidx);
                        SchurMtx.AccSubMatrixTo(1.0, P, default(int[]), Pidx, default(int[]), Pidx);
                        return;
                    }
                case SchurOptions.decoupledApprox:
                    {
                        // Do assembly for approximate Schur inverse
                        invVelMassMatrix = velMassMatrix.CloneAs();
                        invVelMassMatrix.Clear();
                        invVelMassMatrixSqrt = invVelMassMatrix.CloneAs();
                        for (int i = 0; i < velMassMatrix.NoOfCols; i++)
                        {
                            if (ApproxScaling)
                            {
                                invVelMassMatrix.SetDiagonalElement(i, 1 / (velMassMatrix[i, i]));
                                invVelMassMatrixSqrt.SetDiagonalElement(i, 1 / (Math.Sqrt(velMassMatrix[i, i])));
                            }
                            else
                            {
                                invVelMassMatrix.SetDiagonalElement(i, 1);
                                invVelMassMatrixSqrt.SetDiagonalElement(i, 1);
                            }
                        }

                        //invVelMassMatrix.SaveToTextFileSparse("invVelMassMatrix");
                        //velMassMatrix.SaveToTextFileSparse("velMassMatrix");


                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(ConvDiff, pGrad);
                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(divVel, ConvDiffPoissonMtx);

                        // Inverse of mass matrix in Matlab
                        //MultidimensionalArray temp = MultidimensionalArray.Create(Uidx.Length, Uidx.Length);
                        //using (BatchmodeConnector bmc = new BatchmodeConnector())
                        //{
                        //    bmc.PutSparseMatrix(velMassMatrix, "velMassMatrix");
                        //    bmc.Cmd("invVelMassMatrix = inv(full(velMassMatrix))");
                        //    bmc.GetMatrix(temp, "invVelMassMatrix");
                        //    bmc.Execute(false);
                        //}
                        //invVelMassMatrix = temp.ToMsrMatrix();

                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(ConvDiffPoissonMtx, PoissonMtx);
                        //ConvDiffPoissonMtx = MsrMatrix.Multiply(PoissonMtx, ConvDiffPoissonMtx);

                        //ConvDiff.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Uidx);
                        //pGrad.AccSubMatrixTo(1.0, P, default(int[]), Uidx, default(int[]), Pidx);
                        //ConvDiffPoissonMtx.AccSubMatrixTo(1.0, P, default(int[]), Pidx, default(int[]), Pidx);

                        //op.MassMatrix.SaveToTextFileSparse("MassMatrix");
                        //velMassMatrix.SaveToTextFileSparse("velMassMatrix2");


                        // Possion scaled by inverse of the velocity mass matrix 
                        PoissonMtx_T = MsrMatrix.Multiply(invVelMassMatrix, pGrad); 
                        PoissonMtx_T = MsrMatrix.Multiply(divVel, PoissonMtx_T);
                        ////PoissonMtx_T.Acc(PxP, 1); // p.379

                        // Poisson scaled by sqrt of inverse of velocity mass matrix
                        PoissonMtx_H = MsrMatrix.Multiply(invVelMassMatrixSqrt, pGrad);
                        PoissonMtx_H = MsrMatrix.Multiply(divVel, PoissonMtx_H);
                        //PoissonMtx_H.Acc(PxP, 1); // p.379
                        return;
                    }
                case SchurOptions.SIMPLE:
                    {

                        var invdiag_ConvDiff = ConvDiff.CloneAs();
                        invdiag_ConvDiff.Clear();
                        for (int i = 0; i < ConvDiff.NoOfCols; i++)
                        {
                            invdiag_ConvDiff[i, i] = 1 / ConvDiff[i, i];
                        }

                        simpleSchur = MsrMatrix.Multiply(invdiag_ConvDiff, pGrad);
                        simpleSchur = MsrMatrix.Multiply(divVel, simpleSchur);

                        return;
                    }

                default:
                    throw new NotImplementedException("SchurOption");
            }


            //var ConvDiffInvMtx = ConvDiffInv.ToMsrMatrix();


            //// x= inv(P)*b !!!!! To be done with approximate Inverse
            // P.SpMV(1, B, 0, X);
        }

        public void ResetStat()
        {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>
        {

            //using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
            //    solver.DefineMatrix(P);
            //    solver.Solve(X, B);
            //}

            switch (SchurOpt)
            {
                case SchurOptions.exact:
                    {
                        // Directly invert Preconditioning Matrix
                        using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
                        {
                            solver.DefineMatrix(P);
                            solver.Solve(X, B);
                        }
                        return;
                    }

                case SchurOptions.decoupledApprox:
                    {
                        SolveSubproblems(X, B);
                        return;
                    }

                case SchurOptions.SIMPLE:
                    {
                        SolveSIMPLE(X, B);
                        return;
                    }

            }

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
            where V : IList<double>
        {

            var Bu = new double[Uidx.Length];
            var Xu = Bu.CloneAs();
            Bu = B.GetSubVector(Uidx, default(int[]));
            var Bp = new double[Pidx.Length];
            var Xp = Bp.CloneAs();
            Bp = B.GetSubVector(Pidx, default(int[]));

            Xu = X.GetSubVector(Uidx, default(int[]));
            Xp = X.GetSubVector(Pidx, default(int[]));

            ApproxAndSolveSchur(Xp, Bp);

            // Solve ConvDiff*w=v-q*pGrad
            pGrad.SpMVpara(-1, Xp, 1, Bu);

            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(ConvDiff);
                solver.Solve(Xu, Bu);
            }

            var temp = new double[Uidx.Length + Pidx.Length];

            for (int i = 0; i < Uidx.Length; i++)
            {
                temp[Uidx[i]] = Xu[i];
            }

            for (int i = 0; i < Pidx.Length; i++)
            {
                temp[Pidx[i]] = Xp[i];
            }

            X.SetV(temp);
        }


        /// <summary>
        /// Approximate the inverse of the Schur matrix and perform two Poisson solves and Matrix-Vector products. Finite elements and Fast Iterative Solvers p.383
        /// </summary>
        public void ApproxAndSolveSchur<U, V>(U Xp, V Bp)
           where U : IList<double>
           where V : IList<double>
        {

            var temp = new double[Xp.Count];
            var sol = new double[pGrad.NoOfRows];

            // Poisson solve
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(PoissonMtx_T);
                solver.Solve(temp, Bp);
            }

            // Schur Convective part with scaling
            pGrad.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            sol.Clear();
            invVelMassMatrix.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            sol.Clear();
            ConvDiff.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            sol.Clear();
            invVelMassMatrixSqrt.SpMVpara(1, temp, 0, sol);

            temp = sol.ToArray();
            divVel.SpMVpara(1, temp, 0, Xp);

            // Poisson solve
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(PoissonMtx_H);
                solver.Solve(Xp, Xp);
            }
        }

        /// <summary>
        /// Solve with a SIMPLE like approcation of the Schur complement
        /// </summary>
        /// <typeparam name="U"></typeparam>
        /// <typeparam name="V"></typeparam>
        /// <param name="X"></param>
        /// <param name="B"></param>
        public void SolveSIMPLE<U, V>(U X, V B)
           where U : IList<double>
           where V : IList<double>
        {

            var Bu = new double[Uidx.Length];
            var Xu = Bu.CloneAs();
            Bu = B.GetSubVector(Uidx, default(int[]));
            var Bp = new double[Pidx.Length];
            var Xp = Bp.CloneAs();
            Bp = B.GetSubVector(Pidx, default(int[]));

            Xu = X.GetSubVector(Uidx, default(int[]));
            Xp = X.GetSubVector(Pidx, default(int[]));

            // Solve SIMPLE Schur
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(simpleSchur);
                solver.Solve(Xp, Bp);
            }


            // Solve ConvDiff*w=v-q*pGrad
            pGrad.SpMVpara(-1, Xp, 1, Bu);
            using (var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver())
            {
                solver.DefineMatrix(ConvDiff);
                solver.Solve(Xu, Bu);
            }

            var temp = new double[Uidx.Length + Pidx.Length];

            for (int i = 0; i < Uidx.Length; i++)
            {
                temp[Uidx[i]] = Xu[i];
            }

            for (int i = 0; i < Pidx.Length; i++)
            {
                temp[Pidx[i]] = Xp[i];
            }

            X.SetV(temp);
        }
    }
}
