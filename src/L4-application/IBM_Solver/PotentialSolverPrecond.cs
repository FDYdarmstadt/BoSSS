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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Quadrature;
using ilPSP.LinSolvers.PARDISO;
using System.Diagnostics;
using NUnit.Framework;
using System.IO;
using System.Globalization;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using ilPSP;

namespace BoSSS.Application.IBM_Solver {
    
    /// <summary>
    /// A preconditioner, based on the idea of a potential solver;
    /// the idea is to have an preconditioner alternative to <see cref="SIMPLE"/>,
    /// in situations where the SIMPLE lacks good convergence.
    /// </summary>
    class PotentialSolverPrecond : ISolverSmootherTemplate, ISolverWithCallback {


        MultigridOperator m_MgOp;

        public void Init(MultigridOperator op) {
            this.m_MgOp = op;
            int D = this.LsTrk.GridDat.SpatialDimension;

            int[] VelVarIdx = D.ForLoop(d => d);

            this.USubMatrixIdx_Row = this.m_MgOp.Mapping.GetSubvectorIndices(VelVarIdx);
            this.PSubMatrixIdx_Row = this.m_MgOp.Mapping.GetSubvectorIndices(new int[] { D });
            
            this.ExtractMatrices();
        }



        void DoCallBack<V>(double[] Velocity, double[] Pressure, V RHS)
            where V : IList<double> //
        {
            if(this.IterationCallback != null) {
                int LL = this.m_MgOp.Mapping.LocalLength;
                double[] sol = new double[LL];
                double[] res = new double[LL];

                res.SetV(RHS);
                sol.AccV(1.0, Velocity, USubMatrixIdx_Row, default(int[]));
                sol.AccV(1.0, Pressure, PSubMatrixIdx_Row, default(int[]));

                this.m_MgOp.OperatorMatrix.SpMV(-1.0, sol, 1.0, res);

                this.IterationCallback(this.NoOfIterations, sol, res, this.m_MgOp);
            }
        }

        /// <summary>
        /// Computes a potential solution (i.e. neglects convective/diffusive part) to the current residual
        /// </summary>
        /// <param name="X">input/output: solution guess</param>
        /// <param name="RHS">RHS of the saddle point problem</param>
        public void Solve<U, V>(U X, V RHS)
            where U : IList<double>
            where V : IList<double> //
        {

            double[] RESI = RHS.ToArray();
            this.m_MgOp.OperatorMatrix.SpMV(-1.0, X, 1.0, RESI);

            var MM = this.m_MgOp.MassMatrix.CloneAs();
            MM.AccEyeSp(-1.0);
            double nrm = MM.InfNorm();

            double[] R1 = new double[USubMatrixIdx_Row.Length];
            double[] R2 = new double[PSubMatrixIdx_Row.Length];
            RESI.GetSubVector(R1, USubMatrixIdx_Row, default(int[]));
            RESI.GetSubVector(R2, PSubMatrixIdx_Row, default(int[]));


            double gamma;
            if(double.IsInfinity(this.m_SIMPLEOptions.dt))
                gamma = 1.0;
            else
                gamma = (1.0 + this.m_SIMPLEOptions.dt) / this.m_SIMPLEOptions.dt;
            
            
            // RHS of Poisson equation
            double[] Poisson_RHS;
            {
                Poisson_RHS = R2;
                if(this.invMM != null) {
                    double[] tmp = new double[R1.Length];
                    this.invMM.SpMVpara(1.0, R1, 0.0, tmp);
                    this.VelocityDiv.SpMVpara(1.0, tmp, -gamma, Poisson_RHS);
                } else {
                    this.VelocityDiv.SpMVpara(1.0, R1, -gamma, Poisson_RHS);
                }
            }

            // LHS of Poisson equation
            if(m_PressureSolver == null) {
                m_PressureSolver = m_SIMPLEOptions.PressureSolver;
                if(m_PressureSolver is PARDISOSolver) {
                    ((PARDISOSolver)m_PressureSolver).CacheFactorization = true;
                }

                MsrMatrix MXCorrector;
                if(invMM != null) {
                    MXCorrector = MsrMatrix.Multiply(this.VelocityDiv, MsrMatrix.Multiply(this.invMM, this.PressureGrad));
                } else {
                    MXCorrector = MsrMatrix.Multiply(this.VelocityDiv, this.PressureGrad);
                }
                MXCorrector.Acc(-1.0, Stab);

                m_PressureSolver.DefineMatrix(MXCorrector);
            }
            
            // solve Poisson equation
            double[] PressureCorr = new double[R2.Length];
            m_PressureSolver.Solve(PressureCorr, Poisson_RHS);
            
            //// compute velocity correction
            //double[] VelocityCorr = R1;
            //this.PressureGrad.SpMV(-1.0 / gamma, PressureCorr, 1.0, VelocityCorr);

            // apply corrections & return
            //X.AccV(1.0, VelocityCorr, USubMatrixIdx_Row, default(int[])); // the velocity correction is crap; better not to add it.
            X.AccV(1.0, PressureCorr, PSubMatrixIdx_Row, default(int[]));

            NoOfIterations++;
        }

        ISparseSolver m_PressureSolver;


        /// <summary>
        /// iteration counter -- equal to number of solver calls since last counter reset.
        /// </summary>
        private int NoOfIterations = 0;

        public int IterationsInNested {
            get {
                return 0;
            }
        }

        public int ThisLevelIterations {
            get {
                return this.NoOfIterations;
            }
        }

        public bool Converged {
            get {
                return m_Converged;
            }
        }

        bool m_Converged = false;

        public void ResetStat() {
            this.m_Converged = false;
            this.NoOfIterations = 0;
        }
        
        LevelSetTracker LsTrk;

        public SIMPLEOptions m_SIMPLEOptions = null;
        
        MsrMatrix VelocityDiv;
        MsrMatrix PressureGrad;
        MsrMatrix Stab;
        MsrMatrix invMM;


        //double PCorrL2Norm;
        private int[] PSubMatrixIdx_Row;
        private int[] USubMatrixIdx_Row;
        
        public PotentialSolverPrecond(LevelSetTracker _LsTrk) {
            this.LsTrk = _LsTrk;
        }

        void ExtractMatrices() {

            // dispose old solver, if required
            // ===============================

            if(this.m_PressureSolver != null) {
                this.m_PressureSolver.Dispose();
                this.m_PressureSolver = null;
            }
            

            // sub-matrices for the Potential Solver
            // =====================================
            int VelocityLength = this.USubMatrixIdx_Row.Length;
            int PressureLength = this.PSubMatrixIdx_Row.Length;


            this.PressureGrad = new MsrMatrix(VelocityLength, PressureLength, 1, 1);
            this.VelocityDiv = new MsrMatrix(PressureLength, VelocityLength, 1, 1);
            this.Stab = new MsrMatrix(PressureLength, PressureLength, 1, 1);

            var WholeSystemMatrix = this.m_MgOp.OperatorMatrix;

            WholeSystemMatrix.WriteSubMatrixTo(PressureGrad, USubMatrixIdx_Row, default(int[]), PSubMatrixIdx_Row, default(int[]));
            WholeSystemMatrix.WriteSubMatrixTo(VelocityDiv, PSubMatrixIdx_Row, default(int[]), USubMatrixIdx_Row, default(int[]));
            WholeSystemMatrix.WriteSubMatrixTo(Stab, PSubMatrixIdx_Row, default(int[]), PSubMatrixIdx_Row, default(int[]));

            

            // inverse mass matrix 
            // ===================

            if(this.m_SIMPLEOptions.PotentialSolver_UseMassMatrix) {

                // extract the mass-matrix block for the velocity part
                // ---------------------------------------------------
                MsrMatrix MM = new MsrMatrix(this.PressureGrad.RowPartitioning, this.VelocityDiv.ColPartition);
                this.m_MgOp.MassMatrix.WriteSubMatrixTo(MM, this.USubMatrixIdx_Row, default(int[]), this.USubMatrixIdx_Row, default(int[]));

                // invert mass matrix
                // ------------------

                int D = this.LsTrk.GridDat.SpatialDimension;
                this.invMM = new MsrMatrix(MM.RowPartitioning, MM.ColPartition);

                MultidimensionalArray Block = new MultidimensionalArray(2);
                MultidimensionalArray InvBlock = new MultidimensionalArray(2);

                int iRow0 =  MM.RowPartitioning.i0;
                int JAGG = this.m_MgOp.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                int[] DegreeS = this.m_MgOp.Mapping.DgDegree;
                for(int jagg = 0; jagg < JAGG; jagg++) { // loop over aggregate cells...
                    for(int d = 0; d < D; d++) { // loop over velocity components...
                        int N = this.m_MgOp.Mapping.AggBasis[d].GetLength(jagg, DegreeS[d]);
                        if(Block.GetLength(0) != N) {
                            Block.Allocate(N, N);
                            InvBlock.Allocate(N, N);
                        }

                        for(int n = 0; n < N; n++) {
                            CheckMatrix(MM, iRow0, N, iRow0 + n);
                            for (int m = 0; m < N; m++) {
                                Block[n, m] = MM[iRow0 + n, iRow0 + m];
                            }
                        }

                        Block.InvertTo(InvBlock);

                        for(int n = 0; n < N; n++) {
                            for(int m = 0; m < N; m++) {
                                this.invMM[iRow0 + n, iRow0 + m] = InvBlock[n, m];
                            }
                        }

                        iRow0 += N;
                    }
                }

                Debug.Assert(iRow0 == MM.RowPartitioning.iE);

#if DEBUG
                var CheckMX = MM * invMM;
                double TRESH = Math.Max(MM.InfNorm(), invMM.InfNorm()) * 1.0e-10;
                for(int iRow = CheckMX.RowPartitioning.i0; iRow < CheckMX.RowPartitioning.iE; iRow++) {
                    if(Math.Abs(CheckMX.GetDiagonalElement(iRow) - 1.0) > TRESH)
                        throw new ArithmeticException("AapproxInverse is not the Inverse of the Aapprox-Matrix");
                }
#endif

            }
        }

        [Conditional("DEBUG")]
        private static void CheckMatrix(MsrMatrix MM, int iRow0, int N, int iRow) {
            int[] Cols = null;
            double[] Vals = null;
            int LR = MM.GetRow(iRow, ref Cols, ref Vals);
            int cMin = int.MaxValue;
            int cMax = int.MinValue;
            for (int lr = 0; lr < LR; lr++) {
                if (Vals[lr] != 0.0) {
                    cMin = Math.Min(cMin, Cols[lr]);
                    cMax = Math.Max(cMax, Cols[lr]);
                }
            }
            Debug.Assert(cMin >= iRow0);
            Debug.Assert(cMax < iRow0 + N);
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
    }
}
