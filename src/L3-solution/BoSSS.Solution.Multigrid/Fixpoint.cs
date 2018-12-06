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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using ilPSP.Utils;
using BoSSS.Solution.Multigrid;
using BoSSS.Platform;
using BoSSS.Foundation.XDG;
using MPI.Wrappers;
using System.Diagnostics;
using ilPSP.Tracing;

namespace BoSSS.Solution.Multigrid {


    /// <summary>
    /// Fix-point iteration, should have linear convergence.
    /// </summary>
    public class FixpointIterator : NonlinearSolver {


        public FixpointIterator(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq,
            MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig)
            : base(__AssembleMatrix, __AggBasisSeq, __MultigridOperatorConfig) //
        {
        }


        public int MaxIter = 400;
        public int MinIter = 4;
        public double ConvCrit = 1e-9;


        public ISolverSmootherTemplate m_LinearSolver;

        public string m_SessionPath;

        public double UnderRelax = 1;


        ///// <summary>
        ///// delays the evaluation of the coupled instance (e.g. Level-set) during update
        ///// </summary>
        //public event Action DelayCoupledIteration;

        /// <summary>
        /// delegate for checking the convergence criteria of the coupled iteration
        /// </summary>
        /// <returns></returns>
        public delegate bool CoupledConvergenceReached();

        public CoupledConvergenceReached CoupledIteration_Converged;

        public int MaxCoupledIter = 1000;


        public delegate int IterationCounter(int NoIter, ref int coupledIter);

        public IterationCounter Iteration_Count;


        //bool solveVelocity = true;

        //double VelocitySolver_ConvergenceCriterion = 1e-5;

        //double StressSolver_ConvergenceCriterion = 1e-5;


        override public void SolverDriver<S>(CoordinateVector SolutionVec, S RHS) {
            using (var tr = new FuncTrace()) {

                // initial guess and its residual
                // ==============================
                double[] Solution, Residual;

                using (new BlockTrace("Slv Init", tr)) {
                    base.Init(SolutionVec, RHS, out Solution, out Residual);
                }
                double[] Correction = new double[Solution.Length];
                double ResidualNorm = Residual.L2NormPow2().MPISum().Sqrt();
                int NoOfIterations = 0;
                if (m_LinearSolver.GetType() == typeof(SoftGMRES))
                    ((SoftGMRES)m_LinearSolver).m_SessionPath = m_SessionPath;

                OnIterationCallback(NoOfIterations, Solution.CloneAs(), Residual.CloneAs(), this.CurrentLin);

                if (CoupledIteration_Converged == null) 
                    CoupledIteration_Converged = delegate () {
                        return true;
                    };

                int NoOfCoupledIteration = 0;
                if(Iteration_Count == null)
                    Iteration_Count = delegate (int NoIter, ref int coupledIter) {
                        return NoIter + 1;
                    };

                //int[] Velocity_idx = SolutionVec.Mapping.GetSubvectorIndices(false, 0, 1, 2);
                //int[] Stresses_idx = SolutionVec.Mapping.GetSubvectorIndices(false, 3, 4, 5);

                //int[] Velocity_fields = new int[] { 0, 1, 2 };
                //int[] Stress_fields = new int[] { 3, 4, 5 };

                //int NoCoupledIterations = 10;

                // iterate...
                // ==========
                //int NoOfMainIterations = 0;
                using (new BlockTrace("Slv Iter", tr)) {
                    while ((!(ResidualNorm < ConvCrit && CoupledIteration_Converged()) && NoOfIterations < MaxIter && NoOfCoupledIteration < MaxCoupledIter) || (NoOfIterations < MinIter)) {
                        NoOfIterations = Iteration_Count(NoOfIterations, ref NoOfCoupledIteration);
                        //Console.WriteLine("NoOfIterations = {0}", NoOfIterations);

                        //DirectSolver ds = new DirectSolver();
                        //ds.Init(this.CurrentLin);
                        //double L2_Res = Residual.L2Norm();
                        this.m_LinearSolver.Init(this.CurrentLin);
                        Correction.ClearEntries();
                        if (Correction.Length != Residual.Length)
                            Correction = new double[Residual.Length];
                        this.m_LinearSolver.Solve(Correction, Residual);

                        //if (NoOfIterations > NoCoupledIterations)
                        //{
                        //    if (solveVelocity)
                        //    {
                        //        Console.WriteLine("stress correction = 0");
                        //        foreach (int idx in Stresses_idx)
                        //        {
                        //            Correction[idx] = 0.0;
                        //        }
                        //    }
                        //    else
                        //    {
                        //        Console.WriteLine("velocity correction = 0");
                        //        foreach (int idx in Velocity_idx)
                        //        {
                        //            Correction[idx] = 0.0;
                        //        }
                        //    }
                        //}

                        // Residual may be invalid from now on...
                        Solution.AccV(UnderRelax, Correction);

                        // transform solution back to 'original domain'
                        // to perform the linearization at the new point...
                        // (and for Level-Set-Updates ...)
                        this.CurrentLin.TransformSolFrom(SolutionVec, Solution);

                        // update linearization
                        base.Update(SolutionVec.Mapping.Fields, ref Solution);

                        // residual evaluation & callback
                        base.EvalResidual(Solution, ref Residual);
                        ResidualNorm = Residual.L2NormPow2().MPISum().Sqrt();

                        //if (NoOfIterations > NoCoupledIterations)
                        //{

                        //    double coupledL2Res = 0.0;
                        //    if (solveVelocity)
                        //    {
                        //        foreach (int idx in Velocity_idx)
                        //        {
                        //            coupledL2Res += Residual[idx].Pow2();
                        //        }
                        //    }
                        //    else
                        //    {
                        //        foreach (int idx in Stresses_idx)
                        //        {
                        //            coupledL2Res += Residual[idx].Pow2();
                        //        }
                        //    }
                        //    coupledL2Res = coupledL2Res.Sqrt();

                        //    Console.WriteLine("coupled residual = {0}", coupledL2Res);

                        //    if (solveVelocity && coupledL2Res < this.VelocitySolver_ConvergenceCriterion)
                        //    {
                        //        Console.WriteLine("SolveVelocity = false");
                        //        this.solveVelocity = false;
                        //    }
                        //    else if (!solveVelocity && coupledL2Res < this.StressSolver_ConvergenceCriterion)
                        //    {
                        //        Console.WriteLine("SolveVelocity = true");
                        //        this.solveVelocity = true;
                        //    }
                        //}

                        OnIterationCallback(NoOfIterations, Solution.CloneAs(), Residual.CloneAs(), this.CurrentLin);
                    }
                }
            }

        }
    }


}
