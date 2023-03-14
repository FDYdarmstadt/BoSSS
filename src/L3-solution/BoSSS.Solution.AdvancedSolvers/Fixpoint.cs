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
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Platform;
using BoSSS.Foundation.XDG;
using MPI.Wrappers;
using System.Diagnostics;
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Fix-point iteration, should have linear convergence.
    /// </summary>
    public class FixpointIterator : NonlinearSolver {

        /// <summary>
        /// 
        /// </summary>
        public FixpointIterator(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq,
            MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig)
            : base(__AssembleMatrix, __AggBasisSeq, __MultigridOperatorConfig) //
        {
        }


        public int MaxIter = 400;
        public int MinIter = 4;
        public double ConvCrit = 1e-9;
        

        public ISolverFactory m_LinearSolver;

        public string m_SessionPath;

        public double UnderRelax = 1;


        /// <summary>
        /// delegate for checking the convergence criteria of the coupled iteration
        /// </summary>
        /// <returns></returns>
        public delegate bool CoupledConvergenceReached();

        public CoupledConvergenceReached CoupledIteration_Converged;

        public int MinCoupledIterations = 1;
        public int MaxCoupledIter = 10;


        public delegate int IterationCounter(int NoIter, ref int coupledIter);

        public IterationCounter Iteration_Count;


        override public bool SolverDriver<S>(CoordinateVector SolutionVec, S RHS) {
            using (var tr = new FuncTrace()) {

                if (this.ConvCrit <= 0)
                    throw new ArgumentException($"Fixpoint iterator: convergence criterion set to {ConvCrit}: solver will most likely not converge!");

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

                // iterate...
                // ==========
                using (new BlockTrace("Slv Iter", tr)) {
                    bool success = false;
                    
                    while (true) {
                        if(NoOfIterations >= MinIter || NoOfCoupledIteration > MinCoupledIterations) {
                            if(ResidualNorm < ConvCrit && CoupledIteration_Converged()) {
                                success = true;
                                break;
                            }
                            if(NoOfIterations > MaxIter || NoOfCoupledIteration > MaxCoupledIter) {
                                break;
                            }
                        }

                        NoOfIterations = Iteration_Count(NoOfIterations, ref NoOfCoupledIteration);

                        using(var linSlv = this.m_LinearSolver.CreateInstance(this.CurrentLin)) {
                            Correction.ClearEntries();
                            if(Correction.Length != Residual.Length)
                                Correction = new double[Residual.Length];
                            linSlv.Solve(Correction, Residual);
                        }

                        // Residual may be invalid from now on...
                        Solution.AccV(UnderRelax, Correction);

                        // transform solution back to 'original domain'
                        // to perform the linearization at the new point...
                        // (and for Level-Set-Updates ...)
                        this.CurrentLin.TransformSolFrom(SolutionVec, Solution);

                        // update linearization
                        base.Update(SolutionVec.Fields, Solution, 1.0); // SolutionVec -> Solution

                        // residual evaluation & callback
                        base.EvalLinearizedResidual(Solution, ref Residual);
                        ResidualNorm = Residual.L2NormPow2().MPISum().Sqrt();

                        OnIterationCallback(NoOfIterations, Solution.CloneAs(), Residual.CloneAs(), this.CurrentLin);
                    }

                    base.NoOfNonlinearIter = NoOfIterations;
                    return success;
                }
            }

        }
    }


}
