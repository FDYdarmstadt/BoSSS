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
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.LinSolvers.MUMPS;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Platform;
using System.IO;
using ilPSP.Tracing;
using ilPSP.Connectors.Matlab;
using System.Diagnostics;
using ilPSP.LinSolvers.monkey;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.AdvancedSolvers {

    

    /// <summary>
    /// Wrapper around the monkey solver (supports GPU acceleration).
    /// </summary>
    public class MonkeySolver : ISolverSmootherTemplate {

        /// <summary>
        /// 
        /// </summary>
        public enum _whichSolver {

            
            /// <summary>
            /// conjugate gradient solver, see <see cref="ilPSP.LinSolvers.monkey.CG"/>
            /// </summary>
            CG,

            /// <summary>
            /// preconditioned conjugate gradient solver, see <see cref="ilPSP.LinSolvers.monkey.PCG"/>
            /// </summary>
            PCG
            
        }

   

        /// <summary>
        /// Switch between PARDISO and MUMPS.
        /// </summary>
        public _whichSolver WhichSolver = _whichSolver.CG;

        public void Init(MultigridOperator op) {
            using (var tr = new FuncTrace()) {
                var Mtx = op.OperatorMatrix;
                var MgMap = op.Mapping;
                m_MultigridOp = op;

                if (!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");

                m_Mtx = Mtx;
            }
        }

        MultigridOperator m_MultigridOp;

        
        ISparseSolverExt GetSolver(IMutableMatrixEx Mtx) {
            ISparseSolverExt solver;

            
            switch (WhichSolver) {

                case _whichSolver.CG: {
                    var _solver = new CG(); solver = _solver;
                    _solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.Cuda;
                    _solver.MaxIterations = Switcher<int>(_solver.MaxIterations, LinConfig.MaxSolverIterations);
                    _solver.Tolerance = Switcher<double>(_solver.Tolerance, LinConfig.ConvergenceCriterion);
                    break;
                }
                case _whichSolver.PCG: {
                    var _solver = new PCG();solver = _solver;
                    _solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.Cuda;
                    _solver.MaxIterations = Switcher<int>(_solver.MaxIterations, LinConfig.MaxSolverIterations);
                    _solver.Tolerance = Switcher<double>(_solver.Tolerance, LinConfig.ConvergenceCriterion);
                    break;
                }

                default:
                    throw new NotImplementedException();

            }
            

            //solver = new ilPSP.LinSolvers.HYPRE.PCG();

            solver.DefineMatrix(Mtx);

            return solver;
        }

        BlockMsrMatrix m_Mtx;
       

        /// <summary>
        /// %
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (var tr = new FuncTrace()) {
                
                using (var solver = GetSolver(m_Mtx)) {
                    var result = solver.Solve(X, B);
                    this.Converged = result.Converged;
                    this.ThisLevelIterations += result.NoOfIterations;
                }

                
            }
        }


        /// <summary>
        /// %
        /// </summary>
        public int IterationsInNested {
            get { return 0; }
        }

        /// <summary>
        /// 
        /// </summary>
        public int ThisLevelIterations {
            get;
            private set;
        }

        /// <summary>
        /// 
        /// </summary>
        public bool Converged {
            get;
            private set;
        }


        /// <summary>
        /// 
        /// </summary>
        public void ResetStat() {
            ThisLevelIterations = 0;
            Converged = false;
        }

        private static T Switcher<T>(T origin,T setter) {
            T thisreturn;
            if (setter != null) {
                thisreturn = setter;
            } else {
                thisreturn = origin;
            }
            return thisreturn;
        }

        public LinearSolverConfig LinConfig {
            get;
            set;
        }

        /// <summary>
        /// 
        /// </summary>
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        /// <summary>
        /// Release internal memory
        /// </summary>
        public void Dispose() {
            this.m_Mtx = null;
        }

    }
}
