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
using ilPSP.Tracing;
using ilPSP.Utils;
using BoSSS.Platform;

namespace NSE_SIMPLE {

    /// <summary>
    /// Utility class which helps implementing solvers for substeps of SIMPLE iterations.
    /// </summary>
    public abstract class SIMPLESolver {

        ISparseSolver m_Solver;

        /// <summary>
        /// Configuration of SIMPLE algorithm
        /// </summary>
        protected SolverConfiguration m_solverConf;

        //The coefficient matrix and the right-hand side of this solver.
        //One SIMPLESolver can be used to solve for different right-hand sides,
        //but there is always only _one_ coefficient matrix per SIMPLESolver.
        MsrMatrix m_SolverMatrix;
        IList<double> m_Rhs;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConf">Configuration SIMPLE</param>
        /// <param name="_sparseSolver">Configuration sparse solver</param>        
        protected SIMPLESolver(SolverConfiguration solverConf, ISparseSolver _sparseSolver) {

            m_solverConf = solverConf;

            //Create sparse solver
            m_Solver = _sparseSolver;
            if (m_Solver == null)
                throw new ApplicationException("Sorry - unknown error creating solver.");
        }

        /// <summary>
        /// The solver matrix can be defined only once during liftime of this object.
        /// </summary>
        bool MatrixIsDefined = false;

        /// <summary>
        /// Implement this method to define the coefficient matrix of this solver.
        /// </summary>
        /// <param name="dt">time step size</param>
        /// <returns>The coefficient matrix</returns>
        protected abstract MsrMatrix DefineMatrix(double dt);

        /// <summary>
        /// Implement this method to define the right-hand side of this solver.
        /// One <see cref="SIMPLESolver"/> can be solved for different right-hand sides.
        /// </summary>
        /// <param name="dt">time step size</param>
        /// <param name="SpatialComponent">Spatial component - if needed.</param>
        /// <returns>
        /// The right hand side of this solver,
        /// which may depend on <paramref name="SpatialComponent"/>.
        /// </returns>
        protected abstract IList<double> DefineRhs(double dt, int SpatialComponent);

        /// <summary>
        /// Solves the system for <paramref name="Unknowns"/>.
        /// </summary>
        /// <param name="Unknowns">The unknowns to be solved for.</param>
        /// <param name="dt">time step size</param>
        /// <param name="SpatialComponent">Spatial component - this parameter is optional.</param>
        public virtual void Solve(IList<double> Unknowns, double dt, int SpatialComponent = -1) {
            using (new FuncTrace()) {
                //Define and set matrix - can be called only once during lifetime of this object.
                if (!MatrixIsDefined) {
                    m_SolverMatrix = DefineMatrix(dt);
                    m_Solver.DefineMatrix(m_SolverMatrix);
                    MatrixIsDefined = true;
                }

                //Define and set Rhs                        
                m_Rhs = DefineRhs(dt, SpatialComponent);


                //Some checks
                if ((m_SolverMatrix.RowPartitioning.LocalLength != Unknowns.Count) || (m_SolverMatrix.RowPartitioning.LocalLength != m_Rhs.Count))
                    throw new ArgumentException("Mismatch of dimensions!");

                //Solve            
                SolverResult res = m_Solver.Solve(Unknowns, m_Rhs);

                if (m_solverConf.Control.PrintLinerSolverResults && (m_solverConf.MPIRank == 0)) {
                    Console.WriteLine("Linear solver results:");
                    Console.WriteLine("Converged: {0}, NoOfIterations: {1}, Runtime: {2} sec", res.Converged.ToString(), res.NoOfIterations, res.RunTime.TotalSeconds);


                }
            }
        }

        /// <summary>
        /// Dispose _sparseSolver object.
        /// </summary>
        public void Dispose() {
            IDisposable dsolver = m_Solver as IDisposable;
            if (dsolver != null)
                dsolver.Dispose();
        }
    }
}
