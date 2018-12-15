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

using System.Collections.Generic;
using System.Linq;

namespace ilPSP.LinSolvers {

    /// <summary>
    /// Simplified interfaces to solvers.
    /// </summary>
    static public class SimpleSolversInterface {
        
        /// <summary>
        /// Solves a symmetric, definite system using the conjugate gradient algorithm.
        /// </summary>
        /// <param name="Matrix">the matrix of the linear problem.</param>
        /// <param name="X">Output, (hopefully) the solution.</param>
        /// <param name="B">Input, the right-hand-side.</param>
        /// <param name="MaxIterations"></param>
        /// <param name="Tolerance"></param>
        /// <returns>Actual number of iterations.</returns>
        static public int Solve_CG<V, W>(this IMutableMatrixEx Matrix, V X, W B, int MaxIterations = 100000, double Tolerance = 1.0e-10)
            where V : IList<double>
            where W : IList<double> //
        {
            using (var slv = new ilPSP.LinSolvers.monkey.CG()) {
                slv.MaxIterations = MaxIterations;
                slv.Tolerance = Tolerance;
                slv.DevType = monkey.DeviceType.CPU;

                slv.DefineMatrix(Matrix);

                var SolRes = slv.Solve(X, B.ToArray());
                return SolRes.NoOfIterations;
            }
        }

        /// <summary>
        /// Soves a linear system using a direct solver.
        /// </summary>
        /// <param name="Matrix">the matrix of the linear problem.</param>
        /// <param name="X">Output, (hopefully) the solution.</param>
        /// <param name="B">Input, the right-hand-side.</param>
        /// <returns>Actual number of iterations.</returns>
        static public void Solve_Direct<V, W>(this IMutableMatrixEx Matrix, V X, W B)
            where V : IList<double>
            where W : IList<double> //
        {
            /*
            using (var slv = new ilPSP.LinSolvers.PARDISO.PARDISOSolver()) {
                slv.DefineMatrix(Matrix);
                var SolRes = slv.Solve(X, B.ToArray());
            }
            */
            using(var slv = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
                slv.DefineMatrix(Matrix);
                var SolRes = slv.Solve(X, B.ToArray());
            }
        }
    }
}
