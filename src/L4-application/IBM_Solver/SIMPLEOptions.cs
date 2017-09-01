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

using BoSSS.Solution.Control;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;

namespace BoSSS.Application.IBM_Solver {
    class SIMPLEOptions {

        public SIMPLEOptions() {
            this.relax_v = 1.0;
            this.relax_p = 1.0;
            this.Option_Approximation_Predictor = ApproxPredictor.MassMatrix;
            this.CorrrectionMode = SimpleCorrectionMode.classic;
            this.PotentialSolver_UseMassMatrix = true;
        }
        
        public double relax_p {
            get;
            set;
        }

        public double relax_v {
            get;
            set;
        }

        public ApproxPredictor Option_Approximation_Predictor {
            get;
            set;
        }

        /// <summary>
        /// +inf means steady solver!
        /// </summary>
        public double dt = double.PositiveInfinity;



        //private ISparseSolver m_ViscousSolver;
        public ISparseSolver ViscousSolver {
            get {
                return new PARDISOSolver();
            }
        }

        public ISparseSolver PressureSolver {
            get {
                return new PARDISOSolver();
            }
        }


        public SimpleCorrectionMode CorrrectionMode {
            get;
            set;
        }

        /// <summary>
        /// If true, an additional pressure solution is performed.
        /// </summary>
        public bool RunSIMPLER {
            get;
            set;
        }

        /// <summary>
        /// Only used by the <see cref="PotentialSolverPrecond"/>;
        /// </summary>
        public bool PotentialSolver_UseMassMatrix {
            get;
            set;
        }
    }




    public enum SimpleCorrectionMode {
        classic,

        ResidualMinimization,

        ResidualMinimization_perSpecies
    }
}
