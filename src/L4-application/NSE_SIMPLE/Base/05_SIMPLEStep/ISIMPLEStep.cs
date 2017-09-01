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
using BoSSS.Solution;
using BoSSS.Foundation.IO;
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// Status of current SIMPLEStep.
    /// </summary>
    public struct SIMPLEStepStatus {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="SolverConf"></param>
        public SIMPLEStepStatus(SIMPLEControl SolverConf) : this() {
            switch (SolverConf.PhysicsMode) {
                case PhysicsMode.Incompressible:
                case PhysicsMode.Multiphase:
                    this.Timestep = new TimestepNumber(0, 0);
                    break;
                case PhysicsMode.LowMach:
                    this.Timestep = new TimestepNumber(0, 1, 0);
                    break;
                default:
                    throw new NotImplementedException();
            }           
        }

        public TimestepNumber Timestep {
            get;
            private set;
        }

        /// <summary>
        /// Number of SIMPLE iteration of current time step.
        /// </summary>
        public int SIMPLEStepNo {
            get {
                if (Timestep.Length == 2)
                    return Timestep[1];
                else if (Timestep.Length == 3)
                    return Timestep[2];
                else
                    throw new ApplicationException();
            }
        }

        /// <summary>
        /// Consecutively counts all inner iterations (for the unsteady solver) for a nicer display of the residua in the logfile.
        /// </summary>
        public int SIMPLEStepNoTotal {
            get;
            private set;
        }

        /// <summary>
        /// Counter for time steps, which reached <see cref="SIMPLEControl.MaxNoSIMPLEsteps"/>.
        /// </summary>
        public int CntMaxNoSIMPLEsteps {
            get;
            set;
        }

        public bool IsConverged {
            get;
            set;
        }

        public bool TerminateSIMPLE {
            get;
            set;
        }        

        public bool SaveStep {
            get;
            set;
        }

        public void NextSIMPLEIteration() {
            if (Timestep.Length == 2)
                Timestep = Timestep.NextIteration();
            else if (Timestep.Length == 3)
                Timestep = Timestep.NextSubIteration();
            else
                throw new ApplicationException();

            SIMPLEStepNoTotal++;

            IsConverged = false;
            TerminateSIMPLE = false;            
            SaveStep = false;
        }

        public void NextSegregatedStep() {
            if (Timestep.Length == 3) {
                Timestep = Timestep.NextIteration();
            } else {
                throw new ApplicationException();
            }
        }

        public void NextTimestep() {
            Timestep = Timestep.NextTimestep();
        }
    }

    /// <summary>
    /// Interface for SIMPLEStep, which includes all operations for one SIMPLE iteration.
    /// </summary>
    public interface ISIMPLEStep : IDisposable {

        /// <summary>
        /// Performs one SIMPLE iteration.
        /// </summary>
        /// <param name="SIMPLEStatus"></param>
        /// <param name="dt">time step size</param>
        /// <param name="ResLogger"></param>
        void OverallIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger);
    }
}
