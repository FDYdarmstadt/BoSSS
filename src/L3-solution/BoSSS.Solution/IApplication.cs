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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution {

    /// <summary>
    /// Public interface of all applications.
    /// </summary>
    public interface IApplication : IDisposable {
        /// <summary>
        /// Information about the currently active session.
        /// </summary>
        SessionInfo CurrentSessionInfo {
            get;
        }

        /// <summary>
        /// Interface to the database driver
        /// </summary>
        IDatabaseDriver DatabaseDriver {
            get;
        }

        /// <summary>
        /// BoSSS grid
        /// </summary>
        IGrid Grid {
            get;
        }

        /// <summary>
        /// extended grid information
        /// </summary>
        IGridData GridData {
            get;
        }

        /// <summary>
        /// All fields, for which IO should be performed.
        /// </summary>
        ICollection<DGField> IOFields {
            get;
        }

        /// <summary>
        /// rank of this process within the MPI world communicator
        /// </summary>
        int MPIRank {
            get;
        }

        /// <summary>
        /// Size of the MPI world communicator
        /// </summary>
        int MPISize {
            get;
        }

        /// <summary>
        /// The query handler holding all queries to be performed during a run.
        /// </summary>
        QueryHandler QueryHandler {
            get;
        }

        /// <summary>
        /// New 'table' to log query results
        /// </summary>
        QueryResultTable QueryResultTable {
            get;
        }

        /// <summary>
        /// the residual logger for this context
        /// </summary>
        ResidualLogger ResLogger {
            get;
        }

        /// <summary>
        /// User configuration input.
        /// </summary>
        AppControl ControlBase {
            get;
        }

        /// <summary>
        /// Initializes the environment of the application
        /// </summary>
        /// <param name="control">
        /// control object
        /// </param>
        void Init(AppControl control);


        /// <summary>
        /// Runs the application in the "solver"-mode. 
        /// </summary>
        void RunSolverMode();

        /// <summary>
        /// This method should be overridden to support automatic numerical stability analysis of the PDE's operator
        /// </summary>
        /// <returns>
        /// Pairs of property name and value, e.g. ConditionNumber and the respective value of the operators Jacobian matrix condition number.
        /// </returns>
        IDictionary<string, double> OperatorAnalysis(OperatorAnalysisConfig config);


        /// <summary>
        /// Assess to optional level-set tracker; typically null for pure DG-apps (non-XDG).
        /// </summary>
        LevelSetTracker LsTrk {
            get;
        }
    }

    /// <summary>
    /// custom config of the operator analysis
    /// </summary>
    public class OperatorAnalysisConfig {

        /// <summary>
        /// Stencil condition number: works also for higher resolutions than <see cref="CalculateGlobalConditionNumbers"/>, scales linearly with number of cells.
        /// </summary>
        public bool CalculateStencilConditionNumbers = true;

        /// <summary>
        /// Global condition number using MATLAB: very expensive, much more than actual solution of the system.
        /// Only attainable for very small systems, maybe up to 10000 DOFs.
        /// </summary>
        public bool CalculateGlobalConditionNumbers = true;
    }


    /// <summary>
    /// Public interface of all applications.
    /// </summary>
    /// <typeparam name="T">
    /// Type of the control file object
    /// </typeparam>
    public interface IApplication<out T> : IApplication
        where T : AppControl, new() {

        /// <summary>
        /// User configuration input.
        /// </summary>
        T Control {
            get;
        }

        
    }
}
