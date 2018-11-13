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
using BoSSS.Solution;
using CNS.EquationSystem;
using BoSSS.Foundation.IO;
using BoSSS.Solution.CompressibleFlowCommon;

namespace CNS {

    /// <summary>
    /// Generic interface of <see cref="Program"/>.
    /// </summary>
    /// <typeparam name="T">
    /// The type of the control file to be used
    /// </typeparam>
    public interface IProgram<out T> : IApplication<T> where T : CNSControl, new() {

        /// <summary>
        /// A map that determines the active species in some point in the
        /// domain.
        /// </summary>
        ISpeciesMap SpeciesMap {
            get;
        }

        /// <summary>
        /// The employed time stepper
        /// </summary>
        ITimeStepper TimeStepper {
            get;
        }

        /// <summary>
        /// The storage of all current variable values of the current flow
        /// field
        /// </summary>
        CNSFieldSet WorkingSet {
            get;
        }

        /// <summary>
        /// The full operator to be evaluated
        /// </summary>
        Operator FullOperator {
            get;
        }

        /// <summary>
        /// Save the given time-step to the databse
        /// </summary>
        /// <param name="ts"></param>
        /// <param name="phystime"></param>
        void SaveToDatabase(TimestepNumber ts, double phystime);

        /// <summary>
        /// The current (major) time-step number
        /// </summary>
        int TimestepNumber {
            get;
        }
    }
}
