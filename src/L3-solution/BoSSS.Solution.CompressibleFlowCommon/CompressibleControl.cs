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

using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.CompressibleFlowCommon {
    public class CompressibleControl : AppControl, ICloneable {

        /// <summary>
        /// Material parameters to be used
        /// </summary>
        [NotNull]
        public IEquationOfState EquationOfState = IdealGas.Air;

        /// <summary>
        /// The viscosity law to be used, i.e. the variation of the viscosity
        /// due to temperature changes.
        /// </summary>
        public IViscosityLaw ViscosityLaw = new ConstantViscosity();

        /// <summary>
        /// The configured Mach Number in the far field.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double MachNumber;

        /// <summary>
        /// The configured Reynolds number in the far field.
        /// </summary>
        /// <remarks>
        /// This option is ignored if <see cref="DomainType"/> is equal
        /// to "Euler"
        /// </remarks>
        public double ReynoldsNumber;

        /// <summary>
        /// The configured Prandtl number in the far field.
        /// </summary>
        /// <remarks>
        /// This option is ignored if <see cref="DomainType"/> is equal
        /// to "Euler"
        /// </remarks>
        public double PrandtlNumber;

        /// <summary>
        /// The ratio of a characteristic flow velocity to the velocity of a
        /// gravitational wave.
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public double FroudeNumber;

        /// <summary>
        /// The ratio of the bulk viscosity to the shear viscosity.
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public double ViscosityRatio = 0.0;

        /// <summary>
        /// If set to a positive value, defines the interval (in terms of the
        /// time-step number) between log messages on the console (e.g., to
        /// keep file sizes smaller for long runs)
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public int PrintInterval = 1;

        /// <summary>
        /// %
        /// </summary>
        virtual public Material GetMaterial() {
            return new Material(EquationOfState, ViscosityLaw, MachNumber, ReynoldsNumber, PrandtlNumber, FroudeNumber, ViscosityRatio);
        }

        /// <summary>
        /// Utility function to add a new variable to the solver.
        /// </summary>
        /// <param name="variable">
        /// The variable to be added
        /// </param>
        /// <param name="degree">
        /// The desired polynomial degree of the variable
        /// </param>
        /// <param name="saveToDB">
        /// Bool indicating whether the given variable shall be saved to the
        /// database in each saved time-step
        /// </param>
        public void AddVariable(Variable variable, int degree, bool saveToDB = true) {
            variableFields.Add(variable, degree);

            FieldOpts.SaveToDBOpt option;
            if (saveToDB) {
                option = FieldOpts.SaveToDBOpt.TRUE;
            } else {
                option = FieldOpts.SaveToDBOpt.FALSE;
            }

            FieldOptions.Add(variable, new FieldOpts() {
                Degree = degree,
                SaveToDB = option
            });
        }

        /// <summary>
        /// Backing field for <see cref="VariableFields"/>
        /// </summary>
        protected Dictionary<Variable, int> variableFields = new Dictionary<Variable, int>();

        /// <summary>
        /// Dictionary linking field variables (including derived ones) to
        /// the desired polynomial degree
        /// </summary>
        public IReadOnlyDictionary<Variable, int> VariableFields {
            get {
                return variableFields;
            }
        }

        public virtual object Clone() {
            throw new NotImplementedException();
        }
    }
}