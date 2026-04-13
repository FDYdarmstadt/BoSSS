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
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Solution.CompressibleFlowCommon {

    /// <summary>
    /// Control file for compressible flow simulations
    /// </summary>
    [Serializable]
    public class CompressibleControl : AppControl, ICompressibleControl, ICloneable {

        /// <summary>
        /// Material parameters to be used
        /// </summary>
        [NotNull]
        public IEquationOfState EquationOfState {
            get => CompressibleConfiguration.EquationOfState;
            set => CompressibleConfiguration.EquationOfState = value;
        }

        /// <summary>
        /// The viscosity law to be used, i.e. the variation of the viscosity
        /// due to temperature changes.
        /// </summary>
        public IViscosityLaw ViscosityLaw {
            get => CompressibleConfiguration.ViscosityLaw;
            set => CompressibleConfiguration.ViscosityLaw = value;
        }

        /// <summary>
        /// The configured Mach Number in the far field.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double MachNumber {
            get => CompressibleConfiguration.MachNumber;
            set => CompressibleConfiguration.MachNumber = value;
        }

        /// <summary>
        /// The configured Reynolds number in the far field.
        /// </summary>
        /// <remarks>
        /// This option is ignored if <see cref="DomainType"/> is equal to "Euler"
        /// </remarks>
        public double ReynoldsNumber {
            get => CompressibleConfiguration.ReynoldsNumber;
            set => CompressibleConfiguration.ReynoldsNumber = value;
        }

        /// <summary>
        /// The configured Prandtl number in the far field.
        /// </summary>
        /// <remarks>
        /// This option is ignored if <see cref="DomainType"/> is equal to "Euler"
        /// </remarks>
        public double PrandtlNumber {
            get => CompressibleConfiguration.PrandtlNumber;
            set => CompressibleConfiguration.PrandtlNumber = value;
        }

        /// <summary>
        /// The ratio of a characteristic flow velocity to the velocity of a
        /// gravitational wave.
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public double FroudeNumber {
            get => CompressibleConfiguration.FroudeNumber;
            set => CompressibleConfiguration.FroudeNumber = value;
        }

        /// <summary>
        /// The ratio of the bulk viscosity to the shear viscosity.
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public double ViscosityRatio {
            get => CompressibleConfiguration.ViscosityRatio;
            set => CompressibleConfiguration.ViscosityRatio = value;
        }

        /// <summary>
        /// If set to a positive value, defines the interval (in terms of the
        /// time-step number) between log messages on the console (e.g., to
        /// keep file sizes smaller for long runs)
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public int PrintInterval {
            get => CompressibleConfiguration.PrintInterval;
            set => CompressibleConfiguration.PrintInterval = value;
        }

        [DataMember]
        public CompressibleConfiguration CompressibleConfiguration { get; private set; } = new CompressibleConfiguration();

        ICompressibleConfiguration ICompressibleControl.CompressibleConfiguration => CompressibleConfiguration;

        /// <summary>
        /// %
        /// </summary>
        virtual public Material GetMaterial() {
            return CompressibleConfiguration.GetMaterial();
        }

        /// <summary>
        /// The configured polynomial DG degree of the density
        /// </summary>
        public int DensityDegree {
            get {
                return CompressibleConfiguration.DensityDegree;
            }
        }

        /// <summary>
        /// The configured polynomial DG degree of the momentum components.
        /// </summary>
        public int MomentumDegree {
            get {
                return CompressibleConfiguration.MomentumDegree;
            }
        }

        /// <summary>
        /// The configured polynomial DG degree of the energy
        /// </summary>
        public int EnergyDegree {
            get {
                return CompressibleConfiguration.EnergyDegree;
            }
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
            CompressibleConfiguration.SetVariableDegree(variable, degree);

            FieldOpts.SaveToDBOpt option;
            if (saveToDB) {
                option = FieldOpts.SaveToDBOpt.TRUE;
            } else {
                option = FieldOpts.SaveToDBOpt.FALSE;
            }

            var fieldOpts = new FieldOpts() {
                Degree = degree,
                SaveToDB = option
            };

            if (FieldOptions.ContainsKey(variable)) {
                FieldOptions[variable] = fieldOpts;
            } else {
                FieldOptions.Add(variable, fieldOpts);
            }
        }

        /// <summary>
        /// Dictionary linking field variables (including derived ones) to
        /// the desired polynomial degree
        /// </summary>
        public IReadOnlyDictionary<Variable, int> VariableToDegreeMap {
            get {
                return CompressibleConfiguration.VariableToDegreeMap;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override object Clone() {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Indicates that the residual should be calculated (and saved) every
        /// n-th time-step. If zero, no residual calculation should take place
        /// </summary>
        [InclusiveLowerBound(0)]
        public int ResidualInterval {
            get => CompressibleConfiguration.ResidualInterval;
            set => CompressibleConfiguration.ResidualInterval = value;
        }

        /// <summary>
        /// The type of residual logger to be used, see
        /// <see cref="ResidualLoggerTypes"/>.
        /// </summary>
        public ResidualLoggerTypes ResidualLoggerType {
            get => CompressibleConfiguration.ResidualLoggerType;
            set => CompressibleConfiguration.ResidualLoggerType = value;
        }

        /// <summary>
        /// A mapping between residual variables and the corresponding
        /// termination criteria
        /// </summary>
        public IDictionary<string, double> ResidualBasedTerminationCriteria {
            get => CompressibleConfiguration.ResidualBasedTerminationCriteria;
            set => CompressibleConfiguration.ResidualBasedTerminationCriteria = value;
        }
    }
}
