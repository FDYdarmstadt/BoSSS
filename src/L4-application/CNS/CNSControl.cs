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

using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.LoadBalancing;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.ShockCapturing;
using CNS.Source;
using ilPSP;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace CNS {

    /// <summary>
    /// Specialized control file for CNS
    /// </summary>
    [Serializable]
    public class CNSControl : AppControl, ICloneable {

        /// <summary>
        /// Verifies the configuration
        /// </summary>
        public override void Verify() {
            {
                int degree = FieldOptions[Variables.Momentum.xComponent].Degree;

                if (FieldOptions.ContainsKey(Variables.Momentum.yComponent)
                    && FieldOptions[Variables.Momentum.yComponent].Degree != degree) {
                    throw new Exception(
                        "All momentum components must have the same polynomial degree");
                }

                if (FieldOptions.ContainsKey(Variables.Momentum.zComponent)
                    && FieldOptions[Variables.Momentum.zComponent].Degree != degree) {
                    throw new Exception(
                        "All momentum components must have the same polynomial degree");
                }
            }

            if (MachNumber <= 0.0) {
                throw new Exception("Illegal Mach number");
            }
        }

        /// <summary>
        /// Fraction of the theoretically maximum admissible time-step to take.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double CFLFraction;

        /// <summary>
        /// The configured polynomial DG degree of the density
        /// </summary>
        public int DensityDegree {
            get {
                return base.FieldOptions[Variables.Density].Degree;
            }
        }

        /// <summary>
        /// The configured polynomial DG degree of the momentum components.
        /// </summary>
        public int MomentumDegree {
            get {
                return base.FieldOptions[Variables.Momentum.xComponent].Degree;
            }
        }

        /// <summary>
        /// The configured polynomial DG degree of the energy
        /// </summary>
        public int EnergyDegree {
            get {
                return base.FieldOptions[Variables.Energy].Degree;
            }
        }

        /// <summary>
        /// Or-combination of all variables for which initial value definitions
        /// exist
        /// </summary>
        public VariableTypes GetInitialValueVariables() {
            bool conservative = InitialValues_Evaluators.ContainsKey(Variables.Density)
                && InitialValues_Evaluators.ContainsKey(Variables.Energy);
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                conservative &= InitialValues_Evaluators.ContainsKey(Variables.Momentum[d]);
            }

            if (conservative) {
                return VariableTypes.ConservativeVariables;
            }

            bool primitive = InitialValues_Evaluators.ContainsKey(Variables.Density)
                && InitialValues_Evaluators.ContainsKey(Variables.Pressure);
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                primitive &= InitialValues_Evaluators.ContainsKey(Variables.Velocity[d]);
            }

            if (primitive) {
                return VariableTypes.PrimitiveVariables;
            }

            return VariableTypes.None;
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
        private Dictionary<Variable, int> variableFields = new Dictionary<Variable, int>();

        /// <summary>
        /// Dictionary linking field variables (including derived ones) to
        /// the desired polynomial degree
        /// </summary>
        public IReadOnlyDictionary<Variable, int> VariableFields {
            get {
                return variableFields;
            }
        }

        /// <summary>
        /// The type of the domain on which the equations are formulated. For
        /// example, they can be formulated on standard mesh-conformal domain
        /// or truncated by mean of an immersed boundary
        /// </summary>
        public DomainTypes DomainType = DomainTypes.Standard;

        /// <summary>
        /// The active terms of the Navier-Stokes equations.
        /// </summary>
        public Operators ActiveOperators = Operators.Convection;

        /// <summary>
        /// The type of convective fluxes to be used. Currently, only
        /// "HLLC" is supported.
        /// </summary>
        /// <remarks>
        /// This option is ignored if <see cref="DomainType"/> is equal
        /// to "Stokes"
        /// </remarks>
        public ConvectiveFluxTypes ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

        /// <summary>
        /// The type of diffusive fluxes to be used. Currently, only "SIPG"
        /// is supported.
        /// </summary>
        /// <remarks>
        /// This option is ignored if <see cref="DomainType"/> is equal
        /// to "Euler"
        /// </remarks>
        public DiffusiveFluxTypes DiffusiveFluxType = DiffusiveFluxTypes.OptimizedSIPG;

        /// <summary>
        /// The explicit time integration scheme to be used.
        /// </summary>
        public ExplicitSchemes ExplicitScheme = ExplicitSchemes.RungeKutta;

        /// <summary>
        /// The order of the explicit time integration scheme.
        /// </summary>
        [InclusiveLowerBound(1)]
        [InclusiveUpperBound(4)]
        public int ExplicitOrder = 1;

        /// <summary>
        /// The Number of sub-grids for Local Time Stepping
        /// </summary>
        /// This option is only used if <see cref="ExplicitScheme"/> is
        /// equal to LTS
        [InclusiveLowerBound(1)]
        public int NumberOfSubGrids = 1;
        
        /// <summary>
        /// The amount of time steps after which a reclustering is performed if <see cref="ExplicitSchemes.LTS"/>
        /// is used in adaptive mode (<see cref="ReclusteringInterval"/> != 0).
        /// </summary>
        /// <remarks>
        /// This option is only used if <see cref="ExplicitScheme"/> is equal to LTS.
        /// </remarks>
        public int ReclusteringInterval = 0;

        /// <summary>
        /// Enables/Disables the flux correction to obtain a (non-)conservative
        /// <see cref="ExplicitSchemes.LTS"/> scheme. Not avaiblable for adaptive LTS.
        /// </summary>
        /// <remarks>
        /// This option is only used if <see cref="ExplicitScheme"/> is equal to LTS.
        /// </remarks>
        public bool FluxCorrection = true;

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
        /// The viscosity law to be used, i.e. the variation of the viscosity
        /// due to temperature changes.
        /// </summary>
        public IViscosityLaw ViscosityLaw = new ConstantViscosity();

        /// <summary>
        /// Required for the SIPG method. Should be > 1
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double SIPGPenaltyScaling = 1.0;

        /// <summary>
        /// Material parameters to be used
        /// </summary>
        [NotNull]
        public IEquationOfState EquationOfState = IdealGas.Air;

        /// <summary>
        /// The type of residual logger to be used, see
        /// <see cref="ResidualLoggerTypes"/>.
        /// </summary>
        public ResidualLoggerTypes ResidualLoggerType = ResidualLoggerTypes.None;

        /// <summary>
        /// Indicates that the residual should be calculated (and saved) every
        /// n-th time-step. If zero, no residual calculation should take place
        /// </summary>
        [InclusiveLowerBound(0)]
        public int ResidualInterval = 1;

        /// <summary>
        /// If set to a positive value, defines the interval (in terms of the
        /// time-step number) between log messages on the console (e.g., to
        /// keep file sizes smaller for long runs)
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public int PrintInterval = 1;

        /// <summary>
        /// A mapping between residual variables and the corresponding
        /// termination criteria
        /// </summary>
        public IDictionary<string, double> ResidualBasedTerminationCriteria = new Dictionary<string, double>();

        /// <summary>
        /// Custom source terms for the continuity equation
        /// </summary>
        public IList<Func<ISpeciesMap, INonlinearSource>> CustomContinuitySources =
            new List<Func<ISpeciesMap, INonlinearSource>>();

        /// <summary>
        /// Custom source terms for the momentum equations
        /// </summary>
        public IList<Func<ISpeciesMap, INonlinearSource>>[] CustomMomentumSources =
            new IList<Func<ISpeciesMap, INonlinearSource>>[3] {
                new List<Func<ISpeciesMap, INonlinearSource>>(),
                new List<Func<ISpeciesMap, INonlinearSource>>(),
                new List<Func<ISpeciesMap, INonlinearSource>>() };

        /// <summary>
        /// Custom source terms for the energy equation
        /// </summary>
        public IList<Func<ISpeciesMap, INonlinearSource>> CustomEnergySources =
            new List<Func<ISpeciesMap, INonlinearSource>>();

        /// <summary>
        /// Under-relaxation factor for fixed point iterations
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        [InclusiveUpperBound(1.0)]
        public double FixedPointDamping = 1.0;

        /// <summary>
        /// An optional limiter to avoid oscillations near discontinuities
        /// </summary>
        public ILimiter Limiter = null;

        /// <summary>
        /// An optional sensor to detect shocks
        /// </summary>
        public IShockSensor ShockSensor = null;

        /// <summary>
        /// An optional viscosity law to determine the magnitude of the
        /// artificial viscosity if <see cref="ActiveOperators"/> includes
        /// <see cref="Operators.ArtificialViscosity"/>
        /// </summary>
        public IArtificialViscosityLaw ArtificialViscosityLaw = null;

        /// <summary>
        /// Configuration for an option sponge layer defining a non-reflecting
        /// boundary condition.
        /// </summary>
        public SpongeLayerConfig SpongeLayerConfig = null;

        /// <summary>
        /// A classifier that decides which performance class a cell
        /// (currently) belongs to
        /// </summary>
        public ICellClassifier DynamicLoadBalancing_CellClassifier = new IndifferentCellClassifier();

        /// <summary>
        /// Clones this object, but beware: I'm not sure (yet) that I've
        /// covered all cases
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            var clone = (CNSControl)this.MemberwiseClone();
            clone.CustomContinuitySources = new List<Func<ISpeciesMap, INonlinearSource>>();
            clone.CustomContinuitySources.AddRange(this.CustomContinuitySources);

            clone.CustomMomentumSources = new List<Func<ISpeciesMap, INonlinearSource>>[3] {
                new List<Func<ISpeciesMap, INonlinearSource>>(),
                new List<Func<ISpeciesMap, INonlinearSource>>(),
                new List<Func<ISpeciesMap, INonlinearSource>>() };
            for (int d = 0; d < this.CustomMomentumSources.Length; d++) {
                clone.CustomMomentumSources[d].AddRange(this.CustomMomentumSources[d]);
            }

            clone.CustomEnergySources = new List<Func<ISpeciesMap, INonlinearSource>>();
            clone.CustomEnergySources.AddRange(this.CustomEnergySources);

            clone.DynamicLoadBalancing_CellCostEstimatorFactories = new List<Func<IApplication<AppControl>, int, ICellCostEstimator>>();
            clone.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(this.DynamicLoadBalancing_CellCostEstimatorFactories);

            return clone;
        }
    }
}
