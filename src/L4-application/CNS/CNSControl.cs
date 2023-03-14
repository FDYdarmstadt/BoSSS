﻿/* =======================================================================
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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Control;
using BoSSS.Solution.LoadBalancing;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.LoadBalancing;
using CNS.ShockCapturing;
using CNS.Source;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CNS {

    /// <summary>
    /// Specialized control file for CNS
    /// </summary>
    [Serializable]
    public class CNSControl : CompressibleControl {

        /// <summary>
        /// Verifies the configuration
        /// </summary>
        public override void Verify() {
            {
                int degree = FieldOptions[BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Momentum.xComponent].Degree;

                if (FieldOptions.ContainsKey(BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Momentum.yComponent)
                    && FieldOptions[BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Momentum.yComponent].Degree != degree) {
                    throw new Exception(
                        "All momentum components must have the same polynomial degree");
                }

                if (FieldOptions.ContainsKey(BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Momentum.zComponent)
                    && FieldOptions[BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Momentum.zComponent].Degree != degree) {
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
        /// Or-combination of all variables for which initial value definitions
        /// exist
        /// </summary>
        public VariableTypes GetInitialValueVariables() {
            bool conservative = InitialValues_Evaluators.ContainsKey(BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Density)
                && InitialValues_Evaluators.ContainsKey(BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Energy);
            for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                conservative &= InitialValues_Evaluators.ContainsKey(BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Momentum[d]);
            }

            if (conservative) {
                return VariableTypes.ConservativeVariables;
            }

            bool primitive = InitialValues_Evaluators.ContainsKey(BoSSS.Solution.CompressibleFlowCommon.CompressibleVariables.Density)
                && InitialValues_Evaluators.ContainsKey(CNSVariables.Pressure);
            for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                primitive &= InitialValues_Evaluators.ContainsKey(CNSVariables.Velocity[d]);
            }

            if (primitive) {
                return VariableTypes.PrimitiveVariables;
            }

            return VariableTypes.None;
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
        /// Forces the local time stepping algorithm to recluster in the
        /// specified <see cref="ReclusteringInterval"/>
        /// </summary>
        public bool forceReclustering = false;

        /// <summary>
        /// Enables/Disables the flux correction to obtain a (non-)conservative
        /// <see cref="ExplicitSchemes.LTS"/> scheme. Not avaiblable for adaptive LTS.
        /// </summary>
        /// <remarks>
        /// This option is only used if <see cref="ExplicitScheme"/> is equal to LTS.
        /// </remarks>
        public bool FluxCorrection = true;

        /// <summary>
        /// Required for the SIPG method. Should be > 1
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double SIPGPenaltyScaling = 1.0;

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
        /// Configuration for an option sponge layer defining a non-reflecting
        /// boundary condition.
        /// </summary>
        public SpongeLayerConfig SpongeLayerConfig = null;

        /*
        /// <summary>
        /// A classifier that decides which performance class a cell
        /// (currently) belongs to
        /// </summary>
        public ICellClassifier DynamicLoadBalancing_CellClassifier = new IndifferentCellClassifier();
        */

        /// <summary>
        /// The maximum number of sub-steps in the smallest cluster in a LTS clustering
        /// </summary>
        [InclusiveLowerBound(0)]
        public int maxNumOfSubSteps = 0;

        /// <summary>
        /// Writes additional LTS information (number of sub-steps, dt per cluster, etc.) to a text file
        /// </summary>
        public bool WriteLTSLog = false;

        /// <summary>
        /// Enable console output for LTS time stepper
        /// </summary>
        public bool WriteLTSConsoleOutput = false;

        /// <summary>
        /// An optional sensor to detect shocks
        /// </summary>
        public ICNSShockSensor CNSShockSensor = null;

        /// <summary>
        /// An optional viscosity law to determine the magnitude of the
        /// artificial viscosity if <see cref="ActiveOperators"/> includes
        /// <see cref="Operators.ArtificialViscosity"/>
        /// </summary>
        public IArtificialViscosityLaw ArtificialViscosityLaw = null;

        /// <summary>
        /// Clones this object, but beware: I'm not sure (yet) that I've
        /// covered all cases
        /// </summary>
        /// <returns></returns>
        public override object Clone() {
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

            clone.DynamicLoadBalancing_CellCostEstimators = new List<ICellCostEstimator>();
            clone.DynamicLoadBalancing_CellCostEstimators.AddRange(this.DynamicLoadBalancing_CellCostEstimators.Select(e => e.CloneAs()));
            
            
            return clone;
        }

        /// <summary>
        /// To launch CNS.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(CNS.CNSProgram);
        }
    }
}
