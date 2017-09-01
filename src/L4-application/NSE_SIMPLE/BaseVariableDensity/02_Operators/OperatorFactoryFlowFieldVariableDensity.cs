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
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using NSE_SIMPLE.BaseVariableDensity;

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of operators for variable density / variable viscosity flow field, i.e. for momentum and continuity equation.
    /// </summary>
    public abstract class OperatorFactoryFlowFieldVariableDensity {

        /// <summary>
        /// Ctor
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        /// <param name="PressureMapping"></param>
        /// <param name="SolverConf"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>
        /// <param name="Phi0"></param>
        /// <param name="Phi0Mean"></param>
        /// <param name="eta">
        /// Dynamic viscosity for viscous (Laplace) operator.
        /// </param>
        public OperatorFactoryFlowFieldVariableDensity(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping,
            UnsetteledCoordinateMapping PressureMapping,
            SolverConfiguration SolverConf,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SinglePhaseField Phi0, SinglePhaseField Phi0Mean,
            SinglePhaseField eta, SinglePhaseField[] MassFractionsCurrent = null, SinglePhaseField[] MassFractionsMean = null) {

            this.Convection = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.Visc = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.Swip2 = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.Swip3 = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.PressureGradient = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.DivergenceConti = new SIMPLEOperator[SolverConf.SpatialDimension];

            for (int d = 0; d < SolverConf.SpatialDimension; d++) {
                switch (SolverConf.Control.PhysicsMode) {
                    case PhysicsMode.Incompressible:
                        throw new ApplicationException("Wrong OperatorFactory");
                    case PhysicsMode.LowMach:
                    case PhysicsMode.Multiphase:
                        this.Convection[d] = new MomentumConvectionVariableDensityOperator(VelocityMapping, Velocity0, Velocity0Mean, Phi0, Phi0Mean, SolverConf, d);
                        break;
                    default:
                        throw new NotImplementedException("PhysicsMode not implemented");
                }
                this.Visc[d] = new MomentumViscousVariableDensityOperator(VelocityMapping, Phi0, SolverConf, d);
                this.Swip2[d] = new MomentumSwip2(VelocityMapping, VelocityVectorMapping, Phi0, SolverConf, d);
                switch (SolverConf.Control.PhysicsMode) {
                    case PhysicsMode.LowMach:
                        this.Swip3[d] = new MomentumSwip3(VelocityMapping, VelocityVectorMapping, Phi0, SolverConf, d);
                        this.DivergenceConti[d] = GetDivergenceContiOperator(PressureMapping, VelocityMapping, Phi0, SolverConf, d);
                        break;
                    case PhysicsMode.Multiphase:
                        // Divergence of velocity is zero for smooth interface multiphase flows
                        this.Swip3[d] = null;
                        this.DivergenceConti[d] = GetDivergenceContiOperator(PressureMapping, VelocityMapping, Phi0, SolverConf, d);
                        break;
                    default:
                        throw new ApplicationException("PhysicsMode is not compatible with OperatorFactory.");
                }
                this.PressureGradient[d] = new MomentumPressureGradientOperator(VelocityMapping, PressureMapping, SolverConf, d);

            }

            VariableDensitySIMPLEControl varDensConf = SolverConf.Control as VariableDensitySIMPLEControl;
            if (varDensConf.Froude.HasValue) {
                this.BuoyantForce = new SpatialOperator.Evaluator[SolverConf.SpatialDimension];
                for (int d = 0; d < SolverConf.SpatialDimension; d++) {
                    SpatialOperator BuoyancyOperator = (new Buoyancy(varDensConf.GravityDirection,
                        d,
                        varDensConf.Froude.Value,
                        varDensConf.EoS)).Operator();
                    this.BuoyantForce[d] = BuoyancyOperator.GetEvaluator(Phi0.Mapping, VelocityMapping);
                }
            }

            //this.PressureStabilization = new DivergencePressureStabilization(ctx, PressureMapping, SolverConf);
        }

        /// <summary>
        /// Get divergence operator of continuity equation.
        /// </summary>
        /// <param name="PressureMapping"></param>
        /// <param name="VelocityMapping"></param>
        /// <param name="Phi0"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatComp"></param>
        /// <param name="MassFractionsCurrent"></param>
        /// <returns></returns>
        protected abstract SIMPLEOperator GetDivergenceContiOperator(UnsetteledCoordinateMapping PressureMapping, UnsetteledCoordinateMapping VelocityMapping,
            SinglePhaseField Phi0, SolverConfiguration SolverConf, int SpatComp, SinglePhaseField[] MassFractionsCurrent = null);

        /// <summary>
        /// Convection operator of momentum equation.
        /// </summary>
        public SIMPLEOperator[] Convection {
            get;
            private set;
        }

        /// <summary>
        /// Viscous (Laplace) operator of momentum equation.
        /// </summary>
        public SIMPLEOperator[] Visc {
            get;
            private set;
        }

        /// <summary>
        /// Swip2 viscous term.
        /// </summary>
        public SIMPLEOperator[] Swip2 {
            get;
            private set;
        }

        /// <summary>
        /// Swip3 viscous term.
        /// </summary>
        public SIMPLEOperator[] Swip3 {
            get;
            private set;
        }

        /// <summary>
        /// Pressure gradient operator of momentum equation.
        /// </summary>
        public SIMPLEOperator[] PressureGradient {
            get;
            private set;
        }

        /// <summary>
        /// Divergence operator continuity equation.
        /// </summary>
        public SIMPLEOperator[] DivergenceConti {
            get;
            private set;
        }

        /// <summary>
        /// Buoyant force.
        /// </summary>
        public SpatialOperator.Evaluator[] BuoyantForce {
            get;
            private set;
        }

        ///// <summary>
        ///// Pressure stabilization for equal-order formulation.
        ///// </summary>
        //public SIMPLEOperator PressureStabilization {
        //    get;
        //    private set;
        //}
    }
}
