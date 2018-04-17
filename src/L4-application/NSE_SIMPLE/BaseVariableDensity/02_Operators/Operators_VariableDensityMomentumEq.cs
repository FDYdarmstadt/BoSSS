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
using ilPSP.Utils;
using NSE_SIMPLE.BaseVariableDensity;

namespace NSE_SIMPLE {

    /// <summary>
    /// Convection operator for variable density momentum equation.
    /// </summary>
    public class MomentumConvectionVariableDensityOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor for low Mach number flows.
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>
        /// <param name="Phi0"></param>
        /// <param name="Phi0Mean"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialComponent">
        /// Spatial component index of velocity (i.e. u_i) for affine part of operator.
        /// The operator matrix is the same for all spatial directions.
        /// Therefore, if SpatialComponent!=0, only the affine part of the operator is calculated.
        /// </param>    
        public MomentumConvectionVariableDensityOperator(UnsetteledCoordinateMapping VelocityMapping,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SinglePhaseField Phi0, SinglePhaseField Phi0Mean, SolverConfiguration SolverConf, int SpatialComponent)
            : base(VelocityMapping, VelocityMapping, ArrayTools.Cat<SinglePhaseField>(Velocity0.ToArray(), Velocity0Mean.ToArray(), Phi0, Phi0Mean),
                SolverConf, false, SpatialComponent,
                OnlyAffine: (SpatialComponent != 0), OnlyBoundaryEdges: true, MaxUsePerIterMatrix: SolverConf.SpatialDimension) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            VariableDensitySIMPLEControl varDensConf = SolverConf.Control as VariableDensitySIMPLEControl;

            //To be used with OnlyBoundaryEdges: false
            //return (new CoupledLaxFriedrichsVelocity(SpatialComponent, SolverConf.SpatialDimension, SolverConf.EoS, SolverConf.BcMap)).Operator();
            //To be used with OnlyBoundaryEdges: true
            return (new LinearizedConvection(SolverConf.SpatialDimension, SolverConf.BcMap, SpatialComponent, varDensConf.EoS)).Operator(2);
        }
    }

    /// <summary>
    /// Laplace operator for variable density / variable viscosity flows.
    /// </summary>
    public class MomentumViscousVariableDensityOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor for low Mach number flows.
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="Phi"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialComponent">
        /// Spatial component index of velocity (i.e. u_i) for affine part of operator.
        /// The operator matrix is the same for all spatial directions.
        /// Therefore, if SpatialComponent!=0, only the affine part of the operator is calculated.
        /// </param>
        public MomentumViscousVariableDensityOperator(UnsetteledCoordinateMapping VelocityMapping, SinglePhaseField Phi,
            SolverConfiguration SolverConf, int SpatialComponent)
            : base(VelocityMapping, VelocityMapping, new SinglePhaseField[] { Phi }, SolverConf, false, SpatialComponent,
                OnlyAffine: (SpatialComponent != 0), MaxUsePerIterMatrix: SolverConf.SpatialDimension) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            VariableDensitySIMPLEControl varDensConf = SolverConf.Control as VariableDensitySIMPLEControl;

            //return (new Viscosity(SolverConf.PenaltyViscMomentum, SolverConf.reynolds, SolverConf.BcMap, SpatialComponent, SolverConf.ConfigVariableDensity.EoS)).Operator();
            return (new swipViscosity_Term1_variante(SolverConf.PenaltyViscMomentum,
                base.GridData.Cells.cj,
                SpatialComponent,
                SolverConf.SpatialDimension,
                SolverConf.BcMap,
                ViscosityImplementation.H,
                ViscosityOption.VariableViscosityDimensionless,
                reynolds: varDensConf.Reynolds,
                EoS: varDensConf.EoS)).Operator(2);
        }
    }

    /// <summary>
    /// Swip2 viscosity term for low Mach number and smooth interface.
    /// </summary>
    public class MomentumSwip2 : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        /// <param name="Phi"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatDirection"></param>
        public MomentumSwip2(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping, SinglePhaseField Phi,
            SolverConfiguration SolverConf, int SpatDirection)
            : base(VelocityMapping, VelocityVectorMapping, new SinglePhaseField[] { Phi }, SolverConf, false, SpatialDirection: SpatDirection,
            MaxUsePerIterMatrix: SolverConf.SpatialDimension, MaxUsePerIterAffine: SolverConf.SpatialDimension) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            VariableDensitySIMPLEControl varDensConf = SolverConf.Control as VariableDensitySIMPLEControl;

            return (new swipViscosity_Term2(SolverConf.PenaltyViscMomentum,
                base.GridData.Cells.cj,
                SpatialDirection,
                SolverConf.SpatialDimension,
                SolverConf.BcMap,
                ViscosityImplementation.H,
                ViscosityOption.VariableViscosityDimensionless,
                ViscositySolverMode.Segregated,
                reynolds: varDensConf.Reynolds,
                EoS: varDensConf.EoS)).Operator(2);
        }
    }

    /// <summary>
    /// Swip3 viscosity term for low Mach number flows.
    /// </summary>
    public class MomentumSwip3 : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        /// <param name="Phi"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatDirection"></param>
        public MomentumSwip3(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping, SinglePhaseField Phi, SolverConfiguration SolverConf, int SpatDirection)
            : base(VelocityMapping, VelocityVectorMapping, new SinglePhaseField[] { Phi }, SolverConf, false, SpatialDirection: SpatDirection,
            MaxUsePerIterMatrix: SolverConf.SpatialDimension, MaxUsePerIterAffine: SolverConf.SpatialDimension) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            VariableDensitySIMPLEControl varDensConf = SolverConf.Control as VariableDensitySIMPLEControl;

            return (new swipViscosity_Term3(SolverConf.PenaltyViscMomentum,
                base.GridData.Cells.cj,
                SpatialDirection,
                SolverConf.SpatialDimension,
                SolverConf.BcMap,
                ViscosityImplementation.H,
                ViscosityOption.VariableViscosityDimensionless,
                ViscositySolverMode.Segregated,
                reynolds: SolverConf.Control.Reynolds,
                EoS: varDensConf.EoS)).Operator(2);
        }
    }

}
