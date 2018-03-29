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

namespace NSE_SIMPLE {

    /// <summary>
    /// Convection operator for incompressible momentum equation.
    /// </summary>
    public class MomentumConvectionOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="Velocity"></param>
        /// <param name="VelocityMean"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialComponent">
        /// Spatial component index of velocity (i.e. u_i) for affine part of operator.
        /// The operator matrix is the same for all spatial directions.
        /// Therefore, if SpatialComponent!=0, only the affine part of the operator is calculated.
        /// </param>        
        public MomentumConvectionOperator(UnsetteledCoordinateMapping VelocityMapping,
            VectorField<SinglePhaseField> Velocity, VectorField<SinglePhaseField> VelocityMean,
            SolverConfiguration SolverConf, int SpatialComponent)
            : base(VelocityMapping, VelocityMapping, ArrayTools.Cat<SinglePhaseField>(Velocity.ToArray(), VelocityMean.ToArray()),
                SolverConf, false, SpatialComponent, OnlyAffine: (SpatialComponent != 0)) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            return (new LinearizedConvection(SolverConf.SpatialDimension, SolverConf.BcMap, SpatialComponent)).Operator(2);
        }
    }

    /// <summary>
    /// Viscous (Laplace) operator for incompressible momentum equation.
    /// </summary>
    public class MomentumViscousOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor for incompressible flows.
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialComponent">
        /// Spatial component index of velocity (i.e. u_i) for affine part of operator.
        /// The operator matrix is the same for all spatial directions.
        /// Therefore, if SpatialComponent!=0, only the affine part of the operator is calculated.
        /// </param>
        public MomentumViscousOperator(UnsetteledCoordinateMapping VelocityMapping, SolverConfiguration SolverConf, int SpatialComponent)
            : base(VelocityMapping, VelocityMapping, null, SolverConf, true, SpatialComponent, OnlyAffine: (SpatialComponent != 0)) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            //return (new Viscosity(SolverConf.PenaltyViscMomentum, SolverConf.reynolds, SolverConf.BcMap, SpatialComponent, false)).Operator();            
            return (new swipViscosity_Term1_variante(
                SolverConf.PenaltyViscMomentum,
                this.GridData.Cells.cj,
                SpatialComponent,
                SolverConf.SpatialDimension,
                SolverConf.BcMap,
                ViscosityImplementation.H,
                ViscosityOption.ConstantViscosityDimensionless,
                reynolds: SolverConf.Control.Reynolds)).Operator();
        }
    }    
}
