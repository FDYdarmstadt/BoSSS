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

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of operators for incompressible flow field, i.e. for momentum and continuity equation.
    /// </summary>
    public class OperatorFactoryFlowFieldIncompressible {

        /// <summary>
        /// Ctor
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="PressureMapping"></param>
        /// <param name="SolverConf"></param>
        /// <param name="Velocity"></param>
        /// <param name="VelocityMean">
        /// VelocityMean for convective operator.
        /// </param>        
        public OperatorFactoryFlowFieldIncompressible(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping PressureMapping, 
            SolverConfiguration SolverConf, 
            VectorField<SinglePhaseField> Velocity, VectorField<SinglePhaseField> VelocityMean) {

            // Standard operators, which are always used
            this.Convection = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.Visc = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.PressureGradient = new SIMPLEOperator[SolverConf.SpatialDimension];
            this.VelocityDivergence = new SIMPLEOperator[SolverConf.SpatialDimension];            

            for (int d = 0; d < SolverConf.SpatialDimension; d++) {

                this.Convection[d] = new MomentumConvectionOperator(VelocityMapping, Velocity, VelocityMean, SolverConf, d);
                this.Visc[d] = new MomentumViscousOperator(VelocityMapping, SolverConf, d);
                this.PressureGradient[d] = new MomentumPressureGradientOperator(VelocityMapping, PressureMapping, SolverConf, d);                                
                this.VelocityDivergence[d] = new DivergenceVelocityOperator(PressureMapping, VelocityMapping, SolverConf, d);                                                
            }

            // Optional pressure stabilization (needed for equal-order discretization)
            if (SolverConf.Control.PressureStabilizationScaling > 0.0) {
                this.PressureStabilization = new DivergencePressureStabilization(PressureMapping, SolverConf);
            } else {
                this.PressureStabilization = null;
            }

            // Alternative pressure correction replacing M_Div*M_Grad with SIP-discretization
            switch (SolverConf.Control.PredictorApproximation) {
                case PredictorApproximations.Identity:
                case PredictorApproximations.Diagonal:
                case PredictorApproximations.BlockDiagonal:
                    this.IP1PressureCorrection = null;
                    break;
                case PredictorApproximations.Identity_IP1:
                    this.IP1PressureCorrection = new IP1_PressureCorrectionOperator(PressureMapping, SolverConf);
                    break;
                default:
                    throw new NotImplementedException();
            }
        }        

        /// <summary>
        /// Convection operator of momentum equation.
        /// </summary>
        public SIMPLEOperator[] Convection {
            get;
            private set;
        }        

        /// <summary>
        /// Viscous operator of momentum equation.
        /// </summary>
        public SIMPLEOperator[] Visc {
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
        /// Divergence of velocity.
        /// </summary>
        public SIMPLEOperator[] VelocityDivergence {
            get;
            private set;
        }        

        /// <summary>
        /// Pressure stabilization for equal-order formulation.
        /// </summary>
        public SIMPLEOperator PressureStabilization {
            get;
            private set;
        }

        /// <summary>
        /// SIP pressure correction.
        /// </summary>
        public SIMPLEOperator IP1PressureCorrection {
            get;
            private set;
        }
    }
}
