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

namespace NSE_SIMPLE {

    /// <summary>
    /// Pressure gradient for momentum equation.
    /// </summary>
    public class MomentumPressureGradientOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="VelocityMapping"></param>
        /// <param name="PressureMapping"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialDirection">
        /// Spatial index for the direction of the gradient.
        /// </param>
        public MomentumPressureGradientOperator(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping PressureMapping,
            SolverConfiguration SolverConf, int SpatialDirection)
            : base(VelocityMapping, PressureMapping, null, SolverConf, true, SpatialDirection: SpatialDirection, OnlyBoundaryEdges: false) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            SpatialOperator PressureOp = new SpatialOperator(new string[] { VariableNames.Pressure }, new string[] { "p1" },QuadOrderFunc.Linear());
            PressureOp.EquationComponents["p1"].Add(new PressureGradientLin_d(SpatialDirection, SolverConf.BcMap));

            if (SolverConf.Control.PressureGradientSource != null) {
                PressureOp.EquationComponents["p1"].Add(
                    new SrcPressureGradientLin_d(SolverConf.Control.PressureGradientSource[SpatialDirection]));
            }

            PressureOp.Commit();
            return PressureOp;
        }
    }

    /// <summary>
    /// Velocity divergence operator.
    /// For incompressible flows used to calculate divergence of velocity in continuity equation.
    /// For low Mach number flows used in the 'analytical compatibility constraint pressure correction'.
    /// </summary>
    public class DivergenceVelocityOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="RowMapping"></param>
        /// <param name="ColMapping"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialIndex"></param>
        public DivergenceVelocityOperator(UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping,
            SolverConfiguration SolverConf, int SpatialIndex)
            : base(RowMapping, ColMapping, null, SolverConf, true, SpatialIndex, SpatialIndex) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            SpatialOperator DivergenceOp = new SpatialOperator(new string[] { VariableNames.Velocity_d(SpatialComponent) }, new string[] { "div_d" }, QuadOrderFunc.Linear());
            DivergenceOp.EquationComponents["div_d"].Add(new Divergence_DerivativeSource(SpatialComponent, SolverConf.SpatialDimension));
            DivergenceOp.EquationComponents["div_d"].Add(new Divergence_DerivativeSource_Flux(SpatialComponent, SolverConf.BcMap));
            DivergenceOp.Commit();
            return DivergenceOp;
        }
    }

    /// <summary>
    /// Pressure stabilization used in continuity equation for equal-order formulations.
    /// </summary>
    public class DivergencePressureStabilization : SIMPLEOperator {

        /// <summary>
        /// Ctor
        /// </summary>        
        /// <param name="PressureMapping"></param>
        /// <param name="SolverConf"></param>
        public DivergencePressureStabilization(UnsetteledCoordinateMapping PressureMapping, SolverConfiguration SolverConf)
            : base(PressureMapping, PressureMapping, null, SolverConf, true) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            return (new PressureStabilization(
                SolverConf.Control.PressureStabilizationScaling,
                base.GridData.Edges.h_max_Edge,
                SolverConf.Control.Reynolds)).Operator();
        }
    }
}
