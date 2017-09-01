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
    /// Convection operator for scalar equation.
    /// </summary>
    public class LevelSetAdvectionOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="LevelSetMapping"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>        
        /// <param name="SolverConf"></param>
        public LevelSetAdvectionOperator(UnsetteledCoordinateMapping LevelSetMapping,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SolverConfiguration SolverConf)
            : base(LevelSetMapping, LevelSetMapping, ArrayTools.Cat<SinglePhaseField>(Velocity0.ToArray(), Velocity0Mean.ToArray()),
                SolverConf, IsConstant: false) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            return (new LinearizedScalarConvection(SolverConf.SpatialDimension, SolverConf.BcMap, null)).Operator();
            //return (new CoupledLaxFriedrichsScalar(SolverConf.SpatialDimension, SolverConf.EoS, SolverConf.BcMap)).Operator();
            //return (new LinearizedScalarUpwind(SolverConf.SpatialDimension, SolverConf.BcMap, SolverConf.EoS)).Operator();
        }
    }
}