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
using System.Runtime.Serialization;
using System.Linq;
using System.Text;

using ilPSP;
using ilPSP.Utils;

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;


namespace BoSSS.Solution.XheatCommon {


    public class ThermalBoundaryCondMap : BoundaryCondMap<ThermalBcType> {

        static string[] BndFunctions(IGridData g) {
            int D = g.SpatialDimension;

            //return ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Temperature, "HeatFlux");
            return ArrayTools.Cat(VariableNames.HeatFluxVector(D), VariableNames.Temperature, VariableNames.VelocityVector(D));

        }

        /// <summary>
        /// ctor
        /// </summary>
        protected ThermalBoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b, string[] BndFuncName)
            : base(f, b, BndFuncName) {
        }

        /// <summary>
        /// ctor
        /// </summary>
        public ThermalBoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b)
            : base(f, b, BndFunctions(f)) {
        }

    }


    public class ThermalMultiphaseBoundaryCondMap : ThermalBoundaryCondMap {

        static string[] BndFunctions(IGridData g, string[] SpeciesNames) {
            int D = g.SpatialDimension;
            List<string> scalarFields = new List<string>();

            foreach(var S in SpeciesNames) {
                for(int d = 0; d < D; d++) {
                    scalarFields.Add(VariableNames.HeatFluxVectorComponent(d) + "#" + S);
                    scalarFields.Add(VariableNames.Velocity_d(d) + "#" + S);
                }
                scalarFields.Add(VariableNames.Temperature + "#" + S);
                //scalarFields.Add("HeatFlux" + "#" + S);
            }

            return scalarFields.ToArray();
        }


        public ThermalMultiphaseBoundaryCondMap(IGridData f, IDictionary<string, BoSSS.Solution.Control.AppControl.BoundaryValueCollection> b, string[] SpeciesNames)
           : base(f, b, BndFunctions(f, SpeciesNames)) //
        {
            string S0 = "#" + SpeciesNames[0];

            int D = f.SpatialDimension;
            for(int d = 0; d < D; d++) {
                base.bndFunction.Add(VariableNames.HeatFluxVectorComponent(d), base.bndFunction[VariableNames.HeatFluxVectorComponent(d) + S0]);
                base.bndFunction.Add(VariableNames.Velocity_d(d), base.bndFunction[VariableNames.Velocity_d(d) + S0]);
            }
            base.bndFunction.Add(VariableNames.Temperature, base.bndFunction[VariableNames.Temperature + S0]);
            //base.bndFunction.Add("HeatFlux", base.bndFunction["HeatFlux" + S0]);

        }

    }


}
