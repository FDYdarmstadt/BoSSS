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
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Boundary condition mapping for incompressible XDG multiphase methods.
    /// </summary>
    public class IncompressibleMultiphaseBoundaryCondMap : BoSSS.Solution.NSECommon.IncompressibleBoundaryCondMap {


         static string[] BndFunctions(GridData g, string[] SpeciesNames) {
            int D = g.SpatialDimension;
            List<string> scalarFields = new List<string>();

            foreach(var S in SpeciesNames) {
                for(int d = 0; d < D; d++) {
                    scalarFields.Add(VariableNames.Velocity_d(d) + "#" + S);
                }
                scalarFields.Add(VariableNames.Pressure + "#" + S);
            }

            scalarFields.Add(VariableNames.LevelSet);

            return scalarFields.ToArray();
        }


         public IncompressibleMultiphaseBoundaryCondMap(GridData f, IDictionary<string, BoSSS.Solution.Control.AppControl.BoundaryValueCollection> b, string[] SpeciesNames)
            : base(f, b, PhysicsMode.Incompressible, BndFunctions(f, SpeciesNames)) //
         {
             string S0 = "#" + SpeciesNames[0];

             int D = f.SpatialDimension;
             for(int d = 0; d < D; d++) {
                 base.bndFunction.Add(VariableNames.Velocity_d(d), base.bndFunction[VariableNames.Velocity_d(d) + S0]);
             }
             base.bndFunction.Add(VariableNames.Pressure, base.bndFunction[VariableNames.Pressure + S0]);
        }
    }
}
