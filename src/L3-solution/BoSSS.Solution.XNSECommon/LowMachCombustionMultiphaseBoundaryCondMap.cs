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
using ilPSP;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Boundary condition mapping for LowMachCombustion XDG multiphase methods.
    /// </summary>
    public class LowMachCombustionMultiphaseBoundaryCondMap : BoSSS.Solution.NSECommon.IncompressibleBoundaryCondMap {

        static string[] BndFunctions(IGridData g, string[] SpeciesNames, int NoOfChemicalSpecies) {
            int D = g.SpatialDimension;
            List<string> scalarFields = new List<string>();

            foreach (var S in SpeciesNames) {
                for (int d = 0; d < D; d++) {
                    scalarFields.Add(VariableNames.Velocity_d(d) + "#" + S);
                }
                scalarFields.Add(VariableNames.Pressure + "#" + S);
                scalarFields.Add(VariableNames.Temperature + '#' + S);
                for (int s = 0; s < NoOfChemicalSpecies; s++) {
                    scalarFields.Add(VariableNames.MassFractions(NoOfChemicalSpecies)[s] + '#' + S);
                }
            }
            scalarFields.Add(VariableNames.LevelSet);

            return scalarFields.ToArray();
        }

        /// <summary>
        /// Loops over all boundary conditions:
        /// If e.g. `VelocityX` is defined, but not `VelocityX#A`, the value for `VelocityX#A` is inferred from `VelocityX`.
        /// </summary>
        static IDictionary<string, AppControl.BoundaryValueCollection> BndyModify(IDictionary<string, AppControl.BoundaryValueCollection> b, string[] SpeciesNames) {

            bool isNonXname(string s) {
                foreach (var sn in SpeciesNames) {
                    string end = "#" + sn;
                    if (s.Length > 2 && s.EndsWith(end))
                        return false;
                }
                return true;
            }

            var ret = new Dictionary<string, AppControl.BoundaryValueCollection>();
            foreach (var kv in b) {
                var coll = kv.Value.CloneAs();
                ret.Add(kv.Key, coll);

                string[] definedKeys = coll.Evaluators.Keys.ToArray();
                foreach (var varName in definedKeys) {
                    if (isNonXname(varName)) {
                        foreach (var spc in SpeciesNames) {
                            string XvarName = varName + "#" + spc;
                            if (!kv.Value.Evaluators.ContainsKey(XvarName)) {
                                coll.Evaluators.Add(XvarName, coll.Evaluators[varName]);
                            }
                        }
                    }
                }
            }
            return ret;
        }


        public LowMachCombustionMultiphaseBoundaryCondMap(IGridData f, IDictionary<string, BoSSS.Solution.Control.AppControl.BoundaryValueCollection> b, string[] SpeciesNames, int NoOfChemicalSpecies)
           : base(f, BndyModify(b, SpeciesNames), PhysicsMode.Combustion, BndFunctions(f, SpeciesNames, NoOfChemicalSpecies)) {
            string S0 = "#" + SpeciesNames[0];

            int D = f.SpatialDimension;
            for (int d = 0; d < D; d++) {
                base.bndFunction.Add(VariableNames.Velocity_d(d), base.bndFunction[VariableNames.Velocity_d(d) + S0]);
            }
            base.bndFunction.Add(VariableNames.Pressure, base.bndFunction[VariableNames.Pressure + S0]);
            base.bndFunction.Add(VariableNames.Temperature, base.bndFunction[VariableNames.Temperature + S0]);

            for (int s = 0; s < NoOfChemicalSpecies; s++) {
                string varname = VariableNames.MassFractions(NoOfChemicalSpecies)[s] ;
                base.bndFunction.Add(varname, base.bndFunction[varname + S0]);
            }
        }
    }


}