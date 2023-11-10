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

using ApplicationWithIDT;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using System;
using System.Linq;

namespace XESF {

    /// <summary>
    /// Queries for the XESF solver
    /// </summary>
    public static class XESFQueries {

        /// <summary>
        /// L2 error with respect to given reference solution.
        /// </summary>
        /// <returns></returns>
        public static Query L2Error(string fieldName, string speciesName, Func<double[], double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is XESFMain program)) {
                    throw new NotSupportedException("Only valid for applications of type XESFControl");
                }

                XESFControl control = program.Control;
                LevelSetTracker lsTrk = program.LsTrk;

                var field = program.IOFields.Where(delegate (DGField fld) {
                    if(fld.Identification == fieldName) {
                        return true;
                    } else {
                        return false;
                    }
                }).Single();

                if(!(field is XDGField XDGField)) {
                    throw new NotSupportedException("Only valid for XDGFields");
                }
                XDGField sol_exact = XDGField.CloneAs();
                sol_exact.GetSpeciesShadowField(speciesName).ProjectField(x => referenceSolution(x,time));

                return XDGField.L2Error(sol_exact);
            };
        }
        public static Query LHSInfNorm() {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is ApplicationWithIDT<XESFControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type XESFControl");
                }
                return program.LHS.InfNorm();
            };
        }
        public static Query JrUInfNorm() {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is ApplicationWithIDT<XESFControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type XESFControl");
                }
                return program.Jr_U.InfNorm();
            };
        }
        public static Query JrPhiInfNorm() {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is ApplicationWithIDT<XESFControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type XESFControl");
                }
                return program.Jr_phi.InfNorm();
            };
        }
        public static Query JRUInfNorm() {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is ApplicationWithIDT<XESFControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type XESFControl");
                }
                return program.Jobj_U.InfNorm();
            };
        }
        public static Query JRPhiInfNorm() {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is ApplicationWithIDT<XESFControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type XESFControl");
                }
                return program.Jobj_phi.InfNorm();
            };
        }
    }
}