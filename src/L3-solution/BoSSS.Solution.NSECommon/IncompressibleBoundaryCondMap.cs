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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Default boundary condition map that should be used for all kinds of incompressible Navier-Stokes solvers.
    /// </summary>
    public class IncompressibleBoundaryCondMap : BoundaryCondMap<IncompressibleBcType> {

        static string[] BndFunctions(IGridData g, PhysicsMode _PhysicsMode) {
            int D = g.SpatialDimension;
            string[] ScalarFields;

            switch (_PhysicsMode) {
                case PhysicsMode.Incompressible:
                    ScalarFields = new string[] { VariableNames.Pressure };
                    break;
                case PhysicsMode.LowMach:
                    ScalarFields = new string[] { VariableNames.Pressure, VariableNames.Temperature };
                    break;
                case PhysicsMode.Multiphase:
                    ScalarFields = new string[] { VariableNames.Pressure, VariableNames.LevelSet };
                    break;
                case PhysicsMode.Combustion:
                    ScalarFields = new string[] { VariableNames.Pressure, VariableNames.Temperature, VariableNames.MassFraction0, VariableNames.MassFraction1, VariableNames.MassFraction2, VariableNames.MassFraction3 };
                    break;

                case PhysicsMode.Helical:
                    ScalarFields = new string[] { VariableNames.u, VariableNames.v, VariableNames.w, VariableNames.Pressure };
                    return ScalarFields;

                case PhysicsMode.Viscoelastic:
                    ScalarFields = new string[] { VariableNames.Pressure, VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY };
                    break;

                default:
                    throw new ArgumentException();
            }

            return ArrayTools.Cat(VariableNames.VelocityVector(D), ScalarFields);
        }

        public PhysicsMode PhysMode {
            get;
            private set;
        }

        /// <summary>
        /// ctor
        /// </summary>
        protected IncompressibleBoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b, PhysicsMode _PhysicsMode, string[] BndFuncName)
            : base(f, b, BndFuncName) {
            this.PhysMode = _PhysicsMode;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public IncompressibleBoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b, PhysicsMode _PhysicsMode)
            : base(f, b, BndFunctions(f, _PhysicsMode)) {
            this.PhysMode = _PhysicsMode;
        }


        /// <summary>
        /// true if there is at least one pressure-outlet edge:
        /// <see cref="IncompressibleBcType.Pressure_Outlet"/> or <see cref="IncompressibleBcType.Pressure_Dirichlet"/>;
        /// </summary>
        public bool DirichletPressureBoundary {
            get {
                return ((BCTypeUseCount[IncompressibleBcType.Pressure_Outlet] > 0) || (BCTypeUseCount[IncompressibleBcType.Pressure_Dirichlet] > 0));
            }
        }

        /// <summary>
        /// true if there is at least one outflow edge: <see cref="IncompressibleBcType.Outflow"/>;
        /// </summary>
        public bool _OutflowBoundary {
            get {
                return (BCTypeUseCount[IncompressibleBcType.Outflow] > 0);
            }
        }

    }
}
