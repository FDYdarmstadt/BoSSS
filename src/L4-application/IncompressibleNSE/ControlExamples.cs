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
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Utils;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Application.IncompressibleNSE {
    
    /// <summary>
    /// Some (actally, only one) build-in example configuration
    /// </summary>
    static public class ControlExamples {

        static public IncompressibleControl ChannelFlow(int k = 2, int GridRes = 10, bool transient = false) {
            //BoSSS.Application.IncompressibleNSE.ControlExamples.ChannelFlow(k:2,GridRes:4)
            IncompressibleControl C = new IncompressibleControl();

            // database and saving options
            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";



            // Solver Options
            if(transient) {
                C.NoOfTimesteps = 10;
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
                C.dtFixed = 0.01;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            }

            // Re
            C.Density = 1.0;
            C.Viscosity = 1.0/10.0;

            // degree
            C.SetDGdegree(k);

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 10, GridRes * 5 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, GridRes + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, Foundation.Grid.RefElements.CellType.Square_Linear);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    double x = _X[0];
                    double y = _X[1];

                    if(Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom
                        return "Wall_bottom";

                    if(Math.Abs(y - (+1)) < 1.0e-6)
                        // top
                        return "Wall_top";


                    if(Math.Abs(x - (0.0)) < 1.0e-6)
                        // left
                        return "Velocity_inlet";

                    if(Math.Abs(x - (10)) < 1.0e-6)
                        // right
                        return "Pressure_Outlet";

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // set boundary conditions
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", "X => 1 - X[1] * X[1]", TimeDependent: false);
            //C.AddBoundaryValue("Pressure_Outlet");
            //C.AddBoundaryValue("Wall_bottom");
            //C.AddBoundaryValue("Wall_top"); // optional; if not specified, 0.0 is assumed

            // Set Initial Conditions
            C.AddInitialValue("VelocityX", "X => 0.0", TimeDependent: false);

            return C;
        }


    }
}
