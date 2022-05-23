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

namespace BoSSS.Application.IBM_Solver {
     public class HardcodedTestExamples {

        public static IBM_Control IBMCylinderFlow(string _DbPath = null, int k = 2, bool xPeriodic = false, double VelXBase = 0.0) {
            IBM_Control C = new IBM_Control();

            const double BaseSize = 1.0;

            // basic database options
            // ====================

            C.savetodb = false;
            C.ProjectDescription = "Cylinder";
            C.Tags.Add("with immersed boundary method");

            // DG degrees
            // ==========

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                var _xNodes1 = Grid1D.TanhSpacing(0, 1, 5, 1, false);
                _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
                var _xNodes2 = GenericBlas.Linspace(1, 5.5, 35);
                _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
                var _xNodes3 = Grid1D.TanhSpacing(5.5, 22, 20, 1.3, true);

                var xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3);


                var _yNodes1 = Grid1D.TanhSpacing(0, 1, 5, 1.2, false);
                _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                var _yNodes2 = GenericBlas.Linspace(1, 3, 20);
                _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
                var _yNodes3 = Grid1D.TanhSpacing(3, 4.1, 5, 1.2, true);

                var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: xPeriodic);
                grd.EdgeTagNames.Add(1, "Velocity_Inlet_upper");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_lower");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - (0 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] + (-4.1 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic && Math.Abs(X[0] - (0 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (!xPeriodic && Math.Abs(X[0] + (-22 * BaseSize)) <= 1.0e-8)
                        et = 4;


                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };

            C.AddBoundaryValue("Velocity_Inlet_upper", "VelocityX", X => 0);
            C.AddBoundaryValue("Velocity_Inlet_lower", "VelocityX", X => 0);
            if (!xPeriodic) {
                C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX", X => (4 * 1.5 * X[1] * (4.1 - X[1]) / (4.1 * 4.1)));
            }
            C.AddBoundaryValue("Pressure_Outlet_right");

            // Initial Values
            // ==============

            C.particleRadius = 0.5;

            C.InitialValues_Evaluators.Add("Phi", X => -(X[0] - 2).Pow2() + -(X[1] - 2).Pow2() + C.particleRadius.Pow2());


            C.InitialValues_Evaluators.Add("VelocityX", X => 0);


            // Physical Parameters
            // ===================

            C.PhysicalParameters.mu_A = 0.05;
            C.PhysicalParameters.rho_A = 1;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            // misc. solver options
            // ====================

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.LevelSetSmoothing = false;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;

            // Timestepping
            // ============
   
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1E20;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 200;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }

    }
}
