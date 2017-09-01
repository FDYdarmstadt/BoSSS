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
using BoSSS.Solution.Multigrid;
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

            C.AddBoundaryCondition("Velocity_Inlet_upper", "VelocityX", X => 0);
            C.AddBoundaryCondition("Velocity_Inlet_lower", "VelocityX", X => 0);
            if (!xPeriodic) {
                C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX", X => (4 * 1.5 * X[1] * (4.1 - X[1]) / (4.1 * 4.1)));
            }
            C.AddBoundaryCondition("Pressure_Outlet_right");

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

            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 20;
            C.MaxSolverIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NoOfMultigridLevels = 1;

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

        static public IBM_Control ChannelFlow(int k = 2, bool periodic = false, bool pardiso = true ,int xCells = 10, int yCells = 10) {
            IBM_Control C = new IBM_Control();

            // Solver Options
            C.NoOfTimesteps = 100;
            C.MaxSolverIterations = 50;
            C.savetodb = false;
            C.ProjectName = "ChannelFlow";

            // Calculate Navier-Stokes? 
            C.PhysicalParameters.IncludeConvection = true;

            // Timestepper
            C.FixedStreamwisePeriodicBC = periodic;
            C.SrcPressureGrad = new double[] { -0.1, 0 };

            if (pardiso) {
                C.whichSolver = DirectSolver._whichSolver.PARDISO;
            } else {
                C.whichSolver = DirectSolver._whichSolver.MUMPS;
            }
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
            double dt = 1E20;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 60;
            C.NoOfTimesteps = 1;

            // Physical values
            C.PhysicalParameters.rho_A = 1;

            // 1/Re
            C.PhysicalParameters.mu_A = 1.0 / 20;


            // Create Fields
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

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 10, 21);
                var _yNodes = GenericBlas.Linspace(-1, 1, 11);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, Foundation.Grid.RefElements.CellType.Square_Linear, periodicX: periodic);

                if (!periodic) {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                }

                grd.EdgeTagNames.Add(2, "Wall_bottom");
                grd.EdgeTagNames.Add(3, "Wall_top");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom
                        return 2;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top
                        return 3;

                    if (!periodic) {
                        if (Math.Abs(x - (0)) < 1.0e-6)
                            // left
                            return 1;

                        if (Math.Abs(x - (10)) < 1.0e-6)
                            // right
                            return 4;
                    }
                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("2D Channel Flow");

                return grd;
            };


            Func<double[], double, double> VelocityX = (X,t) => (1 - (X[1] * X[1]));
            Func<double[], double, double> VelocityY = (X, t) => 0;
            Func<double[], double, double> Pressure = (X, t) => 0;


            C.ExSol_Velocity_Evaluator = new Func<double[], double, double>[] { VelocityX, VelocityY };

            //Func<double[], double, double> Velocity = new  Func<double[], double, double>;
            //C.ExSol_Velocity[0] = (X, t) => (1 - (X[1] * X[1]));
            //C.ExSol_Velocity["A"][1] = (X, t) => (1 - (X[1] * X[1]));
            //C.ExSol_Pressure["A"] = (X, t) => 0;

            if (!periodic) {
                C.AddBoundaryCondition("Velocity_inlet", "VelocityX", X => 1 - X[1] * X[1]);
                C.AddBoundaryCondition("Pressure_Outlet");
            }

            C.AddBoundaryCondition("Wall_bottom");
            C.AddBoundaryCondition("Wall_top");

            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            C.InitialValues_Evaluators.Add("Phi", X => -1);

            return C;
        }


    }
}
