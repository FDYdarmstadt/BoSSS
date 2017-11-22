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
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.Multigrid;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Solution;
using BoSSS.Platform.Utils.Geom;

namespace BoSSS.Application.IBM_Solver {
    /// <summary>
    /// A few example configurations.
    /// </summary>
    public class HardcodedControl {

        static public IBM_Control ChannelFlow(int k = 2, bool periodic = false, int xCells = 10, int yCells = 10, string dbpath = null) {
            IBM_Control C = new IBM_Control();

            // Solver Options
            C.NoOfTimesteps = 100;
            C.MaxSolverIterations = 100;
            C.MinSolverIterations = 1;
            C.savetodb = false;
            C.DbPath = null;
            C.ProjectName = "ChannelFlow";
            C.SessionName = "GasGebn";

            // Calculate Navier-Stokes? 
            C.PhysicalParameters.IncludeConvection = true;

            // Timestepper
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
            double dt = 1E30;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 60;
            C.NoOfTimesteps = 1;

            // Physical values
            C.PhysicalParameters.rho_A = 1;

            // 1/Re
            C.PhysicalParameters.mu_A = 2.0 / 200;


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
                var _xNodes = GenericBlas.Linspace(0, 10, xCells + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, yCells + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, periodic);

                if (!periodic) {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    //grd.EdgeTagNames.Add(4, "Pressure_Outlet");
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
                            return 1;
                    }
                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("2D Channel Flow");

                return grd;
            };

            Func<double[], double, double> VelocityXex, VelocityYex, Pressure;
            VelocityXex = (X, t) => (1 - (X[1] * X[1]));
            VelocityYex = (X, t) => (0);
            Pressure = (X, t) => (0);



            if (!periodic) {
                C.AddBoundaryCondition("Velocity_inlet", "VelocityX", X => 1 - X[1] * X[1]);
                //C.AddBoundaryCondition("Pressure_Outlet");
            }

            C.AddBoundaryCondition("Wall_bottom");
            C.AddBoundaryCondition("Wall_top");

            // Set Initial Conditions
            //C.InitialValues_Evaluators.Add("VelocityX", X => 1 - X[1] * X[1]);
            //C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues_Evaluators.Add("Pressure", X => 2.0*C.PhysicalParameters.mu_A*(-X[0] + 10));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0.0);
            C.InitialValues_Evaluators.Add("Phi", X => -1);

            return C;
        }

        public static IBM_Control nonIBMCylinderFlow(string _DbPath = null, int k = 2, bool xPeriodic = false, double VelXBase = 0.0) {
            IBM_Control C = new IBM_Control();

            //const double BaseSize = 1.0;

            // basic database options
            // ======================

            C.DbPath = @"\\fdyprime\userspace\krause\BoSSS_DBs\nonIBM_cylinder";
            C.savetodb = true;
            C.ProjectName = "nonIBM/Cylinder/RE100/dt003/DIST";
            C.ProjectDescription = "Cylinder";
            C.Tags.Add("not with immersed boundary method");

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

            C.GridGuid = new Guid("e08e1e84-ceb7-4d5c-a945-d96a5650b66a");


            C.AddBoundaryCondition("Velocity_Inlet", "VelocityX#A", X => 1);

            C.AddBoundaryCondition("Pressure_Outlet");

            C.AddBoundaryCondition("Wall-Upper");
            C.AddBoundaryCondition("Wall-Rear");
            C.AddBoundaryCondition("Wall-Front");
            C.AddBoundaryCondition("Wall-Lower");



            // Initial Values
            // ==============

            //double radius = 0.15;

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            //C.InitialValues.Add("VelocityX", X => 1)
            C.InitialValues_Evaluators.Add("VelocityX", delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                double R = Math.Sqrt((x + 1).Pow2() + y.Pow2());

                double xVel = 1;

                if (R < 0.75) {
                    xVel = 1;
                }
                return xVel;
            });

            C.InitialValues_Evaluators.Add("VelocityY", delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                double R = Math.Sqrt((x + 1).Pow2() + y.Pow2());

                double yVel = 0;

                if (R < 0.75) {
                    yVel = 1;
                }
                return yVel;
            });

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("d789fe46-42ac-4f52-9fd1-bcdb5c23f665"), -1);
            //C.GridGuid = new Guid("0d4abeda-b2ea-40c3-8c3f-345f746569f8");

            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;

            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 20;
            C.MaxSolverIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NoOfMultigridLevels = 0;

            // Timestepping
            // ============

            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
            double dt = 0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = 10000;

            // haben fertig...
            // ===============

            return C;
        }

        public static IBM_Control[] IBMCylinderFlow(string _DbPath = null, int k = 2, bool xPeriodic = false, double VelXBase = 0.0) {
            List<IBM_Control> R = new List<IBM_Control>();

            foreach (int i in new int[] { k }) {

                IBM_Control C = new IBM_Control();

                C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("k", i),
                            };

                k = i;

                const double BaseSize = 1.0;

                // basic database options
                // ======================

                C.DbPath = _DbPath;
                C.savetodb = false;
                C.ProjectName = "FixedCylinderRe100_k" + i + "_CellAgglo02_penalty4_newMesh2";

                switch (i) {
                    case 1:
                        C.MeshFactor = 1.33; // was 1.33
                        break;

                    case 2:
                        C.MeshFactor = 0.92;
                        break;

                    case 3:
                        C.MeshFactor = 0.7; // was 07
                        break;

                    default:

                        throw new ApplicationException();
                }

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
                //Console.WriteLine("Achtung: equal order!!!!");
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

                //grid and boundary conditions
                // ============================

                C.GridFunc = delegate {

                    var _xNodes1 = Grid1D.TanhSpacing(-2, -1, Convert.ToInt32(10 * C.MeshFactor), 0.5, false); //10
                    _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
                    var _xNodes2 = GenericBlas.Linspace(-1, 2, Convert.ToInt32(35 * C.MeshFactor)); //35
                    _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
                    var _xNodes3 = Grid1D.TanhSpacing(2, 20, Convert.ToInt32(60 * C.MeshFactor), 1.5, true); //60

                    var xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3);


                    var _yNodes1 = Grid1D.TanhSpacing(-2, -1, Convert.ToInt32(7 * C.MeshFactor), 0.9, false); //7
                    _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                    var _yNodes2 = GenericBlas.Linspace(-1, 1, Convert.ToInt32(25 * C.MeshFactor)); //25
                    _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
                    var _yNodes3 = Grid1D.TanhSpacing(1, 2.1, Convert.ToInt32(7 * C.MeshFactor), 1.1, true); //7
                    var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);



                    //double[] xNodes = GenericBlas.Linspace(0 * BaseSize, 22 * BaseSize, 25);
                    //double[] yNodes = GenericBlas.Linspace(0 * BaseSize, 4.1 * BaseSize, 25);
                    var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: xPeriodic);
                    grd.EdgeTagNames.Add(1, "Velocity_Inlet_upper");
                    grd.EdgeTagNames.Add(2, "Velocity_Inlet_lower");
                    if (!xPeriodic) {
                        grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                        grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                    }

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1] - (-2 * BaseSize)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - (+2.1 * BaseSize)) <= 1.0e-8)
                            et = 2;
                        if (!xPeriodic && Math.Abs(X[0] - (-2 * BaseSize)) <= 1.0e-8)
                            et = 3;
                        if (!xPeriodic && Math.Abs(X[0] - (+20.0 * BaseSize)) <= 1.0e-8)
                            et = 4;


                        Debug.Assert(et != 0);
                        return et;
                    });

                    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                    return grd;
                };

                //C.GridFunc = delegate {

                //    // Box1
                //    var box1_p1 = new double[2] { -2, -2 };
                //    var box1_p2 = new double[2] { 20, 2.1 };
                //    var box1 = new GridBox(box1_p1, box1_p2, 46, 20); //k1: 70,25 ; k2: 46,20 ; k3: 35,15

                //    // Box2
                //    var box2_p1 = new double[2] { -2, -2 };
                //    var box2_p2 = new double[2] { 3, 2.1 };
                //    var box2 = new GridBox(box2_p1, box2_p2, 26, 40); //k1: 40,50 ; k2: 26,40; k3: 20, 30

                //    // Box3
                //    var box3_p1 = new double[2] { -2, -1 };
                //    var box3_p2 = new double[2] { 1, 1 };
                //    var box3 = new GridBox(box3_p1, box3_p2, 32, 38); //k1: 48,58  ; k2: 32,38; k3: 24, 30

                //    // Box4
                //    var box4_p1 = new double[2] { -0.7, -0.72 };
                //    var box4_p2 = new double[2] { 0.7, 0.7 };
                //    var box4 = new GridBox(box4_p1, box4_p2, 30, 56); //k1: 44,84  ; k2: 30,56; k3: 22, 42

                //    var grd = Grid2D.HangingNodes2D(box1, box2, box3,box4);

                //    grd.EdgeTagNames.Add(1, "Velocity_Inlet_upper");
                //    grd.EdgeTagNames.Add(2, "Velocity_Inlet_lower");
                //    if (!xPeriodic) {
                //        grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                //        grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                //    }

                //    grd.DefineEdgeTags(delegate (double[] X) {
                //        byte et = 0;
                //        if (Math.Abs(X[1] - (-2 * BaseSize)) <= 1.0e-8)
                //            et = 1;
                //        if (Math.Abs(X[1] - (+2.1 * BaseSize)) <= 1.0e-8)
                //            et = 2;
                //        if (!xPeriodic && Math.Abs(X[0] - (-2 * BaseSize)) <= 1.0e-8)
                //            et = 3;
                //        if (!xPeriodic && Math.Abs(X[0] - (+20.0 * BaseSize)) <= 1.0e-8)
                //            et = 4;


                //        Debug.Assert(et != 0);
                //        return et;
                //    });

                //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                //    return grd;
                //};

                C.AddBoundaryCondition("Velocity_Inlet_upper", "VelocityX", X => 0);
                C.AddBoundaryCondition("Velocity_Inlet_lower", "VelocityX", X => 0); //-(4 * 1.5 * X[1] * (4.1 - X[1]) / (4.1 * 4.1))
                if (!xPeriodic) {
                    C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX", X => (4 * 1.5 * (X[1] + 2) * (4.1 - (X[1] + 2)) / (4.1 * 4.1)));
                    //C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX#A", X => 1);   
                }
                C.AddBoundaryCondition("Pressure_Outlet_right");


                // Initial Values
                // ==============

                double radius = 0.5;
                C.PhysicalParameters.rho_A = 1;
                C.PhysicalParameters.mu_A = 1.0 / 20;

                //C.InitialValues.Add("Phi", X => phi(X, 0));

                //C.InitialValues.Add("Phi", X => ((X[0] / (radius * BaseSize)) - mPx) * (X[0] / (radius * BaseSize)) - mPx) + ((X[1]) / (radius * BaseSize)) - 2.)Pow2() - radius.Pow2()));  // quadratic form
                //    );
                C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + radius.Pow2());
                //C.InitialValues.Add("Phi", X => -1);

                C.InitialValues_Evaluators.Add("VelocityX", X => 4 * 1.5 * (X[1] + 2) * (4.1 - (X[1] + 2)) / (4.1 * 4.1));
                //C.InitialValues.Add("VelocityX", delegate (double[] X)
                //{
                //    double x = X[0];
                //    double y = X[1];

                //    double R = Math.Sqrt((x + 1).Pow2() + y.Pow2());

                //    double xVel = 0;

                //    if (R < 0.75)
                //    {
                //        xVel = 1;
                //    }
                //    return xVel;
                //});

                //C.InitialValues.Add("VelocityY", delegate (double[] X) {
                //    double x = X[0];
                //    double y = X[1];

                //    double R = Math.Sqrt((x + 1).Pow2() + (y).Pow2());

                //    double yVel = 0;

                //    if (R < 0.75) {
                //        yVel = 1;
                //    }
                //    return yVel;
                //});

                // For restart
                //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("8f5cfed9-31c7-4be8-aa56-e92e5348e08b"), 95);
                //C.GridGuid = new Guid("71ffc0c4-66aa-4762-b07e-45385f34b03f");

                // Physical Parameters
                // ===================


                C.PhysicalParameters.IncludeConvection = true;


                // misc. solver options
                // ====================

                C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
                C.AdvancedDiscretizationOptions.PenaltySafety = 4;
                C.LevelSetSmoothing = false;
                //C.option_solver = "direct";
                C.MaxKrylovDim = 20;
                C.MaxSolverIterations = 50;
                C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
                C.NoOfMultigridLevels = 0;

                // Timestepping
                // ============

                C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
                double dt = 1E14;
                C.dtMax = dt;
                C.dtMin = dt;
                C.Endtime = 70;
                C.NoOfTimesteps = 1;

                // haben fertig...
                // ===============

                R.Add(C);
            }

            return R.ToArray();
        }

        /// <summary>
        /// Supporting method for the 2D hill flow which describes the hill topology
        /// </summary>
        /// <param name="rPoint"></param>
        /// <param name="sPoint"></param>
        /// <returns></returns>
        public static double[] HillTopology(double rPoint, double sPoint, double h) {

            // Calculate Coordinates
            double[] xyPoint = new double[] { 0, 0 };

            double x = rPoint * 9 * h;
            double y;
            double fx;

            //// HILL TOPOLOGY, Input x, Output y

            // HILL INFLOW
            if (x <= (0.3214 * h)) {
                y = 1 + 2.420E-4 * x * x - 7.588E-5 * x * x * x;
                if (1 < y) {
                    fx = 1;
                } else {
                    fx = y;
                }
            } else if (x <= (0.5 * h)) {
                fx = 0.8955 + 3.484E-2 * x - 3.629E-3 * x * x + 6.749E-5 * x * x * x;
            } else if (x <= (0.7143 * h)) {
                fx = 0.9213 + 2.931E-2 * x - 3.234E-3 * x * x + 5.809E-5 * x * x * x;
            } else if (x <= (1.071 * h)) {
                fx = 1.445 - 4.927E-2 * x + 6.950E-4 * x * x - 7.394E-6 * x * x * x;
            } else if (x <= (1.429 * h)) {
                fx = 0.6401 + 3.123E-2 * x - 1.988E-3 * x * x + 2.242E-5 * x * x * x;
            } else if (x <= (1.929 * h)) {
                y = 2.0139 - 7.180E-2 * x + 5.875E-4 * x * x + 9.553E-7 * x * x * x;
                if (0 > y) {
                    fx = 0;
                } else {
                    fx = y;
                }
            }

            // HILL OUTFLOW
            else if (x >= (9 * h - 0.3214 * h)) {
                double x1 = (9 * h) - x;
                y = 1 + 2.420E-4 * x1 * x1 - 7.588E-5 * x1 * x1 * x1;
                if (1 < y) {
                    fx = 1;
                } else {
                    fx = y;
                }
            } else if (x >= (9 * h - 0.5 * h)) {
                double x1 = (9 * h) - x;
                fx = 0.8955 + 3.484E-2 * x1 - 3.629E-3 * x1 * x1 + 6.749E-5 * x1 * x1 * x1;
            } else if (x >= (9 * h - 0.7143 * h)) {
                double x1 = (9 * h) - x;
                fx = 0.9213 + 2.931E-2 * x1 - 3.234E-3 * x1 * x1 + 5.809E-5 * x1 * x1 * x1;
            } else if (x >= (9 * h - 1.071 * h)) {
                double x1 = (9 * h) - x;
                fx = 1.445 - 4.927E-2 * x1 + 6.950E-4 * x1 * x1 - 7.394E-6 * x1 * x1 * x1;
            } else if (x >= (9 * h - 1.429 * h)) {
                double x1 = (9 * h) - x;
                fx = 0.6401 + 3.123E-2 * x1 - 1.988E-3 * x1 * x1 + 2.242E-5 * x1 * x1 * x1;
            } else if (x >= (9 * h - 1.929 * h)) {
                double x1 = (9 * h) - x;
                y = 2.0139 - 7.180E-2 * x1 + 5.875E-4 * x1 * x1 + 9.553E-7 * x1 * x1 * x1;
                if (0 > y) {
                    fx = 0;
                } else {
                    fx = y;
                }
            } else {
                fx = 0;
            }

            if (fx > 1) {
                Console.WriteLine(fx);
            };

            y = (sPoint * 3.035 + (1 - sPoint) * fx);

            xyPoint[0] = x / h; // x/h
            xyPoint[1] = y; // y/h

            if (xyPoint[1] < 0) {
                Console.WriteLine(xyPoint[1]);
            };

            return xyPoint;
        }


        static public IBM_Control PeriodicHill(int k = 2, double HillHeight = 28) {
            IBM_Control C = new IBM_Control();

            if (HillHeight != 28) {
                throw new NotImplementedException("The Testcase is currently only implemented for a hill height of 28. This is due to the fact that the topology-polynomials are only given in an dimensionless manner for a height of 28.");
            }
            // Solver Options
            C.NoOfTimesteps = 100;
            C.MaxSolverIterations = 50;
            C.savetodb = false;
            C.DbPath = @"P:\BoSSS_DBs\ChannelFlow";
            C.ProjectName = "ChannelFlow";

            C.PhysicalParameters.IncludeConvection = true;

            // Physical Paramter
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1.0 / 200;

            // Setting the timestep size according to a CFL condition (trying to at least)
            // this is a rough estimate (cf. Bjoern Mueller's Dissertation, eq. 4.19)
            // here, twice the average velocity is used to estimate the CFL number
            //double CFLFactor = 2 * C.Rey * (2 * C.DegreeVelocity + 1) * 41 / (9 * HillHeight * HillHeight);
            //double CFLDesired = 0.005;
            // CFL = CFLFactor * dt  <==>  dt = CFL / CFLFactor

            // Setting the timestep
            //C.dtMax = CFLDesired / CFLFactor;
            //C.dtMin = C.dtMax;

            C.ProjectDescription = "dt = " + C.dtMax;


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
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // Create Grid
            // Note: The domain is coordinates are dimensionless
            C.GridFunc = delegate {


                var _rNodes = GenericBlas.Linspace(0, 1, 41);
                var _sNodes = GenericBlas.Linspace(0, 1, 16);

                var grd = Grid2D.CurvedSquareGridChannel(_rNodes, _sNodes, CellType.Square_9, true, (r, s) => HillTopology(r, s, HillHeight));

                grd.EdgeTagNames.Add(2, "Wall_bottom");
                grd.EdgeTagNames.Add(3, "Wall_top");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;

                    #region "Calculation of r and s coordinates out of the physical coordinates ("Reversed" Hill Topology)
                    //---------------------------------------------------------------------
                    // HILL TOPOLOGY, Input x, Output y
                    double x = X[0] * HillHeight;
                    double y = X[1];
                    double r, s;
                    double fx;

                    // HILL INFLOW
                    if (x <= (0.3214 * HillHeight)) {
                        y = 1 + 2.420E-4 * x * x - 7.588E-5 * x * x * x;
                        if (1 < y) {
                            fx = 1;
                        } else {
                            fx = y;
                        }
                    } else if (x <= (0.5 * HillHeight)) {
                        fx = 0.8955 + 3.484E-2 * x - 3.629E-3 * x * x + 6.749E-5 * x * x * x;
                    } else if (x <= (0.7143 * HillHeight)) {
                        fx = 0.9213 + 2.931E-2 * x - 3.234E-3 * x * x + 5.809E-5 * x * x * x;
                    } else if (x <= (1.071 * HillHeight)) {
                        fx = 1.445 - 4.927E-2 * x + 6.950E-4 * x * x - 7.394E-6 * x * x * x;
                    } else if (x <= (1.429 * HillHeight)) {
                        fx = 0.6401 + 3.123E-2 * x - 1.988E-3 * x * x + 2.242E-5 * x * x * x;
                    } else if (x <= (1.929 * HillHeight)) {
                        y = 2.0139 - 7.180E-2 * x + 5.875E-4 * x * x + 9.553E-7 * x * x * x;
                        if (0 > y) {
                            fx = 0;
                        } else {
                            fx = y;
                        }
                    }

                    // HILL OUTFLOW
                    else if (x >= (9 * HillHeight - 0.3214 * HillHeight)) {
                        double x1 = (9 * HillHeight) - x;
                        y = 1 + 2.420E-4 * x1 * x1 - 7.588E-5 * x1 * x1 * x1;
                        if (1 < y) {
                            fx = 1;
                        } else {
                            fx = y;
                        }
                    } else if (x >= (9 * HillHeight - 0.5 * HillHeight)) {
                        double x1 = (9 * HillHeight) - x;
                        fx = 0.8955 + 3.484E-2 * x1 - 3.629E-3 * x1 * x1 + 6.749E-5 * x1 * x1 * x1;
                    } else if (x >= (9 * HillHeight - 0.7143 * HillHeight)) {
                        double x1 = (9 * HillHeight) - x;
                        fx = 0.9213 + 2.931E-2 * x1 - 3.234E-3 * x1 * x1 + 5.809E-5 * x1 * x1 * x1;
                    } else if (x >= (9 * HillHeight - 1.071 * HillHeight)) {
                        double x1 = (9 * HillHeight) - x;
                        ;
                        fx = 1.445 - 4.927E-2 * x1 + 6.950E-4 * x1 * x1 - 7.394E-6 * x1 * x1 * x1;
                    } else if (x >= (9 * HillHeight - 1.429 * HillHeight)) {
                        double x1 = (9 * HillHeight) - x;
                        fx = 0.6401 + 3.123E-2 * x1 - 1.988E-3 * x1 * x1 + 2.242E-5 * x1 * x1 * x1;
                    } else if (x >= (9 * HillHeight - 1.929 * HillHeight)) {
                        double x1 = (9 * HillHeight) - x;
                        y = 2.0139 - 7.180E-2 * x1 + 5.875E-4 * x1 * x1 + 9.553E-7 * x1 * x1 * x1;
                        if (0 > y) {
                            fx = 0;
                        } else {
                            fx = y;
                        }
                    } else {
                        fx = 0;
                    }

                    r = X[0] / 9;
                    s = (X[1] - fx) / (3.035 - fx);

                    //----------------------------------------------------------------------

                    #endregion

                    if (Math.Abs(s - 0) < 1.0e-6)
                        // bottom
                        return 2;

                    if (Math.Abs(s - 1) < 1.0e-6)
                        // top
                        return 3;

                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("PeriodicHill flow");

                return grd;
            };

            Func<double[], double> VelocityXex, VelocityYex, Pressure;
            VelocityXex = X => (1 - (X[1] * X[1]));
            VelocityYex = X => (0);
            Pressure = X => (0);

            //C.AddBoundaryCondition("Velocity_inlet", "VelocityX", VelocityXex);
            //C.AddBoundaryCondition("Velocity_inlet", "VelocityY", VelocityYex);
            C.AddBoundaryCondition("Wall_bottom");
            C.AddBoundaryCondition("Wall_top");
            //C.AddBoundaryCondition("Pressure_Outlet");

            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues_Evaluators.Add("VelocityX", X => (1/C.PhysicalParameters.mu_A) / HillHeight);
            C.InitialValues_Evaluators.Add("VelocityX", X => 1);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Pressure", X => 0);


            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;

            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 20;
            C.MaxSolverIterations = 20;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NoOfMultigridLevels = 0;

            // Timestepping
            // ============

            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 70;
            C.NoOfTimesteps = 100000;

            return C;
        }

        static public IBM_Control Cylinder3D(int k = 2) {
            IBM_Control C = new IBM_Control();

            // Solver options
            C.MaxSolverIterations = 1;
            C.ProjectName = "Cylinder3D HPCCLUSTER 18/06/15";
            C.NoOfTimesteps = 1;
            C.dtMax = 0.1;
            C.dtMin = 0.1;
            // Create Fields
            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = k - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            // Load Grid
            Console.WriteLine("...loading grid");
            C.GridGuid = new Guid("0e9d9f95-ab58-42a2-b5b1-40db65857372");

            #region Creates grid () and sets BC
            //// Create Grid
            //Console.WriteLine("...generating grid");
            //C.GridFunc = delegate {

            //    // x-direction
            //    var _xNodes1 = Grid1D.ExponentialSpaceing(-9.5, -3, 9, 0.98);
            //    _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
            //    var _xNodes2 = Grid1D.ExponentialSpaceing(-3, -1, 9, 0.95);
            //    _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
            //    var _xNodes3 = Grid1D.ExponentialSpaceing(-1, 0, 9, 1);
            //    _xNodes3 = _xNodes3.GetSubVector(0, (_xNodes3.Length - 1));
            //    var _xNodes4 = Grid1D.ExponentialSpaceing(0, 2, 9, 1.05);
            //    _xNodes4 = _xNodes4.GetSubVector(0, (_xNodes4.Length - 1));
            //    var _xNodes5 = Grid1D.ExponentialSpaceing(2, 8.5, 14, 1.02);
            //    _xNodes5 = _xNodes5.GetSubVector(0, (_xNodes5.Length - 1));
            //    var _xNodes6 = Grid1D.ExponentialSpaceing(8.5, 12.5, 9, 1);

            //    var _xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3, _xNodes4, _xNodes5, _xNodes6);

            //    // y-direction
            //    var _yNodes1 = Grid1D.ExponentialSpaceing(-9, -2.5, 6, 0.91);
            //    _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
            //    var _yNodes2 = Grid1D.ExponentialSpaceing(-2.5, -0.5, 6, 0.95);
            //    _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
            //    var _yNodes3 = Grid1D.ExponentialSpaceing(-0.5, 0.5, 6, 1.0);
            //    _yNodes3 = _yNodes3.GetSubVector(0, (_yNodes3.Length - 1));
            //    var _yNodes4 = Grid1D.ExponentialSpaceing(0.5, 2.5, 6, 1.05);
            //    _yNodes4 = _yNodes4.GetSubVector(0, (_yNodes4.Length - 1));
            //    var _yNodes5 = Grid1D.ExponentialSpaceing(2.5, 9, 6, 1.1);

            //    var _yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3, _yNodes4, _yNodes5);

            //    // z-direction
            //    var _zNodes = GenericBlas.Linspace(-3, 3, 14);

            //    // Cut Out
            //    double[] CutOutPoint1 = new double[3];
            //    CutOutPoint1[0] = -1;
            //    CutOutPoint1[1] = -0.5;
            //    CutOutPoint1[2] = -3;

            //    double[] CutOutPoint2 = new double[3];
            //    CutOutPoint2[0] = 0;
            //    CutOutPoint2[1] = 0.5;
            //    CutOutPoint2[2] = 3;

            //   var CutOut = new BoundingBox(3);
            //    CutOut.AddPoint(CutOutPoint1);
            //    CutOut.AddPoint(CutOutPoint2);

            //    var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, false, true, CellType.Cube_Linear, CutOut);

            //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
            //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");
            //    grd.EdgeTagNames.Add(3, "Wall");

            //    grd.DefineEdgeTags(delegate(double[] _X) {
            //        var X = _X;
            //        double x = X[0];
            //        double y = X[1];
            //        double z = X[2];

            //        if (Math.Abs(x - (-9.5)) < 1.0e-6)
            //            // inlet
            //            return 1;

            //        if (Math.Abs(x - (12.5)) < 1.0e-6)
            //            // outlet
            //            return 2;

            //        if (Math.Abs(y - (-9)) < 1.0e-6)
            //            // left
            //            return 2;

            //        if (Math.Abs(y - (9)) < 1.0e-6)
            //            // right
            //            return 2;

            //        if (Math.Abs(x - (-1)) < 1.0e-6)
            //            // Cube front
            //            return 3;

            //        if (Math.Abs(x - (0)) < 1.0e-6)
            //            // cube back
            //            return 3;

            //        if (Math.Abs(y - (-0.5)) < 1.0e-6)
            //            // cube left
            //            return 3;

            //        if (Math.Abs(y - (0.5)) < 1.0e-6)
            //            // cube right
            //            return 3;

            //        throw new ArgumentOutOfRangeException();
            //    });

            //    return grd;
            //};
            #endregion


            #region Creates grid (17710 Cells) and sets BC
            //// Create Grid
            //Console.WriteLine("...generating grid");
            //C.GridFunc = delegate {

            //    // x-direction
            //    var _xNodes1 = Grid1D.ExponentialSpaceing(-9.5, -3, 11, 0.98);
            //    _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
            //    var _xNodes2 = Grid1D.ExponentialSpaceing(-3, -1, 9, 0.95);
            //    _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
            //    var _xNodes3 = Grid1D.ExponentialSpaceing(-1, 0, 8, 1);
            //    _xNodes3 = _xNodes3.GetSubVector(0, (_xNodes3.Length - 1));
            //    var _xNodes4 = Grid1D.ExponentialSpaceing(0, 2, 9, 1.05);
            //    _xNodes4 = _xNodes4.GetSubVector(0, (_xNodes4.Length - 1));
            //    var _xNodes5 = Grid1D.ExponentialSpaceing(2, 8.5, 16, 1.02);
            //    _xNodes5 = _xNodes5.GetSubVector(0, (_xNodes5.Length - 1));
            //    var _xNodes6 = Grid1D.ExponentialSpaceing(8.5, 12.5, 5, 1);

            //    var _xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3, _xNodes4, _xNodes5, _xNodes6);

            //    // y-direction
            //    var _yNodes1 = Grid1D.ExponentialSpaceing(-9, -2.5, 8, 0.91);
            //    _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
            //    var _yNodes2 = Grid1D.ExponentialSpaceing(-2.5, -0.5, 8, 0.95);
            //    _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
            //    var _yNodes3 = Grid1D.ExponentialSpaceing(-0.5, 0.5, 8, 1.0);
            //    _yNodes3 = _yNodes3.GetSubVector(0, (_yNodes3.Length - 1));
            //    var _yNodes4 = Grid1D.ExponentialSpaceing(0.5, 2.5, 8, 1.05);
            //    _yNodes4 = _yNodes4.GetSubVector(0, (_yNodes4.Length - 1));
            //    var _yNodes5 = Grid1D.ExponentialSpaceing(2.5, 9, 8, 1.1);

            //    var _yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3, _yNodes4, _yNodes5);

            //    // z-direction
            //    var _zNodes = GenericBlas.Linspace(-3, 3, 11);

            //    // Cut Out
            //    double[] CutOutPoint1 = new double[3];
            //    CutOutPoint1[0] = -1;
            //    CutOutPoint1[1] = -0.5;
            //    CutOutPoint1[2] = -3;

            //    double[] CutOutPoint2 = new double[3];
            //    CutOutPoint2[0] = 0;
            //    CutOutPoint2[1] = 0.5;
            //    CutOutPoint2[2] = 3;

            //    var CutOut = new BoundingBox(3);
            //    CutOut.AddPoint(CutOutPoint1);
            //    CutOut.AddPoint(CutOutPoint2);

            //    var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, false, true, CellType.Cube_Linear, CutOut);

            //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
            //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");
            //    grd.EdgeTagNames.Add(3, "Wall");

            //    grd.DefineEdgeTags(delegate(double[] _X) {
            //        var X = _X;
            //        double x = X[0];
            //        double y = X[1];
            //        double z = X[2];

            //        if (Math.Abs(x - (-9.5)) < 1.0e-6)
            //            // inlet
            //            return 1;

            //        if (Math.Abs(x - (12.5)) < 1.0e-6)
            //            // outlet
            //            return 2;

            //        if (Math.Abs(z - (-3)) < 1.0e-6)
            //            // left
            //            return 2;

            //        if (Math.Abs(z - (3)) < 1.0e-6)
            //            // right
            //            return 2;

            //        if (Math.Abs(x - (-1)) < 1.0e-6)
            //            // Cube front
            //            return 3;

            //        if (Math.Abs(x - (0)) < 1.0e-6)
            //            // cube back
            //            return 3;

            //        if (Math.Abs(y - (-0.5)) < 1.0e-6)
            //            // cube left
            //            return 3;

            //        if (Math.Abs(y - (0.5)) < 1.0e-6)
            //            // cube right
            //            return 3;

            //        throw new ArgumentOutOfRangeException();
            //    });

            //    return grd;
            //};
            #endregion

            Console.WriteLine("...starting calculation of Cylinder3D");

            // Initial Solution

            // Boundary conditions
            C.AddBoundaryCondition("Velocity_inlet", "VelocityX", (X, t) => 1);
            C.AddBoundaryCondition("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryCondition("Velocity_inlet", "VelocityZ", (X, t) => 0);
            C.AddBoundaryCondition("Wall");
            C.AddBoundaryCondition("Pressure_Outlet");

            // Set Initial Conditions
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(
            //    new Guid("510aac88-fbb3-467c-bdbd-4d16ebf51b88"), 50);
            C.InitialValues_Evaluators.Add("VelocityX", X => 1);
            C.InitialValues_Evaluators.Add("VelocityY", X => 1);
            C.InitialValues_Evaluators.Add("VelocityZ", X => 1);
            C.InitialValues_Evaluators.Add("Pressure", X => 1);

            return C;
        }

        static public IBM_Control KovasnayFlow(int k = 2, int Cells = 8) {
            IBM_Control C = new IBM_Control();

            // Solver Options
            C.MaxSolverIterations = 1000;
            C.savetodb = false;
            C.DbPath = @"P:\BoSSS_DBs\Kovasznay";
            C.ProjectName = "KovasnayFlow";

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
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // Load Grid
            //C.GridGuid = new Guid("379028b7-d27f-4587-8b18-9b8c8f31a512");

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-0.5, 1.5, Cells + 1);
                var _yNodes = GenericBlas.Linspace(-0.5, 1.5, Cells + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);

                grd.EdgeTagNames.Add(1, "Velocity_inlet_all");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-0.5)) < 1.0e-6)
                        // bottom
                        return 1;

                    if (Math.Abs(y - (+1.5)) < 1.0e-6)
                        // top
                        return 1;

                    if (Math.Abs(x - (-0.5)) < 1.0e-6)
                        // left
                        return 1;

                    if (Math.Abs(x - (+1.5)) < 1.0e-6)
                        // right
                        return 1;
                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("KovasnayFlow");

                return grd;
            };

            //C.Stokes = false;

            //double pressureConst = 0;
            //double lambda = ((C.Rey / 2) - (Math.Sqrt((C.Rey * C.Rey / 4) + 4 * Math.PI * Math.PI)));

            //Func<double[], double, double> VelocityXex, VelocityYex, Pressureex;
            //VelocityXex = (X, t) => (1 - (Math.Exp(lambda * X[0]) * Math.Cos(2 * Math.PI * X[1])));
            //VelocityYex = (X, t) => ((lambda / (2 * Math.PI)) * Math.Exp(lambda * X[0]) * Math.Sin(2 * Math.PI * X[1]));
            //Pressureex = (X, t) => ((-0.5 * Math.Exp(2 * lambda * X[0])) + pressureConst);

            //C.AddBoundaryCondition("Velocity_inlet_all", "VelocityX", VelocityXex);
            // C.AddBoundaryCondition("Velocity_inlet_all", "VelocityY", VelocityYex);

            Func<double[], double> InitialGuess = X => 0;

            C.InitialValues_Evaluators.Add("VelocityX", InitialGuess);
            C.InitialValues_Evaluators.Add("VelocityY", InitialGuess);
            C.InitialValues_Evaluators.Add("Pressure", InitialGuess);

            return C;
        }

        static public IBM_Control DrivenCavity(int k = 2, int Cells = 10, string dbpath = null) {
            IBM_Control C = new IBM_Control();


            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication<AppControl> app, int noOfPerformanceClasses) {
                Console.WriteLine("i was called");
                int[] map = new int[] { 1, 5, 100 };
                return new StaticCellCostEstimator(map);
            });


            // Solver Options
            C.MaxSolverIterations = 100;
            C.MinSolverIterations = 1;
            C.savetodb = false;
            C.DbPath = null;
            C.ProjectName = "ChannelFlow";
            C.SessionName = "Channel";
            C.NoOfMultigridLevels = 1;

            // Calculate Navier-Stokes? 
            C.PhysicalParameters.IncludeConvection = true;

            // Timestepper
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
            double dt = 1E50;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 60;
            C.NoOfTimesteps = 1;

            // Physical values
            C.PhysicalParameters.rho_A = 1;

            // 1/Re
            C.PhysicalParameters.mu_A = 2.0 / 200;


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
                var _xNodes = GenericBlas.Linspace(-1, 1, Cells + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, Cells + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear);


                grd.EdgeTagNames.Add(1, "Velocity_inlet");

                grd.EdgeTagNames.Add(2, "Wall");


                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom
                        return 2;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top
                        return 1;

                    if (Math.Abs(x - (-1)) < 1.0e-6)
                        // left
                        return 2;

                    if (Math.Abs(x - (1)) < 1.0e-6)
                        // right
                        return 2;
                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("DrivenCavity Flow");

                return grd;
            };

            Func<double[], double, double> VelocityXex, VelocityYex, Pressure;
            VelocityXex = (X, t) => (1 - (X[1] * X[1]));
            VelocityYex = (X, t) => (0);
            Pressure = (X, t) => (0);



            C.AddBoundaryCondition("Velocity_inlet", "VelocityX", X => 1);
            //C.AddBoundaryCondition("Pressure_Outlet");

            C.AddBoundaryCondition("Wall");

            // Set Initial Conditions
            //C.InitialValues_Evaluators.Add("VelocityX", X => (2*(2*X[1] -1)*(1-((2*X[0] -1)* (2 * X[0] - 1)))));
            //C.InitialValues_Evaluators.Add("VelocityY", X => (-2 * (2 * X[0] - 1) * (1 - ((2 *X[1] - 1)) * (2 *X[1] - 1))));
            //C.InitialValues_Evaluators.Add("Pressure", X => 2.0*C.PhysicalParameters.mu_A*(-X[0] + 10));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0.0);
            C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + 0.5.Pow2());

            return C;
        }

        static public IBM_Control BackwardStep(int k = 2, int cellsX = 20, int cellsY = 10, string dbpath = null) {
            IBM_Control C = new IBM_Control();


            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication<AppControl> app, int noOfPerformanceClasses) {
                Console.WriteLine("i was called");
                int[] map = new int[] { 1, 5, 100 };
                return new StaticCellCostEstimator(map);
            });

            C.DbPath = @"P:\BoSSS_DBs\Bug";

            // Solver Options
            C.MaxSolverIterations = 100;
            C.MinSolverIterations = 1;
            C.savetodb = false;
            C.ProjectName = "BackwardStep";
            C.SessionName = "BackwardStep";
            C.NoOfMultigridLevels = 5;

            // Calculate Navier-Stokes? 
            C.PhysicalParameters.IncludeConvection = true;

            // Timestepper
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 60;
            C.NoOfTimesteps = 5;

            // Physical values
            C.PhysicalParameters.rho_A = 1;

            // 1/Re
            C.PhysicalParameters.mu_A = 2.0 / 10;


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
                var _xNodes = GenericBlas.Linspace(-1, 5, cellsX + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cellsY + 1);

                double[] CutOutPoint1 = new double[2] { -1, -1 };
                double[] CutOutPoint2 = new double[2] { 0, 0 };

                var CutOut = new BoundingBox(2);
                CutOut.AddPoint(CutOutPoint1);
                CutOut.AddPoint(CutOutPoint2);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear,CutOuts:CutOut);


                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if ((Math.Abs(y - (0)) < 1.0e-6) && (x< 1.0e-6))
                        // bottom step
                        return 2;

                    if (Math.Abs(y - (-1)) < 1.0e-6) 
                        // bottom
                        return 2;

                    if ((Math.Abs(x - (0)) < 1.0e-6) && y< 1.0e-6)
                        // bottom step
                        return 2;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top
                        return 2;

                    if ((Math.Abs(x - (-1)) < 1.0e-6) && (y> 1.0e-6))
                        // left
                        return 1;

                    if (Math.Abs(x - (5)) < 1.0e-6)
                        // right
                        return 3;
                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("Backward Step Flow");

                return grd;
            };

            Func<double[], double, double> VelocityXex, VelocityYex, Pressure;
            VelocityXex = (X, t) => (1 - (X[1] * X[1]));
            VelocityYex = (X, t) => (0);
            Pressure = (X, t) => (0);



            C.AddBoundaryCondition("Velocity_inlet", "VelocityX", X => -4*X[1]*(X[1]+4));
            C.AddBoundaryCondition("Pressure_Outlet");

            C.AddBoundaryCondition("Wall");

            // Set Initial Conditions
            //C.InitialValues_Evaluators.Add("VelocityX", X => (2*(2*X[1] -1)*(1-((2*X[0] -1)* (2 * X[0] - 1)))));
            //C.InitialValues_Evaluators.Add("VelocityY", X => (-2 * (2 * X[0] - 1) * (1 - ((2 *X[1] - 1)) * (2 *X[1] - 1))));
            //C.InitialValues_Evaluators.Add("Pressure", X => 2.0*C.PhysicalParameters.mu_A*(-X[0] + 10));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0.0);
            C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + 0.5.Pow2());

            return C;
        }

    }
}
