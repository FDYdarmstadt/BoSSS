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

using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid.RefElements;
using System.Diagnostics;
using BoSSS.Solution.Multigrid;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.IBM_Solver {
    public class HardcodedControl3D {


        public static IBM_Control[] IBMCylinderFlow(string _DbPath = null, int k = 2, bool xPeriodic = false, double VelXBase = 0.0) {
            List<IBM_Control> R = new List<IBM_Control>();

            foreach (int i in new int[] { 3 }) {

                IBM_Control C = new IBM_Control();

                C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("k", i),
                            };

                k = i;

                const double BaseSize = 1.0;

                // basic database options
                // ======================

                C.DbPath = @"\\fdyprime\userspace\krause\BoSSS_DBs\Paper_CellAgglo01_Penalty4";
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

                C.AddBoundaryValue("Velocity_Inlet_upper", "VelocityX", X => 0);
                C.AddBoundaryValue("Velocity_Inlet_lower", "VelocityX", X => 0); //-(4 * 1.5 * X[1] * (4.1 - X[1]) / (4.1 * 4.1))
                if (!xPeriodic) {
                    C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX", X => (4 * 1.5 * (X[1] + 2) * (4.1 - (X[1] + 2)) / (4.1 * 4.1)));
                    //C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX#A", X => 1);   
                }
                C.AddBoundaryValue("Pressure_Outlet_right");


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

                C.LevelSetSmoothing = false;
                //C.option_solver = "direct";
                C.MaxKrylovDim = 20;
                C.MaxSolverIterations = 50;
                C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
                C.NoOfMultigridLevels = 0;

                // Timestepping
                // ============

                C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.ImplicitEuler;
                double dt = 0.05;
                C.dtMax = dt;
                C.dtMin = dt;
                C.Endtime = 70;
                C.NoOfTimesteps = 1000000;

                // haben fertig...
                // ===============

                R.Add(C);
            }

            return R.ToArray();
        }


        static public IBM_Control SphereFlow(int k = 2, int h = 1) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Sphere3D";
            C.savetodb = true;
            C.ProjectName = "Sphere3D";
            C.SessionName = "Sphere3D_" + k + "_Re350";
            C.Tags.Add("with immersed boundary method");

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
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            C.AdaptiveMeshRefinement = false;

            //C.GridPartType = GridPartType.Hilbert;


            C.TimeStepper_Init = Solution.Timestepping.TimeStepperInit.MultiInit;

            // Load Grid
            Console.WriteLine("...loading grid");
            C.GridGuid = new Guid("1a672505-e301-4271-9c7d-050770f48abc");

            C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("d1e5a259-d969-4832-a22e-f8a48b2b7a36"), -1);

            //#region Creates grid () and sets BC
            ////// Create Grid
            //Console.WriteLine("...generating grid");
            //C.GridFunc = delegate {

            //    // x-direction
            //    var _xNodes = GenericBlas.Linspace(-10, 30, (10 * h) + 1);

            //    // y-direction
            //    var _yNodes = GenericBlas.Linspace(-10, 10, (5 * h) + 1);

            //    // z-direction
            //    var _zNodes = GenericBlas.Linspace(-10, 10, (5 * h) + 1);

            //    // Cut Out
            //    var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, false, false, CellType.Cube_Linear);

            //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
            //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");

            //    grd.DefineEdgeTags(delegate (double[] _X) {
            //        var X = _X;
            //        double x = X[0];
            //        double y = X[1];
            //        double z = X[2];

            //        if (Math.Abs(x - (-10)) < 1.0e-6)
            //            // inlet
            //            return 1;

            //        if (Math.Abs(x - (30)) < 1.0e-6)
            //            // outlet
            //            return 2;

            //        if (Math.Abs(y - (-10)) < 1.0e-6)
            //            // left
            //            return 2;

            //        if (Math.Abs(y - (10)) < 1.0e-6)
            //            // right
            //            return 2;

            //        if (Math.Abs(z - (-10)) < 1.0e-6)
            //            // top left
            //            return 2;

            //        if (Math.Abs(z - (10)) < 1.0e-6)
            //            // top right
            //            return 2;

            //        throw new ArgumentOutOfRangeException();
            //    });

            //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

            //    return grd;
            //};
            //#endregion


            // Create Grid with HANGING NODES
            //Console.WriteLine("...generating grid");
            //C.GridFunc = delegate {

            //    // Box1
            //    var box1_p1 = new double[3] { -2, -3, -3 };
            //    var box1_p2 = new double[3] { 10, 3, 3 };
            //    var box1 = new GridCommons.GridBox(box1_p1, box1_p2, 12 * h, 6 * h, 6 * h);

            //    // Box2
            //    var box2_p1 = new double[3] { -1.5, -1, -1 };
            //    var box2_p2 = new double[3] { 10, 1, 1 };
            //    var box2 = new GridCommons.GridBox(box2_p1, box2_p2, 10 * (h + 3), 4 * (h + 3), 4 * (h + 3));

            //    // Cut Out
            //    //var box1 = new GridCommons.GridBox(box1_p1, box1_p2, 10, 5, 5);
            //    //var box2 = new GridCommons.GridBox(box2_p1, box2_p2, 10, 5, 5);
            //    var grd = Grid3D.HangingNodes3D(false, false, false, box1, box2);

            //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
            //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");
            //    //grd.EdgeTagNames.Add(3, "Wall");

            //    grd.DefineEdgeTags(delegate (double[] _X) {
            //        var X = _X;
            //        double x = X[0];
            //        double y = X[1];
            //        double z = X[2];

            //        if (Math.Abs(x - (-2)) < 1.0e-6)
            //            // inlet
            //            return 1;

            //        if (Math.Abs(x - (10)) < 1.0e-6)
            //            // outlet
            //            return 2;

            //        if (Math.Abs(y - (-3)) < 1.0e-6)
            //            // left
            //            return 2;

            //        if (Math.Abs(y - (3)) < 1.0e-6)
            //            // right
            //            return 2;

            //        if (Math.Abs(z - (-3)) < 1.0e-6)
            //            // top left
            //            return 2;

            //        if (Math.Abs(z - (3)) < 1.0e-6)
            //            // top right
            //            return 2;

            //        throw new ArgumentOutOfRangeException();
            //    });

            //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

            //    return grd;
            //};

            //Console.WriteLine("Loading Grid...");
            //C.GridGuid = new Guid("9aca8eed-e98a-4dca-8024-e02edc4d3edc");
            //C.GridGuid = new Guid("099cffa4-238d-42ad-8fbd-85228bfc6b1e");

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

            //    grd.DefineEdgeTags(delegate (double[] _X) {
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

            Console.WriteLine("...starting calculation of Sphere3D");

            // Initial Solution

            // Physical values
            C.particleRadius = 0.5;
            C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.mu_A = (2 * C.particleRadius) / 100;
            C.PhysicalParameters.mu_A = (1.0 * C.particleRadius * 1.0) / 350;

            // Boundary conditions
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", (X, t) => 1);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryValue("Velocity_inlet", "VelocityZ", (X, t) => 0);
            // C.AddBoundaryValue("Wall");
            C.AddBoundaryValue("Pressure_Outlet");

            // Set Initial Conditions
            //C.InitialValues_Evaluators.Add("VelocityX", X => 1);
            ////C.InitialValues_Evaluators.Add("VelocityY", X => 5);
            //C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            //C.InitialValues_Evaluators.Add("Pressure", X => 0);
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + -(X[2]).Pow2() + C.particleRadius.Pow2());

            //C.InitialValues_Evaluators.Add("VelocityY", delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];
            //    double z = X[2];

            //    double R = Math.Sqrt((x + 4).Pow2() + y.Pow2()+z.Pow2());

            //    double yVel = 0;

            //    if (R < 3) {
            //        yVel = 5;
            //    }
            //    return yVel;
            //});
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + C.particleRadius.Pow2());
            //C.InitialValues_Evaluators.Add("Phi", X => -1);




            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 50;
            C.MaxSolverIterations = 50;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NonlinearSolve = NonlinearSolverCodes.NewtonGMRES;
            C.LinearSolve = LinearSolverCodes.exp_schwarz_directcoarse_overlap;
            C.Solver_ConvergenceCriterion = 1E-6;
            C.NoOfMultigridLevels = 2;

            // Timestepping
            // ============


            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000;
            C.NoOfTimesteps = 100000;

            return C;
        }

        static public IBM_Control MittalSphereFlow(int k = 2, int h = 1) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Sphere3D";
            C.savetodb = false;
            C.ProjectName = "Sphere3D";
            C.SessionName = "Sphere3D_" + k + "_Re100";
            C.Tags.Add("with immersed boundary method");

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
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            C.AdaptiveMeshRefinement = false;

            // Load Grid
            Console.WriteLine("...loading grid");
            C.GridGuid = new Guid("202033c5-7c2a-4c29-aaea-86d7a2f7a8a3");

            //#region Creates grid () and sets BC
            ////// Create Grid
            //Console.WriteLine("...generating grid");
            //C.GridFunc = delegate {

            //    // x-direction
            //    var _xNodes = GenericBlas.Linspace(-10, 30, (10 * h) + 1);

            //    // y-direction
            //    var _yNodes = GenericBlas.Linspace(-10, 10, (5 * h) + 1);

            //    // z-direction
            //    var _zNodes = GenericBlas.Linspace(-10, 10, (5 * h) + 1);

            //    // Cut Out
            //    var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, false, false, CellType.Cube_Linear);

            //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
            //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");

            //    grd.DefineEdgeTags(delegate (double[] _X) {
            //        var X = _X;
            //        double x = X[0];
            //        double y = X[1];
            //        double z = X[2];

            //        if (Math.Abs(x - (-10)) < 1.0e-6)
            //            // inlet
            //            return 1;

            //        if (Math.Abs(x - (30)) < 1.0e-6)
            //            // outlet
            //            return 2;

            //        if (Math.Abs(y - (-10)) < 1.0e-6)
            //            // left
            //            return 2;

            //        if (Math.Abs(y - (10)) < 1.0e-6)
            //            // right
            //            return 2;

            //        if (Math.Abs(z - (-10)) < 1.0e-6)
            //            // top left
            //            return 2;

            //        if (Math.Abs(z - (10)) < 1.0e-6)
            //            // top right
            //            return 2;

            //        throw new ArgumentOutOfRangeException();
            //    });

            //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

            //    return grd;
            //};
            //#endregion


            // Create Grid with HANGING NODES
            //Console.WriteLine("...generating grid");
            //C.GridFunc = delegate {

            //    // Box1
            //    var box1_p1 = new double[3] { -6, -7.5, -7.5 };
            //    var box1_p2 = new double[3] { 10, 7.5, 7.5 };
            //    var box1 = new GridCommons.GridBox(box1_p1, box1_p2, 16 * h, 15 * h, 15 * h);

            //    // Box2
            //    var box2_p1 = new double[3] { -2.5, -1.5, -1.5 };
            //    var box2_p2 = new double[3] { 7.5, 1.5, 1.5 };
            //    var box2 = new GridCommons.GridBox(box2_p1, box2_p2, 10 * (h + 3), 4 * (h + 3), 4 * (h + 3));

            //    // Cut Out
            //    //var box1 = new GridCommons.GridBox(box1_p1, box1_p2, 10, 5, 5);
            //    //var box2 = new GridCommons.GridBox(box2_p1, box2_p2, 10, 5, 5);
            //    var grd = Grid3D.HangingNodes3D(false, false, false, box1, box2);

            //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
            //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");
            //    //grd.EdgeTagNames.Add(3, "Wall");

            //    grd.DefineEdgeTags(delegate (double[] _X) {
            //        var X = _X;
            //        double x = X[0];
            //        double y = X[1];
            //        double z = X[2];

            //        if (Math.Abs(x - (-6)) < 1.0e-6)
            //            // inlet
            //            return 1;

            //        if (Math.Abs(x - (10)) < 1.0e-6)
            //            // outlet
            //            return 2;

            //        if (Math.Abs(y - (-7.5)) < 1.0e-6)
            //            // left
            //            return 2;

            //        if (Math.Abs(y - (7.5)) < 1.0e-6)
            //            // right
            //            return 2;

            //        if (Math.Abs(z - (-7.5)) < 1.0e-6)
            //            // top left
            //            return 2;

            //        if (Math.Abs(z - (7.5)) < 1.0e-6)
            //            // top right
            //            return 2;

            //        throw new ArgumentOutOfRangeException();
            //    });

            //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

            //    return grd;
            //};

            //Console.WriteLine("Loading Grid...");
            //C.GridGuid = new Guid("9aca8eed-e98a-4dca-8024-e02edc4d3edc");
            //C.GridGuid = new Guid("099cffa4-238d-42ad-8fbd-85228bfc6b1e");

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

            //    grd.DefineEdgeTags(delegate (double[] _X) {
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

            Console.WriteLine("...starting calculation of Sphere3D");

            // Initial Solution

            // Physical values
            C.particleRadius = 0.5;
            C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.mu_A = (2 * C.particleRadius) / 100;
            C.PhysicalParameters.mu_A = (2.0 * C.particleRadius*1.0) / 100;

            // Boundary conditions
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", (X, t) => 1);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryValue("Velocity_inlet", "VelocityZ", (X, t) => 0);
            // C.AddBoundaryCondition("Wall");
            C.AddBoundaryValue("Pressure_Outlet");

            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("VelocityX", X => 1);
            //C.InitialValues_Evaluators.Add("VelocityY", X => 5);
            C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + -(X[2]).Pow2() + C.particleRadius.Pow2());

            //C.InitialValues_Evaluators.Add("VelocityY", delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];
            //    double z = X[2];

            //    double R = Math.Sqrt((x + 4).Pow2() + y.Pow2()+z.Pow2());

            //    double yVel = 0;

            //    if (R < 3) {
            //        yVel = 5;
            //    }
            //    return yVel;
            //});
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + C.particleRadius.Pow2());
            //C.InitialValues_Evaluators.Add("Phi", X => -1);




            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 30;
            C.MaxSolverIterations = 50;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NonlinearSolve = NonlinearSolverCodes.NewtonGMRES;
            C.LinearSolve = LinearSolverCodes.exp_schwarz_Kcycle_directcoarse_overlap;
            C.Solver_ConvergenceCriterion = 1E-6;
            C.NoOfMultigridLevels = 3;

            // Timestepping
            // ============


            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1E20;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000;
            C.NoOfTimesteps = 1;

            return C;
        }

        static public IBM_Control Simple3DTest(int k = 2, int h = 1) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Sphere3D";
            C.savetodb = false;
            C.ProjectName = "Sphere3D";
            C.SessionName = "Sphere3D_" + k + "_Re350";
            C.Tags.Add("with immersed boundary method");

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
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            C.AdaptiveMeshRefinement = false;

            //C.GridPartType = GridPartType.Hilbert;


            //#region Creates grid () and sets BC
            //// Create Grid
            Console.WriteLine("...generating grid");
            C.GridFunc = delegate {

                // x-direction
                var _xNodes = GenericBlas.Linspace(-1, 1, (3 * h) + 1);

                // y-direction
                var _yNodes = GenericBlas.Linspace(-1, 1, (3 * h) + 1);

                // z-direction
                var _zNodes = GenericBlas.Linspace(-1, 1, (3 * h) + 1);

                // Cut Out
                var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, false, false, CellType.Cube_Linear);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];
                    double z = X[2];

                    if (Math.Abs(x - (-1)) < 1.0e-6)
                        // inlet
                        return 1;

                    if (Math.Abs(x - (1)) < 1.0e-6)
                        // outlet
                        return 2;

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // left
                        return 2;

                    if (Math.Abs(y - (1)) < 1.0e-6)
                        // right
                        return 2;

                    if (Math.Abs(z - (-1)) < 1.0e-6)
                        // top left
                        return 2;

                    if (Math.Abs(z - (1)) < 1.0e-6)
                        // top right
                        return 2;

                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                return grd;
            };

            // Create Grid

            Console.WriteLine("...starting calculation of Sphere3D");

            // Initial Solution

            // Physical values
            C.particleRadius = 0.5;
            C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.mu_A = (2 * C.particleRadius) / 100;
            C.PhysicalParameters.mu_A = (2.0 * C.particleRadius * 1.0) / 100;

            // Boundary conditions
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", (X, t) => 1);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryValue("Velocity_inlet", "VelocityZ", (X, t) => 0);
            // C.AddBoundaryCondition("Wall");
            C.AddBoundaryValue("Pressure_Outlet");

            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("VelocityX", X => 0.2);
            //C.InitialValues_Evaluators.Add("VelocityY", X => 5);
            C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + -(X[2]).Pow2() + C.particleRadius.Pow2());
            C.InitialValues_Evaluators.Add("Phi", X => -1);

            //C.InitialValues_Evaluators.Add("VelocityY", delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];
            //    double z = X[2];

            //    double R = Math.Sqrt((x + 4).Pow2() + y.Pow2()+z.Pow2());

            //    double yVel = 0;

            //    if (R < 3) {
            //        yVel = 5;
            //    }
            //    return yVel;
            //});
            //C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + C.particleRadius.Pow2());
            //C.InitialValues_Evaluators.Add("Phi", X => -1);




            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 4;
            C.MaxSolverIterations = 50;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NonlinearSolve = NonlinearSolverCodes.NewtonGMRES;
            C.LinearSolve = LinearSolverCodes.automatic;
            C.Solver_ConvergenceCriterion = 1E-6;
            C.NoOfMultigridLevels = 2;

            // Timestepping
            // ============


            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000;
            C.NoOfTimesteps = 5;

            return C;
        }

        static public IBM_Control DrivenCavity3D(int k = 2, int cellsXYZ = 1) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Sphere3D";
            C.savetodb = false;
            C.ProjectName = "Sphere3D";
            C.SessionName = "Sphere3D_" + k + "_Re350";
            C.Tags.Add("with immersed boundary method");

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
                Degree = k-1,
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

            C.AdaptiveMeshRefinement = false;

            //C.GridPartType = GridPartType.Hilbert;


            //#region Creates grid () and sets BC
            //// Create Grid
            Console.WriteLine("...generating grid");
            C.GridFunc = delegate {

                // x-direction
                var _xNodes = GenericBlas.Linspace(-0.5, 0.5, cellsXYZ + 1);

                // y-direction
                var _yNodes = GenericBlas.Linspace(-0.5, 0.5, cellsXYZ + 1);

                // z-direction
                var _zNodes = GenericBlas.Linspace(-0.5, 0.5, cellsXYZ + 1);

                // Cut Out
                var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, true, false, CellType.Cube_Linear);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];
                    double z = X[2];

                    if (Math.Abs(x - (-0.5)) < 1.0e-6)
                        // inlet
                        return 2;

                    if (Math.Abs(x - (0.5)) < 1.0e-6)
                        // outlet
                        return 2;

                    if (Math.Abs(y - (-0.5)) < 1.0e-6)
                        // left
                        return 2;

                    if (Math.Abs(y - (0.5)) < 1.0e-6)
                        // right
                        return 2;

                    if (Math.Abs(z - (-0.5)) < 1.0e-6)
                        // top left
                        return 2;

                    if (Math.Abs(z - (0.5)) < 1.0e-6)
                        // top right
                        return 1;

                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                return grd;
            };

            // Create Grid

            Console.WriteLine("...starting calculation of Sphere3D");

            // Initial Solution

            // Physical values
            C.particleRadius = 1;
            C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.mu_A = (2 * C.particleRadius) / 100;
            C.PhysicalParameters.mu_A = (1.0 * C.particleRadius * 1.0) / 400.0;

            // Boundary conditions
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", (X, t) => 1);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (X, t) => 0);
            C.AddBoundaryValue("Velocity_inlet", "VelocityZ", (X, t) => 0);
            C.AddBoundaryValue("Wall");

            // Set Initial Conditions
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            //C.InitialValues_Evaluators.Add("VelocityY", X => 5);
            C.InitialValues_Evaluators.Add("VelocityZ", X => 0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            C.InitialValues_Evaluators.Add("Phi", X => -1);          
            

            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 100;
            C.MaxSolverIterations = 10;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NonlinearSolve = NonlinearSolverCodes.NewtonGMRES;
            C.LinearSolve = LinearSolverCodes.automatic;
            C.Solver_ConvergenceCriterion = 1E-5;
            C.NoOfMultigridLevels = 2;

            // Timestepping
            // ============
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100000;
            C.NoOfTimesteps = 1;

            return C;
        }

    }
}

