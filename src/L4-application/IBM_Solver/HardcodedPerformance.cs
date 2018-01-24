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
using BoSSS.Foundation.IO;
using BoSSS.Solution;

namespace BoSSS.Application.IBM_Solver {
    public class HardcodedPerformance {

        static public IBM_Control SphereFlow(string _DbPath = null, int k = 2, int cells_x = 11, int cells_yz = 9, bool only_channel = true, bool pardiso = true, int no_p = 1, int no_it = 1, bool restart = false, bool load_Grid = false, string _GridGuid = null) {
            IBM_Control C = new IBM_Control();

            // basic database options
            // ======================
            //C.DbPath = _DbPath;
            C.savetodb = true;
            //C.savetodb = false;

            //C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Bug";
            C.DbPath = @"/home/ws35kire/test_db/";

            //string restartSession = "727da287-1b6a-463e-b7c9-7cc19093b5b3";
            //string restartGrid = "3f8f3445-46f1-47ed-ac0e-8f0260f64d8f";

            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication app, int noOfPerformanceClasses) {
                Console.WriteLine("i was called");
                int[] map = new int[] { 1, 5, 100 };
                return new StaticCellCostEstimator(map);
            });



            // Assign correct names

            if (pardiso)
            {
                if (only_channel)
                {
                    C.SessionName = "Channel_Pardiso_k" + k + "_" + cells_x + "x" + cells_yz + "x" + cells_yz + "_no_p" + no_p + "_run" + no_it;

                }
                else
                {
                    C.SessionName = "Sphere_Pardiso_k" + k + cells_x + "x" + cells_yz + "x" + cells_yz + "_no_p" + no_p + "_run" + no_it;
                }
            }
            else
            {
                if (only_channel)
                {
                    C.SessionName = "Channel_Mumps_k" + k + cells_x + "x" + cells_yz + "x" + cells_yz + "_no_p" + no_p + "_run" + no_it;
                }
                else
                {
                    C.SessionName = "Sphere_Mumps_k" + k + cells_x + "x" + cells_yz + "x" + cells_yz + "_no_p" + no_p + "_run" + no_it;
                }
            }
            C.saveperiod = 1;
            //C.SessionName = "Sphere_k" + k + "_h" + h+"Re100";
            C.ProjectName = "Sphere3D_Stokes";
            C.ProjectDescription = "Sphere_k"+k + cells_x + "x" + cells_yz + "x" + cells_yz;
            C.Tags.Add("with immersed boundary method");
            C.Tags.Add("Pardiso " + pardiso);
            C.Tags.Add("only channel " + only_channel);
            C.Tags.Add("k " + k);
            C.Tags.Add("no_p" + no_p);
            C.Tags.Add("run " + no_it);
            C.Tags.Add("cells_x " + cells_x);
            C.Tags.Add("cells_yz " + cells_yz);
            C.Tags.Add("restart " + restart);

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

            //if (restart)
            //{
            //    C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(restartSession), -1);
            //    C.GridGuid = new Guid(restartGrid);
            //}
            // Load Grid
            if (!restart)
            {

            if (load_Grid == true) {
                Console.WriteLine("...loading grid");
                C.GridGuid = new Guid(_GridGuid);
            } else {
                #region Creates grid () and sets BC
                //// Create Grid
                Console.WriteLine("...generating grid");
                C.GridFunc = delegate {

                    // x-direction
                    var _xNodes = GenericBlas.Linspace(-10, 30, cells_x + 1);

                    // y-direction
                    var _yNodes = GenericBlas.Linspace(-10, 10, cells_yz + 1);

                    // z-direction
                    var _zNodes = GenericBlas.Linspace(-10, 10, cells_yz + 1);

                    // Cut Out
                    var grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, false, false, false, CellType.Cube_Linear);

                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(2, "Pressure_Outlet");
                    // grd.EdgeTagNames.Add(3, "Wall");

                    grd.DefineEdgeTags(delegate (double[] _X) {
                        var X = _X;
                        double x = X[0];
                        double y = X[1];
                        double z = X[2];

                        if (Math.Abs(x - (-10)) < 1.0e-6)
                            // inlet
                            return 1;

                        if (Math.Abs(x - (30)) < 1.0e-6)
                            // outlet
                            return 2;

                        if (Math.Abs(y - (-10)) < 1.0e-6)
                            // left
                            return 2;

                        if (Math.Abs(y - (10)) < 1.0e-6)
                            // right
                            return 2;

                        if (Math.Abs(z - (-10)) < 1.0e-6)
                            // top left
                            return 2;

                        if (Math.Abs(z - (10)) < 1.0e-6)
                            // top right
                            return 2;

                        throw new ArgumentOutOfRangeException();
                    });

                    return grd;
                };
            }
                #endregion

                //// Create Grid with HANGING NODES
                //Console.WriteLine("...generating grid");
                //C.GridFunc = delegate {

                //    // Box1
                //    var box1_p1 = new double[3] { -10, -10, -10 };
                //    var box1_p2 = new double[3] { 30, 10, 10 };
                //    var box1 = new GridCommons.GridBox(box1_p1, box1_p2,10,5,5);

                //    // Box2
                //    var box2_p1 = new double[3] { 0, -5, -5 };
                //    var box2_p2 = new double[3] { 20, 5, 5 };
                //    var box2 = new GridCommons.GridBox(box2_p1, box2_p2, 10, 6, 6);

                //    // Cut Out
                //    var grd = Grid3D.HangingNodes3D(false, true, true, box1, box2);

                //    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                //    grd.EdgeTagNames.Add(2, "Pressure_Outlet");
                //    grd.EdgeTagNames.Add(3, "Wall");

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
                //            return 1;

                //        if (Math.Abs(y - (10)) < 1.0e-6)
                //            // right
                //            return 1;

                //        if (Math.Abs(z - (-10)) < 1.0e-6)
                //            // top left
                //            return 1;

                //        if (Math.Abs(z - (10)) < 1.0e-6)
                //            // top right
                //            return 1;

                //        throw new ArgumentOutOfRangeException();
                //    });

                //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                //    return grd;
                //};

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

                // Set Initial Conditions
                C.InitialValues_Evaluators.Add("VelocityX", X => 0.5);
                C.InitialValues_Evaluators.Add("VelocityY", X => 0);
                C.InitialValues_Evaluators.Add("VelocityZ", X => 0.5);
                C.InitialValues_Evaluators.Add("Pressure", X => 0);

                if (only_channel)
                {
                    C.InitialValues_Evaluators.Add("Phi", X => -1);
                }
                else
                {
                    C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + -(X[2]).Pow2() + C.particleRadius.Pow2());
                }

            }
            Console.WriteLine("...starting calculation of Sphere3D");

            // Initial Solution

            // Physical values
            C.particleRadius = 2.5;
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 2.5*2 / 100;

            // Boundary conditions
            C.AddBoundaryCondition("Velocity_inlet", "VelocityX", (X, t) => 1);
            C.AddBoundaryCondition("Velocity_inlet", "VelocityY", (X, t) => 0);
            //C.AddBoundaryCondition("Velocity_inlet", "VelocityZ", (X, t) => 0);
           // C.AddBoundaryCondition("Wall");
            C.AddBoundaryCondition("Pressure_Outlet");

            
            // misc. solver options
            // ====================
            C.PhysicalParameters.IncludeConvection = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 20;
            C.MaxSolverIterations = 50;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;

            // Timestepping
            // ============

            if (pardiso)
            {
                C.whichSolver = DirectSolver._whichSolver.PARDISO;
            }
            else
            {
                C.whichSolver = DirectSolver._whichSolver.MUMPS;
            }
            //C.whichSolver = DirectSolver._whichSolver.MUMPS;
            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10000000;
            //C.NoOfTimesteps = 10;
            C.NoOfTimesteps = 1;
            C.NoOfMultigridLevels = 3;

            return C;
        }
        

        static public IBM_Control IBMCylinderFlow(string _DbPath = null, int k = 2, bool only_channel = true, bool pardiso = true, int no_p = 1, int no_it = 1, bool load_Grid = false, string _GridGuid = null)
        {
            // int cells_x, int cells_yz
            IBM_Control C = new IBM_Control();
            bool xPeriodic = false;
            int i = 2;
            const double BaseSize = 1.0;

            // basic database options
            // ======================

            //C.DbPath = _DbPath;

            C.DbPath = @"\\dc1\userspace\stange\HiWi_database\PerformanceTests";
            C.savetodb = true;

            bool restart = false;
            string restartSession = "67a29dcc-ade9-4704-b198-b3380e774f5a";
            string restartGrid = "42e1ede0-40fc-4267-9d48-94c0397ac9a5";
            bool startFromGivenGrid = true;
            string startGrid = "42e1ede0-40fc-4267-9d48-94c0397ac9a5";
            switch (i)
            {
                case 1:
                    C.MeshFactor = 1.258; // was 1.33
                    break;

                case 2:
                    C.MeshFactor = 3.0; //1.77; //0.92;
                    break;

                case 3:
                    C.MeshFactor = 0.7; // was 07
                    break;

                default:

                    throw new ApplicationException();
            }

            if (pardiso)
            {
                if (only_channel)
                {
                    C.SessionName = "2DChannel_Pardiso_k" + k + "_MeshFactor" + C.MeshFactor + "_no_p" + no_p + "_run" + no_it;

                }
                else
                {
                    C.SessionName = "Cylinder_Pardiso_k" + k + "_MeshFactor" + C.MeshFactor + "_no_p" + no_p + "_run" + no_it;
                }
            }
            else
            {
                if (only_channel)
                {
                    C.SessionName = "2DChannel_Mumps_k" + k + "_MeshFactor" + C.MeshFactor + "_no_p" + no_p + "_run" + no_it;
                }
                else
                {
                    C.SessionName = "Cylinder_Mumps_k" + k + "_MeshFactor" + C.MeshFactor + "_no_p" + no_p + "_run" + no_it;
                }
            }
            C.saveperiod = 1;
            //C.SessionName = "Sphere_k" + k + "_h" + h+"Re100";


            C.Tags.Add("Pardiso " + pardiso);
            C.Tags.Add("only channel " + only_channel);
            C.Tags.Add("k " + k);
            C.Tags.Add("no_p" + no_p);
            C.Tags.Add("run " + no_it);
            C.Tags.Add("MeshFactor " + C.MeshFactor);

            C.ProjectName = "FixedCylinderRe100_k" + i + "_CellAgglo02_penalty4_newMesh2";

            C.ProjectDescription = "Cylinder";

            // DG degrees
            // ==========

            C.FieldOptions.Add("VelocityX", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts()
            {
                Degree = k,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            //Console.WriteLine("Achtung: equal order!!!!");
            C.FieldOptions.Add("Pressure", new FieldOpts()
            {
                Degree = k - 1,

                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // restart options
            // ===============
            if (restart)
            {
                C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(restartSession), -1);
                C.GridGuid = new Guid(restartGrid);
            }
            //grid and boundary conditions
            // ============================

            // Initial Values
            // ==============

            double radius = 0.5;
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 1.0 / 100.0;

            if (!restart)
            {
                if (!startFromGivenGrid)
                {
                    C.GridFunc = delegate
                    {

                        var _xNodes1 = Grid1D.TanhSpacing(-2.0, -1.0, Convert.ToInt32(10.0 * C.MeshFactor), 0.5, false); //10
                        _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
                        var _xNodes2 = GenericBlas.Linspace(-1.0, 2.0, Convert.ToInt32(35.0 * C.MeshFactor)); //35
                        _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
                        var _xNodes3 = Grid1D.TanhSpacing(2.0, 20.0, Convert.ToInt32(60.0 * C.MeshFactor), 1.5, true); //60

                        var xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3);


                        var _yNodes1 = Grid1D.TanhSpacing(-2.0, -1.0, Convert.ToInt32(7.0 * C.MeshFactor), 0.9, false); //7
                        _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                        var _yNodes2 = GenericBlas.Linspace(-1.0, 1.0, Convert.ToInt32(25.0 * C.MeshFactor)); //25
                        _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
                        var _yNodes3 = Grid1D.TanhSpacing(1.0, 2.1, Convert.ToInt32(7.0 * C.MeshFactor), 1.1, true); //7
                        var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);



                        //double[] xNodes = GenericBlas.Linspace(0 * BaseSize, 22 * BaseSize, 25);
                        //double[] yNodes = GenericBlas.Linspace(0 * BaseSize, 4.1 * BaseSize, 25);
                        var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: xPeriodic);
                        grd.EdgeTagNames.Add(1, "Velocity_Inlet_upper");
                        grd.EdgeTagNames.Add(2, "Velocity_Inlet_lower");
                        if (!xPeriodic)
                        {
                            grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                            grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                        }

                        grd.DefineEdgeTags(delegate (double[] X)
                        {
                            byte et = 0;
                            if (Math.Abs(X[1] - (-2.0 * BaseSize)) <= 1.0e-8)
                                et = 1;
                            if (Math.Abs(X[1] - (+2.1 * BaseSize)) <= 1.0e-8)
                                et = 2;
                            if (!xPeriodic && Math.Abs(X[0] - (-2.0 * BaseSize)) <= 1.0e-8)
                                et = 3;
                            if (!xPeriodic && Math.Abs(X[0] - (+20.0 * BaseSize)) <= 1.0e-8)
                                et = 4;


                            Debug.Assert(et != 0);
                            return et;
                        });

                        Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                        return grd;
                    };
                }
                else { C.GridGuid = new Guid(startGrid); }
                if (only_channel)
                {
                    C.InitialValues_Evaluators.Add("Phi", X => -1);
                }
                else
                {
                    C.InitialValues_Evaluators.Add("Phi", X => -(X[0]).Pow2() + -(X[1]).Pow2() + radius.Pow2());
                }


                //C.InitialValues.Add("Phi", X => -1);

                C.InitialValues_Evaluators.Add("VelocityX", X => 4.0 * 1.5 * (X[1] + 2.0) * (4.1 - (X[1] + 2.0)) / (4.1 * 4.1));
            }

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

            C.AddBoundaryCondition("Velocity_Inlet_upper", "VelocityX", X => 0.0);
            C.AddBoundaryCondition("Velocity_Inlet_lower", "VelocityX", X => 0.0); //-(4 * 1.5 * X[1] * (4.1 - X[1]) / (4.1 * 4.1))
            if (!xPeriodic)
            {
                C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX", X => (4.0 * 1.5 * (X[1] + 2.0) * (4.1 - (X[1] + 2.0)) / (4.1 * 4.1)));
                //C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX#A", X => 1);   
            }
            C.AddBoundaryCondition("Pressure_Outlet_right");


            

            //C.InitialValues.Add("Phi", X => phi(X, 0));

            //C.InitialValues.Add("Phi", X => ((X[0] / (radius * BaseSize)) - mPx) * (X[0] / (radius * BaseSize)) - mPx) + ((X[1]) / (radius * BaseSize)) - 2.)Pow2() - radius.Pow2()));  // quadratic form
            //    );

            
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
            //C.NoOfMultigridLevels = 0;

            if (pardiso)
            {
                C.whichSolver = DirectSolver._whichSolver.PARDISO;
            }
            else
            {
                C.whichSolver = DirectSolver._whichSolver.MUMPS;
            }
            // Timestepping
            // ============

            C.Timestepper_Scheme = IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 70;
            C.NoOfTimesteps = 10;

            // haben fertig...
            // ===============
            return C;
            
        }


    }
}
