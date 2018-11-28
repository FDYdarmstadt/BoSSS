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
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.FSI_Solver
{
    public class IBMPaperReproduce : IBM_Solver.HardcodedTestExamples
    {
        public static FSI_Control IBMCylinderFlow(string _DbPath = null, int k = 2, double Re = 20, bool xPeriodic = false)
        {

            FSI_Control C = new FSI_Control();

            C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("k", k),
                            };


            const double BaseSize = 1.0;

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = true;
            C.ProjectDescription = "IBMCylinder_k" + k + "_Re" + Re;
            C.ProjectName = "IBMCylinder_k"+k+"_Re"+Re;
            C.SessionName = "IBMCylinder_k" + k + "_Re" + Re;

            switch (k)
            {
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
            
            C.Tags.Add("IBMCylinderFlow");
            C.Tags.Add("k"+k);

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

            //grid and boundary conditions
            // ============================

            C.GridFunc = delegate
            {

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
                if (!xPeriodic)
                {
                    grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X)
                {
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
            if (!xPeriodic)
            {
                C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX", X => (4 * 1.5 * (X[1] + 2) * (4.1 - (X[1] + 2)) / (4.1 * 4.1)));
                //C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX#A", X => 1);   
            }
            C.AddBoundaryValue("Pressure_Outlet_right");


            // Initial Values
            // ==============

            double radius = 0.5;
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1.0 / Re;

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
            double dt = new double();
            if(Re==20)
            {
                dt = 10E20;
                C.NoOfTimesteps = 1;
            }
            else if (Re == 100)
            {
                dt = 0.05;
                C.NoOfTimesteps = 1000000;
            }
            else throw new ApplicationException();

            C.Timestepper_Scheme = FSI_Control.TimesteppingScheme.BDF2;
            C.Timestepper_Mode = FSI_Control.TimesteppingMode.None;
           
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 70;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control ParticleInShearFlow(string _DbPath = null, int k = 2, double VelXBase = 0.0, double particleRadius = 1)
        {
            FSI_Control C = new FSI_Control();

            const double BaseSize = 1.0;

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = true;
            C.ProjectName = "ShearFlow_k" + k + "_particleRadius" + particleRadius;
            C.SessionName = "ShearFlow_k" + k + "_particleRadius" + particleRadius;
            C.ProjectDescription = "ShearFlow_k" + k + "_particleRadius" + particleRadius;

            C.Tags.Add("ParticleInShearFlow");
            C.Tags.Add("k"+k);

            // Timesteps
            // ==========
            double dt;
            if (particleRadius == 1) dt = 1;
            else if (particleRadius == 0.4) dt = 0.5;
            else if (particleRadius == 0.2) dt = 0.25;
            else throw new ApplicationException();

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

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate
            {
                double[] Xnodes = GenericBlas.Linspace(-2 * BaseSize, 2 * BaseSize, 21);
                double[] Ynodes = GenericBlas.Linspace(-3 * BaseSize, 3 * BaseSize, 31);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: true);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-2 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-2 * BaseSize)) <= 1.0e-8)
                        et = 2;

                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };




            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY", X => 0.02);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY", X => -0.02);

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.LevelSetMovement = "coupled";
            C.includeTranslation = false;
            C.includeRotation = true;

            // Particle Properties
            double particleDensity = 1;
            //C.particleRho = 1;
            C.particleRadius = particleRadius;
            //C.particleMass = Math.PI * C.particleRadius * C.particleRadius * C.particleRho;
            //C.particleMass = Math.PI * C.particleRadius * C.particleRadius * particleDensity;

            Func<double, double> yLevSet = t => (t * t);
            Func<double[], double, double> phi = (X, t) => -(X[0]).Pow2() + -(X[1]).Pow2() + C.particleRadius.Pow2();
            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;

            C.InitialValues_Evaluators.Add("Phi", X => phi(X, 0));
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("fec14187-4e12-43b6-af1e-e9d535c78668"), -1);


            // Physical Parameters
            // ===================

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.01;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxSolverIterations = 100;
            C.NoOfMultigridLevels = 1;

            // Timestepping
            // ============

            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            C.Timestepper_Mode = FSI_Control.TimesteppingMode.None;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 500;
            C.NoOfTimesteps = 2500;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control IBMCylinderFlowUhlmann(string _DbPath = null, int k = 2, bool xPeriodic = false, double VelXBase = 0.0, bool movingMesh = true)
        {
            FSI_Control C = new FSI_Control();

            //const double BaseSize = 1.0;
            //const double MeshFactor = 0.43;

            // basic database options
            // ======================



            C.DbPath = _DbPath;
            C.savetodb = true;
            C.saveperiod = 300;
   
            C.Tags.Add("OscillatingCylinder");
            C.Tags.Add("k"+k);

            if (movingMesh)
            {
                C.ProjectDescription = "OscillatingCylinder_k" + k + "_MM";
                C.ProjectName = "OscillatingCylinder_k" + k + "_MM";
                C.SessionName = "OscillatingCylinder_k" + k + "_MM";
                C.Tags.Add("MM");
            }
            else
            {
                C.ProjectDescription = "OscillatingCylinder_k" + k + "_SP";
                C.ProjectName = "OscillatingCylinder_k" + k + "_SP";
                C.SessionName = "OscillatingCylinder_k" + k + "_SP";
                C.Tags.Add("SP");
            }

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

            // grid and boundary conditions
            // ============================

            //C.GridFunc = delegate {



            //    var _xNodes1 = Grid1D.TanhSpacing(-6.17, -1, Convert.ToInt32(15 * MeshFactor), 2.1, false); //15
            //    _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
            //    var _xNodes2 = GenericBlas.Linspace(-1, 1.5, Convert.ToInt32(50 * MeshFactor)); //50
            //    _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
            //    var _xNodes3 = Grid1D.TanhSpacing(1.5, 20.5, Convert.ToInt32(50 * MeshFactor), 2, true); //50

            //    var xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3);


            //    var _yNodes1 = Grid1D.TanhSpacing(-13.3, -1.2, Convert.ToInt32(15 * MeshFactor), 2.5, false); //15
            //    _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
            //    var _yNodes2 = GenericBlas.Linspace(-1.2, 1.2, Convert.ToInt32(40 * MeshFactor)); //40
            //    _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
            //    var _yNodes3 = Grid1D.TanhSpacing(1.2, 13.3, Convert.ToInt32(15 * MeshFactor), 2.5, true); //15

            //    var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);

            //    //double[] xNodes = GenericBlas.Linspace(-6.17, 20.5, 50);
            //    //double[] yNodes = GenericBlas.Linspace(-13.3, 13.3, 50);
            //    var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: xPeriodic);
            //    grd.EdgeTagNames.Add(1, "Velocity_Inlet_lower");
            //    grd.EdgeTagNames.Add(2, "Velocity_Inlet_upper");
            //    if (!xPeriodic) {
            //        grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
            //        grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
            //    }
            //    //grd.EdgeTagNames.Add(1, "Outflow_lower");
            //    //grd.EdgeTagNames.Add(2, "Outflow_upper");
            //    //if (!xPeriodic)
            //    //{
            //    //    grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
            //    //    grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
            //    //}

            //    grd.DefineEdgeTags(delegate (double[] X) {
            //        byte et = 0;
            //        if (Math.Abs(X[1] - (-13.3)) <= 1.0e-8)
            //            et = 1;
            //        if (Math.Abs(X[1] + (-13.3)) <= 1.0e-8)
            //            et = 2;
            //        if (!xPeriodic && Math.Abs(X[0] - (-6.17)) <= 1.0e-8)
            //            et = 3;
            //        if (!xPeriodic && Math.Abs(X[0] + (-20.5)) <= 1.0e-8)
            //            et = 4;


            //        Debug.Assert(et != 0);
            //        return et;
            //    });

            //    Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

            //    return grd;
            //};

            C.GridFunc = delegate
            {
                int q2 = new int();
                int r2 = new int();
                int q3 = new int();
                int r3 = new int();

                switch (k)
                {
                   case 1:
                        q2 = 35;
                        r2 = 25;
                        q3 = 105;
                        r3 = 75;
                        break;

                    case 2:
                        q2 = 28;
                        r2 = 20;
                        q3 = 64;
                        r3 = 48;
                        break;

                    case 3:
                        q2 = 21;
                        r2 = 15;
                        q3 = 39;
                        r3 = 36;
                        break;

                    default:

                        throw new ApplicationException();
                    }

                // Box1
                var box1_p1 = new double[2] { -6.17, -13.3 };
                var box1_p2 = new double[2] { 20.5, 13.3 };
                var box1 = new GridCommons.GridBox(box1_p1, box1_p2, 15, 15); //k1: ; k2: 15,15; k3: 15,15

                // Box2
                var box2_p1 = new double[2] { -2, -4 };
                var box2_p2 = new double[2] { 10, 4 };
                var box2 = new GridCommons.GridBox(box2_p1, box2_p2, q2, r2); //k1: 35,25 ; k2: 28,20; k3: 21, 15


                // Box3
                var box3_p1 = new double[2] { -1.5, -2.5 };
                var box3_p2 = new double[2] { 6, 2.5 };
                var box3 = new GridCommons.GridBox(box3_p1, box3_p2, q3, r3); //k1: 105, 75 ; k2: 64,48; k3: 39, 36

                var grd = Grid2D.HangingNodes2D(box1, box2, box3);

                grd.EdgeTagNames.Add(1, "Pressure_Outlet_lower");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_upper");
                if (!xPeriodic)
                {
                    grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[1] - (-13.3)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] + (-13.3)) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic && Math.Abs(X[0] - (-6.17)) <= 1.0e-8)
                        et = 3;
                    if (!xPeriodic && Math.Abs(X[0] + (-20.5)) <= 1.0e-8)
                        et = 4;


                    // Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:    {0}", grd.NumberOfCells);

                return grd;
            };


            C.AddBoundaryValue("Pressure_Outlet_lower");
            C.AddBoundaryValue("Pressure_Outlet_upper");
            if (!xPeriodic)
            {
                C.AddBoundaryValue("Velocity_Inlet_left", "VelocityX", X => 1);
            }
            C.AddBoundaryValue("Pressure_Outlet_right");

            //C.AddBoundaryCondition("Outflow_lower");
            //C.AddBoundaryCondition("Outflow_upper");
            //if (!xPeriodic)
            //{
            //    C.AddBoundaryCondition("Velocity_Inlet_left", "VelocityX#A", X => 1);
            //}
            //C.AddBoundaryCondition("Pressure_Outlet_right");

            // Level-Set Movement
            // ===================

            double radius = 0.5;
            C.includeRotation = false;
            C.includeTranslation = false;
            C.LevelSetMovement = "fixed";
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1.0 / 185;


            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) =>0 };

            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0.25 * Math.PI * 2 * 0.166 * -Math.Sin(Math.PI * 2 * 0.195 * t) };

            Func<double, double> yLevSet = t => (0.2 * Math.Cos(Math.PI * 2 * 0.156 * t));
            Func<double[], double, double> phi = (X, t) => -(X[0]).Pow2() + -(X[1] - yLevSet(t)).Pow2() + radius.Pow2();
            //Func<double[], double, double> phi = (X, t) => -(X[0]).Pow2() + -(X[1]-1).Pow2() + radius.Pow2();
            C.MovementFunc = phi;

            Func<double, double> xVelocity = t => 0;
            Func<double, double> yVelocity = t => (0.2 * Math.PI * 2 * 0.156 * -Math.Sin(Math.PI * 2 * 0.156 * t));
            Func<double, double>[] particleTransVelocity = { xVelocity, yVelocity };
            Func<double, double>[] particleAnglVelocity = { xVelocity, xVelocity };

            C.transVelocityFunc = particleTransVelocity;
            C.anglVelocityFunc = particleAnglVelocity;

            // Initial Values
            // ==============           

            C.InitialValues_Evaluators.Add("Phi", X => phi(X, 0));
            ////C.InitialValues.Add("Phi", X => -1);

            C.InitialValues_Evaluators.Add("VelocityX", X => 1);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => 0);
            ////C.InitialValues.Add("VelocityY#A", X => osciVelocity(X, 0));
            //C.InitialValues.Add("VelocityX", delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];

            //    double R = Math.Sqrt((x + 1).Pow2() + y.Pow2());

            //    double xVel = 1;

            //    if (R < 0.75) {
            //        xVel = 1;
            //    }
            //    return xVel;
            //});

            ////C.InitialValues.Add("VelocityY", delegate (double[] X)
            ////{
            ////    double x = X[0];
            ////    double y = X[1];

            ////    double R = Math.Sqrt((x + 2).Pow2() + y.Pow2());

            ////    double yVel = 0;

            ////    if (R < 0.75)
            ////    {
            ////        yVel = 1;
            ////    }
            ////    return yVel;
            ////});

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("8bef674a-7e37-45a0-9ee2-81a7f0cd02eb"), 1500);
            //C.GridGuid = new Guid("be76c5d5-010c-41a5-b342-e87b42d9734e");

            // Physical Parameters
            // ===================
            C.PhysicalParameters.IncludeConvection = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxKrylovDim = 20;
            C.MaxSolverIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NoOfMultigridLevels = 0;

            // Timestepping
            // ============

            if(movingMesh) {
                C.Timestepper_Mode = FSI_Control.TimesteppingMode.MovingMesh;
            } else
            {
                C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            }
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 300;
            C.NoOfTimesteps = 1000000000;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control ParticleUnderGravity(string _DbPath = null, int k = 2, double VelXBase = 0.0, bool movingMesh = true, bool restart = false)
        {
            //List<FSI_Control> R = new List<FSI_Control>();

            // foreach (int i in new int[] {1,2, 3 }) {
            string restartSession = "e73d770a-d26f-412b-b4f2-c68421898e9e";
            string restartGrid = "dff0fdc4-fc46-4e94-acc3-9ad99e7be5cf";

            FSI_Control C = new FSI_Control();

            const double BaseSize = 1.0;

            C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("k", k),
                            };

            // k = i;

            // basic database options
            // ======================

            //C.DbPath = _DbPath;
            C.DbPath = @"\\dc1\userspace\stange\HiWi_database\ParticleUnderGravity";
            C.savetodb = true;
            C.saveperiod = 1;

            C.Tags.Add("ParticleUnderGravity");
            C.Tags.Add("k"+k);
            C.Tags.Add("restart_" + restart);

            if (movingMesh)
            {

                C.ProjectDescription = "ParticleUnderGravity_dt0.001_" + k + "_MM";
                C.ProjectName = "ParticleUnderGravity_dt0.0001_k" + k + "_MM";
                C.SessionName = "ParticleUnderGravity_dt0.001_k" + k + "_MM_MFVOneStepGaussAndStokes"; //_MFVOneStepGaussAndStokes
                C.Tags.Add("MM");
            }
            else
            {
                C.ProjectDescription = "ParticleUnderGravity_dt0.0001_k" + k + "_SP";
                C.ProjectName = "ParticleUnderGravity_dt0.0001_k" + k + "_SP";
                C.SessionName = "ParticleUnderGravity_dt0.001_k" + k + "_SP_MFVOneStepGaussAndStokes_DoF150000"; 
                C.Tags.Add("SP");
            }
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

            // grid and boundary conditions
            // ============================
            if (!restart)
            {
            C.GridFunc = delegate
            {

                int q = new int();
                int r = new int();

                switch (k)
                {
                    case 1:
                        q = 113; //60; DoF 150000: 113-191
                        r = 191; //178;
                        break;

                    case 2:
                        q = 41; //41; DoF 150000: 71-141
                        r = 121; //121;
                        break;

                    case 3:
                        q = 31;//45;//31;
                        r = 91;//129;//91;
                        break;

                    default:

                        throw new ApplicationException();
                }

                double[] Xnodes = GenericBlas.Linspace(-1 * BaseSize, 1 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(0 * BaseSize, 6 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X)
                {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (0 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-6 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };

                C.Particles[0] = new Particle(2, 4) {
                    radius_P = 0.125,
                    rho_P = 1.25,
                
                };
                C.Particles[0].currentPos_P[0] = new double[] { 0.0, 4.0 };

                //Func<double[], double, double> phi = (X, t) => -(X[0] - C.initialPos[0][0]).Pow2() + -(X[1] - C.initialPos[0][1]).Pow2() + C.particleRadius.Pow2();
                //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
                //C.MovementFunc = phi;

                C.InitialValues_Evaluators.Add("Phi", X => C.Particles[0].phi_P(X, 0));
                //C.InitialValues.Add("VelocityX#B", X => 1);
                C.InitialValues_Evaluators.Add("VelocityX", X => 0);
                C.InitialValues_Evaluators.Add("VelocityY", X => 0);
                //C.InitialValues.Add("Phi", X => -1);
                //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

            }


            C.AddBoundaryValue("Velocity_Inlet_left", "VelocityY", X => 0);
            C.AddBoundaryValue("Velocity_Inlet_right", "VelocityY", X => 0);
            C.AddBoundaryValue("Wall_lower");
            C.AddBoundaryValue("Pressure_Outlet");

            // Boundary values for level-set
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0.1 * 2 * Math.PI * -Math.Sin(Math.PI * 2 * 1 * t), (t) =>  0};
            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0 };

            // Initial Values
            // ==============

            // Coupling Properties
            C.LevelSetMovement = "coupled";
            C.includeRotation = false;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;

            // Particle Properties
            //C.particleRho = 1.25; // 1.25;
            //C.PhysicalParameters.mu_B = 0.1;
            C.particleRadius = 0.125;
            //C.particleMass = Math.PI * C.particleRadius * C.particleRadius * C.particleRho;
            //C.particleMass = 1;



            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("42c82f3c-bdf1-4531-8472-b65feb713326"), 400);
            //C.GridGuid = new Guid("f1659eb6-b249-47dc-9384-7ee9452d05df");


            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.MaxSolverIterations = 100;


            // Timestepping
            // ============

            if (movingMesh)
            {
                C.Timestepper_Mode = FSI_Control.TimesteppingMode.MovingMesh;
            }
            else
            {
                C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            }
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1;
            C.NoOfTimesteps = 1000000;

            // haben fertig...
            // ===============

            return C;
           // R.Add(C);
            //}

            //return R.ToArray();
        }

    }


}
