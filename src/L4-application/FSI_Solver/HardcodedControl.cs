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
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Foundation.IO;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.FSI_Solver {
    public class HardcodedControl : IBM_Solver.HardcodedTestExamples {

        public static FSI_Control ParticleInShearFlow(string _DbPath = null, int k = 2, double VelXBase = 0.0) {
            FSI_Control C = new FSI_Control();

            const double BaseSize = 1.0;

            // basic database options
            // ======================
            
            //C.DbPath = @"\\dc1\userspace\rieckmann\local\FSI\Test_db";
            C.savetodb = false;

            C.ProjectName = "ShearFlow_Test";
            C.ProjectDescription = "ShearFlow";
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
            C.FieldOptions.Add("Curvature", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-2 * BaseSize, 2 * BaseSize, 21);
                double[] Ynodes = GenericBlas.Linspace(-3 * BaseSize, 3 * BaseSize, 31);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: true);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");


                grd.DefineEdgeTags(delegate (double[] X) {
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
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeTranslation = false;
            C.includeRotation = true;

            // Particle Properties
            //double particleDensity = 1;
            C.particleRadius = 0.4;

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
            C.PhysicalParameters.mu_A = 0.25;
            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.NoOfMultigridLevels = 1;

            // Timestepping
            // ============

            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 120;
            C.NoOfTimesteps = 2;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control IBMCylinderFlowUhlmann(string _DbPath = null, int k = 2, bool xPeriodic = false, double VelXBase = 0.0) {
            FSI_Control C = new FSI_Control();

            //const double BaseSize = 1.0;
            //const double MeshFactor = 0.43;

            // basic database options
            // ======================



            C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\osciCylinder";
            C.savetodb = true;
            C.saveperiod = 10;
            C.ProjectName = "Uhlmann_k2_Re185_dt01_CellAgglo02_penalty4_SinglePhase";
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

            C.GridFunc = delegate {

                // Box1
                var box1_p1 = new double[2] { -6.17, -13.3 };
                var box1_p2 = new double[2] { 20.5, 13.3 };
                var box1 = new GridCommons.GridBox(box1_p1, box1_p2, 15, 15); //k1: ; k2: 15,15; k3: 15,15

                // Box2
                var box2_p1 = new double[2] { -2, -4 };
                var box2_p2 = new double[2] { 10, 4 };
                var box2 = new GridCommons.GridBox(box2_p1, box2_p2, 28, 20); //k1: 35,25 ; k2: 28,20; k3: 21, 15


                // Box3
                var box3_p1 = new double[2] { -1.5, -2.5 };
                var box3_p2 = new double[2] { 6, 2.5 };
                var box3 = new GridCommons.GridBox(box3_p1, box3_p2, 64, 48); //k1: 105, 75 ; k2: 64,48; k3: 39, 36 



                var grd = Grid2D.HangingNodes2D(box1, box2, box3);

                grd.EdgeTagNames.Add(1, "Pressure_Outlet_lower");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_upper");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "Velocity_Inlet_left");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet_right");
                }


                grd.DefineEdgeTags(delegate (double[] X) {
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
            if (!xPeriodic) {
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

            C.particleRadius = 0.5;
            C.includeRotation = false;
            C.includeTranslation = false;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 1.0 / 185;


            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) =>0 };

            //C.BoundaryFunc = new Func<double, double>[] { (t) => 0, (t) => 0.25 * Math.PI * 2 * 0.166 * -Math.Sin(Math.PI * 2 * 0.195 * t) };

            Func<double, double> yLevSet = t => (0.2 * Math.Cos(Math.PI * 2 * 0.156 * t));
            Func<double[], double, double> phi = (X, t) => -(X[0]).Pow2() + -(X[1] - yLevSet(t)).Pow2() + C.particleRadius.Pow2();
            //Func<double[], double, double> phi = (X, t) => -(X[0]).Pow2() + -(X[1]).Pow2() + C.particleRadius.Pow2();
            //Func<double[], double, double> phi = (X, t) => -(X[0]).Pow2() + -(X[1]-1).Pow2() + radius.Pow2();
            C.MovementFunc = phi;

            Func<double, double> xVelocity = t => 0;
            Func<double, double> yVelocity = t => (0.2 * Math.PI * 2 * 0.156 * -Math.Sin(Math.PI * 2 * 0.156 * t));
            Func<double, double>[] particleTransVelocity = { xVelocity, yVelocity };
            Func<double, double>[] particleAnglVelocity = { xVelocity, xVelocity };

            C.transVelocityFunc = particleTransVelocity;
            C.anglVelocityFunc = particleAnglVelocity;
            //C.VelocityFunc[1] = t => 0;// (0.2 * Math.PI * 2 * 0.156 * -Math.Sin(Math.PI * 2 * 0.156 * t));

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

            //C.InitialValues_Evaluators.Add("VelocityY", delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];

            //    double R = Math.Sqrt((x + 2).Pow2() + y.Pow2());

            //    double yVel = 0;

            //    if (R < 0.75) {
            //        yVel = 1;
            //    }
            //    return yVel;
            //});

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
            C.LinearSolver.MaxKrylovDim = 20;
            C.LinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.LinearSolver.NoOfMultigridLevels = 0;

            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.1;
            C.dtFixed = dt;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 300;
            C.NoOfTimesteps = 3000;

            // haben fertig...
            // ===============

            return C;
        }

        public static FSI_Control[] ParticleUnderGravity(string _DbPath = null, int k = 2, double VelXBase = 0.0) {
            List<FSI_Control> R = new List<FSI_Control>();

            // foreach (int i in new int[] {1,2, 3 }) {

            FSI_Control C = new FSI_Control();

            const double BaseSize = 1.0;

            //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
             //                   new Tuple<string,object>("k", k),
              //              };

            // k = i;

            // basic database options
            // ======================

            //C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Bug";
            C.savetodb = false;
            C.saveperiod = 1;
            C.ProjectName = "ParticleUnderGravity_k" + k + "_CellAgglo02_penalty4";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


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
            C.FieldOptions.Add("Curvature", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                switch (k) {
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
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


                grd.DefineEdgeTags(delegate (double[] X) {
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
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = false;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleRadius = 0.125;
            //C.particleMass = 1;

            

            C.Particles.Add(new Particle_Sphere(2, new double[] { 0.0, 1.0 }) {
                radius_P = (0.125/2.0),
                particleDensity = 1.25
            });

            //Func<double[], double, double> phiComplete = (X, t) => -1 * (C.Particles[0].phi_P(X, t));

            //Func<double, double> yLevSet = t => (t * t);
            //C.initialPos[0] = new double[] { 0.0, 4.0 };
            //Func<double[], double, double> phi = (X, t) => -(X[0] - C.initialPos[0][0]).Pow2() + -(X[1] - C.initialPos[0][1]).Pow2() + C.particleRadius.Pow2();
            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;

            //C.InitialValues_Evaluators.Add("Phi", X => phi(X, 0));
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);

            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

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
            C.LinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.NoOfMultigridLevels = 1;

            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = IBM_Solver.IBM_Control.TimesteppingScheme.BDF2;
            double dt = 0.0005;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1;
            C.NoOfTimesteps = 10000;

            // haben fertig...
            // ===============

            R.Add(C);
            //}

            return R.ToArray();
        }

        public static FSI_Control ParticleCollision(string _DbPath = null, int k = 2, double VelXBase = 0.0) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;

            //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
            //                    new Tuple<string,object>("k", k),
            //                };

            // k = i;

            // basic database options
            // ======================

            C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\Collisions";
            C.savetodb = false;
            C.saveperiod = 1;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


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

                int q = new int();
                int r = new int();

                switch (99) {
                    /*
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;
                        */
                    case 99:
                        q = 21;
                        r = 31;
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-1 * BaseSize, 1 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(0 * BaseSize, 3 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (0 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-3 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };




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
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleRadius = 0.25;
            //C.particleMass = 1;


            
            C.Particles.Add(new Particle_Ellipsoid(new double[] { 0.4, 1.0 }) {
                particleDensity = 1.0
            });

            C.Particles.Add(new Particle_Sphere(2, new double[] { 0.2, 0.5 }) {
                radius_P = 0.2,
                particleDensity = 1.0,        
            });
            C.Particles[1].transVelocityAtTimestep[0][0] = 0.5;
            C.Particles[1].transVelocityAtTimestep[0][1] = 1.0;

            C.Particles.Add(new Particle_Sphere(2, new double[] { 0.5, 2.0 }) {
                radius_P = 0.2,
                particleDensity = 1.0,
            });
            //C.Particles[2].transVelocityAtTimestep[0][0] = -1.5;
            //C.Particles[2].transVelocityAtTimestep[0][1] = -0.5;



            //C.Particles[1].transVelocityAtTimestep[1][1] = 1.0;
            //C.Particles[1].transVelocityAtTimestep[2][1] = 1.0;

            //C.Particles.Add(new Particle(2, 4, new double[] { -0.25, 1.0 }) {
            //    radius_P = 0.25,
            //    particleDensity = 0.75
            //});


            //Func<double[], double, double> phiComplete = (X, t) => 1 * (C.Particles[0].phi_P(X, t)* C.Particles[1].phi_P(X, t)* C.Particles[2].phi_P(X, t));

            //for (int i = 0;i<C.Particles.Count; i++) {
            //    phiComplete = (X,t) => phiComplete(X,t)*C.Particles[i].phi_P(X,t);
            //}


            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;         

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

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
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }

        /// <summary>
        /// Dry Collisions between particles of different shape
        /// </summary>
        public static FSI_Control CollisionDiss(string _DbPath = null, double VelXBase = 0.0) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;

            //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
            //                    new Tuple<string,object>("k", k),
            //                };

            // k = i;

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.saveperiod = 1;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


            // DG degrees
            // ==========

            C.SetDGdegree(1);
            C.FieldOptions["PhiDG"].Degree = 4;
            C.FieldOptions["Phi"].Degree = 4;

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                switch (99) {
                    /*
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;
                        */
                    case 99:
                        q = 31;
                        r = 46;
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(-2 * BaseSize, 2 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-2 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-2 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };




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
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.1;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleRadius = 0.25;
            //C.particleMass = 1;


            
            C.Particles.Add(new Particle_Ellipsoid(new double[] { -1.2, 0.9 }, startAngl: 90.0) {
                particleDensity = 1.0,
                thickness_P = 0.05,
                length_P = 0.1
            });
            C.Particles[0].transVelocityAtTimestep[0][0] = -5.0;
            C.Particles[0].rotationalVelocityAtTimestep[0] = -0.1;

            C.Particles.Add(new Particle_Sphere(2, new double[] { -0.6, 0.3},startAngl:-90.0) {
                radius_P = 0.25,
                particleDensity = 1.0,
            });
            C.Particles[1].transVelocityAtTimestep[0] = new double[2] { -5.0, 5.0 };
            C.Particles[1].rotationalVelocityAtTimestep[0] = -10;

 
            /*
            C.Particles.Add(new Particle_Hippopede(2, new double[] { -0.2, -0.5 }, startAngl:-45) {
                radius_P = 0.15,
                particleDensity = 1.0,
            });

            C.Particles[2].transVelocityAtTimestep[0] = new double[2] { -5.0,0.0};

            C.Particles.Add(new Particle_Squircle(2, new double[] { 1.0, 1.0 }, startAngl: -20.0) {
                radius_P = 0.25,
                particleDensity = 1.0,
            });
            C.Particles[3].transVelocityAtTimestep[0] = new double[2] { -5.0, -5.0 };

            C.Particles.Add(new Particle_Bean(2, new double[] { 1.0, -1.0 }, startAngl: -20.0) {
                radius_P = 0.25,
                particleDensity = 1.0,
            });
            
*/
            C.pureDryCollisions = true;
            C.collisionModel = FSI_Control.CollisionModel.MomentumConservation;



            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            //Func<double[], double, double> phiComplete = (X, t) => (1 * (C.Particles[0].phi_P(X, t) * C.Particles[1].phi_P(X, t)) * (C.Particles[2].phi_P(X, t) * C.Particles[3].phi_P(X, t))* C.Particles[4].phi_P(X, t));

            //for (int i = 0;i<C.Particles.Count; i++) {
            //    phiComplete = (X,t) => phiComplete(X,t)*C.Particles[i].phi_P(X,t);
            //}


            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;         

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

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
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.000005;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1.2;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }

        public static FSI_Control FiveRandomParticles(string _DbPath = null, int k = 2, double dt = 0.001, double VelXBase = 0.0, int collisionModelInt = 1) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;

            //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
            //                    new Tuple<string,object>("k", k),
            //                };

            // k = i;

            C.collisionModel = (FSI_Solver.FSI_Control.CollisionModel)collisionModelInt;

            // basic database options
            // ======================     
            
            C.DbPath = @"\\hpccluster\hpccluster-scratch\krause\cluster_db";
            C.savetodb = true;
            C.saveperiod = (int)(0.01/dt);
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


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

                int q = new int();
                int r = new int();

                switch (99) {
                    /*
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;
                        */
                    case 99:
                        q = 21; //21/31
                        r = 81; //81/121
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-1.0 * BaseSize, 1.0 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(0.0 * BaseSize, 8 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.0 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.0 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (0.0 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-8.0 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };




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
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.01;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleRadius = 0.25;
            //C.particleMass = 1;


            C.Particles.Add(new Particle_Sphere(2, new double[] { -0.2, 7.5 }, startAngl: 45.0) {
                radius_P = 0.1,
                particleDensity = 3.0
            });

            C.Particles.Add(new Particle_Ellipsoid(new double[] { 0.2, 7.3 }, startAngl: 30.0) {
                particleDensity = 3.0,
                thickness_P = 0.05,
                length_P = 0.1
            });

            C.Particles.Add(new Particle_Squircle(2, new double[] { -0.2, 6.95 }, startAngl: -20.0) {
                radius_P = 0.1,
                particleDensity = 3.0,
            });

            C.Particles.Add(new Particle_Sphere(2, new double[] { -0.5, 7.2 }, startAngl: -45.0) {
                radius_P = 0.15,
                particleDensity = 3.0,
            });

            C.Particles.Add(new Particle_Squircle(2, new double[] { 0.2, 6.5 }, startAngl: -45.0) {
                radius_P = 0.15,
                particleDensity = 3.0,
            });


            C.collisionModel = FSI_Control.CollisionModel.MomentumConservation;

            C.pureDryCollisions = false;

            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            ////Define level-set
            //Func<double[], double, double> phiComplete = delegate (double[] X, double t) {
            //    int exp = C.Particles.Count - 1;
            //    double ret = Math.Pow(-1, exp);
            //    for (int i = 0; i < C.Particles.Count; i++) {
            //        ret *= C.Particles[i].phi_P(X, t);
            //    }
            //    return ret;
            //};


            //Func<double[], double, double> phiComplete = (X, t) => (1 * (C.Particles[0].phi_P(X, t) * C.Particles[1].phi_P(X, t) * C.Particles[2].phi_P(X, t)* C.Particles[3].phi_P(X, t) * C.Particles[4].phi_P(X, t)));

            //for (int i = 0;i<C.Particles.Count; i++) {
            //    phiComplete = (X,t) => phiComplete(X,t)*C.Particles[i].phi_P(X,t);
            //}


            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;         

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("42c82f3c-bdf1-4531-8472-b65feb713326"), 400);
            //C.GridGuid = new Guid("f1659eb6-b249-47dc-9384-7ee9452d05df");


            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = true;


            // misc. solver options
            // ====================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            //double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }

        public static FSI_Control FallingEllipse(string _DbPath = null, int k = 2, double VelXBase = 0.0, double angle= 0.0) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;


            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.saveperiod = 100;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


            // DG degrees
            // ==========

            //C.FieldOptions.Add("VelocityX", new FieldOpts() {
            //    Degree = k,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("VelocityY", new FieldOpts() {
            //    Degree = k,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("Pressure", new FieldOpts() {
            //    Degree = k - 1,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("PhiDG", new FieldOpts() {
            //    Degree = 2,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("Phi", new FieldOpts() {
            //    Degree = 2,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            C.SetDGdegree(k);

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                switch (99) {
                    /*
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;
*/
                    case 99:
                        q = 61; //61 //31
                        r = 81; //81 //41
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(-2 * BaseSize, 2 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-2 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-2 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };

            C.GridPartType = GridPartType.Hilbert;


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
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleMass = 1;


            C.Particles.Add(new Particle_Ellipsoid(new double[] { 0.0*BaseSize, 1.0*BaseSize }, startAngl: angle) {
                particleDensity = 10.0,
                length_P = 0.1*BaseSize,
                thickness_P = 0.2*BaseSize,
                gravityVertical = 9.81
            });

            //C.Particles[0].rotationalVelocityAtTimestep[0] = 10;

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            //Func<double[], double, double> phiComplete = delegate (double[] X, double t) {
            //    double r = 1 * (C.Particles[0].phi_P(X, t));
            //    if (double.IsNaN(r) || double.IsInfinity(r))
            //        throw new ArithmeticException();
            //    return r;
            //};

            //for (int i = 0;i<C.Particles.Count; i++) {
            //    phiComplete = (X,t) => phiComplete(X,t)*C.Particles[i].phi_P(X,t);
            //}


            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;         

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

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
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 2.0;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }

        public static FSI_Control DraftKissingTumbling(string _DbPath = null, int k = 2, double VelXBase = 0.0,int collisionModelInt = 0) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;

            //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
            //                    new Tuple<string,object>("k", k),
            //                };

            // k = i;

            C.collisionModel = (FSI_Solver.FSI_Control.CollisionModel)collisionModelInt;

            // basic database options
            // ======================

            //C.DbPath = @"\\hpccluster\hpccluster-scratch\krause\DraftKissing_db";
            C.savetodb = false;
            C.saveperiod = 10;
            C.ProjectName = "DraftKissingTumbling";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


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
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {

                int q = new int();
                int r = new int();

                switch (99) {
                    /*
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;
                    */
                    case 99:
                        q = 61; //31/41
                        r = 241; //121/161
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-1.0 * BaseSize, 1.0 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(0.0 * BaseSize, 8 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.0 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.0 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (0.0 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-8.0 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };           
          
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
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.mu_A = 0.01;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleMass = 1;


            C.Particles.Add(new Particle_Sphere(2, new double[] { 0.0, 7.2 }) {
                radius_P = 0.1,
                particleDensity = 1.01
            });

            //C.Particles[0].transVelocityAtTimestep[0][1] = -0.5;

            C.Particles.Add(new Particle_Sphere(2, new double[] { 0.0, 6.8 }) {
                radius_P = 0.1,
                particleDensity = 1.01,
            });

            //C.Particles[1].transVelocityAtTimestep[0][1] = -0.5;

            //Func<double[], double, double> phiComplete = (X, t) => -1 * (C.Particles[0].phi_P(X, t) * C.Particles[1].phi_P(X, t));
      

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);


            // For restart
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("8e89aedc-f1bf-44b0-8e9a-8d028057b4ce"), 900);
            //C.GridGuid = new Guid("bdacff75-dd60-400c-8327-24aabff4495f");


            // Physical Parameters
            // ===================

            C.PhysicalParameters.IncludeConvection = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.LevelSetSmoothing = false;
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;
            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }

        public static FSI_Control WallCollisionTest(string _DbPath = null, int k = 2, double VelXBase = 0.0, double angle = 0.0, int collisionModelInt = 0) {
            FSI_Control C = new FSI_Control();


            const double BaseSize = 1.0;

            C.collisionModel = (FSI_Solver.FSI_Control.CollisionModel)collisionModelInt;


            // basic database options
            // ======================

            //C.DbPath = @"\\hpccluster\hpccluster-scratch\krause\cluster_db";
            //C.DbPath = @"\\dc1\userspace\krause\BoSSS_DBs\FallingEllipse";
            C.savetodb = false;
            C.saveperiod = 100;
            C.ProjectName = "ParticleCollisionTest";
            C.ProjectDescription = "Gravity";
            C.SessionName = C.ProjectName;
            C.Tags.Add("with immersed boundary method");
            C.AdaptiveMeshRefinement = true;


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

                int q = new int();
                int r = new int();
                int iCase = 99; // this construction prevents compile warning 

                switch (iCase) {
                    case 1:
                        q = 60;
                        r = 178;
                        break;

                    case 2:
                        q = 41;
                        r = 121;
                        break;

                    case 3:
                        q = 31;
                        r = 91;
                        break;
                        
                    case 99:
                        q = 31; //61 //31
                        r = 31; //81 //41
                        break;

                    default:

                        throw new ApplicationException();
                }

                //q = 16;
                //r = 46;

                double[] Xnodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, q); //k1: 71; k2:41; k3: 31
                double[] Ynodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, r); //k1: 211; k2:121; k3: 91

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet_left");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(3, "Wall_lower");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] + (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 4;

                    Debug.Assert(et != 0);
                    return et;
                });

                Console.WriteLine("Cells:" + grd.NumberOfCells);

                return grd;
            };




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
            //C.LevelSetMovement = "coupled";
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.includeRotation = true;
            C.includeTranslation = true;

            // Fluid Properties
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 0.1;

            // Particle Properties
            //C.PhysicalParameters.mu_B = 0.1;
            //C.particleMass = 1;

            C.Particles.Add(new Particle_Sphere(2, new double[] { -0.5, -1.35 }, startAngl: 0.0) {
                radius_P = 0.1,
                particleDensity = 1.25,
            });

            C.Particles.Add(new Particle_Sphere(2, new double[] { 0.8, -1.35 }, startAngl: 0.0) {
                radius_P = 0.1,
                particleDensity = 1.25,
            });

            //C.Particles[0].rotationalVelocityAtTimestep[0] = 10;

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            //Func<double[], double, double> phiComplete = (X, t) => -1 * (C.Particles[0].phi_P(X, t)* C.Particles[1].phi_P(X, t));

            //for (int i = 0;i<C.Particles.Count; i++) {
            //    phiComplete = (X,t) => phiComplete(X,t)*C.Particles[i].phi_P(X,t);
            //}


            //Func<double[], double, double> phi = (X, t) => -(X[0] - t+X[1]);
            //C.MovementFunc = phi;         

            //C.InitialValues_Evaluators.Add("Phi", X => phiComplete(X, 0));
            //C.InitialValues_Evaluators.Add("Phi", X => -1);
            //C.InitialValues.Add("VelocityX#B", X => 1);
            C.InitialValues_Evaluators.Add("VelocityX", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0);
            //C.InitialValues.Add("Phi", X => -1);
            //C.InitialValues.Add("Phi", X => (X[0] - 0.41));

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
            C.LinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.NoOfMultigridLevels = 1;


            // Timestepping
            // ============

            //C.Timestepper_Mode = FSI_Control.TimesteppingMode.Splitting;
            C.Timestepper_Scheme = FSI_Solver.FSI_Control.TimesteppingScheme.BDF2;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 2.0;
            C.NoOfTimesteps = 100000;

            // haben fertig...
            // ===============

            return C;

        }
    }


}
