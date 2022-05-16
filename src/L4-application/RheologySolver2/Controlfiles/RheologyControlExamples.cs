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
using System.Text;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation;
using BoSSS.Solution.Tecplot;
using System.IO;
using BoSSS.Solution.GridImport;
using System.Configuration;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Control File for specified examples for calculation with viscoelastic fluid
    /// </summary>
    static public class RheologyControlExamples {

        /// <summary>
        /// 4:1 Contraction Flow
        /// </summary>
        static public RheologyControl Contraction(string path = null, int degree = 2, int GridLevel = 3) { //int kelem = 4
            RheologyControl C = new RheologyControl();

            //Path für cluster
            //\\dc1\userspace\kikker\cluster\cluster_db\ContractionNYC

            //Path für lokale DB
            //C:\AnnesBoSSSdb\Contraction

            //Solver Options          
            C.savetodb = false;
            C.DbPath = path;
            C.ProjectName = "Contraction";

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1E-7;

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.useFDJacobianForOperatorMatrix = true;

            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.NoOfTimesteps = 1;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;

            C.ObjectiveParam = 1.0;

            //Debugging and Solver Analysis
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = true;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = false;
            //C.UsePerssonSensor = true;
            //C.SensorLimit = 1e-4;
            C.AdaptiveMeshRefinement = true;
            C.RefinementLevel = 4;
            //C.AMR_startUpSweeps = 1;
            //C.UseArtificialDiffusion = false;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.beta = 1;// 0.11;
            C.Reynolds = 1;
            C.Weissenberg = 1.0;
            C.RaiseWeissenberg = true;
            C.giesekusfactor = 0.0;

            //Grid Params
            double L = 20;
            double H = 2;
            int cellsX = 40;
            int cellsY = 8;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 1;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 0;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;


            //Exact Solution Contraction

            // Set Initial Conditions
            Func<double[], double, double> VelocityXfunction = (X, t) => (1.0 - (X[1] * X[1]) / (2 * H));
            Func<double[], double, double> VelocityYfunction = (X, t) => (0.0);
            Func<double[], double, double> Pressurefunction = (X, t) => 2 / (2 * H) * C.Reynolds * (L - X[0]);
            Func<double[], double, double> StressXXfunction = (X, t) => 2 * C.Weissenberg * (1 - C.beta) * ((-2 / (2 * H) * X[1]) * (-2 / (2 * H) * X[1])); //aim Weissenberg!!!
            Func<double[], double, double> StressXYfunction = (X, t) => (1 - C.beta) * (-2 / (2 * H) * X[1]);
            Func<double[], double, double> StressYYfunction = (X, t) => (0.0);

            // Insert Exact Solution
            //C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            //C.ExSol_Pressure = Pressurefunction;
            //C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            // Create Fields
            C.SetDGdegree(degree);

            // Create Grid

            double[] pt1a = new double[] { L / 2, H / 4 };
            double[] pt1b = new double[] { L, H };

            BoundingBox boundingBox1;
            boundingBox1 = new BoundingBox(pt1a, pt1b);

            double[] pt2a = new double[] { L / 2, -H / 4 };
            double[] pt2b = new double[] { L, -H };

            BoundingBox boundingBox2;
            boundingBox2 = new BoundingBox(pt2a, pt2b);

            BoundingBox[] BoundingBox = new BoundingBox[] { boundingBox1, boundingBox2 };

            C.GridFunc = delegate {
                // UNIFORM CARTESIAN GRID
                var _xNodes = GenericBlas.Linspace(0, L, cellsX + 1);
                var _yNodes = GenericBlas.Linspace(-H, H, (cellsY + 1));
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC, false, null, boundingBox1, boundingBox2);
                //var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_9, C.FixedStreamwisePeriodicBC, false, null, boundingBox1, boundingBox2);

                #region grids
                //var _xNodes = GenericBlas.Linspace(0, 10, cells2 + 1);// 10 * GridLevel + 1); //(10 * kelem + 1));
                //var _yNodes = GenericBlas.Linspace(-2, 2, (cells2 / 4) + 1);// (int)(2 * 1.5 * GridLevel) + GridLevel + 1); //(int)((2 * 1.5 * kelem) + kelem + 1));
                //var grd1 = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                //var _xNodes2 = GenericBlas.Linspace(10, 20, cells2 + 1);// 10 * GridLevel + 1); //(10 * kelem + 1));
                //var _yNodes2 = GenericBlas.Linspace(-0.5, 0.5, (cells2 / 4) + 1);// (int)(2 * 1.5 * GridLevel) + GridLevel + 1); //(int)((2 * 1.5 * kelem) + kelem + 1));
                //var grd2 = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                //var grdM = GridCommons.MergeLogically(grd1, grd2);
                //var grd = GridCommons.Seal(grdM);

                // NON_UNIFORM CARTESIAN GRID
                //var _xNodes1 = Grid1D.TanhSpacing(0, 10, (cells2 / 4) + 1, 2, false);
                //_xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
                //var _xNodes2 = Grid1D.TanhSpacing(10, 20, (cells2 / 4) + 1, 2, true);
                //var _xNodes = ArrayTools.Cat(_xNodes1, _xNodes2);

                //var _yNodes1 = Grid1D.TanhSpacing(0, 0.5, (cells2 / 6) + 1, 1.5, false);
                //_yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                //var _yNodes2 = Grid1D.TanhSpacing(0.5, 2, (6 * cells2 / 16) + 1, 1.5, true);
                ////var _yNodes = ArrayTools.Cat(_yNodes1, _yNodes2);

                //var _yNodes3 = Grid1D.TanhSpacing(-0.5, 0, (cells2 / 6) + 1, 1.5, true);
                //_yNodes3 = _yNodes3.GetSubVector(0, (_yNodes3.Length - 1));
                //var _yNodes4 = Grid1D.TanhSpacing(-2, -0.5, (6 * cells2 / 16) + 1, 1.5, false);
                //_yNodes4 = _yNodes4.GetSubVector(0, (_yNodes4.Length - 1));
                //var _yNodes = ArrayTools.Cat(_yNodes4, _yNodes3, _yNodes1, _yNodes2);



                //// CARTESIAN GRID WITH HANGING NODES REFINEMENT
                //double[] ecke = new double[] { 10, 0.5 };

                //var boxA1_p1 = new double[2] { 0, -2 };
                //var boxA1_p2 = new double[2] { 10, 2 };
                //var boxA1 = new GridCommons.GridBox(boxA1_p1, boxA1_p2, 2 * 40, 2 * 16);

                //var boxA2_p1 = new double[2] { 9, -1 };
                //var boxA2_p2 = new double[2] { 10, 1 };
                ////var boxA2_p1 = new double[2] { ecke[0] - 0.5, ecke[1] - 0.25 };
                ////var boxA2_p2 = new double[2] { ecke[0], ecke[1] + 0.25 };
                //var boxA2 = new GridCommons.GridBox(boxA2_p1, boxA2_p2, 2 * 8, 2 * 16);

                //var boxA3_p1 = new double[2] { 9.5, -0.75 };
                //var boxA3_p2 = new double[2] { 10, 0.75 };
                //var boxA3 = new GridCommons.GridBox(boxA3_p1, boxA3_p2, 2 * 8, 2 * 24);

                //var grdA = Grid2D.HangingNodes2D(boxA1, boxA2, boxA3);

                //var boxB1_p1 = new double[2] { 10, -0.5 };
                //var boxB1_p2 = new double[2] { 20, 0.5 };
                //var boxB1 = new GridCommons.GridBox(boxB1_p1, boxB1_p2, 2 * 40, 2 * 4);

                //var boxB2_p1 = new double[2] { 10, -0.5 };
                //var boxB2_p2 = new double[2] { 11, 0.5 };
                //var boxB2 = new GridCommons.GridBox(boxB2_p1, boxB2_p2, 2 * 8, 2 * 8);

                //var boxB3_p1 = new double[2] { 10, -0.5 };
                //var boxB3_p2 = new double[2] { 10.5, 0.5 };
                //var boxB3 = new GridCommons.GridBox(boxB3_p1, boxB3_p2, 2 * 8, 2 * 16);

                //var grdB = Grid2D.HangingNodes2D(boxB1, boxB2, boxB3);

                //var grdM = GridCommons.MergeLogically(grdA, grdB);
                //var grd = GridCommons.Seal(grdM);

                // COARSE CARTESIAN GRID WITH HANGING NODES REFINEMENT - FOR DEBUGGING!
                //double[] ecke = new double[] { 10, 0.5 };

                //var boxA1_p1 = new double[2] { 0, -2 };
                //var boxA1_p2 = new double[2] { 10, 2 };
                //var boxA1 = new GridCommons.GridBox(boxA1_p1, boxA1_p2, 2 * 10, 2 * 4);

                //var boxA2_p1 = new double[2] { 9, -1 };
                //var boxA2_p2 = new double[2] { 10, 1 };
                ////var boxA2_p1 = new double[2] { ecke[0] - 0.5, ecke[1] - 0.25 };
                ////var boxA2_p2 = new double[2] { ecke[0], ecke[1] + 0.25 };
                //var boxA2 = new GridCommons.GridBox(boxA2_p1, boxA2_p2, 2 * 2, 2 * 4);

                ////var boxA3_p1 = new double[2] { 9.5, -0.75 };
                ////var boxA3_p2 = new double[2] { 10, 0.75 };
                ////var boxA3 = new GridCommons.GridBox(boxA3_p1, boxA3_p2, 2 * 4, 2 * 12);

                //var grdA = Grid2D.HangingNodes2D(boxA1, boxA2);

                //var boxB1_p1 = new double[2] { 10, -0.5 };
                //var boxB1_p2 = new double[2] { 20, 0.5 };
                //var boxB1 = new GridCommons.GridBox(boxB1_p1, boxB1_p2, 2 * 10, 2 * 1);

                //var boxB2_p1 = new double[2] { 10, -0.5 };
                //var boxB2_p2 = new double[2] { 11, 0.5 };
                //var boxB2 = new GridCommons.GridBox(boxB2_p1, boxB2_p2, 2 * 2, 2 * 2);

                ////var boxB3_p1 = new double[2] { 10, -0.5 };
                ////var boxB3_p2 = new double[2] { 10.5, 0.5 };
                ////var boxB3 = new GridCommons.GridBox(boxB3_p1, boxB3_p2, 2 * 4, 2 * 8);

                //var grdB = Grid2D.HangingNodes2D(boxB1, boxB2);

                //var grdM = GridCommons.MergeLogically(grdA, grdB);
                //var grd = GridCommons.Seal(grdM);
                #endregion

                if (!C.FixedStreamwisePeriodicBC) {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                }
                               

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (!C.FixedStreamwisePeriodicBC) {
                        if (Math.Abs(x - (0)) < 1.0e-6) {
                            //left
                            return "Velocity_inlet";
                        }

                        if (Math.Abs(x - (L)) < 1.0e-6) {
                            //right
                            return "Pressure_Outlet";
                        }
                    }

                    if (Math.Abs(y - (-H)) < 1.0e-6 && x < L / 2 + 1.0e-6) {
                        //bottom front
                        return "Wall_bottom";
                    }

                    //if (Math.Abs(y - (0)) < 1.0e-6)
                    //{
                    //    //symmetry line
                    //    return "Wall_bottom";
                    //}

                    if (Math.Abs(y - (+H)) < 1.0e-6 && x < L / 2 + 1.0e-6) {
                        //top front
                        return "Wall_top";
                    }

                    if (Math.Abs(y - (-H / 4)) < 1.0e-6 && x > L / 2 - 1.0e-6) {
                        // bottom back
                        return "Wall_bottom";
                    }

                    if (Math.Abs(y - (H / 4)) < 1.0e-6 && x > L / 2 - 1.0e-6) {
                        // top back
                        return "Wall_top";
                    }

                    if (Math.Abs(x - (L / 2)) < 1.0e-6 && y < -H / 4 - 1.0e-6) {
                        // bottom contraction
                        return "Wall_Contraction_bottom";
                    }

                    if (Math.Abs(x - (L / 2)) < 1.0e-6 && y > H / 4 - 1.0e-6) {
                        //top contraction
                        return "Wall_Contraction_top";
                    }

                    throw new ArgumentOutOfRangeException("at x = " + x + "and y = " + y);
                });


                return grd;
            };


            // Analytical Sol for Params
            //if (C.SetParamsAnalyticalSol == true) {
            //    C.VelFunctionU = X => VelocityXfunction(X, 0);
            //    C.VelFunctionV = X => VelocityYfunction(X, 0);
            //    C.PresFunction = X => Pressurefunction(X, 0);
            //}

            // Set Initial Conditions
            if (C.SetInitialConditions == true) {

                C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXX", X => StressXXfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));

                if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true) {
                    C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));
                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            C.AddBoundaryValue("Wall_bottom");
            C.AddBoundaryValue("Wall_top");
            //C.AddBoundaryValue("Wall_Contraction_bottom");
            //C.AddBoundaryValue("Wall_Contraction_top");
            //C.AddBoundaryCondition("FreeSlip");


            if (!C.FixedStreamwisePeriodicBC) {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", VelocityXfunction);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", VelocityYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", StressXXfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", StressXYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", StressYYfunction);
                //C.AddBoundaryCondition("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Pressure_Outlet");

            }
            return C;
        }
        //__________________________________________________________________________________________________________________

        /// <summary>
        /// Confined cylinder in a channel flow
        /// </summary>
        static public RheologyControl ConfinedCylinder(
            string path = null,
            int degree = 1) {

            RheologyControl C = new RheologyControl();

                                          
            //Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            C.DbPath = path;
            C.ProjectName = "Cylinder";

            C.AlternateDbPaths = new[] {
                (@"\\dc1\userspace\kikker\cluster\cluster_db\ConfinedCylinder_Drag", "hpccluster"), // Anne's Account
                (@"d:\Users\kummer\default_bosss_db", "terminal03"),     // Florian's account on Terminal03
                (@"c:\Users\florian\default_bosss_db", "rennmaschin"),   // Florian's old Laptop
                (@"c:\Users\flori\default_bosss_db", "stormbreaker")     // Florian's new Laptop
            };

            C.savetodb = true;
            C.DbPath = path;
            C.ProjectName = "Cylinder";
            C.SessionName = "with continuation";

            // stationary calculation
            // ======================
            
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            

            /*
            // transient calculation
            // =====================
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            C.dt = 0.1;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            */
            // linear solver config 
            // ====================
            //C.LinearSolver = LinearSolverCode.classic_mumps.GetConfig();
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();//.exp_Kcycle_schwarz_4Rheology;
          

            // nonlinear solver config
            // =======================
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1E-5;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();


            // more solver options
            // ===================
            
            C.ObjectiveParam = 1.0;
            C.useFDJacobianForOperatorMatrix = false;

            //C.UsePerssonSensor = false;
            //C.SensorLimit = 1e-4;

            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 10;

            //C.UseArtificialDiffusion = false;

            C.Bodyforces = true;

            // Homotopy:
            C.RaiseWeissenberg = false;
            //C.WeissenbergIncrement = 0.1;

            //Debugging and Solver Analysis
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = false;

            //Physical Params
            // ==============
            double u0 = 1.5; // 0.375;// 0.66;// 3 / 2;
            double h = 4;

            C.Stokes = false;
            C.StokesConvection = false;
            C.FixedStreamwisePeriodicBC = false;
            C.beta = 0.59;
            C.Reynolds = 0.1;
            C.Weissenberg = 0.00000; //aim Weissenberg number!
            C.giesekusfactor = 0.0;

            //Penalties
            C.ViscousPenaltyScaling = 1.0;
            C.Penalty2 = 1;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 1;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;
            C.alpha = 1;
            C.StressPenalty = 1.0;

            //Exact Solution Confined Cylinder

            // Set Initial Conditions / Boundary Conditions   
            
            Func<double[], double, double> VelocityXfunction = delegate (double[] X, double t) {
                return u0 * (1 - (X[1] * X[1]) / h);
            };
            Func<double[], double, double> VelocityYfunction = (X, t) => 0.0;
            Func<double[], double, double> Pressurefunction = (X, t) => u0 * 0.5 * C.Reynolds * (35 - X[0]);
            Func<double[], double, double> StressXXfunction = (X, t) =>  2 * C.Weissenberg * (1 - C.beta) * u0 * (-2 / h) * X[1] * u0 * (-2 / h) * X[1];
            Func<double[], double, double> StressXYfunction = (X, t) => (1 - C.beta) * u0 * (-2 / h) * X[1];
            Func<double[], double, double> StressYYfunction = (X, t) => (0.0);

            // Insert Exact Solution
            //C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            //C.ExSol_Pressure = Pressurefunction;
            //C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            C.SetDGdegree(degree);

            // Create Grid
            // grids used by florian
            string grid = "445bf1e8-5082-4100-83db-a53f2d2aeb97"; // Florian new laptop (full, level 0); Cell Count = 640; Dim = 2 }
            //string grid = "b87b95c1-10d8-4b99-bc82-9c5fec70219d"; // Florian new laptop (full, level 1); Cell Count = 2416; Dim = 2 }
            //string grid = "3cfced02-dbe6-4449-a3d0-d9e258813608"; // Florian new laptop (full, level 2); Cell Count = 8504; Dim = 2 }
            //string grid = "70c0393d-62b9-4c32-9c86-1173251617da"; // Florian new laptop (full, level 3); Cell Count = 29424; Dim = 2 }
            //string grid = "c88c914b-c387-4894-9697-a78bad31f2da"; // florian terminal03 (full, level 0)
            //string grid = "061e7cfb-7ffe-4540-bc74-bfffce824fef"; // florian terminal03 (full, level 1)
            //string grid = "51aadb49-e3d5-4e88-897e-13b6b329995b"; // florian terminal03 (full, level 2)
            //string grid = "ab9b8f9a-e1aa-469d-8eea-b56a88e672a4"; // florian terminal03 (full, level 3)


            Guid gridGuid;
            if(Guid.TryParse(grid, out gridGuid)) {
                C.GridGuid = gridGuid;
            } else {
                C.GridFunc = delegate () {
                    GridCommons _grid;

                    _grid = GridImporter.Import(grid);

                    return _grid;
                };
            }
            

            // Analytical Sol for Params
            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
                C.PresFunction = X => Pressurefunction(X, 0);
            }

            // restart (florian, terminal03)
            //Guid restartID = new Guid("45c813f2-8be5-43ab-9e41-7abbca99cc99"); 
            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, new TimestepNumber(1, 4)); // Weissenberg 0.4, deg = 2
            
            

            if (C.RestartInfo == null) {

                C.SetInitialConditions = false;

                //Set Initial Conditions
                if (C.SetInitialConditions == true)
                {

                    C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                    C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                    C.InitialValues_Evaluators.Add("StressXX", X => 0);// StressXXfunction(X, 0));
                    C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                    C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));

                    if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true)
                    {
                        C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));
                    }
                }
            }

            // Set Boundary Conditions
            //C.AddBoundaryValue("Wall_bottom", "VelocityX", X => 0);
            C.AddBoundaryValue("Wall_top", "VelocityX", "X => 0", TimeDependent:false);
            //C.AddBoundaryValue("Wall_bottom", "VelocityY", X => 0);
            //C.AddBoundaryValue("Wall_top", "VelocityY", X => 0);
            //C.AddBoundaryValue("Wall_cylinder", "VelocityX", X => 0);
            //C.AddBoundaryValue("Wall_cylinder", "VelocityY", X => 0);
            //C.AddBoundaryValue("Freeslip");


            if(!C.FixedStreamwisePeriodicBC) {
                
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", $"X => {u0} * (1 - (X[1] * X[1]) / {h})", TimeDependent:false);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", $"X => 0.0", TimeDependent:false);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", $"X => 2 * {C.Weissenberg} * (1 - {C.beta}) * {u0} * (-2 / {h}) * X[1] * {u0} * (-2 / {h}) * X[1]", TimeDependent:false);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", $"X => (1 - {C.beta}) * {u0} * (-2 / {h}) * X[1]", TimeDependent:false);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", $"X => 0.0", TimeDependent:false);
                //C.AddBoundaryCondition("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Pressure_Outlet");

            }
            return C;
        }
        //__________________________________________________________________________________________________________________

    }
}