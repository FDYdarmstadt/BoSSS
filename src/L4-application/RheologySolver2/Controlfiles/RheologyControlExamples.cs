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
        static public RheologyControl Contraction(string path = @"\\dc1\userspace\kikker\cluster\cluster_db\ContractionNYC", int degree = 1, int GridLevel = 3) { //int kelem = 4
            RheologyControl C = new RheologyControl();

            //Path für cluster
            //\\dc1\userspace\kikker\cluster\cluster_db\ContractionNYC

            //Path für lokale DB
            //C:\AnnesBoSSSdb\Contraction

            //Solver Options
            C.NoOfTimesteps = 15;
            C.savetodb = true;
            C.DbPath = path;
            C.ProjectName = "Contration";
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1E-7;
            C.LinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.ConvergenceCriterion = 1E-7;

            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.NewtonGMRES;
            C.ObjectiveParam = 1.0;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = false;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.beta = 0.11;
            C.Reynolds = 1;
            C.Weissenberg = 0.1;
            C.RaiseWeissenberg = true;

            //Grid Params
            //double GridLevel = 5;
            double h = Math.Pow(2, -GridLevel + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 1;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 1;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;
            

            //Exact Solution Contraction

            // Set Initial Conditions
            Func<double[], double, double> VelocityXfunction = (X, t) => ((3/2) *4* (1 - (X[1] * X[1])/4)); //(4 - X[1] * X[1]);
            Func <double[], double, double> VelocityYfunction = (X, t) => (0.0);
            Func<double[], double, double> Pressurefunction = (X, t) => 2 / C.Reynolds * (20 -  X[0]);
            Func<double[], double, double> StressXXfunction = (X, t) => 2 * C.Weissenberg * (1 - C.beta) * (((-2 * X[1])) * ((-2 * X[1])));
            Func<double[], double, double> StressXYfunction = (X, t) => (1 - C.beta) * ((-2 * X[1]));
            Func<double[], double, double> StressYYfunction = (X, t) => (0.0);

            // Insert Exact Solution
            C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            C.ExSol_Pressure = Pressurefunction;
            C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            // Create Fields
            //int degree = 2;
            C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            // Create Grid

            double[] pt1a = new double[] { 10, 0.5 };
            double[] pt1b = new double[] { 20, 2 };

            BoundingBox boundingBox1;
            boundingBox1 = new BoundingBox(pt1a, pt1b);

            double[] pt2a = new double[] { 10, -0.5 };
            double[] pt2b = new double[] { 20, -2 };

            BoundingBox boundingBox2;
            boundingBox2 = new BoundingBox(pt2a, pt2b);

            C.GridFunc = delegate
            {
                // UNIFORM CARTESIAN GRID
                //var _xNodes = GenericBlas.Linspace(0, 20, cells2 + 1);// 10 * GridLevel + 1); //(10 * kelem + 1));
                //var _yNodes = GenericBlas.Linspace(-2, 2, (cells2 / 4) + 1);// (int)(2 * 1.5 * GridLevel) + GridLevel + 1); //(int)((2 * 1.5 * kelem) + kelem + 1));

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

                //var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC, false, boundingBox1, boundingBox2);

                //// CARTESIAN GRID WITH HANGING NODES REFINEMENT
                double[] ecke = new double[] { 10, 0.5 };

                var boxA1_p1 = new double[2] { 0, -2 };
                var boxA1_p2 = new double[2] { 10, 2 };
                var boxA1 = new GridCommons.GridBox(boxA1_p1, boxA1_p2, 2 * 40, 2 * 16);

                var boxA2_p1 = new double[2] { 9, -1 };
                var boxA2_p2 = new double[2] { 10, 1 };
                //var boxA2_p1 = new double[2] { ecke[0] - 0.5, ecke[1] - 0.25 };
                //var boxA2_p2 = new double[2] { ecke[0], ecke[1] + 0.25 };
                var boxA2 = new GridCommons.GridBox(boxA2_p1, boxA2_p2, 2 * 8, 2 * 16);

                var boxA3_p1 = new double[2] { 9.5, -0.75 };
                var boxA3_p2 = new double[2] { 10, 0.75 };
                var boxA3 = new GridCommons.GridBox(boxA3_p1, boxA3_p2, 2 * 8, 2 * 24);

                var grdA = Grid2D.HangingNodes2D(boxA1, boxA2, boxA3);

                var boxB1_p1 = new double[2] { 10, -0.5 };
                var boxB1_p2 = new double[2] { 20, 0.5 };
                var boxB1 = new GridCommons.GridBox(boxB1_p1, boxB1_p2, 2 * 40, 2 * 4);

                var boxB2_p1 = new double[2] { 10, -0.5 };
                var boxB2_p2 = new double[2] { 11, 0.5 };
                var boxB2 = new GridCommons.GridBox(boxB2_p1, boxB2_p2, 2 * 8, 2 * 8);

                var boxB3_p1 = new double[2] { 10, -0.5 };
                var boxB3_p2 = new double[2] { 10.5, 0.5 };
                var boxB3 = new GridCommons.GridBox(boxB3_p1, boxB3_p2, 2 * 8, 2 * 16);

                var grdB = Grid2D.HangingNodes2D(boxB1, boxB2, boxB3);

                var grdM = GridCommons.MergeLogically(grdA, grdB);
                var grd = GridCommons.Seal(grdM);

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

                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                }

                //grd.EdgeTagNames.Add(2, "FreeSlip");
                grd.EdgeTagNames.Add(2, "Wall_bottom");
                grd.EdgeTagNames.Add(3, "Wall_top");


                grd.EdgeTagNames.Add(5, "Wall_Contraction_bottom");
                grd.EdgeTagNames.Add(6, "Wall_Contraction_top");

                grd.DefineEdgeTags(delegate (double[] _X)
                {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (0)) < 1.0e-6)
                        {
                            //left
                            return 1;
                        }

                        if (Math.Abs(x - (20)) < 1.0e-6)
                        {
                            //right
                            return 4;
                        }
                    }

                    if (Math.Abs(y - (-2)) < 1.0e-6 && x < 10 + 1.0e-6)
                    {
                        //bottom front
                        return 2;
                    }

                    //if (Math.Abs(y - (0)) < 1.0e-6)
                    //{
                    //    //symmetry line
                    //    return 2;
                    //}

                    if (Math.Abs(y - (+2)) < 1.0e-6 && x < 10 + 1.0e-6)
                    {
                        //top front
                        return 3;
                    }

                    if (Math.Abs(y - (-0.5)) < 1.0e-6 && x > 10 - 1.0e-6)
                    {
                        // bottom back
                        return 2;
                    }

                    if (Math.Abs(y - (+0.5)) < 1.0e-6 && x > 10 - 1.0e-6)
                    {
                        // top back
                        return 3;
                    }

                    if (Math.Abs(x - (10)) < 1.0e-6 && y < -0.5 - 1.0e-6)
                    {
                        // bottom contraction
                        return 5;
                    }

                    if (Math.Abs(x - (10)) < 1.0e-6 && y > 0.5 - 1.0e-6)
                    {
                        //top contraction
                        return 6;
                    }

                    throw new ArgumentOutOfRangeException("at x = " + x + "and y = " + y);
                });


                return grd;
            };


            // Analytical Sol for Params
            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
                C.PresFunction = X => Pressurefunction(X, 0);
            }

            // Set Initial Conditions
            if (C.SetInitialConditions == true)
            {

                C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXX", X => StressXXfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));

                if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true)
                {
                    C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));
                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            C.AddBoundaryValue("Wall_bottom");
            C.AddBoundaryValue("Wall_top");
            C.AddBoundaryValue("Wall_Contraction_bottom");
            C.AddBoundaryValue("Wall_Contraction_top");
            //C.AddBoundaryCondition("FreeSlip");
            

            if (!C.FixedStreamwisePeriodicBC)
            {
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
        static public RheologyControl ConfinedCylinder(string path = @"\\dc1\userspace\kikker\cluster\cluster_db\ConfinedCylinder", int degree = 2) {
            RheologyControl C = new RheologyControl();

            //Path für cluster
            //\\dc1\userspace\kikker\cluster\cluster_db\ConfinedCylinder

            //Path für lokale DB
            //C:\AnnesBoSSSdb\ConfinedCylinder

            //Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            C.DbPath = path;
            C.ProjectName = "Cylinder";

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1E-6;

            C.LinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MinSolverIterations = 1;          
            C.LinearSolver.ConvergenceCriterion = 1E-6;

            //C.UnderRelax = 1.0;
            C.dt = 1e6;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Newton;
            C.LinearSolver.SolverCode = LinearSolverConfig.Code.exp_gmres_levelpmg;
            C.ObjectiveParam = 1.0;
            C.useJacobianForOperatorMatrix = false;

            C.UsePerssonSensor = true;
            C.SensorLimit = 1e-4;

            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 10;

            C.UseArtificialDiffusion = false;

            C.Bodyforces = true;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = false;

            //Physical Params
            double u0 = 1.5; // 0.375;// 0.66;// 3 / 2;
            double h = 4;

            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.beta = 0.59;
            C.Reynolds = 1;
            C.Weissenberg = 0.1; //aim Weissenberg number!
            C.RaiseWeissenberg = false;
            C.WeissenbergIncrement = 0.1;

            //Penalties
            C.ViscousPenaltyScaling = 1;
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
            Func<double[], double, double> VelocityXfunction = (X, t) => u0 * (1  - (X[1] *  X[1])/h);
            Func<double[], double, double> VelocityYfunction = (X, t) => 0.0;
            Func<double[], double, double> Pressurefunction = (X, t) => u0 * 0.5 * C.Reynolds * (35 - X[0]);
            Func<double[], double, double> StressXXfunction = (X, t) =>  2 * C.Weissenberg * (1 - C.beta) * u0 * (-2 / h) * X[1] * u0 * (-2 / h) * X[1];
            Func<double[], double, double> StressXYfunction = (X, t) => (1 - C.beta) * u0 * (-2 / h) * X[1];
            Func<double[], double, double> StressYYfunction = (X, t) => (0.0);

            // Insert Exact Solution
            C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            C.ExSol_Pressure = Pressurefunction;
            C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            // Create Fields
            //int degree = 2;
            C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("ResidualMomentumX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("ResidualMomentumY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("ResidualConti", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("ResidualStressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("ResidualStressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("ResidualStressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            // Create Grid
            //fine grid - only on cluster!           
            //string grid = "70797022-eba0-4c77-b179-334c665044b5";

            //more refined in wake of cylinder - only on cluster!
            //string grid = "3637610b-bcdf-4cdd-a647-cd7f91e373e8";


            //coarser grid - works without cluster!
            string grid = "f9aa12dc-53bb-4e2c-81b3-ffccc251a3f7";

            //very coarse grid as starting point for refinement
            //string grid = "e296a1b2-98f9-4fdf-8a32-04e0954ff369";

            //Dennis Zylinder for drag validation
            //string grid = "a67192f5-6b59-4caf-a95a-0a08730c3365";


            Guid gridGuid;
            if (Guid.TryParse(grid, out gridGuid))
            {
                C.GridGuid = gridGuid;
            }
            else
            {
                C.GridFunc = delegate ()
                {
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

            //restart
            //var database = new DatabaseInfo(path);
            //Guid restartID = new Guid("9ae08191-ee15-4803-9e3f-566f119c9de4");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

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

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            C.AddBoundaryValue("Wall_bottom", "VelocityX", X => 0);
            C.AddBoundaryValue("Wall_top", "VelocityX", X => 0);
            C.AddBoundaryValue("Wall_bottom", "VelocityY", X => 0);
            C.AddBoundaryValue("Wall_top", "VelocityY", X => 0);
            C.AddBoundaryValue("Wall_cylinder", "VelocityX", X => 0);
            C.AddBoundaryValue("Wall_cylinder", "VelocityY", X => 0);


            if (!C.FixedStreamwisePeriodicBC)
            {
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

    }
}