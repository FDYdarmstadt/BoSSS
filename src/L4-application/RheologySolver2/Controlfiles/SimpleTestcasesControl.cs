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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.Rheology
{
    /// <summary>
    /// Control file of all simple testcases for debugging purpose, e.g. channel or cnsistency checks
    /// </summary>
    static public class SimpleTestcasesControl
    {
        /// <summary>
        /// Channel Flow
        /// </summary>
        static public RheologyControl Channel(string path = @"C:\AnnesBoSSSdb\Channel", int degree = 2, int GridLevel = 5)
        {
            RheologyControl C = new RheologyControl();

            //Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            C.DbPath = path;
            C.ProjectName = "Channel";
            C.NonLinearSolver.MaxSolverIterations = 30;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.NonLinearSolver.ConvergenceCriterion = 1E-10;
            C.LinearSolver.MaxSolverIterations = 30;
            C.LinearSolver.MinSolverIterations = 3;
            C.LinearSolver.ConvergenceCriterion = 1E-10;
            C.dt = 0.1;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Newton;
            C.ObjectiveParam = 1.0;

            C.UsePerssonSensor = false;
            C.SensorLimit = 1e-4;

            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 10;

            C.UseArtificialDiffusion = false;

            C.Bodyforces = false;
            //C.WhichWall = "Wall_Cylinder";

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = true;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = false;
            C.GravitySource = false;
            C.GravityX = (X, t) => 1;
            C.GravityY = (X, t) => 0;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.beta = 0.59;
            C.Reynolds = 1;
            C.Weissenberg = 1.0; //aim Weissenberg number!
            C.RaiseWeissenberg = false;
            C.WeissenbergIncrement = 0.1;

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
            C.alpha = 1;
            C.StressPenalty = 1.0;

            //Exact Solution Channel
            Func<double[], double, double> VelocityXfunction = (X, t) => 1 - (X[1] * X[1]);
            Func<double[], double, double> VelocityYfunction = (X, t) => 0;
            Func<double[], double, double> Pressurefunction = (X, t) => 2* C.Reynolds * (20 - X[0]);
            Func<double[], double, double> StressXXfunction = (X, t) => 0;// 2  * (1 - C.beta) * (((-2 * X[1])) * ((-2 * X[1])));
            Func<double[], double, double> StressXYfunction = (X, t) => (1 - C.beta) * (-2 * X[1]);
            Func<double[], double, double> StressYYfunction = (X, t) => (0.0);

            // Insert Exact Solution
            C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            C.ExSol_Pressure = Pressurefunction;
            C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            // Create Fields
            C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 20, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, (cells2 / 4) + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                }

                grd.EdgeTagNames.Add(2, "Wall_bottom");
                grd.EdgeTagNames.Add(3, "Wall_top");
                //grd.EdgeTagNames.Add(2, "FreeSlip");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom
                        return 2;

                    if (Math.Abs(y - (1)) < 1.0e-6)
                        // top
                        return 3;

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (0)) < 1.0e-6)
                            // left
                            return 1;

                        if (Math.Abs(x - (20)) < 1.0e-6)
                            // right
                            return 4;
                    }
                    throw new ArgumentOutOfRangeException();
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
            C.AddBoundaryValue("Wall_bottom", "VelocityX", VelocityXfunction);
            C.AddBoundaryValue("Wall_top", "VelocityX", VelocityXfunction);
            //C.AddBoundaryValue("Wall_bottom", "VelocityY", VelocityYfunction);
            //C.AddBoundaryValue("Wall_top", "VelocityY", VelocityYfunction);
            //C.AddBoundaryValue("Wall_bottom", "VelocityX", X => 0);
            //C.AddBoundaryValue("Wall_top", "VelocityX", X => 0);
            //C.AddBoundaryValue("Wall_bottom", "VelocityY", X => 0);
            //C.AddBoundaryValue("Wall_top", "VelocityY", X => 0);
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
        /// Consistency Constitutive equation
        /// </summary>
        static public RheologyControl ConsistencyConstitutive(string path = @"C:\AnnesBoSSSdb\ConsistencyConstitutive_withDiv", int degree = 1, int GridLevel = 5)
        {


            RheologyControl C = new RheologyControl();

            // Solver Options
            C.NoOfTimesteps = 10;
            C.savetodb = false;
            C.DbPath = path;
            C.SessionName = "Degree" + degree + ", GridLevel" + GridLevel;
            C.ProjectName = "ConsistencyStudyConstitutive";
            C.NonLinearSolver.MaxSolverIterations = 30;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion= 1E-7;
            C.LinearSolver.MaxSolverIterations = 30;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.ConvergenceCriterion = 1E-7;
            //C.ConvCritGMRES = 1E-13;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.NewtonGMRES;
            C.ObjectiveParam = 1.0;

            //Grid Params
            //double GridLevel = 5;
            double h = Math.Pow(2, -GridLevel + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = false;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = true;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.GravitySource = true;
            C.beta = 1;
            C.Reynolds = 1;
            C.Weissenberg = 0;
            C.RaiseWeissenberg = true;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 1;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 1;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;
            C.StressPenalty = 1.0;

            // Exact Solution manufactured Solution
            Func<double[], double, double> VelocityXfunction = (X, t) => X[0] * X[0];
            Func<double[], double, double> VelocityYfunction = (X, t) => -X[1];
            Func<double[], double, double> Pressurefunction = (X, t) => 0;
            Func<double[], double, double> StressXXfunction = (X, t) => X[0] * X[0];
            Func<double[], double, double> StressXYfunction = (X, t) => X[0] * X[0] + X[1] * X[1];
            Func<double[], double, double> StressYYfunction = (X, t) => X[1] * X[1];

            // Exact Solution manufactured Solution turned 90°
            //Func<double[], double, double> VelocityXfunction = (X, t) => -X[1];
            //Func<double[], double, double> VelocityYfunction = (X, t) => X[0] * X[0];
            //Func<double[], double, double> Pressurefunction = (X, t) => 0;
            //Func<double[], double, double> StressXXfunction = (X, t) => X[1] * X[1];
            //Func<double[], double, double> StressXYfunction = (X, t) => X[0] * X[0] + X[1] * X[1];
            //Func<double[], double, double> StressYYfunction = (X, t) => X[0] * X[0];

            //Gravity sources
            // Weissenberg = 1
            //C.GravityX = (X, t) => 2 * X[0] * X[0] * X[0] + X[0] * X[0] * (2 * X[0] - 1) - 2 * X[0] - 2 * X[1];
            //C.GravityY = (X, t) => -X[1] - X[1] * (2 * X[0] - 1) - 2 * X[0];
            //C.GravityXX = (X, t) => X[0] * X[0] + 2 * X[0] * X[0] * X[0] - 4 * X[0];
            //C.GravityXY = (X, t) => 2 * X[0] * X[0] * X[0] + X[0] * X[0] - X[1] * X[1];
            //C.GravityYY = (X, t) => -X[1] * X[1] + 2;
            //C.GravityDiv = (X, t) => 2 * X[0] - 1;

            ////Weissenberg = 0
            //C.GravityX = (X, t) => 2 * X[0] * X[0] * X[0] + X[0] * X[0] * (2 * X[0] - 1) - 2 * X[0] - 2 * X[1];
            //C.GravityY = (X, t) => -X[1] - X[1] * (2 * X[0] - 1) - 2 * X[0];
            //C.GravityXX = (X, t) => X[0] * X[0] - 4 * X[0];
            //C.GravityXY = (X, t) => X[0] * X[0] + X[1] * X[1];
            //C.GravityYY = (X, t) => X[1] * X[1] + 2;
            //C.GravityDiv = (X, t) => 2 * X[0] - 1;

            //Weissenberg variable including Objective Terms and beta!
            C.GravityX = (X, t) => 2 * X[0] * X[0] * X[0] + X[0] * X[0] * (2 * X[0] - 1) - 2 * X[0] - 2 * X[1] - 2 * C.beta;
            C.GravityY = (X, t) => -X[1] - X[1] * (2 * X[0] - 1) - 2 * X[0];
            C.GravityXX = (X, t) => X[0] * X[0] - 2 * C.Weissenberg * X[0] * X[0] * X[0] - 4 * (1 - C.beta) * X[0];
            C.GravityXY = (X, t) => (X[0] * X[0] + X[1] * X[1] + C.Weissenberg * (2 * X[0] * X[0] * X[0] - 2 * X[1] * X[1] - (2 * X[0] - 1) * (X[0] * X[0] + X[1] * X[1])));
            C.GravityYY = (X, t) => X[1] * X[1] + 2 - 2 * C.beta;
            C.GravityDiv = (X, t) => 2 * X[0] - 1;

            //Weissenberg variable including Objective Terms and beta! ohne grad(u)^T
            //C.GravityX = (X, t) => 2 * X[0] * X[0] * X[0] + X[0] * X[0] * (2 * X[0] - 1) - 2 * X[0] - 2 * X[1] - 2 * C.beta;
            //C.GravityY = (X, t) => -X[1] - X[1] * (2 * X[0] - 1) - 2 * X[0];
            //C.GravityXX = (X, t) => X[0] * X[0] - 2 * C.Weissenberg * X[0] * X[0] * X[0] - 2 * (1 - C.beta) * X[0];
            //C.GravityXY = (X, t) => (X[0] * X[0] + X[1] * X[1] + C.Weissenberg * (2 * X[0] * X[0] * X[0] - 2 * X[1] * X[1] - (2 * X[0] - 1) * (X[0] * X[0] + X[1] * X[1])));
            //C.GravityYY = (X, t) => X[1] * X[1] + 1 - 1 * C.beta;
            //C.GravityDiv = (X, t) => 2 * X[0] - 1;

            //Weissenberg variable including Objective Terms and beta turned 90°
            //C.GravityX = (X, t) => -X[0] * X[0] - 2 * X[1];
            //C.GravityY = (X, t) => -2 * X[0] * X[1] - 2 * X[0] - 2 * C.beta;
            //C.GravityXX = (X, t) => X[1] * X[1] + C.Weissenberg * (2 * X[0] * X[0] * X[1] + 2 * X[0] * X[0] + 2 * X[1] * X[1]);
            //C.GravityXY = (X, t) => (X[0] * X[0]) + (X[1] * X[1]) + (C.Weissenberg * (2 * X[0] * X[0] * X[1] - 2 * X[0] * X[1] * X[1]
            //                        + X[0] * X[0] - 2 * X[1] * X[0])) - (1 - C.beta) * (-1 + 2 * X[0]);
            //C.GravityYY = (X, t) => X[0] * X[0] + C.Weissenberg * (-2 * X[1] * X[0] - 4 * X[0] * (X[0] * X[0] + X[1] * X[1]));
            //C.GravityDiv = (X, t) => 0;

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

            //C.FieldOptions.Add("VelocityXGradientX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            //C.FieldOptions.Add("VelocityXGradientY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            //C.FieldOptions.Add("VelocityYGradientX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            //C.FieldOptions.Add("VelocityYGradientY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });


            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC); //C.FixedStreamwisePeriodicBC periodicX:true

                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                }

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom
                        return 1;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top
                        return 1;

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (-1)) < 1.0e-6)
                            // left
                            return 1;

                        if (Math.Abs(x - (+1)) < 1.0e-6)
                            // right
                            return 1;
                    }

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // Analytical Sol for Params
            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
            }

            // Set Initial Conditions
            if (C.SetInitialConditions == true)
            {

                C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXX", X => StressXXfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));
                C.InitialValues_Evaluators.Add(VariableNames.GravityX, X => C.GravityX(X, 0));
                C.InitialValues_Evaluators.Add(VariableNames.GravityY, X => C.GravityY(X, 0));
                C.InitialValues_Evaluators.Add("GravityXX", X => C.GravityXX(X, 0));
                C.InitialValues_Evaluators.Add("GravityXY", X => C.GravityXY(X, 0));
                C.InitialValues_Evaluators.Add("GravityYY", X => C.GravityYY(X, 0));
                C.InitialValues_Evaluators.Add("GravityDiv", X => C.GravityDiv(X, 0));

                if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true)
                {
                    C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));

                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions

            if (!C.FixedStreamwisePeriodicBC)
            {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", VelocityXfunction);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", VelocityYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", StressXXfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", StressXYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", StressYYfunction);
                C.AddBoundaryValue("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Velocity_inlet", "GravityX", C.GravityX);
                C.AddBoundaryValue("Velocity_inlet", "GravityY", C.GravityY);
                C.AddBoundaryValue("Velocity_inlet", "GravityXX", C.GravityXX);
                C.AddBoundaryValue("Velocity_inlet", "GravityXY", C.GravityXY);
                C.AddBoundaryValue("Velocity_inlet", "GravityYY", C.GravityYY);

                //C.AddBoundaryCondition("Pressure_Outlet", "Pressure", Pressurefunction);

            }
            return C;
        }

        //__________________________________________________________________________________________________________________     
        /// <summary>
        /// Stagnation point flow (test of the wall BC for constitutive equations
        /// </summary>
        static public RheologyControl Staupunkt(string path = @"C:\AnnesBoSSSdb\Staupunkt", int degree = 2, int GridLevel = 5)
        {


            RheologyControl C = new RheologyControl();

            // Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            C.DbPath = path;
            C.SessionName = "Degree" + degree + ", GridLevel" + GridLevel;
            C.ProjectName = "Staupunkt";
            //C.MaxIter = 20;
            //C.MinIter = 3;
            //C.ConvCrit = 1E-14;
            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.NonLinearSolver.ConvergenceCriterion = 1E-14;
            C.LinearSolver.MaxSolverIterations = 20;
            C.LinearSolver.MinSolverIterations = 3;
            C.LinearSolver.ConvergenceCriterion = 1E-14;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;

            //Grid Params
            //double GridLevel = 5;
            double h = Math.Pow(2, -GridLevel + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = true;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.GravitySource = false;
            C.beta = 1;
            C.Reynolds = 1;
            C.Weissenberg = 0;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 20;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 1;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;

            // Exact Solution manufactured Solution
            Func<double[], double, double> VelocityXfunction = (X, t) => 0;
            Func<double[], double, double> VelocityYfunction = (X, t) => -X[1];
            Func<double[], double, double> Pressurefunction = (X, t) => 0;
            Func<double[], double, double> StressXXfunction = (X, t) => 0;
            Func<double[], double, double> StressXYfunction = (X, t) => 0;
            Func<double[], double, double> StressYYfunction = (X, t) => 0;

            //Gravity sources
            C.GravityX = (X, t) => 0;
            C.GravityY = (X, t) => 0;
            C.GravityXX = (X, t) => 0;
            C.GravityXY = (X, t) => 0;
            C.GravityYY = (X, t) => 0;
            C.GravityDiv = (X, t) => 0;

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
            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });


            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(2, "Wall");
                    grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                }

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

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (-1)) < 1.0e-6)
                            // left
                            return 3;

                        if (Math.Abs(x - (+1)) < 1.0e-6)
                            // right
                            return 3;
                    }

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // Analytical Sol for Params
            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
            }

            // Set Initial Conditions
            if (C.SetInitialConditions == true)
            {

                C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXX", X => StressXXfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));
                C.InitialValues_Evaluators.Add(VariableNames.GravityX, X => C.GravityX(X, 0));
                C.InitialValues_Evaluators.Add(VariableNames.GravityY, X => C.GravityY(X, 0));
                C.InitialValues_Evaluators.Add("GravityXX", X => C.GravityXX(X, 0));
                C.InitialValues_Evaluators.Add("GravityXY", X => C.GravityXY(X, 0));
                C.InitialValues_Evaluators.Add("GravityYY", X => C.GravityYY(X, 0));

                if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true)
                {
                    C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));

                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            C.AddBoundaryValue("Wall");

            if (!C.FixedStreamwisePeriodicBC)
            {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", VelocityXfunction);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", VelocityYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", StressXXfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", StressXYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", StressYYfunction);
                //C.AddBoundaryCondition("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Velocity_inlet", "GravityX", C.GravityX);
                C.AddBoundaryValue("Velocity_inlet", "GravityY", C.GravityY);
                C.AddBoundaryValue("Velocity_inlet", "GravityXX", C.GravityXX);
                C.AddBoundaryValue("Velocity_inlet", "GravityXY", C.GravityXY);
                C.AddBoundaryValue("Velocity_inlet", "GravityYY", C.GravityYY);

                C.AddBoundaryValue("Pressure_Outlet", "Pressure", Pressurefunction);

            }
            return C;
        }
        //__________________________________________________________________________________________________________________     
        /// <summary>
        /// Convergence of the Stokes system with LDG
        /// </summary>
        static public RheologyControl ConvergenceStokesLDG(string path = @"C:\AnnesBoSSSdb\ConvergenceStokesLDG_Paper", int degree = 2, int GridLevel = 2)
        {
            // Path wenn lokal gerechnet wird: C:\AnnesBoSSSdb\ConvergenceStokesLDG 
            // Path für lokal 2. Versuch ohne penalty in RB: C:\AnnesBoSSSdb\ConvergenceStokesLDG2exclBEpen
            // Path wenn mit Cluster gerechnet wird: \\dc1\userspace\kikker\cluster\cluster_db\ConvergenceStudyLDG

            RheologyControl C = new RheologyControl();

            // Solver Options
            C.NoOfTimesteps = 5;
            C.savetodb = false;
            C.DbPath = path;
            C.SessionName = "Degree" + degree + ", GridLevel" + GridLevel;
            C.ProjectName = "ConvStudyLDG";
            //C.MaxIter = 3;
            //C.MinIter = 1;
            //C.ConvCrit = 5E-7;
            C.NonLinearSolver.MaxSolverIterations = 3;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 5E-7;
            C.LinearSolver.MaxSolverIterations = 3;
            C.LinearSolver.MinSolverIterations = 1;
            C.LinearSolver.ConvergenceCriterion = 5E-7;
            //C.ConvCritGMRES = 1E-08;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            //C.NonlinearMethod = NonlinearSolverMethod.Newton;
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Newton;
            C.ObjectiveParam = 1.0;
            C.UsePerssonSensor = true;
            C.AdaptiveMeshRefinement = true;
            C.RefinementLevel = 3;

            //Grid Params
            //double GridLevel = 5;
            double h = Math.Pow(2, -GridLevel + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;

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
            C.GravitySource = true;
            C.beta = 0;
            C.Reynolds = 1;
            C.Weissenberg = 0.0;
            C.RaiseWeissenberg = false;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 1;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 1;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;
            C.StressPenalty = 1.0;

            // Exact Solution manufactured Solution
            Func<double[], double, double> VelocityXfunction = (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> VelocityYfunction = (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]);
            Func<double[], double, double> Pressurefunction = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]);
            Func<double[], double, double> StressXXfunction = (X, t) => -2 * (1 - C.beta) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> StressXYfunction = (X, t) => -2 * (1 - C.beta) * Math.Exp(X[0]) * (Math.Cos(X[1]) - X[1] * Math.Sin(X[1]));
            Func<double[], double, double> StressYYfunction = (X, t) => 2 * (1 - C.beta) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));

            //Gravity Stokes system with We = 0
            //C.GravityX = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]) * ((C.Reynolds + C.beta - 1) / C.Reynolds);
            //C.GravityY = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.Reynolds + C.beta - 1) / C.Reynolds);
            //C.GravityXX = (X, t) => 0;
            //C.GravityXY = (X, t) => 0;
            //C.GravityYY = (X, t) => 0;
            //C.GravityDiv = (X, t) => 0;


            //Gravity Stokes system without objective terms
            //C.GravityX = (X, t) => -1 / C.Reynolds * Math.Exp(X[0]) * (Math.Exp(X[0]) * Math.Cos(X[1]) * Math.Cos(X[1]) * C.Reynolds
            //- Math.Exp(X[0]) * C.Reynolds * X[1] * X[1] - Math.Exp(X[0]) * C.Reynolds - 2 * Math.Sin(X[1]) * C.Reynolds + 2 * Math.Sin(X[1]));
            //C.GravityY = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.Reynolds - 1) / C.Reynolds);
            //C.GravityXX = (X, t) => 2 * Math.Exp((2 * X[0])) * C.Weissenberg * (Math.Pow(Math.Cos(X[1]), 2) - X[1] * X[1] - 1) * (C.beta - 1);
            //C.GravityXY = (X, t) => -2 * Math.Exp((2 * X[0])) * C.Weissenberg * (Math.Cos(X[1]) * Math.Sin(X[1]) * C.beta - Math.Cos(X[1]) * Math.Sin(X[1]) + C.beta * X[1] - X[1]);
            //C.GravityYY = (X, t) => -2 * Math.Exp((2 * X[0])) * C.Weissenberg * (Math.Pow(Math.Cos(X[1]), 2) - X[1] * X[1] - 1) * (C.beta - 1);
            //C.GravityDiv = (X, t) => 0;

            //Gravity for full system
            C.GravityX = (X, t) => -1 / C.Reynolds * Math.Exp(X[0]) * (Math.Exp(X[0]) * Math.Cos(X[1]) * Math.Cos(X[1]) * C.Reynolds
                        - Math.Exp(X[0]) * C.Reynolds * X[1] * X[1] - Math.Exp(X[0]) * C.Reynolds - 2 * Math.Sin(X[1]) * C.Reynolds + 2 * Math.Sin(X[1]));

            C.GravityY = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.Reynolds - 1) / C.Reynolds);

            C.GravityXX = (X, t) => 2 * Math.Exp(2 * X[0]) * C.Weissenberg * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.beta * X[1] * C.Weissenberg
                        + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.beta + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] + C.beta * X[1] * X[1]
                        - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) - X[1] * X[1] + C.beta - 1);

            C.GravityXY = (X, t) => -2 * Math.Exp(2 * X[0]) *(C.beta - 1) * C.Weissenberg * (2 * Math.Pow(Math.Cos(X[1]), 2) * X[1] + 3 * Math.Cos(X[1]) * Math.Sin(X[1]) + X[1]);

            C.GravityYY = (X, t) => -2 * Math.Exp(2 * X[0]) * C.Weissenberg * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.beta * X[1]
                        + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.beta + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] - 3 * C.beta * X[1] * X[1]
                        - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) + 3 * X[1] * X[1] - 3 * C.beta + 3);

            C.GravityDiv = (X, t) => 0;

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

            C.FieldOptions.Add("VelocityXGradientX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityXGradientY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityYGradientX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityYGradientY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });


            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC); //C.FixedStreamwisePeriodicBC

                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    //grd.EdgeTagNames.Add(2, "Pressure_Outlet");
                }

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom
                        return 1;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top
                        return 1;

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (-1)) < 1.0e-6)
                            // left
                            return 1;

                        if (Math.Abs(x - (+1)) < 1.0e-6)
                            // right
                            return 1;
                    }

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // Analytical Sol for Params
            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
            }

            // Set Initial Conditions
            if (C.SetInitialConditions == true)
            {

                C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXX", X => StressXXfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));
                C.InitialValues_Evaluators.Add(VariableNames.GravityX, X => C.GravityX(X, 0));
                C.InitialValues_Evaluators.Add(VariableNames.GravityY, X => C.GravityY(X, 0));

                if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true)
                {
                    C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));

                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            if (!C.FixedStreamwisePeriodicBC)
            {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", VelocityXfunction);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", VelocityYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", StressXXfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", StressXYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", StressYYfunction);
                C.AddBoundaryValue("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Velocity_inlet", "GravityX", C.GravityX);
                C.AddBoundaryValue("Velocity_inlet", "GravityY", C.GravityY);
                //C.AddBoundaryCondition("Pressure_Outlet", "Pressure", Pressurefunction);

            }
            return C;
        }
        //__________________________________________________________________________________________________________________       
        /// <summary>
        /// Channel Flow with moving wall
        /// </summary>
        static public RheologyControl MovingWallChannel()
        {
            RheologyControl C = new RheologyControl();

            // Solver Options
            C.ViscousPenaltyScaling = 1;
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            C.DbPath = @"C:\AnnesBoSSSdb\Channel";
            C.ProjectName = "Channel";
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = true;
            //C.MaxIter = 7;
            //C.MinIter = 3;
            //C.ConvCrit = 1E-13;
            C.NonLinearSolver.MaxSolverIterations = 7;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.NonLinearSolver.ConvergenceCriterion = 1E-13;
            C.LinearSolver.MaxSolverIterations = 7;
            C.LinearSolver.MinSolverIterations = 3;
            C.LinearSolver.ConvergenceCriterion = 1E-13;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;


            // Create Fields
            int degree = 2;
            C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = degree - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });


            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 10, 41);
                var _yNodes = GenericBlas.Linspace(-1, 1, 11);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                }

                grd.EdgeTagNames.Add(2, "Wall_bottom");
                grd.EdgeTagNames.Add(3, "Velocity_inlet_oben");

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

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (0)) < 1.0e-6)
                            // left
                            return 1;

                        if (Math.Abs(x - (10)) < 1.0e-6)
                            // right
                            return 4;
                    }
                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => (1 - (X[1] * X[1]));
                C.VelFunctionV = X => 0;
            }

            Func<double[], double> stressXXfunction = X => 2 * C.Weissenberg * (1 - C.beta) * (((-2 * X[1])) * ((-2 * X[1])));
            Func<double[], double> stressXYfunction = X => (1 - C.beta) * ((-2 * X[1]));

            // Set Initial Conditions
            if (C.SetInitialConditions == true)
            {

                C.InitialValues_Evaluators.Add("VelocityX", X => 0.5 * X[1] + 0.5);
                C.InitialValues_Evaluators.Add("VelocityY", X => 0);
                C.InitialValues_Evaluators.Add("StressXX", X => 2 * C.Weissenberg * (1 - C.beta) * 0.25);
                C.InitialValues_Evaluators.Add("StressXY", X => (1 - C.beta) * 0.5);
                C.InitialValues_Evaluators.Add("StressYY", X => 0);
                C.InitialValues_Evaluators.Add("Phi", X => -1);

                if (C.SetInitialPressure == true)
                {
                    C.InitialValues_Evaluators.Add("Pressure", X => 10 * (2 * C.beta / C.Reynolds - (-2 + 2 * C.beta) / C.Reynolds) - (2 * C.beta / C.Reynolds - (-2 + 2 * C.beta) / C.Reynolds) * X[0]);

                    //only Stokes without extra stresses
                    //C.InitialValues_Evaluators.Add("Pressure", X => 10* (2 * C.beta / C.Reynolds) - (2 * C.beta / C.Reynolds) * X[0]);
                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            C.AddBoundaryValue("Wall_bottom");
            C.AddBoundaryValue("Wall_top");
            C.AddBoundaryValue("Velocity_inlet_oben", "VelocityX", X => 0.5 * X[1] + 0.5);
            C.AddBoundaryValue("Velocity_inlet_oben", "StressXX", X => 2 * C.Weissenberg * (1 - C.beta) * 0.25);
            C.AddBoundaryValue("Velocity_inlet_oben", "StressXY", X => (1 - C.beta) * 0.5);
            C.AddBoundaryValue("Velocity_inlet_oben", "StressYY", X => 0);

            if (!C.FixedStreamwisePeriodicBC)
            {

                C.AddBoundaryValue("Velocity_inlet", "VelocityX", X => 0.5 * X[1] + 0.5);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", X => 2 * C.Weissenberg * (1 - C.beta) * 0.25);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", X => (1 - C.beta) * 0.5);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", X => 0);
                C.AddBoundaryValue("Pressure_Outlet");

            }
            return C;
        }
        //__________________________________________________________________________________________________________________
        /// <summary>
        /// Parameter Study of Channel Flow
        /// </summary>
        static public RheologyControl[] ChannelParameterStudy()
        {

            List<RheologyControl> All = new List<RheologyControl>();

            //int[] DGdegree = new int[] {1, 2, 3, 4 };
            //int[] relDGdegree = new int[] { -1, 0, 1, 2 };
            //int[] GrdRes = new int[] { 2, 5, 10, 15 };

            int[] DGdegree = new int[] { 1, 2, 3, 4 };
            int[] relDGdegree = new int[] { 0, 1, 2 };
            int[] GrdRes = new int[] { 2, 5, 10 };

            foreach (int p in DGdegree)
            {
                foreach (int q in relDGdegree)
                {
                    foreach (int h in GrdRes)
                    {
                        var C = new RheologyControl();

                        // Solver Options
                        C.ViscousPenaltyScaling = 1;
                        C.NoOfTimesteps = 1;
                        C.savetodb = false;
                        C.DbPath = @"C:\AnnesBoSSSdb\Channel";
                        C.ProjectName = "Channel";
                        C.Stokes = false;
                        C.FixedStreamwisePeriodicBC = true;
                        //C.MaxIter = 2;
                        //C.MinIter = 1;
                        //C.ConvCrit = 1E-10;
                        C.NonLinearSolver.MaxSolverIterations = 2;
                        C.NonLinearSolver.MinSolverIterations = 1;
                        C.NonLinearSolver.ConvergenceCriterion = 1E-10;
                        C.LinearSolver.MaxSolverIterations = 2;
                        C.LinearSolver.MinSolverIterations = 1;
                        C.LinearSolver.ConvergenceCriterion = 1E-10;


                        // Create Fields
                        C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = p + q, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
                        C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = p + q, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
                        C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = p, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

                        C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = p + q, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
                        C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = p + q, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
                        C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = p + q, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
                        C.FieldOptions.Add("PhiDG", new FieldOpts()
                        {
                            Degree = 0,
                            SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                        });
                        C.FieldOptions.Add("Phi", new FieldOpts()
                        {
                            Degree = 0,
                            SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                        });

                        // Create Grid
                        C.GridFunc = delegate {
                            var _xNodes = GenericBlas.Linspace(0, 10, 2 * h + 1);
                            var _yNodes = GenericBlas.Linspace(-1, 1, h + 1);

                            var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                            if (!C.FixedStreamwisePeriodicBC)
                            {
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

                                if (!C.FixedStreamwisePeriodicBC)
                                {
                                    if (Math.Abs(x - (0)) < 1.0e-6)
                                        // left
                                        return 1;

                                    if (Math.Abs(x - (10)) < 1.0e-6)
                                        // right
                                        return 4;
                                }
                                throw new ArgumentOutOfRangeException();
                            });

                            return grd;
                        };

                        // Set Initial Conditions
                        C.InitialValues_Evaluators.Add("VelocityX", X => 1 - (X[1] * X[1]));
                        C.InitialValues_Evaluators.Add("VelocityY", X => 0);
                        C.InitialValues_Evaluators.Add("ExternalStressXX", X => 0);
                        C.InitialValues_Evaluators.Add("ExternalStressXY", X => -X[1]);
                        C.InitialValues_Evaluators.Add("ExternalStressYY", X => 0);


                        // Set Boundary Conditions
                        C.AddBoundaryValue("Wall_bottom");
                        C.AddBoundaryValue("Wall_top");

                        if (!C.FixedStreamwisePeriodicBC)
                        {
                            C.AddBoundaryValue("Velocity_inlet", "VelocityX", X => (1 - (X[1] * X[1])));
                            C.AddBoundaryValue("Pressure_Outlet");
                        }


                        C.Paramstudy_CaseIdentification = new Tuple<string, object>[] {

                                    new Tuple<string, object>("DGdegree", p),

                                    new Tuple<string, object>("relativeDGdegree", q),

                                    new Tuple<string, object>("GridResolutionFactor", h),

                                    };
                        All.Add(C);

                    }
                }
            }

            return All.ToArray();
        }
        //__________________________________________________________________________________________________________________
        /// <summary>
        /// Channel test for Flow
        /// </summary>
        static public RheologyControl ChannelLDG(string path = @"\\dc1\userspace\kikker\cluster\cluster_db\ConvergenceStudyLDG", int degree = 2, int GridLevel = 5)
        {
            // Path wenn lokal gerechnet wird: C:\AnnesBoSSSdb\ConvergenceStokesLDG 
            // Path für lokal 2. Versuch ohne penalty in RB: C:\AnnesBoSSSdb\ConvergenceStokesLDG2exclBEpen
            // Path wenn mit Cluster gerechnet wird: \\dc1\userspace\kikker\cluster\cluster_db\ConvergenceStudyLDG

            RheologyControl C = new RheologyControl();

            // Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            C.DbPath = path;
            C.SessionName = "Degree" + degree + ", GridLevel" + GridLevel;
            C.ProjectName = "ConvStudyLDG";
            //C.MaxIter = 15;
            //C.MinIter = 3;
            //C.ConvCrit = 1E-14;
            C.NonLinearSolver.MaxSolverIterations = 15;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.NonLinearSolver.ConvergenceCriterion = 1E-14;
            C.LinearSolver.MaxSolverIterations = 15;
            C.LinearSolver.MinSolverIterations = 3;
            C.LinearSolver.ConvergenceCriterion = 1E-14;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;

            //Grid Params
            //double GridLevel = 5;
            double h = Math.Pow(2, -GridLevel + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = true;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.beta = 0;
            C.Reynolds = 1;
            C.Weissenberg = 0;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 1 / h;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = h;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;

            //Exact Solution Channel
            Func<double[], double, double> VelocityXfunction = (X, t) => 4 - X[1] * X[1];
            Func<double[], double, double> VelocityYfunction = (X, t) => (0.0);
            Func<double[], double, double> Pressurefunction = (X, t) => 2 / C.Reynolds * (20 - X[0]);
            Func<double[], double, double> StressXXfunction = (X, t) => 2 * C.Weissenberg * (1 - C.beta) * (((-2 * X[1])) * ((-2 * X[1])));
            Func<double[], double, double> StressXYfunction = (X, t) => (1 - C.beta) * (-2 * X[1]);
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
            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });


            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 20, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-2, 2, cells2 + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);


                if (!C.FixedStreamwisePeriodicBC)
                {
                    grd.EdgeTagNames.Add(1, "Velocity_inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                }

                grd.EdgeTagNames.Add(2, "Wall_bottom");
                grd.EdgeTagNames.Add(3, "Wall_top");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-2)) < 1.0e-6)
                        // bottom
                        return 2;

                    if (Math.Abs(y - (+2)) < 1.0e-6)
                        // top
                        return 3;

                    if (!C.FixedStreamwisePeriodicBC)
                    {
                        if (Math.Abs(x - (0)) < 1.0e-6)
                            // left
                            return 1;

                        if (Math.Abs(x - (20)) < 1.0e-6)
                            // right
                            return 4;
                    }

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            // Analytical Sol for Params
            if (C.SetParamsAnalyticalSol == true)
            {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
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


            if (!C.FixedStreamwisePeriodicBC)
            {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", VelocityXfunction);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", VelocityYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", StressXXfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", StressXYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", StressYYfunction);
                //C.AddBoundaryCondition("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Pressure_Outlet", "Pressure", Pressurefunction);

            }
            return C;
        }
        //___________________________________________________________________________________________________________________________________________

    }
}
