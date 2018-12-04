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

namespace BoSSS.Application.Rheology {
    static public class RheologyLocalDGTest {

        static public RheologyControl LocalDGComputeRes(int GridRes = 2, int PolyDeg = 1, double beta = 0) {
            RheologyControl C = LocalDGGeneric(GridRes, PolyDeg, beta);

            Console.WriteLine("Test 1: Insert exact solution and only compute residual.");
            C.Stokes = true;
            C.SkipSolveAndEvaluateResidual = true;
            C.SetInitialConditions = true;
            C.SetInitialPressure = true;
            C.SetParamsAnalyticalSol = false;
            C.beta = beta;

            Func<double[], double, double> Pressurefunction = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]);
            C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));

            return C;
        }

        static public RheologyControl LocalDGStokes(int GridRes = 3, int PolyDeg = 1, double beta = 0.11) {
            RheologyControl C = LocalDGGeneric(GridRes, PolyDeg, beta);

            Console.WriteLine("Test 2: Compute Stokes system.");
            C.Stokes = true;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.beta = beta;

            return C;
        }

        static public RheologyControl LocalDGLinearized(int GridRes = 4, int PolyDeg = 2, double beta = 0)
        {
            RheologyControl C = LocalDGGeneric(GridRes, PolyDeg, beta);

            Console.WriteLine("Test 3: Insert exact solution as Parameter and solve linearized Problem.");
            C.Stokes = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = false;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = true;
            C.beta = beta;
            C.Weissenberg = 1.0;

            return C;
        }

        static public RheologyControl LocalDGNavierStokes(int GridRes = 2, int PolyDeg = 2, double beta = 0)
        {
            RheologyControl C = LocalDGGeneric(GridRes, PolyDeg, beta);

            Console.WriteLine("Test 4: Compute whole system and see if solution converges.");
            C.Stokes = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = false;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.beta = beta;
            C.Weissenberg = 1.0;

            return C;
        }

        // Convergence Manufactured Solution
        static public RheologyControl LocalDGGeneric(int GridLevel = 2, int PolyDeg = 2, double beta = 0) {
            RheologyControl C = new RheologyControl();

            // Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            //C.DbPath = @"C:\AnnesBoSSSdb\ConvergenceStokesLDG";
            C.ProjectName = "ConvStudyLDG";
            C.MaxIter = 10;
            C.MinIter = 10;
            C.ConvCrit = 1E-14;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;

            //Grid Params
            double h = Math.Pow(2, -GridLevel + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;
            C.grd = cells2;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.ComputeL2Error = true;

            //Physical Params
            C.Stokes = true;
            C.FixedStreamwisePeriodicBC = false;
            C.GravitySource = true;
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

            // Exact Solution manufactured Solution
            Func<double[], double, double> VelocityXfunction = (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> VelocityYfunction = (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]);
            Func<double[], double, double> Pressurefunction = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]);
            Func<double[], double, double> StressXXfunction = (X, t) => -2 * (1 - C.beta) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> StressXYfunction = (X, t) => -2 * (1 - C.beta) * Math.Exp(X[0]) * (Math.Cos(X[1]) - X[1] * Math.Sin(X[1]));
            Func<double[], double, double> StressYYfunction = (X, t) => 2 * (1 - C.beta) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));

            //PROBABLY ALWAYS ZERO?
            C.GravityX = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]) * ((C.Reynolds + C.beta - 1) / C.Reynolds);
            C.GravityY = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.Reynolds + C.beta - 1) / C.Reynolds);

            // Insert Exact Solution
            C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            C.ExSol_Pressure = Pressurefunction;
            C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            // Create Fields
            C.deg = PolyDeg;
            C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = PolyDeg - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Phi", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });


            // Create Grid
            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CellType.Square_Linear, C.FixedStreamwisePeriodicBC);

                if (!C.FixedStreamwisePeriodicBC) {
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

                    if (!C.FixedStreamwisePeriodicBC) {
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
            if (C.SetParamsAnalyticalSol == true) {
                C.VelFunctionU = X => VelocityXfunction(X, 0);
                C.VelFunctionV = X => VelocityYfunction(X, 0);
            }

            // Set Initial Conditions
            if (C.SetInitialConditions == true) {

                C.InitialValues_Evaluators.Add("VelocityX", X => VelocityXfunction(X, 0));
                C.InitialValues_Evaluators.Add("VelocityY", X => VelocityYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXX", X => StressXXfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressXY", X => StressXYfunction(X, 0));
                C.InitialValues_Evaluators.Add("StressYY", X => StressYYfunction(X, 0));
                C.InitialValues_Evaluators.Add("GravityX", X => C.GravityX(X, 0));
                C.InitialValues_Evaluators.Add("GravityY", X => C.GravityY(X, 0));

                if (C.SetInitialPressure == true || C.SkipSolveAndEvaluateResidual == true) {
                    C.InitialValues_Evaluators.Add("Pressure", X => Pressurefunction(X, 0));

                }
            }

            C.InitialValues_Evaluators.Add("Phi", X => -1);

            // Set Boundary Conditions
            if (!C.FixedStreamwisePeriodicBC) {
                C.AddBoundaryValue("Velocity_inlet", "VelocityX", VelocityXfunction);
                C.AddBoundaryValue("Velocity_inlet", "VelocityY", VelocityYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXX", StressXXfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressXY", StressXYfunction);
                C.AddBoundaryValue("Velocity_inlet", "StressYY", StressYYfunction);
                C.AddBoundaryValue("Velocity_inlet", "Pressure", Pressurefunction);
                C.AddBoundaryValue("Velocity_inlet", "GravityX", C.GravityX);
                C.AddBoundaryValue("Velocity_inlet", "GravityY", C.GravityY);
                //C.AddBoundaryCondition("Pressure_Outlet");

            }
            return C;
        }       
    }
}
