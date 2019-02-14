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
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Application.Rheology {
    static public class RheologyConsistencyTest {
 
        //Consistency constitutive equation with convective term

        static public RheologyControl ConsistencyConstitutiveComputeRes(int GridRes = 3, int PolyDeg = 2, double beta = 0)
        {
            RheologyControl C = ConsistencyConstitutiveGeneric(GridRes, PolyDeg, beta);

            Console.WriteLine("Test 1: Insert exact polynomial solution and only compute residual.");
            C.Stokes = false;
            C.SkipSolveAndEvaluateResidual = true;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.beta = beta;

            return C;
        }


        static public RheologyControl ConsistencyConstitutiveStability(int GridRes = 2, int PolyDeg = 2, double beta = 0)
        {
            RheologyControl C = ConsistencyConstitutiveGeneric(GridRes, PolyDeg, beta);

            Console.WriteLine("Test 2: Insert exact polynomial solution and test if nonlinear solver (NEWTON) stays stable.");
            C.Stokes = false;
            C.SkipSolveAndEvaluateResidual = false;
            C.SetInitialConditions = true;
            C.SetInitialPressure = false;
            C.SetParamsAnalyticalSol = false;
            C.beta = beta;


            return C;
        }


        static public RheologyControl ConsistencyConstitutiveGeneric(int GridRes = 2, int PolyDeg = 2, double beta = 0)
        {


            RheologyControl C = new RheologyControl();

            // Solver Options
            C.NoOfTimesteps = 1;
            C.savetodb = false;
            //C.DbPath = "C:\AnnesBoSSSdb\ConsistencyConstitutive_withDiv";
            C.SessionName = "Degree" + PolyDeg + ", GridLevel" + GridRes;
            C.ProjectName = "ConsistencyStudyConstitutive";
            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.NonLinearSolver.ConvergenceCriterion = 1E-20;
            C.LinearSolver.MaxSolverIterations = 20;
            C.LinearSolver.MinSolverIterations = 3;
            C.LinearSolver.ConvergenceCriterion = 1E-13;


            //C.MaxIter = 20;
            //C.MinIter = 3;
            //C.ConvCrit = 1E-20;
            //C.ConvCritGMRES = 1E-13;
            C.dt = 1E20;
            C.dtMax = C.dt;
            C.dtMin = C.dt;
            C.Timestepper_Scheme = RheologyControl.TimesteppingScheme.ImplicitEuler;
            C.NonLinearSolver.SolverCode = NonLinearSolverConfig.Code.Newton;//C.NonlinearMethod = NonlinearSolverMethod.Newton;

            //Grid Params
            //double GridLevel = 5;
            double h = Math.Pow(2, -GridRes + 1);
            double cells = 1 / h;
            int cells2 = (int)cells;

            //Debugging and Solver Analysis
            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = true;
            C.SetInitialConditions = true;
            C.SetInitialPressure = true;
            C.SetParamsAnalyticalSol = false;
            C.ComputeL2Error = true;

            //Physical Params
            C.Stokes = false;
            C.FixedStreamwisePeriodicBC = false;
            C.GravitySource = true;
            C.beta = 0;
            C.Reynolds = 1;
            C.Weissenberg = 1;

            //Penalties
            C.ViscousPenaltyScaling = 1;
            C.Penalty2 = 1;
            C.Penalty1[0] = 0.0;
            C.Penalty1[1] = 0.0;
            C.PresPenalty2 = 1;
            C.PresPenalty1[0] = 0.0;
            C.PresPenalty1[1] = 0.0;

            // Exact Solution manufactured Solution
            Func<double[], double, double> VelocityXfunction = (X, t) => X[0] * X[0];
            Func<double[], double, double> VelocityYfunction = (X, t) => -X[1];
            Func<double[], double, double> Pressurefunction = (X, t) => 0;
            Func<double[], double, double> StressXXfunction = (X, t) => X[0] * X[0];
            Func<double[], double, double> StressXYfunction = (X, t) => X[0] * X[0] + X[1] * X[1];
            Func<double[], double, double> StressYYfunction = (X, t) => X[1] * X[1];

            //Gravity sources
            //Weissenberg = 1 including Objective Terms!
            C.GravityX = (X, t) => 2 * X[0] * X[0] * X[0] + X[0] * X[0] * (2 * X[0] - 1) - 2 * X[0] - 2 * X[1];
            C.GravityY = (X, t) => -X[1] - X[1] * (2 * X[0] - 1) - 2 * X[0];
            C.GravityXX = (X, t) => X[0] * X[0] - 2 * X[0] * X[0] * X[0] - 4 * X[0];
            C.GravityXY = (X, t) => X[0] * X[0] - X[1] * X[1] + 2 * X[0] * X[0] * X[0] - (2 * X[0] - 1) * (X[0] * X[0] + X[1] * X[1]);
            C.GravityYY = (X, t) => X[1] * X[1] + 2;
            C.GravityDiv = (X, t) => 2 * X[0] - 1;

            // Insert Exact Solution
            C.ExSol_Velocity = new Func<double[], double, double>[] { VelocityXfunction, VelocityYfunction };
            C.ExSol_Pressure = Pressurefunction;
            C.ExSol_Stress = new Func<double[], double, double>[] { StressXXfunction, StressXYfunction, StressYYfunction };

            // Create Fields
            //int degree = 2;
            C.FieldOptions.Add("VelocityX", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("VelocityY", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("Pressure", new FieldOpts() { Degree = PolyDeg - 1, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });

            C.FieldOptions.Add("StressXX", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressXY", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            C.FieldOptions.Add("StressYY", new FieldOpts() { Degree = PolyDeg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
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
                C.InitialValues_Evaluators.Add("GravityX", X => C.GravityX(X, 0));
                C.InitialValues_Evaluators.Add("GravityY", X => C.GravityY(X, 0));
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
    }
}
