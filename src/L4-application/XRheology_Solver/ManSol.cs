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
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools.TestCases;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XRheology_Solver {

    /// <summary>
    /// class providing Controls for the channel flow type testcases
    /// </summary>
    public static class ManSol {


        /// <summary>
        /// control object for various testing
        /// </summary>
        /// <returns></returns>
        public static XRheology_Control ManSol_Consistency(int p = 2, int kelem = 16, int wallBC = 0) {

            XRheology_Control C = new XRheology_Control();

            string _DbPath = null; // @"D:\local\local_test_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XRheology/Channel";
            C.ProjectDescription = "Consistency test with non-polynomial manufactured solution";

            C.ContinueOnIoError = false;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            }); C.FieldOptions.Add("StressXX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("StressXY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("StressYY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Gravity", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            //C.FieldOptions.Add("GravityXX", new FieldOpts() {
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("GravityXY", new FieldOpts() {
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("GravityYY", new FieldOpts() {
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            //C.FieldOptions.Add("KineticEnergy", new FieldOpts() {
            //    Degree = 2*p,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.reynolds_A = 1.0;
            C.PhysicalParameters.reynolds_B = 1.0;
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Weissenberg_a = 0.0;
            C.PhysicalParameters.Weissenberg_b = 0.0;
            C.PhysicalParameters.beta_a = 1.0;
            C.PhysicalParameters.beta_b = 1.0;
            double sigma = 0.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.beta_S = 0.05;
            //C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-1, 1, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);// periodicX:true, periodicY:true);

                //grd.EdgeTagNames.Add(2, "Freeslip");
                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                //grd.EdgeTagNames.Add(2, "Velocity_inlet_right");
                //grd.EdgeTagNames.Add(3, "Pressure_Outlet");


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

                    if (Math.Abs(x - (-1)) < 1.0e-6)
                        // left
                        return 1;

                    if (Math.Abs(x - (+1)) < 1.0e-6)
                        // right
                        return 1;

                    throw new ArgumentOutOfRangeException();
                });

                return grd;
            };

            #endregion

            // Functions for exact solution
            // ==============
            #region func

            Func<double[], double, double> PhiFunc = (X,t) => X[1] - (1 / 2.0) + 0.2; // + (H/20)*Math.Cos(8 * Math.PI * X[0] / L)); //-1);

            Func<double[], double, double> VelX_A_Func = (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1])); //X[0] * X[0];
            Func<double[], double, double> VelY_A_Func = (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]); //-2 * X[0] * X[1];

            Func<double[], double, double> VelX_B_Func = (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));// X[0] * X[0];
            Func<double[], double, double> VelY_B_Func = (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]); //-2 * X[0] * X[1];

            Func<double[], double, double> Pres_A_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]); //2 * X[0];
            Func<double[], double, double> Pres_B_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]); //2 * X[0]; 

            Func<double[], double, double> StressXX_A_Func = (X, t) => -2 * (1 - C.PhysicalParameters.beta_a) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> StressXY_A_Func = (X, t) => -2 * (1 - C.PhysicalParameters.beta_a) * Math.Exp(X[0]) * (Math.Cos(X[1]) - X[1] * Math.Sin(X[1]));
            Func<double[], double, double> StressYY_A_Func = (X, t) => 2 * (1 - C.PhysicalParameters.beta_a) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));

            Func<double[], double, double> StressXX_B_Func = (X, t) => -2 * (1 - C.PhysicalParameters.beta_b) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> StressXY_B_Func = (X, t) => -2 * (1 - C.PhysicalParameters.beta_b) * Math.Exp(X[0]) * (Math.Cos(X[1]) - X[1] * Math.Sin(X[1]));
            Func<double[], double, double> StressYY_B_Func = (X, t) => 2 * (1 - C.PhysicalParameters.beta_b) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));

            Func<double[], double, double> GravityX_A_Func = (X, t) => -1 / C.PhysicalParameters.reynolds_A * Math.Exp(X[0]) * (Math.Exp(X[0]) * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.reynolds_A - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_A * X[1] * X[1] - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_A - 2 * Math.Sin(X[1]) * C.PhysicalParameters.reynolds_A + 2 * Math.Sin(X[1]));
            Func<double[], double, double> GravityY_A_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.PhysicalParameters.reynolds_A - 1) / C.PhysicalParameters.reynolds_A);

            Func<double[], double, double> GravityX_B_Func = (X, t) => -1 / C.PhysicalParameters.reynolds_B * Math.Exp(X[0]) * (Math.Exp(X[0]) * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.reynolds_B - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_B * X[1] * X[1] - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_B - 2 * Math.Sin(X[1]) * C.PhysicalParameters.reynolds_B + 2 * Math.Sin(X[1]));
            Func<double[], double, double> GravityY_B_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.PhysicalParameters.reynolds_B - 1) / C.PhysicalParameters.reynolds_B);

            Func<double[], double, double> GravityXX_A_Func = (X, t) => 2 * Math.Exp(2 * X[0]) * C.PhysicalParameters.Weissenberg_a * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParameters.beta_a * X[1] * C.PhysicalParameters.Weissenberg_a
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.beta_a + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] + C.PhysicalParameters.beta_a * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) - X[1] * X[1] + C.PhysicalParameters.beta_a - 1);
            Func<double[], double, double> GravityXY_A_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * (C.PhysicalParameters.beta_a - 1) * C.PhysicalParameters.Weissenberg_a * (2 * Math.Pow(Math.Cos(X[1]), 2) * X[1] 
                + 3 * Math.Cos(X[1]) * Math.Sin(X[1]) + X[1]);
            Func<double[], double, double> GravityYY_A_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * C.PhysicalParameters.Weissenberg_a * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParameters.beta_a * X[1]
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.beta_a + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] - 3 * C.PhysicalParameters.beta_a * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) 
                + 3 * X[1] * X[1] - 3 * C.PhysicalParameters.beta_a + 3);

            Func<double[], double, double> GravityXX_B_Func = (X, t) => 2 * Math.Exp(2 * X[0]) * C.PhysicalParameters.Weissenberg_b * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParameters.beta_b * X[1] * C.PhysicalParameters.Weissenberg_b
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.beta_b + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] + C.PhysicalParameters.beta_b * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) - X[1] * X[1] + C.PhysicalParameters.beta_b - 1);
            Func<double[], double, double> GravityXY_B_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * (C.PhysicalParameters.beta_b - 1) * C.PhysicalParameters.Weissenberg_b * (2 * Math.Pow(Math.Cos(X[1]), 2) * X[1]
                + 3 * Math.Cos(X[1]) * Math.Sin(X[1]) + X[1]);
            Func<double[], double, double> GravityYY_B_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * C.PhysicalParameters.Weissenberg_b * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParameters.beta_b * X[1]
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.beta_b + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] - 3 * C.PhysicalParameters.beta_b * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1])
                + 3 * X[1] * X[1] - 3 * C.PhysicalParameters.beta_b + 3);

            #endregion

            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X,0));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => VelX_A_Func(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX#B", X => VelX_B_Func(X, 0));

            C.InitialValues_Evaluators.Add("VelocityY#A", X => VelY_A_Func(X, 0));
            C.InitialValues_Evaluators.Add("VelocityY#B", X => VelY_B_Func(X, 0));

            C.InitialValues_Evaluators.Add("Pressure#A", X => Pres_A_Func(X, 0));
            C.InitialValues_Evaluators.Add("Pressure#B", X => Pres_B_Func(X, 0));
            //double Pjump = sigma / radius;
            //C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);

            C.InitialValues_Evaluators.Add("GravityX#A", X => GravityX_A_Func(X, 0));
            C.InitialValues_Evaluators.Add("GravityX#B", X => GravityX_B_Func(X, 0));

            C.InitialValues_Evaluators.Add("GravityY#A", X => GravityY_A_Func(X, 0));
            C.InitialValues_Evaluators.Add("GravityY#B", X => GravityY_B_Func(X, 0));

            //C.InitialValues_Evaluators.Add("GravityXX#A", X => GravityXX_A_Func(X, 0));
            //C.InitialValues_Evaluators.Add("GravityXX#B", X => GravityXX_B_Func(X, 0));

            //C.InitialValues_Evaluators.Add("GravityXY#A", X => GravityXY_A_Func(X, 0));
            //C.InitialValues_Evaluators.Add("GravityXY#B", X => GravityXY_B_Func(X, 0));

            //C.InitialValues_Evaluators.Add("GravityYY#A", X => GravityYY_A_Func(X, 0));
            //C.InitialValues_Evaluators.Add("GravityYY#B", X => GravityYY_B_Func(X, 0));

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("cf6bd7bf-a19f-409e-b8c2-0b89388daad6");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 10);

            #endregion

            // exact solution
            // ==============
            #region exact

            C.Phi = PhiFunc;

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { VelX_A_Func, VelY_A_Func });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { VelX_B_Func, VelY_B_Func });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", Pres_A_Func);
            //C.ExactSolutionPressure.Add("B", Pres_B_Func);

            //C.ExactSolutionStressXX = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionStressXX.Add("A", StressXX_A_Func);
            //C.ExactSolutionStressXX.Add("B", StressXX_B_Func);

            //C.ExactSolutionStressXY = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionStressXY.Add("A", StressXY_A_Func);
            //C.ExactSolutionStressXY.Add("B", StressXY_B_Func);

            //C.ExactSolutionStressYY = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionStressYY.Add("A", StressYY_A_Func);
            //C.ExactSolutionStressYY.Add("B", StressYY_B_Func);

            #endregion

            // gravity source
            // ===================
            #region gravity
            C.GravityX = new Dictionary<string, Func<double[], double, double>>();
            C.GravityX.Add("A", GravityX_A_Func);
            C.GravityX.Add("B", GravityX_B_Func);

            C.GravityY = new Dictionary<string, Func<double[], double, double>>();
            C.GravityY.Add("A", GravityY_A_Func);
            C.GravityY.Add("B", GravityY_B_Func);

            //C.GravityXX = new Dictionary<string, Func<double[], double, double>>();
            //C.GravityXX.Add("A", GravityXX_A_Func);
            //C.GravityXX.Add("B", GravityXX_B_Func);

            //C.GravityXY = new Dictionary<string, Func<double[], double, double>>();
            //C.GravityXY.Add("A", GravityXY_A_Func);
            //C.GravityXY.Add("B", GravityXY_B_Func);

            //C.GravityYY = new Dictionary<string, Func<double[], double, double>>();
            //C.GravityYY.Add("A", GravityYY_A_Func);
            //C.GravityYY.Add("B", GravityYY_B_Func);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            //C.AddBoundaryValue("Freeslip");
            C.AddBoundaryValue("Velocity_inlet", "VelocityX#A", VelX_A_Func);
            C.AddBoundaryValue("Velocity_inlet", "VelocityX#B", VelX_B_Func);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY#A", VelY_A_Func);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY#B", VelY_B_Func);

            C.AddBoundaryValue("Velocity_inlet", "Pressure#A", Pres_A_Func);
            C.AddBoundaryValue("Velocity_inlet", "Pressure#B", Pres_B_Func);

            //C.AddBoundaryValue("Velocity_inlet_right", "VelocityX#A", (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1])));
            //C.AddBoundaryValue("Velocity_inlet_right", "VelocityX#B", (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1])));
            //C.AddBoundaryValue("Velocity_inlet_right", "VelocityY#A", (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]));
            //C.AddBoundaryValue("Velocity_inlet_right", "VelocityY#B", (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]));

            //C.AddBoundaryValue("Velocity_inlet_right", "Pressure#A", (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]));
            //C.AddBoundaryValue("Velocity_inlet_right", "Pressure#B", (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]));

            //C.AddBoundaryValue("Pressure_Outlet");

            //C.AddBoundaryValue("Velocity_inlet", "GravityX#A", GravityX_A_Func);
            //C.AddBoundaryValue("Velocity_inlet", "GravityX#B", GravityX_B_Func);
            //C.AddBoundaryValue("Velocity_inlet", "GravityY#A", GravityY_A_Func);
            //C.AddBoundaryValue("Velocity_inlet", "GravityY#B", GravityY_B_Func);


            //C.AddBoundaryValue("Velocity_inlet_right", "GravityX#A", GravityX_A_Func);
            //C.AddBoundaryValue("Velocity_inlet_right", "GravityX#B", GravityX_B_Func);
            //C.AddBoundaryValue("Velocity_inlet_right", "GravityY#A", GravityY_A_Func);
            //C.AddBoundaryValue("Velocity_inlet_right", "GravityY#B", GravityY_B_Func);

            //C.AddBoundaryValue("Velocity_inlet", "GravityXX#A", GravityXX_A_Func);
            //C.AddBoundaryValue("Velocity_inlet", "GravityXX#B", GravityXX_B_Func);

            //C.AddBoundaryValue("Velocity_inlet", "GravityXY#A", GravityXY_A_Func);
            //C.AddBoundaryValue("Velocity_inlet", "GravityXY#B", GravityXY_B_Func);

            //C.AddBoundaryValue("Velocity_inlet", "GravityYY#A", GravityYY_A_Func);
            //C.AddBoundaryValue("Velocity_inlet", "GravityYY#B", GravityYY_B_Func);

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_mumps;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.Penalty2 = 1;
            C.SkipSolveAndEvaluateResidual = false;


            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XRheology_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.CompMode = AppControl._CompMode.Steady;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 300;
            C.saveperiod = 10;

            #endregion


            return C;
        }



        // ==========================================
        // Control-objects for elementalTestProgramm
        // ==========================================


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XRheology_Control CF_BoundaryTest(int bc = 3, bool Xperiodic = true, bool init_exact = false, bool slip = false) {

            int p = 2;
            int kelem = 16;


            XRheology_Control C = new XRheology_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null; //_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Channel flow for BC testing";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 2;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: Xperiodic);

                switch (bc) {
                    case 1: {
                            grd.EdgeTagNames.Add(1, "freeslip_lower");
                            grd.EdgeTagNames.Add(2, "freeslip_upper");
                            break;
                        }
                    case 2: {
                            grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                            grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                            break;
                        }
                    case 3: {
                            grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                            grd.EdgeTagNames.Add(2, "freeslip_upper");
                            break;
                        }
                    case 4: {
                            grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                            grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                            break;
                        }
                    default: {
                            throw new NotImplementedException("No such testcase available");
                        }
                }

                if (!Xperiodic) {
                    grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                    grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (!Xperiodic) {
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 1.0;

            if (init_exact) 
                C.InitialValues_Evaluators.Add("VelocityX#A", X => U);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            switch (bc) {
                case 1: {
                        C.AddBoundaryValue("freeslip_lower");
                        if (slip)
                            C.AddBoundaryValue("freeslip_upper", "VelocityX#A", X => -U);
                        else
                            C.AddBoundaryValue("freeslip_upper");

                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                case 2: {
                        //C.AddBoundaryCondition("navierslip_linear_lower");
                        //if (slip)
                        //    C.AddBoundaryCondition("navierslip_linear_upper", "VelocityX#A", X => -U);
                        //else
                        //    C.AddBoundaryCondition("navierslip_linear_upper");

                        C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", X => -U);
                        C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#A", X => U);

                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                case 3: {
                        C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                        C.AddBoundaryValue("freeslip_upper");
                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                case 4: {
                        C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                        C.AddBoundaryValue("navierslip_linear_upper");
                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                default: {
                        break;
                    }
            }

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations=50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XRheology_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.CompMode = AppControl._CompMode.Steady;

            double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 100;
            C.saveperiod = 1;

            #endregion


            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XRheology_Control CF_LevelSetMovementTest(int boundarySetup = 2, double characteristicLength = 1.0, LevelSetEvolution lsEvo = LevelSetEvolution.FastMarching, 
            LevelSetHandling lsHandl = LevelSetHandling.Coupled_Once, XRheology_Control.TimesteppingScheme tsScheme = XRheology_Control.TimesteppingScheme.ImplicitEuler) {

            int p = 2;
            int kelem = 16;
            double cLength = characteristicLength;

            XRheology_Control C = new XRheology_Control();
            
            // basic database options
            // ======================
            #region db

            C.DbPath = null; //_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Two-phase Channel flow for testing the level set movement";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 2;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);

                bool xPeriodic = (boundarySetup == 1) ? true : false;
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);

                grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");

                switch (boundarySetup) {
                    case 1:
                        grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                        grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                        break;
                    case 2:
                        grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                        grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                        break;
                    default:
                        throw new ArgumentException("invalid boundary setup");

                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic) {
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc;

            switch (boundarySetup) {
                case 1: {
                        // horizontal interface
                        PhiFunc = (X => ((X[0]- cLength).Pow2()).Sqrt() - cLength/2);
                        break;
                    }
                case 2: {
                        // radial interface
                        double[] center = new double[] { L / 4.0, H / 2.0 };
                        double radius = cLength;
                        PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius);
                        break;
                    }
                default:
                    PhiFunc = (X => -1);
                    break;
                
            }
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 1.0;

            switch (boundarySetup) {
                case 1:
                    //C.InitialValues_Evaluators.Add("VelocityY#A", X => U);
                    //C.InitialValues_Evaluators.Add("VelocityY#B", X => U);
                    C.InitialValues_Evaluators.Add("VelocityX#A", X => U);
                    C.InitialValues_Evaluators.Add("VelocityX#B", X => U);
                    break;
                case 2:
                    C.InitialValues_Evaluators.Add("VelocityX#A", X => U);
                    C.InitialValues_Evaluators.Add("VelocityX#B", X => U);
                    break;
                default:
                    throw new ArgumentException("invalid boundary setup");

            }


            #endregion


            // boundary conditions
            // ===================
            #region BC

            switch (boundarySetup) {
                case 1:
                    //C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", X => U);
                    //C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#B", X => U);
                    //C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#A", X => U);
                    //C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => U);
                    C.AddBoundaryValue("pressure_outlet_right");
                    break;
                case 2:
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => U);
                    C.AddBoundaryValue("pressure_outlet_right");
                    break;
                default:
                    break;
            }

            #endregion

            // advanced settings for Fourier-Level-Set
            // ======================
            #region Fourier level-set

            switch (lsEvo)
            {
                case LevelSetEvolution.Fourier:
                    {
                        switch (boundarySetup)
                        {
                            case 1:
                                {
                                    throw new ArgumentException("Fourier Level-Set not implemented in Line Movement Test");
                                    break;
                                }
                            case 2:
                                {
                                    int numSp = 640;
                                    double[] FourierP = new double[numSp];
                                    double[] samplP = new double[numSp];
                                    double[] center = new double[] { L / 4.0, H / 2.0 };
                                    double radius = cLength;
                                    for (int sp = 0; sp < numSp; sp++)
                                    {
                                        FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                                        samplP[sp] = radius;
                                    }

                                    C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem)
                                    {
                                        center = center,
                                        FourierEvolve = Fourier_Evolution.MaterialPoints,
                                        centerMove = CenterMovement.Reconstructed,
                                    };

                                    C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
                                    break;
                                }                               
                            default:
                                break;
                               
                        }
                        break;
                    }
                default:
                    break;
            }

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 4;
            C.LinearSolver.MinSolverIterations = 4;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.Option_LevelSetEvolution = lsEvo;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            C.Timestepper_LevelSetHandling = lsHandl;
            C.Timestepper_Scheme = tsScheme;

            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;
            C.saveperiod = 1;

            #endregion


            return C;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XRheology_Control CF_LevelSetRotationTest(int boundarySetup = 1, double characteristicLength = 1.0, LevelSetEvolution lsEvo = LevelSetEvolution.FastMarching,
            LevelSetHandling lsHandl = LevelSetHandling.Coupled_Once, XRheology_Control.TimesteppingScheme tsScheme = XRheology_Control.TimesteppingScheme.ImplicitEuler)
        {

            int p = 2;
            int kelem = 16;
            double cLength = characteristicLength;

            XRheology_Control C = new XRheology_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null; //_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Two-phase flow for testing the level set movement in solid body rotation";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts()
            {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts()
            {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 1;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L/2, L/2, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-H/2, H/2, kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");                   
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "velocity_inlet_right");               

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + H/2) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H/2) <= 1.0e-8)
                        et = 2;                    
                    if (Math.Abs(X[0] + L/2) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L/2) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc;

            switch (boundarySetup)
            {
                case 1:
                    {
                        // elliptoid
                        double[] center = new double[] { 0.15, 0.0 };
                        double[] shape = new double[] { 1, 0.36 };
                        double radius = cLength;
                        PhiFunc = (X => ((X[0] - center[0]).Pow2() / shape[0] + (X[1] - center[1]).Pow2() / shape[1]).Sqrt() - radius);
                        break;
                    }
                case 2:
                    {
                        // slotted disk
                        double[] xCutout = new double[] { -0.1, 0.1 };
                        double yCutout = -0.1;
                        double radius = cLength;
                        ZalesaksDisk disk = new ZalesaksDisk(xCutout, yCutout, radius);
                        PhiFunc = (X => disk.SignedDistanceLevelSet(X));
                        break;
                    }
                default:
                    PhiFunc = (X => -1);
                    break;

            }
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("VelocityX#A", X => -X[1]);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => -X[1]);
            C.InitialValues_Evaluators.Add("VelocityY#A", X => X[0]);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => X[0]);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#B", X => -X[1]);

            C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#B", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#B", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#B", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityY#B", X => X[0]);

            #endregion

            // advanced settings for Fourier-Level-Set
            // ======================
            #region Fourier level-set

            switch (lsEvo)
            {
                case LevelSetEvolution.Fourier:
                    {
                        switch (boundarySetup)
                        {
                            case 1:
                                {
                                    int numSp = 640;
                                    double[] FourierP = new double[numSp];
                                    double[] samplP = new double[numSp];
                                    double[] center = new double[] { 0.15, 0.0 };
                                    double radius = cLength;
                                    for (int sp = 0; sp < numSp; sp++)
                                    {
                                        FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                                        samplP[sp] = radius / (Math.Cos(FourierP[sp]).Pow2() + Math.Sin(FourierP[sp]).Pow2() / 0.36).Sqrt();
                                    }

                                    C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem)
                                    {
                                        center = center,
                                        FourierEvolve = Fourier_Evolution.MaterialPoints,
                                        centerMove = CenterMovement.Reconstructed,
                                        PeriodicFunc = (X => radius / (Math.Cos(X).Pow2() + Math.Sin(X).Pow2() / 0.36).Sqrt())
                                    };

                                    C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
                                    break;
                                }
                            case 2:
                                {
                                    throw new ArgumentException("Fourier Level-Set is not suitable for Slotted Disk");
                                }
                            default:
                                break;

                        }
                        break;
                    }
                default:
                    break;
            }

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 4;
            C.LinearSolver.MinSolverIterations = 4;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.Option_LevelSetEvolution = lsEvo;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            C.Timestepper_LevelSetHandling = lsHandl;
            C.Timestepper_Scheme = tsScheme;

            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;
            C.saveperiod = 1;

            #endregion


            return C;
        }

    }
}
