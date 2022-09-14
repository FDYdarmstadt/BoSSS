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
using BoSSS.Solution.LevelSetTools;

namespace BoSSS.Application.XRheology_Solver {

    /// <summary>
    /// class providing Controls for unphysical tests
    /// </summary>
    public static class GenericTestcases {


        /// <summary>
        /// control object for quasi 1d test where one element border is replaced by levelset. The pointwise error for each flux can be compared.
        /// </summary>
        /// <returns></returns>
        public static XRheology_Control InterfaceTest(int p = 2) {

            XRheology_Control C = new XRheology_Control();

            string _DbPath = null; // @"D:\local\local_test_db";

            // basic options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XRheology/InterfaceTest";
            C.ProjectDescription = "Cell border replaced by vertical interface";
            
            C.OperatorMatrixAnalysis = true;
            C.SkipSolveAndEvaluateResidual = true;
            bool PHI = true;
            C.InterfaceTest = true;

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
            });
            C.FieldOptions.Add("StressXX", new FieldOpts() {
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

            C.PhysicalParametersRheology.reynolds_A = 1.0;
            C.PhysicalParametersRheology.reynolds_B = 1.0;
            C.PhysicalParametersRheology.beta_a = 0.0;
            C.PhysicalParametersRheology.beta_b = 0.0;

            C.RaiseWeissenberg = false;
            C.PhysicalParametersRheology.Weissenberg_a = 1.0;
            C.PhysicalParametersRheology.Weissenberg_b = 1.0;
            C.WeissenbergIncrement = 0.1;

            double sigma = 0.0;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            C.GridFunc = delegate () {

                double[] Xnodes;

                if (PHI == true) {
                    //with interface at 0
                    Xnodes = new double[] { -2, -1, 1, 2 };
                    Console.WriteLine("PHI is inserted at x = 0 instead of inner edge between cells");
                } else {
                    //without interface at 0
                    Xnodes = new double[] { -2, -1, 0, 1, 2 };
                }

                double[] Ynodes = new double[] { 0, 1};
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);
                //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "velocity_inlet");
                grd.EdgeTagNames.Add(4, "pressure_outlet");
                grd.EdgeTagNames.Add(2, "wall");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - 1) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[1] ) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] - 2) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + 2) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            // functions
            // ===========
            #region function definition

            Func<double[], double> PhiFunc;

            if (PHI == true) {
                //with interface at 0
                PhiFunc = (X => X[0]);
            } else {
                //without interface at 0
                PhiFunc = (X => -1);
            }

            Func<double[], double, double> VelocityXfunction_A = (X, t) => (X[0] * X[0]);  
            Func<double[], double, double> VelocityXfunction_B = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> VelocityYfunction_A = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> VelocityYfunction_B = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> Pressurefunction_A = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> Pressurefunction_B = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> StressXXfunction_A = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> StressXXfunction_B = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> StressXYfunction_A = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> StressXYfunction_B = (X, t) => (X[0] * X[0]);
            Func<double[], double, double> StressYYfunction_A = (X, t) => (X[0] * X[0]); 
            Func<double[], double, double> StressYYfunction_B = (X, t) => (X[0] * X[0]);
            #endregion

            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            #endregion

            // exact solution
            // ==============
            #region exact

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { VelocityXfunction_A, VelocityYfunction_A });
            C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { VelocityXfunction_B, VelocityYfunction_B });

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", Pressurefunction_A);
            C.ExactSolutionPressure.Add("B", Pressurefunction_B);

            C.ExactSolutionStressXX = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionStressXX.Add("A", (X, t) => StressXXfunction_A(X, t));// * C.PhysicalParameters.Weissenberg_a);
            C.ExactSolutionStressXX.Add("B", (X, t) => StressXXfunction_B(X, t));// * C.PhysicalParameters.Weissenberg_b);

            C.ExactSolutionStressXY = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionStressXY.Add("A", StressXYfunction_A);
            C.ExactSolutionStressXY.Add("B", StressXYfunction_B);

            C.ExactSolutionStressYY = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionStressYY.Add("A", StressYYfunction_A);
            C.ExactSolutionStressYY.Add("B", StressYYfunction_B);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("Velocity_inlet", "VelocityX#A", VelocityXfunction_A);
            C.AddBoundaryValue("Velocity_inlet", "VelocityX#B", VelocityXfunction_B);

            C.AddBoundaryValue("Velocity_inlet", "VelocityY#A", VelocityYfunction_A);
            C.AddBoundaryValue("Velocity_inlet", "VelocityY#B", VelocityYfunction_B);

            C.AddBoundaryValue("Velocity_inlet", "StressXX#A", StressXXfunction_A);
            C.AddBoundaryValue("Velocity_inlet", "StressXX#B", StressXXfunction_B);

            C.AddBoundaryValue("Velocity_inlet", "StressXY#A", StressXYfunction_A);
            C.AddBoundaryValue("Velocity_inlet", "StressXY#B", StressXYfunction_B);

            C.AddBoundaryValue("Velocity_inlet", "StressYY#A", StressYYfunction_A);
            C.AddBoundaryValue("Velocity_inlet", "StressYY#B", StressYYfunction_B);

            C.AddBoundaryValue("pressure_outlet");
            C.AddBoundaryValue("wall");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 3;
            C.NonLinearSolver.MinSolverIterations = 1;     
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.Penalty2 = 1;
            C.AdvancedDiscretizationOptions.Penalty1[0] = 0;
            C.AdvancedDiscretizationOptions.Penalty1[1] = 0;
            //C.AdvancedDiscretizationOptions.PresPenalty2 = 0.0;

            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 300;
            C.saveperiod = 10;

            #endregion


            return C;
        }

        /// <summary>
        /// control object for various testing
        /// </summary>
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

            C.PhysicalParametersRheology.reynolds_A = 1.0;
            C.PhysicalParametersRheology.reynolds_B = 1.0;
            C.PhysicalParametersRheology.rho_A = 1;
            C.PhysicalParametersRheology.rho_B = 1;
            C.PhysicalParametersRheology.mu_A = 1;
            C.PhysicalParametersRheology.mu_B = 1;
            C.PhysicalParametersRheology.Weissenberg_a = 0.0;
            C.PhysicalParametersRheology.Weissenberg_b = 0.0;
            C.PhysicalParametersRheology.beta_a = 1.0;
            C.PhysicalParametersRheology.beta_b = 1.0;
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

            Func<double[], double, double> PhiFunc = (X, t) => X[1] - (1 / 2.0) + 0.2; // + (H/20)*Math.Cos(8 * Math.PI * X[0] / L)); //-1);

            Func<double[], double, double> VelX_A_Func = (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1])); //X[0] * X[0];
            Func<double[], double, double> VelY_A_Func = (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]); //-2 * X[0] * X[1];

            Func<double[], double, double> VelX_B_Func = (X, t) => -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));// X[0] * X[0];
            Func<double[], double, double> VelY_B_Func = (X, t) => Math.Exp(X[0]) * X[1] * Math.Sin(X[1]); //-2 * X[0] * X[1];

            Func<double[], double, double> Pres_A_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]); //2 * X[0];
            Func<double[], double, double> Pres_B_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Sin(X[1]); //2 * X[0]; 

            Func<double[], double, double> StressXX_A_Func = (X, t) => -2 * (1 - C.PhysicalParametersRheology.beta_a) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> StressXY_A_Func = (X, t) => -2 * (1 - C.PhysicalParametersRheology.beta_a) * Math.Exp(X[0]) * (Math.Cos(X[1]) - X[1] * Math.Sin(X[1]));
            Func<double[], double, double> StressYY_A_Func = (X, t) => 2 * (1 - C.PhysicalParametersRheology.beta_a) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));

            Func<double[], double, double> StressXX_B_Func = (X, t) => -2 * (1 - C.PhysicalParametersRheology.beta_b) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            Func<double[], double, double> StressXY_B_Func = (X, t) => -2 * (1 - C.PhysicalParametersRheology.beta_b) * Math.Exp(X[0]) * (Math.Cos(X[1]) - X[1] * Math.Sin(X[1]));
            Func<double[], double, double> StressYY_B_Func = (X, t) => 2 * (1 - C.PhysicalParametersRheology.beta_b) * Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));

            Func<double[], double, double> GravityX_A_Func = (X, t) => -1 / C.PhysicalParameters.reynolds_A * Math.Exp(X[0]) * (Math.Exp(X[0]) * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.reynolds_A - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_A * X[1] * X[1] - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_A - 2 * Math.Sin(X[1]) * C.PhysicalParameters.reynolds_A + 2 * Math.Sin(X[1]));
            Func<double[], double, double> GravityY_A_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.PhysicalParameters.reynolds_A - 1) / C.PhysicalParameters.reynolds_A);

            Func<double[], double, double> GravityX_B_Func = (X, t) => -1 / C.PhysicalParameters.reynolds_B * Math.Exp(X[0]) * (Math.Exp(X[0]) * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParameters.reynolds_B - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_B * X[1] * X[1] - Math.Exp(X[0]) * C.PhysicalParameters.reynolds_B - 2 * Math.Sin(X[1]) * C.PhysicalParameters.reynolds_B + 2 * Math.Sin(X[1]));
            Func<double[], double, double> GravityY_B_Func = (X, t) => 2 * Math.Exp(X[0]) * Math.Cos(X[1]) * ((C.PhysicalParameters.reynolds_B - 1) / C.PhysicalParameters.reynolds_B);

            Func<double[], double, double> GravityXX_A_Func = (X, t) => 2 * Math.Exp(2 * X[0]) * C.PhysicalParametersRheology.Weissenberg_a * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParametersRheology.beta_a * X[1] * C.PhysicalParametersRheology.Weissenberg_a
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParametersRheology.beta_a + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] + C.PhysicalParametersRheology.beta_a * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) - X[1] * X[1] + C.PhysicalParametersRheology.beta_a - 1);
            Func<double[], double, double> GravityXY_A_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * (C.PhysicalParametersRheology.beta_a - 1) * C.PhysicalParametersRheology.Weissenberg_a * (2 * Math.Pow(Math.Cos(X[1]), 2) * X[1]
                + 3 * Math.Cos(X[1]) * Math.Sin(X[1]) + X[1]);
            Func<double[], double, double> GravityYY_A_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * C.PhysicalParametersRheology.Weissenberg_a * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParametersRheology.beta_a * X[1]
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParametersRheology.beta_a + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] - 3 * C.PhysicalParametersRheology.beta_a * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1])
                + 3 * X[1] * X[1] - 3 * C.PhysicalParametersRheology.beta_a + 3);

            Func<double[], double, double> GravityXX_B_Func = (X, t) => 2 * Math.Exp(2 * X[0]) * C.PhysicalParametersRheology.Weissenberg_b * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParametersRheology.beta_b * X[1] * C.PhysicalParametersRheology.Weissenberg_b
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParametersRheology.beta_b + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] + C.PhysicalParametersRheology.beta_b * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1]) - X[1] * X[1] + C.PhysicalParametersRheology.beta_b - 1);
            Func<double[], double, double> GravityXY_B_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * (C.PhysicalParametersRheology.beta_b - 1) * C.PhysicalParametersRheology.Weissenberg_b * (2 * Math.Pow(Math.Cos(X[1]), 2) * X[1]
                + 3 * Math.Cos(X[1]) * Math.Sin(X[1]) + X[1]);
            Func<double[], double, double> GravityYY_B_Func = (X, t) => -2 * Math.Exp(2 * X[0]) * C.PhysicalParametersRheology.Weissenberg_b * (-2 * Math.Cos(X[1]) * Math.Sin(X[1]) * C.PhysicalParametersRheology.beta_b * X[1]
                + 3 * Math.Cos(X[1]) * Math.Cos(X[1]) * C.PhysicalParametersRheology.beta_b + 2 * Math.Cos(X[1]) * Math.Sin(X[1]) * X[1] - 3 * C.PhysicalParametersRheology.beta_b * X[1] * X[1] - 3 * Math.Cos(X[1]) * Math.Cos(X[1])
                + 3 * X[1] * X[1] - 3 * C.PhysicalParametersRheology.beta_b + 3);

            #endregion

            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0));

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
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 3;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
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

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 300;
            C.saveperiod = 10;

            #endregion


            return C;
        }
    }
}
