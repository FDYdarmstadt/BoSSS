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
    /// class providing Controls for unphysical tests
    /// </summary>
    public static class GenericTestcases {


        /// <summary>
        /// control object for various testing
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

            C.PhysicalParameters.reynolds_A = 1.0;
            C.PhysicalParameters.reynolds_B = 1.0;
            C.PhysicalParameters.beta_a = 0.0;
            C.PhysicalParameters.beta_b = 0.0;

            C.RaiseWeissenberg = false;
            C.PhysicalParameters.Weissenberg_a = 1.0;
            C.PhysicalParameters.Weissenberg_b = 1.0;
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
            C.LinearSolver.SolverCode = LinearSolverCode.classic_mumps;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.LinearSolver.MaxSolverIterations = 10;
            C.LinearSolver.ConvergenceCriterion = 1e-8;

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

            C.Timestepper_Scheme = XRheology_Control.TimesteppingScheme.ImplicitEuler;
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
