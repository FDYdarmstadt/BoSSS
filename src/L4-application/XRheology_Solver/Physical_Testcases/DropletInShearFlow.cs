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
    /// class providing Controls for the channel flow type testcases
    /// </summary>
    public static class DropletInShearFlow {


        /// <summary>
        /// control object for various testing
        /// </summary>
        public static XRheology_Control DropletInShearFlow_1 (int p = 2, int kelem = 8) {

            XRheology_Control C = new XRheology_Control();

            string _DbPath = null;

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XRheology/Channel";
            C.ProjectDescription = "Channel flow with droplet inside";

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
            C.PhysicalParametersRheology.beta_a = 1.0;
            C.PhysicalParametersRheology.beta_b = 0.0;

            C.RaiseWeissenberg = false;
            C.PhysicalParametersRheology.Weissenberg_a = 0.0;// .3;
            C.PhysicalParametersRheology.Weissenberg_b = 0.0;
            C.WeissenbergIncrement = 0.1;

            double sigma = 0.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.beta_S = 0.05;
            //C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            C.FixedStreamwisePeriodicBC = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 10;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-H, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: C.FixedStreamwisePeriodicBC);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                if (!C.FixedStreamwisePeriodicBC) {
                    grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                    grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] + H) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;

                    if (!C.FixedStreamwisePeriodicBC) {
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
            double r = 0.25;

            Func<double[], double> PhiFunc = (X => ((X[0] - 5).Pow2() + (X[1]).Pow2()).Sqrt() - r);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => X[1]);

            #endregion

            // exact solution
            // ==============
            #region exact

            //Exact Solution Channel
            Func<double[], double, double> VelocityXfunction_A = (X, t) => X[1];   
            Func<double[], double, double> VelocityXfunction_B = (X, t) => X[1];
            Func<double[], double, double> VelocityYfunction_A = (X, t) => 0;
            Func<double[], double, double> VelocityYfunction_B = (X, t) => 0;
            Func<double[], double, double> Pressurefunction_A = (X, t) => 0;
            Func<double[], double, double> Pressurefunction_B = (X, t) => 0;
            Func<double[], double, double> StressXXfunction_A = (X, t) => 0; // WEISSENBERG IS MULTIPLIED IN BC IN FLUX! 
            Func<double[], double, double> StressXXfunction_B = (X, t) => 0;
            Func<double[], double, double> StressXYfunction_A = (X, t) => 0;
            Func<double[], double, double> StressXYfunction_B = (X, t) => 0;
            Func<double[], double, double> StressYYfunction_A = (X, t) => (0.0);
            Func<double[], double, double> StressYYfunction_B = (X, t) => (0.0);

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { VelocityXfunction_A, VelocityYfunction_A });
            C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { VelocityXfunction_B, VelocityYfunction_B });

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", Pressurefunction_A);
            C.ExactSolutionPressure.Add("B", Pressurefunction_B);

            C.ExactSolutionStressXX = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionStressXX.Add("A", (X, t) => StressXXfunction_A(X,t) * C.PhysicalParametersRheology.Weissenberg_a);
            C.ExactSolutionStressXX.Add("B", (X, t) => StressXXfunction_B(X, t) * C.PhysicalParametersRheology.Weissenberg_b);

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

            //C.AddBoundaryValue("wall_upper", "VelocityX#A", X => 1);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", X => 1);

            //C.AddBoundaryValue("wall_lower", "VelocityX#A", X => -1);
            C.AddBoundaryValue("wall_lower", "VelocityX#B", X => -1);

            if (!C.FixedStreamwisePeriodicBC) {
                C.AddBoundaryValue("Velocity_inlet_left", "VelocityX#A", VelocityXfunction_A);
                C.AddBoundaryValue("Velocity_inlet_left", "VelocityX#B", VelocityXfunction_B);

                C.AddBoundaryValue("Velocity_inlet_left", "VelocityY#A", VelocityYfunction_A);
                C.AddBoundaryValue("Velocity_inlet_left", "VelocityY#B", VelocityYfunction_B);

                C.AddBoundaryValue("Velocity_inlet_left", "StressXX#A", StressXXfunction_A);
                C.AddBoundaryValue("Velocity_inlet_left", "StressXX#B", StressXXfunction_B);

                C.AddBoundaryValue("Velocity_inlet_left", "StressXY#A", StressXYfunction_A);
                C.AddBoundaryValue("Velocity_inlet_left", "StressXY#B", StressXYfunction_B);

                C.AddBoundaryValue("Velocity_inlet_left", "StressYY#A", StressYYfunction_A);
                C.AddBoundaryValue("Velocity_inlet_left", "StressYY#B", StressYYfunction_B);

                C.AddBoundaryValue("pressure_outlet_right");
            }


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;     
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.Penalty2 = 1;
            C.AdvancedDiscretizationOptions.Penalty1[0] = 0;
            C.AdvancedDiscretizationOptions.Penalty1[1] = 0;
            //C.AdvancedDiscretizationOptions.PresPenalty2 = 0.0;
            //C.AdvancedDiscretizationOptions.UseWeightedAverages = true;

            C.OperatorMatrixAnalysis = false;
            C.SkipSolveAndEvaluateResidual = false;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XRheology_Control.RefinementStrategy.CurvatureRefined;
            C.RefinementLevel = 3;
            C.BaseRefinementLevel = 3;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 2000;
            C.saveperiod = 1;

            #endregion


            return C;
        }
    }
}
