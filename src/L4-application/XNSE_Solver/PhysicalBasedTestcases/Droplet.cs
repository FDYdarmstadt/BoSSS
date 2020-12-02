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

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the droplet testcases
    /// </summary>
    public static class Droplet {

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control StaticDroplet_Free(int p = 2, int kelem = 16, int AMRlvl = 0) {

            XNSE_Control C = new XNSE_Control();

            int D = 3;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;
            bool steadyInterface = true;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_test_db";
            string _DbPath = null; // @"\\HPCCLUSTER\hpccluster-scratch\smuda\XNSE_studyDB";
            //string _DbPath = @"\\terminal03\Users\smuda\local\terminal03_XNSE_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "StaticDroplet";
            //C.ProjectDescription = "Static droplet";
            //C.SessionName = "SD_meshStudy_Hysing_mesh" + kelem; // "_AMR"+AMRlvl;

            C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.Dropletlike;
            C.PostprocessingModules.Add(new Dropletlike() { LogPeriod = 10 });

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
            if(D == 3) {
                C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }
            //C.FieldOptions.Add("GravityY", new FieldOpts() {
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = Math.Max(2, p),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = Math.Max(2, p),
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

            //C.Tags.Add("Hysing");
            C.Tags.Add("La = 5000");
            C.PhysicalParameters.rho_A = 1e4;
            C.PhysicalParameters.rho_B = 1e4;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            double sigma = 1.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.Tags.Add("La = 0.005");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 10;
            //double sigma = 1;
            //C.PhysicalParameters.Sigma = sigma;

            //Air - Water(lenght scale == centimeters, 3D space)
            //C.PhysicalParameters.rho_A = 1e3;      // kg / cm^3
            //C.PhysicalParameters.rho_B = 1.2;    // kg / cm^3
            //C.PhysicalParameters.mu_A = 1e-3;       // kg / cm * sec
            //C.PhysicalParameters.mu_B = 17.1e-6;    // kg / cm * sec
            //double sigma = 72.75e-3;                // kg / sec^2 
            //C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double Lscl = 1.0;
            double xSize = Lscl * 1.0;
            double ySize = Lscl * 1.0;
            double zSize = Lscl * 1.0;

            if(D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize/2.0, xSize/2.0, kelem + 0);
                    double[] Ynodes = GenericBlas.Linspace(-ySize/2.0, ySize/2.0, kelem + 0);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                    //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "wall_lower");
                    grd.EdgeTagNames.Add(2, "wall_upper");
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[1] + ySize/2.0) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[1] - ySize/2.0) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0] + xSize/2.0) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize/2.0) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            if(D == 3) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize / 2.0, xSize / 2.0, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-ySize / 2.0, ySize / 2.0, kelem + 1);
                    double[] Znodes = GenericBlas.Linspace(-zSize / 2.0, zSize / 2.0, kelem + 1);
                    var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

                    grd.EdgeTagNames.Add(1, "wall_lower");
                    grd.EdgeTagNames.Add(2, "wall_upper");
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");
                    grd.EdgeTagNames.Add(5, "wall_front");
                    grd.EdgeTagNames.Add(6, "wall_back");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1] + ySize / 2.0) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - ySize / 2.0) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0] + xSize / 2.0) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - xSize / 2.0) <= 1.0e-8)
                            et = 4;
                        if (Math.Abs(X[2] + zSize / 2.0) <= 1.0e-8)
                            et = 5;
                        if (Math.Abs(X[2] - zSize / 2.0) <= 1.0e-8)
                            et = 6;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double r = Lscl * 0.25;

            Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - r);         // signed distance
            //Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()) - r.Pow2());         // quadratic

            if(D == 3) {
                double a = 1.25;
                double b = 1.25;
                double c = 0.64;
                PhiFunc = (X => (((X[0] - 0.0).Pow2() / a.Pow2()) + ((X[1] - 0.0).Pow2() / b.Pow2() ) + ((X[2] - 0.0).Pow2() / c.Pow2())).Sqrt() - r);
            }       

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            double Pjump = sigma / r;
            if(D == 3)
                Pjump *= 2.0;

            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => 0.0);


            // restart
            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("c95b2833-288e-4014-ba08-affa65a2398e");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // exact solution
            // ==============
            #region exact

            C.Phi = ((X, t) => PhiFunc(X));

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            if(D == 2) {
                C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
                C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            }

            if(D == 3) {
                C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0, (X, t) => 0.0 });
                C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0, (X, t) => 0.0 });
            }

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");
            if(D == 3) {
                C.AddBoundaryValue("wall_front");
                C.AddBoundaryValue("wall_back");
            }


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;
            C.solveKineticEnergyEquation = false;

            C.CheckJumpConditions = false;
            C.CheckInterfaceProps = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //if (AMRlvl > 0) {
            //    C.AdaptiveMeshRefinement = true;
            //    C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            //    C.BaseRefinementLevel = AMRlvl;
            //}
            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.BaseRefinementLevel = 1;
            C.RefinementLevel = 1;
            //C.AMR_startUpSweeps = 2;

            //C.InitSignedDistance = false;
            C.adaptiveReInit = false;

            //C.LinearSolver.SolverCode = LinearSolverCode.exp_Kcycle_schwarz;
            C.LinearSolver.SolverCode = LinearSolverCode.automatic;
            //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.LinearSolver.NoOfMultigridLevels = 2;
            C.NonLinearSolver.MaxSolverIterations = 50; 
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-9;
            C.LevelSet_ConvergenceCriterion = 1e-7;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };


            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.Option_LevelSetEvolution = (steadyInterface) ? LevelSetEvolution.None : LevelSetEvolution.Fourier;
            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.None;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = 1.0 * sigma;
            //C.PhysicalParameters.lambda_I = 2.0 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.UseWeightedAverages = false;
            //C.InterAverage = XNSE_Control.InterfaceAveraging.viscosity;


            //C.SessionName = "SD_meshStudy_Hysing_mesh" + kelem + "_Fourier";


            if (C.Option_LevelSetEvolution == LevelSetEvolution.Fourier) {

                int numSp = 1800;
                double[] FourierP = new double[numSp];
                double[] samplP = new double[numSp];
                for (int sp = 0; sp < numSp; sp++) {
                    FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                    samplP[sp] = r;
                }

                C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
                    center = new double[] { 0.0, 0.0},
                    FourierEvolve = Fourier_Evolution.MaterialPoints,
                    centerMove = CenterMovement.Reconstructed,
                };
            }


            #endregion


            // Timestepping
            // ============
            #region time

            //switch(p) {
            //    case 1: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //            break;
            //        }
            //    case 2: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //            C.Timestepper_BDFinit = TimeStepperInit.MultiInit;
            //            break;
            //        }
            //    default:
            //        C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //        C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //        break;
            //}

            //if(D == 3) {
            //    C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //    C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //}

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;


            //C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            //C.LSunderrelax = 0.05;
            C.Timestepper_LevelSetHandling = (steadyInterface) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = 0.1; //0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1; // 12500; // (int)(125.0 / dt);
            C.saveperiod = 10;
            
            #endregion


            return C;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control OscillatingDroplet(int p = 2, int kelem = 21, int method = 0, double mu_scl = 1.0, bool onlyB = false, double tScale = 1.0) {

            XNSE_Control C = new XNSE_Control();

            bool hysing = true;
            bool kummer = false;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"\\dc1\userspace\yotov\bosss-db";
            //string _DbPath = @"\\HPCCLUSTER\hpccluster-scratch\smuda\XNSE_studyDB";
            //string _DbPath = @"\\terminal03\Users\smuda\local\terminal03_XNSE_studyDB";
            //string _DbPath = @"D:\local\local_Testcase_databases\Testcase_OscillatingDroplet";
            string _DbPath = null;

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "OscillatingDroplet";
            if (hysing) {
                C.SessionName = "OD_meshStudy_Hysing2_mesh" + kelem + "_rerun";
            } else {
                C.SessionName = "OD_AirWater_";
            }

            C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.Dropletlike;
            //C.LogPeriod = 10;
            C.PostprocessingModules.Add(new Dropletlike() { LogPeriod = 4 });

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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
                Degree = Math.Max(p,2),
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

            if (hysing) {
                double rho = 1e4;
                double mu = 1.0;
                double sigma = 0.1;

                C.Tags.Add("Hysing");
                C.Tags.Add("La = 500");
                C.PhysicalParameters.rho_A = rho;
                C.PhysicalParameters.rho_B = rho;
                C.PhysicalParameters.mu_A = mu;
                C.PhysicalParameters.mu_B = mu;
                C.PhysicalParameters.Sigma = sigma;

            }

            if (kummer) {

                C.PhysicalParameters.rho_A = 1;
                C.PhysicalParameters.rho_B = 1;
                C.PhysicalParameters.mu_A = 0.5;
                C.PhysicalParameters.mu_B = 0.05;
                C.PhysicalParameters.Sigma = 0.1;
            }

            //C.Tags.Add("La = 0.005");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 10;
            //double sigma = 1;
            //C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.beta_S = 0.05;
            //C.PhysicalParameters.Theta_e = Math.PI / 3.0;

            if (!hysing) {
                // Air - Water: 
                C.PhysicalParameters.rho_A = 1e3;
                C.PhysicalParameters.rho_B = 1.2;
                C.PhysicalParameters.mu_A = (onlyB) ? 1e-3 : 1e-3 * mu_scl;
                C.PhysicalParameters.mu_B = 17.1e-6 * mu_scl;
                C.PhysicalParameters.Sigma = 72.75e-3;
            }


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // grid generation
            // ===============
            #region grid

            double Lscale = hysing ? 1.0 : 0.01;
            double L = 1.0 * Lscale;
            //double xSize = 1.0;
            //double ySize = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-(L / 2.0), (L / 2.0), kelem + 0);
                double[] Ynodes = GenericBlas.Linspace(-(L / 2.0), (L / 2.0), kelem + 0);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double r = 0.25 * Lscale;

            //Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - r);         // signed distance
            //Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()) - r.Pow2());         // quadratic

            double a = r;
            double b = r;
            if (hysing) {
                a = 1.25 * r;
                b = 0.8 * r;
            } else {
                a = 1.1 * r;
                b = 0.91 * r;
            }
            Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() / a.Pow2() + (X[1] - 0.0).Pow2() / b.Pow2()) - 1);          // ellipse                     

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            //double Pjump = sigma / r;
            //C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            //C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => 0.0);


            //// restart
            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("c95b2833-288e-4014-ba08-affa65a2398e");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // exact solution
            // ==============
            #region exact

            //C.Phi = ((X, t) => PhiFunc(X));

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.solveKineticEnergyEquation = false;

            //C.ComputeEnergyProperties = true;
            //C.FieldOptions.Add("KineticEnergy", new FieldOpts() {
            //    Degree = 2 * p,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});

            C.CheckJumpConditions = false;
            C.CheckInterfaceProps = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.EnforceLevelSetConservation = true;

            //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //if (!hysing) {
            //    C.AdaptiveMeshRefinement = true;
            //    C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            //    C.BaseRefinementLevel = 3;
            //    C.AMR_startUpSweeps = 4;
            //    C.RefinementLevel = 1;
            //}

            //C.ReInitPeriod = 100;
            //C.useFiltLevSetGradientForEvolution = true;

            //C.AdvancedDiscretizationOptions.LFFA = 0.9;
            //C.AdvancedDiscretizationOptions.LFFB = 0.9;

            //C.AdaptiveMeshRefinement = true;
            //C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            //C.BaseRefinementLevel = 2;
            //C.AMR_startUpSweeps = 2;
            //C.SessionName = C.SessionName + "_AMR1";

            int numSp = 1440;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                double angle = sp * (2 * Math.PI / (double)numSp);
                FourierP[sp] = angle;
                samplP[sp] = a * b / Math.Sqrt((a * Math.Cos(angle + Math.PI / 2)).Pow2() + (b * Math.Sin(angle + Math.PI / 2)).Pow2());
            }

            FourierLevSetControl FourierCntrl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (Math.Pow(2,C.BaseRefinementLevel)*(double)kelem)) {
                center = new double[] { 0.0, 0.0 },
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                centerMove = CenterMovement.Reconstructed,
            };


            C.SetLevelSetMethod(method, FourierCntrl);
            //C.SessionName = "OscillatingDroplet_setup3_muScl"+mu_scl+"_methodStudy_k2_" + C.methodTagLS;
            //C.Option_LevelSetEvolution = LevelSetEvolution.None;
            //C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.InitSignedDistance = true;


            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = 1.0;
            //C.PhysicalParameters.lambda_I = 2.0;
            //C.SessionName = C.SessionName + "_curvature";

            C.AdvancedDiscretizationOptions.SetSurfaceTensionMaxValue = false;
            //C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.surfaceDeformationRateLocal;

            #endregion


            // Timestepping
            // ============
            #region time

            //switch (p) {
            //    case 1: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //            break;
            //        }
            //    case 2: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //            C.Timestepper_BDFinit = TimeStepperInit.MultiInit;
            //            break;
            //        }
            //    default:
            //        C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //        C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //        break;

            //}
            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;
            //C.LSunderrelax = 0.5;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = (hysing) ? 0.5 : 1.5e-5 * tScale;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = (hysing) ? 1000 : 0.15;
            C.NoOfTimesteps = (hysing) ? 2000 : (int)(0.15 / dt);
            C.saveperiod = 10;

            #endregion

            return C;
        }


        public static XNSE_Control OscillatingDroplet_Kummer(int p = 2, int kelem = 27, int method = 0, double mu_scl = 1.0, bool onlyB = false) {

            XNSE_Control C = new XNSE_Control();

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            string _DbPath = null; // @"D:\local\local_Testcase_databases\Testcase_OscillatingDroplet";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "OscillatingDroplet";
            C.SessionName = "OD_AirWater_Kummer";
            
            C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.Dropletlike;
            //C.LogPeriod = 10;
            C.PostprocessingModules.Add(new Dropletlike());
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

            // Air - Water: 
            C.PhysicalParameters.rho_A = 1e3;
            C.PhysicalParameters.rho_B = 1.2;
            C.PhysicalParameters.mu_A = (onlyB) ? 1e-3 : 1e-3 * mu_scl;
            C.PhysicalParameters.mu_B = 17.1e-6 * mu_scl;
            C.PhysicalParameters.Sigma = 72.75e-3;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // grid generation
            // ===============
            #region grid

            double L = 9e-2;
            //double xSize = 1.0;
            //double ySize = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-(L / 2.0), (L / 2.0), kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-(L / 2.0), (L / 2.0), kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - (L / 2.0)) <= 1.0e-8)
                        et = 1;
                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init


            double a = 0.9185e-2;
            double b = 0.7515e-2;
            Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() / a.Pow2() + (X[1] - 0.0).Pow2() / b.Pow2()) - 1);          // ellipse                     

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            //double Pjump = sigma / r;
            //C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            //C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => 0.0);


            //// restart
            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("c95b2833-288e-4014-ba08-affa65a2398e");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // exact solution
            // ==============
            #region exact

            //C.Phi = ((X, t) => PhiFunc(X));

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.solveKineticEnergyEquation = false;

            C.ComputeEnergyProperties = true;
            C.FieldOptions.Add("KineticEnergy", new FieldOpts() {
                Degree = 2 * p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            C.CheckJumpConditions = false;
            C.CheckInterfaceProps = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.EnforceLevelSetConservation = true;

            //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.BaseRefinementLevel = 1;
            //C.RefinementLevel = 4;

            //C.ReInitPeriod = 100;
            //C.useFiltLevSetGradientForEvolution = true;

            //C.AdvancedDiscretizationOptions.LFFA = 0.9;
            //C.AdvancedDiscretizationOptions.LFFB = 0.9;

            //C.AdaptiveMeshRefinement = true;
            //C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            //C.BaseRefinementLevel = 1;
            //C.SessionName = C.SessionName + "_AMR1";

            int numSp = 1440;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                double angle = sp * (2 * Math.PI / (double)numSp);
                FourierP[sp] = angle;
                samplP[sp] = a * b / Math.Sqrt((a * Math.Cos(angle + Math.PI / 2)).Pow2() + (b * Math.Sin(angle + Math.PI / 2)).Pow2());
            }

            FourierLevSetControl FourierCntrl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (Math.Pow(2, C.BaseRefinementLevel) * (double)kelem)) {
                center = new double[] { 0.0, 0.0 },
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                centerMove = CenterMovement.Reconstructed,
            };


            C.SetLevelSetMethod(method, FourierCntrl);
            //C.SessionName = "OscillatingDroplet_setup3_muScl"+mu_scl+"_methodStudy_k2_" + C.methodTagLS;


            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = 1.0;
            //C.PhysicalParameters.lambda_I = 2.0;
            //C.SessionName = C.SessionName + "_curvature";

            C.AdvancedDiscretizationOptions.SetSurfaceTensionMaxValue = false;

            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;
            //C.LSunderrelax = 0.5;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = 0.0005;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1.5;
            C.NoOfTimesteps = (int)(1.5 / dt);
            C.saveperiod = 10;
           
            #endregion

            return C;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control TheKartoffel(int p = 2, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            _DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_OscillatingDroplet";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Static droplet";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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
            //C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
            //    Degree = p - 1,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("EnergyJump", new FieldOpts() {
            //    Degree = 2 * p - 1,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("SurfEnergy_Changerate", new FieldOpts() {
            //    Degree = 2 * p - 1,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("CurvEnergy", new FieldOpts() {
            //    Degree = 2 * p - 1,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Hysing");
            C.Tags.Add("La = 500");
            C.PhysicalParameters.rho_A = 1e4;
            C.PhysicalParameters.rho_B = 1e4;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            double sigma = 0.1;
            C.PhysicalParameters.Sigma = sigma;

            //C.Tags.Add("La = 0.005");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 1;

            // Air-Water (lenght scale == centimeters, 3D space)
            //C.PhysicalParameters.rho_A = 1e-3;      // kg / cm^3
            //C.PhysicalParameters.rho_B = 1.2e-6;    // kg / cm^3
            //C.PhysicalParameters.mu_A = 1e-5;       // kg / cm * sec
            //C.PhysicalParameters.mu_B = 17.1e-8;    // kg / cm * sec
            //double sigma = 72.75e-3;                // kg / sec^2 
            //C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            double sizeFactor = 1.0 / 3.0;
            double xSize = sizeFactor * 4.5;
            double ySize = sizeFactor * 4.5;

            int xkelem = 40 * 1;
            int ykelem = 40 * 1;

            double hMin = Math.Min(2 * xSize / (xkelem), 2 * ySize / (ykelem));


            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, ykelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + ySize) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;
                    return et;
                });

                return grd;
            };


            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => (X[0].Pow2() / (0.8 * 0.8) * 1.25 + X[1].Pow2() / (0.8 * 0.8) * 0.8) - 1.0 + 0.2 * Math.Sin(10 * X[0] * X[1]))); // Kartoffel      

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //// restart
            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("73d5810e-4076-4dfd-8615-dc7e3a60bdc8");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.solveKineticEnergyEquation = false;
            C.ComputeEnergyProperties = false;

            C.CheckJumpConditions = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.AdaptiveMeshRefinement = true;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 80;
            C.LinearSolver.MaxSolverIterations = 80;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };

            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.ExtensionVelocity;
            //C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.FullBoussinesqScriven;
            C.PhysicalParameters.mu_I = 1 * sigma;
            C.PhysicalParameters.lambda_I = 2 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            //if (C.Option_LevelSetEvolution == LevelSetEvolution.Fourier) {

            //    int numSp = 640;
            //    double[] FourierP = new double[numSp];
            //    double[] samplP = new double[numSp];
            //    for (int sp = 0; sp < numSp; sp++) {
            //        FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
            //        samplP[sp] = r;
            //    }

            //    C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
            //        center = new double[] { 0.5, 0.5 },
            //        FourierEvolve = Fourier_Evolution.MaterialPoints,
            //        centerMove = CenterMovement.Reconstructed,
            //    };
            //}


            #endregion


            // Timestepping
            // ============
            #region time

            //switch (p) {
            //    case 1: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //            break;
            //        }
            //    case 2: {
            //            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //            C.Timestepper_BDFinit = TimeStepperInit.MultiInit;
            //            break;
            //        }
            //    default:
            //        C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //        C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //        break;

            //}
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;


            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;
            //C.LSunderrelax = 0.5;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10000; // (int)(125.0 / dt);
            C.saveperiod = 10;

            #endregion


            return C;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CollidingDroplets(int p = 1, int kelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_CollidingDroplets";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplets";
            C.ProjectDescription = "colliding droplets";

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
                Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // grid genration
            // ==============
            #region grid


            double xSize = 3.0;
            double ySize = 3.0;

            C.GridFunc = delegate () {
                //double[] Xnodes = GenericBlas.Linspace(-xSize / 2.0, xSize / 2.0, (int)xSize * kelem + 1);     // + 1 collision at cell boundary
                //double[] Ynodes = GenericBlas.Linspace(-ySize / 2.0, ySize / 2.0, (int)ySize * kelem + 1);

                int scl = kelem / 2;

                var _xNodes1 = Grid1D.TanhSpacing(-1.5, -0.1, 2 * scl + 1, 1.5, false);
                _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
                var _xNodes2 = GenericBlas.Linspace(-0.1, 0.1, scl + 1);
                _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
                var _xNodes3 = Grid1D.TanhSpacing(0.1, 1.5, 2 * scl + 1, 1.5, true);

                var xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3);


                var _yNodes1 = Grid1D.TanhSpacing(-1.5, 0.0, (5 / 2) * scl + 1, 1.5, false);
                _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                //var _yNodes2 = GenericBlas.Linspace(-1.2, 1.2, Convert.ToInt32(40 * MeshFactor)); //40
                //_yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
                var _yNodes3 = Grid1D.TanhSpacing(0.0, 1.5, (5 / 2) * scl + 1, 1.5, true);

                var yNodes = ArrayTools.Cat(_yNodes1, _yNodes3);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false);

                //int scl = kelem;
                //double xSizeC = 0.6;
                //double ySizeC = 0.3;

                //var _xNodesL = Grid1D.TanhSpacing(-1.5, -xSizeC, scl/2, 0.5, false);
                //var _xNodesR = Grid1D.TanhSpacing(xSizeC, 1.5, scl/2, 0.5, true);
                //var _xNodesC = GenericBlas.Linspace(-xSizeC, xSizeC, 2*scl);
                //_xNodesC = _xNodesC.GetSubVector(1, (_xNodesC.Length - 2));
                //var xNodes = ArrayTools.Cat(_xNodesL, _xNodesC, _xNodesR);

                //var _yNodesB = Grid1D.TanhSpacing(-1.5, -ySizeC, scl/2, 0.5, false);
                //var _yNodesU = Grid1D.TanhSpacing(ySizeC, 1.5, scl/2, 0.5, true);
                //var _yNodesC = GenericBlas.Linspace(-ySizeC, ySizeC, scl);
                //_yNodesC = _yNodesC.GetSubVector(1, (_yNodesC.Length - 2));
                //var yNodes = ArrayTools.Cat(_yNodesB, _yNodesC, _yNodesU);

                //var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                grd.EdgeTagNames.Add(1, "pressure_outlet_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + ySize / 2.0) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize / 2.0) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + xSize / 2.0) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xSize / 2.0) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("pressure_outlet_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");


            #endregion


            // Initial Values
            // ==============
            #region init

            // left droplet
            double[] center_l = new double[] { -0.3, 0.0 };
            double radius_l = 0.25;
            Func<double[], double> bubble_l = (X => ((X[0] - center_l[0]).Pow2() + (X[1] - center_l[1]).Pow2()).Sqrt() - radius_l); // signed-distance form

            // right droplet
            double[] center_r = new double[] { 0.3, 0.0 };
            double radius_r = 0.25;
            Func<double[], double> bubble_r = (X => ((X[0] - center_r[0]).Pow2() + (X[1] - center_r[1]).Pow2()).Sqrt() - radius_r); // signed-distance form

            Func<double[], double> PhiFunc = (X => Math.Min(bubble_l(X), bubble_r(X)));

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double vel_collision = 2.0;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => (X[0] < 0.0) ? vel_collision / 2.0 : -vel_collision / 2.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            #endregion


            // Physical Parameters
            // ===================
            #region physics


            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1e-3;
            C.PhysicalParameters.Sigma = 0.0;


            //C.PhysicalParameters.rho_A = 1e-2;
            //C.PhysicalParameters.rho_B = 1.5e-5;
            //C.PhysicalParameters.mu_A = 7e-4;
            //C.PhysicalParameters.mu_B = 6e-6;
            //C.PhysicalParameters.Sigma = 0.05;

            //// tetradecane(A) in nitrogen(B): in m 
            //C.PhysicalParameters.rho_A = 764;
            //C.PhysicalParameters.rho_B = 1.25;
            //C.PhysicalParameters.mu_A = 2e-3;
            //C.PhysicalParameters.mu_B = 16.6e-6;
            //C.PhysicalParameters.Sigma = 26.56e-3;

            //// tetradecane(A) in nitrogen(B): in mm 
            //C.PhysicalParameters.rho_A = 7.64e-7;
            //C.PhysicalParameters.rho_B = 1.25e-9;
            //C.PhysicalParameters.mu_A = 2e-6;
            //C.PhysicalParameters.mu_B = 16.6e-9;
            //C.PhysicalParameters.Sigma = 26.56e-3;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            //C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 25;
            C.LinearSolver.MaxSolverIterations = 25;
            //C.Solver_MaxIterations = 25;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            double dt = 5e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 10;

            #endregion

            return C;

        }


        /// <summary>
        /// control object for the use in BoSSSPad Worksheet mode
        /// additional parameters need to be set during job creation
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control Droplet_forWorksheet(bool steadyInterface) {

            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================
            #region db

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = true;
            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.Dropletlike;
            C.PostprocessingModules.Add(new Dropletlike());
            #endregion


            // DG degrees
            // ==========
            #region degrees

            // need to be set by user via setDGdegree() in worksheet

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

            // need to be set during job creation 

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid genration
            // ==============
            #region grid

            // need to be set by user via setGrid() in worksheet

            #endregion


            // boundary conditions
            // ===================
            #region BC

            // need to be set during job creation 

            #endregion


            // Initial Values
            // ==============
            #region init

            // need to be set during job creation 

            #endregion


            // additional parameters (evaluation)
            // ==================================

            // need to be set during job creation 


            // misc. solver options
            // ====================
            #region solver

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;

            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;

            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion


            // Level-Set options (AMR)
            // =======================
            #region levset

            C.LSContiProjectionMethod = (steadyInterface) ? Solution.LevelSetTools.ContinuityProjectionOption.None : Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.Option_LevelSetEvolution = (steadyInterface) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = (steadyInterface) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            //C.dtMax = dt; // need to be set according to grid and DG degree
            //C.dtMin = dt;
            C.Endtime = 1000;
            //C.NoOfTimesteps = 0; 

            C.saveperiod = 1;
            
            #endregion


            return C;

        }


        // ==========================================
        // control-objects for elementalTestProgramm
        // ==========================================

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control SD_SurfaceTensionTest(int p = 2, int kelem = 20, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_test_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Static droplet";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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

            //C.Tags.Add("Hysing");
            //C.Tags.Add("La = 5000");
            //C.PhysicalParameters.rho_A = 1e4;
            //C.PhysicalParameters.rho_B = 1e4;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 1;
            //double sigma = 1;
            //C.PhysicalParameters.Sigma = sigma;

            //C.Tags.Add("La = 0.005");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 1;

            // Air-Water (lenght scale == centimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-3;      // kg / cm^3
            C.PhysicalParameters.rho_B = 1.2e-6;    // kg / cm^3
            C.PhysicalParameters.mu_A = 1e-5;       // kg / cm * sec
            C.PhysicalParameters.mu_B = 17.1e-8;    // kg / cm * sec
            double sigma = 72.75e-3;                // kg / sec^2 
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double r = 0.25;

            Func<double[], double> PhiFunc = (X => ((X[0] - 0.5).Pow2() + (X[1] - 0.5).Pow2()).Sqrt() - r);         // signed distance
            ////Func<double[], double> PhiFunc = (X => ((X[0] - 0.5).Pow2() + (X[1] - 0.5).Pow2()) - r.Pow2());         // quadratic
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            double Pjump = sigma / r;
            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => 0.0);


            //// restart
            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("73d5810e-4076-4dfd-8615-dc7e3a60bdc8");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // exact solution
            // ==============
            #region exact

            C.Phi = ((X, t) => PhiFunc(X));

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.solveKineticEnergyEquation = true;
            C.ComputeEnergyProperties = false;

            C.CheckJumpConditions = true;
            C.CheckInterfaceProps = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.AdaptiveMeshRefinement = true;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 80;
            C.LinearSolver.MaxSolverIterations = 80;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-9;
            C.LevelSet_ConvergenceCriterion = 1e-7;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            //C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };

            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.PhysicalParameters.mu_I = 10 * sigma;
            C.PhysicalParameters.lambda_I = 20 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            if (C.Option_LevelSetEvolution == LevelSetEvolution.Fourier) {

                int numSp = 640;
                double[] FourierP = new double[numSp];
                double[] samplP = new double[numSp];
                for (int sp = 0; sp < numSp; sp++) {
                    FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                    samplP[sp] = r;
                }

                C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
                    center = new double[] { 0.5, 0.5 },
                    FourierEvolve = Fourier_Evolution.MaterialPoints,
                    centerMove = CenterMovement.Reconstructed,
                };
            }


            #endregion


            // Timestepping
            // ============
            #region time

            switch (p) {
                case 1: {
                        C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
                        C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
                        break;
                    }
                case 2: {
                        C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
                        C.Timestepper_BDFinit = TimeStepperInit.MultiInit;
                        break;
                    }
                default:
                    C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
                    C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
                    break;

            }
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;


            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Iterative;
            //C.LSunderrelax = 0.05;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = 12.5e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 125; // (int)(125.0 / dt);
            C.saveperiod = 1;

            #endregion


            return C;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control OscillatingDroplet_AMRtest(int p = 2, int kelem = 20) {

            XNSE_Control C = new XNSE_Control();

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "AMR Test";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 0);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 0);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // prescribed level-set movement
            // =============================

            double r = 0.25;

            double a = 0.6;
            double b = -0.6;
            double period = 1.0;
            Func<double[], double, double> PhiFunc = ((X, t) => ((X[0] - (xSize / 2.0)).Pow2() / (r * (1.0 + a * Math.Sin(2.0 * Math.PI * t / period))).Pow2()
                                                                    + (X[1] - (xSize / 2.0)).Pow2() / (r * (1.0 + b * Math.Sin(2.0 * Math.PI * t / period))).Pow2()) - 1);          // ellipse     


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0.0));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            #endregion


            // exact solution
            // ==============
            #region exact

            C.Phi = (PhiFunc);

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("Wall_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.solveKineticEnergyEquation = false;
            C.ComputeEnergyProperties = false;

            C.CheckJumpConditions = false;
            C.CheckInterfaceProps = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.EnforceLevelSetConservation = true;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 80;
            C.LinearSolver.MaxSolverIterations = 80;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.Prescribed;
            //C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.CurvatureRefined;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.LieSplitting;

            C.TimesteppingMode = compMode;

            double dt = period / 100.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 5 * period;
            C.NoOfTimesteps = 500;
            C.saveperiod = 1;

            #endregion


            return C;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control DropletOnPlate_AMRtest(int p = 2, int kelem = 16) {

            XNSE_Control C = new XNSE_Control();

            int D = 2;

            if(D == 3)
                C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "AMR test";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;

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
            if(D == 3) {
                C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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

            //C.PhysicalParameters.betaS_A = 0.05;
            //C.PhysicalParameters.betaS_B = 0.05;

            //C.PhysicalParameters.betaL = 0;
            //C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            double xSize = 0.3;
            double ySize = 0.3;
            double zSize = 0.3;

            if(D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 0);
                    double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 0);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            if(D == 3) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 1);
                    double[] Znodes = GenericBlas.Linspace(0, zSize, kelem + 1);
                    var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    grd.EdgeTagNames.Add(5, "navierslip_linear_front");
                    grd.EdgeTagNames.Add(6, "navierslip_linear_back");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[2]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[2] - zSize) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;
                        if(Math.Abs(X[1]) <= 1.0e-8)
                            et = 5;
                        if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                            et = 6;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.1;
            //double Theta_e = Math.PI / 2.0;
            //double s = 2 * R * Math.Sin(Theta_e);
            //double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            Func<double[], double, double> PhiFunc = ((X, t) => -1.0);
            double prescribedVel = 0.0;

            if(D == 2) {
                double[] center_0 = new double[] { xSize / 2.0, ySize / 2.0 };
                prescribedVel = ySize / 2.1;
                PhiFunc = ((X, t) => ((X[0] - (center_0[0] + (0.0 * prescribedVel))).Pow2() + (X[1] - (center_0[1] - (t * prescribedVel))).Pow2()).Sqrt() - R);
            }

            if(D == 3) {
                double[] center_0 = new double[] { xSize / 2.0, ySize / 2.0, zSize / 2.0 };
                prescribedVel = zSize / 2.0;
                PhiFunc = ((X, t) => ((X[0] - center_0[0]).Pow2() + (X[1] - center_0[1]).Pow2() + (X[2] - (center_0[2] - (t * prescribedVel))).Pow2()).Sqrt() - R);
            }

            C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0.0));

            C.InitialValues_Evaluators.Add("Pressure#A", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("navierslip_linear_lower");
            C.AddBoundaryValue("navierslip_linear_upper");
            C.AddBoundaryValue("navierslip_linear_left");
            C.AddBoundaryValue("navierslip_linear_right");

            if(D == 3) {
                C.AddBoundaryValue("navierslip_linear_front");
                C.AddBoundaryValue("navierslip_linear_back");
            }

            #endregion


            // exact solution
            // ==============

            C.Phi = PhiFunc;

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);



            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.Prescribed;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.ContactLineRefined;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.LieSplitting;

            C.TimesteppingMode = compMode;

            double dt = prescribedVel / 50.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 5;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }


    }

}

