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
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the Capillary rise testcases
    /// </summary>
    public static class CapillaryRise {


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CapillaryRise_Tube(int p = 2, int kelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            _DbPath = @"\\fdyprime\userspace\smuda\cluster\CapillaryRise\CapillaryRise_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryRise";
            C.ProjectDescription = "Capillary rise in tube";

            C.ContinueOnIoError = false;

            C.LogValues = XNSE_Control.LoggingValues.CapillaryHeight;
            C.LogPeriod = 10;

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

            C.PhysicalParameters.rho_A = 9.97e-4;    
            C.PhysicalParameters.rho_B = 1.204e-6;     
            C.PhysicalParameters.mu_A = 1e-5;        
            C.PhysicalParameters.mu_B = 17.1e-8;     
            C.PhysicalParameters.Sigma = 72.75e-3;

            C.PhysicalParameters.betaS_A = 1e-3;
            C.PhysicalParameters.betaS_B = 1e-5;

            C.PhysicalParameters.betaL = 1e-4;
            C.PhysicalParameters.theta_e = Math.PI / 3.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 0.75e-1;
            double H = 10 * R;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-R, R, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 5 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "pressure_outlet_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                grd.EdgeTagNames.Add(4, "navierslip_linear_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + R) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - R) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double h0 = 1.0e-1;

            Func<double[], double> PhiFunc = (X => X[1] - h0);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("GravityY", X => -9.81e2);

            C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("ae0490a8-7ae3-474d-978e-ccdd1ed8726a"); 
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("pressure_outlet_lower");
            C.AddBoundaryCondition("pressure_outlet_upper");
            C.AddBoundaryCondition("navierslip_linear_left");
            C.AddBoundaryCondition("navierslip_linear_right");


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 50000;
            C.saveperiod = 10;

            #endregion


            return C;
        }


    }
}
