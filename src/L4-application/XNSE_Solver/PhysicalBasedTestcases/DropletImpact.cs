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
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Application.XNSE_Solver.Logging;
using System.Configuration;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for droplet impact testcases
    /// </summary>
    public static class DropletImpact {


        public static XNSE_Control DropletImpactTest_hydrophilicSurface(int k = 3, int numCells = 8) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/DropletImpact";
            C.ProjectDescription = "droplet impact on solid wall (hydrophilic)";

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(k);


            // physical parameters
            // ===================
            #region physics

            C.Tags.Add("Water/Air at 1atm");
            C.PhysicalParameters.rho_A = 1e3;
            C.PhysicalParameters.rho_B = 1.0;
            double mu_scl = 1.0;
            C.PhysicalParameters.mu_A = 1e-3 * mu_scl;
            C.PhysicalParameters.mu_B = 1.48e-5 * mu_scl;
            double sigma = 0.07;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = (60.0 / 180.0) * Math.PI;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double domainSize = 3.0e-3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "slipsymmetry_left");
                grd.EdgeTagNames.Add(2, "pressure_outlet_right");
                grd.EdgeTagNames.Add(3, "navierslip_linear_surface");
                grd.EdgeTagNames.Add(4, "pressure_outlet_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - domainSize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - domainSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // initial values
            // ==============
            #region init

            double radius = 1.0e-3;
            double offset = 0.0; 

            Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()).Sqrt() - radius);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / radius);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            double U0 = 0.2;
            C.InitialValues_Evaluators.Add("VelocityY#A", X => -U0);

            double g = 9.81;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("slipsymmetry_left");
            C.AddBoundaryValue("pressure_outlet_right");
            C.AddBoundaryValue("navierslip_linear_surface");
            C.AddBoundaryValue("pressure_outlet_upper");


            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;


            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }


        public static XNSE_Control DropletImpactTest_hydrophobicSurface(int k = 3, int numCells = 8) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/DropletImpact";
            C.ProjectDescription = "droplet impact on solid wall (hydrophobic)";

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(k);


            // physical parameters
            // ===================
            #region physics

            C.Tags.Add("Water/Air at 1atm");
            C.PhysicalParameters.rho_A = 1e3;
            C.PhysicalParameters.rho_B = 1.0;
            double mu_scl = 1.0;
            C.PhysicalParameters.mu_A = 1e-3 * mu_scl;
            C.PhysicalParameters.mu_B = 1.48e-5 * mu_scl;
            double sigma = 0.07;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = (150.0 / 180.0) * Math.PI;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double domainSize = 3.15e-3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "slipsymmetry_left");
                grd.EdgeTagNames.Add(2, "pressure_outlet_right");
                grd.EdgeTagNames.Add(3, "navierslip_linear_surface");
                grd.EdgeTagNames.Add(4, "pressure_outlet_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - domainSize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - domainSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // initial values
            // ==============
            #region init

            double radius = 1.05e-3;
            double offset = 0.0;

            Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()).Sqrt() - radius);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / radius);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            double U0 = 0.5;
            C.InitialValues_Evaluators.Add("VelocityY#A", X => -U0);

            double g = 9.81;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("slipsymmetry_left");
            C.AddBoundaryValue("pressure_outlet_right");
            C.AddBoundaryValue("navierslip_linear_surface");
            C.AddBoundaryValue("pressure_outlet_upper");


            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;


            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 2e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }

    }
}
