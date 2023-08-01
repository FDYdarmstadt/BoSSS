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
using static BoSSS.Solution.AMRLevelIndicatorLibrary;
using BoSSS.Solution;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for droplet impact testcases
    /// </summary>
    public static class DropletImpact {


        public static XNSE_Control MovingWallBoundaryLayer(int k = 3, int numCells = 4) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/MovingWallBundaryLayer";
            C.ProjectDescription = "test setup for resolving the boundary layer at a moving wall";

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(k);


            // physical parameters
            // ===================
            #region physics

            //C.Tags.Add("Water at 1atm");
            //C.PhysicalParameters.rho_A = 1e3;
            //C.PhysicalParameters.mu_A = 1e-3;

            C.Tags.Add("Air at 1atm");
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 1.48e-5;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double domainSize = 3.0e-3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 10 * domainSize, 10 * numCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                //grd.EdgeTagNames.Add(2, "freestream_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    //if (Math.Abs(X[0]) <= 1.0e-8)
                    //    et = 1;
                    //if (Math.Abs(X[0] - domainSize) <= 1.0e-8)
                    //    et = 2;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - 10 * domainSize) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            #endregion


            // initial values
            // ==============
            #region init

            Func<double[], double> PhiFunc = X => -1.0;

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //double g = 9.81;
            //C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double Uwall = -0.1;
            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => Uwall);
            C.AddBoundaryValue("velocity_inlet_upper");
            //C.AddBoundaryValue("freestream_upper");

            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;

            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonBoundary(new byte[] { 1 }) { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;

            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            #endregion


            return C;

        }


        static public Grid3D RotatingDiskSector_Linearized(double radiusOP, double l_azimuthal, double l_radial, double h_axial, int res_azimuthal, int res_radial, int res_axial) {

            double[] xNodes = GenericBlas.Linspace(-l_azimuthal / 2.0, l_azimuthal / 2.0, res_azimuthal + 1);    // azimuthal direction
            double[] yNodes = GenericBlas.Linspace(-l_radial / 2.0, l_radial / 2.0, res_radial + 1);    // radial direction
            double[] zNodes = GenericBlas.Linspace(0.0, h_axial, res_axial + 1);    // axial direction

            var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes, periodicX: true);
            grd.Name = $"RotatingDiskSector3D_Linearized_{res_azimuthal}x{res_radial}x{res_axial}";

            grd.EdgeTagNames.Add(1, "wall_rotatingDisk");
            grd.EdgeTagNames.Add(2, "freestream_top");
            grd.EdgeTagNames.Add(3, "pressure_outlet_front");
            grd.EdgeTagNames.Add(4, "pressure_outlet_back");
            // grd.EdgeTagNames.Add(5, "pressure_outlet_upstream");
            // grd.EdgeTagNames.Add(6, "pressure_outlet_downstream");

            grd.DefineEdgeTags(delegate (Vector X) {
                byte et = 0;
                if (X.z.Abs() <= 1e-8)
                    et = 1;
                if ((X.z - h_axial).Abs() <= 1e-8)
                    et = 2;
                if ((X.y + (l_radial / 2.0)).Abs() <= 1e-8)
                    et = 3;
                if ((X.y - (l_radial / 2.0)).Abs() <= 1e-8)
                    et = 4;
                // if((X.x + (l_azimuthal / 2.0)).Abs() <= 1e-8)
                //     et = 5;
                // if((X.x - (l_azimuthal / 2.0)).Abs() <= 1e-8)
                //     et = 6;

                return et;
            });

            return grd;
        }


        public static XNSE_Control RotatingDiskBoundaryLayer(int k = 3, int gridRes = 2) {

            var C = new XNSE_Control();

            C.SetDGdegree(k);


            double radiusOP = 100; // operating point -> Re = radiusOp / Lstar
            double viscosity = 1.0; // kinematic viscosity
            double omega = 1.0; // rotation rate
            double Lstar = Math.Sqrt(viscosity / omega);

            double z = 10;
            double zstar = z * Lstar;


            // physical parameters
            double density = 1.0;
            C.PhysicalParameters.rho_A = density;
            C.PhysicalParameters.mu_A = density * viscosity;

            C.PhysicalParameters.IncludeConvection = true;


            // set grid
            double l_azimuthal = zstar / 5.0;
            double l_radial = zstar / 10.0;
            double h_axial = zstar;

            int res_azimuthal = 2 * gridRes;
            int res_radial = 1 * gridRes;
            int res_axial = 10 * gridRes;

            Grid3D grd = RotatingDiskSector_Linearized(radiusOP, l_azimuthal, l_radial, h_axial, res_azimuthal, res_radial, res_axial);
            C.GridFunc = delegate () { return grd; };



            // boundary condition
            C.AddBoundaryValue("wall_rotatingDisk", "VelocityX#A", X => -(X[1] + radiusOP) * omega);


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            // C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            // C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            // C.dtFixed = Math.PI * 1.0e-3;
            // C.NoOfTimesteps = 2000;

            //{
            //    C.AdaptiveMeshRefinement = true;
            //    int AMRlevel = 2;
            //    C.activeAMRlevelIndicators.Add(new AMRLevelIndicatorLibrary.AMRonBoundary(new byte[] { 1 }) { maxRefinementLevel = AMRlevel });
            //    C.AMR_startUpSweeps = AMRlevel;
            //}

            Func<double[], double> PhiFunc = X => -1.0;
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.Option_LevelSetEvolution = LevelSetEvolution.None;


            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-9;

            return C;
        }



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
            C.activeAMRlevelIndicators.Add(new AMRonBoundary(new byte[] { 3 }) { maxRefinementLevel = 1 });
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

            //C.PostprocessingModules.Add(new MovingContactLineLogging() { LogPeriod = 1, printToConsole = true });

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

            C.PhysicalParameters.betaS_A = 5.0;
            C.PhysicalParameters.betaS_B = 5.0;

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
            double offset = -1.0e-6;

            Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()).Sqrt() - radius);
            //Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()) - radius.Pow2());

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / radius);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            double U0 = 0.0;
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


            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
            //C.activeAMRlevelIndicators.Add(new AMRonBoundary(new byte[] { 3 }) { maxRefinementLevel = 3 });
            C.AMR_startUpSweeps = 3;

            //C.ReInitPeriod = 1;
            //C.InitSignedDistance = true;  // Deprecated


            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 5e-6;
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
