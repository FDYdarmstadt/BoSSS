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
using BoSSS.Solution.NSECommon;
using ilPSP;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.Multigrid;
using ilPSP.Utils;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// A few example configurations.
    /// </summary>
    public static class HardcodedControl {
       public static XNSE_Control SloshingTank(string _DbPath = null, int k = 3) {
            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = false; // _DbPath != null;
            C.ProjectName = "XNSE/SloshingTank";
            C.ProjectDescription = "Multiphase Sloshing tank";

            // DG degrees
            // ==========

            C.AddFieldOption("VelocityX", k, FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("VelocityY", k, FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("SurfaceForceDiagnosticX", SaveOpt:FieldOpts.SaveToDBOpt.FALSE);
            C.AddFieldOption("SurfaceForceDiagnosticY", SaveOpt: FieldOpts.SaveToDBOpt.FALSE);
            C.AddFieldOption("Pressure", k - 1, SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("PhiDG", SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("Phi", Degree: 2, SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("Curvature", Degree: 6, SaveOpt: FieldOpts.SaveToDBOpt.TRUE);



            // grid and boundary conditions
            // ============================

            const bool xPeriodic = false;

            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-2, 2, 20);
                double[] Ynodes = GenericBlas.Linspace(-2, 2, 20);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);
                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - (-2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (+2.0)) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic && Math.Abs(X[0] - (-2.0)) <= 1.0e-8)
                        et = 3;
                    if (!xPeriodic && Math.Abs(X[0] - (+2.0)) <= 1.0e-8)
                        et = 4;


                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (x, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (x, t) => 0.0);
            C.AddBoundaryCondition("wall_left", "VelocityX#A", (x, t) => 0.0);
            C.AddBoundaryCondition("wall_right", "VelocityX#A", (x, t) => 0.0);

            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (x, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (x, t) => 0.0);
            C.AddBoundaryCondition("wall_left", "VelocityX#B", (x, t) => 0.0);
            C.AddBoundaryCondition("wall_right", "VelocityX#B", (x, t) => 0.0);



            // Initial Values
            // ==============

            /*
            C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(
                new Guid("..."),
                950);
            */


            C.InitialValues_Evaluators.Add("Phi",
                (X => (X[1] + Math.Sin(X[0] * Math.PI / 2.0) * 0.1))   // quadratic form
                );
            C.InitialValues_Evaluators.Add("GravityX", X => 0.0);
            C.InitialValues_Evaluators.Add("GravityY", X => -1000.0); // cm / sec^2


            // Physical Parameters
            // ===================


            /*
            // Air-Water (lenght scale == meters, 3D space)
            C.PhysicalParameters.rho_A = 1000; //     kg / m³
            C.PhysicalParameters.rho_B = 1.2; //      kg / m³
            C.PhysicalParameters.mu_A = 1.0e-3; //    kg / m / sec
            C.PhysicalParameters.mu_B = 17.1e-6; //   kg / m / sec
            C.PhysicalParameters.Sigma = 72.75e-3; // kg / sec²     
             */

            /*
            // Water (A) vs. Air (B) // (lenght scale == centimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-3; //     kg / cm³
            C.PhysicalParameters.rho_B = 1.2e-6; //   kg / cm³
            C.PhysicalParameters.mu_A = 1e-5; //      kg / cm / sec
            C.PhysicalParameters.mu_B = 17.1e-8; //   kg / cm / sec
            C.PhysicalParameters.Sigma =  72.75e-3; // kg / sec²     
            //*/

            // Water (A) vs. Octane (B) // (lenght scale == centimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-3; //     kg / cm³
            C.PhysicalParameters.rho_B = 0.7e-3; //   kg / cm³
            C.PhysicalParameters.mu_A = 1e-5; //      kg / cm / sec
            C.PhysicalParameters.mu_B = 0.5e-5; //   kg / cm / sec
            C.PhysicalParameters.Sigma = 0;  // no surface tension  




            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.option_solver = "direct";
            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.PressureBlockPrecondMode = MultigridOperator.Mode.IdMass_DropIndefinite;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            C.ComputeEnergy = true;

            // Timestepping
            // ============

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 0.5;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e12;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }

        /// <summary>
        /// Maintainer: kummer
        /// </summary>
        public static XNSE_Control TransientDroplet(
            //string _DbPath = @"\\fdyprime\userspace\kummer\BoSSS-db-XNSE",
            string _DbPath = null,
            int degree = 2,
            double dt = 2e-4,
            double elipsDelta = 0.1,
            int NoOfTs = 100000) {


            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Multiphase Droplet";
            C.Tags.Add("oscillating");
            C.Tags.Add("fourier");

            // DG degrees
            // ==========

            C.SetFieldOptions(degree, 2);

            // grid and boundary conditions
            // ============================

            const double BaseSize = 1.0;
            const bool xPeriodic = false;
            const double VelXBase = 0.0;

            int xkelem = 54;
            int ykelem = 54;
            double xSize = -4.5 * BaseSize;
            double ySize = -4.5 * BaseSize;

            double hMin = Math.Min(2 * xSize / (xkelem), 2 * ySize / (ykelem));


            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-4.5 * BaseSize, 4.5 * BaseSize, 55);
                double[] Ynodes = GenericBlas.Linspace(-4.5 * BaseSize, 4.5 * BaseSize, 55);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);
                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - (-4.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (+4.5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic && Math.Abs(X[0] - (-4.5 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (!xPeriodic && Math.Abs(X[0] - (+4.5 * BaseSize)) <= 1.0e-8)
                        et = 4;


                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (x, t) => VelXBase);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (x, t) => VelXBase);
            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (x, t) => VelXBase);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (x, t) => VelXBase);
            if (!xPeriodic) {
                C.AddBoundaryCondition("wall_left", "VelocityX#A", (x, t) => VelXBase);
                C.AddBoundaryCondition("wall_right", "VelocityX#A", (x, t) => VelXBase);
                C.AddBoundaryCondition("wall_left", "VelocityX#B", (x, t) => VelXBase);
                C.AddBoundaryCondition("wall_right", "VelocityX#B", (x, t) => VelXBase);

            }

            // Initial Values
            // ==============

            //var database = new DatabaseInfo(_DbPath);

            //var latestSession = database.Sessions.OrderByDescending(e => e.CreationTime)
            //    .First(sess => sess.ProjectName == "XNSE/Droplet" && (sess.Timesteps.Last().TimeStepNumber.MajorNumber - sess.Timesteps.First().TimeStepNumber.MajorNumber) > 50);

            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(latestSession.ID, null);

            //ISessionInfo latestSession = database.Sessions.OrderByDescending(e => e.CreationTime)
            //    .FirstOrDefault(sess => sess.ProjectName == "XNSE/Droplet" && sess.Timesteps.Count > 50 && sess.Tags.Contains("highPenalty"));

            //if (latestSession == null) {
            //    C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(new Guid("6c567939-b310-44a8-b45c-d94880e04cbf"), new TimestepNumber(500));
            //} else {
            //    C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(latestSession.ID, null);
            //}





            double radius = 0.835;

            C.InitialValues_Evaluators.Add("Phi",
                (X => (X[0] / (radius * BaseSize * (1.0 + elipsDelta))).Pow2() + (X[1] / (radius * BaseSize * (1.0 - elipsDelta))).Pow2() - 1.0)   // quadratic form
                );
            C.InitialValues_Evaluators.Add("VelocityX", X => VelXBase);

            // Physical Parameters
            // ===================


            // Air-Water (lenght scale == centimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-3; //     kg / cm³
            C.PhysicalParameters.rho_B = 1.2e-6; //   kg / cm³
            C.PhysicalParameters.mu_A = 1e-5; //      kg / cm / sec
            C.PhysicalParameters.mu_B = 17.1e-8; //   kg / cm / sec
            C.PhysicalParameters.Sigma = 72.75e-3; // kg / sec²     
            //*/

            /*
            // Air-Water (lenght scale == centimeters, 2D space i.e. pressure = Force/Len, density = mass/Len/Len, etc.)
            // Dimensions are different, therefore different scaling. hovever, results scale the same way.
            C.PhysicalParameters.rho_A = 1e-1; //     kg / cm²
            C.PhysicalParameters.rho_B = 1.2e-4; //   kg / cm²
            C.PhysicalParameters.mu_A = 1e-3; //      kg / sec
            C.PhysicalParameters.mu_B = 17.1e-6; //   kg / sec
            C.PhysicalParameters.Sigma = 72.75e-1; // kg cm / sec²     
            */

            /*
            // Air-Water (lenght scale == millimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-6; //     kg / mm³
            C.PhysicalParameters.rho_B = 1.2e-9; //   kg / mm³
            C.PhysicalParameters.mu_A = 1e-6; //      kg / mm / sec
            C.PhysicalParameters.mu_B = 17.1e-9; //   kg / mm / sec
            C.PhysicalParameters.Sigma = 72.75e-3; // kg / sec²     
            */

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            // misc. solver options
            // ====================

            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;

            C.Solver_ConvergenceCriterion = 1.0e-6;
            C.LevelSet_ConvergenceCriterion = 1.0e-6;


            bool useFourierLevelSet = false;
            if (useFourierLevelSet) {
                // Fourier -- level-set

                Func<double, double> radius_of_alpha = delegate (double alpha) {
                    double ret = radius + Math.Cos(2.0 * alpha) * radius * elipsDelta;
                    return ret;
                };

                C.FourierLevSetControl = new FourierLevSetControl() { 
                    FType = FourierType.Polar,
                    numSp = 1024,
                    DomainSize = 2 * Math.PI,
                    PeriodicFunc = radius_of_alpha,
                    //FilterWidth = 0.5,
                    UnderRelax = 0.5,
                    InterpolationType = Interpolationtype.LinearSplineInterpolation};

                C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
                C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
            } else {

                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
                C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
                C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 2;
            }

            C.ComputeEnergy = true;

            // Timestepping
            // ============

            C.CompMode = AppControl._CompMode.Transient;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = NoOfTs;

            // haben fertig...
            // ===============

            return C;
        }

        /// <summary>
        /// An ellipsoidal droplet without level-set motion.
        /// </summary>
        /// <param name="_DbPath"></param>
        /// <param name="sizeFactor">scaling of the droplet.</param>
        /// <param name="degree"></param>
        /// <param name="elipsDelta"></param>
        /// <param name="solver"></param>
        /// <param name="IncludeConvection">true: Navier-Stokes; false: Stokes.</param>
        public static XNSE_Control StaticDroplet(
            string _DbPath = null,
            int sizeFactor = 2, int degree = 3, double elipsDelta = 0.0, string solver = null, bool IncludeConvection = false) {
            XNSE_Control C = new XNSE_Control();

            if (sizeFactor < 1)
                throw new ArgumentOutOfRangeException();
            if (degree < 1)
                throw new ArgumentOutOfRangeException();

            const double VelXBase = 0.0;
            const double BaseSize = 1.0;
            bool xPeriodic = false;
            C.FakePoisson = true;

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Multiphase Droplet";
            C.Tags.Add("h/p conv_study");
            C.Tags.Add("Stokes");
            C.Tags.Add(string.Format("k={0}", degree));
            C.Tags.Add(string.Format("sizeFactor={0}", sizeFactor));

            #region DG Degrees

            C.AddFieldOption("VelocityX", degree, FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("VelocityY", degree, FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("SurfaceForceDiagnosticX", SaveOpt: FieldOpts.SaveToDBOpt.FALSE);
            C.AddFieldOption("SurfaceForceDiagnosticY", SaveOpt: FieldOpts.SaveToDBOpt.FALSE);
            C.AddFieldOption("Pressure", degree - 1, SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("PhiDG", SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("Phi", Degree: 2, SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            C.AddFieldOption("Curvature", Degree: 8, SaveOpt: FieldOpts.SaveToDBOpt.TRUE);
            

            #endregion
            // grid and boundary conditions
            // ============================
            #region Grid and BC's
            double h = 3.0 * BaseSize / sizeFactor;
            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, 9 * sizeFactor + 1);
                double[] Ynodes = GenericBlas.Linspace(-1.5 * BaseSize, 1.5 * BaseSize, 9 * sizeFactor + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);
                //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (+1.5 * BaseSize)) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic && Math.Abs(X[0] - (-1.5 * BaseSize)) <= 1.0e-8)
                        et = 3;
                    if (!xPeriodic && Math.Abs(X[0] - (+1.5 * BaseSize)) <= 1.0e-8)
                        et = 4;


                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (X, t) => VelXBase);
            if (!xPeriodic) {
                C.AddBoundaryCondition("wall_left", "VelocityX#A", (X, t) => VelXBase);
                C.AddBoundaryCondition("wall_right", "VelocityX#A", (X, t) => VelXBase);
                C.AddBoundaryCondition("wall_left", "VelocityX#B", (X, t) => VelXBase);
                C.AddBoundaryCondition("wall_right", "VelocityX#B", (X, t) => VelXBase);

#pragma warning disable 162
                if (VelXBase != 0.0) {
                    C.BoundaryValues["wall_left"].type = IncompressibleBcType.Velocity_Inlet.ToString();
                    C.BoundaryValues["wall_right"].type = IncompressibleBcType.Velocity_Inlet.ToString();
                }
#pragma warning restore 162
            }

            #endregion

            #region Initial Values
            // Initial Values
            // ==============


            double xShift = h * 0.0;
            double yShift = 0.0 * h;

            C.InitialValues_Evaluators.Add("Phi",
                //(X => ((X[0] - xShift) / (0.8 * BaseSize * (1.0 + elipsDelta))).Pow2() + ((X[1] - yShift) / (0.8 * BaseSize * (1.0 - elipsDelta))).Pow2() - 1.0)   // quadratic form
                (X => (0.8 - ((X[0] - xShift).Pow2() + (X[1] - yShift).Pow2()).Sqrt()))
                //(X => (0.8*0.8 - ((X[0] - xShift).Pow2() + (X[1] - yShift).Pow2()))) 
                );
            C.InitialValues_Evaluators.Add("VelocityX", X => VelXBase);


            // Physical Parameters
            // ===================


            // Air-Water
            C.PhysicalParameters.rho_A = 1000;
            C.PhysicalParameters.rho_B = 1.2;
            C.PhysicalParameters.mu_A = 1.0e-3;
            C.PhysicalParameters.mu_B = 17.1e-6;
            C.PhysicalParameters.Sigma = 72.75e-3;
            //*/

            /*
            // Lame
            C.PhysicalParameters.rho_A = 10;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.5;
            C.PhysicalParameters.mu_B = 0.005;
            C.PhysicalParameters.Sigma = 0.1;
            //*/


            /*
            // convective not stable
            C.PhysicalParameters.rho_A = 10;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.05;
            C.PhysicalParameters.mu_B = 0.001;
            C.PhysicalParameters.Sigma = 0.1;
             */

            /*
            // Lamer
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.1;
            */

            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;
            #endregion
            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.UseXDG4Velocity = true;

            C.NoOfMultigridLevels = 3;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.PressureBlockPrecondMode = MultigridOperator.Mode.IdMass_DropIndefinite;
            C.Solver_MaxKrylovDim = 100;

            C.FakePoisson = false;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            if (solver == null) {
                C.option_solver =
                    //"iterativesimple_resmini";
                    //"nonlingmres2+simple";
                    //"nonlingmres+schwarz";
                    //"fixpointiterator";
                    //"flexgmres+simple";
                    //"orthonormalization";
                    "direct";
                //"gmres+schwarz";
                //"gmres+schwarz+coarse";
                //"gmres+multigrid";
                //"nonlingmres+simple";
                //"directsimple";
                //"gmres+simple";
            } else {
                C.option_solver = solver;
            }

            C.Solver_MaxIterations = 1000;

            // NO Timestepping
            // ===============

            C.CompMode = AppControl._CompMode.Steady;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            // haben fertig...
            // ===============

            return C;
        }

        /// <summary>
        /// A very simple example for debugging purposes
        /// </summary>
        /// <param name="_DbPath"></param>
        /// <param name="sizeFactor"></param>
        /// <param name="degree"></param>
        /// <param name="elipsDelta"></param>
        /// <param name="solver"></param>
        /// <param name="IncludeConvection"></param>
        /// <returns></returns>
        public static XNSE_Control Pseudo1D(
            string _DbPath = null,
            int degree = 1, string solver = null, bool IncludeConvection = false) {
            XNSE_Control C = new XNSE_Control();

            if (degree < 1)
                throw new ArgumentOutOfRangeException();

            const double VelXBase = 1.0;

            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/Pseudo1D";
            C.ProjectDescription = "Pseudo 1D calculation";
            C.Tags.Add(string.Format("k={0}", degree));


            // DG degrees
            // ==========

            C.SetFieldOptions(degree, 2);
            
            // grid and boundary conditions
            // ============================

            C.GridFunc = delegate {
                double[] Xnodes = new double[] { -5, -3, -1, 1, 3, 7 };
                double[] Ynodes = new double[] { -1, 1 };
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - (-1)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (+1)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] - (-5)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - (+7)) <= 1.0e-8)
                        et = 4;


                    Debug.Assert(et != 0);
                    return et;
                });

                return grd;
            };

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_left", "VelocityX#A", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_right", "VelocityX#A", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_left", "VelocityX#B", (X, t) => VelXBase);
            C.AddBoundaryCondition("wall_right", "VelocityX#B", (X, t) => VelXBase);

#pragma warning disable 162
            if (VelXBase != 0.0) {
                C.BoundaryValues["wall_left"].type = IncompressibleBcType.Velocity_Inlet.ToString();
                C.BoundaryValues["wall_right"].type = IncompressibleBcType.Velocity_Inlet.ToString();
            }
#pragma warning restore 162


            // Initial Values
            // ==============


            C.InitialValues_Evaluators.Add("Phi", (X => (X[0] - 0.9)));
            C.InitialValues_Evaluators.Add("VelocityX", X => VelXBase);


            // Physical Parameters
            // ===================

            // Lame
            C.PhysicalParameters.rho_A = 10;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.5;
            C.PhysicalParameters.mu_B = 0.005;
            C.PhysicalParameters.Sigma = 0.1;
            //*/


            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;

            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.UseXDG4Velocity = false;


            C.NoOfMultigridLevels = 3;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.Solver_MaxKrylovDim = 100;
            C.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource = Solution.XNSECommon.CurvatureAlgorithms.LevelSetSource.fromDG;

            if (solver == null) {
                C.option_solver =
                    //"iterativesimple_resmini";
                    //"nonlingmres2+simple";
                    //"nonlingmres+schwarz";
                    //"fixpointiterator";
                    //"flexgmres+simple";
                    //"orthonormalization";
                    "direct";
                //"gmres+schwarz";
                //"gmres+schwarz+coarse";
                //"gmres+multigrid";
                //"nonlingmres+simple";
                //"directsimple";
                //"gmres+simple";
            } else {
                C.option_solver = solver;
            }

            C.Solver_MaxIterations = 1000;

            // NO Timestepping
            // ===============

            C.CompMode = AppControl._CompMode.Steady;
            C.dtMax = 0.1;
            C.dtMin = 0.1;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1;

            // haben fertig...
            // ===============

            return C;
        }



        static public XNSE_Control[] StaticDroplet_LinearSolver_ParamStudy() {
            List<XNSE_Control> R = new List<XNSE_Control>();

            foreach (int sizefactor in new int[] { 1, 2, 4, 8 }) {
                foreach (int dgdeg in new int[] { 1, 2, 3 }) {

                    if (dgdeg >= 3 && sizefactor >= 8)
                        continue;

                    foreach (var pcMode in new MultigridOperator.Mode[] { MultigridOperator.Mode.IdMass, MultigridOperator.Mode.SymPart_DiagBlockEquilib }) {

                        var BasicDesc = new Tuple<string, object>[] {
                            new Tuple<string,object>("Grid", 18*sizefactor),
                            new Tuple<string,object>("DG-Degree", dgdeg),
                            new Tuple<string,object>("VelBlockPrecond", pcMode)
                        };


                        /*
                        // "iterativesimple_classic"
                        // =================
                        {
                            var C1 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                            C1.VelocityBlockPrecondMode = pcMode;

                            C1.option_solver = "iterativesimple_classic";
                            C1.Paramstudy_ContinueOnError = true;
                            C1.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("option_solver", C1.option_solver)
                            }.Cat(BasicDesc);
                            R.Add(C1);
                        }

                        // "iterativesimple_resmini"
                        // =================
                        {
                            var C5 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                            C5.VelocityBlockPrecondMode = pcMode;

                            C5.option_solver = "iterativesimple_resmini";
                            C5.Paramstudy_ContinueOnError = true;
                            C5.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("option_solver", C5.option_solver)
                            }.Cat(BasicDesc);
                            R.Add(C5);
                        }

                        // "iterativesimpler"
                        // ==================
                        {
                            var C6 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                            C6.VelocityBlockPrecondMode = pcMode;

                            C6.option_solver = "iterativesimpler";
                            C6.Paramstudy_ContinueOnError = true;
                            C6.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                new Tuple<string,object>("option_solver", C6.option_solver)
                            }.Cat(BasicDesc);
                            R.Add(C6);
                        }
                        */

                        foreach (int kdim in new int[] { 10, 20, 100 }) {

                            /*
                            // orthonormalization
                            // ==================
                            {
                                var C2 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                                C2.VelocityBlockPrecondMode = pcMode;
                                C2.option_solver = "orthonormalization";
                                C2.MaxKrylovDim = kdim;
                                C2.Paramstudy_ContinueOnError = true;
                                C2.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string,object>("option_solver", C2.option_solver),
                                    new Tuple<string,object>("KrylovDim", kdim)
                                }.Cat(BasicDesc);
                                R.Add(C2);
                            }
                             */

                            /*
                            // gmres+schwarz
                            // =============
                            {
                                var C3 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                                C3.VelocityBlockPrecondMode = pcMode;
                                C3.option_solver = "gmres+schwarz";
                                C3.MaxKrylovDim = kdim;
                                C3.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string,object>("option_solver", C3.option_solver),
                                    new Tuple<string,object>("KrylovDim", kdim)
                                }.Cat(BasicDesc);


                                R.Add(C3);
                            }*/


                            // gmres+schwarz
                            // =============
                            {
                                var C7 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                                C7.VelocityBlockPrecondMode = pcMode;
                                C7.option_solver = "gmres+schwarz+coarse";
                                C7.Solver_MaxKrylovDim = kdim;
                                C7.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string,object>("option_solver", C7.option_solver),
                                    new Tuple<string,object>("KrylovDim", kdim)
                                }.Cat(BasicDesc);
                                R.Add(C7);
                            }

                            /*
                            // gmres+simple
                            // ============
                            {
                                var C4 = StaticDroplet(sizeFactor: sizefactor, degree: dgdeg);
                                C4.VelocityBlockPrecondMode = pcMode;
                                C4.option_solver = "gmres+simple";
                                C4.MaxKrylovDim = kdim;
                                C4.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                                    new Tuple<string,object>("option_solver", C4.option_solver),
                                    new Tuple<string,object>("KrylovDim", kdim)
                                }.Cat(BasicDesc);
                                R.Add(C4);
                            }
                             */
                        }

                    }
                }
            }



            return R.ToArray();
        }


        public static XNSE_Control ManufacturedDroplet() {
            XNSE_Control C = new XNSE_Control();


            C.DbPath = null;
            C.savetodb = false;

            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Multiphase Droplet";

            C.SetFieldOptions(3, 4);

            C.GridFunc = delegate {
                //double[] Xnodes = GenericBlas.Linspace(-1.5, 1.5, 18);
                //double[] Ynodes = GenericBlas.Linspace(-1.5, 1.5, 18);
                double[] Xnodes = GenericBlas.Linspace(-2, 2, 7);
                double[] Ynodes = GenericBlas.Linspace(-2, 2, 8);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                grd.EdgeTagNames.Add(1, "velocity_inlet");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });

                return grd;
            };

            const double dt = 0.3;
            const double CC_A = 1.0;
            const double CC_B = 0.0;
            const double RHO_A = 1.2;
            const double RHO_B = 0.1;
            const double MU_A = 0.1;
            const double MU_B = 0.01;
            const double a0 = 1.0; // Parameter fuer Ellipse
            const double b0 = 0.5; // Parameter fuer Ellipse


            C.InitialValues_Evaluators.Add("Phi",
                X => -1.0 + (X[0] / a0).Pow2() + (X[1] / b0).Pow2()
                //X => -1.0 + ((X[0] / a0).Pow2() + (X[1] / b0).Pow2()).Sqrt()
                );
            C.InitialValues_Evaluators.Add("VelocityX", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY", X => 0.0);
            C.InitialValues_Evaluators.Add("GravityX#A", X => -(-2.0 * X[0] * CC_A / RHO_A - (1.0 / dt) * (-X[0])));
            C.InitialValues_Evaluators.Add("GravityX#B", X => -(-2.0 * X[0] * CC_B / RHO_B - (1.0 / dt) * (-X[0])));
            C.InitialValues_Evaluators.Add("GravityY#A", X => +((1.0 / dt) * (+X[1])));
            C.InitialValues_Evaluators.Add("GravityY#B", X => +((1.0 / dt) * (+X[1])));



            C.InitialValues_Evaluators.Add("SurfaceForceX", X => -((CC_A - CC_B) * (1 + X[0].Pow2()) + 2.0 * (MU_A - MU_B)));
            C.InitialValues_Evaluators.Add("SurfaceForceY", X => -((CC_A - CC_B) * (1 + X[0].Pow2()) - 2.0 * (MU_A - MU_B)));

            C.AddBoundaryCondition("velocity_inlet", "VelocityX#A", (X, t) => -X[0]);
            C.AddBoundaryCondition("velocity_inlet", "VelocityY#A", (X, t) => X[1]);
            C.AddBoundaryCondition("velocity_inlet", "VelocityX#B", (X, t) => -X[0]);
            C.AddBoundaryCondition("velocity_inlet", "VelocityY#B", (X, t) => X[1]);


            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;

            C.PhysicalParameters.useArtificialSurfaceForce = true;
            C.PhysicalParameters.rho_A = RHO_A;
            C.PhysicalParameters.rho_B = RHO_B;
            C.PhysicalParameters.mu_A = MU_A;
            C.PhysicalParameters.mu_B = MU_B;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.Standard;


            //C.LevelSetSmoothing = false;
            C.CompMode = AppControl._CompMode.Transient;
            //C.LevelSetOptions.CutCellVelocityProjectiontype = Solution.LevelSetTools.Advection.NonconservativeAdvection.CutCellVelocityProjectiontype.L2_plain;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 2 * dt;
            C.NoOfTimesteps = 1;


            return C;
        }


        public static XNSE_Control TaylorCouette(string _DbPath = null, int k = 3, int sizeFactor = 4) {
            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Multiphase Droplet";
            C.Tags.Add("oscillating");

            // DG degrees
            // ==========

            C.SetFieldOptions(k, 2);

            // grid and boundary conditions
            // ============================

            string innerWallTag = IncompressibleBcType.Velocity_Inlet.ToString() + "_inner";
            string outerWallTag = IncompressibleBcType.Velocity_Inlet.ToString() + "_outer";


            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-2, 2, 8 * sizeFactor + 1);
                double[] Ynodes = GenericBlas.Linspace(-2, 2, 8 * sizeFactor + 1);
                var cutOut = new BoundingBox(new double[] { -0.5, -0.5 }, new double[] { +0.5, +0.5 });
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, CutOuts: cutOut);
                grd.EdgeTagNames.Add(1, innerWallTag);
                grd.EdgeTagNames.Add(2, outerWallTag);

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - (-0.5)) <= 1.0e-8 || Math.Abs(X[0] - (+0.5)) <= 1.0e-8
                        || Math.Abs(X[1] - (-0.5)) <= 1.0e-8 || Math.Abs(X[1] - (+0.5)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - (-2)) <= 1.0e-8 || Math.Abs(X[0] - (+2)) <= 1.0e-8
                        || Math.Abs(X[1] - (-2)) <= 1.0e-8 || Math.Abs(X[1] - (+2)) <= 1.0e-8)
                        et = 2;

                    if (et == 0)
                        throw new ApplicationException("error in DefineEdgeTags");
                    return et;
                });

                return grd;
            };


            // Physical Parameters
            // ===================

            const double Ui = 2;
            const double Ua = 1;
            const double rhoA = 0.1;
            const double rhoB = 1.3;
            const double muA = 0.01;
            const double muB = 0.2;
            const double sigma = 0.9;
            double Ri = Math.Sqrt(2) / 2, Ra = 2, Rm = (Ri + Ra) / 2;



            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.rho_B = rhoB;
            C.PhysicalParameters.mu_A = muA;
            C.PhysicalParameters.mu_B = muB;
            C.PhysicalParameters.Sigma = sigma;

            // Exact solution
            // ==============

            //const double _C1A = 1.119266055, _C1B = 1.328440367, _C2A = .2201834863, _C2B = 0.1100917431e-1, _C3A = .9221025166, _C3B = 0;

            //double _C1A = (2*(12*Ua*muB+5*Ui*muA-9*Ui*muB))/(5*muA+27*muB),
            //    _C1B = (2*(3*Ua*muA+9*Ua*muB-4*Ui*muA))/(5*muA+27*muB),
            //    _C2A = -6*muB*(Ua-3*Ui)/(5*muA+27*muB),
            //    _C2B = -6*muA*(Ua-3*Ui)/(5*muA+27*muB),
            //    _C3A = (108*Ua.Pow2()*muA*muB*rhoB-270*Ua.Pow2()*muB.Pow2()*rhoA+162*Ua.Pow2()*muB.Pow2()*rhoB+60*Ua*Ui*muA.Pow2()*rhoB-240*Ua*Ui*muA*muB*rhoA-144*Ua*Ui*muA*muB*rhoB+324*Ua*Ui*muB.Pow2()*rhoA-50*Ui.Pow2()*muA.Pow2()*rhoA-130*Ui.Pow2()*muA.Pow2()*rhoB+180*Ui.Pow2()*muA*muB*rhoA+25*muA.Pow2()*sigma+270*muA*muB*sigma+729*muB.Pow2()*sigma)/(5*muA+27*muB).Pow2(),
            //    _C3B = 0;

            double
               _C1A = (Ra.Pow2() * Ri * Ui * muA - Ra.Pow2() * Ri * Ui * muB + Ra * Rm.Pow2() * Ua * muB - Ri * Rm.Pow2() * Ui * muA) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA),
               _C1B = (Ra * Ri.Pow2() * Ua * muA - Ra * Ri.Pow2() * Ua * muB + Ra * Rm.Pow2() * Ua * muB - Ri * Rm.Pow2() * Ui * muA) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA),
               _C2A = Ri * Ra * Rm.Pow2() * muB * (Ra * Ui - Ri * Ua) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA),
               _C2B = Ra * Ri * Rm.Pow2() * muA * (Ra * Ui - Ri * Ua) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA),
               _C3A = (1.0 / 2.0) * (-Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA.Pow2() * rhoA - Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA.Pow2() * rhoB + Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muB.Pow2() * rhoA + Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muB.Pow2() * rhoB - 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muB.Pow2() * rhoB + 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA.Pow2() * rhoA + 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow2() * muA * muB * sigma + 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow2() * muA * muB * sigma - 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(4) * muA * muB * sigma - Ra.Pow2() * Rm.Pow(7) * Ua.Pow2() * muB.Pow2() * rhoA + Ra.Pow2() * Rm.Pow(7) * Ua.Pow2() * muB.Pow2() * rhoB - Ri.Pow2() * Rm.Pow(7) * Ui.Pow2() * muA.Pow2() * rhoA + Ri.Pow2() * Rm.Pow(7) * Ui.Pow2() * muA.Pow2() * rhoB - 4 * Ra.Pow(4) * Ri.Pow(4) * muA * muB * sigma - 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow2() * muB.Pow2() * sigma - 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow2() * muA.Pow2() * sigma + 2 * Ra.Pow(4) * Ri.Pow(4) * muA.Pow2() * sigma + 2 * Ra.Pow(4) * Ri.Pow(4) * muB.Pow2() * sigma + 2 * Ra.Pow(4) * Rm.Pow(4) * muB.Pow2() * sigma + 2 * Ri.Pow(4) * Rm.Pow(4) * muA.Pow2() * sigma + 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA.Pow2() * rhoB * Math.Log(Rm) - 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muB.Pow2() * rhoA * Math.Log(Rm) - 4 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muB.Pow2() * rhoA * Math.Log(Rm) + 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muA * muB * rhoB * Math.Log(Rm) - 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muA * muB * rhoB * Math.Log(Rm) + 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA * muB * rhoA * Math.Log(Rm) + 4 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA.Pow2() * rhoB * Math.Log(Rm) - 2 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muA * muB * rhoA + 2 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA * muB * rhoB + 2 * Ra * Ri * Rm.Pow(7) * Ua * Ui * muA * muB * rhoA - 2 * Ra * Ri * Rm.Pow(7) * Ua * Ui * muA * muB * rhoB - 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA * muB * rhoA * Math.Log(Rm) + 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA * muB * rhoA * Math.Log(Rm) - 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA * muB * rhoB * Math.Log(Rm) + 4 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muA * muB * rhoB * Math.Log(Rm) - 4 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA * muB * rhoA * Math.Log(Rm) + 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muB.Pow2() * rhoA * Math.Log(Rm) - 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muA.Pow2() * rhoB * Math.Log(Rm) + 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muB.Pow2() * rhoA * Math.Log(Rm) - 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA.Pow2() * rhoB * Math.Log(Rm) + 2 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA * muB * rhoA + 2 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA.Pow2() * rhoB - 2 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muB.Pow2() * rhoA + 2 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muB.Pow2() * rhoA - 2 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muA * muB * rhoB + 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muA * muB * rhoB - 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA * muB * rhoA - 2 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA.Pow2() * rhoB) / (Rm * (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA).Pow2()),
               _C3B = 0;


            Func<double, double> vA = r => _C1A * r + _C2A / r;
            Func<double, double> vB = r => _C1B * r + _C2B / r;
            Func<double, double> psiA = r => 0.5 * rhoA * _C1A.Pow2() * r.Pow2() - 0.5 * rhoA * _C2A.Pow2() / (r.Pow2()) + 2 * rhoA * _C1A * _C2A * Math.Log(r) + _C3A;
            Func<double, double> psiB = r => 0.5 * rhoB * _C1B.Pow2() * r.Pow2() - 0.5 * rhoB * _C2B.Pow2() / (r.Pow2()) + 2 * rhoB * _C1B * _C2B * Math.Log(r) + _C3B;
            //Func<double,double> _psiA = r => 0.5838609245e-2 * r.Pow2() - .1256228666 / r.Pow2() - .1083300757 * Math.Log(r) + .4668609891;
            //Func<double,double> _psiB = r => .1498764510 * r.Pow2() - 0.4082743166e-2 / r.Pow2() + 0.9894702064e-1 * Math.Log(r);

            //Console.WriteLine("Errors: {0}, {1}", psiA(0.7) - _psiA(0.7), psiB(1.7) - _psiB(1.7));
            //Console.WriteLine("Drücke: {0}, {1}", psiA(Rm), psiB(Rm));


            Func<double[], double, double> UA1 = (X, t) => (-X[1] / X.L2Norm()) * vA(X.L2Norm());
            Func<double[], double, double> UA2 = (X, t) => (+X[0] / X.L2Norm()) * vA(X.L2Norm());
            Func<double[], double, double> UB1 = (X, t) => (-X[1] / X.L2Norm()) * vB(X.L2Norm());
            Func<double[], double, double> UB2 = (X, t) => (+X[0] / X.L2Norm()) * vB(X.L2Norm());

            Func<double[], double, double> PsiA = (X, t) => psiA(X.L2Norm());
            Func<double[], double, double> PsiB = (X, t) => psiB(X.L2Norm());


            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { UA1, UA2 });
            C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { UB1, UB2 });

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", PsiA);
            C.ExactSolutionPressure.Add("B", PsiB);


            // Boundary condition
            // ==================

            C.AddBoundaryCondition(innerWallTag, "VelocityX#A", UA1);
            C.AddBoundaryCondition(innerWallTag, "VelocityY#A", UA2);
            C.AddBoundaryCondition(innerWallTag, "VelocityX#B", (X, t) => double.NaN);
            C.AddBoundaryCondition(innerWallTag, "VelocityY#B", (X, t) => double.NaN);

            C.AddBoundaryCondition(outerWallTag, "VelocityX#A", (X, t) => double.NaN);
            C.AddBoundaryCondition(outerWallTag, "VelocityY#A", (X, t) => double.NaN);
            C.AddBoundaryCondition(outerWallTag, "VelocityX#B", UB1);
            C.AddBoundaryCondition(outerWallTag, "VelocityY#B", UB2);


            // Initial Values
            // ==============


            C.InitialValues_Evaluators.Add("Phi",
                (X => X.L2NormPow2() - Rm.Pow2())  // signed-distance form
                );

            C.InitialValues_Evaluators.Add("VelocityX#A", x => UA1(x, 0));
            C.InitialValues_Evaluators.Add("VelocityY#A", x => UA2(x, 0));
            C.InitialValues_Evaluators.Add("VelocityX#B", x => UB1(x, 0));
            C.InitialValues_Evaluators.Add("VelocityY#B", x => UB2(x, 0));

            C.InitialValues_Evaluators.Add("Pressure#A", x => PsiA(x, 0));
            C.InitialValues_Evaluators.Add("Pressure#B", x => PsiB(x, 0));



            // misc. solver options
            // ====================

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.option_solver = "fixpointiterator";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 3;
            C.Solver_MaxIterations = 20;

            // Timestepping
            // ============

            C.CompMode = AppControl._CompMode.Steady;

            // haben fertig...
            // ===============

            return C;
        }


        public static XNSE_Control CouettePoiseuille(string _DbPath = null, int p = 2) {

            XNSE_Control C = new XNSE_Control();

            //const double pSize = 1.0;
            bool xPeriodic = true;

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/CouettePoiseuille";
            C.ProjectDescription = "Static Multiphase";

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
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // physical parameters
            // ===================
            #region physics

            const double rhoA = 20;
            const double rhoB = 10;
            const double muA = 3;
            const double muB = 1;
            const double sigma = 0.0;

            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.rho_B = rhoB;
            C.PhysicalParameters.mu_A = muA;
            C.PhysicalParameters.mu_B = muB;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid genration
            // ==============
            #region grid

            double H = 1;
            double L = 2 * H;

            int ykelem = 5;
            int xkelem = 2 * ykelem;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L / 2, L / 2, xkelem);
                double[] Ynodes = GenericBlas.Linspace(0, H, ykelem);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);
                grd.EdgeTagNames.Add(1, "Velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "Velocity_inlet_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;

                    return et;
                });
                return grd;
            };

            #endregion


            // boundary conditions
            // ===================
            #region BC

            const double u_w = 0.5;

            C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => u_w);
            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#B", (X, t) => u_w);

            #endregion


            // initial values
            // ==================
            #region init

            double y_s = 3.0 / 5.0;

            C.InitialValues_Evaluators.Add("Phi",
                (X => X[1] - y_s)
                );

            double eps = muA / muB;
            double Re = rhoA * H * u_w / muA;

            Func<double[], double> u0_A = X => u_w * (Math.Exp(Re * X[1]) - 1.0) / (Math.Exp(Re * y_s) - 1.0 - Math.Exp(Re * y_s * (1.0 - eps)) * (Math.Exp(eps * Re * y_s) - Math.Exp(eps * Re)));
            Func<double[], double> u0_B = X => u_w + u_w * Math.Exp(Re * y_s * (1.0 - eps)) * (Math.Exp(eps * Re * X[1]) - Math.Exp(eps * Re)) / (Math.Exp(Re * y_s) - 1.0 - Math.Exp(Re * y_s * (1.0 - eps)) * (Math.Exp(eps * Re * y_s) - Math.Exp(eps * Re)));

            C.InitialValues_Evaluators.Add("VelocityX#A", u0_A);
            C.InitialValues_Evaluators.Add("VelocityX#B", u0_B);

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.ComputeEnergy = false;

            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;
            C.Solver_ConvergenceCriterion = 1e-8;

            #endregion


            // Timestepping
            // ===============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.CompMode = AppControl._CompMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 15;

            #endregion

            return C;
        }


        public static XNSE_Control TranspiratingChannel(string _DbPath = null, int p = 2) {

            XNSE_Control C = new XNSE_Control();

            const int kelem = 9;
            const double pSize = 1;
            bool xPeriodic = true;
            // const double y_0 = 0.1; // has to be also changed in XNSE_SolverMain

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/TranspiratingChannel";
            C.ProjectDescription = "Manufactured Solution";

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
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            double xSize = 2 * pSize;
            double ySize = pSize;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);
                grd.EdgeTagNames.Add(1, "Velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "Velocity_inlet_upper");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "Pressure_Dirichlet_left");
                    grd.EdgeTagNames.Add(4, "Pressure_Dirichlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + ySize) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic) {
                        if (Math.Abs(X[0] + xSize) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion

            // physical parameters
            // ===================
            #region physics

            const double rhoA = 1;
            const double rhoB = 2;
            const double muA = 1;
            const double muB = 2;
            const double sigma = 1;

            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.rho_B = rhoB;
            C.PhysicalParameters.mu_A = muA;
            C.PhysicalParameters.mu_B = muB;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion

            // exact solution
            // ==============
            #region solution

            const double v0 = 0.1;
            const double psi0 = -0.1;

            Func<double, double> a_A = t => (-1 / rhoA) * psi0 * t;
            Func<double, double> a_B = t => (-1 / rhoB) * psi0 * t;
            double b_A = -(psi0 / v0) * (((1 / rhoB) - (1 / rhoA)) / (1 - (muA / muB)));
            double b_B = (muA / muB) * b_A;

            Func<double[], double, double> u_A = (X, t) => a_A(t) + b_A * X[1];
            Func<double[], double, double> u_B = (X, t) => a_B(t) + b_B * X[1];
            Func<double[], double, double> v_0 = (X, t) => v0;

            Func<double[], double, double> psi_0 = (X, t) => psi0 * X[0] + (xSize * psi0 + 1);

            //Func<double, double> phi = t => t * v0;

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { u_A, v_0 });
            C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { u_B, v_0 });

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", psi_0);
            C.ExactSolutionPressure.Add("B", psi_0);

            //double[] X_test = new double[] {0,1};
            //double u_B_test = u_B(X_test,0);
            //Console.WriteLine("u_B = {0}", u_B_test);
            //u_B_test = u_B(X_test, 1);
            //Console.WriteLine("u_B = {0}", u_B_test);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", u_A);
            C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => double.NaN);
            C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityY#A", v_0);
            C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityY#B", (X, t) => double.NaN);


            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => double.NaN);
            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#B", u_B);
            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityY#A", (X, t) => double.NaN);
            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityY#B", v_0);


            if (!xPeriodic) {
                C.AddBoundaryCondition("Pressure_Dirichlet_left", "Pressure#A", psi_0);
                C.AddBoundaryCondition("Pressure_Dirichlet_left", "Pressure#B", psi_0);
                C.AddBoundaryCondition("Pressure_Dirichlet_right", "Pressure#A", psi_0);
                C.AddBoundaryCondition("Pressure_Dirichlet_right", "Pressure#B", psi_0);
            }


            #endregion

            // initial values
            // ==================
            #region init

            C.InitialValues_Evaluators.Add("Phi",
                (X => X[1])
                );
            C.InitialValues_Evaluators.Add("VelocityX#A", X => u_A(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX#B", X => u_B(X, 0));
            C.InitialValues_Evaluators.Add("VelocityY#A", X => v_0(X, 0));
            C.InitialValues_Evaluators.Add("VelocityY#B", X => v_0(X, 0));


            C.InitialValues_Evaluators.Add("Pressure#A", X => psi_0(X, 0));
            C.InitialValues_Evaluators.Add("Pressure#B", X => psi_0(X, 0));


            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;


            C.NoOfMultigridLevels = 1;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            //C.MaxKrylovDim = 100;

            C.option_solver = "direct";
            //C.MaxSolverIterations = 20;

            #endregion

            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;

            #endregion

            return C;
        }

        /*
        /// <summary>
        /// 
        /// </summary>
        /// <remarks>
        /// Maintainer: Florian 
        /// </remarks>
        public static XNSE_Control OscillatingDroplet(int p = 2, int kelem = 54) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = @"\\fdyprime\userspace\kummer\BoSSS-db-XNSE";
            C.savetodb = false;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Oscillating Droplet";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new AppControl.BaseConfig.FieldOpts() {
                Degree = p,
                SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new AppControl.BaseConfig.FieldOpts() {
                Degree = p,
                SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new AppControl.BaseConfig.FieldOpts() {
                Degree = p - 1,
                SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new AppControl.BaseConfig.FieldOpts() {
                SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new AppControl.BaseConfig.FieldOpts() {
                Degree = 4,
                SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new AppControl.BaseConfig.FieldOpts() {
                Degree = 8,
                SaveToDB = AppControl.BaseConfig.FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // grid genration
            // ==============
            #region grid

            double sizeFactor = 1.0;
                        double xSize = sizeFactor * 4.5;
            double ySize = sizeFactor * 4.5;

            
            int xkelem = kelem * 1;
            int ykelem = kelem * 1;

            double hMin = Math.Min(2 * xSize / (xkelem), 2 * ySize / (ykelem));


            C.GridFunc = delegate() {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, ykelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");
                
                grd.DefineEdgeTags(delegate(double[] X) {
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


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_left", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_left", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_right", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_right", "VelocityX#B", (X, t) => 0.0);
 
            #endregion


            // Initial Values
            // ==============
            #region init

            double baseSize = 1.0;
            double elipsDelta = 0.1;
            double radius = 0.835;

            C.InitialValues_Evaluators.Add("Phi",
                //(X => (X[0].Pow2() / (0.8 * 0.8) * 1.25 + X[1].Pow2() / (0.8 * 0.8) * 0.8) - 1.0 + 0.2 * Math.Sin(10 * X[0] * X[1])) // Kartoffel
                (X => (X[0] / (radius * baseSize * (1.0 + elipsDelta))).Pow2() + (X[1] / (radius * baseSize * (1.0 - elipsDelta))).Pow2() - 1.0)   // quadratic form
                //(X => ((X[0] / (radius * baseSize * (1.0 + elipsDelta))).Pow2() + (X[1] / (radius * baseSize * (1.0 - elipsDelta))).Pow2()).Sqrt() - 1.0)  // signed-distance form
                //(X => ((X[0] / (0.8 * BaseSize)).Abs().Pow(1.2) + (X[1] / (0.8 * BaseSize)).Abs().Pow(1.2)) - 1.0.Abs().Pow(1.2))
                );
            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            // Air-Water (lenght scale == centimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-3;      // kg / cm^3
            C.PhysicalParameters.rho_B = 1.2e-6;    // kg / cm^3
            C.PhysicalParameters.mu_A = 1e-5;       // kg / cm * sec
            C.PhysicalParameters.mu_B = 17.1e-8;    // kg / cm * sec
            C.PhysicalParameters.Sigma = 72.75e-3;  // kg / sec^2     
            

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityImplementation = ViscosityImplementation.H;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.EnforceLevelSetContinuity  = true;
            //C.Option_LevelSetEvolution = LevelSetEvolution.Outer_Loop;
            C.Option_LevelSetEvolution = LevelSetEvolution.Inner_Loop;
            C.option_solver = "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 300;
            C.Solver_ConvergenceCriterion = 1.0e-6;
            C.LevelSet_ConvergenceCriterion = 1.0e-6;

            // extension-velocity -- level-set
            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_ClosestPoint;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            // Fourier -- level-set
            C.FourierLevSetControl = new FourierLevSetControl(
                FourierType.Cylindrical,
                1024,
                2 * Math.PI,
                hMin,
                alpha => radius + Math.Cos(alpha)*radius*elipsDelta,
                Interpolationtype.LinearSplineInterpolation,
                0.5);

            C.ComputeEnergy = true;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 2e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1e12;
            C.NoOfTimesteps = 3;

            #endregion

            return C;
        }
        */

        public static XNSE_Control SloshingTank2(string _DbPath = @"\\fdyprime\userspace\smuda\Databases\test_db", int p = 2) {

            XNSE_Control C = new XNSE_Control();

            _DbPath = @"D:\local\local_test_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Tank";
            C.ProjectDescription = "Sloshing Tank";
            C.Tags.Add("freeslip boundary condition");

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
            C.FieldOptions.Add("FilteredVelocityX", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("FilteredVelocityY", new FieldOpts() {
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
                Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            #endregion

            // grid genration
            // ==============
            #region grid

            double L = 1;
            double H = 0.5 * L;

            bool xPeriodic = false;

            int xkelem = 16;
            int ykelem_Interface = 10 * xkelem + 1;
            int ykelem_outer = 2 * xkelem;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L / 2, L / 2, xkelem + 1);
                double[] Ynodes_Interface = GenericBlas.Linspace(-L, L, ykelem_Interface + 1);
                Ynodes_Interface = Ynodes_Interface.GetSubVector(1, Ynodes_Interface.GetLength(0) - 2);
                double[] Ynodes_lower = GenericBlas.Linspace(-(H + L), -L, ykelem_outer + 1);
                double[] Ynodes_upper = GenericBlas.Linspace(L, (H + L), ykelem_outer + 1);
                double[] Ynodes = Ynodes_lower.Concat(Ynodes_Interface).ToArray().Concat(Ynodes_upper).ToArray();

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "freeslip_left");
                    grd.EdgeTagNames.Add(4, "freeslip_right");
                    //grd.EdgeTagNames.Add(3, "wall_left");
                    //grd.EdgeTagNames.Add(4, "wall_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (H + L)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (H + L)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + L / 2) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L / 2) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            if (!xPeriodic) {
                C.AddBoundaryCondition("freeslip_left");
                C.AddBoundaryCondition("freeslip_right");
                //C.AddBoundaryCondition("wall_left");
                //C.AddBoundaryCondition("wall_right");
            }


            #endregion

            // Initial Values
            // ==============
            #region init

            //double sigma = 0.1;
            //Func<double, double> gaussbump = x => ((1 / Math.Sqrt(2 * Math.PI * sigma.Pow2())) * Math.Exp(-x.Pow2() / (2 * sigma.Pow2())));
            Func<double, double> wave = x => 0.01 * Math.Cos(x * 2 * Math.PI);

            C.InitialValues_Evaluators.Add("Phi",
                (X => (X[1] - wave(X[0]))));

            C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e2);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e2);


            //var database = new DatabaseInfo(_DbPath);
            //var sessTank = database.Sessions.Where(s => s.Name.ToLower().Contains("tank"));
            //var latestSession = sessTank.OrderByDescending(e => e.CreationTime).First();
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(latestSession.ID, null);


            #endregion

            // Physical Parameters
            // ===================
            #region physics


            // Air-Water (length scale: meters)
            //C.PhysicalParameters.rho_A = 1000;      // kg / m^3
            //C.PhysicalParameters.rho_B = 1.2;       // kg / m^3
            //C.PhysicalParameters.mu_A = 1.0e-3;     // kg / m * s
            //C.PhysicalParameters.mu_B = 17.1e-6;    // kg / m * s
            //C.PhysicalParameters.Sigma = 72.75e-3;  // kg / s^2     


            // Air-Water (length scale: centimeters)
            C.PhysicalParameters.rho_A = 1e-3;      // kg / cm^3
            C.PhysicalParameters.rho_B = 1.2e-6;    // kg / cm^3
            C.PhysicalParameters.mu_A = 1e-5;       // kg / cm * s
            C.PhysicalParameters.mu_B = 17.1e-8;    // kg / cm * s
            C.PhysicalParameters.Sigma = 72.75e-3;  // kg / s^2     


            // 
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 0.1;
            //C.PhysicalParameters.mu_A = 1e-3;
            //C.PhysicalParameters.mu_B = 1e-4;

            C.PhysicalParameters.Sigma = 0.0;   // free surface boundary condition

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.Solver_MaxIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.PressureBlockPrecondMode = MultigridOperator.Mode.IdMass;
            C.NoOfMultigridLevels = 1;

            C.ComputeEnergy = false;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;

            #endregion


            return C;
        }


        public static XNSE_Control ForcedSloshingTank(string _DbPath = null, int p = 2) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/Tank";
            C.ProjectDescription = "Forced Sloshing Tank";

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
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            double xSize = 100;
            double ySize = 120;

            int xkelem = 20;
            int ykelem = 25;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize / 2, xSize / 2, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, ykelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);


                // oscillation induced by wall motion
                //grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                //grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                //grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                //grd.EdgeTagNames.Add(4, "velocity_inlet_right");

                // oscillation induced by body force
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
                    if (Math.Abs(X[0] + xSize / 2) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xSize / 2) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            // boundary conditions
            // ===================
            #region BC


            //double A = 0.93;
            //double omega = 1 / 1.183;

            // oscillation induced by wall motion
            //C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", (X, t) => A * omega * Math.Sin(omega * t));
            //C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#B", (X, t) => A * omega * Math.Sin(omega * t));
            //C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#A", (X, t) => A * omega * Math.Sin(omega * t));
            //C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#B", (X, t) => A * omega * Math.Sin(omega * t));

            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", (X, t) => A * omega * Math.Sin(omega * t));
            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", (X, t) => A * omega * Math.Sin(omega * t));
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#A", (X, t) => A * omega * Math.Sin(omega * t));
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#B", (X, t) => A * omega * Math.Sin(omega * t));


            // oscillation induced by body force
            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (X, t) => 0.0);

            C.AddBoundaryCondition("wall_left", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_right", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_left", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_right", "VelocityX#B", (X, t) => 0.0);


            #endregion

            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi",
                (X => (X[1] - 50))
                );

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e2);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e2);


            #endregion

            // Physical Parameters
            // ===================
            #region physics

            // Air-Water (length scale: centimeters)
            C.PhysicalParameters.rho_A = 1e-3;      // kg / cm^3
            C.PhysicalParameters.rho_B = 1.2e-6;    // kg / cm^3
            C.PhysicalParameters.mu_A = 1e-5;       // kg / cm * s
            C.PhysicalParameters.mu_B = 17.1e-8;    // kg / cm * s
            C.PhysicalParameters.Sigma = 72.75e-3;  // kg / s^2   

            C.PhysicalParameters.Sigma = 0.0;       // free surface condition

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            //C.option_solver = "GMRES+schwarz+coarse";
            C.option_solver = "direct";
            //C.option_solver = "fixpointiterator";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;

            #endregion


            return C;
        }


        public static XNSE_Control Test_PlanarFourierLS(string _DbPath = @"\\fdyprime\userspace\smuda\Databases\test_db", int p = 2) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null; // _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/PlanarFourierLS_test";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ==============
            #region grid

            double lambda = 1;

            double L = lambda;
            double H = 2 * L;

            int xkelem = 9;
            int ykelem_Interface = 10 * xkelem + 1;
            int ykelem_outer = 1 * xkelem;

            double grdSize = L / (double)xkelem;

            C.GridFunc = delegate () {
                //double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
                //double[] Ynodes = GenericBlas.Linspace(-L / 2, L / 2, xkelem);
                //var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                //grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                //grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                //grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                //grd.EdgeTagNames.Add(4, "velocity_inlet_right");

                //grd.DefineEdgeTags(delegate (double[] X) {
                //    byte et = 0;
                //    if (Math.Abs(X[1] + (L/2)) <= 1.0e-8)
                //        et = 1;
                //    if (Math.Abs(X[1] - (L/2)) <= 1.0e-8)
                //        et = 2;
                //    if (Math.Abs(X[0]) <= 1.0e-8)
                //        et = 3;
                //    if (Math.Abs(X[0] - L) <= 1.0e-8)
                //        et = 4;

                //    return et;
                //});

                double[] Xnodes = GenericBlas.Linspace(0, 6, 2 * xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 3, xkelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "velocity_inlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + 0) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - 3) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - 6) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            double Yvel = 1.0;

            C.AddBoundaryCondition("velocity_inlet_lower", "VelocityY#A", X => Yvel);
            C.AddBoundaryCondition("velocity_inlet_lower", "VelocityY#B", X => Yvel);
            C.AddBoundaryCondition("velocity_inlet_upper", "VelocityY#A", X => Yvel);
            C.AddBoundaryCondition("velocity_inlet_upper", "VelocityY#B", X => Yvel);

            C.AddBoundaryCondition("velocity_inlet_left", "VelocityY#A", X => Yvel);
            C.AddBoundaryCondition("velocity_inlet_left", "VelocityY#B", X => Yvel);
            C.AddBoundaryCondition("velocity_inlet_right", "VelocityY#A", X => Yvel);
            C.AddBoundaryCondition("velocity_inlet_right", "VelocityY#B", X => Yvel);

            //double yVel_max = 1.0;
            //Func<double[], double, double> yVel_seesaw = (X, t) => yVel_max * Math.Sin(2*Math.PI*t);

            //C.AddBoundaryCondition("velocity_inlet_lower", "VelocityY#A", (X,t) => yVel_seesaw(X, -t));
            //C.AddBoundaryCondition("velocity_inlet_upper", "VelocityY#B", yVel_seesaw);

            #endregion


            // Initial Values
            // ==============
            #region init

            double h0 = 0.0;

            double A0 = 0.1; // lambda / 10;
            Func<double, double> PeriodicFunc = x => h0 + A0 * Math.Sin(x * 2 * Math.PI / lambda);

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - 1)); // Math.Sin(X[0]) + Math.Cos(X[0]) + X[0] - (X[1] + 1))); 

            C.InitialValues_Evaluators.Add("VelocityY#A", X => Yvel);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => Yvel);


            //var database = new DatabaseInfo(_DbPath);
            //var sessProjName = database.Sessions.Where(s => s.Name.ToLower().Contains("test"));
            //var latestSession = sessProjName.OrderByDescending(e => e.CreationTime).First();
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(latestSession.ID, null);

            #endregion


            // exact solution
            // ==============

            //C.Phi = (X, t) => (X[1] - (h0 + Yvel * t));

            //C.ExSol_Velocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExSol_Velocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 1.0 });
            //C.ExSol_Velocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 1.0 });

            //C.ExSol_Pressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExSol_Pressure.Add("A", (X, t) => 0.0);
            //C.ExSol_Pressure.Add("B", (X, t) => 0.0);


            // Fourier Level-Set
            // =================
            #region Fourier

            //int numSp = 64;
            //double[] FourierP = new double[numSp];
            //double[] samplP = new double[numSp];
            //for (int sp = 0; sp < numSp; sp++) {
            //    FourierP[sp] = sp * (L / (double)numSp);
            //    samplP[sp] = PeriodicFunc(FourierP[sp]);
            //}

            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, L, FourierP, samplP, 1.0 / (double)xkelem) {
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //};


            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.3;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;

            C.LSContiProjectionMethod = ContinuityProjectionOption.ContinuousDG;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-8;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            //C.AdvancedDiscretizationOptions.surfTensionMode = SurfaceTensionMode.Curvature_Fourier;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Iterative;
            C.LSunderrelax = 0.7;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 50;

            #endregion


            return C;
        }


        public static XNSE_Control Test_PolarFourierLS(string _DbPath = null, int p = 2) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;  //@"D:\local\local_test_db";
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/PolarFourierLS_test";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            //C.Tags.Add("Bubble");
            //C.PhysicalParameters.rho_A = 100;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 24.5;

            //C.Tags.Add("Droplet");
            //C.PhysicalParameters.rho_A = 1000;
            //C.PhysicalParameters.rho_B = 100;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 1;
            //C.PhysicalParameters.Sigma = 0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ==============
            #region grid

            double L = 1.0;
            double ch_fac = 3;
            int k = 20;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L, ch_fac * L, ((int)ch_fac + 1) * k + 1);
                double[] Ynodes = GenericBlas.Linspace(-L, L, 2 * k + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                //grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                //grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                //grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                //grd.EdgeTagNames.Add(4, "velocity_inlet_right");
                grd.EdgeTagNames.Add(1, "freeslip_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (L)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (L)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + (L)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - (ch_fac * L)) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            double Xvel = 1.0;

            //C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#B", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#A", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#B", X => Xvel);

            C.AddBoundaryCondition("freeslip_lower");
            C.AddBoundaryCondition("freeslip_upper");

            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#A", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#B", X => Xvel);

            C.AddBoundaryCondition("pressure_outlet_right");
            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => Xvel);


            #endregion


            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.0, 0.0 };
            //double a = 0.3;
            //double b = 0.15;
            //Func<double, double> radius = phi => a * b / Math.Sqrt(a.Pow2() * Math.Sin(phi).Pow2() + b.Pow2() * Math.Cos(phi).Pow2());
            double radius = 0.25;


            C.InitialValues_Evaluators.Add("Phi",
                //(X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2())   // quadratic form
                (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius)  // signed-distance form
                                                                                                //(X => (X[0].Pow2() / a.Pow2() + X[1].Pow2() / b.Pow2()) - 1)                    // ellipsoid
                );

            C.InitialValues_Evaluators.Add("VelocityX#A", X => Xvel);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => Xvel);

            C.InitialValues_Evaluators.Add("Pressure#A", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityX#A", X => 3); //0.1
            //C.InitialValues_Evaluators.Add("GravityX#B", X => 3);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("7b94348f-5133-48aa-88c0-be625a70ff92");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // exact solution
            // ==============


            C.Phi = ((X, t) => ((X[0] - (center[0] + Xvel * t)).Pow2() + (X[1] - (center[1])).Pow2()).Sqrt() - radius);

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => Xvel, (X, t) => 0.0 });
            C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => Xvel, (X, t) => 0.0 });

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", (X, t) => 0.0);
            C.ExactSolutionPressure.Add("B", (X, t) => 0.0);


            // Fourier Level-Set
            // =================
            #region Fourier

            int numSp = 640;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                samplP[sp] = radius;
            }

            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)k) {
            //    center = center,
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    centerMove = CenterMovement.Reconstructed,
            //};


            #endregion

            // misc. solver options
            // ====================
            #region solver

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;

            C.LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            //C.AdvancedDiscretizationOptions.surfTensionMode = SurfaceTensionMode.Curvature_Fourier;
            //C.FourierLevSetControl.Timestepper = FourierLevelSet_Timestepper.ExplicitEuler;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 4;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
 
            C.CompMode = AppControl._CompMode.Transient;

            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;

            #endregion


            return C;

        }


        public static XNSE_Control Test_ChannelFlow(int degree = 2) {

            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================
            #region db
            C.DbPath = null;
            C.savetodb = false;

            #endregion

            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = degree - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = degree + 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = degree + 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            double L = 2;
            double H = 1;
            bool xPeriodic = false;

            int wall_bc = 3; // 1 = wall, 2 = freeslip, 3 = velocity_inlet, 4 = mixed

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 22);
                double[] Ynodes = GenericBlas.Linspace(0, H, 15);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);
                switch (wall_bc) {
                    case 1: {
                            grd.EdgeTagNames.Add(1, "wall_lower");
                            grd.EdgeTagNames.Add(2, "wall_upper");
                            break;
                        }
                    case 2: {
                            grd.EdgeTagNames.Add(1, "freeslip_lower");
                            grd.EdgeTagNames.Add(2, "freeslip_upper");
                            break;
                        }
                    case 3: {
                            grd.EdgeTagNames.Add(1, "Velocity_inlet_lower");
                            grd.EdgeTagNames.Add(2, "Velocity_inlet_upper");
                            break;
                        }
                    case 4:
                        {
                            grd.EdgeTagNames.Add(1, "freeslip_lower");
                            grd.EdgeTagNames.Add(2, "Velocity_inlet_upper");
                            //grd.EdgeTagNames.Add(1, "Velocity_inlet_lower");
                            //grd.EdgeTagNames.Add(2, "freeslip_upper");
                            break;
                        }
                }

                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "Velocity_inlet_left");
                    grd.EdgeTagNames.Add(4, "pressure_outlet_right");
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

            // physical parameters
            // ===================
            #region physics

            const double rhoA = 1;
            const double rhoB = 1;
            const double muA = 1;
            const double muB = 1;
            const double sigma = 1;

            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.rho_B = rhoB;
            C.PhysicalParameters.mu_A = muA;
            C.PhysicalParameters.mu_B = muB;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double velX = 1.0;

            switch (wall_bc) {
                case 1:
                    {
                        C.AddBoundaryCondition("wall_lower");
                        C.AddBoundaryCondition("wall_upper");
                        break;
                    }
                case 2:
                    {
                        C.AddBoundaryCondition("freeslip_lower");
                        C.AddBoundaryCondition("freeslip_upper");
                        break;
                    }
                case 3:
                    {
                        C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", (X, t) => velX);
                        C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => velX);
                        C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => velX);
                        C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#B", (X, t) => velX);
                        //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", (X, t) => 0.0);
                        //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => 0.0);
                        //C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => 0.0);
                        //C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#B", (X, t) => 0.0);
                        break;
                    }
                case 4:
                    {
                        C.AddBoundaryCondition("freeslip_lower");
                        C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => velX);
                        C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#B", (X, t) => velX);
                        //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", (X, t) => velX);
                        //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => velX);
                        //C.AddBoundaryCondition("freeslip_upper");
                        break;
                    }
            }

            if (!xPeriodic) {
                C.AddBoundaryCondition("Velocity_inlet_left", "VelocityX#A", (X, t) => velX);
                C.AddBoundaryCondition("Velocity_inlet_left", "VelocityX#B", (X, t) => velX);
                C.AddBoundaryCondition("Pressure_outlet_right");
            }


            #endregion


            // initial values
            // ==============
            #region init


            double fx = 1.0;

            if(xPeriodic)
                C.InitialValues_Evaluators.Add("GravityX", X => fx);


            // exact solution for periodic testcase with lower freeslip and upper velocity inlet
            //Func<double[], double, double> u = (X, t) => 0.5 * fx * (H.Pow2() - X[1].Pow2()) + velX;

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { u, (X, t) => 0.0 });

            C.InitialValues_Evaluators.Add("Phi",
                (X => X[0] - 1)
                );
            C.InitialValues_Evaluators.Add("VelocityX#A", X => velX);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => velX);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct"; 
            C.Solver_MaxIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;

            C.ComputeEnergy = false;

            #endregion

            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 100;

            #endregion

            return C;
        }


        public static XNSE_Control Test_StagnationPointFlow(int degree = 2) {

            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================
            #region db
            C.DbPath = null;
            C.savetodb = false;

            #endregion

            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = degree - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            double L = 2;
            double H = 1;

            int kelem = 10;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L/2.0, L/2.0, 2 * kelem +1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "freeslip_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + (L / 2.0)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - (L / 2.0)) <= 1.0e-8)
                        et = 4;
                    return et;
                });

                return grd;
            };

            #endregion

            // physical parameters
            // ===================
            #region physics

            const double rhoA = 1;
            const double rhoB = 1;
            const double muA = 1;
            const double muB = 1;
            const double sigma = 1;

            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.rho_B = rhoB;
            C.PhysicalParameters.mu_A = muA;
            C.PhysicalParameters.mu_B = muB;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double a = 1.0;

            C.AddBoundaryCondition("freeslip_lower");

            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => a * X[0]);
            C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityY#A", (X, t) => - a * X[1]);

            C.AddBoundaryCondition("pressure_outlet_left");
            C.AddBoundaryCondition("Pressure_outlet_right");

            #endregion


            // initial values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi",
                (X => -1)
                );
            C.InitialValues_Evaluators.Add("VelocityX#A", X => a * X[0]);
            C.InitialValues_Evaluators.Add("VelocityY#A", X => - a * X[1]);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.Solver_MaxIterations = 100;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;

            C.ComputeEnergy = false;

            #endregion

            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Steady;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;

            #endregion

            return C;
        }


        public static XNSE_Control Test_channelFlow_fk(int degree = 2) {

            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================
            #region db
            C.DbPath = null;
            C.savetodb = false;

            #endregion

            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = degree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = degree - 1,
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
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            IncompressibleBcType outBc = IncompressibleBcType.Pressure_Outlet;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 10, 33);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, 9);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false, periodicY: false);

                grd.EdgeTagNames.Add(1, IncompressibleBcType.Wall.ToString() + "_upper");
                grd.EdgeTagNames.Add(2, IncompressibleBcType.Wall.ToString() + "_lower");
                grd.EdgeTagNames.Add(3, outBc.ToString());
                grd.EdgeTagNames.Add(4, IncompressibleBcType.Velocity_Inlet.ToString());
               
                grd.DefineEdgeTags(delegate (double[] X) {
                    double x = X[0], y = X[1];

                    if (Math.Abs(y - (+1)) <= 1.0e-8)
                        return 1;
                    if(Math.Abs(y - (-1)) <= 1.0e-8)
                        return 2;
                    if(Math.Abs(x - (10)) <= 1.0e-8)
                        return 3;
                    if (Math.Abs(x - (0)) <= 1.0e-8)
                        return 4;
                    throw new ApplicationException();
                });

                return grd;
            };

            #endregion

            // physical parameters
            // ===================
            #region physics

            const double rhoA = 1;
            const double rhoB = 1;
            const double muA = 0.1;
            const double muB = 0.1;
            const double sigma = 0;

            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.rho_B = rhoB;
            C.PhysicalParameters.mu_A = muA;
            C.PhysicalParameters.mu_B = muB;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityX#A", (X, t) => (1.0 - X[1] * X[1]));
            C.AddBoundaryCondition(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityX#B", (X, t) => (1.0 - X[1] * X[1]));

            C.AddBoundaryCondition(outBc.ToString());

            C.AddBoundaryCondition(IncompressibleBcType.Wall.ToString() + "_upper");

            C.AddBoundaryCondition(IncompressibleBcType.Wall.ToString() + "_lower");

            #endregion


            // initial values
            // ==================
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => X[0] - 2));
            //C.InitialValues_Evaluators.Add("VelocityX#A", (X) => (2.0 / 3.0));
            //C.InitialValues_Evaluators.Add("VelocityX#B", (X) => (2.0 / 3.0));


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.TransposeTermMissing;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;

            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.Solver_MaxIterations = 3;
            C.Solver_ConvergenceCriterion = 1.0e-10;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;

            C.ComputeEnergy = false;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion

            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            C.dtFixed = 0.01;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1;

            #endregion

            return C;
        }


        public static XNSE_Control FallingDropletOnSurface(int p = 2, int kelem = 40, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            _DbPath = @"D:\local\local_test_db";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "falling droplet";

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
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Testcase 1");
            C.PhysicalParameters.rho_A = 1000;
            C.PhysicalParameters.rho_B = 100;
            C.PhysicalParameters.mu_A = 10;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 24.5;


            //C.Tags.Add("Testcase 1 - higher parameters");
            //C.PhysicalParameters.rho_A = 1000;
            //C.PhysicalParameters.rho_B = 10000;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            //C.Tags.Add("Testcase 2");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 1.96;

            // Re = 3.5 ; Bo(Eo) = 1
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            //// Re = 35 ; Bo(Eo) = 100
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 2.45;

            //// Re = 70 ; Bo(Eo) = 10
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.05;
            //C.PhysicalParameters.mu_B = 5;
            //C.PhysicalParameters.Sigma = 24.5;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 2.0;

            //int kelem = 160;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

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

            double[] center = new double[] { 0.5, 1.5 }; //0.5,0.5
            double radius = 0.25;
            Func<double[], double> droplet = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form

            Func<double[], double> surface = (X => (X[1] - 1.0));

            Func<double[], double> PhiFunc = X => Math.Min( droplet(X), surface(X));
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("58745416-3320-4e0c-a5fa-fc3a2c5203c7");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;


            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-3; // (1.0 / (double)kelem) / 16.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 3000; // (int)(3 / dt);
            C.saveperiod = 30;

            #endregion

            return C;

        }


        public static XNSE_Control OscillatingSphere(int p = 1, int kelem = 19, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_OscillatingSphere";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/OscillatingSphere";
            C.ProjectDescription = "static droplet";
            C.Tags.Add("hysing");

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


            // grid genration
            // ==============
            #region grid

            //double L = 4.5;

            //double h_min = L / (double)kelem;

            //C.GridFunc = delegate () {
            //    double[] Xnodes = GenericBlas.Linspace(-L, L, kelem + 1);
            //    double[] Ynodes = GenericBlas.Linspace(-L, L, kelem + 1);
            //    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

            //    grd.EdgeTagNames.Add(1, "wall_lower");
            //    grd.EdgeTagNames.Add(2, "wall_upper");
            //    grd.EdgeTagNames.Add(3, "wall_left");
            //    grd.EdgeTagNames.Add(4, "wall_right");

            //    grd.DefineEdgeTags(delegate (double[] X) {
            //        byte et = 0;
            //        if (Math.Abs(X[1] + L) <= 1.0e-8)
            //            et = 1;
            //        if (Math.Abs(X[1] - L) <= 1.0e-8)
            //            et = 2;
            //        if (Math.Abs(X[0] + L) <= 1.0e-8)
            //            et = 3;
            //        if (Math.Abs(X[0] - L) <= 1.0e-8)
            //            et = 4;
            //        return et;
            //    });
            //    return grd;
            //};

            // smolianski
            double L = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;
                    return et;
                });
                return grd;
            };


            #endregion


            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Bubble");
            //C.PhysicalParameters.rho_A = 100;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 24.5;

            //C.Tags.Add("Droplet");
            //C.PhysicalParameters.rho_A = 1000;
            //C.PhysicalParameters.rho_B = 100;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 1;
            //C.PhysicalParameters.Sigma = 24.5;

            //C.PhysicalParameters.rho_A = 1.0;
            //C.PhysicalParameters.rho_B = 1.0;
            //C.PhysicalParameters.mu_A = 1.0;
            //C.PhysicalParameters.mu_B = 1.0;
            //C.PhysicalParameters.Sigma = 1.0;

            //// Air-Water (lenght scale == centimeters, 3D space)
            //C.PhysicalParameters.rho_A = 1e-3; //     kg / cm³
            //C.PhysicalParameters.rho_B = 1.2e-6; //   kg / cm³
            //C.PhysicalParameters.mu_A = 1e-5; //      kg / cm / sec
            //C.PhysicalParameters.mu_B = 17.1e-8; //   kg / cm / sec
            //C.PhysicalParameters.Sigma = 72.75e-3; // kg / sec²   

            // La = 5000
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0e-4;
            C.PhysicalParameters.Sigma = 1.0;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial Values
            // ==============
            #region init


            double[] center = new double[] { 0.5, 0.5 };
            double a = 2.4;
            double b = 2.4;
            //Func<double, double> radius = phi => a * b / Math.Sqrt(a.Pow2() * Math.Sin(phi).Pow2() + b.Pow2() * Math.Cos(phi).Pow2());
            double radius = 0.25;
            Func<double, double> radiusFunc = phi => radius;

            double delta = 0.0;
            C.InitialValues_Evaluators.Add("Phi",
                //(X => (X[0].Pow2() / a.Pow2() + X[1].Pow2() / b.Pow2()) - 1)
                (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius)  // signed distance
                //(X => ((1 + delta) * (X[0] - center[0])).Pow2() + ((1.0 - delta) * (X[1] - center[1])).Pow2() - radius.Pow2())   // quadratic form
                );

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //double pressureJump = C.PhysicalParameters.Sigma / radius;
            //C.InitialValues_Evaluators.Add("Pressure#A", X => pressureJump);
            //C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("fb857c4c-c060-4d10-a86a-e4ef6a93f5c8");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 180);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("wall_left");
            C.AddBoundaryCondition("wall_right");

            #endregion


            // exact solution
            // ==============
            #region exact

            //C.Phi = ((X, t) => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius);

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => pressureJump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // Fourier Level-Set
            // =================
            #region Fourier

            int numSp = 640;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                samplP[sp] = radius;
            }

            C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, L / (double)kelem) {
                center = center,
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                centerMove = CenterMovement.Reconstructed,
            };


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeInterfaceEnergy = true;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            //C.ContiField = XNSE_Control.ContinuityProjection.ContinuousDG;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 2;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.LaplaceBeltrami_Local;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.CompMode = AppControl._CompMode.Transient;

            double dt = 5e-5;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 125;
            C.NoOfTimesteps = 20;
            C.saveperiod = 10;

            #endregion

            return C;

        }


        public static XNSE_Control CollidingDroplets(int p = 2, int kelem = 40, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            _DbPath = @"D:\local\local_Testcase_databases\Testcase_CollidingDroplets";

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
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
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


                var _xNodes1 = Grid1D.TanhSpacing(-1.5, -0.1, 20 + 1, 1.5, false); 
                _xNodes1 = _xNodes1.GetSubVector(0, (_xNodes1.Length - 1));
                var _xNodes2 = GenericBlas.Linspace(-0.1, 0.1, 10 + 1); 
                _xNodes2 = _xNodes2.GetSubVector(0, (_xNodes2.Length - 1));
                var _xNodes3 = Grid1D.TanhSpacing(0.1, 1.5, 20 + 1, 1.5, true);

                var xNodes = ArrayTools.Cat(_xNodes1, _xNodes2, _xNodes3);


                var _yNodes1 = Grid1D.TanhSpacing(-1.5, 0.0, 24 + 1, 1.5, false);
                _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                //var _yNodes2 = GenericBlas.Linspace(-1.2, 1.2, Convert.ToInt32(40 * MeshFactor)); //40
                //_yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
                var _yNodes3 = Grid1D.TanhSpacing(0.0, 1.5, 24 + 1, 1.5, true);

                var yNodes = ArrayTools.Cat(_yNodes1, _yNodes3);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false);


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

            C.AddBoundaryCondition("pressure_outlet_lower");
            C.AddBoundaryCondition("pressure_outlet_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");


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


            C.PhysicalParameters.rho_A = 1e-2;
            C.PhysicalParameters.rho_B = 1.5e-5;
            C.PhysicalParameters.mu_A = 7e-4;
            C.PhysicalParameters.mu_B = 6e-6;
            C.PhysicalParameters.Sigma = 0.5;

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

            C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.CompMode = AppControl._CompMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            double dt = 1e-4; 
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000; 
            C.saveperiod = 10;

            #endregion

            return C;

        }


        public static XNSE_Control FreeSlipBCTest(string _DbPath = null, int p = 2) {

            XNSE_Control C = new XNSE_Control();


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false;

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
            C.PhysicalParameters.Sigma = 1;
            

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid genration
            // ==============
            #region grid

            double L = 1.0;
            int kelem = 20;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-L, L, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "freeslip_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                //grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                //grd.EdgeTagNames.Add(4, "velocity_inlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (L)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (L)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] + (L)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - (L)) <= 1.0e-8)
                        et = 4;
                    return et;
                });
                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("freeslip_lower");
            C.AddBoundaryCondition("freeslip_upper");
            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => 1.0);
            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", X => 1.0);
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#A", X => 1.0);
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#B", X => 1.0);

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi",
                (X => -1)
                );

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => 2.0 * X[1] * (1 - X[0].Pow2()));
            //C.InitialValues_Evaluators.Add("VelocityY#A", X => -2.0 * X[0] * (1 - X[1].Pow2()));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityX#A", X => 0.1);
            C.InitialValues_Evaluators.Add("GravityX#B", X => 0.1);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.ComputeEnergy = false;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpoint+levelset" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-2; 
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000; 

            #endregion

            return C;
        }


        public static XNSE_Control KelvinHelmholtzInstability(string _DbPath = @"\\fdyprime\userspace\smuda\Databases\test_db", int p = 2) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = _DbPath != null;
            C.ProjectName = "XNSE/Instability";
            C.ProjectDescription = "Kelvin Helmholtz Instability";

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
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            bool xPeriodic = true;

            double xSize = 4.0;
            double ySize = 2.0;

            int xkelem = 41;
            int ykelem = 21;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, ykelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);


                grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                    grd.EdgeTagNames.Add(4, "velocity_inlet_right");
                }


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + ySize) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic) {
                        if (Math.Abs(X[0] + xSize) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion

            // exact solution (viscous potential flow)
            // =======================================
            #region exact

            // Air-Water (length scale: centimeters)
            double rho_l = 1e-3;      // kg / cm^3
            double rho_a = 1.2e-6;    // kg / cm^3
            double mu_l = 1e-5;       // kg / cm * s
            double mu_a = 17.1e-8;    // kg / cm * s
            double sigma = 72.75e-3;  // kg / s^2     

            double U_l = 100.0;     //undisturbed x-velocity for water phase
            double U_a = 600.0;     //undisturbed x-velocity for air phase

            double h_l = ySize;
            double h_a = ySize;

            double A_0 = 0.01;           //(complex) amplitude of inital disturbance   
            double k = 2 * Math.PI;       //wavenumber of disturbance

            double beta_R = -11.476;        //growth rate: beta = beta_R + i*beta_I
            double beta_I = -520.368;

            double A_lR = (A_0 * beta_R) / (k * Math.Sinh(k * h_l));          //complex amplitude for water potential
            double A_lI = (A_0 * (beta_I + k * U_l)) / (k * Math.Sinh(k * h_l));
            double A_aR = -(A_0 * beta_R) / (k * Math.Sinh(k * h_a));          //complex amplitude for air potential
            double A_aI = -(A_0 * (beta_I + k * U_a)) / (k * Math.Sinh(k * h_a));


            Func<double[], double, double> h = (X, t) => 2 * A_0 * Math.Exp(beta_R * t) * Math.Cos(beta_I * t + k * X[0]);
            Func<double[], double, double> u_l = (X, t) => U_l - 2 * k * Math.Exp(beta_R * t) * Math.Cosh(k * (X[1] + h_l)) * (A_lR * Math.Sin(beta_I * t + k * X[0]) + A_lI * Math.Cos(beta_I * t + k * X[0]));
            Func<double[], double, double> w_l = (X, t) => 2 * Math.Exp(beta_R * t) * Math.Sinh(k * (X[1] + h_l)) * (A_aR * Math.Cos(beta_I * t + k * X[0]) - A_aI * Math.Sin(beta_I * t + k * X[0]));
            Func<double[], double, double> u_a = (X, t) => U_a - 2 * k * Math.Exp(beta_R * t) * Math.Cosh(k * (X[1] - h_a)) * (A_aR * Math.Sin(beta_I * t + k * X[0]) + A_aI * Math.Cos(beta_I * t + k * X[0]));
            Func<double[], double, double> w_a = (X, t) => 2 * Math.Exp(beta_R * t) * Math.Sinh(k * (X[1] - h_a)) * (A_aR * Math.Cos(beta_I * t + k * X[0]) - A_aI * Math.Sin(beta_I * t + k * X[0]));

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", u_l);
            C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#B", u_a);
            C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#A", u_l);
            C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#B", u_a);
            if (!xPeriodic) {
                C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", u_l);
                C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", u_a);
                C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#A", u_l);
                C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#B", u_a);

                //C.AddBoundaryCondition("velocity_inlet_left", "Phi", X => (X[1] - h(X, 0)));
                //C.AddBoundaryCondition("velocity_inlet_right", "Phi", X => (X[1] - h(X, 0)));
            }


            #endregion

            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi",
                (X => (X[1] - h(X, 0)))
                );

            C.InitialValues_Evaluators.Add("VelocityX#A", X => u_l(X, 0));
            C.InitialValues_Evaluators.Add("VelocityX#B", X => u_a(X, 0));

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e2);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e2);


            #endregion

            // Physical Parameters
            // ===================
            #region physics


            // Air-Water (length scale: meters)
            //C.PhysicalParameters.rho_A = 1000;      // kg / m^3
            //C.PhysicalParameters.rho_B = 1.2;       // kg / m^3
            //C.PhysicalParameters.mu_A = 1.0e-3;     // kg / m * s
            //C.PhysicalParameters.mu_B = 17.1e-6;    // kg / m * s
            //C.PhysicalParameters.Sigma = 72.75e-3;  // kg / s^2     


            // Air-Water (length scale: centimeters)
            C.PhysicalParameters.rho_A = rho_l;     // kg / cm^3
            C.PhysicalParameters.rho_B = rho_a;     // kg / cm^3
            C.PhysicalParameters.mu_A = mu_l;       // kg / cm * s
            C.PhysicalParameters.mu_B = mu_a;       // kg / cm * s
            C.PhysicalParameters.Sigma = sigma;     // kg / s^2     


            // Air-Oil (length scale: centimeters)
            //C.PhysicalParameters.rho_A = 8.63e-4;
            //C.PhysicalParameters.rho_B = 1.2e-6;
            //C.PhysicalParameters.mu_A = 2e-4;
            //C.PhysicalParameters.mu_B = 17.1e-8;

            //C.PhysicalParameters.Sigma = 0.0;   // free surface boundary condition

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-5;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;

            #endregion

            return C;

        }


        public static XNSE_Control KH_Instability(string _DbPath = @"\\fdyprime\userspace\smuda\Databases\test_db", int p = 2) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Instability";
            C.ProjectDescription = "Kelvin Helmholtz Instability";
            C.Tags.Add("specialized LevelSet");

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            // grid genration
            // ==============
            #region grid

            double L = 2 * Math.PI;

            double h_l = 0.1;
            double h_a = 5 * h_l;

            int xkelem = 60;
            int ykelem = 10 * 6 + 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-h_l, h_a, ykelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");

                //grd.EdgeTagNames.Add(1, "wall_lower");
                //grd.EdgeTagNames.Add(2, "wall_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + h_l) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - h_a) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            // Physical Parameters
            // ===================
            #region physics  

            // Air-Water (length scale: centimeters)
            double rho_l = 1e-3;          // kg / cm^3
            double rho_a = 1.2e-6;        // kg / cm^3
            double mu_l = 1e-5;           // kg / cm * s
            double mu_a = 17.1e-8;        // kg / cm * s
            double sigma = 72.75e-3;      // kg / s^2     


            // Air-Oil (length scale: centimeters)
            //double rho_A = 8.63e-4;
            //double rho_B = 1.2e-6;
            //double mu_A = 2e-4;
            //double mu_B = 17.1e-8;


            C.PhysicalParameters.rho_A = rho_l;          // kg / cm^3
            C.PhysicalParameters.rho_B = rho_a;        // kg / cm^3
            C.PhysicalParameters.mu_A = mu_l;           // kg / cm * s
            C.PhysicalParameters.mu_B = mu_a;       // kg / cm * s
            C.PhysicalParameters.Sigma = sigma;      // kg / s^2     
            //C.PhysicalParameters.Sigma = 0.0;   // free surface boundary condition


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // boundary conditions
            // ===================
            #region BC

            double U_l = 0.0;
            double U_a = 0.0;

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", X => U_l);
            C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#B", X => U_a);

            //C.AddBoundaryCondition("wall_lower", "VelocityX#A", X => 0.0);
            //C.AddBoundaryCondition("wall_upper", "VelocityX#B", X => 0.0);

            #endregion

            // Initial Values
            // ==============
            #region init

            double A0 = 0.005;
            double k = 2;
            Func<double, double> h = x => A0 * Math.Sin(k * x);

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - h(X[0])));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => U_l);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => U_a);

            double g = 9.81e2;
            double alpha = Math.PI / 6;

            C.InitialValues_Evaluators.Add("GravityY#A", X => -g * Math.Cos(alpha));
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g * Math.Cos(alpha));

            C.InitialValues_Evaluators.Add("GravityX#A", X => g * Math.Sin(alpha));
            C.InitialValues_Evaluators.Add("GravityX#B", X => g * Math.Sin(alpha));

            //var database = new DatabaseInfo(_DbPath);
            //var sessTank = database.Sessions.Where(s => s.Name.ToLower().Contains("instability"));
            //var latestSession = sessTank.OrderByDescending(e => e.CreationTime).First();
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(latestSession.ID, null);


            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            //C.PressureBlockPrecondMode = MultigridOperator.Mode.IdMass;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;

            C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Fourier;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;

            #endregion

            return C;
        }


        public static XNSE_Control RisingBubble_fn_BenchmarkTest()
        {

            XNSE_Control C = new XNSE_Control();

            C.LogValues = XNSE_Control.LoggingValues.RisingBubble;

            C.savetodb = false;
            int p = 2;

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
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts()
            {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // grid genration
            // ==============
            #region grid

            bool xPeriodic = false;

            double size = 1.5;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-size, size, 19);
                double[] Ynodes = GenericBlas.Linspace(-size, size, 19);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                if (!xPeriodic)
                {
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + size) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - size) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic)
                    {
                        if (Math.Abs(X[0] + size) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - size) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", (X, t) => 0.0);
            if (!xPeriodic)
            {
                C.AddBoundaryCondition("wall_left", "VelocityX#A", (X, t) => 0.0);
                C.AddBoundaryCondition("wall_right", "VelocityX#A", (X, t) => 0.0);
                C.AddBoundaryCondition("wall_left", "VelocityX#B", (X, t) => 0.0);
                C.AddBoundaryCondition("wall_right", "VelocityX#B", (X, t) => 0.0);
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double x0_c = 0.1;
            double x1_c = -0.3;
            double radius = 0.8;
            Func<double[], double> circle_ls_signd = X => Math.Sqrt((X[0] - x0_c).Pow2() + (X[1] - x1_c).Pow2()) - radius;
            Func<double[], double> circle_ls_quadr = X => ((X[1] - x0_c) / radius).Pow2() + ((X[1] - x1_c) / radius).Pow2() - 1.0;

            C.InitialValues_Evaluators.Add("Phi", (circle_ls_signd));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            #endregion



            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.UseGhostPenalties = false;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpointiterator" : "direct";
            C.Solver_MaxIterations = 50;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;

            C.Solver_ConvergenceCriterion = 1e-8;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion


            C.CompMode = AppControl._CompMode.Steady;


            return C;

        }


        public static XNSE_Control HMF_hangingNodesTest() {

            XNSE_Control C = new XNSE_Control();

            C.savetodb = false;

            int p = 1;

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
                Degree = p + 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(0, 3, 4);
                double[] Ynodes = GenericBlas.Linspace(0, 3, 4);
                var grid = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                //var box_outer_p1 = new double[2] { 0, 0 };
                //var box_outer_p2 = new double[2] { 3, 3 };
                //var box_outer = new GridCommons.GridBox(box_outer_p1, box_outer_p2, 3, 3);

                //var box_inner_p1 = new double[2] { 1, 1 };
                //var box_inner_p2 = new double[2] { 2, 2 };
                //var box_inner = new GridCommons.GridBox(box_inner_p1, box_inner_p2, 2, 2);

                //var grid = Grid2D.HangingNodes2D(box_outer, box_inner);

                grid.EdgeTagNames.Add(1, "wall_lower");
                grid.EdgeTagNames.Add(2, "wall_upper");
                grid.EdgeTagNames.Add(3, "wall_left");
                grid.EdgeTagNames.Add(4, "wall_right");

                grid.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - 3) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - 3) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grid;
            };


            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("wall_left");
            C.AddBoundaryCondition("wall_right");

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - X[0] + 0.2));

            C.ComputeEnergy = false;

            C.CompMode = AppControl._CompMode.Steady;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;


            return C;

        }

    }
}
