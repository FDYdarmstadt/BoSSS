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
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Utils;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.XDG;
using BoSSS.Application.XNSE_Solver.Loadbalancing;
using BoSSS.Application.XNSE_Solver.LoadBalancing;
using MathNet.Numerics.LinearAlgebra.Factorization;
using BoSSS.Solution;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// A few example configurations.
    /// </summary>
    public static class HardcodedControl {


        /// <summary>
        /// Maintainer: kummer
        /// </summary>
        public static XNSE_Control TransientDroplet(
            //string _DbPath = @"\\fdyprime\userspace\kummer\BoSSS-db-XNSE",
            string _DbPath = null,
            int degree = 2,
            double dt = 2e-4,
            double elipsDelta = 0.1,
            int resolution = 54,
            int NoOfTs = 100000) {


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

            C.SetFieldOptions(degree, 2);

            // grid and boundary conditions
            // ============================

            const double BaseSize = 1.0;
            const bool xPeriodic = false;
            const double VelXBase = 0.0;

            int xkelem = resolution;
            int ykelem = resolution;
            double xSize = -4.5 * BaseSize;
            double ySize = -4.5 * BaseSize;

            double hMin = Math.Min(2 * xSize / (xkelem), 2 * ySize / (ykelem));


            C.GridFunc = delegate {
                double[] Xnodes = GenericBlas.Linspace(-4.5 * BaseSize, 4.5 * BaseSize, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-4.5 * BaseSize, 4.5 * BaseSize, ykelem + 1);
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

            C.AddBoundaryValue("wall_lower", "VelocityX#A", (x, t) => VelXBase);
            C.AddBoundaryValue("wall_upper", "VelocityX#A", (x, t) => VelXBase);
            C.AddBoundaryValue("wall_lower", "VelocityX#B", (x, t) => VelXBase);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", (x, t) => VelXBase);
            if (!xPeriodic) {
                C.AddBoundaryValue("wall_left", "VelocityX#A", (x, t) => VelXBase);
                C.AddBoundaryValue("wall_right", "VelocityX#A", (x, t) => VelXBase);
                C.AddBoundaryValue("wall_left", "VelocityX#B", (x, t) => VelXBase);
                C.AddBoundaryValue("wall_right", "VelocityX#B", (x, t) => VelXBase);

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

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            // misc. solver options
            // ====================

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
                    InterpolationType = Interpolationtype.LinearSplineInterpolation
                };

                C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
                C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
            } else {

                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
                C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
                C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 2;
            }

            C.ComputeEnergyProperties = true;

            // Timestepping
            // ============

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = NoOfTs;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
            C.AMR_startUpSweeps = 0;

            // haben fertig...
            // ===============

            return C;
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

            C.AddBoundaryValue("velocity_inlet", "VelocityX#A", (X, t) => -X[0]);
            C.AddBoundaryValue("velocity_inlet", "VelocityY#A", (X, t) => X[1]);
            C.AddBoundaryValue("velocity_inlet", "VelocityX#B", (X, t) => -X[0]);
            C.AddBoundaryValue("velocity_inlet", "VelocityY#B", (X, t) => X[1]);


            C.AgglomerationThreshold = 0.0;

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
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //C.LevelSetOptions.CutCellVelocityProjectiontype = Solution.LevelSetTools.Advection.NonconservativeAdvection.CutCellVelocityProjectiontype.L2_plain;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 2 * dt;
            C.NoOfTimesteps = 1;


            return C;
        }


        /// <summary>
        /// See:
        /// Extended discontinuous Galerkin methods for two-phase flows: the spatial discretization, F. Kummer, IJNME 109 (2), 2017, section 6.3.
        /// </summary>
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

            C.AddBoundaryValue(innerWallTag, "VelocityX#A", UA1);
            C.AddBoundaryValue(innerWallTag, "VelocityY#A", UA2);
            C.AddBoundaryValue(innerWallTag, "VelocityX#B", (X, t) => double.NaN);
            C.AddBoundaryValue(innerWallTag, "VelocityY#B", (X, t) => double.NaN);

            C.AddBoundaryValue(outerWallTag, "VelocityX#A", (X, t) => double.NaN);
            C.AddBoundaryValue(outerWallTag, "VelocityY#A", (X, t) => double.NaN);
            C.AddBoundaryValue(outerWallTag, "VelocityX#B", UB1);
            C.AddBoundaryValue(outerWallTag, "VelocityY#B", UB2);


            // Initial Values
            // ==============


            C.InitialValues_Evaluators.Add("Phi",
                (X => X.L2NormPow2() - Rm.Pow2())  // quadratic form
                );

            C.InitialValues_Evaluators.Add("VelocityX#A", x => UA1(x, 0));
            C.InitialValues_Evaluators.Add("VelocityY#A", x => UA2(x, 0));
            C.InitialValues_Evaluators.Add("VelocityX#B", x => UB1(x, 0));
            C.InitialValues_Evaluators.Add("VelocityY#B", x => UB2(x, 0));

            C.InitialValues_Evaluators.Add("Pressure#A", x => PsiA(x, 0));
            C.InitialValues_Evaluators.Add("Pressure#B", x => PsiB(x, 0));



            // misc. solver options
            // ====================

            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.NonLinearSolver.MaxSolverIterations = 20;
            //C.Solver_MaxIterations = 20;

            // Timestepping
            // ============

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            // haben fertig...
            // ===============

            return C;
        }

        public static XNSE_Control testcube() {
            XNSE_Control C = new XNSE_Control();

            int NoOfTimeSteps = 10;
            int k = 3;
            bool IncludeConvection = true;
            int Res = 15;
            int SpaceDim = 3;
            bool Steady = false;
            bool useLoadBal = true;
            bool useAMR = false;
            var Gshape = Shape.Cube;

            C.GridFunc = delegate {

                int xMin = -2, yMin = -1, zMin = -1;
                int xMax = 4, yMax = 1, zMax = 1;

                // x-direction
                var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                // y-direction
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
                // z-direction
                var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);

                GridCommons grd;
                switch (SpaceDim) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }

                //grd.AddPredefinedPartitioning("debug", MakeDebugPart);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x, y, z;
                    x = X[0];
                    y = X[1];
                    z = X[2];
                    if (Math.Abs(x - xMin) < 1E-8)
                        return 1;
                    else
                        return 3;
                });
                return grd;
            };


            // basic database options
            // ======================
            C.savetodb = false;
            C.SessionName = "test";
            var MachineName = System.Environment.MachineName;
            if (MachineName == "PCMIT32")
                C.DbPath = @"D:\trash_db";
            if (IncludeConvection) {
                C.SessionName += "_NSE";
                C.Tags.Add("NSE");
            } else {
                C.SessionName += "_Stokes";
                C.Tags.Add("Stokes");
            }
            C.Tags.Add(SpaceDim + "D");
            if (Steady) C.Tags.Add("steady");
            else C.Tags.Add("transient");

            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));
            C.saveperiod = 1;
            //C.TracingNamespaces = "*";

            C.GridPartType = GridPartType.clusterHilbert;
            //C.DynamicLoadbalancing_ClassifierType = ClassifierType.CutCells;
            C.DynamicLoadBalancing_On = useLoadBal;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0.1;


            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double Re = 100;
            double muA = 1E-2;

            double partRad = 0.3666;
            double anglev = Re * muA / rhoA / (2 * partRad);
            double d_hyd = 2 * partRad;
            double VelocityIn = Re * muA / rhoA / d_hyd;
            double[] pos = new double[SpaceDim];
            double ts = 2 * Math.PI / anglev / 10;
            double inletdelay = 5 * ts;

            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            // movement of IBM
            C.Rigidbody.SetParameters(pos, anglev, partRad, SpaceDim);
            C.Rigidbody.SpecifyShape(Gshape);
            C.Rigidbody.SetRotationAxis("x");
            C.AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));

            C.AddInitialValue("Pressure", new Formula(@"X => 0"));
            C.AddBoundaryValue("Pressure_Outlet");
            //C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => {VelocityIn}*(double)(t<={inletdelay}?(t/{inletdelay}):1)", true));
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => {VelocityIn}", false));
            C.AddInitialValue("VelocityX", new Formula($"(X) => {VelocityIn}"));

            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1E-6;
            C.NonLinearSolver.MaxSolverIterations = 6;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.SetMaximalRefinementLevel(2);
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============
            double dt = -1;
            if (Steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                dt = 1000;
                C.NoOfTimesteps = 1;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                dt = ts;
                C.NoOfTimesteps = NoOfTimeSteps;
            }
            C.TimeSteppingScheme = TimeSteppingScheme.BDF4;
            C.dtFixed = dt;
            return C;
        }

        public static XNSE_Control testcube_ArtithmeticError() {
            XNSE_Control C = new XNSE_Control();

            int NoOfTimeSteps = 10;
            int k = 3;
            bool IncludeConvection = true;
            int Res = 15;
            int SpaceDim = 3;
            bool Steady = false;
            bool useLoadBal = true;
            bool useAMR = false;
            var Gshape = Shape.Cube;

            C.GridFunc = delegate {

                int xMin = -1, yMin = -1, zMin = -1;
                int xMax = 1, yMax = 1, zMax = 1;
                int Stretch = (xMax - xMin) / (yMax - yMin) * Res;

                // x-direction
                var _xNodes = GenericBlas.Linspace(xMin, xMax, Stretch + 1);
                // y-direction
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
                // z-direction
                var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);

                GridCommons grd;
                switch (SpaceDim) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }

                //grd.AddPredefinedPartitioning("debug", MakeDebugPart);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x, y, z;
                    x = X[0];
                    y = X[1];
                    if (SpaceDim > 2)
                        z = X[2];
                    if (Math.Abs(x - xMin) < 1E-8)
                        return 1;
                    else
                        return 3;
                });
                return grd;
            };


            // basic database options
            // ======================
            C.savetodb = false;
            C.SessionName = "test";
            var MachineName = System.Environment.MachineName;
            if (MachineName == "PCMIT32")
                C.DbPath = @"D:\trash_db";
            if (IncludeConvection) {
                C.SessionName += "_NSE";
                C.Tags.Add("NSE");
            } else {
                C.SessionName += "_Stokes";
                C.Tags.Add("Stokes");
            }
            C.Tags.Add(SpaceDim + "D");
            if (Steady) C.Tags.Add("steady");
            else C.Tags.Add("transient");

            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));
            C.saveperiod = 1;
            //C.TracingNamespaces = "*";

            C.GridPartType = GridPartType.clusterHilbert;
            //C.DynamicLoadbalancing_ClassifierType = ClassifierType.CutCells;
            C.DynamicLoadBalancing_On = useLoadBal;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0.1;


            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double Re = 100;
            double muA = 1E-2;

            double partRad = 0.3666;
            double anglev = Re * muA / rhoA / (2 * partRad);
            double d_hyd = 2 * partRad;
            double VelocityIn = Re * muA / rhoA / d_hyd;
            double[] pos = new double[SpaceDim];
            double ts = 2 * Math.PI / anglev;
            double inletdelay = 5 * ts;

            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            //C.Rigidbody.SetParameters(pos, 0.0, partRad, SpaceDim);
            //C.Rigidbody.SpecifyShape(Gshape);
            //C.Rigidbody.SetRotationAxis("z");
            C.SetOptionsResFields(k);
            C.AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));

            C.AddInitialValue("Pressure", new Formula(@"X => 0"));
            C.AddBoundaryValue("Pressure_Outlet");
            //C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => {VelocityIn}*(double)(t<={inletdelay}?(t/{inletdelay}):1)", true));
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => {VelocityIn}", true));

            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1E-6;
            //C.NonLinearSolver.MaxSolverIterations = 6;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.SetMaximalRefinementLevel(2);
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============
            double dt = -1;
            if (Steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                dt = int.MaxValue;
                C.NoOfTimesteps = 1;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                dt = ts;
                C.NoOfTimesteps = NoOfTimeSteps;
            }
            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.dtFixed = dt;
            return C;
        }

        public static XNSE_Control SphereFlowKrause() {
            XNSE_Control C = new XNSE_Control();

            int NoOfTimeSteps = 10;
            int k = 3;
            bool IncludeConvection = true;
            int SpaceDim = 2;
            bool Steady = false;
            bool useLoadBal = true;
            bool useAMR = false;
            var Gshape = Shape.Sphere;

            C.GridFunc = delegate {

                double xMin = -0.5, yMin = -0.5, zMin = -0.5;
                double xMax = 1.5, yMax = 0.5, zMax = 0.5;

                // x-direction
                var _xNodes = GenericBlas.Linspace(xMin, xMax, 128 + 1);
                // y-direction
                var _yNodes = GenericBlas.Linspace(yMin, yMax, 64 + 1);
                // z-direction
                var _zNodes = GenericBlas.Linspace(zMin, zMax, 13 + 1);

                GridCommons grd;
                switch (SpaceDim) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }

                //grd.AddPredefinedPartitioning("debug", MakeDebugPart);

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x, y, z;
                    x = X[0];
                    y = X[1];
                    if (SpaceDim > 2)
                        z = X[2];
                    if (Math.Abs(x - xMin) < 1E-8)
                        return 1;
                    if (Math.Abs(x - xMax) < 1E-8)
                        return 3;
                    return 2;
                });
                return grd;
            };


            // basic database options
            // ======================
            C.savetodb = false;
            C.SessionName = "test";
            var MachineName = System.Environment.MachineName;
            if (MachineName == "PCMIT32")
                C.DbPath = @"D:\trash_db";
            if (IncludeConvection) {
                C.SessionName += "_NSE";
                C.Tags.Add("NSE");
            } else {
                C.SessionName += "_Stokes";
                C.Tags.Add("Stokes");
            }
            C.Tags.Add(SpaceDim + "D");
            if (Steady) C.Tags.Add("steady");
            else C.Tags.Add("transient");

            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));
            C.saveperiod = 1;
            //C.TracingNamespaces = "*";

            C.GridPartType = GridPartType.clusterHilbert;
            //C.DynamicLoadbalancing_ClassifierType = ClassifierType.CutCells;
            C.DynamicLoadBalancing_On = useLoadBal;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0.1;


            // Physical Parameters
            // ===================
            const double rhoA = 1;
            double muA = 1E-2;

            double partRad = 0.1;
            double[] pos = new double[SpaceDim];

            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            C.Rigidbody.SetParameters(pos, 0.0, partRad, SpaceDim);
            C.Rigidbody.SpecifyShape(Gshape);
            //C.Rigidbody.SetRotationAxis("y");
            C.SetOptionsResFields(k);
            C.AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));

            C.AddInitialValue("Pressure", new Formula(@"X => 0"));
            C.AddBoundaryValue("Pressure_Outlet");
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => 1-X[1]*X[1]", true));
            //C.AddBoundaryValue("Velocity_inlet", "VelocityY", new Formula($"(X,t) => 1-x*x", true));

            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            //C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 1E-6;
            //C.NonLinearSolver.MaxSolverIterations = 6;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.SetMaximalRefinementLevel(2);
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============
            double dt = -1;
            if (Steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                dt = int.MaxValue;
                C.NoOfTimesteps = 1;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                dt = 0.1;
                C.NoOfTimesteps = NoOfTimeSteps;
            }
            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.dtFixed = dt;
            return C;
        }

        public static XNSE_Control RestartTest(
            int NoOfTimeSteps = 10,
            int k = 2,
            bool IncludeConvection = false,
            int Res = 15,
            int SpaceDim = 3,
            bool Steady = false,
            bool useLoadBal = true,
            bool useAMR = false,
            Shape Gshape = Shape.Cube
        ) {
            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            C.savetodb = false;
            var thisOS = System.Environment.OSVersion.Platform;
            var MachineName = System.Environment.MachineName;
            if (MachineName == "PCMIT32")
                C.DbPath = @"D:\trash_db";

            Guid SID = new Guid("7edc2ae0-fb80-4226-b9e2-0fc97f1d4931");
            C.RestartInfo = new Tuple<Guid, TimestepNumber>(SID, new TimestepNumber("1"));

            C.SessionName = "test";
            if (IncludeConvection) {
                C.SessionName += "_NSE";
                C.Tags.Add("NSE");
            } else {
                C.SessionName += "_Stokes";
                C.Tags.Add("Stokes");
            }
            C.Tags.Add(SpaceDim + "D");
            if (Steady) C.Tags.Add("steady");
            else C.Tags.Add("transient");

            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));
            C.saveperiod = 1;
            //C.TracingNamespaces = "*";

            C.GridPartType = GridPartType.clusterHilbert;
            //C.DynamicLoadbalancing_ClassifierType = ClassifierType.CutCells;
            C.DynamicLoadBalancing_On = useLoadBal;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0.1;


            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double Re = 1000;
            double muA = 1E-2;

            double partRad = 0.4333;
            double anglev = Re * muA / rhoA / (2 * partRad);
            //double anglev = 0;
            double d_hyd = 2 * partRad;
            double VelocityIn = Re * muA / rhoA / d_hyd;
            double[] pos = new double[SpaceDim];
            double ts = 2 * Math.PI / anglev / (double)NoOfTimeSteps;
            //double ts = 0.1;
            double inletdelay = 5 * ts;

            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            C.Rigidbody.SetParameters(pos, anglev, partRad, SpaceDim);
            C.Rigidbody.SpecifyShape(Gshape);
            C.Rigidbody.SetRotationAxis("z");
            C.AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));
            C.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi", (X, t) => 0);
            C.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi", (X, t) => 0);

            C.UseImmersedBoundary = true;

            C.AddInitialValue("Pressure", new Formula(@"X => 0"));
            C.AddBoundaryValue("Pressure_Outlet");
            //C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => {VelocityIn}*(double)(t<={inletdelay}?(t/{inletdelay}):1)", true));
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X) => {VelocityIn}"));

            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1E-3;
            C.NonLinearSolver.MaxSolverIterations = 5;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.SetMaximalRefinementLevel(2);
                C.AMR_startUpSweeps = 0;
            }

            // Timestepping
            // ============
            double dt = -1;
            if (Steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                dt = 1000;
                C.NoOfTimesteps = 1;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                dt = ts;
                C.NoOfTimesteps = NoOfTimeSteps;
            }
            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.dtFixed = dt;
            return C;
        }

        public static XNSE_Control Rotating_Cube(int k = 4, int Res = 20, int SpaceDim = 2, bool useAMR = true, int NoOfTimesteps = 2, bool writeToDB = false, bool tracing = false, bool loadbalancing = false, bool IncludeConv = false) {
            double Re = 1000;
            double particleRad = 0.261;

            var C = Rotating_Something(k, Res, SpaceDim, useAMR, NoOfTimesteps, writeToDB, tracing, loadbalancing, Re, particleRad);
            //C.LS_TrackerWidth = 6;
            C.Rigidbody.SpecifyShape(Shape.Cube);
            C.Rigidbody.SetRotationAxis("z");
            C.PhysicalParameters.IncludeConvection = IncludeConv;
            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();

            return C;
        }

        public static XNSE_Control Rotating_Sphere(int k = 4, int Res = 10, int SpaceDim = 3, bool useAMR = false, int NoOfTimesteps = 3, bool writeToDB = true, bool tracing = false, bool loadbalancing = false, bool IncludeConv = false) {
            // --control 'cs:BoSSS.Application.XNSE_Solver.HardcodedControl.Rotating_Sphere(k: 2, Res: 10, SpaceDim: 2, useAMR: false, NoOfTimesteps: 1, writeToDB: false, loadbalancing: false, IncludeConv: false)'


            double particleRad = 0.48;
            double Re = 1000;
            //double anglev = muA/(rhoA*(1- particleRad))*41.3 * Math.Sqrt((particleRad + 1) / (2* (1 - particleRad)));

            var C = Rotating_Something(k, Res, SpaceDim, useAMR, NoOfTimesteps, writeToDB, tracing, loadbalancing, Re, particleRad);
            C.Rigidbody.SpecifyShape(Shape.Sphere);
            C.Rigidbody.SetRotationAxis("z");
            C.PhysicalParameters.IncludeConvection = IncludeConv;

            return C;
        }

        public static XNSE_Control Rotating_Something(int k, int Res, int SpaceDim, bool useAMR, int NoOfTimesteps, bool writeToDB, bool tracing, bool loadbalancing, double Reynoldsnumber, double particleRad) {


            XNSE_Control C = new XNSE_Control();
            // basic database options
            // ======================

            if (writeToDB) {
                var thisOS = System.Environment.OSVersion.Platform;
                var MachineName = System.Environment.MachineName;
                switch (thisOS) {
                    case PlatformID.Unix:
                        C.AlternateDbPaths = new[] {
                            (@"/work/scratch/jw52xeqa/DB_IBM_test", ""),
                            (@"W:\work\scratch\jw52xeqa\DB_IBM_test","")};
                        break;
                    case PlatformID.Win32NT:
                        if (MachineName == "PCMIT32") {
                            C.DbPath = @"D:\trash_db";
                            //C.DbPath = @"D:\2D_Partitioning_samples";
                        } else {
                            C.DbPath = @"\\dc1\userspace\kummer\bosss-db-sept21";
                        }
                        break;
                    default:
                        throw new Exception("No Db-path specified. You stupid?");
                }
                (@"C:\Users\flori\default_bosss_db", "stormbreaker").AddToArray(ref C.AlternateDbPaths);
                (@"C:\Users\flori\default_bosss_db", "stormbreaker").AddToArray(ref C.AlternateDbPaths);
                (@"\\dc1\userspace\kummer\bosss-db-sept21", default(string)).AddToArray(ref C.AlternateDbPaths);
            }
            C.savetodb = writeToDB;
            C.ProjectName = "XNSE/IBM_benchmark";
            C.ProjectDescription = "rotating cube";
            C.Tags.Add(SpaceDim + "D");

            // DG degrees
            // ==========

            //C.SetFieldOptions(k, Math.Max(6, k * 2));
            C.SetFieldOptions(k, Math.Max(k, 2));
            C.SessionName = "XNSE_rotsphere";
            C.saveperiod = 1;
            if (tracing)
                C.TracingNamespaces = "*";

            // grid and boundary conditions
            // ============================

            //// Create Grid
            double xMin = -2, yMin = -2, zMin = -1;
            double xMax = 3, yMax = 2, zMax = 1;

            Func<double[], int> MakeDebugPart = delegate (double[] X) {
                double x = X[0];
                double range = xMax - xMin;
                double interval = range / ilPSP.Environment.MPIEnv.MPI_Size;
                return (int)((x - xMin) / interval);
            };

            C.GridFunc = delegate {

                // x-direction

                var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                //var _xNodes = GenericBlas.Logspace(0, 3, cells_x + 1);
                // y-direction
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
                // z-direction
                var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);


                GridCommons grd;
                switch (SpaceDim) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }

                //grd.AddPredefinedPartitioning("debug", MakeDebugPart);

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x, y, z;
                    x = X[0];
                    y = X[1];
                    if (SpaceDim == 3)
                        z = X[2];
                    if (Math.Abs(x - xMin) < 1E-8)
                        return "Velocity_inlet";
                    else
                        return "Pressure_Outlet";
                });

                return grd;

            };


            //C.GridPartType = GridPartType.Predefined;
            //C.GridPartOptions = "debug";
            C.GridPartType = GridPartType.clusterHilbert;
            C.DynamicLoadBalancing_On = loadbalancing;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimators.Clear();
            C.DynamicLoadBalancing_CellCostEstimators.Add(new Loadbalancing.XNSECellCostEstimator());
            C.DynamicLoadBalancing_ImbalanceThreshold = -0.1;


            //Func<double[], double, double[]> VelocityAtIB = delegate (double[] X, double time) {

            //    if (pos.Length != X.Length)
            //        throw new ArgumentException("check dimension of center of mass");

            //    Vector angVelo = new Vector(new double[] { 0, 0, anglev });
            //    Vector CenterofMass = new Vector(pos);
            //    Vector radialVector = new Vector(X) - CenterofMass;
            //    Vector transVelocity = new Vector(new double[SpaceDim]);
            //    Vector pointVelocity;

            //    switch (SpaceDim) {
            //        case 2:
            //        pointVelocity = new Vector(transVelocity[0] - angVelo[2] * radialVector[1], transVelocity[1] + angVelo[2] * radialVector[0]);
            //        break;
            //        case 3:
            //        pointVelocity = transVelocity + angVelo.CrossProduct(radialVector);
            //        break;
            //        default:
            //        throw new NotImplementedException("this number of dimensions is not supported");
            //    }

            //    return pointVelocity;
            //};

            //Func<double[], double, double> VelocityX = delegate (double[] X, double time) { return VelocityAtIB(X, time)[0]; };
            //Func<double[], double, double> VelocityY = delegate (double[] X, double time) { return VelocityAtIB(X, time)[1]; };
            //Func<double[], double, double> VelocityZ = delegate (double[] X, double time) { return VelocityAtIB(X, time)[2]; };

            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double muA = 1E-3;
            double d_hyd = 2 * particleRad;
            double anglev = Reynoldsnumber * muA / rhoA / d_hyd;
            double VelocityIn = Reynoldsnumber * muA / rhoA / d_hyd;
            double ts = 2 * Math.PI / anglev;
            double inletdelay = 5 * ts;
            Console.WriteLine("Reminder: angular velocity set to zero");
            anglev = 0.0;
            inletdelay = 0.0;
            double[] pos = new double[SpaceDim];
            C.Rigidbody.SetParameters(pos, anglev, particleRad, SpaceDim);

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), X => -1);
            C.UseImmersedBoundary = true;
            //if (C.UseImmersedBoundary) {
            //    C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
            //    //C.InitialValues_EvaluatorsVec.Add(VariableNames.LevelSetCGidx(1), PhiFuncDelegate);
            //    //C.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi2", VelocityX);
            //    //C.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi2", VelocityY);
            //    //if (SpaceDim == 3)
            //    //    C.InitialValues_Evaluators_TimeDep.Add("VelocityZ@Phi2", VelocityZ);
            //}
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", new Formula($"(X,t) => {VelocityIn}*(double)(t<={inletdelay}?(t/{inletdelay}):1)", true));
            C.AddBoundaryValue("Pressure_Outlet");

            //C.OperatorMatrixAnalysis = false;

            // misc. solver options
            // ====================

            //C.EqualOrder = false;
            //C.PressureStabilizationFactor = 1;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            //C.UseSchurBlockPrec = true;
            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            //C.PressureBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.TransposeTermMissing;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1E-6;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;

            C.Tags.Add("convection:" + C.PhysicalParameters.IncludeConvection);

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.SetMaximalRefinementLevel(1);
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            //C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //double dt = 0.01;
            //C.dtMax = dt;
            //C.dtMin = dt*1E-2;
            //C.dtFixed = dt;
            //C.NoOfTimesteps = NoOfTimesteps;


            // haben fertig...
            // ===============

            return C;

        }

        public static XNSE_Control KarmanVortexStreet(int k = 2, int Res = 20, int SpaceDim = 2, int NoOfTimeSteps = 100, bool UseAMR = true, bool writeToDB = false, bool loadbalancing = true) {
            XNSE_Control C = new XNSE_Control();

            // Session Options
            // ===================
            if (writeToDB) {
                var thisOS = System.Environment.OSVersion.Platform;
                var MachineName = System.Environment.MachineName;
                switch (thisOS) {
                    case PlatformID.Unix:
                        C.AlternateDbPaths = new[] {
                            (@"/work/scratch/jw52xeqa/DB_IBM_test", ""),
                            (@"W:\work\scratch\jw52xeqa\DB_IBM_test","")};
                        break;
                    case PlatformID.Win32NT:
                        if (MachineName == "PCMIT32") {
                            C.DbPath = @"D:\trash_db";
                            //C.DbPath = @"D:\2D_Partitioning_samples";
                        } else {
                            C.DbPath = @"\\hpccluster\hpccluster-scratch\weber\DB_IBM_test";
                        }
                        break;
                    default:
                        throw new Exception("No Db-path specified. You stupid?");
                }
                (@"C:\Users\flori\default_bosss_db", "stormbreaker").AddToArray(ref C.AlternateDbPaths);
            }

            C.savetodb = writeToDB;
            C.saveperiod = 1;
            C.ProjectName = "XNSE/IBM_benchmark";
            C.ProjectDescription = "rotating cube";
            C.SessionName = "XNSE_rotsphere";
            C.Tags.Add(SpaceDim + "D");


            // Create Grid (spacial discretization)
            // ===================
            C.SetFieldOptions(k, Math.Max(k, 2));

            double xMin = -1, yMin = -1, zMin = -1;
            double xMax = 3, yMax = 1, zMax = 1;
            double corr = ((xMax - xMin) / (yMax - yMin));

            C.GridFunc = delegate {
                // x-direction
                var _xNodes = GenericBlas.Linspace(xMin, xMax, (int)(Res * corr + 1));
                // y-direction
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
                // z-direction
                var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);
                GridCommons grd;

                switch (SpaceDim) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "Velocity_inlet");
                grd.EdgeTagNames.Add(2, "Wall");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                grd.EdgeTagNames.Add(4, "FreeSlip");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    byte edgeID = 0;
                    double x, y, z;
                    x = X[0];
                    y = X[1];

                    if (Math.Abs(x - xMin) < 1.0e-8)
                        edgeID = 1;
                    if (Math.Abs(x - xMax) < 1.0e-8)
                        edgeID = 3;
                    if (Math.Abs(y - yMax) < 1.0e-8 || Math.Abs(y - yMin) < 1.0e-8)
                        edgeID = 3;
                    if (SpaceDim == 3) {
                        z = X[2];
                        if (Math.Abs(z - zMin) < 1.0e-8 || Math.Abs(z - zMax) < 1.0e-8)
                            edgeID = 3;
                    }
                    return edgeID;
                });
                return grd;
            };

            // Physical Parameters
            // ===================
            const double rhoA = 1; // kg / m3
            const double muA = 1E-3; // Pa*s
            //const double muA = 1; // Pa*s
            double[] pos = new double[SpaceDim];
            const double radius = 0.5;
            const double Re = 700;
            double VeloMag = Re * muA / (2 * radius * rhoA);
            double M = VeloMag / (SpaceDim == 3 ? 1.0 / 3.0 : 2.0 / 3.0);
            C.Tags.Add("Re=" + Re);

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            // Functions and Delegates
            // ====================
            Func<double[], double, double> SphereFunc = delegate (double[] X, double t) {
                switch (SpaceDim) {
                    case 2:
                        // circle
                        return -X[0] * X[0] - X[1] * X[1] + radius * radius;

                    case 3:
                        // sphere
                        return -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + radius * radius;

                    default:
                        throw new NotImplementedException();
                }
            };

            Func<double[], double, double[]> VelocityAtIB = delegate (double[] X, double time) {

                if (pos.Length != X.Length)
                    throw new ArgumentException("check dimension of center of mass");

                Vector angVelo = new Vector(new double[] { 0, 0, 0 });
                Vector CenterofMass = new Vector(pos);
                Vector radialVector = new Vector(X) - CenterofMass;
                Vector transVelocity = new Vector(new double[SpaceDim]);
                Vector pointVelocity;

                switch (SpaceDim) {
                    case 2:
                        pointVelocity = new Vector(transVelocity[0] - angVelo[2] * radialVector[1], transVelocity[1] + angVelo[2] * radialVector[0]);
                        break;
                    case 3:
                        pointVelocity = transVelocity + angVelo.CrossProduct(radialVector);
                        break;
                    default:
                        throw new NotImplementedException("this number of dimensions is not supported");
                }

                return pointVelocity;
            };
            Func<double[], double, double> VelocityX = delegate (double[] X, double time) { return VelocityAtIB(X, time)[0]; };
            Func<double[], double, double> VelocityY = delegate (double[] X, double time) { return VelocityAtIB(X, time)[1]; };
            Func<double[], double, double> VelocityZ = delegate (double[] X, double time) { return VelocityAtIB(X, time)[2]; };

            Func<double, double> timeramp = (double t) => (t <= 1.5 ? (Math.Sin(t - 1.5 / 2) + 1) : 2) / 2;

            // Initial Conditions
            // ====================
            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), X => -1);
            C.UseImmersedBoundary = true;
            if (C.UseImmersedBoundary) {
                C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), SphereFunc);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi2", VelocityX);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi2", VelocityY);
                if (SpaceDim == 3)
                    C.InitialValues_Evaluators_TimeDep.Add("VelocityZ@Phi2", VelocityZ);
            }
            C.InitialValues_Evaluators.Add("Pressure", X => 0);

            // Boundary Conditions
            // ====================
            //C.AddBoundaryValue("FreeSlip");
            //if(SpaceDim==2)
            //    C.AddBoundaryValue("Velocity_inlet","VelocityX",(double[] X)=> M*(1-1/(yMax*yMax)*X[1]*X[1]));
            //if (SpaceDim == 3)
            //    C.AddBoundaryValue("Velocity_inlet", "VelocityX", (double[] X, double t) => M*(1 - 1 / (yMax * yMax) * X[1] * X[1] - 1 / (zMax * zMax) * X[2] * X[2])*(t<=1?t:1));
            C.AddBoundaryValue("Velocity_inlet", "VelocityX", (double[] X, double t) => VeloMag * timeramp(t));
            C.AddBoundaryValue("Velocity_inlet", "VelocityY", (double[] X) => 0);
            C.AddBoundaryValue("Velocity_inlet", "VelocityZ", (double[] X) => 0);
            C.AddBoundaryValue("Pressure_Outlet");

            // misc. solver options
            C.LSContiProjectionMethod = ContinuityProjectionOption.None;


            // ====================
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.None;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1.1123E-6;
            C.NonLinearSolver.MaxSolverIterations = 5;
            C.NonLinearSolver.verbose = true;

            C.Tags.Add("convection:" + C.PhysicalParameters.IncludeConvection);

            // AMR
            // ============
            C.AdaptiveMeshRefinement = UseAMR;
            if (UseAMR) {
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 4 });
                C.AMR_startUpSweeps = 1;
            }

            // Loadbalancing
            // ============
            C.GridPartType = GridPartType.clusterHilbert;
            //C.DynamicLoadbalancing_ClassifierType = ClassifierType.CutCells;
            C.DynamicLoadBalancing_On = loadbalancing;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_CellCostEstimators.Clear();
            C.DynamicLoadBalancing_CellCostEstimators.Add(new Loadbalancing.XNSECellCostEstimator());
            C.DynamicLoadBalancing_ImbalanceThreshold = -0.1;

            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            double dt = 0.1;
            C.dtFixed = dt;
            C.NoOfTimesteps = NoOfTimeSteps;

            return C;
        }


        public static XNSE_Control Rotating_Cube2(int dim = 3, int p = 2, int kelem = 20, bool useAMR = true) {

            XNSE_Control C = new XNSE_Control();

            bool useIB = true;

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            string _DbPath = null;

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "RotatingCube3D";
            C.SessionName = "SetupTest";

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
            if (dim == 3) {
                C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = Math.Max(p, 2),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            if (useIB) {
                C.FieldOptions.Add("PhiDG2", new FieldOpts() {
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
                C.FieldOptions.Add("Phi2", new FieldOpts() {
                    Degree = Math.Max(p, 2),
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }
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
            C.PhysicalParameters.Sigma = 0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            C.GridFunc = delegate () {
                double[] cube = GenericBlas.Linspace(-1.0, 1.0, kelem + 1);
                GridCommons grd;
                if (dim == 3)
                    grd = Grid3D.Cartesian3DGrid(cube, cube, cube);
                else
                    grd = Grid2D.Cartesian2DGrid(cube, cube);

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (1.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (1.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (1.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - (1.0)) <= 1.0e-8)
                        et = 1;
                    if (dim == 3) {
                        if (Math.Abs(X[2] + (1.0)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[2] - (1.0)) <= 1.0e-8)
                            et = 1;
                    }
                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                double[] pos = new double[3];
                double anglev = 10;
                //double t = 0;
                double angle = -(anglev * t) % (2 * Math.PI);
                double particleRad = 0.261;

                if (dim == 3) {
                    return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                            Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)),
                                            Math.Abs(X[2] - pos[2])))
                                            + particleRad;
                } else {
                    return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                        Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                        + particleRad;
                }

            };


            if (useIB) {
                C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), X => -1);
                C.UseImmersedBoundary = true;
                C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
            } else {
                C.InitialValues_Evaluators_TimeDep.Add("Phi", PhiFunc);
            }


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

            C.CheckJumpConditions = false;
            C.CheckInterfaceProps = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 2;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            #endregion


            // level set options
            // ====================
            #region solver

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            if (useIB)
                C.Option_LevelSetEvolution2 = LevelSetEvolution.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            //C.InitSignedDistance = true;


            C.AdaptiveMeshRefinement = useAMR;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
            C.AMR_startUpSweeps = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = compMode;

            if (compMode == AppControl._TimesteppingMode.Transient) {
                double dt = 1;
                C.dtMax = dt;
                C.dtMin = dt;
                C.Endtime = 1000;
                C.NoOfTimesteps = 2;
            }

            #endregion

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

            C.SetDGdegree(p);

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

            C.AddBoundaryValue("Velocity_inlet_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_inlet_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#A", (X, t) => u_w);
            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#B", (X, t) => u_w);

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
            C.ComputeEnergyProperties = false;

            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            #endregion


            // Timestepping
            // ===============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 15;

            #endregion

            return C;
        }

        /// <summary>
        /// 2D manufactured solution from:
        /// 
        /// p‑Multilevel Preconditioners for HHO Discretizations of the Stokes Equations with Static Condensation,
        /// Lorenzo Botti, Daniele A. Di Pietro; https://doi.org/10.1007/s42967-021-00142-5
        /// </summary>
        public static XNSE_Control BottiDiPietro2D(int Res = 20, int p = 2) {
            // --control cs: BoSSS.Application.XNSE_Solver.HardcodedControl.BottiDiPietro2D()
            var C = new XNSE_Control();

            C.GridFunc = delegate () {
                GridCommons g;
                double[] xNodes = GenericBlas.Linspace(-1, +1, Res + 1);
                double[] yNodes = xNodes;
                g = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                g.DefineEdgeTags(delegate (double[] X) {
                    double x = X[0];
                    if (Math.Abs(x - (-1)) < 1e-8)
                        return "pressure_outlet";
                    return "wall";
                });

                return g;
            };


            C.SetDGdegree(p);


            C.PhysicalParameters.rho_A = 1; // not relevant, since density is not present in steady-state Stokes.
            C.PhysicalParameters.rho_B = 1; // not relevant, since density is not present in steady-state Stokes.
            C.PhysicalParameters.mu_A = 1; // dimensionless
            C.PhysicalParameters.mu_B = 1; // dimensionless
            C.PhysicalParameters.Sigma = 0; // not relevant, since single phase
            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            double VelocityXex(double[] X) {
                return -Math.Exp(X[0]) * (X[1] * Math.Cos(X[1]) + Math.Sin(X[1]));
            }

            double VelocityYex(double[] X) {
                return Math.Exp(X[0]) * X[1] * Math.Sin(X[1]);
            }


            C.AddBoundaryValue("wall", "VelocityX", VelocityXex);
            C.AddBoundaryValue("wall", "VelocityY", VelocityYex);

            //C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
            //C.LinearSolver = new PmgConfig() {
            //    ConvergenceCriterion = 1e-9
            //};
            C.LinearSolver = new OrthoMGSchwarzConfig() {
                ConvergenceCriterion = 1e-9
            };


            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.TransposeTermMissing;


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

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

            C.AddBoundaryValue("Velocity_inlet_lower", "VelocityX#A", u_A);
            C.AddBoundaryValue("Velocity_inlet_lower", "VelocityX#B", (X, t) => double.NaN);
            C.AddBoundaryValue("Velocity_inlet_lower", "VelocityY#A", v_0);
            C.AddBoundaryValue("Velocity_inlet_lower", "VelocityY#B", (X, t) => double.NaN);


            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#A", (X, t) => double.NaN);
            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#B", u_B);
            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityY#A", (X, t) => double.NaN);
            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityY#B", v_0);


            if (!xPeriodic) {
                C.AddBoundaryValue("Pressure_Dirichlet_left", "Pressure#A", psi_0);
                C.AddBoundaryValue("Pressure_Dirichlet_left", "Pressure#B", psi_0);
                C.AddBoundaryValue("Pressure_Dirichlet_right", "Pressure#A", psi_0);
                C.AddBoundaryValue("Pressure_Dirichlet_right", "Pressure#B", psi_0);
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

            C.AgglomerationThreshold = 0.2;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;


            #endregion

            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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

            C.AddBoundaryValue("wall_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("wall_left", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_left", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("wall_right", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_right", "VelocityX#B", (X, t) => 0.0);
 
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
            C.ConvergenceCriterion = 1.0e-6;
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

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            if (!xPeriodic) {
                C.AddBoundaryValue("freeslip_left");
                C.AddBoundaryValue("freeslip_right");
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

            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.NonLinearSolver.MaxSolverIterations = 100;

            C.ComputeEnergyProperties = false;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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
            C.AddBoundaryValue("wall_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", (X, t) => 0.0);

            C.AddBoundaryValue("wall_left", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_right", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_left", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("wall_right", "VelocityX#B", (X, t) => 0.0);


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

            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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

            C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", X => Yvel);
            C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#B", X => Yvel);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#A", X => Yvel);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#B", X => Yvel);

            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#A", X => Yvel);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#B", X => Yvel);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityY#A", X => Yvel);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityY#B", X => Yvel);

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

            C.ComputeEnergyProperties = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.3;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;

            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
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

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Iterative;
            C.LSunderrelax = 0.7;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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

            C.AddBoundaryValue("freeslip_lower");
            C.AddBoundaryValue("freeslip_upper");

            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#A", X => Xvel);
            //C.AddBoundaryCondition("velocity_inlet_right", "VelocityX#B", X => Xvel);

            C.AddBoundaryValue("pressure_outlet_right");
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => Xvel);


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
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
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

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 4;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

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
                    case 4: {
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
                case 1: {
                    C.AddBoundaryValue("wall_lower");
                    C.AddBoundaryValue("wall_upper");
                    break;
                }
                case 2: {
                    C.AddBoundaryValue("freeslip_lower");
                    C.AddBoundaryValue("freeslip_upper");
                    break;
                }
                case 3: {
                    C.AddBoundaryValue("Velocity_inlet_lower", "VelocityX#A", (X, t) => velX);
                    C.AddBoundaryValue("Velocity_inlet_lower", "VelocityX#B", (X, t) => velX);
                    C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#A", (X, t) => velX);
                    C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#B", (X, t) => velX);
                    //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", (X, t) => 0.0);
                    //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => 0.0);
                    //C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#A", (X, t) => 0.0);
                    //C.AddBoundaryCondition("Velocity_inlet_upper", "VelocityX#B", (X, t) => 0.0);
                    break;
                }
                case 4: {
                    C.AddBoundaryValue("freeslip_lower");
                    C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#A", (X, t) => velX);
                    C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#B", (X, t) => velX);
                    //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#A", (X, t) => velX);
                    //C.AddBoundaryCondition("Velocity_inlet_lower", "VelocityX#B", (X, t) => velX);
                    //C.AddBoundaryCondition("freeslip_upper");
                    break;
                }
            }

            if (!xPeriodic) {
                C.AddBoundaryValue("Velocity_inlet_left", "VelocityX#A", (X, t) => velX);
                C.AddBoundaryValue("Velocity_inlet_left", "VelocityX#B", (X, t) => velX);
                C.AddBoundaryValue("Pressure_outlet_right");
            }


            #endregion


            // initial values
            // ==============
            #region init


            double fx = 1.0;

            if (xPeriodic)
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
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.ComputeEnergyProperties = false;

            #endregion

            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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
                double[] Xnodes = GenericBlas.Linspace(-L / 2.0, L / 2.0, 2 * kelem + 1);
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

            C.AddBoundaryValue("freeslip_lower");

            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityX#A", (X, t) => a * X[0]);
            C.AddBoundaryValue("Velocity_inlet_upper", "VelocityY#A", (X, t) => -a * X[1]);

            C.AddBoundaryValue("pressure_outlet_left");
            C.AddBoundaryValue("Pressure_outlet_right");

            #endregion


            // initial values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi",
                (X => -1)
                );
            C.InitialValues_Evaluators.Add("VelocityX#A", X => a * X[0]);
            C.InitialValues_Evaluators.Add("VelocityY#A", X => -a * X[1]);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.NonLinearSolver.MaxSolverIterations = 100;
            C.ComputeEnergyProperties = false;

            #endregion

            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            double dt = 0.1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;

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

            Func<double[], double> PhiFunc = X => Math.Min(droplet(X), surface(X));
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

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.NonLinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
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

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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
            //double a = 2.4;
            //double b = 2.4;
            //Func<double, double> radius = phi => a * b / Math.Sqrt(a.Pow2() * Math.Sin(phi).Pow2() + b.Pow2() * Math.Cos(phi).Pow2());
            double radius = 0.25;
            Func<double, double> radiusFunc = phi => radius;

            //double delta = 0.0;
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

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");

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

            C.ComputeEnergyProperties = true;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            //C.ContiField = XNSE_Control.ContinuityProjection.ContinuousDG;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.LaplaceBeltrami_Local;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            double dt = 5e-5;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 125;
            C.NoOfTimesteps = 20;
            C.saveperiod = 10;

            #endregion

            return C;

        }

        /// <summary>
        /// Benchmark. Do not change!
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <param name="D">2D or 3D</param>
        /// <returns></returns>
        public static XNSE_Control StokesSphere(int p = 4, int kelem = 32, string _DbPath = null, int D = 2) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db
            //_DbPath = @"D:\trash_DB";
            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/StokesSphere";
            C.ProjectDescription = "static droplet";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("Velocity*", new FieldOpts() {
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
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = Math.Max(4, 2 * p + 2),
                //Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // grid generation
            // ===============
            #region grid



            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-1, 1, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, kelem + 1);
                double[] Znodes = GenericBlas.Linspace(-1, 1, kelem + 1);

                GridCommons grd;
                switch (D) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }


                grd.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[0] - (-1)) <= 1.0e-8)
                        return "wall_left";
                    if (Math.Abs(X[0] - (+1)) <= 1.0e-8)
                        return "wall_right";
                    if (Math.Abs(X[1] - (-1)) <= 1.0e-8)
                        return "wall_front";
                    if (Math.Abs(X[1] - (+1)) <= 1.0e-8)
                        return "wall_back";
                    if (D > 2) {
                        if (Math.Abs(X[2] - (-1)) <= 1.0e-8)
                            return "wall_top";
                        if (Math.Abs(X[2] - (+1)) <= 1.0e-8)
                            return "wall_bottom";
                    }

                    throw new ArgumentException("unknown wall");
                });
                return grd;
            };


            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Air-Water (lenght scale == centimeters, 3D space)
            C.PhysicalParameters.rho_A = 1e-3; //     kg / cm³
            C.PhysicalParameters.rho_B = 1.2e-6; //   kg / cm³
            C.PhysicalParameters.mu_A = 1e-5; //      kg / cm / sec
            C.PhysicalParameters.mu_B = 17.1e-8; //   kg / cm / sec
            C.PhysicalParameters.Sigma = 72.75e-3; // kg / sec²   



            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial Values
            // ==============
            #region init

            //double[] center = new double[] { 0.5, 0.5 };
            ////double a = 2.4;
            ////double b = 2.4;
            ////Func<double, double> radius = phi => a * b / Math.Sqrt(a.Pow2() * Math.Sin(phi).Pow2() + b.Pow2() * Math.Cos(phi).Pow2());
            //double radius = 0.25;
            //Func<double, double> radiusFunc = phi => radius;

            //double delta = 0.0;

            double r = 0.49;
            double nonsp = 0.5;

            if (D == 2)
                C.AddInitialValue("Phi", new Formula($"X => (X[0]/{r * nonsp}).Pow2() + (X[1]/{r}).Pow2() - 1.0", false));
            else if (D == 3)
                C.AddInitialValue("Phi", new Formula($"X => (X[0]/{r * nonsp}).Pow2() + (X[1]/{r}).Pow2() + (X[2]/{r}).Pow2() - 1.0", false));
            else
                throw new ArgumentOutOfRangeException();

            C.LSContiProjectionMethod = ContinuityProjectionOption.None;


            #endregion


            // boundary conditions
            // ===================
            #region BC

            //C.AddBoundaryValue("wall_lower");
            //C.AddBoundaryValue("wall_upper");
            //C.AddBoundaryValue("wall_left");
            //C.AddBoundaryValue("wall_right");
            //C.AddBoundaryValue("wall_left");
            //C.AddBoundaryValue("wall_right");

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.ComputeEnergyProperties = false;
            C.AgglomerationThreshold = 0.1;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            // Solver related Stuff
            //C.ContiField = XNSE_Control.ContinuityProjection.ContinuousDG;
            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.IdMass_DropIndefinite;
            //C.PressureBlockPrecondMode = MultigridOperator.Mode.IdMass_DropIndefinite;
            C.UseSchurBlockPrec = true;
            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            // Levelset & Co
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.LaplaceBeltrami_Local;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            //C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //C.dtFixed = 1E-4;


            #endregion

            C.SessionName = String.Format("J{0}_p{1}_{2}", Math.Pow(kelem, D), p, C.LinearSolver.Shortname);

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

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

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

            C.AddBoundaryValue("freeslip_lower");
            C.AddBoundaryValue("freeslip_upper");
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


            C.ComputeEnergyProperties = false;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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

            C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", u_l);
            C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", u_a);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", u_l);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", u_a);
            if (!xPeriodic) {
                C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", u_l);
                C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", u_a);
                C.AddBoundaryValue("velocity_inlet_right", "VelocityX#A", u_l);
                C.AddBoundaryValue("velocity_inlet_right", "VelocityX#B", u_a);

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

            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
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

            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => U_l);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => U_a);

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

            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = Solution.XNSECommon.ViscosityMode.FullySymmetric;
            C.NonLinearSolver.MaxSolverIterations = 100;

            C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Fourier;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;

            #endregion

            return C;
        }


        public static XNSE_Control RisingBubble_fn_BenchmarkTest() {

            XNSE_Control C = new XNSE_Control();

            //C.LogValues = XNSE_Control.LoggingValues.RisingBubble;
            C.PostprocessingModules.Add(new PhysicalBasedTestcases.RisingBubble2DBenchmarkQuantities());

            C.savetodb = false;
            int p = 2;

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
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
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
                if (!xPeriodic) {
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + size) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - size) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic) {
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

            C.AddBoundaryValue("wall_lower", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#A", (X, t) => 0.0);
            C.AddBoundaryValue("wall_lower", "VelocityX#B", (X, t) => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", (X, t) => 0.0);
            if (!xPeriodic) {
                C.AddBoundaryValue("wall_left", "VelocityX#A", (X, t) => 0.0);
                C.AddBoundaryValue("wall_right", "VelocityX#A", (X, t) => 0.0);
                C.AddBoundaryValue("wall_left", "VelocityX#B", (X, t) => 0.0);
                C.AddBoundaryValue("wall_right", "VelocityX#B", (X, t) => 0.0);
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

            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 0;

            #endregion


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;


            return C;

        }


        public static XNSE_Control HMF_hangingNodesTest() {

            XNSE_Control C = new XNSE_Control();

            C.savetodb = false;

            int p = 1;

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


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - X[0] + 0.2));

            C.ComputeEnergyProperties = false;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;


            return C;

        }


        public static XNSE_Control HMF_3DcontactlineTest() {

            XNSE_Control C = new XNSE_Control();

            int D = 3;

            if (D == 3)
                C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Steady;

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/HMF3D";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            int p = 2;

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            if (D == 3) {
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

            C.Tags.Add("Reusken");
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            double sigma = 0.0;
            C.PhysicalParameters.Sigma = sigma;

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

            double scale = 0.2;

            int xkelem = 8;
            int ykelem = 8;
            int zkelem = 1;

            double xSize = (scale * (double)xkelem) / 2.0;
            double ySize = (scale * (double)ykelem) / 2.0;
            double zSize = 2.0 * scale * (double)zkelem;

            if (D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, xkelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, ykelem + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");

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
            }

            if (D == 3) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, xkelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, ykelem + 1);
                    double[] Znodes = GenericBlas.Linspace(0, zSize, zkelem + 1);
                    var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    grd.EdgeTagNames.Add(5, "navierslip_linear_front");
                    grd.EdgeTagNames.Add(6, "navierslip_linear_back");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[2]) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[2] - zSize) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0] + xSize) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;
                        if (Math.Abs(X[1] + ySize) <= 1.0e-8)
                            et = 5;
                        if (Math.Abs(X[1] - ySize) <= 1.0e-8)
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

            double R = (scale * 5.0) / 2.0;

            Func<double[], double> PhiFunc = X => -1.0;
            if (D == 2) {
                PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - R);
            }

            if (D == 3) {
                PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - R);  //zylinder
                //PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2() + (X[2] - 0.0).Pow2()).Sqrt() - R);  //sphere
            }

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double pJump = sigma / R;
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

            if (D == 3) {
                C.AddBoundaryValue("navierslip_linear_front");
                C.AddBoundaryValue("navierslip_linear_back");
            }

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.LieSplitting;

            C.TimesteppingMode = compMode;
            double dt = 3e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }



        public static XNSE_Control LateralAdhesionForceGrid() {

            XNSE_Control C = new XNSE_Control();

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Steady;

            // basic database options
            // ======================
            #region db

            C.savetodb = false;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            int p = 2;

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

            C.PhysicalParameters.IncludeConvection = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 1.0;
            double H = 1.0;

            double h = 0.3;
            double b = 0.01;

            C.GridFunc = delegate () {

                var TestSecLeft_p1 = new double[2] { -L, 0 };
                var TestSecLeft_p2 = new double[2] { -b / 2, H };
                var TestSecLeft = new GridCommons.GridBox(TestSecLeft_p1, TestSecLeft_p2, 10, 10);

                var TestSecCenter_p1 = new double[2] { -b / 2, 0 };
                var TestSecCenter_p2 = new double[2] { b / 2, h };
                var TestSecCenter = new GridCommons.GridBox(TestSecCenter_p1, TestSecCenter_p2, 1, 6);

                var TestSecRight_p1 = new double[2] { b / 2, 0 };
                var TestSecRight_p2 = new double[2] { L, H };
                var TestSecRight = new GridCommons.GridBox(TestSecRight_p1, TestSecRight_p2, 10, 10);


                var grd = Grid2D.HangingNodes2D(TestSecLeft, TestSecCenter, TestSecRight);

                grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                grd.EdgeTagNames.Add(5, "navierslip_linear_sensor");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + 0) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8 && ((X[0] >= b / 2) || (X[0] <= -b / 2)))
                        et = 2;
                    if (Math.Abs(X[0] + L) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;
                    if (Math.Abs(X[1] - h) <= 1.0e-8 && (X[0] >= -b / 2) && (X[0] <= b / 2))
                        et = 5;
                    if (Math.Abs(X[0] - b / 2) <= 1.0e-8 && (X[1] >= h))
                        et = 5;
                    if (Math.Abs(X[0] + b / 2) <= 1.0e-8 && (X[1] >= h))
                        et = 5;

                    return et;
                });

                return grd;
            };


            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = X => -1.0;

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("navierslip_linear_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("pressure_outlet_left");
            C.AddBoundaryValue("pressure_outlet_right");

            C.AddBoundaryValue("navierslip_linear_sensor");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.LSContiProjectionMethod = ContinuityProjectionOption.None;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = compMode;
            double dt = 1.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }
        public static XNSE_Control RotCubeDomainDecompoitionError() {
            //gives a domain decomposition error on 43th time step with np=3
            //default parameter set with AMR
            //cs:BoSSS.Application.XNSE_Solver.HardcodedControl.RotCubeDomainDecompoitionError()

            var C = Rotating_Something_Unsteady(4, 30, 2, true);
            return C;
        }

        public static XNSE_Control RotCubeMemoryError() {
            var C = Rotating_Something_Unsteady(4, 20, 3, true);
            return C;
        }

        public static XNSE_Control RotCubePicardConvergenceError() {
            //this error does happen when number of processsors =2 but not with np=4
            var C = Rotating_Something_Unsteady(4, 20, 2, false);
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 0.000001; //it can also happen with 0.001 depending on parameters
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            C.PhysicalParameters.IncludeConvection = true;
            C.NoOfTimesteps = 300;
            return C;
        }

        public static XNSE_Control RotCubeFreeMeanValueError() {
            // cs:BoSSS.Application.XNSE_Solver.HardcodedControl.RotCubeFreeMeanValueError()

            //this error does happen when number of processors =2 but not with np=4
            var C = Rotating_Something_Unsteady(k: 4, Res: 30, SpaceDim: 2, useAMR: false, Gshape: Shape.Cube, OuterBcType: IncompressibleBcType.Wall);
            
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            //C.NonLinearSolver.ConvergenceCriterion = 1.0e-8;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.Globalization = Newton.GlobalizationOption.None;
            C.NonLinearSolver.ConvergenceCriterion = 0.0; // As accurate as possible
            C.NonLinearSolver.MinSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 20;
            
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.UseSchurBlockPrec = false;

            C.PhysicalParameters.IncludeConvection = true;
            C.NoOfTimesteps = 300;
            C.NoOfTimesteps = 1;

         

            
            //C.UseSchurBlockPrec = true;
            
            C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
                        
            C.TracingNamespaces = "BoSSS";

            return C;
        }

        public enum TestCase {
            Bubble,
            Droplet
        }

        public static XNSE_Control MergingDroplet3DLevelSetError() {
            var C = MergingBubble(1, 20, 3, TestCase.Droplet,g: 0);

            // Update relevant parameters"
            // ==============
            string LargeDrop, SmallDrop;
            double[] center_l, center_s;
            double radius_l = 0.25;
            double radius_s = 0.2;


            // Lvl Set 3D
            center_l = new double[] { 0.5, 1.0008, 0.5 };
            center_s = new double[] { 0.5, 0.5402, 0.5 };

            LargeDrop = $"(Math.Sqrt((X[0] - {center_l[0]}).Pow2() + (X[1] - {center_l[1]}).Pow2() + (X[2] - {center_l[2]}).Pow2()) - {radius_l}) ";
            SmallDrop = $"(Math.Sqrt((X[0] - {center_s[0]}).Pow2() + (X[1] - {center_s[1]}).Pow2() + (X[2] - {center_s[2]}).Pow2()) - {radius_s} )";


            string code = $"(X) => 0";
            code = $"(X) => -Math.Min(" + LargeDrop + " , " + SmallDrop + " ) ";
                
            C.InitialValues.Add("VelocityY#B", new Formula( $"X =>  X[1] > {center_l[1]-radius_l}  ?  0.1 : 1.5")  );
            C.SessionName = "Droplet_" + C.SessionName;

            var my_formula = new Formula(code);

            C.InitialValues.Remove("Phi");
            //Phi
            C.InitialValues.Add("Phi",
                        new Formula(code)
                        );
            C.NoOfTimesteps = 100;
            return C;
        }

        //test if source-to-source  agglomeration groups can be formed (test case does not have a physical meaning)
        public static XNSE_Control TwoTorusesAggTestCase(int k = 2, int Res = 20, LevelSetHandling LSMethod = LevelSetHandling.LieSplitting, bool AMR = true, double g = -9.81) {
            XNSE_Control C = new XNSE_Control();
            int SpaceDim = 2;
            //C.DbPath = @"C:\debug_db";
            C.ProjectName = "XNSE-Torus";
            C.ProjectDescription = "agglomeration test";
            C.Tags.Add("level set");
            C.Tags.Add(String.Format("{0}D", SpaceDim));
            C.savetodb = false;

            // Simulation settings
            // ==============
            int NoOfTimeSteps = 10;
            bool Steady = false;
            double bigR = 0.26;
            double smallR = 0.22;
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 0;

            //bool IncludeConvection = true;
            var LeftCenter = new double[] { -0.5, 0 };
            var RightCenter = new double[] { 0.5, 0.0 };

            C.SessionName = $"2DTorus_k{k}_Res{Res}_AMR{AMR}_LS{LSMethod}";
            C.GridFunc = GridFuncFactory(SpaceDim, Res, false, IncompressibleBcType.Pressure_Outlet);
            C.GridPartType = GridPartType.Hilbert;
            C.DynamicLoadBalancing_On = false;

            // Physical Parameters
            // =================== 
            // Bo = 250?, Re = 35?
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 100;
            C.PhysicalParameters.mu_A = 0.01;
            C.PhysicalParameters.mu_B = 0.1;
            C.PhysicalParameters.Sigma = 0.097;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;


            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));

            //level-set for two toruses
            string LeftTorus = $"-Math.Sqrt( (Math.Sqrt((X[0] -({LeftCenter[0]}))*(X[0] -({LeftCenter[0]})) + (X[1]- ({LeftCenter[1]}))*(X[1]- ({LeftCenter[1]}))  ) - {bigR}) * (Math.Sqrt((X[0] -({LeftCenter[0]}))*(X[0] -({LeftCenter[0]})) + (X[1]- ({LeftCenter[1]}))*(X[1]- ({LeftCenter[1]}))  ) - {bigR}) ) + {smallR}";
            string RightTorus = $"-Math.Sqrt( (Math.Sqrt((X[0] -({RightCenter[0]}))*(X[0] -({RightCenter[0]})) + (X[1]- ({RightCenter[1]}))*(X[1]- ({RightCenter[1]}))  ) - {bigR}) * (Math.Sqrt((X[0] -({RightCenter[0]}))*(X[0] -({RightCenter[0]})) + (X[1]- ({RightCenter[1]}))*(X[1]- ({RightCenter[1]}))  ) - {bigR}) ) + {smallR}";

            string code = $"(X) => Math.Max(" + LeftTorus + " , " + RightTorus + "  ) ";

            var my_formula = new Formula(code);

            //Phi
            //.UseImmersedBoundary = true;
            C.AddInitialValue("Phi", new Formula(code)); //


            C.UseImmersedBoundary = false;
            C.PlotAgglomeration = true;
            C.AddInitialValue("Pressure", new Formula(@"X => 0"));
            C.AddBoundaryValue("Pressure_Outlet");

            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.95;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.None;
            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() {
                NoOfMultigridLevels = 5,
                ConvergenceCriterion = 1E-8,
                MaxSolverIterations = 100,
                //MaxKrylovDim = 30,
                TargetBlockSize = 10000,
                //verbose = true
            };


            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1E-3;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = AMR;
            if (AMR) {
                C.SetMaximalRefinementLevel(1);
                C.AMR_startUpSweeps = 1;
            }

            // Timestepping
            // ============
            double dt = -1;
            if (Steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                dt = 1000;
                C.NoOfTimesteps = 1;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                dt = 1;
                C.NoOfTimesteps = NoOfTimeSteps;
            }
            C.TimeSteppingScheme = TimeSteppingScheme.ExplicitEuler; //BD4
            C.dtFixed = dt;
            C.SkipSolveAndEvaluateResidual = true;
            //C.SessionName = "SolverOn_" + C.SessionName;
            return C;
        }


        public static XNSE_Control MergingBubble(int k = 2, int Res = 20, int SpaceDim = 2, TestCase myTestCase = TestCase.Bubble, LevelSetHandling LSMethod = LevelSetHandling.LieSplitting, bool AMR = true, double g = -9.81) {
            XNSE_Control C = new XNSE_Control();
            //C.DbPath = @"C:\debug_db";
            C.ProjectName = "XNSE-Bubble";
            C.ProjectDescription = "merging bubble";
            C.Tags.Add("level set");
            C.Tags.Add(String.Format("{0}D", SpaceDim));
            C.savetodb = false;

            C.SessionName = $"{myTestCase}Merger_k{k}_Res{Res}_AMR{AMR}_LS{LSMethod}";
            C.GridFunc = LongGridFuncFactory(SpaceDim, Res);
            C.GridPartType = GridPartType.Hilbert;
            C.DynamicLoadBalancing_On = false;

            // Physical Parameters
            // =================== 
            // Bo = 250?, Re = 35?
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 100;
            C.PhysicalParameters.mu_A = 0.01;
            C.PhysicalParameters.mu_B = 0.1;
            C.PhysicalParameters.Sigma = 0.097;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            // boundary conditions
            // ===================
            C.AddBoundaryValue("Wall");
            C.AddBoundaryValue("FreeSlip");

            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));

            // Initial Values
            // ==============
            string LargeDrop, SmallDrop;
            double[] center_l, center_s;
            double radius_l = 0.25;
            double radius_s = 0.2;

            if (SpaceDim == 3) {
                // Lvl Set 3D
                center_l = new double[] { 0.5, 1.0002, 0.5 };
                center_s = new double[] { 0.5, 0.5409, 0.5 };

                LargeDrop = $"(Math.Sqrt((X[0] - {center_l[0]}).Pow2() + (X[1] - {center_l[1]}).Pow2() + (X[2] - {center_l[2]}).Pow2()) - {radius_l}) ";
                SmallDrop = $"(Math.Sqrt((X[0] - {center_s[0]}).Pow2() + (X[1] - {center_s[1]}).Pow2() + (X[2] - {center_s[2]}).Pow2()) - {radius_s} )";
                C.SessionName += "_3D";

            } else {
                // Lvl Set 2D
                center_l = new double[] { 0.5, 1.0002 };
                center_s = new double[] { 0.5, 0.5409 };

                LargeDrop = $"(Math.Sqrt((X[0] - {center_l[0]}).Pow2() + (X[1] - {center_l[1]}).Pow2()) - {radius_l}) ";
                SmallDrop = $"(Math.Sqrt((X[0] - {center_s[0]}).Pow2() + (X[1] - {center_s[1]}).Pow2()) - {radius_s} )";
            }

            string code = $"(X) => 0";
            if (myTestCase == TestCase.Droplet) {
                code = $"(X) => -Math.Min(" + LargeDrop + " , " + SmallDrop + " ) ";
                //C.InitialValues.Add("VelocityY#B", new Formula( $"X =>  X[1] > {center_l[1]-radius_l}  ?  -0.1 : -1.5")  );
                C.SessionName = "Droplet_" + C.SessionName;
            } else if (myTestCase == TestCase.Bubble) {
                code = $"(X) => Math.Min(" + LargeDrop + " , " + SmallDrop + " ) ";
                //C.InitialValues.Add("VelocityY#A", new Formula( $"X =>  X[1] > {center_l[1]-radius_l}  ?  -0.1 : -1.5")  );
                C.SessionName = "Bubble_" + C.SessionName;
            }
            var my_formula = new Formula(code);


            //Phi
            C.InitialValues.Add("Phi",
                        new Formula(code)
                        );


            C.UseImmersedBoundary = false;


            double G = g;

            C.InitialValues.Add("GravityY#A",
                        new Formula($"X => {G}")
                        );
            C.InitialValues.Add("GravityY#B",
                    new Formula($"X => {G}")
                    );

            C.SessionName += $"_g{Math.Abs(G):f3}";

            // misc. solver options
            // ====================

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            //C.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            //C.LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;
            //C.SessionName += "_ConstDG";

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.SessionName += "_Picard";

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;


            // Timestepping
            // ============
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = XdgTimestepping.TimeStepperInit.SingleInit;

            string LSname = "";
            switch (LSMethod) {
                case LevelSetHandling.Coupled_Once:
                    LSname = "Moving";
                    break;
                case LevelSetHandling.LieSplitting:
                    LSname = "Lie";
                    break;
                case LevelSetHandling.Coupled_Iterative:
                    LSname = "Iterative";
                    break;
            }

            C.SessionName += "_" + LSname;
            C.Timestepper_LevelSetHandling = LSMethod;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            double dt = 2e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 500;
            C.saveperiod = 10;

            C.AdaptiveMeshRefinement = AMR;
            if (AMR) {
                int AMRlvl = 1;
                C.SetMaximalRefinementLevel(AMRlvl);
                C.AMR_startUpSweeps = AMRlvl;
            }


            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye; //default

            return C;
        }

        public static XNSE_Control RotatingTilted3DTorusAgg0() {
            // throws an error due to the SymPart_DiagBlockEquilib_DropIndefinite
            var C = RotatingTiltedXRigid(2, 16, 3, true, AMRLevel: 1, TiltAngle: Math.PI / 4, SolverOn: false);
            C.AgglomerationThreshold = 0.0;
            C.NoOfTimesteps = 1;
            return C;
        }

        public static XNSE_Control RotatingTilted3DTorus() {
            // throws an error due to the SymPart_DiagBlockEquilib_DropIndefinite
            var C = RotatingTiltedXRigid(1, 20, 3, true, AMRLevel: 2, TiltAngle: Math.PI/6, SolverOn: true);
            C.NoOfTimesteps = 100;
            return C;
        }

        public static XNSE_Control VanishingRotatingCircle() {
            // throws an error due to the SymPart_DiagBlockEquilib_DropIndefinite
            var C = RotatingTiltedXRigid(1, 20, 2, false, shape: Shape.Sphere,  RotAxis: "z", SolverOn: true, rateOfRadius: -4);
            C.PlotAgglomeration = true;
            C.NoOfTimesteps = 100;
            return C;
        }

        public static XNSE_Control RotatingPopcorn() {
            var C = RotatingTiltedXRigid(1, 32, 2, false, shape: Shape.Popcorn, RotAxis: "z", SolverOn: true, rateOfRadius: 0.0, TiltAngle: 0.0, partRad: 0.6);
            C.PlotAgglomeration = true;
            C.NoOfTimesteps = 10;
            var config = new OperatorAnalysisConfig();
            config.CalculateMassMatrix = true;
            C.PostprocessingModules.Add(new Logging.CondLogger(config));
            return C;
        }

        public static XNSE_Control InfiniteConditionNumberTorus() {
            // this test case was resulting in infinite condition numbers
            var C = RotatingTiltedXRigid(2, 128, 2, false, shape: Shape.Torus, RotAxis: "z", SolverOn: false, rateOfRadius: 0.0, TiltAngle: 0.0, partRad: 0.6);
            C.PlotAgglomeration = true;
            C.NoOfTimesteps = 10;
            C.AgglomerationThreshold = 0;
            //C.CalculateConditionNumber = AppControl._ConditionStudy.OperatorandMass;

            return C;
        }

        public static XNSE_Control RotatingTiltedXRigid(int k = 3, int Res = 20, int SpaceDim = 2, bool AMR = true, int AMRLevel = 1, bool LoadBalance = false, Shape shape = Shape.Torus, double TiltAngle = Math.PI/4, string RotAxis = "y", IncompressibleBcType OuterBcType = IncompressibleBcType.Pressure_Outlet, bool SolverOn = true, double rateOfRadius = 0.0, double partRad = 0.39) {
            XNSE_Control C = new XNSE_Control();

            // Simulation Settings
            // ===================
            bool IncludeConvection = true;
            int NoOfTimeSteps = 200;

            C.SessionName = "Solver" + SolverOn + "_" + C.SessionName;
            C.savetodb = false;
            //C.DbPath = @"D:\trash_db";
            C.ProjectName = "XNSE/IBM_test";
            C.ProjectDescription = "rotating" + shape.ToString();
            C.Tags.Add("rotating");
            C.Tags.Add("level set");
            C.Tags.Add(String.Format("{0}D", SpaceDim));

            switch (OuterBcType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Pressure_Outlet:
                    // ok;
                    break;

                default:
                    throw new ArgumentException("not recommended to use boundary condition: " + OuterBcType);

            }

            C.GridFunc = GridFuncFactory(SpaceDim, Res, false, OuterBcType);

            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double Re = 1000;
            double muA = 0.01;

            //double partRad = 0.39;
            double d_hyd = 2 * partRad;
            double VelocityMax = Re * muA / rhoA / d_hyd;
            double anglev = VelocityMax / partRad;

            // The longest arm from the center of rotation is larger than the given parameter partRad for the popcorn case.
            if (shape == Shape.Popcorn)
                anglev = anglev / 1.39;

            double[] pos = new double[SpaceDim];
            double ts = 2 * Math.PI / anglev / NoOfTimeSteps; //   1 revolution around its rot. axis
            Console.WriteLine("VelocityMax: {0} Angular Velocity: {1}", VelocityMax, anglev);

            C.PhysicalParameters.IncludeConvection = IncludeConvection;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;

            string rotAxis;
            if (RotAxis == null)
                rotAxis = "y";
            else
                rotAxis = RotAxis;

            C.SessionName = string.Format("k{0}_Re{1}_t{2}_ti{3:f3}at{4}_r{5}_AMR{6}_LoadBalance{7}", k, Re, NoOfTimeSteps, TiltAngle, rotAxis.ToUpper(), partRad, AMR, LoadBalance);
            if (IncludeConvection) {
                C.SessionName += "_NSE";
                C.Tags.Add("NSE");
            } else {
                C.SessionName += "_Stokes";
                C.Tags.Add("Stokes");
            }
            C.Tags.Add(SpaceDim + "D");
            C.Tags.Add("transient");

            // DG degrees
            // ==========
            C.SetFieldOptions(k, Math.Max(k, 2));
            C.saveperiod = 5;

            C.GridPartType = GridPartType.Hilbert;
            //C.DynamicLoadbalancing_ClassifierType = ClassifierType.CutCells;
            C.DynamicLoadBalancing_On = LoadBalance;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 10;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0.1;

            C.ImmediatePlotPeriod = 10;
            C.SuperSampling = 0;

            //Set xRigid 
            double ringRad = partRad / 2;
            C.Rigidbody.SetParameters(pos, anglev, partRad, SpaceDim, ringRad, rateOfRadius);
            C.Rigidbody.SpecifyShape(shape);
            C.Rigidbody.SetRotationAxis(rotAxis);

            var tiltAxis = new Vector(1, 0, 0);
            C.SessionName += $"_W{anglev:f2}";
            C.Rigidbody.SetTilt(tiltAxis, TiltAngle);
            C.AddInitialValue(VariableNames.LevelSetCGidx(0), new Formula("X => -1"));
            C.UseImmersedBoundary = true;

            C.AddInitialValue("Pressure", new Formula(@"X => 0"));
            C.AddBoundaryValue(OuterBcType.ToString());

            //double inletdelay = 5*ts;
            //C.AddBoundaryValue("Velocity_inlet","VelocityX",new Formula($"(X,t) => {VelocityIn}*(double)(t<={inletdelay}?(t/{inletdelay}):1)",true));
            //C.AddBoundaryValue("Velocity_inlet","VelocityX",new Formula($"(X) => {VelocityIn}"));
            //C.AddInitialValue("VelocityX", new Formula($"(X) => {VelocityIn}"));


            // discretization settings
            C.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.2;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() {
                NoOfMultigridLevels = 5,
                ConvergenceCriterion = 1E-8,
                MaxSolverIterations = 100,
                //MaxKrylovDim = 30,
                TargetBlockSize = 10000,
                //verbose = true
            };
            //C.LinearSolver.NoOfMultigridLevels = 5;
            //C.LinearSolver.ConvergenceCriterion = 1E-6;
            //C.LinearSolver.MaxSolverIterations = 200;
            //C.LinearSolver.MaxKrylovDim = 50;
            //C.LinearSolver.TargetBlockSize = 10000;
            //C.LinearSolver.verbose = true;
            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1E-3;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = AMR;
            if (AMR) {
                C.SetMaximalRefinementLevel(AMRLevel);
                C.AMR_startUpSweeps = AMRLevel;
            }


            // Timestepping
            // ============
            C.NoOfTimesteps = 1; //NoOfTimeSteps
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler; //BD4
            C.dtFixed = ts;
            C.SkipSolveAndEvaluateResidual = !SolverOn;
            Console.WriteLine(C.SessionName); 

            return C;

        }


        //A copy form XNSE_Solver_MPItest.cs
        public static XNSE_Control Rotating_Something_Unsteady(int k = 4, int Res = 30, int SpaceDim = 2, bool useAMR = true, bool useLoadBal = false, Shape Gshape = Shape.Cube, bool UsePredefPartitioning = false, IncompressibleBcType OuterBcType = IncompressibleBcType.Wall) {
            XNSE_Control C = new XNSE_Control();

            switch(OuterBcType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Pressure_Outlet:
                    // ok;
                    break;

                default:
                    throw new ArgumentException("not recommended to use boundary condition: " + OuterBcType);

            }


            // basic database options
            // ======================

            C.savetodb = false;
            //C.DbPath = @"D:\trash_db";
            C.ProjectName = "XNSE/IBM_test";
            C.ProjectDescription = "rotating cube";
            C.Tags.Add("rotating");
            C.Tags.Add("level set");
            C.Tags.Add(String.Format("{0}D", SpaceDim));

            // DG degrees
            // ==========

            C.SetFieldOptions(k, Math.Max(2, k * 2));
            if (UsePredefPartitioning) {
                C.GridPartType = GridPartType.Predefined;
                C.GridPartOptions = "testgrid";
            } else
                C.GridPartType = GridPartType.METIS;

            C.SessionName = "XNSE_rotcube_test";
            C.saveperiod = 1;


            // grid and boundary conditions
            // ============================

            //// Create Grid
            C.GridFunc = GridFuncFactory(SpaceDim, Res, UsePredefPartitioning, OuterBcType);

            // Physical Parameters
            // ===================
            const double rhoA = 1;
            const double muA = 1;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;
            C.PhysicalParameters.rho_A = rhoA;
            C.PhysicalParameters.mu_A = muA;
            double anglev = 10;
            double[] pos = new double[SpaceDim];
            double particleRad = 0.261;

            var PhiFunc = GenPhiFunc(anglev, Gshape, SpaceDim, particleRad, pos);

            Func<double[], double, double[]> VelocityAtIB = delegate (double[] X, double time) {

                if (pos.Length != X.Length)
                    throw new ArgumentException("check dimension of center of mass");

                Vector angVelo = new Vector(new double[] { 0, 0, anglev });
                Vector CenterofMass = new Vector(pos);
                Vector radialVector = new Vector(X) - CenterofMass;
                Vector transVelocity = new Vector(new double[SpaceDim]);
                Vector pointVelocity;

                switch (SpaceDim) {
                    case 2:
                        pointVelocity = new Vector(transVelocity[0] - angVelo[2] * radialVector[1], transVelocity[1] + angVelo[2] * radialVector[0]);
                        break;
                    case 3:
                        pointVelocity = transVelocity + angVelo.CrossProduct(radialVector);
                        break;
                    default:
                        throw new NotImplementedException("this number of dimensions is not supported");
                }

                return pointVelocity;
            };

            Func<double[], double, double> VelocityX = delegate (double[] X, double time) { return VelocityAtIB(X, time)[0]; };
            Func<double[], double, double> VelocityY = delegate (double[] X, double time) { return VelocityAtIB(X, time)[1]; };
            Func<double[], double, double> VelocityZ = delegate (double[] X, double time) { return VelocityAtIB(X, time)[2]; };

            var PhiFuncDelegate = BoSSS.Solution.Utils.NonVectorizedScalarFunction.Vectorize(PhiFunc);

            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(0), X => -1);
            C.UseImmersedBoundary = true;
            if (C.UseImmersedBoundary) {
                //C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
                C.InitialValues_EvaluatorsVec.Add(VariableNames.LevelSetCGidx(1), PhiFuncDelegate);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi2", VelocityX);
                C.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi2", VelocityY);
                if (SpaceDim == 3)
                    C.InitialValues_Evaluators_TimeDep.Add("VelocityZ@Phi2", VelocityZ);
            }
            C.InitialValues_Evaluators.Add("Pressure", X => 0);
            
            // misc. solver options
            // ====================

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.UseSchurBlockPrec = true;
            C.AgglomerationThreshold = 0.1;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.Option_LevelSetEvolution2 = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.verbose = true;

            C.AdaptiveMeshRefinement = useAMR;
            if (useAMR) {
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
                C.AMR_startUpSweeps = 1;
            }

            C.DynamicLoadBalancing_On = useLoadBal;
            C.DynamicLoadBalancing_RedistributeAtStartup = true;
            C.DynamicLoadBalancing_Period = 1;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0;

            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            double dt = 0.005;
            C.dtMax = dt;
            C.dtMin = dt;
            C.dtFixed = dt;
            C.NoOfTimesteps = 50;

            return C;
        }

        private static Func<double[], double, double> GenPhiFunc(double anglev, Shape Gshape, int SpaceDim, double particleRad, double[] pos) {
            return delegate (double[] X, double t) {
                double angle = -(anglev * t) % (2 * Math.PI);
                switch (Gshape) {
                    case Shape.Cube:
                        switch (SpaceDim) {
                            case 2:
                                return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                    Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                                    + particleRad;
                            case 3:
                                return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                                        Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)), Math.Abs(X[2] - pos[2])))
                                                        + particleRad;
                            default:
                                throw new NotImplementedException("Dimension not supported");
                        };
                    case Shape.Sphere:
                        switch (SpaceDim) {
                            case 2:
                                return -X[0] * X[0] - X[1] * X[1] + particleRad * particleRad;
                            case 3:
                                return -X[0] * X[0] - X[1] * X[1] - X[2] * X[2] + particleRad * particleRad;
                            default:
                                throw new NotImplementedException("Dimension not supported");
                        }
                    //case Geometry.Parted:
                    //    return -X[0] + 0.7;
                    default:
                        throw new NotImplementedException("Shape unknown");
                }
            };
        }

        public static Func<IGrid> LongGridFuncFactory(int SpaceDim, int Res) {
            double xMax = 1, yMax = 2, zMax = 1;
            double xMin = 0, yMin = 0, zMin = 0;

            if (SpaceDim == 2) {
                return delegate {
                    var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                    var _yNodes = GenericBlas.Linspace(yMin, yMax, Res * 2 + 1);

                    GridCommons grd;
                    grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);


                    grd.EdgeTagNames.Add(1, "Wall");
                    //grd.EdgeTagNames.Add(2, "wall_upper");
                    grd.EdgeTagNames.Add(3, "FreeSlip");
                    //grd.EdgeTagNames.Add(4, "freeslip_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (X[1] <= yMin)
                            et = 1;
                        if (X[1] >= yMax)
                            et = 1;
                        if (X[0] <= xMin)
                            et = 3;
                        if (X[0] >= xMax)
                            et = 3;

                        return et;
                    });
                    return grd;
                };
            } else if (SpaceDim == 3) {
                return delegate {

                    var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                    var _yNodes = GenericBlas.Linspace(yMin, yMax, Res * 2 + 1);
                    var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);
                    GridCommons grd;
                    grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes);
                    grd.EdgeTagNames.Add(1, "Wall");
                    //grd.EdgeTagNames.Add(2, "wall_upper");
                    grd.EdgeTagNames.Add(3, "FreeSlip");
                    //grd.EdgeTagNames.Add(4, "freeslip_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (X[1] <= yMin)
                            et = 1;
                        if (X[1] >= yMax)
                            et = 1;
                        if (X[0] <= xMin)
                            et = 3;
                        if (X[0] >= xMax)
                            et = 3;
                        if (X[2] <= zMin)
                            et = 3;
                        if (X[2] >= zMax)
                            et = 3;
                        return et;
                    });
                    return grd;
                };

            } else {
                throw new NotImplementedException();
            }

        }


        public static Func<IGrid> GridFuncFactory(int SpaceDim, int Res, bool UsePredefPartitioning, IncompressibleBcType OuterBcType) {
            double xMin = -1, yMin = -1, zMin = -1;
            double xMax = 1, yMax = 1, zMax = 1;

            // Predefined Partitioning
            Func<double[], int> MakeDebugPart = delegate (double[] X) {
                double x = X[0];
                double range = xMax - xMin;
                double interval = range * 0.9 / (ilPSP.Environment.MPIEnv.MPI_Size - 1); // last rank gets only a 10% stripe of the domain
                return (int)((x - xMin) / interval);
            };

            // The Grid Function
            return delegate {

                var _xNodes = GenericBlas.Linspace(xMin, xMax, Res + 1);
                var _yNodes = GenericBlas.Linspace(yMin, yMax, Res + 1);
                var _zNodes = GenericBlas.Linspace(zMin, zMax, Res + 1);

                GridCommons grd;
                switch (SpaceDim) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                        break;
                    case 3:
                        grd = Grid3D.Cartesian3DGrid(_xNodes, _yNodes, _zNodes, Foundation.Grid.RefElements.CellType.Cube_Linear, false, false, false);
                        break;
                    default:
                        throw new ArgumentOutOfRangeException();
                }
                if (UsePredefPartitioning) grd.AddPredefinedPartitioning("testgrid", MakeDebugPart);
                grd.EdgeTagNames.Add(2, OuterBcType.ToString());
                grd.DefineEdgeTags(delegate (double[] _X) {
                    return 2;
                });
                return grd;
            };
        }
    }
}
