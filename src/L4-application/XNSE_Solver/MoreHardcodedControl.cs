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
using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// A few example configurations.
    /// </summary>
    public static class MoreHardcodedControl {

        public static XNSE_Control MeniskusBenchmark(double mu = 0.01, double rho = 1.0, double g = 10.0, double sigma = 0.01, double slip = 0.001, double theta = Math.PI / 2.0 * 30.0 / 90.0, int level = 2, bool boundary = false) {

            var ctrl = new XNSE_Control();

            ctrl.ProjectName = "MeniskusBenchmark";
            ctrl.SessionName = "MeniskusBenchmark" + (boundary ? "Open" : "Closed") + "_AMR" + level + "Rho" + rho;

            ctrl.DbPath = null;
            ctrl.savetodb = ctrl.DbPath != null;
            ctrl.saveperiod = 10;

            if(ctrl.DbPath == null) {
                ctrl.ImmediatePlotPeriod = 1;
                ctrl.SuperSampling = 2;
            }

            int p = 2;
            ctrl.SetDGdegree(p);

            //grid
            int kelemR = 4;

            string[] Bndy;
            if (boundary) {
                Bndy = new string[] {  "Inner",
                                        "NavierSlip_linear_right",
                                        "pressure_outlet_top",
                                        "Slipsymmetry_left",
                                        "pressure_outlet_bottom"};
            } else {
                Bndy = new string[] {  "Inner",
                                        "NavierSlip_linear_right",
                                        "wall_top",
                                        "Slipsymmetry_left",
                                        "wall_bottom"};
            }

            ctrl.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 10, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 10, kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                for (byte i = 1; i < Bndy.Count(); i++) {
                    grd.EdgeTagNames.Add(i, Bndy[i]);
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Xnodes.Last()) < 1e-8)
                        return 1;
                    if (Math.Abs(X[0] - Xnodes.First()) < 1e-8)
                        return 3;
                    if (Math.Abs(X[1] - Ynodes.Last()) < 1e-8)
                        return 2;
                    if (Math.Abs(X[1] - Ynodes.First()) < 1e-8)
                        return 4;
                    return et;
                });

                return grd;
            };


            //material
            double a = Math.Sqrt(sigma / (rho * g));
            double P = mu * Math.Sqrt(9 * Math.Cos(Math.PI / 2.0 * 30.0 / 90.0) * sigma / (Math.Pow(rho, 3) * Math.Pow(g, 2) * Math.Pow(a, 5)));
            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("AMR", level));
            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("rho", rho));
            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("g", g));
            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("sigma", sigma));
            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("P", P));
            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Ls", slip));

            double g_scaled = g / a;
            PhysicalParameters physParams = new PhysicalParameters() {
                rho_A = rho * Math.Pow(a, 3),
                rho_B = rho * Math.Pow(a, 3) / 1000.0,

                mu_A = mu * Math.Pow(a, 1),
                mu_B = mu * Math.Pow(a, 1) / 1000.0,

                theta_e = theta,
                Sigma = sigma,
                sliplength = slip,

                IncludeConvection = false,
                Material = true
            };
            ctrl.PhysicalParameters = physParams;
            ctrl.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;

            // inital values
            double y0 = 5.0 - 1e-3; // small offset to avoid coinciding with edge
            ctrl.AddInitialValue("Phi", $"(X, t) => -{y0} + X[1]", true);
            ctrl.AddInitialValue("GravityY#A", $"(X, t) => -{g_scaled}", true);
            ctrl.AddInitialValue("GravityY#B", $"(X, t) => -{g_scaled}", true);

            //Boundary Conditions
            if (boundary) {
                ctrl.AddBoundaryValue(Bndy[1]);
                ctrl.AddBoundaryValue(Bndy[3]);
                ctrl.AddBoundaryValue(Bndy[2], "Pressure#B", $"(X, t) => 0.0", true);
                ctrl.AddBoundaryValue(Bndy[4], "Pressure#A", $"(X, t) => {y0} * {ctrl.PhysicalParameters.rho_A} * {g_scaled} + (10 - {y0}) * {ctrl.PhysicalParameters.rho_B} * {g_scaled}", true);
            } else {
                ctrl.AddBoundaryValue(Bndy[1]);
                ctrl.AddBoundaryValue(Bndy[3]);
                ctrl.AddBoundaryValue(Bndy[2]);
                ctrl.AddBoundaryValue(Bndy[4]);
            }

            // AMR
            ctrl.AdaptiveMeshRefinement = level > 0;
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRonNarrowband() { levelSet = 0, maxRefinementLevel = level });
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.AMRLevelIndicatorLibrary.AMRonBoundary(1) { maxRefinementLevel = level });
            ctrl.AMR_startUpSweeps = level;

            //Timestepping
            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting;
            ctrl.FastMarchingReInitPeriod = 150;

            ctrl.SkipSolveAndEvaluateResidual = false;

            ctrl.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 0.01;
            ctrl.Endtime = 10.0;
            ctrl.NoOfTimesteps = 1000; //int.MaxValue; // timesteps can be adapted, simulate until endtime is reached

            // Insitu Postprocessing
            ctrl.PostprocessingModules.Add(new MovingContactLineLogging() { LogPeriod = 1 });


            return ctrl;
        }

    }
}
