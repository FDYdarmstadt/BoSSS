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

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools;


namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// Unit tests for the XNSE-solver checking various properties for the static droplet test case
    /// - exact pressure jump
    /// - spurious velocities
    /// - integral properties (gauss and stokes)
    /// </summary>
    [TestFixture]
    public class StaticDropletElementalTest {

        /// <summary>
        /// 2D static droplet test case
        /// </summary>
        /// <returns></returns>
        static public XNSE_Control StaticDroplet2D_ExactEvalution(bool quadraticLS = true, int BCsetup = 1, bool useAMR = false) { 
        
            var C = new XNSE_Control();

            int k = 1;
            C.SetDGdegree(k);

            // Physical Parameters
            // ===================
            #region physics params

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            double sigma = 1.23;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double Lscl = 1.0;
            double xSize = Lscl * 1.0;
            double ySize = Lscl * 1.0;

            int kelem = 5;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize / 2.0, xSize / 2.0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-ySize / 2.0, ySize / 2.0, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                switch (BCsetup) {
                    case 0: {
                            grd.EdgeTagNames.Add(1, "wall_bottom");
                            grd.EdgeTagNames.Add(2, "wall_top");
                            grd.EdgeTagNames.Add(3, "wall_left");
                            grd.EdgeTagNames.Add(4, "wall_right");
                            break;
                        }
                    case 1: {
                            grd.EdgeTagNames.Add(1, "pressure_outlet_bottom");
                            grd.EdgeTagNames.Add(2, "pressure_outlet_top");
                            grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                            grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                            break;
                        }
                    default: {
                            goto case 0;
                        }
                }


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
            #region boundary conditions

            switch (BCsetup) { 
                case 0: {
                        C.AddBoundaryValue("wall_bottom");
                        C.AddBoundaryValue("wall_top");
                        C.AddBoundaryValue("wall_left");
                        C.AddBoundaryValue("wall_right");
                        break;
                    }
                case 1: {
                        C.AddBoundaryValue("pressure_outlet_bottom");
                        C.AddBoundaryValue("pressure_outlet_top");
                        C.AddBoundaryValue("pressure_outlet_left");
                        C.AddBoundaryValue("pressure_outlet_right");
                        break;
                    }
                default: { 
                    goto case 0;
                    }
            }

            #endregion


            // initial values & level set
            // ==============
            #region initial values

            double[] center = new double[] { 0.0, 0.0 };
            double R = 0.2312;
            if (quadraticLS) {
                C.InitialValues_Evaluators.Add("Phi", (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2()));         // quadratic
                C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;
            } else {
                C.InitialValues_Evaluators.Add("Phi", (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R));         // signed distance
                C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            }

            //Console.WriteLine($"pressure jump: {sigma / R}");
            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / R);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Steady;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;


            #endregion


            // solver
            // ======

            C.SkipSolveAndEvaluateResidual = true;

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            return C;
        }


        /// <summary>
        /// 3D static droplet test case
        /// </summary>
        /// <returns></returns>
        static public XNSE_Control StaticDroplet3D_ExactEvaluation() {

            var C = new XNSE_Control();



            return C;
        }


        [Test]
        public static void TestStaticDroplet2D_ExactSolution() {

            var ctrl = StaticDroplet2D_ExactEvalution();

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 3;

            using (var TestRun = new XNSE()) {

                TestRun.Init(ctrl);
                TestRun.RunSolverMode();
            }

        }

    }
}
