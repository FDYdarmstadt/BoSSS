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

using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.XDG;

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
        static public XNSE_Control StaticDroplet2D_ExactEvalution(int k, bool quadraticLS = true, int BCsetup = 0, bool useAMR = false) { 
        
            var C = new XNSE_Control();

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

            C.CutCellQuadratureType = CutCellQuadratureMethod.OneStepGaussAndStokes;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.LevelSetTools.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


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


        static void CheckIntegralProperties(XNSE solver) {

            var LsTrk = solver.LsTrk;
            var CCmask = LsTrk.Regions.GetCutCellMask();

            Basis testBasis = new Basis(LsTrk.GridDat, 0);

            SinglePhaseField testFieldX = new SinglePhaseField(testBasis);
            Func<double[], double> fx = (X => 1.0);
            testFieldX.ProjectField(fx);
            SinglePhaseField testFieldY = new SinglePhaseField(testBasis);
            Func<double[], double> fy = (X => 1.0);
            testFieldY.ProjectField(fy);
            SinglePhaseField testFieldZ = new SinglePhaseField(testBasis);
            Func<double[], double> fz = (X => 1.0);
            testFieldZ.ProjectField(fz);

            VectorField<SinglePhaseField> testField = (LsTrk.GridDat.SpatialDimension == 2) ?
                testField = new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY }) :
                testField = new VectorField<SinglePhaseField>(new SinglePhaseField[] { testFieldX, testFieldY, testFieldZ });

            int order = 2 * ((LevelSet)LsTrk.LevelSets[0]).Basis.Degree + 1;


            // check gauss
            // ===========
            double totalGaussError = 0.0;
            Dictionary<int, (double error, double gaussInVolume, double gaussAtInterface, double gaussAtEdges)> gaussErrorPerCell = new Dictionary<int, (double error, double gaussInVolume, double gaussAtInterface, double gaussAtEdges)>();
            SinglePhaseField gaussErrorField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "gaussError");

            foreach (int jC in CCmask.ItemEnum) {
                var result = XNSEUtils.CheckGaussInCutCell(jC, LsTrk, LsTrk.GetSpeciesId("A"), testField, order);
                Console.WriteLine($"Gauss error in Cell {jC}: {result.error}");
                totalGaussError += result.error.Abs();
                gaussErrorPerCell.Add(jC, result);
                gaussErrorField.SetMeanValue(jC, result.error);
            }

            Console.WriteLine($"total Gauss error = {totalGaussError}");
            //Assert.LessOrEqual(totalGaussError, 0.0, $"Gauss error above threshold.");

            // check stokes
            // ============
            double totalStokesError = 0.0;
            Dictionary<int, (double error, double stokesAtInterface, double stokesAtEdges)> stokesErrorPerCell = new Dictionary<int, (double error, double stokesAtInterface, double stokesAtEdges)>();
            SinglePhaseField stokesErrorField = new SinglePhaseField(new Basis(LsTrk.GridDat, 0), "stokesError");

            foreach (int jC in CCmask.ItemEnum) {
                var result = XNSEUtils.CheckStokesForCell(jC, LsTrk, testField, order);
                Console.WriteLine($"Stokes error in Cell {jC}: {result.error}");
                stokesErrorPerCell.Add(jC, result);
                totalStokesError += result.error.Abs();
                stokesErrorField.SetMeanValue(jC, result.error);
            }
            
            Console.WriteLine($"total Stokes error = {totalStokesError}");
            //Assert.LessOrEqual(totalStokesError, 0.0, $"Stokes error above threshold.");

        }


        [Test]
        public static void TestStaticDroplet2D_ExactSolution() {

            var ctrl = StaticDroplet2D_ExactEvalution(1);

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 5;

            using (var TestSolver = new XNSE()) {

                TestSolver.Init(ctrl);
                TestSolver.RunSolverMode();

                CheckIntegralProperties(TestSolver);

            }

        }

    }
}
