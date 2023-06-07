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

using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// An all-up NUnit test for the XNSEC application.
    /// </summary>
    [TestFixture]
    static public partial class NUnitTest {

       

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        //[Test]
        public static void LevelSetAdvectionTest2D_fwd(
               [Values(2, 3, 4)] int LSdegree,
               [Values(0, 1, 2)] int AMRlevel,
               [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
               [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {
            LevelSetAdvectionTest2D(LSdegree, AMRlevel, levelSetEvolution, levelSetHandling, false);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        //[Test]
        public static void LevelSetAdvectionTest2D_reverse(
           [Values(2, 3, 4)] int LSdegree,
           [Values(0, 1, 2)] int AMRlevel,
           [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
           [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {
            LevelSetAdvectionTest2D(LSdegree, AMRlevel, levelSetEvolution, levelSetHandling, true);
        }

        public static void LevelSetAdvectionTest2D(

    [Values(2, 3, 4)] int LSdegree,
    [Values(0, 1, 2)] int AMRlevel,
    [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
    [Values(false, true)] bool reversed) {
            int gridResolution;
            switch (LSdegree) {
                case 2: gridResolution = 3; break;
                case 3: gridResolution = 2; break;
                //case 4: gridResolution = 1; break;
                default:
                    gridResolution = 1; break;
            }

            if (LSdegree == 4 && AMRlevel > 0 && levelSetEvolution == LevelSetEvolution.StokesExtension)
                return;

            var Tst = new LevelSetAdvectionTestXNSEC(2, LSdegree, reversed);
            var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling, gridResolution, AMRlevel);
            C.SkipSolveAndEvaluateResidual = false;
            //C.ImmediatePlotPeriod = 1;
            C.NonLinearSolver.verbose = true;
            //string IO = $"LSAdvectionTest2D-deg{LSdegree}-amrLvl{AMRlevel}-lsEvo{levelSetEvolution}-rev{reversed}-grdRes{gridResolution}";
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 3;
            LevelSetTest(Tst, C);

            //Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!
        }

        public class LevelSetAdvectionTestXNSEC : LevelSetAdvectionTest, IXNSEClsTest {

            public LevelSetAdvectionTestXNSEC(int spatDim, int LevelSetDegree, bool reversed) : base(spatDim, LevelSetDegree, reversed) {
            }
            public bool EnableMassFractions => false;
            public bool EnableTemperature => false;
            public int NumberOfChemicalComponents => 1;

            public bool ChemicalReactionTermsActive => false;

            public double[] GravityDirection => new double[] { 0, 0, 0 };

            public Func<double[], double, double> GetMassFractions(string species, int q) {
                switch (base.SpatialDimension) {
                    case 2:
                        if (q == 0) {
                            return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetTemperature(string species) {
                switch (base.SpatialDimension) {
                    case 2:
                        return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();

                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }
        }

        public interface IXNSEClsTest : IXNSElsTest, IXNSECTest {
        }

        public static void LevelSetTest(IXNSETest Tst, XNSEC_Control C, string IO = null) {
            using (var solver = new XNSEC()) {
                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();
                //solver.OperatorAnalysis();

                //-------------------Evaluate Error ----------------------------------------
                var evaluator = new LevelSetErrorEvaluator<XNSEC_Control>(solver);
                double[] LastErrors = evaluator.ComputeL2Error(C.Endtime, C);

                //string[] ErrNames = new string[] { "Phi", "PhiDG", "Gradient PhiDG" };
                string[] ErrNames = new string[] { "Phi", "PhiDG" }; //, "Gradient PhiDG", "Interface size", "Interface Area #A", "Interface Area #B" };
                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length) {
                    Console.WriteLine("LastErrors.Length = {0} not equal to ErrThresh.Length = {1}", LastErrors.Length, ErrThresh.Length);
                    throw new ApplicationException();
                }
                for (int i = 0; i < ErrThresh.Length; i++) {
                    Console.WriteLine("L2 error, '{0}': \t{1}", ErrNames[i], LastErrors[i]);
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                //-------------------check parallel simulations ----------------------------------------
                int RefMPIsize = 1;
                if (IO != null) {
                    LevelSet PhiDG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet;
                    LevelSet PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet;

                    var projCheck = new TestingIO(solver.GridData, $"{IO}.csv", true, RefMPIsize);
                    projCheck.AddDGField(PhiDG);
                    projCheck.AddDGField(PhiCG);
                    projCheck.DoIOnow();

                    Assert.Less(projCheck.AbsError(PhiDG), 1.0e-15, "Mismatch in projected PhiDG between single-core and parallel run.");
                    Assert.Less(projCheck.AbsError(PhiCG), 1.0e-15, "Mismatch in projected PhiCG between single-core and parallel run.");
                }
            }
        }

        private static XNSEC_Control LSTstObj2CtrlObj(IXNSEClsTest tst, int FlowSolverDegree, int maxNoTimesteps,
            LevelSetEvolution levelSetEvolution, LevelSetHandling levelSetHandling,
            int GridResolution = 1, int AMRlevel = 0) {
            XNSEC_Control C = new XNSEC_Control();
            int D = tst.SpatialDimension;

            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "LevelSet/" + tst.GetType().Name;
            C.ProjectDescription = "Test";
            C.NumberOfChemicalSpecies = 1;
            C.EnableMassFractions = false;
            C.EnableTemperature = false;
            C.rhoOne = true;
            C.NumberOfChemicalSpecies = tst.NumberOfChemicalComponents;
            C.ChemicalReactionActive = tst.ChemicalReactionTermsActive;
            // DG degree
            // =========

            C.SetDGdegree(FlowSolverDegree);

            // grid
            // ====

            C.GridFunc = () => tst.CreateGrid(GridResolution);

            // boundary conditions
            // ===================

            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // Physical parameters
            // ====================

            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;
            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.Material = tst.Material;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;

            // initial values and exact solution
            // =================================

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));

                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));

                for (int d = 0; d < D; d++) {
                    Func<double[], double, double> Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            // advanced spatial discretization settings
            // ========================================

            //C.CutCellQuadratureType = CutCellQuadratureType;

            // timestepping and solver
            // =======================

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.Option_LevelSetEvolution = levelSetEvolution;
            //if (levelSetEvolution == LevelSetEvolution.Fourier) {   // ToDo: refactor Fourier initialization!
            //    int numSp = 640;
            //    double[] FourierP = new double[numSp];
            //    double[] samplP = new double[numSp];
            //    for (int sp = 0; sp < numSp; sp++) {
            //        FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
            //        samplP[sp] = ((LevelSetAdvectiontTest)tst).Radius;
            //    }

            //    C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 0.1) {
            //        center = new double[] { 0.0, 0.0},
            //        FourierEvolve = Fourier_Evolution.MaterialPoints,
            //        centerMove = CenterMovement.Reconstructed,
            //    };
            //}

            C.Timestepper_LevelSetHandling = levelSetHandling;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            // adaptive mesh refinement
            if (AMRlevel > 0) {
                C.AdaptiveMeshRefinement = true;
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = AMRlevel });
                C.AMR_startUpSweeps = AMRlevel;
            }

            double dt = tst.ComputeTimestep(GridResolution, FlowSolverDegree, AMRlevel);
            C.dtFixed = dt;
            double T = tst.getEndTime();
            int timesteps = (int)(T / dt);
            C.NoOfTimesteps = (timesteps > maxNoTimesteps) ? maxNoTimesteps : timesteps;
            C.Endtime = dt * C.NoOfTimesteps;

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            // return
            // ======

            return C;
        }

     
    }
}