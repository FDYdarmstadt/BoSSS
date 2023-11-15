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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using ilPSP;
using BoSSS.Solution.AdvancedSolvers.Testing;
using ilPSP.Connectors.Matlab;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Solution.Tecplot;

namespace BoSSS.Application.XNSE_Solver.Tests {


    public interface IXNSElsTest : IXNSETest {

        /// <summary>
        /// compute a suitable timestep for various combinations of grid resolutions and level-set degrees
        /// </summary>
        double ComputeTimestep(int gridResolution, int lsDegree, int AMRlevel);

        /// <summary>
        /// 
        /// </summary>
        double getEndTime();

    }

    /// <summary>
    /// A collection of all-up NUnit tests for the XNSE solver regarding various level set handling
    /// </summary>
    [TestFixture]
    static public partial class LevelSetUnitTests {

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        [Test]
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
        [Test]
        public static void LevelSetAdvectionTest2D_reverse(
           [Values(2, 3, 4)] int LSdegree,
           [Values(0, 1, 2)] int AMRlevel,
           [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
           [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {
            LevelSetAdvectionTest2D(LSdegree, AMRlevel, levelSetEvolution, levelSetHandling, true);
        }



        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>
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

            var Tst = new LevelSetAdvectionTest(2, LSdegree, reversed);
            var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling, gridResolution, AMRlevel);
            C.SkipSolveAndEvaluateResidual = true;

            //string IO = $"LSAdvectionTest2D-deg{LSdegree}-amrLvl{AMRlevel}-lsEvo{levelSetEvolution}-rev{reversed}-grdRes{gridResolution}";
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 3;

            LevelSetTest(Tst, C);

            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!

        }

        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetAdvectionTest3D(
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(0, 1, 2)] int AMRlevel,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling, 
        //    [Values(false, true)] bool reversed)
        //    {

        //    var Tst = new LevelSetAdvectionTest(3, LSdegree, reversed);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 20, levelSetEvolution, levelSetHandling, 1);
        //    C.SkipSolveAndEvaluateResidual = true;

        //    LevelSetTest(Tst, C);

        //}


        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetAdvectionOnWallTest2D(
        //    [Values(Math.PI/4, Math.PI/3, Math.PI/6)] double contactAngle,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(0, 1, 2)] int AMRlevel,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {

        //    var Tst = new LevelSetAdvectionOnWallTest(2, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling, 1, AMRlevel);
        //    C.SkipSolveAndEvaluateResidual = true;

        //    LevelSetTest(Tst, C);

        //}

        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetAdvectionOnWallTest3D(
        //    [Values(Math.PI / 4, Math.PI / 3, Math.PI / 6)] double contactAngle,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(0, 1, 2)] int AMRlevel,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling) {

        //    var Tst = new LevelSetAdvectionOnWallTest(3, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling, 1, AMRlevel);
        //    C.SkipSolveAndEvaluateResidual = true;

        //    LevelSetTest(Tst, C);

        //}

        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetScalingTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetScalingTest(
        //    [Values(2)] int spatialDimension,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling)
        //    //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timeSteppingScheme) 
        //    {
        //    // Todo: singleInit/multiInit, 

        //    var Tst = new LevelSetScalingTest(spatialDimension, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling);

        //    LevelSetTest(Tst, C);

        //}


        ///// <summary>
        ///// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetRotationTest"/>
        ///// </summary>
        //[Test]
        //public static void LevelSetRotationTest(
        //    [Values(2)] int spatialDimension,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling)
        //    //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timeSteppingScheme) 
        //    {
        //    // Todo: singleInit/multiInit, 

        //    var Tst = new LevelSetRotationTest(spatialDimension, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling);

        //    LevelSetTest(Tst, C);

        //}




        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>



        public static void LevelSetTest(IXNSETest Tst, XNSE_Control C, string IO = null) {

            using (var solver = new XNSE()) {

                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                AgglomerationAlgorithm.Katastrophenplot = delegate (DGField[] plotFields,string Tag) {
                    Tag ??= "";
                    List<DGField> allfields = new();
                    allfields.AddRange(plotFields);
                    
                    foreach(var f in solver.RegisteredFields) {
                        if(!allfields.Contains(f, (a, b) => object.ReferenceEquals(a, b)))
                            allfields.Add(f);
                    }

                    Tecplot.PlotFields(allfields, Tag + "AgglomFail", 0.0, 4);
                };

                solver.Init(C);
                solver.RunSolverMode();
                //solver.OperatorAnalysis();

                //-------------------Evaluate Error ---------------------------------------- 
                var evaluator = new LevelSetErrorEvaluator<XNSE_Control>(solver);
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


        class ASLS_XNSE_Control : XNSE_Control {
            public override Type GetSolverType() {
                return typeof(XNSE);
            }
        }


        static ASLS_XNSE_Control LSTstObj2CtrlObj(IXNSElsTest tst, int FlowSolverDegree, int maxNoTimesteps,
            LevelSetEvolution levelSetEvolution, LevelSetHandling levelSetHandling, 
            int GridResolution = 1, int AMRlevel = 0) {

            ASLS_XNSE_Control C = new ASLS_XNSE_Control();
            int D = tst.SpatialDimension;


            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "LevelSet/" + tst.GetType().Name;
            C.ProjectDescription = "Test";


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
