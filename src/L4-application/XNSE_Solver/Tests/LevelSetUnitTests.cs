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

namespace BoSSS.Application.XNSE_Solver.Tests {


    interface IXNSElsTest : IXNSETest {

        /// <summary>
        /// compute a suitable timestep for various combinationa of grid resolutions and level-set degrees
        /// </summary>
        /// <returns></returns>
        double ComputeTimestep(int gridResolution, int lsDegree);

    }

    /// <summary>
    /// A collection of all-up NUnit tests for the XNSE solver regarding various level set handling
    /// </summary>
    [TestFixture]
    static public partial class LevelSetUnitTests {


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetScalingTest"/>
        /// </summary>
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


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetRotationTest"/>
        /// </summary>
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
        [Test]
        public static void LevelSetAdvectionTest(
            [Values(2)] int spatialDimension,
            [Values(2,3,4)] int LSdegree,
            [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling)
            //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timeSteppingScheme) 
            {
                // Todo: singleInit/multiInit, 

            var Tst = new LevelSetAdvectionTest(spatialDimension, LSdegree);
            var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling);

            LevelSetTest(Tst, C);

        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>
        //[Test]
        //public static void LevelSetShearingTest(
        //    [Values(2)] int spatialDimension,
        //    [Values(2, 3, 4)] int LSdegree,
        //    [Values(LevelSetEvolution.FastMarching)] LevelSetEvolution levelSetEvolution,
        //    [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling)
        //    //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timeSteppingScheme) 
        //    {
        //    // Todo: singleInit/multiInit, 

        //    var Tst = new LevelSetShearingTest(spatialDimension, LSdegree);
        //    var C = LSTstObj2CtrlObj(Tst, LSdegree, 40, levelSetEvolution, levelSetHandling);

        //    LevelSetTest(Tst, C);

        //}


        private static void LevelSetTest(IXNSETest Tst, XNSE_Control C) {

            using (var solver = new XNSE()) {

                Console.WriteLine("Warning! - enabled immediate plotting");
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();
                //solver.OperatorAnalysis();

                //-------------------Evaluate Error ---------------------------------------- 
                LevelSetErrorEvaluator evaluator = new LevelSetErrorEvaluator(solver);
                double[] LastErrors = evaluator.ComputeL2Error(C.Endtime, C);

                string[] ErrNames = new string[] { "Phi", "PhiDG", "Gradient PhiDG" };
                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length; i++) {
                    Console.WriteLine("L2 error, '{0}': \t{1}", ErrNames[i], LastErrors[i]);
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

            }
        }


        class ASLS_XNSE_Control : XNSE_Control {
            public override Type GetSolverType() {
                return typeof(XNSE);
            }
        }


        static ASLS_XNSE_Control LSTstObj2CtrlObj(IXNSElsTest tst, int FlowSolverDegree, int NoTimesteps,
            LevelSetEvolution levelSetEvolution, LevelSetHandling levelSetHandling,
            int GridResolution = 2) {

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


            C.NoOfTimesteps = NoTimesteps;
            double dt = tst.ComputeTimestep(GridResolution, FlowSolverDegree);
            C.dtFixed = dt;
            C.Endtime = (double) (NoTimesteps * dt);


            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver.ConvergenceCriterion = 1e-9;

            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            // return
            // ======

            return C;
        }


    }


}
