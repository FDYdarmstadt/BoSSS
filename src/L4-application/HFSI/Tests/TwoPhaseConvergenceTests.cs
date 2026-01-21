using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSE_Solver.Tests;
using HFSISolver.Tests;
using BoSSS.Solution.Control;
using BoSSS.Solution;
using BoSSS.Solution.Statistic;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HFSISolver.Tests {

    /// <summary>
    /// 
    /// </summary>
    [TestFixture]
    static public class TwoPhaseConvergenceTests {


        [Test]
        public static void FluidSolidTaylorCouette_Experimental_P1() {
            Experimental(1);
        }
        [Test]
        public static void FluidSolidTaylorCouette_Experimental_P2() {
            Experimental(2);
        }
        [Test]
        public static void FluidSolidTaylorCouette_Exact_P1() {
            Exact(1);
        }
        [Test]
        public static void FluidSolidTaylorCouette_Exact_P2() {
            Exact(2);
        }
        [Test]
        public static void FluidSolidTaylorCouette_Exact_P3() {
            Exact(3);
        }
        [Test]
        public static void FluidSolidTaylorCouette_Exact_P4() {
            Exact(4);
        }

        public static void Experimental(int p = 2) {
            //double dt = 1.0e200;
            // --test=HFSISolver.Tests.SolidOnlyTests.RotationConvergenceTest

            BoSSS.Solution.Application.DeleteOldPlotFiles();

            var Tst = new FSTaylorCouette();

            List<HFSI_Control> controlFiles = new List<HFSI_Control>();

            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 3));
            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 4));
            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 5));
            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 6));
            if(p < 2)
                controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 7));

            foreach(var c in controlFiles)
                c.ImmediatePlotPeriod = -1;

            //foreach(var c in controlFiles) {
            //    Assert.IsTrue(c.SkipSolveAndEvaluateResidual == false);
            //    c.dtFixed = dt;
            //    Assert.IsTrue(c.TimesteppingMode == BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady);
            //    c.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            //    c.NonLinearSolver.ConvergenceCriterion = 0.0; // as accurate as possible
            //    c.NonLinearSolver.MaxSolverIterations = 100; 
            //}

            controlFiles.SolverConvergenceTest_Experimental("SolverConvP" + p,
            (VariableNames.DisplacementX, NormType.L2_embedded, p - 1.5, 0, 100000),
            (VariableNames.DisplacementY, NormType.L2_embedded, p - 1.5, 0, 100000),
            (BoSSS.Solution.NSECommon.VariableNames.Pressure, NormType.L2noMean_embedded, p - 1.5, 0, 100000),
            (BoSSS.Solution.NSECommon.VariableNames.Temperature, NormType.L2noMean_embedded, p - 1.5, 0, 100000),
            (BoSSS.Solution.NSECommon.VariableNames.VelocityX, NormType.L2_embedded, p - 1.5, 0, 100000),
            (BoSSS.Solution.NSECommon.VariableNames.VelocityY, NormType.L2_embedded, p - 1.5, 0, 100000)
            );

        }

        public static void Exact(int p = 2) {
            //double dt = 1.0e200;
            // --test=HFSISolver.Tests.SolidOnlyTests.RotationConvergenceTest

            var Tst = new FSTaylorCouette();

            BoSSS.Solution.Application.DeleteOldPlotFiles();
            List<HFSI_Control> controlFiles = new List<HFSI_Control>();
            //var cs = new ZLS_Control[4];

            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 3));
            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 4));
            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 5));
            controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 6));
            if(p < 2)
                controlFiles.Add(HFSISolver.Tests.FSTC.SmallCircle(Tst, p, 7));

            foreach(var c in controlFiles)
                c.ImmediatePlotPeriod = -1;

            var cs = controlFiles.ToArray();

            //foreach(var c in controlFiles) {
            //    Assert.IsTrue(c.SkipSolveAndEvaluateResidual == false);
            //    c.dtFixed = dt;
            //    Assert.IsTrue(c.TimesteppingMode == BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady);
            //    c.NonLinearSolver.SolverCode = BoSSS.Solution.Control.NonLinearSolverCode.Newton;
            //    c.NonLinearSolver.ConvergenceCriterion = 0.0; // as accurate as possible
            //    c.NonLinearSolver.MaxSolverIterations = 100; 
            //}

            SolverConvergenceTest_Exact(Tst, cs, true, new[] {
                ("DisplacementX", p + 1 - 0.6, -1.568, 3.0),
                ("DisplacementY", p + 1 - 0.6, -1.568, 3.0),
                ("VelocityX", p + 1 - 0.5, -1.04, 3.0), 
                ("VelocityY", p + 1 - 0.5, -1.04, 3.0),
                ("Pressure", p - 0.5, -0.507, 3.0),
                ("Temperature", p + 1 - 0.5, -0.507, 3.0)}
                );

        }

        private static (string Name, double slope, double Intercept)[] SolverConvergenceTest_Exact(IHFSITest Tst, HFSI_Control[] CS, bool useExactSolution, (string Name, double Slope, double intercept, double interceptTol)[] RegResults) {
            int D = Tst.SpatialDimension;
            //var CS = __CS.ToArray();
            int NoOfMeshes = CS.Length;

            foreach(var c in CS)
                c.ImmediatePlotPeriod = -1;

            //if(RegResults.Length != D + 1)
            //    throw new ArgumentException("Expecting slopes for velocity and pressure.");

            var Ret = new List<(string Name, double slope, double Intercept)>();

            if(useExactSolution) {
                if(NoOfMeshes < 2)
                    throw new ArgumentException("At least two meshes required for convergence against exact solution.");

                MultidimensionalArray errorS = null;
                string[] Names = null;

                double[] hS = new double[NoOfMeshes];
                HFSI[] solvers = new HFSI[NoOfMeshes];
                //IApplication[] solvers = new IApplication[NoOfMeshes];

                for(int k = 0; k < CS.Length; k++) {

                    var C = CS[k];
                    //using(var solver = new XNSE()) {
                    //var solverClass = C.GetSolverType();
                    //IApplication solver = (IApplication)Activator.CreateInstance(solverClass);
                    var solver = new HFSI();
                    solvers[k] = solver;

                    {
                        //Console.WriteLine("Warning! - enabled immediate plotting");
                        //C.ImmediatePlotPeriod = 1;
                        //C.SuperSampling = 3;

                        solver.Init(C);
                        solver.RunSolverMode();

                        //-------------------Evaluate Error ---------------------------------------- 
                        var evaluator = new ZLSErrorEvaluator<HFSI_Control>(solver);
                        double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);
                        double[] ErrThresh = Tst.AcceptableL2Error;


                        if(k == 0) {
                            errorS = MultidimensionalArray.Create(NoOfMeshes, LastErrors.Length);
                            Names = new string[LastErrors.Length];
                            if(RegResults.Length != Names.Length)
                                throw new ArgumentOutOfRangeException();
                        } else {
                            if(LastErrors.Length != Names.Length)
                                throw new ApplicationException();
                        }

                        if(LastErrors.Length != ErrThresh.Length)
                            throw new ApplicationException();
                        for(int i = 0; i < ErrThresh.Length; i++) {
                            Console.WriteLine($"L2 error, '{solver.Operator.DomainVar[i]}': \t{LastErrors[i]}");
                            Names[i] = solver.Operator.DomainVar[i];
                        }

                        errorS.SetRow(k, LastErrors);
                        hS[k] = evaluator.GetGrid_h();
                    }

                }

                for(int i = 0; i < errorS.GetLength(1); i++) {
                    var RegModel = hS.LogLogRegression(errorS.GetColumn(i));
                    var RegRef = RegResults.Single(ttt => ttt.Name == Names[i]);
                    Console.WriteLine($"Convergence slope for Error of '{Names[i]}': \t{RegModel.Slope}\tIntercept: \t{RegModel.Intercept}\t(Expecting: {RegRef.Slope}, {RegRef.intercept}+/-{RegRef.interceptTol})");
                    Ret.Add((Names[i], RegModel.Slope, RegModel.Intercept));
                }

                for(int i = 0; i < errorS.GetLength(1); i++) {
                    var RegModel = hS.LogLogRegression(errorS.GetColumn(i));
                    var RegRef = RegResults.Single(ttt => ttt.Name == Names[i]);


                    Assert.GreaterOrEqual(RegModel.Slope, RegRef.Slope, $"Convergence Slope of {Names[i]} is degenerate.");
                    Assert.GreaterOrEqual(RegModel.Intercept, RegRef.intercept - RegRef.interceptTol, $"Convergence Intercept of {Names[i]} is degenerate.");
                    Assert.LessOrEqual(RegModel.Intercept, RegRef.intercept + RegRef.interceptTol, $"Convergence Intercept of {Names[i]} overshoot.");
                }

                foreach(var s in solvers) {
                    s.Dispose();
                }
            } else {
                if(NoOfMeshes < 3)
                    throw new ArgumentException("At least three meshes required for convergence if finest solution is assumed to be exact.");



                //BoSSS.Solution.Statistic.ConvergenceTest.SolverConvergenceTest_Experimental(
                //    CS,
                //    "Experimental Convergence",
                //    (D + 1).ForLoop(iVar => (iVar < D ? VariableNames.Velocity_d(iVar) : VariableNames.Pressure,
                //                             iVar < D ? NormType.L2_approximate : NormType.L2noMean_approximate,
                //                             RegResults[iVar].Slope,
                //                             RegResults[iVar].intercept, RegResults[iVar].interceptTol)));

            }



            return Ret.ToArray();


        }



    }
}
