using BoSSS.Foundation;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using System;
using System.Linq;
using static BoSSS.Solution.AdvancedSolvers.Testing.ConditionNumberScalingTest;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.XNSECommon;
using System.Security.Policy;

namespace StokesHelical_Ak.TestSpartial {

    [TestFixture]
    static public class TestSpatial {

        [Test]
        static public void SpatialConvergence_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialConvergence_without_R0fix

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] urErrorLx = new double[gridSize.Length];
            double[] uetaErrorLx = new double[gridSize.Length];
            double[] uxiErrorLx = new double[gridSize.Length];
            double[] psiErrorLx = new double[gridSize.Length];
            HelicalControl[] spaceConvergence = new HelicalControl[gridSize.Length];

            for(int ell = 0; ell < r_min.Length; ell++) {
                {
                    for(int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence[i] = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell]);
                    }

                    Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");
                    for(int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(spaceConvergence[i].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                    }

                    for(int i = 0; i < gridSize.Length; i++) {
                        var solver = new HelicalMain();
                        solver.Init(spaceConvergence[i]);
                        solver.RunSolverMode();
                        Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                        Assert.IsTrue(spaceConvergence[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(spaceConvergence[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                        urErrorLx[i] = solver.urErrorLx;
                        uetaErrorLx[i] = solver.uetaErrorLx;
                        uxiErrorLx[i] = solver.uxiErrorLx;
                        psiErrorLx[i] = solver.psiErrorLx;
                        h[i] = 2 * Math.PI / gridSize[i];
                    }

                    using(var gp = new Gnuplot()) {
                        gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                        gp.OutputFile = $"SpatialConvergence_with_R0fix-p{pOrder}.png";
                        gp.SetTitle($"SpatialConvergence_with_R0fix-p{pOrder}");
                        gp.PlotLogXLogY(h, urErrorLx, title: "ur", format: new PlotFormat("-sr"));
                        gp.PlotLogXLogY(h, uetaErrorLx, title: "ueta", format: new PlotFormat("-^b"));
                        gp.PlotLogXLogY(h, uxiErrorLx, title: "uxi", format: new PlotFormat("-xm"));
                        gp.PlotLogXLogY(h, psiErrorLx, title: "psi", format: new PlotFormat("-*k"));
                        gp.RunAndExit();
                        //Console.WriteLine("endless loop");
                        //while(true) ;
                    }

                    (double slope, double intercept) regressionUr = DoubleExtensions.LogLogRegression(h, urErrorLx);
                    (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(h, uetaErrorLx);
                    (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(h, uxiErrorLx);
                    (double slope, double intercept) regressionPsi = DoubleExtensions.LogLogRegression(h, psiErrorLx);

                    Console.WriteLine("Order_" + pOrder);
                    Console.WriteLine("slopeUR_" + regressionUr.slope);
                    Console.WriteLine("slopeUETA_" + regressionUeta.slope);
                    Console.WriteLine("slopeUXI_" + regressionUxi.slope);
                    Console.WriteLine("slopePSI_" + regressionPsi.slope);

                    double tolerance = 0.6;
                    Assert.GreaterOrEqual(regressionUr.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ur lower than expected (was {0} but should be {1} +- {2})", regressionUr.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionPsi.slope, pOrder - tolerance, String.Format("Convergence rate for Pressure lower than expected (was {0} but should be {1} +- {2})", regressionPsi.slope, pOrder, tolerance));
                }

            }




            if(pOrder == 5) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // error thresholds are specified for the finest grid with pOrder == 5
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                double thresholdPsi = 5e-9;
                Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h.First(), r_min.Last(), psiErrorLx.First());
                Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx.First(), thresholdPsi);
                Assert.LessOrEqual(psiErrorLx.First(), thresholdPsi, "Error. psiErrorLx not fulfilled");

                double thresholdUr = 1e-11;
                Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h.First(), r_min.Last(), urErrorLx.First());
                Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx.First(), thresholdUr);
                Assert.LessOrEqual(urErrorLx.First(), thresholdUr, "Error. urErrorL2 not fulfilled");

                double thresholdUeta = 1e-11;
                Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h.First(), r_min.Last(), uetaErrorLx.First());
                Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx.First(), thresholdUeta);
                Assert.LessOrEqual(uetaErrorLx.First(), thresholdUeta, "Error. uetaErrorL2 not fulfilled");

                double thresholdUxi = 1e-11;
                Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h.First(), r_min.Last(), uxiErrorLx.First());
                Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx.First(), thresholdUxi);
                Assert.LessOrEqual(uxiErrorLx.First(), thresholdUxi, "Error. uxiErrorL2 not fulfilled");
            }
        }

        [Test]
        static public void SpatialComparison_Direct_vs_Iterativ_with_R0fix([Values(5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_with_R0fix(5)

            //###########################################################
            // Direct Solver
            //###########################################################

            int gridSize = 64;
            double h;
            double r_min = 0;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_direct = new HelicalControl();
            HelicalControl spaceConvergence_iterative = new HelicalControl();

            spaceConvergence_direct = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min);

            var helical_direct = new HelicalMain();
            helical_direct.Init(spaceConvergence_direct);
            helical_direct.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_direct.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_direct.PressureReferencePoint, "Pressure Reference Point has to be true");
            urErrorLx = helical_direct.urErrorLx;
            uetaErrorLx = helical_direct.uetaErrorLx;
            uxiErrorLx = helical_direct.uxiErrorLx;
            psiErrorLx = helical_direct.psiErrorLx;
            h = 2 * Math.PI / gridSize;

            double thresholdPsi = 5e-9;
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi, "Error. psiErrorLx not fulfilled");

            double thresholdUr = 1e-11;
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr);
            Assert.LessOrEqual(urErrorLx, thresholdUr, "Error. urErrorL2 not fulfilled");

            double thresholdUeta = 1e-11;
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi = 1e-11;
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi, "Error. uxiErrorL2 not fulfilled");

            //###########################################################
            // Iterativ Solver
            //###########################################################

            spaceConvergence_iterative = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min);
            spaceConvergence_iterative.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { };
            var helical_iterativ = new HelicalMain();
            helical_iterativ.Init(spaceConvergence_iterative);
            helical_iterativ.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_iterative.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_iterative.PressureReferencePoint, "Pressure Reference Point has to be true");
            urErrorLx = helical_iterativ.urErrorLx;
            uetaErrorLx = helical_iterativ.uetaErrorLx;
            uxiErrorLx = helical_iterativ.uxiErrorLx;
            psiErrorLx = helical_iterativ.psiErrorLx;
            h = 2 * Math.PI / gridSize;

            double thresholdPsi_ = 5e-9;
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            double thresholdUr_ = 1e-10;
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr_);
            Assert.LessOrEqual(urErrorLx, thresholdUr_, "Error. urErrorL2 not fulfilled");

            double thresholdUeta_ = 2e-8;
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi_ = 1e-8;
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");

            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterativ.ur, helical_iterativ.ueta, helical_iterativ.uxi, helical_iterativ.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for(int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            double diff_threshold = 4e-7;
            Console.WriteLine("L2 norm of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm());
            Console.WriteLine("If L2 norm of the difference is {0} < {1} than good :) ", diff.L2Norm(), diff_threshold);
            Assert.LessOrEqual(diff.L2Norm(), diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }


        [Test]
        static public void HangingNodes_with_R0fix([Values(4)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.HangingNodes_with_R0fix(5)

            int gridSize = 64;

            double r_min = 0;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_hangingNodes = new HelicalControl();
            spaceConvergence_hangingNodes = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Hanging_Nodes(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min);
            var helical_hangingNodes = new HelicalMain();
            helical_hangingNodes.Init(spaceConvergence_hangingNodes);
            helical_hangingNodes.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_hangingNodes.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_hangingNodes.PressureReferencePoint, "Pressure Reference Point has to be true");
            urErrorLx = helical_hangingNodes.urErrorLx;
            uetaErrorLx = helical_hangingNodes.uetaErrorLx;
            uxiErrorLx = helical_hangingNodes.uxiErrorLx;
            psiErrorLx = helical_hangingNodes.psiErrorLx;

            double thresholdPsi_ = 1e-5;
            Console.WriteLine("The psiErrorLx error for rMin = {0} is = {1}", r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            double thresholdUr_ = 5e-8;
            Console.WriteLine("The urErrorLx error for  rMin =  {0} is = {1}", r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr_);
            Assert.LessOrEqual(urErrorLx, thresholdUr_, "Error. urErrorL2 not fulfilled");

            double thresholdUeta_ = 5e-8;
            Console.WriteLine("The uetaErrorLx error for  rMin =  {0} is = {1}", r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi_ = 5e-8;
            Console.WriteLine("The uxiErrorL2 error for rMin =  {0} is = {1}", r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");
        }


        [Test]
        static public void ConditionNumberScaling_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {


            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for(int ell = 0; ell < r_min.Length; ell++) {
                {
                    for(int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell]);
                    }

                    // StencilCondNo-bndyUncut
                    var conf = new ConditionNumberScalingTest.Config() {
                        plot = true,
                        title = "ConditionNumberScaling_with_R0fix-p" + pOrder,
                        ComputeGlobalCondNo = false
                    };

                    // probably because of the R0fix, the stencil numbers at the r=0 - boundary go crazy.
                    conf.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_bndyUncut] = (XAxisDesignation.Grid_1Dres, 0.5, -1);

                    ConditionNumberScalingTest.Perform(gridSeries, conf);
                    for(int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(gridSeries[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(gridSeries[i].steady, "Should be steady");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                    }

                }
            }
        }


        [Test]
        static public void SpatialConvergence_without_R0fix([Values(2, 3, 4, 5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialConvergence_without_R0fix 

            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //Console.WriteLine("remove testcode");
            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

            ilPSP.Environment.NumThreads = 1;

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0.1 };
            double[] urErrorL2 = new double[gridSize.Length];
            double[] uetaErrorL2 = new double[gridSize.Length];
            double[] uxiErrorL2 = new double[gridSize.Length];
            double[] psiErrorL2 = new double[gridSize.Length];
            HelicalControl[] spaceConvergence = new HelicalControl[gridSize.Length];

            for(int ell = 0; ell < r_min.Length; ell++) {
                {
                    for(int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence[i] = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell]);
                    }
                    Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");

                    for(int i = 0; i < gridSize.Length; i++) {
                        var solver = new HelicalMain();
                        solver.Init(spaceConvergence[i]);
                        solver.RunSolverMode();
                        Assert.IsFalse(spaceConvergence[i].R0fixOn, "R0fix should be false");
                        Assert.IsTrue(spaceConvergence[i].PressureReferencePoint, "PressureReferencePoint should be true. Since R0_fix is false");
                        Assert.IsTrue(spaceConvergence[i].steady, "Should be steady");
                        Assert.AreEqual(Globals.activeMult, Globals.Multiplier.one, $"Multiplier expected to be {Globals.Multiplier.one}");
                        urErrorL2[i] = solver.urErrorL2;
                        uetaErrorL2[i] = solver.uetaErrorL2;
                        uxiErrorL2[i] = solver.uxiErrorL2;
                        psiErrorL2[i] = solver.psiErrorL2;
                        h[i] = 2 * Math.PI / gridSize[i];
                    }

                    using(var gp = new Gnuplot()) {
                        gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                        gp.OutputFile = $"SpatialConvergence_without_R0fix-p{pOrder}.png";
                        gp.SetTitle($"SpatialConvergence_without_R0fix-p{pOrder}");
                        gp.PlotLogXLogY(h, urErrorL2, title: "ur", format: new PlotFormat("-sr"));
                        gp.PlotLogXLogY(h, uetaErrorL2, title: "ueta", format: new PlotFormat("-^b"));
                        gp.PlotLogXLogY(h, uxiErrorL2, title: "uxi", format: new PlotFormat("-xm"));
                        gp.PlotLogXLogY(h, psiErrorL2, title: "psi", format: new PlotFormat("-*k"));
                        gp.RunAndExit();
                        //Console.WriteLine("endless loop");
                        //while(true) ;
                    }


                    (double slope, double intercept) regressionUr = DoubleExtensions.LogLogRegression(h, urErrorL2);
                    (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(h, uetaErrorL2);
                    (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(h, uxiErrorL2);
                    (double slope, double intercept) regressionPsi = DoubleExtensions.LogLogRegression(h, psiErrorL2);

                    Console.WriteLine("Order_" + pOrder);
                    Console.WriteLine("slopeUR_" + regressionUr.slope);
                    Console.WriteLine("slopeUETA_" + regressionUeta.slope);
                    Console.WriteLine("slopeUXI_" + regressionUxi.slope);
                    Console.WriteLine("slopePSI_" + regressionPsi.slope);

                    double tolerance = pOrder <= 3 ? 0.5 : 0.8;
                    Console.WriteLine("If rate " + regressionUr.slope + " >= " + pOrder + " minus tolerance " + tolerance + " than good");
                    Assert.GreaterOrEqual(regressionUr.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ur lower than expected (was {0} but should be {1} +- {2})", regressionUr.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionPsi.slope, pOrder - tolerance, String.Format("Convergence rate for Pressure lower than expected (was {0} but should be {1} +- {2})", regressionPsi.slope, pOrder, tolerance));
                }

            }

            if(pOrder == 5) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // error thresholds are specified for the finest grid with pOrder == 5
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                double thresholdPsi = 5e-9;
                Console.WriteLine("The psiErrorL2 error for {0} xi-Cells and rMin = {1} is = {2}", h.First(), r_min.Last(), psiErrorL2.First());
                Console.WriteLine("If the psiErrorL2 error is {0} < {1} than good :) ", psiErrorL2.First(), thresholdPsi);
                Assert.LessOrEqual(psiErrorL2.First(), thresholdPsi, "Error. psiErrorLx not fulfilled");

                double thresholdUr = 1e-11;
                Console.WriteLine("The urErrorL2 error for {0} xi-Cells and {1} is = {2}", h.First(), r_min.Last(), urErrorL2.First());
                Console.WriteLine("If the urErrorL2 error is {0} < {1} than good :) ", urErrorL2.First(), thresholdUr);
                Assert.LessOrEqual(urErrorL2.First(), thresholdUr, "Error. urErrorL2 not fulfilled");

                double thresholdUeta = 1e-11;
                Console.WriteLine("The uetaErrorL2 error for {0} xi-Cells and {1} is = {2}", h.First(), r_min.Last(), uetaErrorL2.First());
                Console.WriteLine("If the uetaErrorL2 error is {0} < {1} than good :) ", uetaErrorL2.First(), thresholdUeta);
                Assert.LessOrEqual(uetaErrorL2.First(), thresholdUeta, "Error. uetaErrorL2 not fulfilled");

                double thresholdUxi = 1e-11;
                Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h.First(), r_min.Last(), uxiErrorL2.First());
                Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorL2.First(), thresholdUxi);
                Assert.LessOrEqual(uxiErrorL2.First(), thresholdUxi, "Error. uxiErrorL2 not fulfilled");
            }
        }

        [Test]
        static public void SpatialComparison_Direct_vs_Iterativ_without_R0fix([Values(5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_without_R0fix(5)

            //###########################################################
            // Direct Solver
            //###########################################################

            int gridSize = 64;
            double h;
            double r_min = 0.1;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_direct = new HelicalControl();
            HelicalControl spaceConvergence_iterative = new HelicalControl();

            spaceConvergence_direct = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min);

            var helical_direct = new HelicalMain();
            helical_direct.Init(spaceConvergence_direct);
            helical_direct.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.one, $"Multiplier expected to be {Globals.Multiplier.one}");
            Assert.IsFalse(spaceConvergence_direct.R0fixOn, "R0fix must be turned off");
            Assert.IsTrue(spaceConvergence_direct.PressureReferencePoint, "Pressure Reference Point has to be true");
            urErrorLx = helical_direct.urErrorLx;
            uetaErrorLx = helical_direct.uetaErrorLx;
            uxiErrorLx = helical_direct.uxiErrorLx;
            psiErrorLx = helical_direct.psiErrorLx;
            h = 2 * Math.PI / gridSize;

            double thresholdPsi = 5e-9;
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi, "Error. psiErrorLx not fulfilled");

            double thresholdUr = 1e-11;
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr);
            Assert.LessOrEqual(urErrorLx, thresholdUr, "Error. urErrorL2 not fulfilled");

            double thresholdUeta = 1e-11;
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi = 1e-11;
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi, "Error. uxiErrorL2 not fulfilled");

            //###########################################################
            // Iterativ Solver
            //###########################################################

            spaceConvergence_iterative = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min);
            spaceConvergence_iterative.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { };
            var helical_iterativ = new HelicalMain();
            helical_iterativ.Init(spaceConvergence_iterative);
            helical_iterativ.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.one, $"Multiplier expected to be {Globals.Multiplier.one}");
            Assert.IsFalse(spaceConvergence_iterative.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_iterative.PressureReferencePoint, "Pressure Reference Point has to be true");
            urErrorLx = helical_iterativ.urErrorLx;
            uetaErrorLx = helical_iterativ.uetaErrorLx;
            uxiErrorLx = helical_iterativ.uxiErrorLx;
            psiErrorLx = helical_iterativ.psiErrorLx;
            h = 2 * Math.PI / gridSize;

            double thresholdPsi_ = 5e-9;
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            double thresholdUr_ = 1e-10;
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr_);
            Assert.LessOrEqual(urErrorLx, thresholdUr_, "Error. urErrorL2 not fulfilled");

            double thresholdUeta_ = 2e-8;
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi_ = 1e-8;
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");

            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterativ.ur, helical_iterativ.ueta, helical_iterativ.uxi, helical_iterativ.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for(int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            double diff_threshold = 4e-7;
            Console.WriteLine("L2 norm of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm());
            Console.WriteLine("If L2 norm of the difference is {0} < {1} than good :) ", diff.L2Norm(), diff_threshold);
            Assert.LessOrEqual(diff.L2Norm(), diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }

        [Test]
        static public void HangingNodes_without_R0fix([Values(4)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.HangingNodes_with_R0fix(5)

            int gridSize = 64;
            double r_min = 0.1;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_hangingNodes = new HelicalControl();
            spaceConvergence_hangingNodes = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Hanging_Nodes(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min);
            var helical_hangingNodes = new HelicalMain();
            helical_hangingNodes.Init(spaceConvergence_hangingNodes);
            helical_hangingNodes.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.one, $"Multiplier expected to be {Globals.Multiplier.one}");
            Assert.IsFalse(spaceConvergence_hangingNodes.R0fixOn, "R0fix must be turned off");
            Assert.IsTrue(spaceConvergence_hangingNodes.PressureReferencePoint, "Pressure Reference Point has to be true");
            urErrorLx = helical_hangingNodes.urErrorLx;
            uetaErrorLx = helical_hangingNodes.uetaErrorLx;
            uxiErrorLx = helical_hangingNodes.uxiErrorLx;
            psiErrorLx = helical_hangingNodes.psiErrorLx;

            double thresholdPsi_ = 1e-5;
            Console.WriteLine("The psiErrorLx error for rMin = {0} is = {1}", r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            double thresholdUr_ = 5e-8;
            Console.WriteLine("The urErrorLx error for {0} is = {1}", r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr_);
            Assert.LessOrEqual(urErrorLx, thresholdUr_, "Error. urErrorL2 not fulfilled");

            double thresholdUeta_ = 5e-8;
            Console.WriteLine("The uetaErrorLx error for {0} is = {1}", r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi_ = 5e-8;
            Console.WriteLine("The uxiErrorL2 error for {0} is = {1}", r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");
        }
        [Test]
        static public void ConditionNumberScaling_without_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            //ilPSP.Environment.NumThreads = 1;

            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0.1 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for(int ell = 0; ell < r_min.Length; ell++) {
                {
                    for(int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = StokesHelical_Ak.HardcodedControl.ManSol_Steady_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell]);
                    }

                    ConditionNumberScalingTest.Perform(gridSeries, new ConditionNumberScalingTest.Config() { plot = true, title = "ConditionNumberScaling_without_R0fix-p" + pOrder, ComputeGlobalCondNo = false });
                    for(int i = 0; i < gridSize.Length; i++) {
                        Assert.IsFalse(gridSeries[i].R0fixOn, "R0fix must be turned off");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "PressureReferencePoint should be true. Since R0_fix is false");
                    }
                }
            }
        }

        /// <summary>
        /// Tests a steady-state Hagen Poiseulle flow (aka. Pipe flow).
        /// The solver computes a stationary Stokes solution, this is achieved by:
        /// - setting the initial values to 0.0, so that the convective terms, which are explicitly treated, in the first time-step are 0
        /// - computing one very long time-step, so that the (1/dt)-contributions and the initial value are negligible
        /// Beside spatial convergence, one also expects the radial velocity and the pressure to be close to 0.0.
        /// </summary>
        /// <remarks>
        /// Note:
        /// - there is no Navier-Stokes version of this, since the explicit treatment of the convective terms with a large time-step is not numerically stable.
        /// </remarks>
        [Test]
        static public void SteadyHagenPoiseulle_Stokes([Values(2, 3, 4, 5)] int pOrder) {

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] ur_L2_Error = new double[gridSize.Length];
            double[] ueta_L2_Error = new double[gridSize.Length];
            double[] uxi_L2_Error = new double[gridSize.Length];
            double[] pressure_L2_Error = new double[gridSize.Length];

            HelicalControl[] spaceConvergence_HP = new HelicalControl[gridSize.Length];
            SinglePhaseField[] exSol_ur = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_ueta = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_uxi = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_pressure = new SinglePhaseField[gridSize.Length];

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            for(int ell = 0; ell < r_min.Length; ell++) {
                //ilPSP.Environment.NumThreads = 1;
                for(int i = 0; i < gridSize.Length; i++) {
                    spaceConvergence_HP[i] = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50);
                }

                Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence_HP.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");
                for(int i = 0; i < gridSize.Length; i++) {
                    Assert.IsTrue(spaceConvergence_HP[i].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                    spaceConvergence_HP[i].InitialValues.Clear();
                    spaceConvergence_HP[i].InitialValues_Evaluators.Clear();
                    spaceConvergence_HP[i].savetodb = false;
                    spaceConvergence_HP[i].TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                    // Exakte Loesung
                    double MaxAmp = spaceConvergence_HP[i].maxAmpli;
                    Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (spaceConvergence_HP[i].rMax * spaceConvergence_HP[i].rMax - X[0] * X[0])) / (4 * nu);
                    Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (spaceConvergence_HP[i].rMax * spaceConvergence_HP[i].rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                    Func<double[], double> ur = (X) => 0;
                    Func<double[], double> pressure = (X) => 0;

                    using(var solver = new HelicalMain()) {
                        solver.Init(spaceConvergence_HP[i]);
                        solver.RunSolverMode();
                        Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                        Assert.IsTrue(spaceConvergence_HP[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(spaceConvergence_HP[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                        h[i] = 2 * Math.PI / gridSize[i];
                        // Deep Copy
                        exSol_ur[i] = solver.ur.CloneAs();
                        exSol_ueta[i] = solver.ueta.CloneAs();
                        exSol_uxi[i] = solver.uxi.CloneAs();
                        exSol_pressure[i] = solver.Pressure.CloneAs();

                        //Exakte Lösung auf SinglePhaseField
                        exSol_ur[i].ProjectField(ur);
                        exSol_ueta[i].ProjectField(ueta);
                        exSol_uxi[i].ProjectField(uxi);
                        exSol_pressure[i].ProjectField(pressure);

                        //L2 Fehler berechnen. Meine Loesung vs richtige Loesung
                        ur_L2_Error[i] = solver.ur.L2Error(exSol_ur[i]);
                        uxi_L2_Error[i] = solver.uxi.L2Error(exSol_uxi[i]);
                        ueta_L2_Error[i] = solver.ueta.L2Error(exSol_ueta[i]);
                        pressure_L2_Error[i] = solver.Pressure.L2Error(exSol_pressure[i]);


                        Console.WriteLine($"ur       L2 Error: {ur_L2_Error[i]:0.###e-00} (should be close to 0.0)");
                        Console.WriteLine($"uxi      L2 Error: {uxi_L2_Error[i]:0.###e-00}");
                        Console.WriteLine($"ueta     L2 Error: {ueta_L2_Error[i]:0.###e-00}");
                        Console.WriteLine($"pressure L2 Error: {pressure_L2_Error[i]:0.###e-00} (should be close to 0.0)");


                        Assert.LessOrEqual(ur_L2_Error[i], 1.0e-11, $"ur L2 Error out of range: {ur_L2_Error[i]:0.###e-00} (should be close to 0.0 )");
                        Assert.LessOrEqual(pressure_L2_Error[i], 1.0e-10, $"pressure L2 Error out of range: {pressure_L2_Error[i]:0.###e-00} (should be close to 0.0)");
                    }
                }

                using(var gp = new Gnuplot()) {
                    gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                    gp.OutputFile = $"SpatialConvergence_with_R0fix-p{pOrder}.png";
                    gp.SetTitle($"SpatialConvergence_Hagen_Poiseulle_with_R0fix-p{pOrder}");
                    gp.PlotLogXLogY(h, ur_L2_Error, title: "ur", format: new PlotFormat("-sr"));
                    gp.PlotLogXLogY(h, ueta_L2_Error, title: "ueta", format: new PlotFormat("-^b"));
                    gp.PlotLogXLogY(h, uxi_L2_Error, title: "uxi", format: new PlotFormat("-xm"));
                    gp.PlotLogXLogY(h, pressure_L2_Error, title: "psi", format: new PlotFormat("-*k"));
                    gp.RunAndExit();
                    //Console.WriteLine("endless loop");
                    //while(true) ;
                }

                (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(h, ueta_L2_Error);
                (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(h, uxi_L2_Error);

                Console.WriteLine("Order_" + pOrder);
                Console.WriteLine("slopeUETA_" + regressionUeta.slope);
                Console.WriteLine("slopeUXI_" + regressionUxi.slope);

                double tolerance = 0.75;

                // Since Ur and pressure are zero, there are no slops
                // Therefore, looking at nslopes for Ueta und Uxi
                if(pOrder != 5) {
                    Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                }
                if(pOrder == 5) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // error thresholds are specified for the SECOUND!!!!! finest grid with pOrder == 5
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    Assert.LessOrEqual(ur_L2_Error[1], 1.0e-10, $"ur L2 Error out of range: {ur_L2_Error[1]:0.###e-00} (should be close to 0.0 )");
                    Assert.LessOrEqual(uxi_L2_Error[1], 1.0e-10, $"uxi L2 Error out of range: {uxi_L2_Error[1]:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2_Error[1], 1.0e-10, $"ueta L2 Error out of range: {ueta_L2_Error[1]:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(pressure_L2_Error[1], 1.0e-10, $"pressure L2 Error out of range: {pressure_L2_Error[1]:0.###e-00} (should be close to 0.0)");
                }
            }
        }


        /// <summary>
        /// Tests a steady-state Centrifuge flow (aka. Centrifuge flow).
        /// The solver computes a stationary Stokes solution, this is achieved by:
        /// - setting the initial values to 0.0, so that the convective terms, which are explicitly treated, in the first time-step are 0
        /// - computing one very long time-step, so that the (1/dt)-contributions and the initial value are negligible
        /// Beside spatial convergence, one also expects the radial velocity and the pressure to be close to 0.0.
        /// </summary>
        /// <remarks>
        /// Note:
        /// -1) There is no Navier-Stokes version of this, since the explicit treatment of the convective terms with a large time-step is not numerically stable.
        /// -2) Attention!! Steady Stokes Solution for pressure (p=0) is not the same as for the Navier Stokes Solution (p=MaxAmp * MaxAmp * r * r * 0.5). 
        /// -3) Here Only the Steady-state Solution is looked at. 
        /// </remarks>
        [Test]
        static public void SteadyCentrifuge_Stokes([Values(2, 3, 4, 5)] int pOrder) {

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] ur_L2_Error = new double[gridSize.Length];
            double[] ueta_L2_Error = new double[gridSize.Length];
            double[] uxi_L2_Error = new double[gridSize.Length];
            double[] pressure_L2_Error = new double[gridSize.Length];

            HelicalControl[] spaceConvergence_HP = new HelicalControl[gridSize.Length];
            SinglePhaseField[] exSol_ur = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_ueta = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_uxi = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_pressure = new SinglePhaseField[gridSize.Length];

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            for(int ell = 0; ell < r_min.Length; ell++) {
                //ilPSP.Environment.NumThreads = 1;
                for(int i = 0; i < gridSize.Length; i++) {
                    spaceConvergence_HP[i] = StokesHelical_Ak.DNS_Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50);
                }

                Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence_HP.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");
                for(int i = 0; i < gridSize.Length; i++) {
                    Assert.IsTrue(spaceConvergence_HP[i].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                    spaceConvergence_HP[i].InitialValues.Clear();
                    spaceConvergence_HP[i].InitialValues_Evaluators.Clear();
                    spaceConvergence_HP[i].savetodb = false;
                    spaceConvergence_HP[i].TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                    // Exact Solution for steady state!
                    double MaxAmp = spaceConvergence_HP[i].maxAmpli;
                    Func<double[], double> ur = (X) => 0;
                    Func<double[], double> pressure = (X) => 0;
                    Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                    Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;

                    using(var solver = new HelicalMain()) {
                        solver.Init(spaceConvergence_HP[i]);
                        solver.RunSolverMode();
                        Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                        Assert.IsTrue(spaceConvergence_HP[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(spaceConvergence_HP[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                        h[i] = 2 * Math.PI / gridSize[i];
                        // Deep Copy
                        exSol_ur[i] = solver.ur.CloneAs();
                        exSol_ueta[i] = solver.ueta.CloneAs();
                        exSol_uxi[i] = solver.uxi.CloneAs();
                        exSol_pressure[i] = solver.Pressure.CloneAs();

                        //Exakte Lösung auf SinglePhaseField
                        exSol_ur[i].ProjectField(ur);
                        exSol_ueta[i].ProjectField(ueta);
                        exSol_uxi[i].ProjectField(uxi);
                        exSol_pressure[i].ProjectField(pressure);

                        //L2 Fehler berechnen. Meine Loesung vs richtige Loesung
                        ur_L2_Error[i] = solver.ur.L2Error(exSol_ur[i]);
                        uxi_L2_Error[i] = solver.uxi.L2Error(exSol_uxi[i]);
                        ueta_L2_Error[i] = solver.ueta.L2Error(exSol_ueta[i]);
                        pressure_L2_Error[i] = solver.Pressure.L2Error(exSol_pressure[i]);


                        Console.WriteLine($"ur       L2 Error: {ur_L2_Error[i]:0.###e-00} (should be close to 0.0)");
                        Console.WriteLine($"uxi      L2 Error: {uxi_L2_Error[i]:0.###e-00}");
                        Console.WriteLine($"ueta     L2 Error: {ueta_L2_Error[i]:0.###e-00}");
                        Console.WriteLine($"pressure L2 Error: {pressure_L2_Error[i]:0.###e-00} (should be close to 0.0)");


                        Assert.LessOrEqual(ur_L2_Error[i], 1.0e-11, $"ur L2 Error out of range: {ur_L2_Error[i]:0.###e-00} (should be close to 0.0 )");
                        Assert.LessOrEqual(pressure_L2_Error[i], 1.0e-10, $"Pressure L2 Error out of range: {pressure_L2_Error[i]:0.###e-00} (should be close to 0.0 )");
                    }
                }

                using(var gp = new Gnuplot()) {
                    gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                    gp.OutputFile = $"SpatialConvergence_with_R0fix-p{pOrder}.png";
                    gp.SetTitle($"SpatialConvergence_Hagen_Poiseulle_with_R0fix-p{pOrder}");
                    gp.PlotLogXLogY(h, ur_L2_Error, title: "ur", format: new PlotFormat("-sr"));
                    gp.PlotLogXLogY(h, ueta_L2_Error, title: "ueta", format: new PlotFormat("-^b"));
                    gp.PlotLogXLogY(h, uxi_L2_Error, title: "uxi", format: new PlotFormat("-xm"));
                    gp.PlotLogXLogY(h, pressure_L2_Error, title: "psi", format: new PlotFormat("-*k"));
                    gp.RunAndExit();
                    //Console.WriteLine("endless loop");
                    //while(true) ;
                }

                (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(h, ueta_L2_Error);
                (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(h, uxi_L2_Error);

                Console.WriteLine("Order_" + pOrder);
                Console.WriteLine("slopeUETA_" + regressionUeta.slope);
                Console.WriteLine("slopeUXI_" + regressionUxi.slope);

                double tolerance = 0.75;

                // Since Ur is zero, there are no slop
                // Therefore, looking at slopes for Ueta, Uxi and Pressure
                if(pOrder != 5) {
                    Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                    Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                }
                if(pOrder == 5) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // error thresholds are specified for the finest grid with pOrder == 5
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                    Assert.LessOrEqual(ur_L2_Error[0], 1.0e-10, $"ur L2 Error out of range: {ur_L2_Error[0]:0.###e-00} (should be close to 0.0 )");
                    Assert.LessOrEqual(uxi_L2_Error[0], 1.0e-10, $"uxi L2 Error out of range: {uxi_L2_Error[0]:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2_Error[0], 1.0e-10, $"ueta L2 Error out of range: {ueta_L2_Error[0]:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(pressure_L2_Error[0], 1.0e-10, $"pressure L2 Error out of range: {pressure_L2_Error[0]:0.###e-00} (should be close to 0.0)");
                }
            }
        }
    }
}
