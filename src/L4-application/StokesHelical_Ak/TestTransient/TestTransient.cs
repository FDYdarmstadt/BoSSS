using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using NUnit.Framework.Constraints;
using NUnitLite;
using StokesHelical_Ak.Hard_Coded_Control;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace StokesHelical_Ak.TestTransient {
    [TestFixture]
    static public class TestTransient {

        #region DDD
        /// <summary>
        /// Transient Tests for the Solution in Dr Dominiks Dierkes Paper (Formular 4.4) for BDF 1 and no R0fix
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_TimeConv_BDF1_DDD_no_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF1" };
            for (int ell = 0; ell < timeSchemes.Length; ell++) {
                for (int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = Man_Sol_DDD.ManSol_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0.1);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == false, "R0_fix should be false");
                    Assert.That(Globals.pressureReferencePoint == true, "We have to calculate with PRP, since R0 is off");
                    Assert.That(Globals.activeMult == Globals.Multiplier.one);
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if (timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if (timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }


                using (var gp = new Gnuplot()) {
                    gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                    gp.OutputFile = $"TemporalConvergence_without_R0fix-BDF{order}.png";
                    gp.SetTitle($"TemporalConvergence_without_R0fix-BDF{order}");
                    gp.PlotLogXLogY(timeStepSize, urErrorL2, title: "ur", format: new PlotFormat("-sr"));
                    gp.PlotLogXLogY(timeStepSize, uetaErrorL2, title: "ueta", format: new PlotFormat("-^b"));
                    gp.PlotLogXLogY(timeStepSize, uxiErrorL2, title: "uxi", format: new PlotFormat("-xm"));
                    gp.PlotLogXLogY(timeStepSize, psiErrorL2, title: "psi", format: new PlotFormat("-*k"));
                    gp.RunAndExit();
                    //Console.WriteLine("endless loop");
                    //while(true) ;
                }

                (double slope, double intercept) regressionUr = DoubleExtensions.LogLogRegression(timeStepSize, urErrorL2);
                (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(timeStepSize, uetaErrorL2);
                (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(timeStepSize, uxiErrorL2);
                (double slope, double intercept) regressionPsi = DoubleExtensions.LogLogRegression(timeStepSize, psiErrorL2);

                Console.WriteLine("Order " + order);
                Console.WriteLine("slopeUR " + regressionUr.slope);
                Console.WriteLine("slopeUETA " + regressionUeta.slope);
                Console.WriteLine("slopeUXI " + regressionUxi.slope);
                Console.WriteLine("slopePSI " + regressionPsi.slope);

                double tolerance = 0.5;
                Console.WriteLine("If rate " + regressionUr.slope + " >= " + order + " minus tolerance " + tolerance + " than good");
                Assert.That(regressionUr.slope >= order - tolerance, String.Format("Convergence rate for Ur lower than expected (was {0} but should be {1} +- {2})", regressionUr.slope, order, tolerance));
                Assert.That(regressionUeta.slope >= order - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, order, tolerance));
                Assert.That(regressionUxi.slope >= order - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, order, tolerance));
                Assert.That(regressionPsi.slope >= order - tolerance, String.Format("Convergence rate for Psi lower than expected (was {0} but should be {1} +- {2})", regressionPsi.slope, order, tolerance));
            }
            Console.WriteLine("Remember: R0_fix is " + timeConvergence.First().R0fixOn);
            Console.WriteLine("Remember: Calculating with PRP is " + Globals.pressureReferencePoint);


            double thresholdPsi = 5e-2;
            Console.WriteLine("The psiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), psiErrorL2.First());
            Console.WriteLine("If the psiErrorL2 error is {0} < {1} than good :) ", psiErrorL2.First(), thresholdPsi);
            Assert.That(psiErrorL2.First() < thresholdPsi, "Error. psiErrorLx not fullfilled");

            double thresholdUr = 2e-3;
            Console.WriteLine("The urErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), urErrorL2.First());
            Console.WriteLine("If the urErrorL2 error is {0} < {1} than good :) ", urErrorL2.First(), thresholdUr);
            Assert.That(urErrorL2.First() < thresholdUr, "Error. urErrorL2 not fullfilled");

            double thresholdUeta = 2e-3;
            Console.WriteLine("The uetaErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uetaErrorL2.First());
            Console.WriteLine("If the uetaErrorL2 error is {0} < {1} than good :) ", uetaErrorL2.First(), thresholdUeta);
            Assert.That(uetaErrorL2.First() < thresholdUeta, "Error. uetaErrorL2 not fullfilled");

            double thresholdUxi = 2e-3;
            Console.WriteLine("The uxiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uxiErrorL2.First());
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorL2.First(), thresholdUxi);
            Assert.That(uxiErrorL2.First() < thresholdUxi, "Error. uxiErrorL2 not fullfilled");
        }
        /// <summary>
        /// Transient Tests for the Solution in Dr Dominiks Dierkes Paper (Formular 4.4) for BDF 3 and no R0fix
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_TimeConv_BDF3_DDD_no_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF3" };
            for (int ell = 0; ell < timeSchemes.Length; ell++) {
                for (int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = Man_Sol_DDD.ManSol_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0.1);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == false, "R0_fix should be false");
                    Assert.That(Globals.pressureReferencePoint == true, "We have to calculate with PRP, since R0 is off");
                    Assert.That(Globals.activeMult == Globals.Multiplier.one);
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if (timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if (timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }

                using (var gp = new Gnuplot()) {
                    gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                    gp.OutputFile = $"TemporalConvergence_without_R0fix-BDF{order}.png";
                    gp.SetTitle($"TemporalConvergence_without_R0fix-BDF{order}");
                    gp.PlotLogXLogY(timeStepSize, urErrorL2, title: "ur", format: new PlotFormat("-sr"));
                    gp.PlotLogXLogY(timeStepSize, uetaErrorL2, title: "ueta", format: new PlotFormat("-^b"));
                    gp.PlotLogXLogY(timeStepSize, uxiErrorL2, title: "uxi", format: new PlotFormat("-xm"));
                    gp.PlotLogXLogY(timeStepSize, psiErrorL2, title: "psi", format: new PlotFormat("-*k"));
                    gp.RunAndExit();
                    //Console.WriteLine("endless loop");
                    //while(true) ;
                }


                (double slope, double intercept) regressionUr = DoubleExtensions.LogLogRegression(timeStepSize, urErrorL2);
                (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(timeStepSize, uetaErrorL2);
                (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(timeStepSize, uxiErrorL2);
                (double slope, double intercept) regressionPsi = DoubleExtensions.LogLogRegression(timeStepSize, psiErrorL2);

                Console.WriteLine("Order " + order);
                Console.WriteLine("slopeUR " + regressionUr.slope);
                Console.WriteLine("slopeUETA " + regressionUeta.slope);
                Console.WriteLine("slopeUXI " + regressionUxi.slope);
                Console.WriteLine("slopePSI " + regressionPsi.slope);

                double tolerance = 0.5;
                Console.WriteLine("If rate " + regressionUr.slope + " >= " + order + " minus tolerance " + tolerance + " than good");
                Assert.That(regressionUr.slope >= order - tolerance, String.Format("Convergence rate for Ur lower than expected (was {0} but should be {1} +- {2})", regressionUr.slope, order, tolerance));
                Assert.That(regressionUeta.slope >= order - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, order, tolerance));
                Assert.That(regressionUxi.slope >= order - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, order, tolerance));
                Assert.That(regressionPsi.slope >= order - tolerance, String.Format("Convergence rate for Psi lower than expected (was {0} but should be {1} +- {2})", regressionPsi.slope, order, tolerance));
            }
            Console.WriteLine("Remember: R0_fix is " + timeConvergence.First().R0fixOn);
            Console.WriteLine("Remember: Calculating with PRP is " + Globals.pressureReferencePoint);

            double thresholdPsi = 5e-4;
            Console.WriteLine("The psiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), psiErrorL2.First());
            Console.WriteLine("If the psiErrorL2 error is {0} < {1} than good :) ", psiErrorL2.First(), thresholdPsi);
            Assert.That(psiErrorL2.First() < thresholdPsi, "Error. psiErrorLx not fullfilled");

            double thresholdUr = 2e-5;
            Console.WriteLine("The urErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), urErrorL2.First());
            Console.WriteLine("If the urErrorL2 error is {0} < {1} than good :) ", urErrorL2.First(), thresholdUr);
            Assert.That(urErrorL2.First() < thresholdUr, "Error. urErrorL2 not fullfilled");

            double thresholdUeta = 1e-5;
            Console.WriteLine("The uetaErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uetaErrorL2.First());
            Console.WriteLine("If the uetaErrorL2 error is {0} < {1} than good :) ", uetaErrorL2.First(), thresholdUeta);
            Assert.That(uetaErrorL2.First() < thresholdUeta, "Error. uetaErrorL2 not fullfilled");

            double thresholdUxi = 1e-5;
            Console.WriteLine("The uxiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uxiErrorL2.First());
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorL2.First(), thresholdUxi);
            Assert.That(uxiErrorL2.First() < thresholdUxi, "Error. uxiErrorL2 not fullfilled");
        }

        /// <summary>
        /// Transient Tests for the Solution in Dr Dominiks Dierkes Paper (Formular 4.4) for BDF 1 and with R0fix
        /// </summary>
        /// <remarks>
        /// ATTENTION MUTIPLIER IS ONE BUT should be BSQ (See in the Helical Main). Don't Know why BSQ is not working in this Testcase
        /// </remarks>
        [Test]
        static public void Transient_TimeConv_BDF1_DDD_with_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF1" };
            for (int ell = 0; ell < timeSchemes.Length; ell++) {
                for (int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = Man_Sol_DDD.ManSol_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == true, "R0_fix should be true");
                    Assert.That(Globals.pressureReferencePoint == true, "We have to calculate with PRP");
                    if (Globals.activeMult == Globals.Multiplier.one && timeConvergence[i].rMin < 10e-6) {
                        Console.WriteLine("Friendly Reminder: Mutiplier One and rMin<10e-6");
                    }
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if (timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if (timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }

                using (var gp = new Gnuplot()) {
                    gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                    gp.OutputFile = $"TemporalConvergence_with_R0fix-BDF{order}.png";
                    gp.SetTitle($"TemporalConvergence_with_R0fix-BDF{order}");
                    gp.PlotLogXLogY(timeStepSize, urErrorL2, title: "ur", format: new PlotFormat("-sr"));
                    gp.PlotLogXLogY(timeStepSize, uetaErrorL2, title: "ueta", format: new PlotFormat("-^b"));
                    gp.PlotLogXLogY(timeStepSize, uxiErrorL2, title: "uxi", format: new PlotFormat("-xm"));
                    gp.PlotLogXLogY(timeStepSize, psiErrorL2, title: "psi", format: new PlotFormat("-*k"));
                    gp.RunAndExit();
                    //Console.WriteLine("endless loop");
                    //while(true) ;
                }

                (double slope, double intercept) regressionUr = DoubleExtensions.LogLogRegression(timeStepSize, urErrorL2);
                (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(timeStepSize, uetaErrorL2);
                (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(timeStepSize, uxiErrorL2);
                (double slope, double intercept) regressionPsi = DoubleExtensions.LogLogRegression(timeStepSize, psiErrorL2);

                Console.WriteLine("Order " + order);
                Console.WriteLine("slopeUR " + regressionUr.slope);
                Console.WriteLine("slopeUETA " + regressionUeta.slope);
                Console.WriteLine("slopeUXI " + regressionUxi.slope);
                Console.WriteLine("slopePSI " + regressionPsi.slope);

                double tolerance = 0.5;
                Console.WriteLine("If rate " + regressionUr.slope + " >= " + order + " minus tolerance " + tolerance + " than good");
                Assert.That(regressionUr.slope >= order - tolerance, String.Format("Convergence rate for Ur lower than expected (was {0} but should be {1} +- {2})", regressionUr.slope, order, tolerance));
                Assert.That(regressionUeta.slope >= order - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, order, tolerance));
                Assert.That(regressionUxi.slope >= order - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, order, tolerance));
                Assert.That(regressionPsi.slope >= order - tolerance, String.Format("Convergence rate for Psi lower than expected (was {0} but should be {1} +- {2})", regressionPsi.slope, order, tolerance));
            }
            Console.WriteLine("Remember: R0_fix is " + timeConvergence.First().R0fixOn);
            Console.WriteLine("Remember: Calculating with PRP is " + Globals.pressureReferencePoint);

            double thresholdPsi = 5e-2;
            Console.WriteLine("The psiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), psiErrorL2.First());
            Console.WriteLine("If the psiErrorL2 error is {0} < {1} than good :) ", psiErrorL2.First(), thresholdPsi);
            Assert.That(psiErrorL2.First() < thresholdPsi, "Error. psiErrorLx not fullfilled");

            double thresholdUr = 2e-3;
            Console.WriteLine("The urErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), urErrorL2.First());
            Console.WriteLine("If the urErrorL2 error is {0} < {1} than good :) ", urErrorL2.First(), thresholdUr);
            Assert.That(urErrorL2.First() < thresholdUr, "Error. urErrorL2 not fullfilled");

            double thresholdUeta = 2e-3;
            Console.WriteLine("The uetaErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uetaErrorL2.First());
            Console.WriteLine("If the uetaErrorL2 error is {0} < {1} than good :) ", uetaErrorL2.First(), thresholdUeta);
            Assert.That(uetaErrorL2.First() < thresholdUeta, "Error. uetaErrorL2 not fullfilled");

            double thresholdUxi = 2e-3; ;
            Console.WriteLine("The uxiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uxiErrorL2.First());
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorL2.First(), thresholdUxi);
            Assert.That(uxiErrorL2.First() < thresholdUxi, "Error. uxiErrorL2 not fullfilled");
        }
        /// <summary>
        /// Transient Tests for the Solution in Dr Dominiks Dierkes Paper (Formular 4.4) for BDF 3 and with R0fix
        /// </summary>
        /// <remarks>
        ///  ATTENTION MUTIPLIER IS ONE BUT should be BSQ (See in the Helical Main). Don't Know why BSQ is not working in this Testcase
        /// </remarks>
        [Test]
        // Convergence for BDF3
        static public void Transient_TimeConv_BDF3_DDD_with_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF3" };
            for (int ell = 0; ell < timeSchemes.Length; ell++) {
                for (int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = Man_Sol_DDD.ManSol_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == true, "R0_fix should be true");
                    Assert.That(Globals.pressureReferencePoint == true, "We have to calculate without PRP, since R0 is on");
                    if (Globals.activeMult == Globals.Multiplier.one && timeConvergence[i].rMin < 10e-6) {
                        Console.WriteLine("Friendly Reminder: Mutiplier One and rMin<10e-6");
                    }
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if (timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if (timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }

                using (var gp = new Gnuplot()) {
                    gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                    gp.OutputFile = $"TemporalConvergence_with_R0fix-BDF{order}.png";
                    gp.SetTitle($"TemporalConvergence_with_R0fix-BDF{order}");
                    gp.PlotLogXLogY(timeStepSize, urErrorL2, title: "ur", format: new PlotFormat("-sr"));
                    gp.PlotLogXLogY(timeStepSize, uetaErrorL2, title: "ueta", format: new PlotFormat("-^b"));
                    gp.PlotLogXLogY(timeStepSize, uxiErrorL2, title: "uxi", format: new PlotFormat("-xm"));
                    gp.PlotLogXLogY(timeStepSize, psiErrorL2, title: "psi", format: new PlotFormat("-*k"));
                    gp.RunAndExit();
                    //Console.WriteLine("endless loop");
                    //while(true) ;
                }

                (double slope, double intercept) regressionUr = DoubleExtensions.LogLogRegression(timeStepSize, urErrorL2);
                (double slope, double intercept) regressionUeta = DoubleExtensions.LogLogRegression(timeStepSize, uetaErrorL2);
                (double slope, double intercept) regressionUxi = DoubleExtensions.LogLogRegression(timeStepSize, uxiErrorL2);
                (double slope, double intercept) regressionPsi = DoubleExtensions.LogLogRegression(timeStepSize, psiErrorL2);

                Console.WriteLine("Order " + order);
                Console.WriteLine("slopeUR " + regressionUr.slope);
                Console.WriteLine("slopeUETA " + regressionUeta.slope);
                Console.WriteLine("slopeUXI " + regressionUxi.slope);
                Console.WriteLine("slopePSI " + regressionPsi.slope);

                double tolerance = 0.5;
                Console.WriteLine("If rate " + regressionUr.slope + " >= " + order + " minus tolerance " + tolerance + " than good");
                Assert.That(regressionUr.slope >= order - tolerance, String.Format("Convergence rate for Ur lower than expected (was {0} but should be {1} +- {2})", regressionUr.slope, order, tolerance));
                Assert.That(regressionUeta.slope >= order - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, order, tolerance));
                Assert.That(regressionUxi.slope >= order - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, order, tolerance));
                Assert.That(regressionPsi.slope >= order - tolerance, String.Format("Convergence rate for Psi lower than expected (was {0} but should be {1} +- {2})", regressionPsi.slope, order, tolerance));
            }
            Console.WriteLine("Remember: R0_fix is " + timeConvergence.First().R0fixOn);
            Console.WriteLine("Remember: Calculating with PRP is " + Globals.pressureReferencePoint);

            double thresholdPsi = 5e-4;
            Console.WriteLine("The psiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), psiErrorL2.First());
            Console.WriteLine("If the psiErrorL2 error is {0} < {1} than good :) ", psiErrorL2.First(), thresholdPsi);
            Assert.That(psiErrorL2.First() < thresholdPsi, "Error. psiErrorLx not fullfilled");

            double thresholdUr = 2e-5;
            Console.WriteLine("The urErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), urErrorL2.First());
            Console.WriteLine("If the urErrorL2 error is {0} < {1} than good :) ", urErrorL2.First(), thresholdUr);
            Assert.That(urErrorL2.First() < thresholdUr, "Error. urErrorL2 not fullfilled");

            double thresholdUeta = 1e-5;
            Console.WriteLine("The uetaErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uetaErrorL2.First());
            Console.WriteLine("If the uetaErrorL2 error is {0} < {1} than good :) ", uetaErrorL2.First(), thresholdUeta);
            Assert.That(uetaErrorL2.First() < thresholdUeta, "Error. uetaErrorL2 not fullfilled");

            double thresholdUxi = 1e-5;
            Console.WriteLine("The uxiErrorL2 error for {0} timesteps and {1} is = {2}", timeSteps.First(), timeSchemes.Last(), uxiErrorL2.First());
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorL2.First(), thresholdUxi);
            Assert.That(uxiErrorL2.First() < thresholdUxi, "Error. uxiErrorL2 not fullfilled");

        }
        #endregion

        #region Hagen Poiseulle
        /// <summary>
        /// Tests if a steady-state, Stokes Hagen Poiseulle flow (aka. Pipe flow),
        /// is also the solution of the instationary solver.
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void PseudoSteady_HP_Re_10_Stokes_with_R0fix([Values(2, 3, 4)] int pOrder = 4, bool NavierStokes = false) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");


            //ilPSP.Environment.NumThreads = 1;
            var ctrlStat = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 1E50, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 40);

            // Initial Values = 0!
            ctrlStat.InitialValues.Clear();
            ctrlStat.InitialValues_Evaluators.Clear();
            ctrlStat.savetodb = true;
            ctrlStat.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            Guid steadyStateSession;
            double[] SteadyStateSolution;
            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrlStat);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverStat.Control.maxAmpli;

                // Exact Solution
                Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (4 * nu);
                Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => 0;


                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverStat.Pressure.L2Error(exSol_pressure);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error out of range: {uxi_L2:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error out of range: {ueta_L2:0.###e-00} (should be close to 0.0)");
                }
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error out of range: {pressure_L2:0.###e-00} (should be close to 0.0)");

                steadyStateSession = solverStat.CurrentSessionInfo.ID;
                SteadyStateSolution = solverStat.CurrentSolution.ToArray();
            }

            // Restart Solution from Exact Solution!!!!!
            // Now dt = 0.0001

            var ctrlTransient = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 40);
            ctrlTransient.GridFunc = null;
            ctrlTransient.RestartInfo = new Tuple<Guid, TimestepNumber>(steadyStateSession, null); // 2nd arg null -> take last timestep;
            ctrlTransient.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrlTransient.NavierStokes = NavierStokes;

            using (var solverTransient = new HelicalMain()) {
                solverTransient.Init(ctrlTransient);
                solverTransient.RunSolverMode();

                SinglePhaseField exSol_ur = solverTransient.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverTransient.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverTransient.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverTransient.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverTransient.Control.maxAmpli;

                Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverTransient.Control.rMax * solverTransient.Control.rMax - X[0] * X[0])) / (4 * nu);
                Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverTransient.Control.rMax * solverTransient.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => 0;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverTransient.ur.L2Error(exSol_ur);
                double uxi_L2 = solverTransient.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverTransient.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverTransient.Pressure.L2Error(exSol_pressure);

                double solDist = solverTransient.CurrentSolution.MPI_L2Dist(SteadyStateSolution);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine();
                Console.WriteLine($"l2 dist to steady sol: {solDist:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error: {uxi_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error: {ueta_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(solDist, 1.0e-10, $"l2 dist to steady sol: {solDist: 0.###e-00} (should be close to 0.0)");
                }
            }

        }


        /// <summary>
        /// Tests if a steady-state, Full Navier Stokes Hagen Poiseulle flow (aka. Pipe flow),
        /// is also the solution of the instationary solver.
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void PseudoSteady_HP_Re_10_Navier_Stokes_with_R0fix([Values(2, 3, 4)] int pOrder = 4, bool NavierStokes = true) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");


            //ilPSP.Environment.NumThreads = 1;
            var ctrlStat = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 1E50, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 40);

            // Initial Values = 0!
            ctrlStat.InitialValues.Clear();
            ctrlStat.InitialValues_Evaluators.Clear();
            ctrlStat.savetodb = true;
            ctrlStat.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            Guid steadyStateSession;
            double[] SteadyStateSolution;
            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrlStat);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverStat.Control.maxAmpli;

                // Exact Solution
                Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (4 * nu);
                Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => 0;


                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverStat.Pressure.L2Error(exSol_pressure);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error out of range: {uxi_L2:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error out of range: {ueta_L2:0.###e-00} (should be close to 0.0)");
                }
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error out of range: {pressure_L2:0.###e-00} (should be close to 0.0)");

                steadyStateSession = solverStat.CurrentSessionInfo.ID;
                SteadyStateSolution = solverStat.CurrentSolution.ToArray();
            }

            // Restart Solution from Exact Solution!!!!!
            // Now dt = 0.0001

            var ctrlTransient = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 40);
            ctrlTransient.GridFunc = null;
            ctrlTransient.RestartInfo = new Tuple<Guid, TimestepNumber>(steadyStateSession, null); // 2nd arg null -> take last timestep;
            ctrlTransient.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrlTransient.NavierStokes = NavierStokes;

            using (var solverTransient = new HelicalMain()) {
                solverTransient.Init(ctrlTransient);
                solverTransient.RunSolverMode();

                SinglePhaseField exSol_ur = solverTransient.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverTransient.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverTransient.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverTransient.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverTransient.Control.maxAmpli;

                Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverTransient.Control.rMax * solverTransient.Control.rMax - X[0] * X[0])) / (4 * nu);
                Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverTransient.Control.rMax * solverTransient.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => 0;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverTransient.ur.L2Error(exSol_ur);
                double uxi_L2 = solverTransient.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverTransient.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverTransient.Pressure.L2Error(exSol_pressure);

                double solDist = solverTransient.CurrentSolution.MPI_L2Dist(SteadyStateSolution);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine();
                Console.WriteLine($"l2 dist to steady sol: {solDist:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error: {uxi_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error: {ueta_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(solDist, 1.0e-10, $"l2 dist to steady sol: {solDist: 0.###e-00} (should be close to 0.0)");
                }
            }

        }

        /// <summary>
        /// Full Navier Stokes Hagen Poiseulle flow (aka. Pipe flow),
        /// With 10% White Noise over laminar Solutin
        /// Reynolds 10
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_HP_Re_10_White_Noise_10_Procent_with_R0fix([Values(3)] int pOrder = 3) {

            double maxAmpitude = 40;
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            var ctrl = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1000, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, MaxAmp: maxAmpitude);

            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;
            ctrl.savetodb = true;
            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();
            ctrl.ImmediatePlotPeriod = 10;


            var random = 0.1 * maxAmpitude;
            // Initial Values
            // ==============
            string InitialValue_ur_p =
            "static class MyInitialValue_ur_p {" +
            " public static double GenerateRandomValue(double[] X, double t) {" +
            "    var random = new Random();" +
            "    double randomValue = (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
            "   return randomValue;" +
            " }" +
            "}";

            string InitialValue_uxi =
            "static class MyInitialValue_uxi {" +
            " public static double GenerateRandomValue(double[] X, double t) {" +
            "    var random = new Random();" +
            "    double randomValue_and_lami = - " + maxAmpitude + " * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
            "     * (" + (a * a) + " * (" + (ctrl.rMax * ctrl.rMax) + " - X[0] * X[0])) / (4 * " + nu + ")" +
            "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
            "   return randomValue_and_lami;" +
            " }" +
            "}";


              string InitialValue_ueta =
               "static class MyInitialValue_ueta {" +
               " public static double GenerateRandomValue(double[] X, double t) {" +
               "    var random = new Random();" +
               "    double randomValue_and_lami = " + maxAmpitude + " * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
               "     * (" + (a * b) + " * (" + (ctrl.rMax * ctrl.rMax) + " - X[0] * X[0])) / (X[0] * 4 * " + nu + ")" +
               "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
               "   return randomValue_and_lami;" +
               " }" +
               "}";


            var ur0_p0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur_p.GenerateRandomValue", true, InitialValue_ur_p);
            var uxi0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var ueta0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            ctrl.AddInitialValue("Pressure", ur0_p0);
            ctrl.AddInitialValue("ur", ur0_p0);
            ctrl.AddInitialValue("ueta", ueta0);
            ctrl.AddInitialValue("uxi", uxi0);

            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrl);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();
                SinglePhaseField pressureError = solverStat.Pressure.CloneAs();

                // Exact Solution
                Func<double[], double> uxi = (X) => -maxAmpitude * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (4 * nu);
                Func<double[], double> ueta = (X) => maxAmpitude * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => 0;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                pressureError.ProjectField(pressure);
                pressureError.Acc(-1, solverStat.Pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);

                double psiErrorMean = pressureError.GetMeanValueTotal(null);
                pressureError.AccConstant(-psiErrorMean);

                double pressure_L2 = pressureError.L2Norm();
                Console.WriteLine($"ur       L2 Error/maxAmpitude: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error/maxAmpitude: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error/maxAmpitude: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error/maxAmpitude: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2 / maxAmpitude, 1.0e-10, $"ur L2 Error/maxAmpitude out of range: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(uxi_L2 / maxAmpitude, 1.0e-10, $"uxi L2 Error/maxAmpitude out of range: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(ueta_L2 / maxAmpitude, 1.0e-10, $"ueta L2 Error/maxAmpitude out of range: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2 / maxAmpitude, 1.0e-8, $"pressure L2 Error/maxAmpitude out of range: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
            }

        }

        /// <summary>
        /// Full Navier Stokes Hagen Poiseulle flow (aka. Pipe flow),
        /// With 10% White Noise over laminar Solutin
        /// Reynolds 10000
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_HP_Re_2500_White_Noise_10_Procent_with_R0fix([Values(3)] int pOrder = 3) {

            double maxAmpitude = 10000;
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            var ctrl = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 10000, deltaT: 0.000005, _DbPath: tempDB.Path, bdfOrder: 1, MaxAmp: maxAmpitude);

            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;
            ctrl.savetodb = true;
            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();
            ctrl.ImmediatePlotPeriod = 10;


            var random = 0.1 * maxAmpitude;
            // Initial Values
            // ==============
            string InitialValue_ur_p =
            "static class MyInitialValue_ur_p {" +
            " public static double GenerateRandomValue(double[] X, double t) {" +
            "    var random = new Random();" +
            "    double randomValue = (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
            "   return randomValue;" +
            " }" +
            "}";
            string InitialValue_uxi =
            "static class MyInitialValue_uxi {" +
            " public static double GenerateRandomValue(double[] X, double t) {" +
            "    var random = new Random();" +
            "    double randomValue_and_lami = - " + maxAmpitude + " * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
            "     * (" + (a * a) + " * (" + (ctrl.rMax * ctrl.rMax) + " - X[0] * X[0])) / (4 * " + nu + ")" +
            "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
            "   return randomValue_and_lami;" +
            " }" +
            "}";
            string InitialValue_ueta =
             "static class MyInitialValue_ueta {" +
             " public static double GenerateRandomValue(double[] X, double t) {" +
             "    var random = new Random();" +
             "    double randomValue_and_lami = " + maxAmpitude + " * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
             "     * (" + (a * b) + " * (" + (ctrl.rMax * ctrl.rMax) + " - X[0] * X[0])) / (X[0] * 4 * " + nu + ")" +
             "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
             "   return randomValue_and_lami;" +
             " }" +
             "}";


            var ur0_p0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur_p.GenerateRandomValue", true, InitialValue_ur_p);
            var uxi0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var ueta0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            ctrl.AddInitialValue("Pressure", ur0_p0);
            ctrl.AddInitialValue("ur", ur0_p0);
            ctrl.AddInitialValue("ueta", ueta0);
            ctrl.AddInitialValue("uxi", uxi0);

            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrl);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();
                SinglePhaseField pressureError = solverStat.Pressure.CloneAs();

                // Exact Solution
                Func<double[], double> uxi = (X) => -maxAmpitude * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (4 * nu);
                Func<double[], double> ueta = (X) => maxAmpitude * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (solverStat.Control.rMax * solverStat.Control.rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => 0;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                pressureError.ProjectField(pressure);
                pressureError.Acc(-1, solverStat.Pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);

                double psiErrorMean = pressureError.GetMeanValueTotal(null);
                pressureError.AccConstant(-psiErrorMean);

                double pressure_L2 = pressureError.L2Norm();
                Console.WriteLine($"ur       L2 Error/maxAmpitude: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error/maxAmpitude: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error/maxAmpitude: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error/maxAmpitude: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2 / maxAmpitude, 1.0e-10, $"ur L2 Error/maxAmpitude out of range: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(uxi_L2 / maxAmpitude, 1.0e-10, $"uxi L2 Error/maxAmpitude out of range: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(ueta_L2 / maxAmpitude, 1.0e-10, $"ueta L2 Error/maxAmpitude out of range: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2 / maxAmpitude, 1.0e-8, $"pressure L2 Error/maxAmpitude out of range: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
            }

        }
        #endregion

        #region Centrifugal Flow
        /// <summary>
        /// Tests if the exact solution for the, Stokes Centrifugal flow (aka. Centrifugal flow),
        /// is also the solution of the instationary solver. Solver should find the exact Solution. 
        /// For BDF1
        /// </summary>
        /// <remarks>
        /// In the Stokes Case for the centrifugal flow the Pressure is Zero!
        /// </remarks>
        /// 
        [Test]
        static public void PseudoSteady_CF_Re_10_Stokes_with_R0fix([Values(2, 3, 4)] int pOrder = 4, bool NavierStokes = false) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            //ilPSP.Environment.NumThreads = 1;
            var ctrlStat = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 1E50, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 10);
            // Initial Values = 0!
            ctrlStat.InitialValues.Clear();
            ctrlStat.InitialValues_Evaluators.Clear();
            ctrlStat.savetodb = true;
            ctrlStat.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            Guid steadyStateSession;
            double[] SteadyStateSolution;
            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrlStat);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverStat.Control.maxAmpli;

                // Exact Solution
                Func<double[], double> ur = (X) => 0;
                //Func<double[], double> pressure = (X) => MaxAmp * MaxAmp * X[0] * X[0] * 0.5;
                Func<double[], double> pressure = (X) => 0;
                Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;


                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverStat.Pressure.L2Error(exSol_pressure);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error out of range: {uxi_L2:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error out of range: {ueta_L2:0.###e-00} (should be close to 0.0)");
                }
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error out of range: {pressure_L2:0.###e-00} (should be close to 0.0)");

                steadyStateSession = solverStat.CurrentSessionInfo.ID;
                SteadyStateSolution = solverStat.CurrentSolution.ToArray();
            }

            // Restart Solution from Exact Solution!!!!!
            // Now dt = 0.0001

            var ctrlTransient = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 10);
            ctrlTransient.GridFunc = null;
            ctrlTransient.RestartInfo = new Tuple<Guid, TimestepNumber>(steadyStateSession, null); // 2nd arg null -> take last timestep;
            ctrlTransient.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrlTransient.NavierStokes = NavierStokes;

            using (var solverTransient = new HelicalMain()) {
                solverTransient.Init(ctrlTransient);
                solverTransient.RunSolverMode();

                SinglePhaseField exSol_ur = solverTransient.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverTransient.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverTransient.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverTransient.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverTransient.Control.maxAmpli;

                // Exact Solution
                Func<double[], double> ur = (X) => 0;
                //Func<double[], double> pressure = (X) => MaxAmp * MaxAmp * X[0] * X[0] * 0.5;
                Func<double[], double> pressure = (X) => 0;
                Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverTransient.ur.L2Error(exSol_ur);
                double uxi_L2 = solverTransient.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverTransient.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverTransient.Pressure.L2Error(exSol_pressure);

                double solDist = solverTransient.CurrentSolution.MPI_L2Dist(SteadyStateSolution);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine();
                Console.WriteLine($"l2 dist to steady sol: {solDist:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error: {uxi_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error: {ueta_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(solDist, 1.0e-10, $"l2 dist to steady sol: {solDist: 0.###e-00} (should be close to 0.0)");
                }
            }

        }

        /// <summary>
        /// Tests if the exact solution for the, Stokes Centrifugal flow (aka. Centrifugal flow),
        /// is also the solution of the instationary solver. Solver should find the exact Solution. 
        /// For BDF1
        /// </summary>
        /// <remarks>
        /// In the Stokes Case for the centrifugal flow the Pressure is Zero!
        /// </remarks>
        /// 
        [Test]
        static public void PseudoSteady_CF_Re_10_Navier_Stokes_with_R0fix([Values(2, 3, 4)] int pOrder = 4, bool NavierStokes = true) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            //ilPSP.Environment.NumThreads = 1;
            var ctrlStat = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 1E50, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 10);
            // clear Inital Values
            ctrlStat.InitialValues.Clear();
            ctrlStat.InitialValues_Evaluators.Clear();
            // Set Initial Values to 1
            ctrlStat.AddInitialValue("Pressure", new Formula("(X) =>1"));
            ctrlStat.AddInitialValue("ur", new Formula("(X) => 1"));
            ctrlStat.AddInitialValue("ueta", new Formula($"(X) => 1"));
            ctrlStat.AddInitialValue("uxi", new Formula($"(X) => 1"));
            ctrlStat.savetodb = true;
            ctrlStat.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
            ctrlStat.ImmediatePlotPeriod = 1;
            Guid steadyStateSession;
            double[] SteadyStateSolution;
            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrlStat);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverStat.Control.maxAmpli;

                // Exact Solution
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => MaxAmp * MaxAmp * X[0] * X[0] * 0.5;
                Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;


                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverStat.Pressure.L2Error(exSol_pressure);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error out of range: {uxi_L2:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error out of range: {ueta_L2:0.###e-00} (should be close to 0.0)");
                }
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error out of range: {pressure_L2:0.###e-00} (should be close to 0.0)");

                steadyStateSession = solverStat.CurrentSessionInfo.ID;
                SteadyStateSolution = solverStat.CurrentSolution.ToArray();
            }

            // Restart Solution from Exact Solution!!!!!
            // Now dt = 0.0001

            var ctrlTransient = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, rMin: 0, MaxAmp: 10);
            ctrlTransient.GridFunc = null;
            ctrlTransient.RestartInfo = new Tuple<Guid, TimestepNumber>(steadyStateSession, null); // 2nd arg null -> take last timestep;
            ctrlTransient.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrlTransient.NavierStokes = NavierStokes;
            ctrlTransient.ImmediatePlotPeriod = 1;
            using (var solverTransient = new HelicalMain()) {
                solverTransient.Init(ctrlTransient);
                solverTransient.RunSolverMode();

                SinglePhaseField exSol_ur = solverTransient.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverTransient.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverTransient.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverTransient.Pressure.CloneAs();

                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;
                double MaxAmp = solverTransient.Control.maxAmpli;

                // Exact Solution
                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => MaxAmp * MaxAmp * X[0] * X[0] * 0.5;
                Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                exSol_pressure.ProjectField(pressure);

                double ur_L2 = solverTransient.ur.L2Error(exSol_ur);
                double uxi_L2 = solverTransient.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverTransient.ueta.L2Error(exSol_ueta);
                double pressure_L2 = solverTransient.Pressure.L2Error(exSol_pressure);

                double solDist = solverTransient.CurrentSolution.MPI_L2Dist(SteadyStateSolution);

                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine();
                Console.WriteLine($"l2 dist to steady sol: {solDist:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");
                if (pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error: {uxi_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error: {ueta_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(solDist, 1.0e-10, $"l2 dist to steady sol: {solDist: 0.###e-00} (should be close to 0.0)");
                }
            }

        }

        /// <summary>
        /// Full Navier Stokes Cylindrical flow 
        /// With 10% White Noise over laminar Solutin
        /// Reynolds 10
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_CF_Re_10_White_Noise_10_Procent_with_R0fix([Values(3)] int pOrder = 3) {

            double maxAmpitude = 10;
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            var ctrl = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1000, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, MaxAmp: maxAmpitude);

            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;
            ctrl.savetodb = true;
            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();
            ctrl.ImmediatePlotPeriod = 10;


            var random = 0.1 * maxAmpitude;
            // Initial Values
            // ==============
            string InitialValue_p =
                "static class MyInitialValue_p {" +
                " public static double GenerateRandomValue(double[] X, double t) {" +
                "    var random = new Random();" +
                "    double randomValue = " + (maxAmpitude * maxAmpitude * 0.5) + " * X[0] * X[0] " +
                "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
                "   return randomValue;" +
                " }" +
                "}";

            string InitialValue_ur =
            "static class MyInitialValue_ur {" +
            " public static double GenerateRandomValue(double[] X, double t) {" +
            "    var random = new Random();" +
            "    double randomValue = (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
            "   return randomValue;" +
            " }" +
            "}";

            string InitialValue_uxi =
                "static class MyInitialValue_uxi {" +
                " public static double GenerateRandomValue(double[] X, double t) {" +
                "    var random = new Random();" +
                "    double randomValue_and_lami = " + (maxAmpitude * b) + " * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
                "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
                "   return randomValue_and_lami;" +
                " }" +
                "}";

            string InitialValue_ueta =
                "static class MyInitialValue_ueta {" +
                " public static double GenerateRandomValue(double[] X, double t) {" +
                "    var random = new Random();" +
                "    double randomValue_and_lami = " + (maxAmpitude * a) + " * X[0] * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
                "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
                "   return randomValue_and_lami;" +
                " }" +
                "}";

            var p0 = new BoSSS.Solution.Control.Formula("MyInitialValue_p.GenerateRandomValue", true, InitialValue_p);
            var ur0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur.GenerateRandomValue", true, InitialValue_ur);
            var uxi0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var ueta0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            ctrl.AddInitialValue("Pressure", p0);
            ctrl.AddInitialValue("ur", ur0);
            ctrl.AddInitialValue("ueta", ueta0);
            ctrl.AddInitialValue("uxi", uxi0);

            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrl);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();
                SinglePhaseField pressureError = solverStat.Pressure.CloneAs();

                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => maxAmpitude * maxAmpitude * X[0] * X[0] * 0.5;
                Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * maxAmpitude * X[0];
                Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * maxAmpitude;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                pressureError.ProjectField(pressure);
                pressureError.Acc(-1, solverStat.Pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);

                double psiErrorMean = pressureError.GetMeanValueTotal(null);
                pressureError.AccConstant(-psiErrorMean);

                double pressure_L2 = pressureError.L2Norm();
                Console.WriteLine($"ur       L2 Error/maxAmpitude: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error/maxAmpitude: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error/maxAmpitude: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error/maxAmpitude: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2 / maxAmpitude, 1.0e-10, $"ur L2 Error/maxAmpitude out of range: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(uxi_L2 / maxAmpitude, 1.0e-10, $"uxi L2 Error/maxAmpitude out of range: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(ueta_L2 / maxAmpitude, 1.0e-10, $"ueta L2 Error/maxAmpitude out of range: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2 / maxAmpitude, 1.0e-8, $"pressure L2 Error/maxAmpitude out of range: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
            }

        }

        /// <summary>
        /// Full Navier Stokes Cylindrical flow 
        /// With 10% White Noise over laminar Solutin
        /// Reynolds 2500
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void Transient_CF_Re_2500_White_Noise_10_Procent_with_R0fix([Values(3)] int pOrder = 3) {

            double maxAmpitude = 2500;
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");
            var ctrl = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, numOfTimesteps: 1000, deltaT: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1, MaxAmp: maxAmpitude);

            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;
            ctrl.savetodb = true;
            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();
            ctrl.ImmediatePlotPeriod = 1;


            var random = 0.1 * maxAmpitude;
            // Initial Values
            // ==============
            string InitialValue_p =
                "static class MyInitialValue_p {" +
                " public static double GenerateRandomValue(double[] X, double t) {" +
                "    var random = new Random();" +
                "    double randomValue = " + (maxAmpitude * maxAmpitude * 0.5) + " * X[0] * X[0] " +
                "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
                "   return randomValue;" +
                " }" +
                "}";

            string InitialValue_ur =
            "static class MyInitialValue_ur {" +
            " public static double GenerateRandomValue(double[] X, double t) {" +
            "    var random = new Random();" +
            "    double randomValue = (random.NextDouble() - 0.5) * Math.Sin(X[0]) * Math.Sin(X[1] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
            "   return randomValue;" +
            " }" +
            "}";

            string InitialValue_uxi =
                "static class MyInitialValue_uxi {" +
                " public static double GenerateRandomValue(double[] X, double t) {" +
                "    var random = new Random();" +
                "    double randomValue_and_lami = " + (maxAmpitude * b) + " * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
                "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
                "   return randomValue_and_lami;" +
                " }" +
                "}";

            string InitialValue_ueta =
                "static class MyInitialValue_ueta {" +
                " public static double GenerateRandomValue(double[] X, double t) {" +
                "    var random = new Random();" +
                "    double randomValue_and_lami = " + (maxAmpitude * a) + " * X[0] * (X[0] / (Math.Sqrt(" + (a * a) + " * X[0] * X[0] + " + (b * b) + ")))" +
                "     + (random.NextDouble() - 0.5) * Math.Sin(X[1]) * Math.Sin(X[0] * Math.PI * 0.5 / " + ctrl.rMax + ");" +
                "   return randomValue_and_lami;" +
                " }" +
                "}";



            var p0 = new BoSSS.Solution.Control.Formula("MyInitialValue_p.GenerateRandomValue", true, InitialValue_p);
            var ur0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ur.GenerateRandomValue", true, InitialValue_ur);
            var uxi0 = new BoSSS.Solution.Control.Formula("MyInitialValue_uxi.GenerateRandomValue", true, InitialValue_uxi);
            var ueta0 = new BoSSS.Solution.Control.Formula("MyInitialValue_ueta.GenerateRandomValue", true, InitialValue_ueta);


            ctrl.AddInitialValue("Pressure", p0);
            ctrl.AddInitialValue("ur", ur0);
            ctrl.AddInitialValue("ueta", ueta0);
            ctrl.AddInitialValue("uxi", uxi0);

            using (var solverStat = new HelicalMain()) {
                solverStat.Init(ctrl);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();
                SinglePhaseField pressureError = solverStat.Pressure.CloneAs();

                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => maxAmpitude * maxAmpitude * X[0] * X[0] * 0.5;
                Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * maxAmpitude * X[0];
                Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * maxAmpitude;

                exSol_ur.ProjectField(ur);
                exSol_ueta.ProjectField(ueta);
                exSol_uxi.ProjectField(uxi);
                pressureError.ProjectField(pressure);
                pressureError.Acc(-1, solverStat.Pressure);

                double ur_L2 = solverStat.ur.L2Error(exSol_ur);
                double uxi_L2 = solverStat.uxi.L2Error(exSol_uxi);
                double ueta_L2 = solverStat.ueta.L2Error(exSol_ueta);

                double psiErrorMean = pressureError.GetMeanValueTotal(null);
                pressureError.AccConstant(-psiErrorMean);

                double pressure_L2 = pressureError.L2Norm();
                Console.WriteLine($"ur       L2 Error/maxAmpitude: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error/maxAmpitude: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error/maxAmpitude: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error/maxAmpitude: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2 / maxAmpitude, 1.0e-10, $"ur L2 Error/maxAmpitude out of range: {ur_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(uxi_L2 / maxAmpitude, 1.0e-10, $"uxi L2 Error/maxAmpitude out of range: {uxi_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(ueta_L2 / maxAmpitude, 1.0e-10, $"ueta L2 Error/maxAmpitude out of range: {ueta_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2 / maxAmpitude, 1.0e-8, $"pressure L2 Error/maxAmpitude out of range: {pressure_L2 / maxAmpitude:0.###e-00} (should be close to 0.0)");
            }

        }
        #endregion

    }
}
