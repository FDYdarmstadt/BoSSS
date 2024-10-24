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


        [Test]
        // Convergence for BDF1
        static public void TimeConvergenceBDF1_no_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF1"};
            for(int ell = 0; ell < timeSchemes.Length; ell++) {
                for(int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = StokesHelical_Ak.HardcodedControl.ManSol_Transient_DDD_Paper(noOfCellsR:32, noOfCellsXi:32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree:5, rMin:0.1);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == false, "R0_fix should be false");
                    Assert.That(timeConvergence[i].PressureReferencePoint == true, "We have to calculate with PRP, since R0 is off");
                    Assert.That(Globals.activeMult == Globals.Multiplier.one);
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if(timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if(timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }


                using(var gp = new Gnuplot()) {
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
            Console.WriteLine("Remember: Calculating with PRP is " + timeConvergence.First().PressureReferencePoint);


            double thresholdPsi = 5e-2;
            Console.WriteLine("The psiErrorL2 error for {0} timesteps and {1} is = {2}" , timeSteps.First(), timeSchemes.Last(), psiErrorL2.First()) ;
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


        [Test]
        // Convergence for BDF3
        static public void TimeConvergenceBDF3_no_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] {"BDF3" };
            for(int ell = 0; ell < timeSchemes.Length; ell++) {
                for(int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = StokesHelical_Ak.HardcodedControl.ManSol_Transient_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0.1);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == false, "R0_fix should be false");
                    Assert.That(timeConvergence[i].PressureReferencePoint == true, "We have to calculate with PRP, since R0 is off");
                    Assert.That(Globals.activeMult == Globals.Multiplier.one);
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if(timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if(timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }

                using(var gp = new Gnuplot()) {
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
            Console.WriteLine("Remember: Calculating with PRP is " + timeConvergence.First().PressureReferencePoint);

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

        [Test]
        // Convergence for BDF1
        static public void TimeConvergenceBDF1_with_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF1"};
            for(int ell = 0; ell < timeSchemes.Length; ell++) {
                for(int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = StokesHelical_Ak.HardcodedControl.ManSol_Transient_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == true, "R0_fix should be true");
                    Assert.That(timeConvergence[i].PressureReferencePoint == true, "We have to calculate with PRP");
                    if(Globals.activeMult == Globals.Multiplier.one && timeConvergence[i].rMin < 10e-6) {
                        Console.WriteLine("Friendly Reminder: Mutiplier One and rMin<10e-6");
                    }
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if(timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if(timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }

                using(var gp = new Gnuplot()) {
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
            Console.WriteLine("Remember: Calculating with PRP is " + timeConvergence.First().PressureReferencePoint);

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

        [Test]
        // Convergence for BDF3
        static public void TimeConvergenceBDF3_with_R0fix() {

            int[] timeSteps = new int[] { 64, 32, 16, 8, 4 };
            double[] urErrorL2 = new double[timeSteps.Length];
            double[] uetaErrorL2 = new double[timeSteps.Length];
            double[] uxiErrorL2 = new double[timeSteps.Length];
            double[] psiErrorL2 = new double[timeSteps.Length];
            double[] timeStepSize = new double[timeSteps.Length];
            HelicalControl[] timeConvergence = new HelicalControl[timeSteps.Length];
            string[] timeSchemes = new string[] { "BDF3" };
            for(int ell = 0; ell < timeSchemes.Length; ell++) {
                for(int i = 0; i < timeSteps.Length; i++) {
                    timeConvergence[i] = StokesHelical_Ak.HardcodedControl.ManSol_Transient_DDD_Paper(noOfCellsR: 32, noOfCellsXi: 32, dtRefining: timeSteps[i], bdfOrder: timeSchemes[ell], degree: 5, rMin: 0);
                    var solver = new HelicalMain();
                    solver.Init(timeConvergence[i]);
                    solver.RunSolverMode();
                    Assert.That(timeConvergence[i].R0fixOn == true, "R0_fix should be true");
                    Assert.That(timeConvergence[i].PressureReferencePoint == true, "We have to calculate without PRP, since R0 is on");
                    if(Globals.activeMult == Globals.Multiplier.one && timeConvergence[i].rMin < 10e-6) {
                        Console.WriteLine("Friendly Reminder: Mutiplier One and rMin<10e-6");
                    }
                    psiErrorL2[i] = solver.psiErrorL2;
                    urErrorL2[i] = solver.urErrorL2;
                    uetaErrorL2[i] = solver.uetaErrorL2;
                    uxiErrorL2[i] = solver.uxiErrorL2;
                    timeStepSize[i] = 2 * Math.PI / timeSteps[i];
                }
                double order;

                if(timeSchemes[ell] == "BDF3") {
                    order = 3;
                } else if(timeSchemes[ell] == "BDF1") {
                    order = 1;
                } else {
                    throw new ArgumentException("Unsupported BDF scheme: " + timeSchemes[ell]);
                }

                using(var gp = new Gnuplot()) {
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
            Console.WriteLine("Remember: Calculating with PRP is " + timeConvergence.First().PressureReferencePoint);

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
        /// Tests if a steady-state, laminar Hagen Poiseulle flow (aka. Pipe flow),
        /// is also the solution of the instationary solver.
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void PseudoSteadyHagenPoiseulle(
            [Values(2, 3, 4)] int pOrder = 4,
            [Values(false, true)] bool NavierStokes = false
            ) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");


            //ilPSP.Environment.NumThreads = 1;
            var ctrlStat = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, dtRefining: 1, Tend: 1E50, _DbPath:tempDB.Path, bdfOrder:1);

            ctrlStat.InitialValues.Clear();
            ctrlStat.InitialValues_Evaluators.Clear();
            ctrlStat.savetodb = true;
            ctrlStat.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            Guid steadyStateSession;
            double[] SteadyStateSolution;
            using(var solverStat = new HelicalMain()) {
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
                if(pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error out of range: {uxi_L2:0.###e-00} (should be close to 0.0)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error out of range: {ueta_L2:0.###e-00} (should be close to 0.0)");
                }
                Assert.LessOrEqual(pressure_L2, 1.0e-10, $"pressure L2 Error out of range: {pressure_L2:0.###e-00} (should be close to 0.0)");

                steadyStateSession = solverStat.CurrentSessionInfo.ID;
                SteadyStateSolution = solverStat.CurrentSolution.ToArray();
            }

            //SplittingTimestepper.Uninfty = SteadyStateSolution.CloneAs();

            var ctrlTransient = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, dtRefining: 1, Tend: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1);
            ctrlTransient.GridFunc = null;
            ctrlTransient.RestartInfo = new Tuple<Guid, TimestepNumber>(steadyStateSession, null); // 2nd arg null -> take last timestep;
            ctrlTransient.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrlTransient.NavierStokes = NavierStokes;

            using(var solverTransient = new HelicalMain()) {
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
                if(pOrder == 4) {
                    Assert.LessOrEqual(uxi_L2, 1.0e-9, $"uxi L2 Error: {uxi_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(ueta_L2, 1.0e-9, $"ueta L2 Error: {ueta_L2:0.###e-00} (should nonzero)");
                    Assert.LessOrEqual(solDist, 1.0e-10, $"l2 dist to steady sol: {solDist: 0.###e-00} (should be close to 0.0)");
                }
            }

        }

        /// <summary>
        /// Tests if the exact solution for the, Centrifugal flow (aka. Centrifugal flow),
        /// is also the solution of the instationary solver. Solver should find the exact Solution. 
        /// For BDF1
        /// </summary>
        /// <remarks>
        /// </remarks>
        [Test]
        static public void PseudoSteadyCentrifuge(
            [Values(4)] int pOrder = 4
            ) {
            var tempDB = DatabaseInfo.CreateOrOpen("tempDB");

            var ctrlStat = StokesHelical_Ak.DNS_Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: 64, noOfCellsXi: 64, dtRefining: 1, Tend: 0.0001, _DbPath: tempDB.Path, bdfOrder: 1);

            double a = Globals.a;
            double b = Globals.b;
            double MaxAmp = ctrlStat.maxAmpli;

            ctrlStat.savetodb = true;
            ctrlStat.NoOfTimesteps = 1;
            ctrlStat.ImmediatePlotPeriod = 1;
            // Clear Initial Values
            ctrlStat.InitialValues.Clear();
            ctrlStat.InitialValues_Evaluators.Clear();
            
            // Exact Solution
            ctrlStat.AddInitialValue("Pressure", new Formula($"(X) => {MaxAmp} * {MaxAmp} * X[0] * X[0] *0.5"));
            ctrlStat.AddInitialValue("ur", new Formula($"(X) => 0"));
            ctrlStat.AddInitialValue("ueta", new Formula($"(X) => (X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} )))*{a}*{MaxAmp}* X[0]"));
            ctrlStat.AddInitialValue("uxi", new Formula($"(X) => (X[0]/(Math.Sqrt({a * a} * X[0] * X[0] + {b * b} )))*{b}*{MaxAmp}"));
            using(var solverStat = new HelicalMain()) {
                solverStat.Init(ctrlStat);
                solverStat.RunSolverMode();

                SinglePhaseField exSol_ur = solverStat.ur.CloneAs();
                SinglePhaseField exSol_ueta = solverStat.ueta.CloneAs();
                SinglePhaseField exSol_uxi = solverStat.uxi.CloneAs();
                SinglePhaseField exSol_pressure = solverStat.Pressure.CloneAs();
                SinglePhaseField pressureError = solverStat.Pressure.CloneAs();

                Func<double[], double> ur = (X) => 0;
                Func<double[], double> pressure = (X) => MaxAmp * MaxAmp * X[0] * X[0] *0.5;
                Func<double[], double> ueta = (X) => (X[0]/(Math.Sqrt(a * a * X[0] * X[0] + b * b )))*a*MaxAmp* X[0];
                Func<double[], double> uxi = (X) => (X[0]/(Math.Sqrt(a * a * X[0] * X[0] + b * b )))*b*MaxAmp;

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
                Console.WriteLine($"ur       L2 Error: {ur_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"uxi      L2 Error: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"ueta     L2 Error: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Console.WriteLine($"pressure L2 Error: {pressure_L2:0.###e-00} (should be close to 0.0)");

                Assert.LessOrEqual(ur_L2, 1.0e-10, $"ur L2 Error out of range: {ur_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(uxi_L2, 1.0e-10, $"uxi L2 Error out of range: {uxi_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(ueta_L2, 1.0e-10, $"ueta L2 Error out of range: {ueta_L2:0.###e-00} (should be close to 0.0)");
                Assert.LessOrEqual(pressure_L2, 1.0e-8, $"pressure L2 Error out of range: {pressure_L2:0.###e-00} (should be close to 0.0)");
            }

        }
    }
}
