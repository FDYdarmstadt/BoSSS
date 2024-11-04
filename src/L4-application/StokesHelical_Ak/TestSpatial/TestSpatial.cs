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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
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
using StokesHelical_Ak.Hard_Coded_Control;

namespace StokesHelical_Ak.TestSpartial
{

    [TestFixture]
    static public class TestSpatial {

        /// Hagen Poiseulle Flow.
        /// See PhD-Draft from Schahin Akbari
        #region HagenPoiseulle 

        #region With R0fix
        /// <summary>
        /// Tests a steady-state Convergence Hagen Poiseulle flow Re=10^5 (aka. Pipe flow).
        /// The Reynolds Number is set to 100.000. The amplitude is 400.000 since Re = Amplitude/4
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
        static public void Steady_SpatialConv_HP_Re_100000_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            int[] amplitude = new int[] { 400000 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] ur_L2_Error = new double[gridSize.Length];
            double[] ueta_L2_Error = new double[gridSize.Length];
            double[] uxi_L2_Error = new double[gridSize.Length];
            double[] pressure_L2_Error = new double[gridSize.Length];

            HelicalControl[,] spaceConvergence_HP = new HelicalControl[gridSize.Length, amplitude.Length];
            SinglePhaseField[] exSol_ur = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_ueta = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_uxi = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_pressure = new SinglePhaseField[gridSize.Length];

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            for (int ell = 0; ell < r_min.Length; ell++) {

                //ilPSP.Environment.NumThreads = 1;
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence_HP[i, k] = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: amplitude[k]);
                    }
                }
                // Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence_HP.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                        spaceConvergence_HP[i, k].InitialValues.Clear();
                        spaceConvergence_HP[i, k].InitialValues_Evaluators.Clear();
                        spaceConvergence_HP[i, k].savetodb = false;
                        spaceConvergence_HP[i, k].TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                        // Exakte Loesung
                        double MaxAmp = spaceConvergence_HP[i, k].maxAmpli;
                        Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (spaceConvergence_HP[i, k].rMax * spaceConvergence_HP[i, k].rMax - X[0] * X[0])) / (4 * nu);
                        Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (spaceConvergence_HP[i, k].rMax * spaceConvergence_HP[i, k].rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                        Func<double[], double> ur = (X) => 0;
                        Func<double[], double> pressure = (X) => 0;

                        using (var solver = new HelicalMain()) {
                            solver.Init(spaceConvergence_HP[i, k]);
                            solver.RunSolverMode();
                            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                            Assert.IsTrue(spaceConvergence_HP[i, k].R0fixOn, "R0fix must be turned on");
                            Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, "Pressure Reference Point has to be true");
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


                            Console.WriteLine($"ur       L2 Error/{amplitude[k]}: {ur_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                            Console.WriteLine($"uxi      L2 Error/{amplitude[k]}: {uxi_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"ueta     L2 Error/{amplitude[k]}: {ueta_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"pressure L2 Error/{amplitude[k]}: {pressure_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");


                            Assert.LessOrEqual(ur_L2_Error[i] / amplitude[k], 1.0e-11, $"ur L2 Error / amplitude out of range: {ur_L2_Error[i]:0.###e-00} (should be close to 0.0 )");
                            Assert.LessOrEqual(pressure_L2_Error[i] / amplitude[k], 1.0e-10, $"pressure L2 Error / amplitude  out of range: {pressure_L2_Error[i]:0.###e-00} (should be close to 0.0)");
                        }
                    }
                    using (var gp = new Gnuplot()) {
                        gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                        gp.OutputFile = $"SpatialConvergence Amplitude {amplitude[k]} with R0fix p{pOrder}.png";
                        gp.SetTitle($"Spatial Convergence Hagen Poiseulle Amplitude {amplitude[k]} with R0fix p{pOrder}");
                        gp.PlotLogXLogY(h, ur_L2_Error, title: $"ur L2error", format: new PlotFormat("-sr"));
                        gp.PlotLogXLogY(h, ueta_L2_Error, title: $"ueta L2error", format: new PlotFormat("-^b"));
                        gp.PlotLogXLogY(h, uxi_L2_Error, title: $"uxi L2error", format: new PlotFormat("-xm"));
                        gp.PlotLogXLogY(h, pressure_L2_Error, title: $"psi L2error", format: new PlotFormat("-*k"));
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
                    if (pOrder != 5) {
                        Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                        Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                    }
                    if (pOrder == 5) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // error thresholds are specified for the SECOND!!!!!! finest grid with pOrder == 5 
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        Assert.LessOrEqual(ur_L2_Error[1] / amplitude[k], 1.0e-10, $"ur L2 Error / amplitude out of range: {ur_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                        Assert.LessOrEqual(uxi_L2_Error[1] / amplitude[k], 1.0e-10, $"uxi L2 Error / amplitude out of range: {uxi_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(ueta_L2_Error[1] / amplitude[k], 1.0e-10, $"ueta L2 Error / amplitude out of range: {ueta_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(pressure_L2_Error[1] / amplitude[k], 1.0e-10, $"pressure L2 Error / amplitude out of range: {pressure_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                    }
                }
            }
        }

        /// <summary>
        /// Tests a steady-state Convergence Hagen Poiseulle flow Re=10 (aka. Pipe flow).
        /// The Reynolds Number is set to 10 The amplitude is 40 since Re = Amplitude/4
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
        static public void Steady_SpatialConv_HP_Re_10_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            int[] amplitude = new int[] { 40 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] ur_L2_Error = new double[gridSize.Length];
            double[] ueta_L2_Error = new double[gridSize.Length];
            double[] uxi_L2_Error = new double[gridSize.Length];
            double[] pressure_L2_Error = new double[gridSize.Length];

            HelicalControl[,] spaceConvergence_HP = new HelicalControl[gridSize.Length, amplitude.Length];
            SinglePhaseField[] exSol_ur = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_ueta = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_uxi = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_pressure = new SinglePhaseField[gridSize.Length];

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            for (int ell = 0; ell < r_min.Length; ell++) {

                //ilPSP.Environment.NumThreads = 1;
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence_HP[i, k] = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: amplitude[k]);
                    }
                }
                // Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence_HP.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                        spaceConvergence_HP[i, k].InitialValues.Clear();
                        spaceConvergence_HP[i, k].InitialValues_Evaluators.Clear();
                        spaceConvergence_HP[i, k].savetodb = false;
                        spaceConvergence_HP[i, k].TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                        // Exakte Loesung
                        double MaxAmp = spaceConvergence_HP[i, k].maxAmpli;
                        Func<double[], double> uxi = (X) => -MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * a * (spaceConvergence_HP[i, k].rMax * spaceConvergence_HP[i, k].rMax - X[0] * X[0])) / (4 * nu);
                        Func<double[], double> ueta = (X) => MaxAmp * (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * (a * b * (spaceConvergence_HP[i, k].rMax * spaceConvergence_HP[i, k].rMax - X[0] * X[0])) / (X[0] * 4 * nu);
                        Func<double[], double> ur = (X) => 0;
                        Func<double[], double> pressure = (X) => 0;

                        using (var solver = new HelicalMain()) {
                            solver.Init(spaceConvergence_HP[i, k]);
                            solver.RunSolverMode();
                            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                            Assert.IsTrue(spaceConvergence_HP[i, k].R0fixOn, "R0fix must be turned on");
                            Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, "Pressure Reference Point has to be true");
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


                            Console.WriteLine($"ur       L2 Error/{amplitude[k]}: {ur_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                            Console.WriteLine($"uxi      L2 Error/{amplitude[k]}: {uxi_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"ueta     L2 Error/{amplitude[k]}: {ueta_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"pressure L2 Error/{amplitude[k]}: {pressure_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");


                            Assert.LessOrEqual(ur_L2_Error[i] / amplitude[k], 1.0e-11, $"ur L2 Error out of range: {ur_L2_Error[i]:0.###e-00} (should be close to 0.0 )");
                            Assert.LessOrEqual(pressure_L2_Error[i] / amplitude[k], 1.0e-10, $"pressure L2 Error out of range: {pressure_L2_Error[i]:0.###e-00} (should be close to 0.0)");
                        }
                    }
                    using (var gp = new Gnuplot()) {
                        gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                        gp.OutputFile = $"SpatialConvergence Amplitude {amplitude[k]} with R0fix p{pOrder}.png";
                        gp.SetTitle($"Spatial Convergence Hagen Poiseulle Amplitude {amplitude[k]} with R0fix p{pOrder}");
                        gp.PlotLogXLogY(h, ur_L2_Error , title: $"ur L2error", format: new PlotFormat("-sr"));
                        gp.PlotLogXLogY(h, ueta_L2_Error, title: $"ueta L2error", format: new PlotFormat("-^b"));
                        gp.PlotLogXLogY(h, uxi_L2_Error, title: $"uxi L2error", format: new PlotFormat("-xm"));
                        gp.PlotLogXLogY(h, pressure_L2_Error, title: $"psi L2error", format: new PlotFormat("-*k"));
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
                    if (pOrder != 5) {
                        Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                        Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                    }
                    if (pOrder == 5) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // error thresholds are specified for the SECOND!!!!!! finest grid with pOrder == 5 
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        Assert.LessOrEqual(ur_L2_Error[1] / amplitude[k], 1.0e-10, $"ur L2 Error / amplitude out of range: {ur_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                        Assert.LessOrEqual(uxi_L2_Error[1] / amplitude[k], 1.0e-10, $"uxi L2 Error / amplitude out of range: {uxi_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(ueta_L2_Error[1] / amplitude[k], 1.0e-10, $"ueta L2 Error / amplitude out of range: {ueta_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(pressure_L2_Error[1] / amplitude[k], 1.0e-10, $"pressure L2 Error / amplitude out of range: {pressure_L2_Error[1] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                    }
                }
            }
        }


        /// <summary>
        /// ConditionNumberScaling Tests a steady-state Hagen Poiseulle flow Re=10^5 (aka. Pipe flow).
        /// The Reynolds Number is set to 10^5. The amplitude is 40^*10^5 since Re = Amplitude/4
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
        static public void ConditionNumberScaling_HP_Re_100000_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {


            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: 400000);
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
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(gridSeries[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(gridSeries[i].steady, "Should be steady");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                    }

                }
            }
        }


        /// <summary>
        /// ConditionNumberScaling Tests a steady-state Hagen Poiseulle flow Re=10 (aka. Pipe flow).
        /// The Reynolds Number is set to 10. The amplitude is 40 since Re = Amplitude/4
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
        static public void ConditionNumberScaling_HP_Re_10_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {


            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: 40);
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
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(gridSeries[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(gridSeries[i].steady, "Should be steady");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                    }

                }
            }
        }

        /// <summary>
        /// Tests a steady-state Direct vs Iterativ Linear Solver for Hagen Poiseulle flow Re 100.000 (aka. Hagenpoiseulle flow).
        /// The Reynolds Number is set to 100.000 The amplitude is 400.000 since Re = Amplitude/4
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
        // [Test]  Iterative Solver cannot minimize Residual up to 10^(-13). Stops at ... I dont know (04.11.2024 Schahin Akbari)
        static public void Direct_vs_Iterativ_HP_Re_100000_with_R0fix([Values(5)] int pOrder) {
            // Polynomorder 5 take too long for our test runners!
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_with_R0fix(4)

            //###########################################################
            // Direct Solver
            //###########################################################

            int gridSize = 64;
            double h;
            double r_min = 0;
            double maxAmplitude = 400000;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_direct = new HelicalControl();
            HelicalControl spaceConvergence_iterative = new HelicalControl();

            spaceConvergence_direct = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);

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

            double thresholdPsi = 5e-9; // Poly Order 5
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx / maxAmplitude);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi);
            Assert.LessOrEqual(psiErrorLx / maxAmplitude, thresholdPsi, "Error. psiErrorLx not fulfilled");

            double thresholdUr = 1e-11; // Poly Order 5
            Console.WriteLine("The urErrorLx/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx / maxAmplitude);
            Console.WriteLine("If the urErrorLx / maxAmplitude error is {0} < {1} than good :) ", urErrorLx / maxAmplitude, thresholdUr);
            Assert.LessOrEqual(urErrorLx / maxAmplitude, thresholdUr, "Error. urErrorL2 not fulfilled");

            double thresholdUeta = 1e-11; // Poly Order 5
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx / maxAmplitude);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx / maxAmplitude, thresholdUeta);
            Assert.LessOrEqual(uetaErrorLx / maxAmplitude, thresholdUeta, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi = 1e-11; // Poly Order 5
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx / maxAmplitude);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx / maxAmplitude, thresholdUxi);
            Assert.LessOrEqual(uxiErrorLx / maxAmplitude, thresholdUxi, "Error. uxiErrorL2 not fulfilled");
            //###########################################################
            // Iterativ Solver
            //###########################################################
            spaceConvergence_iterative = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);
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

            double thresholdPsi_ = 5e-9; // Poly Order 5
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx / maxAmplitude);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx / maxAmplitude, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            double thresholdUr_ = 1e-11; // Poly Order 5
            Console.WriteLine("The urErrorLx/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx / maxAmplitude);
            Console.WriteLine("If the urErrorLx / maxAmplitude error is {0} < {1} than good :) ", urErrorLx / maxAmplitude, thresholdUr_);
            Assert.LessOrEqual(urErrorLx / maxAmplitude, thresholdUr_, "Error. urErrorL2 not fulfilled");

            double thresholdUeta_ = 1e-11; // Poly Order 5
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx / maxAmplitude);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx / maxAmplitude, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx / maxAmplitude, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi_ = 1e-11; // Poly Order 5
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx / maxAmplitude);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx / maxAmplitude, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx / maxAmplitude, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");


            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterativ.ur, helical_iterativ.ueta, helical_iterativ.uxi, helical_iterativ.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for (int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            double diff_threshold = 4e-5; // Poly Order 4
            Console.WriteLine("L2 norm of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm());
            Console.WriteLine("If L2 norm of the difference is {0} < {1} than good :) ", diff.L2Norm(), diff_threshold);
            Assert.LessOrEqual(diff.L2Norm(), diff_threshold, "Error. L2 norm of the difference not fulfilled");
        }

        /// <summary>
        /// Tests a steady-state Direct vs Iterativ Linear Solver for Hagen Poiseulle flow Re 10 (aka. Hagenpoiseulle flow).
        /// The Reynolds Number is set to 10 The amplitude is 40 since Re = Amplitude/4
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
        // [Test]  Iterative Solver cannot minimize Residual up to 10^(-13). Stops at ... I dont know  (04.11.2024 Schahin Akbari)
        static public void Direct_vs_Iterativ_HP_Re_10_with_R0fix([Values(5)] int pOrder) {
            // Polynomorder 5 take too long for our test runners!
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_with_R0fix(4)

            //###########################################################
            // Direct Solver
            //###########################################################

            int gridSize = 64;
            double h;
            double r_min = 0;
            double maxAmplitude = 40;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_direct = new HelicalControl();
            HelicalControl spaceConvergence_iterative = new HelicalControl();

            spaceConvergence_direct = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);

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

            double thresholdPsi = 5e-9; // Poly Order 5
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx / maxAmplitude);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi);
            Assert.LessOrEqual(psiErrorLx / maxAmplitude, thresholdPsi, "Error. psiErrorLx not fulfilled");

            double thresholdUr = 1e-11; // Poly Order 5
            Console.WriteLine("The urErrorLx/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx / maxAmplitude);
            Console.WriteLine("If the urErrorLx / maxAmplitude error is {0} < {1} than good :) ", urErrorLx / maxAmplitude, thresholdUr);
            Assert.LessOrEqual(urErrorLx / maxAmplitude, thresholdUr, "Error. urErrorL2 not fulfilled");

            double thresholdUeta = 1e-11; // Poly Order 5
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx / maxAmplitude);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx / maxAmplitude, thresholdUeta);
            Assert.LessOrEqual(uetaErrorLx / maxAmplitude, thresholdUeta, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi = 1e-11; // Poly Order 5
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx / maxAmplitude);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx / maxAmplitude, thresholdUxi);
            Assert.LessOrEqual(uxiErrorLx / maxAmplitude, thresholdUxi, "Error. uxiErrorL2 not fulfilled");
            //###########################################################
            // Iterativ Solver
            //###########################################################
            spaceConvergence_iterative = StokesHelical_Ak.Hagen_Poiseulle.HagenPoiseulle(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);
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

            double thresholdPsi_ = 5e-9; // Poly Order 5
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx / maxAmplitude);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx / maxAmplitude, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            double thresholdUr_ = 1e-11; // Poly Order 5
            Console.WriteLine("The urErrorLx/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx / maxAmplitude);
            Console.WriteLine("If the urErrorLx / maxAmplitude error is {0} < {1} than good :) ", urErrorLx / maxAmplitude, thresholdUr_);
            Assert.LessOrEqual(urErrorLx / maxAmplitude, thresholdUr_, "Error. urErrorL2 not fulfilled");

            double thresholdUeta_ = 1e-11; // Poly Order 5
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx / maxAmplitude);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx / maxAmplitude, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx / maxAmplitude, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi_ = 1e-11; // Poly Order 5
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx / maxAmplitude);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx / maxAmplitude, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx / maxAmplitude, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");


            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterativ.ur, helical_iterativ.ueta, helical_iterativ.uxi, helical_iterativ.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for (int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            double diff_threshold = 4e-7; // Poly Order 5
            //double diff_threshold = 4e-5; // Poly Order 4
            Console.WriteLine("L2 norm of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm());
            Console.WriteLine("If L2 norm of the difference is {0} < {1} than good :) ", diff.L2Norm(), diff_threshold);
            Assert.LessOrEqual(diff.L2Norm(), diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }

        #endregion
        #endregion

        /// Centrifugal Flow.
        /// See PhD-Draft from Schahin Akbari
        #region Centrifuge
        #region With_R0fix
        /// <summary>
        /// ConditionNumberScaling Tests a steady-state Centrifuge flow Re 100.000 (aka. Centrifuge flow).
        /// The Reynolds Number is set to 100000
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
        static public void ConditionNumberScaling_CF_Re_100000_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {


            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];
            double MaxAmplitude = 100000;
            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: MaxAmplitude);
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
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(gridSeries[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(gridSeries[i].steady, "Should be steady");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                    }

                }
            }
        }
        /// <summary>
        /// ConditionNumberScaling Tests a steady-state Centrifuge flow Re 10 (aka. Centrifuge flow).
        /// The Reynolds Number is set to 10
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
        static public void ConditionNumberScaling_CF_Re_10_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            double MaxAmplitude = 10;
            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: MaxAmplitude);
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
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(gridSeries[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(gridSeries[i].steady, "Should be steady");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                    }

                }
            }
        }
        /// <summary>
        /// Tests a steady-state Convergence Test for Centrifuge flow Re 10^5 (aka. Centrifuge flow).
        /// The Reynolds Number is set to 100000
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
        static public void Steady_SpatialConv_CF_Re_100000_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            int[] amplitude = new int[] { 100000 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] ur_L2_Error = new double[gridSize.Length];
            double[] ueta_L2_Error = new double[gridSize.Length];
            double[] uxi_L2_Error = new double[gridSize.Length];
            double[] pressure_L2_Error = new double[gridSize.Length];

            HelicalControl[,] spaceConvergence_HP = new HelicalControl[gridSize.Length, amplitude.Length];
            SinglePhaseField[] exSol_ur = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_ueta = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_uxi = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_pressure = new SinglePhaseField[gridSize.Length];

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            for (int ell = 0; ell < r_min.Length; ell++) {
                //ilPSP.Environment.NumThreads = 1;
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence_HP[i, k] = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: amplitude[k]);
                    }
                }
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                        spaceConvergence_HP[i, k].InitialValues.Clear();
                        spaceConvergence_HP[i, k].InitialValues_Evaluators.Clear();
                        spaceConvergence_HP[i, k].savetodb = false;
                        spaceConvergence_HP[i, k].TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                        // Exact Solution for steady state!
                        double MaxAmp = spaceConvergence_HP[i, k].maxAmpli;
                        Func<double[], double> ur = (X) => 0;
                        Func<double[], double> pressure = (X) => 0;
                        Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                        Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;

                        using (var solver = new HelicalMain()) {
                            solver.Init(spaceConvergence_HP[i, k]);
                            solver.RunSolverMode();
                            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                            Assert.IsTrue(spaceConvergence_HP[i, k].R0fixOn, "R0fix must be turned on");
                            Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, "Pressure Reference Point has to be true");
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


                            Console.WriteLine($"ur       L2 Error/{amplitude[k]}: {ur_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                            Console.WriteLine($"uxi      L2 Error/{amplitude[k]}: {uxi_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"ueta     L2 Error/{amplitude[k]}: {ueta_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"pressure L2 Error/{amplitude[k]}: {pressure_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");


                            Assert.LessOrEqual(ur_L2_Error[i] / amplitude[k], 1.0e-11, $"ur L2 Error / amplitude out of range: {ur_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                            Assert.LessOrEqual(pressure_L2_Error[i] / amplitude[k], 1.0e-10, $"Pressure L2 Error  / amplitude out of range: {pressure_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                        }
                    }

                    using (var gp = new Gnuplot()) {
                        gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                        gp.OutputFile = $"SpatialConvergence Amplitude {amplitude[k]} with R0fix p{pOrder}.png";
                        gp.SetTitle($"Spatial Convergence Centrifugal Flow Amplitude {amplitude[k]} with R0fix p{pOrder}");
                        gp.PlotLogXLogY(h, ur_L2_Error, title: $"ur L2error", format: new PlotFormat("-sr"));
                        gp.PlotLogXLogY(h, ueta_L2_Error, title: $"ueta L2error", format: new PlotFormat("-^b"));
                        gp.PlotLogXLogY(h, uxi_L2_Error, title: $"uxi L2error", format: new PlotFormat("-xm"));
                        gp.PlotLogXLogY(h, pressure_L2_Error, title: $"psi L2error", format: new PlotFormat("-*k"));
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
                    if (pOrder != 5) {
                        Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                        Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                    }
                    if (pOrder == 5) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // error thresholds are specified for the finest grid with pOrder == 5
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                        Assert.LessOrEqual(ur_L2_Error[0] / amplitude[k], 1.0e-10, $"ur L2 Error out of range: {ur_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                        Assert.LessOrEqual(uxi_L2_Error[0] / amplitude[k], 1.0e-10, $"uxi L2 Error out of range: {uxi_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(ueta_L2_Error[0] / amplitude[k], 1.0e-10, $"ueta L2 Error out of range: {ueta_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(pressure_L2_Error[0] / amplitude[k], 1.0e-10, $"pressure L2 Error out of range: {pressure_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                    }
                }
            }
        }
        /// <summary>
        /// Tests a steady-state Convergence Test for Centrifuge flow Re 10 (aka. Centrifuge flow).
        /// The Reynolds Number is set to 10
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
        static public void Steady_SpatialConv_CF_Re_10_Stokes_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            int[] amplitude = new int[] { 10 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] ur_L2_Error = new double[gridSize.Length];
            double[] ueta_L2_Error = new double[gridSize.Length];
            double[] uxi_L2_Error = new double[gridSize.Length];
            double[] pressure_L2_Error = new double[gridSize.Length];

            HelicalControl[,] spaceConvergence_HP = new HelicalControl[gridSize.Length, amplitude.Length];
            SinglePhaseField[] exSol_ur = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_ueta = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_uxi = new SinglePhaseField[gridSize.Length];
            SinglePhaseField[] exSol_pressure = new SinglePhaseField[gridSize.Length];

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            for (int ell = 0; ell < r_min.Length; ell++) {
                //ilPSP.Environment.NumThreads = 1;
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence_HP[i, k] = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], dtRefining: 1, Tend: 1.0e50, MaxAmp: amplitude[k]);
                    }
                }
                for (int k = 0; k < amplitude.Length; k++) {
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                        spaceConvergence_HP[i, k].InitialValues.Clear();
                        spaceConvergence_HP[i, k].InitialValues_Evaluators.Clear();
                        spaceConvergence_HP[i, k].savetodb = false;
                        spaceConvergence_HP[i, k].TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

                        // Exact Solution for steady state!
                        double MaxAmp = spaceConvergence_HP[i, k].maxAmpli;
                        Func<double[], double> ur = (X) => 0;
                        Func<double[], double> pressure = (X) => 0;
                        Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * MaxAmp * X[0];
                        Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * MaxAmp;

                        using (var solver = new HelicalMain()) {
                            solver.Init(spaceConvergence_HP[i, k]);
                            solver.RunSolverMode();
                            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
                            Assert.IsTrue(spaceConvergence_HP[i, k].R0fixOn, "R0fix must be turned on");
                            Assert.IsTrue(spaceConvergence_HP[i, k].PressureReferencePoint, "Pressure Reference Point has to be true");
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

                            Console.WriteLine($"ur       L2 Error/{amplitude[k]}: {ur_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                            Console.WriteLine($"uxi      L2 Error/{amplitude[k]}: {uxi_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"ueta     L2 Error/{amplitude[k]}: {ueta_L2_Error[i] / amplitude[k]:0.###e-00}");
                            Console.WriteLine($"pressure L2 Error/{amplitude[k]}: {pressure_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0)");


                            Assert.LessOrEqual(ur_L2_Error[i] / amplitude[k], 1.0e-11, $"ur L2 Error / amplitude out of range: {ur_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                            Assert.LessOrEqual(pressure_L2_Error[i] / amplitude[k], 1.0e-10, $"Pressure L2 Error/ amplitude out of range: {pressure_L2_Error[i] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                        }
                    }

                    using (var gp = new Gnuplot()) {
                        gp.Terminal = string.Format("pngcairo size {0},{1}", 1024, 768);
                        gp.OutputFile = $"SpatialConvergence Amplitude {amplitude[k]} with R0fix p{pOrder}.png";
                        gp.SetTitle($"Spatial Convergence Centrifugal Flow Amplitude {amplitude[k]} with R0fix p{pOrder}");
                        gp.PlotLogXLogY(h, ur_L2_Error, title: $"ur L2error", format: new PlotFormat("-sr"));
                        gp.PlotLogXLogY(h, ueta_L2_Error, title: $"ueta L2error", format: new PlotFormat("-^b"));
                        gp.PlotLogXLogY(h, uxi_L2_Error, title: $"uxi L2error", format: new PlotFormat("-xm"));
                        gp.PlotLogXLogY(h, pressure_L2_Error, title: $"psi L2error", format: new PlotFormat("-*k"));
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
                    if (pOrder != 5) {
                        Assert.GreaterOrEqual(regressionUeta.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Ueta lower than expected (was {0} but should be {1} +- {2})", regressionUeta.slope, pOrder + 1, tolerance));
                        Assert.GreaterOrEqual(regressionUxi.slope, pOrder + 1 - tolerance, String.Format("Convergence rate for Uxi lower than expected (was {0} but should be {1} +- {2})", regressionUxi.slope, pOrder + 1, tolerance));
                    }
                    if (pOrder == 5) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // error thresholds are specified for the finest grid with pOrder == 5
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                        Assert.LessOrEqual(ur_L2_Error[0] / amplitude[k], 1.0e-10, $"ur L2 Error / amplitude out of range: {ur_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0 )");
                        Assert.LessOrEqual(uxi_L2_Error[0] / amplitude[k], 1.0e-10, $"uxi L2 Error / amplitude out of range: {uxi_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(ueta_L2_Error[0] / amplitude[k], 1.0e-10, $"ueta L2 Error / amplitude out of range: {ueta_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                        Assert.LessOrEqual(pressure_L2_Error[0] / amplitude[k], 1.0e-9, $"pressure L2 Error / amplitude out of range: {pressure_L2_Error[0] / amplitude[k]:0.###e-00} (should be close to 0.0)");
                    }
                }
            }
        }
        /// <summary>
        /// Tests a steady-state Direct vs Iterativ Linear Solver for Centrifuge flow Re 100000 (aka. Centrifuge flow).
        /// The Reynolds Number is set to 100000
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
        // [Test]  Iterative Solver cannot minimize Residual up to 10^(-13). Stops at ... I dont know  (04.11.2024 Schahin Akbari)
        static public void Direct_vs_Iterativ_CF_Re_100000_with_R0fix([Values(5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_with_R0fix(4)

            //###########################################################
            // Direct Solver
            //###########################################################

            int gridSize = 64;
            double h;
            double r_min = 0;
            double maxAmplitude = 100000;
            double ur_L2_Error;
            double ueta_L2_Error;
            double uxi_L2_Error;
            double pressure_L2_Error;
            HelicalControl spaceConvergence_direct = new HelicalControl();
            HelicalControl spaceConvergence_iterative = new HelicalControl();

            SinglePhaseField exSol_ur;
            SinglePhaseField exSol_ueta;
            SinglePhaseField exSol_uxi;
            SinglePhaseField exSol_pressure;

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            Func<double[], double> ur = (X) => 0;
            Func<double[], double> pressure = (X) => 0;
            Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * maxAmplitude * X[0];
            Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * maxAmplitude;

            spaceConvergence_direct = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);


            var helical_direct = new HelicalMain();
            helical_direct.Init(spaceConvergence_direct);
            helical_direct.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_direct.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_direct.PressureReferencePoint, "Pressure Reference Point has to be true");


            // Deep Copy
            exSol_ur = helical_direct.ur.CloneAs();
            exSol_ueta = helical_direct.ueta.CloneAs();
            exSol_uxi = helical_direct.uxi.CloneAs();
            exSol_pressure = helical_direct.Pressure.CloneAs();

            //Exakte Lösung auf SinglePhaseField
            exSol_ur.ProjectField(ur);
            exSol_ueta.ProjectField(ueta);
            exSol_uxi.ProjectField(uxi);
            exSol_pressure.ProjectField(pressure);

            //L2 Fehler berechnen. Meine Loesung vs richtige Loesung
            ur_L2_Error = helical_direct.ur.L2Error(exSol_ur);
            ueta_L2_Error = helical_direct.ueta.L2Error(exSol_ueta);
            uxi_L2_Error = helical_direct.uxi.L2Error(exSol_uxi);
            pressure_L2_Error = helical_direct.Pressure.L2Error(exSol_pressure);

            h = 2 * Math.PI / gridSize;

            double thresholdUr = 1e-11; // Poly Order 5
            Console.WriteLine("The ur_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ur_L2_Error / maxAmplitude);
            Console.WriteLine("If the ur_L2_Error / maxAmplitude error is {0} < {1} than good :) ", ur_L2_Error / maxAmplitude, thresholdUr);
            Assert.LessOrEqual(ur_L2_Error / maxAmplitude, thresholdUr, "Error. ur_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUeta = 1e-11; // Poly Order 5
            Console.WriteLine("The ueta_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ueta_L2_Error / maxAmplitude);
            Console.WriteLine("If the ueta_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", ueta_L2_Error / maxAmplitude, thresholdUeta);
            Assert.LessOrEqual(ueta_L2_Error / maxAmplitude, thresholdUeta, "Error. ueta_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUxi = 1e-11; // Poly Order 5
            Console.WriteLine("The uxi_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, uxi_L2_Error / maxAmplitude);
            Console.WriteLine("If the uxi_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", uxi_L2_Error / maxAmplitude, thresholdUxi);
            Assert.LessOrEqual(uxi_L2_Error / maxAmplitude, thresholdUxi, "Error. uxi_L2_Error/ maxAmplitude not fulfilled");

            double thresholdPsi = 5e-9; // Poly Order 5
            Console.WriteLine("The pressure_L2_Error / maxAmplitude error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, pressure_L2_Error / maxAmplitude);
            Console.WriteLine("If the pressure_L2_Error / maxAmplitude error is {0} < {1} than good :) ", pressure_L2_Error / maxAmplitude, thresholdPsi);
            Assert.LessOrEqual(pressure_L2_Error / maxAmplitude, thresholdPsi, "Error. pressure_L2_Error not fulfilled");
            //###########################################################
            // Iterativ Solver
            //###########################################################


            spaceConvergence_iterative = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);
            spaceConvergence_iterative.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { };


            var helical_iterative = new HelicalMain();
            helical_iterative.Init(spaceConvergence_iterative);
            helical_iterative.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_iterative.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_iterative.PressureReferencePoint, "Pressure Reference Point has to be true");


            // Deep Copy
            exSol_ur = helical_iterative.ur.CloneAs();
            exSol_ueta = helical_iterative.ueta.CloneAs();
            exSol_uxi = helical_iterative.uxi.CloneAs();
            exSol_pressure = helical_iterative.Pressure.CloneAs();

            //Exakte Lösung auf SinglePhaseField
            exSol_ur.ProjectField(ur);
            exSol_ueta.ProjectField(ueta);
            exSol_uxi.ProjectField(uxi);
            exSol_pressure.ProjectField(pressure);

            //L2 Fehler berechnen. Meine Loesung vs richtige Loesung
            ur_L2_Error = helical_iterative.ur.L2Error(exSol_ur);
            ueta_L2_Error = helical_iterative.ueta.L2Error(exSol_ueta);
            uxi_L2_Error = helical_iterative.uxi.L2Error(exSol_uxi);
            pressure_L2_Error = helical_iterative.Pressure.L2Error(exSol_pressure);

            h = 2 * Math.PI / gridSize;

            double thresholdUr_ = 1e-11; // Poly Order 5
            Console.WriteLine("The ur_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ur_L2_Error / maxAmplitude);
            Console.WriteLine("If the ur_L2_Error / maxAmplitude error is {0} < {1} than good :) ", ur_L2_Error / maxAmplitude, thresholdUr_);
            Assert.LessOrEqual(ur_L2_Error / maxAmplitude, thresholdUr_, "Error. ur_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUeta_ = 1e-11; // Poly Order 5
            Console.WriteLine("The ueta_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ueta_L2_Error / maxAmplitude);
            Console.WriteLine("If the ueta_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", ueta_L2_Error / maxAmplitude, thresholdUeta_);
            Assert.LessOrEqual(ueta_L2_Error / maxAmplitude, thresholdUeta_, "Error. ueta_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUxi_ = 1e-11; // Poly Order 5
            Console.WriteLine("The uxi_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, uxi_L2_Error / maxAmplitude);
            Console.WriteLine("If the uxi_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", uxi_L2_Error / maxAmplitude, thresholdUxi_);
            Assert.LessOrEqual(uxi_L2_Error / maxAmplitude, thresholdUxi_, "Error. uxi_L2_Error/ maxAmplitude not fulfilled");

            double thresholdPsi_ = 5e-9; // Poly Order 5
            Console.WriteLine("The pressure_L2_Error / maxAmplitude error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, pressure_L2_Error / maxAmplitude);
            Console.WriteLine("If the pressure_L2_Error / maxAmplitude error is {0} < {1} than good :) ", pressure_L2_Error / maxAmplitude, thresholdPsi_);
            Assert.LessOrEqual(pressure_L2_Error / maxAmplitude, thresholdPsi_, "Error. pressure_L2_Error not fulfilled");


            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterative.ur, helical_iterative.ueta, helical_iterative.uxi, helical_iterative.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for (int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            double diff_threshold = 4e-7; // Poly Order 5
            //double diff_threshold = 4e-5; // Poly Order 4
            Console.WriteLine("L2 norm / maxAmplitude of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm() / maxAmplitude);
            Console.WriteLine("If L2 norm / maxAmplitude of the difference is {0} < {1} than good :) ", diff.L2Norm() / maxAmplitude, diff_threshold);
            Assert.LessOrEqual(diff.L2Norm() / maxAmplitude, diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }

        /// <summary>
        /// Tests a steady-state Direct vs Iterativ Linear Solver for Centrifuge flow Re 10 (aka. Centrifuge flow).
        /// The Reynolds Number is set to 10
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
        // [Test]  Iterative Solver cannot minimize Residual up to 10^(-13). Stops at 10^(-6) I think (04.11.2024 Schahin Akbari)
        static public void Direct_vs_Iterativ_CF_Re_10_with_R0fix([Values(5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_with_R0fix(5)

            //###########################################################
            // Direct Solver
            //###########################################################

            int gridSize = 64;
            double h;
            double r_min = 0;
            double maxAmplitude = 10;
            double ur_L2_Error;
            double ueta_L2_Error;
            double uxi_L2_Error;
            double pressure_L2_Error;
            HelicalControl spaceConvergence_direct = new HelicalControl();
            HelicalControl spaceConvergence_iterative = new HelicalControl();

            SinglePhaseField exSol_ur;
            SinglePhaseField exSol_ueta;
            SinglePhaseField exSol_uxi;
            SinglePhaseField exSol_pressure;

            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            Func<double[], double> ur = (X) => 0;
            Func<double[], double> pressure = (X) => 0;
            Func<double[], double> ueta = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * a * maxAmplitude * X[0];
            Func<double[], double> uxi = (X) => (X[0] / (Math.Sqrt(a * a * X[0] * X[0] + b * b))) * b * maxAmplitude;

            spaceConvergence_direct = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);


            var helical_direct = new HelicalMain();
            helical_direct.Init(spaceConvergence_direct);
            helical_direct.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_direct.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_direct.PressureReferencePoint, "Pressure Reference Point has to be true");


            // Deep Copy
            exSol_ur = helical_direct.ur.CloneAs();
            exSol_ueta = helical_direct.ueta.CloneAs();
            exSol_uxi = helical_direct.uxi.CloneAs();
            exSol_pressure = helical_direct.Pressure.CloneAs();

            //Exakte Lösung auf SinglePhaseField
            exSol_ur.ProjectField(ur);
            exSol_ueta.ProjectField(ueta);
            exSol_uxi.ProjectField(uxi);
            exSol_pressure.ProjectField(pressure);

            //L2 Fehler berechnen. Meine Loesung vs richtige Loesung
            ur_L2_Error = helical_direct.ur.L2Error(exSol_ur);
            ueta_L2_Error = helical_direct.ueta.L2Error(exSol_ueta);
            uxi_L2_Error = helical_direct.uxi.L2Error(exSol_uxi);
            pressure_L2_Error = helical_direct.Pressure.L2Error(exSol_pressure);

            h = 2 * Math.PI / gridSize;

            double thresholdUr = 1e-11; // Poly Order 5
            Console.WriteLine("The ur_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ur_L2_Error / maxAmplitude);
            Console.WriteLine("If the ur_L2_Error / maxAmplitude error is {0} < {1} than good :) ", ur_L2_Error / maxAmplitude, thresholdUr);
            Assert.LessOrEqual(ur_L2_Error / maxAmplitude, thresholdUr, "Error. ur_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUeta = 1e-11; // Poly Order 5
            Console.WriteLine("The ueta_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ueta_L2_Error / maxAmplitude);
            Console.WriteLine("If the ueta_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", ueta_L2_Error / maxAmplitude, thresholdUeta);
            Assert.LessOrEqual(ueta_L2_Error / maxAmplitude, thresholdUeta, "Error. ueta_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUxi = 1e-11; // Poly Order 5
            Console.WriteLine("The uxi_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, uxi_L2_Error / maxAmplitude);
            Console.WriteLine("If the uxi_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", uxi_L2_Error / maxAmplitude, thresholdUxi);
            Assert.LessOrEqual(uxi_L2_Error / maxAmplitude, thresholdUxi, "Error. uxi_L2_Error/ maxAmplitude not fulfilled");

            double thresholdPsi = 5e-9; // Poly Order 5
            Console.WriteLine("The pressure_L2_Error / maxAmplitude error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, pressure_L2_Error / maxAmplitude);
            Console.WriteLine("If the pressure_L2_Error / maxAmplitude error is {0} < {1} than good :) ", pressure_L2_Error / maxAmplitude, thresholdPsi);
            Assert.LessOrEqual(pressure_L2_Error / maxAmplitude, thresholdPsi, "Error. pressure_L2_Error not fulfilled");
            //###########################################################
            // Iterativ Solver
            //###########################################################


            spaceConvergence_iterative = StokesHelical_Ak.Centrifuge.Centrifuge_Flow(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, dtRefining: 1, Tend: 1.0e50, MaxAmp: maxAmplitude);
            spaceConvergence_iterative.LinearSolver = new BoSSS.Solution.AdvancedSolvers.OrthoMGSchwarzConfig() { };


            var helical_iterative = new HelicalMain();
            helical_iterative.Init(spaceConvergence_iterative);
            helical_iterative.RunSolverMode();
            Assert.AreEqual(Globals.activeMult, Globals.Multiplier.Bsq, $"Multiplier expected to be {Globals.Multiplier.Bsq}");
            Assert.IsTrue(spaceConvergence_iterative.R0fixOn, "R0fix must be turned on");
            Assert.IsTrue(spaceConvergence_iterative.PressureReferencePoint, "Pressure Reference Point has to be true");


            // Deep Copy
            exSol_ur = helical_iterative.ur.CloneAs();
            exSol_ueta = helical_iterative.ueta.CloneAs();
            exSol_uxi = helical_iterative.uxi.CloneAs();
            exSol_pressure = helical_iterative.Pressure.CloneAs();

            //Exakte Lösung auf SinglePhaseField
            exSol_ur.ProjectField(ur);
            exSol_ueta.ProjectField(ueta);
            exSol_uxi.ProjectField(uxi);
            exSol_pressure.ProjectField(pressure);

            //L2 Fehler berechnen. Meine Loesung vs richtige Loesung
            ur_L2_Error = helical_iterative.ur.L2Error(exSol_ur);
            ueta_L2_Error = helical_iterative.ueta.L2Error(exSol_ueta);
            uxi_L2_Error = helical_iterative.uxi.L2Error(exSol_uxi);
            pressure_L2_Error = helical_iterative.Pressure.L2Error(exSol_pressure);

            h = 2 * Math.PI / gridSize;

            double thresholdUr_ = 1e-11; // Poly Order 5
            Console.WriteLine("The ur_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ur_L2_Error / maxAmplitude);
            Console.WriteLine("If the ur_L2_Error / maxAmplitude error is {0} < {1} than good :) ", ur_L2_Error / maxAmplitude, thresholdUr_);
            Assert.LessOrEqual(ur_L2_Error / maxAmplitude, thresholdUr_, "Error. ur_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUeta_ = 1e-11; // Poly Order 5
            Console.WriteLine("The ueta_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, ueta_L2_Error / maxAmplitude);
            Console.WriteLine("If the ueta_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", ueta_L2_Error / maxAmplitude, thresholdUeta_);
            Assert.LessOrEqual(ueta_L2_Error / maxAmplitude, thresholdUeta_, "Error. ueta_L2_Error/ maxAmplitude not fulfilled");

            double thresholdUxi_ = 1e-11; // Poly Order 5
            Console.WriteLine("The uxi_L2_Error/ maxAmplitude error for {0} xi-Cells and {1} is = {2}", h, r_min, uxi_L2_Error / maxAmplitude);
            Console.WriteLine("If the uxi_L2_Error/ maxAmplitude error is {0} < {1} than good :) ", uxi_L2_Error / maxAmplitude, thresholdUxi_);
            Assert.LessOrEqual(uxi_L2_Error / maxAmplitude, thresholdUxi_, "Error. uxi_L2_Error/ maxAmplitude not fulfilled");

            double thresholdPsi_ = 5e-9; // Poly Order 5
            Console.WriteLine("The pressure_L2_Error / maxAmplitude error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, pressure_L2_Error / maxAmplitude);
            Console.WriteLine("If the pressure_L2_Error / maxAmplitude error is {0} < {1} than good :) ", pressure_L2_Error / maxAmplitude, thresholdPsi_);
            Assert.LessOrEqual(pressure_L2_Error / maxAmplitude, thresholdPsi_, "Error. pressure_L2_Error not fulfilled");


            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterative.ur, helical_iterative.ueta, helical_iterative.uxi, helical_iterative.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for (int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            double diff_threshold = 4e-7; // Poly Order 5
            //double diff_threshold = 4e-5; // Poly Order 4
            Console.WriteLine("L2 norm / maxAmplitude of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm() / maxAmplitude);
            Console.WriteLine("If L2 norm / maxAmplitude of the difference is {0} < {1} than good :) ", diff.L2Norm() / maxAmplitude, diff_threshold);
            Assert.LessOrEqual(diff.L2Norm() / maxAmplitude, diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }

        #endregion

        #endregion
        // Here is the Stuff of Doktor Dokinik Dierkes In his Paper Formula 4.4
        #region DoktorDominikDierkes

        // With R0_fix
        #region With_R0fix
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void Steady_SpatialConvergence_DDD_Paper_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialConvergence_without_R0fix(X)

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0 };
            double[] urErrorLx = new double[gridSize.Length];
            double[] uetaErrorLx = new double[gridSize.Length];
            double[] uxiErrorLx = new double[gridSize.Length];
            double[] psiErrorLx = new double[gridSize.Length];
            HelicalControl[] spaceConvergence = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence[i] = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell], steady: true);
                    }

                    Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(spaceConvergence[i].PressureReferencePoint, $"Pressure Reference Point must be true (i = {i})");
                    }

                    for (int i = 0; i < gridSize.Length; i++) {
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

                    using (var gp = new Gnuplot()) {
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




            if (pOrder == 5) {
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
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void Direct_vs_Iterativ_DDD_Paper_with_R0fix([Values(5)] int pOrder) {
            // Polynomorder 5 take too long for our test runners!
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialComparison_Direct_vs_Iterativ_with_R0fix(4)

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

            spaceConvergence_direct = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min, steady: true);

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

            double thresholdPsi = 5e-9; // Poly Order 5
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi, "Error. psiErrorLx not fulfilled");

            double thresholdUr = 1e-11; // Poly Order 5
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr);
            Assert.LessOrEqual(urErrorLx, thresholdUr, "Error. urErrorL2 not fulfilled");

            double thresholdUeta = 1e-11; // Poly Order 5
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta, "Error. uetaErrorL2 not fulfilled");

            double thresholdUxi = 1e-11; // Poly Order 5
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi, "Error. uxiErrorL2 not fulfilled");

            //###########################################################
            // Iterativ Solver
            //###########################################################

            spaceConvergence_iterative = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min, steady: true);
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

            //double thresholdPsi_ = 5e-9; // Poly Order 5
            double thresholdPsi_ = 5e-7; // Poly Order 4
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            //double thresholdUr_ = 1e-10; // Poly Order 5
            double thresholdUr_ = 1e-8; // Poly Order 4
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr_);
            Assert.LessOrEqual(urErrorLx, thresholdUr_, "Error. urErrorL2 not fulfilled");

            //double thresholdUeta_ = 2e-8; // Poly Order 5
            double thresholdUeta_ = 2e-6; // Poly Order 4
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            //double thresholdUxi_ = 1e-8; // Poly Order 5
            double thresholdUxi_ = 1e-6; // Poly Order 4
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");

            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterativ.ur, helical_iterativ.ueta, helical_iterativ.uxi, helical_iterativ.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for (int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            //double diff_threshold = 4e-7; // Poly Order 5
            double diff_threshold = 4e-5; // Poly Order 4
            Console.WriteLine("L2 norm of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm());
            Console.WriteLine("If L2 norm of the difference is {0} < {1} than good :) ", diff.L2Norm(), diff_threshold);
            Assert.LessOrEqual(diff.L2Norm(), diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }
        /// <summary>
        /// FISCH
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void ConditionNumberScaling_DDD_Paper_with_R0fix([Values(2, 3, 4, 5)] int pOrder) {


            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell], steady: true);
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
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsTrue(gridSeries[i].R0fixOn, "R0fix must be turned on");
                        Assert.IsTrue(gridSeries[i].steady, "Should be steady");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "Pressure Reference Point has to be true");
                    }

                }
            }
        }
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void HangingNodes_DDD_Paper_with_R0fix([Values(4)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.HangingNodes_with_R0fix(5)

            int gridSize = 64;

            double r_min = 0;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_hangingNodes = new HelicalControl();
            spaceConvergence_hangingNodes = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min, steady: true);

            // Grid
            #region Grid
            spaceConvergence_hangingNodes.GridFunc = delegate {
                int numberOfSegments = 4;
                // Initialisierung einer Liste für die Gitter
                GridCommons[] grids = new GridCommons[numberOfSegments];
                // Basisanzahl der Zellen im ersten Segment
                int N = (gridSize + 1) / numberOfSegments;
                // Startpunkt für xNodes
                double xStart = r_min;
                // Die Schrittweite für xNodes-Bereiche erhöhen
                double xStep = (spaceConvergence_hangingNodes.rMax / numberOfSegments) - xStart;

                for (int i = 0; i < numberOfSegments; i++) {
                    // Berechnung der Knotenanzahl für das aktuelle Segment
                    int nodesInSegment = N * (int)Math.Pow(2, i) + 1;

                    // Bestimmung der x-Bereichsgrenzen für das aktuelle Segment
                    double xEnd = xStart + xStep;

                    // Erzeugung der x-Knoten und y-Knoten für das aktuelle Segment
                    double[] xNodes = GenericBlas.Linspace(xStart, xEnd, nodesInSegment);
                    double[] yNodes = GenericBlas.Linspace(0, 2 * Math.PI, nodesInSegment);

                    // Erstellung des Gitters für das aktuelle Segment und Hinzufügen in die Liste
                    grids[i] = Grid2D.Cartesian2DGrid(xNodes, yNodes, type: CellType.Square_Linear);
                    // Aktualisierung des Startpunkts für den nächsten Durchlauf
                    xStart = xEnd;
                }

                // Zusammenführen der Gitter aus der Liste
                var grdJ = grids[0];
                for (int i = 1; i < grids.Length; i++) {
                    grdJ = GridCommons.MergeLogically(grids[i], grdJ);
                }
                // Versiegelung des zusammengeführten Gitters
                GridCommons grd = GridCommons.Seal(grdJ, 4);

                // Hinzufügen der Edge Tag Names und Definition der Edge Tags...
                // (Füge hier den restlichen Teil deines Codes ein)

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if (Math.Abs(r - spaceConvergence_hangingNodes.rMax) < 1E-8 || (spaceConvergence_hangingNodes.rMin >= 1E-6 && Math.Abs(r - spaceConvergence_hangingNodes.rMin) < 1E-8))
                        return 1;
                    else
                        return 2;
                });
                return grd;
            };
            #endregion


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
        #endregion
        // Without R0_fix
        #region Without R0fix
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void SpatialConvergence_without_R0fix_DDD_Paper([Values(2, 3, 4, 5)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.SpatialConvergence_without_R0fix 

            ilPSP.Environment.NumThreads = 1;

            int[] gridSize = new int[] { 64, 32, 16, 8, 4 };
            double[] h = new double[gridSize.Length];
            double[] r_min = new double[] { 0.1 };
            double[] urErrorL2 = new double[gridSize.Length];
            double[] uetaErrorL2 = new double[gridSize.Length];
            double[] uxiErrorL2 = new double[gridSize.Length];
            double[] psiErrorL2 = new double[gridSize.Length];
            HelicalControl[] spaceConvergence = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        spaceConvergence[i] = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell], steady: true);
                    }
                    Console.WriteLine($"pOrder = {pOrder}, ell = {ell}, {spaceConvergence.Select(C => C.PressureReferencePoint).ToConcatString("[", ", ", "]")}");

                    for (int i = 0; i < gridSize.Length; i++) {
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

                    using (var gp = new Gnuplot()) {
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

            if (pOrder == 5) {
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
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void Direct_vs_Iterativ_without_R0fix_DDD_Paper([Values(5)] int pOrder) {
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

            spaceConvergence_direct = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min, steady: true);

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

            //double thresholdPsi = 5e-9;   // Poly Order 5
            double thresholdPsi = 5e-7;     // Poly Order 4
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi, "Error. psiErrorLx not fulfilled");

            //double thresholdUr = 1e-11;// Poly Order 5
            double thresholdUr = 1e-9;// Poly Order 4
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr);
            Assert.LessOrEqual(urErrorLx, thresholdUr, "Error. urErrorL2 not fulfilled");

            //double thresholdUeta = 1e-11;// Poly Order 5
            double thresholdUeta = 1e-9;// Poly Order 4
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta, "Error. uetaErrorL2 not fulfilled");

            //double thresholdUxi = 1e-11;// Poly Order 5
            double thresholdUxi = 1e-9;// Poly Order 4
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi, "Error. uxiErrorL2 not fulfilled");

            //###########################################################
            // Iterativ Solver
            //###########################################################

            spaceConvergence_iterative = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min, steady: true);
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

            //double thresholdPsi_ = 5e-9; // Poly Order 5
            double thresholdPsi_ = 5e-7; // Poly Order 4
            Console.WriteLine("The psiErrorLx error for {0} xi-Cells and rMin = {1} is = {2}", h, r_min, psiErrorLx);
            Console.WriteLine("If the psiErrorLx error is {0} < {1} than good :) ", psiErrorLx, thresholdPsi_);
            Assert.LessOrEqual(psiErrorLx, thresholdPsi_, "Error. psiErrorLx not fulfilled");

            //double thresholdUr_ = 1e-10; // Poly Order 5
            double thresholdUr_ = 1e-8; // Poly Order 4
            Console.WriteLine("The urErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, urErrorLx);
            Console.WriteLine("If the urErrorLx error is {0} < {1} than good :) ", urErrorLx, thresholdUr_);
            Assert.LessOrEqual(urErrorLx, thresholdUr_, "Error. urErrorL2 not fulfilled");

            //double thresholdUeta_ = 2e-8; // Poly Order 5
            double thresholdUeta_ = 2e-6; // Poly Order 4
            Console.WriteLine("The uetaErrorLx error for {0} xi-Cells and {1} is = {2}", h, r_min, uetaErrorLx);
            Console.WriteLine("If the uetaErrorLx error is {0} < {1} than good :) ", uetaErrorLx, thresholdUeta_);
            Assert.LessOrEqual(uetaErrorLx, thresholdUeta_, "Error. uetaErrorL2 not fulfilled");

            //double thresholdUxi_ = 1e-8;  // Poly Order 5
            double thresholdUxi_ = 1e-6;  // Poly Order 4
            Console.WriteLine("The uxiErrorL2 error for {0} xi-Cells and {1} is = {2}", h, r_min, uxiErrorLx);
            Console.WriteLine("If the uxiErrorL2 error is {0} < {1} than good :) ", uxiErrorLx, thresholdUxi_);
            Assert.LessOrEqual(uxiErrorLx, thresholdUxi_, "Error. uxiErrorL2 not fulfilled");

            CoordinateMapping helicalSol_direct = new CoordinateMapping(helical_direct.ur, helical_direct.ueta, helical_direct.uxi, helical_direct.Pressure);
            CoordinateVector helicalSolVec_direct = new CoordinateVector(helicalSol_direct);

            CoordinateMapping helicalSol_iterativ = new CoordinateMapping(helical_iterativ.ur, helical_iterativ.ueta, helical_iterativ.uxi, helical_iterativ.Pressure);
            CoordinateVector helicalSolVec_iterativ = new CoordinateVector(helicalSol_iterativ);

            // calculate the difference
            double[] diff = new double[helicalSolVec_iterativ.Length];
            for (int i = 0; i < helicalSolVec_iterativ.Length; i++) {
                diff[i] = helicalSolVec_iterativ[i] - helicalSolVec_direct[i];
            }

            //double diff_threshold = 4e-7; // Poly Order 5
            double diff_threshold = 4e-5; // Poly Order 4
            Console.WriteLine("L2 norm of the difference between iterativ Solver and direct solver is = {0}", diff.L2Norm());
            Console.WriteLine("If L2 norm of the difference is {0} < {1} than good :) ", diff.L2Norm(), diff_threshold);
            Assert.LessOrEqual(diff.L2Norm(), diff_threshold, "Error. L2 norm of the difference not fulfilled");

        }
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void HangNodes_without_R0fix_DDD_Paper_([Values(4)] int pOrder) {
            // --test=StokesHelical_Ak.TestSpartial.TestSpatial.HangingNodes_with_R0fix(5)

            int gridSize = 64;
            double r_min = 0.1;
            double urErrorLx;
            double uetaErrorLx;
            double uxiErrorLx;
            double psiErrorLx;
            HelicalControl spaceConvergence_hangingNodes = new HelicalControl();
            spaceConvergence_hangingNodes = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize, noOfCellsXi: gridSize, rMin: r_min, steady: true);

            // Grid
            #region Grid
            spaceConvergence_hangingNodes.GridFunc = delegate {
                int numberOfSegments = 4;
                // Initialisierung einer Liste für die Gitter
                GridCommons[] grids = new GridCommons[numberOfSegments];
                // Basisanzahl der Zellen im ersten Segment
                int N = (gridSize + 1) / numberOfSegments;
                // Startpunkt für xNodes
                double xStart = r_min;
                // Die Schrittweite für xNodes-Bereiche erhöhen
                double xStep = (spaceConvergence_hangingNodes.rMax / numberOfSegments) - xStart;

                for (int i = 0; i < numberOfSegments; i++) {
                    // Berechnung der Knotenanzahl für das aktuelle Segment
                    int nodesInSegment = N * (int)Math.Pow(2, i) + 1;

                    // Bestimmung der x-Bereichsgrenzen für das aktuelle Segment
                    double xEnd = xStart + xStep;

                    // Erzeugung der x-Knoten und y-Knoten für das aktuelle Segment
                    double[] xNodes = GenericBlas.Linspace(xStart, xEnd, nodesInSegment);
                    double[] yNodes = GenericBlas.Linspace(0, 2 * Math.PI, nodesInSegment);

                    // Erstellung des Gitters für das aktuelle Segment und Hinzufügen in die Liste
                    grids[i] = Grid2D.Cartesian2DGrid(xNodes, yNodes, type: CellType.Square_Linear);
                    // Aktualisierung des Startpunkts für den nächsten Durchlauf
                    xStart = xEnd;
                }

                // Zusammenführen der Gitter aus der Liste
                var grdJ = grids[0];
                for (int i = 1; i < grids.Length; i++) {
                    grdJ = GridCommons.MergeLogically(grids[i], grdJ);
                }
                // Versiegelung des zusammengeführten Gitters
                GridCommons grd = GridCommons.Seal(grdJ, 4);

                // Hinzufügen der Edge Tag Names und Definition der Edge Tags...
                // (Füge hier den restlichen Teil deines Codes ein)

                grd.EdgeTagNames.Add(1, "Dirichlet");
                grd.EdgeTagNames.Add(2, "Stuff");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double r, xi;
                    r = X[0];
                    xi = X[1];
                    if (Math.Abs(r - spaceConvergence_hangingNodes.rMax) < 1E-8 || (spaceConvergence_hangingNodes.rMin >= 1E-6 && Math.Abs(r - spaceConvergence_hangingNodes.rMin) < 1E-8))
                        return 1;
                    else
                        return 2;
                });
                return grd;
            };
            #endregion
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
        /// <summary>
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// Text
        /// </summary>
        /// <remarks>
        /// 
        /// 
        /// </remarks>
        [Test]
        static public void ConditionNumberScaling_DDD_without_R0fix([Values(2, 3, 4, 5)] int pOrder) {

            //ilPSP.Environment.NumThreads = 1;

            int[] gridSize = new int[] { 4, 8, 16, 32, 64 };
            double[] r_min = new double[] { 0.1 };
            HelicalControl[] gridSeries = new HelicalControl[gridSize.Length];

            for (int ell = 0; ell < r_min.Length; ell++) {
                {
                    for (int i = 0; i < gridSize.Length; i++) {
                        gridSeries[i] = Man_Sol_DDD.ManSol_DDD_Paper(degree: pOrder, noOfCellsR: gridSize[i], noOfCellsXi: gridSize[i], rMin: r_min[ell], steady: true);
                    }

                    ConditionNumberScalingTest.Perform(gridSeries, new ConditionNumberScalingTest.Config() { plot = true, title = "ConditionNumberScaling_without_R0fix-p" + pOrder, ComputeGlobalCondNo = false });
                    for (int i = 0; i < gridSize.Length; i++) {
                        Assert.IsFalse(gridSeries[i].R0fixOn, "R0fix must be turned off");
                        Assert.IsTrue(gridSeries[i].PressureReferencePoint, "PressureReferencePoint should be true. Since R0_fix is false");
                    }
                }
            }
        }
        #endregion

        #endregion



    }
}
