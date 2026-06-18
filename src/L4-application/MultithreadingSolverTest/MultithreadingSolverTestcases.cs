using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using NUnit.Framework;
using System;
using System.Diagnostics;
using System.IO;

namespace BoSSS.Application.MultithreadingSolverTest {

    internal sealed class ThreadScalingTestcase {
        public string Name;

        /// <summary>
        /// Creates a fresh solver control object for one run.
        /// </summary>
        public Func<AppControl> CreateControl;

        /// <summary>
        /// Runs the control object and returns the elapsed wall time in seconds.
        /// </summary>
        public Func<AppControl, double> RunControl;

        /// <summary>
        /// Minimum expected speedup from 1 to 2 threads: T1 / T2.
        /// </summary>
        public double MinSpeedup12 = 1.1;

        /// <summary>
        /// Minimum expected speedup from 1 to 4 threads: T1 / T4.
        /// </summary>
        public double MinSpeedup14 = 1.2;

        /// <summary>
        /// Minimum expected speedup from 2 to 4 threads: T2 / T4.
        /// Use values below 1.0 if 4 threads may be slightly slower than 2.
        /// </summary>
        public double MinSpeedup24 = 1.1;

        /// <summary>
        /// Number of measured runs after the warm-up.
        /// </summary>
        public int Runs = 3;
    }

    internal static class MultithreadingSolverTestcases {

        public static readonly ThreadScalingTestcase[] All = new[] {
            new ThreadScalingTestcase() {
                Name = "ipPossion",
                CreateControl = () => {
                    var C =  BoSSS.Application.SipPoisson.SipHardcodedControl.Square(64,64);
                    C.ProjectName = "MultithreadingScalingTest/ipPoisson";
                    C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
                    PrepareControl(C, "ipPoisson");
                    return C;
                },
                RunControl = control => RunControl<BoSSS.Application.SipPoisson.SipControl, BoSSS.Application.SipPoisson.SipPoissonMain>((BoSSS.Application.SipPoisson.SipControl)control),

            },
            new ThreadScalingTestcase {
                Name = "XdgPoisson3",
                CreateControl = () => {
                    var C = BoSSS.Application.XdgPoisson3.HardCodedControl.Circle(Resolution: 64, p: 2);
                    C.ProjectName = "MultithreadingScalingTest/XdgPoisson3";
                    C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
                    PrepareControl(C, "XdgPoisson3");
                    return C;
                },
                RunControl = control => RunControl<BoSSS.Application.XdgPoisson3.XdgPoisson3Control, BoSSS.Application.XdgPoisson3.XdgPoisson3Main>((BoSSS.Application.XdgPoisson3.XdgPoisson3Control)control),
            },
            new ThreadScalingTestcase {
                Name = "XNSERisingBubble",
                CreateControl = () => {
                    var C = XNSE_Solver_MPItest.RisingBubble(p: 2, kelem: 18);
                    C.NoOfTimesteps = 4;
                    C.ProjectName = "MultithreadingScalingTest/Bubble";
                    PrepareControl(C, "RisingBubble");
                    return C;
                },
                RunControl = control => RunControl<XNSE_Control, XNSE_Solver.XNSE>((XNSE_Control)control),
            },
        };

        private static void PrepareControl<AppControlType>(AppControlType C, string testcase) where AppControlType : AppControl {
            C.ImmediatePlotPeriod = -1;
            C.SuperSampling = 0;
            C.savetodb = false;
            C.SessionName = $"ThreadPerf_{testcase}";
        }


        private static double RunControl<AppControlType, ApplicationType>(AppControlType control) where AppControlType : AppControl, new() where ApplicationType : Application<AppControlType>, new() {
            ApplicationType solver = null;
            double totalSeconds = 0;

            try {
                ilPSP.Tracing.Tracer.Reset();
                solver = new ApplicationType();
                solver.Init(control);
         

                Stopwatch total = Stopwatch.StartNew();
                solver.RunSolverMode();
                total.Stop();

                var json = solver.OnlineProfiling.Serialize();
                File.WriteAllText($"profiling_bin.json", json);

                totalSeconds = total.Elapsed.TotalSeconds;
            } catch(Exception e) {
                Assert.Fail($"Exception during solver run: {e}");
            } finally {
                solver?.Dispose();
            }

            return totalSeconds;
        }
    }
}