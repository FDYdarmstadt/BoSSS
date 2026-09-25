using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using NUnit.Framework;
using System;
using System.Diagnostics;
using System.IO;
using System.Threading;
using static ilPSP.Tracing.Tracer;

namespace BoSSS.Application.MultithreadingSolverTest {

    internal sealed class ThreadScalingRequirement {
        /// <summary>
        /// Thread count used for the baseline measurement.
        /// </summary>
        public int FromThreads;

        /// <summary>
        /// Thread count used for the scaled measurement.
        /// </summary>
        public int ToThreads;

        /// <summary>
        /// Minimum required speedup: T(FromThreads) / T(ToThreads).
        /// Zero or a negative value disables this requirement.
        /// </summary>
        public double MinSpeedup;
    }

    internal sealed class ThreadScalingMethod {
        /// <summary>
        /// Extra-P method-name suffix or full region name. All matching invocation
        /// paths are aggregated before calculating the speedup.
        /// </summary>
        public string Name;

        /// <summary>
        /// Explicit thread-pair requirements for this method. An empty list disables
        /// all method checks.
        /// </summary>
        public ThreadScalingRequirement[] Requirements = Array.Empty<ThreadScalingRequirement>();
    }

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
        /// Explicit thread-pair requirements for total runtime. An empty list disables
        /// all total-runtime checks.
        /// </summary>
        public ThreadScalingRequirement[] Requirements = Array.Empty<ThreadScalingRequirement>();

        /// <summary>
        /// Extra-P methods to validate in addition to total runtime. Each method has
        /// independent speedup requirements.
        /// </summary>
        public ThreadScalingMethod[] ScalingMethods = Array.Empty<ThreadScalingMethod>();

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
                    var C =  BoSSS.Application.SipPoisson.SipHardcodedControl.Square(96,96);
                    C.ProjectName = "MultithreadingScalingTest/ipPoisson";
                    C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
                    PrepareControl(C, "ipPoisson");
                    
                    return C;
                },
                Requirements = new[] {
                    new ThreadScalingRequirement { FromThreads = 1, ToThreads = 2, MinSpeedup = 1.2 },
                    new ThreadScalingRequirement { FromThreads = 1, ToThreads = 4, MinSpeedup = 1.3 },
                    new ThreadScalingRequirement { FromThreads = 2, ToThreads = 4, MinSpeedup = 1.0 }
                },
                ScalingMethods = new[] {
                    new ThreadScalingMethod {
                        Name = "AssembleMatrix",
                        Requirements = new[] {
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 2, MinSpeedup = 1.4 },
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 4, MinSpeedup = 1.3 }
                        }
                    },
                    new ThreadScalingMethod {
                        Name = "BoSSS.Foundation.DifferentialOperator+EvaluatorNonLin.Evaluate",
                        Requirements = new[] {
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 2, MinSpeedup = 1.3 }
                        }
                    }
                },
                RunControl = control => RunControl<BoSSS.Application.SipPoisson.SipControl, BoSSS.Application.SipPoisson.SipPoissonMain>((BoSSS.Application.SipPoisson.SipControl)control),
            },
            new ThreadScalingTestcase {
                // this is not really good now
                Name = "XdgPoisson3",
                CreateControl = () => {
                    var C = BoSSS.Application.XdgPoisson3.HardCodedControl.Circle(Resolution: 192, p: 4);
                    C.ProjectName = "MultithreadingScalingTest/XdgPoisson3";
                    C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
                    PrepareControl(C, "XdgPoisson3");
                    return C;
                },
                ScalingMethods = new[] {
                    new ThreadScalingMethod {
                        Name = "AssembleMatrix",
                        Requirements = new[] {
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 2, MinSpeedup = 1.1 }
                        }
                    }
                },
                RunControl = control => RunControl<BoSSS.Application.XdgPoisson3.XdgPoisson3Control, BoSSS.Application.XdgPoisson3.XdgPoisson3Main>((BoSSS.Application.XdgPoisson3.XdgPoisson3Control)control),
            },
            new ThreadScalingTestcase {
                Name = "XNSERisingBubble",
                CreateControl = () => {
                    var C = XNSE_Solver_MPItest.RisingBubble(p: 2, kelem: 64);
                    C.NoOfTimesteps = 1;
                    C.ProjectName = "MultithreadingScalingTest/Bubble";
                    PrepareControl(C, "RisingBubble");
                    return C;
                },
                ScalingMethods = new[] {
                    new ThreadScalingMethod {
                        Name = "AssembleMatrixCallback",
                        Requirements = new[] {
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 2, MinSpeedup = 1.3 },
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 4, MinSpeedup = 1.4 }
                        }
                    },
                    new ThreadScalingMethod {
                        Name = "EvaluateOperator",
                        Requirements = new[] {
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 2, MinSpeedup = 1.3 },
                            new ThreadScalingRequirement { FromThreads = 1, ToThreads = 4, MinSpeedup = 1.4 }
                        }
                    }
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
                solver = new ApplicationType();
                solver.Init(control);

                string profilingJson;
    
                using(IsolatedTraceRootScope traceScope = ilPSP.Tracing.Tracer.IsolateRoot("root_frame")) {
                    Stopwatch total = Stopwatch.StartNew();
                    solver.RunSolverMode();
                    total.Stop();
                    traceScope.Root.UpdateTime();
                    profilingJson = solver.OnlineProfiling.SerializeWithCustomMCR(traceScope.Root);
                    totalSeconds = total.Elapsed.TotalSeconds;
                }

                File.WriteAllText("profiling_bin.json", profilingJson);

            } catch(Exception e) {
                Assert.Fail($"Exception during solver run: {e}");
            } finally {
                solver?.Dispose();
            }

            return totalSeconds;
        }
    }
}
