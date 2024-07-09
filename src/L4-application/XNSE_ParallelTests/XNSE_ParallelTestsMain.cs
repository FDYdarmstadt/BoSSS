using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;
using System.Drawing;
using BoSSS.Application.XNSE_Solver;
using XNSE_ParallelTets;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;

namespace XNSE_ParallelTests {

    [TestFixture]
    public class XNSE_ParallelTests {

        [Test]
        public static void Test_ChannelFlow2D() {
            var C = Controls.Test_ChannelFlow2D(true, false);
            RunTest(C, "ChannelFlow2D_baseCase");
        }

        [Test]
        public static void Test_ChannelFlow3D() {
            var C = Controls.Test_ChannelFlow3D(true, false);
            RunTest(C, "ChannelFlow3D_baseCase");
        }


        static void Main(string[] args) {

            // to test individual setups
            //var C = Controls.Test_ChannelFlow2D(true, false);
            var C = Controls.Test_ChannelFlow3D(true, false);

            RunTest(C, "localTestcase");

        }


        private static void RunTest(XNSE_Control control, string testcaseName) {

            //System.Environment.SetEnvironmentVariable("OMP_NUM_THREADS", "1");
            BoSSS.Solution.Application.InitMPI(num_threads: 1);

            csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int procs);
            csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int rank);

            Console.WriteLine($"Running {testcaseName} on {procs} procs");

            List<double> MomentumRes = new List<double>();

            using (var solver = new XNSE()) {
                try {
                    solver.Init(control);
                    solver.RunSolverMode();
                    int D = solver.Grid.SpatialDimension;
                    MomentumRes.Add(solver.CurrentResidual.Fields.Take(D).Sum(f => f.L2Norm()).MPISum());
                    CompareParallelRun(solver, testcaseName);
                } catch (Exception e) {
                    Console.WriteLine($"run on {procs} procs failed");
                    Console.WriteLine(e.Message);
                    Console.WriteLine(e.StackTrace);
                    MomentumRes.Add(-1.0);
                }
            }

            for (int i = 0; i < MomentumRes.Count; i++) {
                Assert.Less(MomentumRes[i].Abs(), 1e-6, "Momentum Residual too high.");
            }
        }

        private static void CompareParallelRun(XNSE solver, string testcaseName) {

            int RefMPIsize = 1;

            var FieldChecker = new TestingIO(solver.GridData, $"SolutionsFields_{testcaseName}.csv", true, RefMPIsize);
            int D = solver.GridData.SpatialDimension;
            for (int d = 0; d < D; d++) {
                FieldChecker.AddDGField(solver.Velocity[d]);
            }
            FieldChecker.AddDGField(solver.Pressure);
            FieldChecker.DoIOnow();

            for (int d = 0; d < D; d++) {
                Console.WriteLine($"absolute L2 error for {solver.Velocity[d].Identification} field: {FieldChecker.AbsError(solver.Velocity[d])}");
                Assert.Less(FieldChecker.AbsError(solver.Velocity[d]), 1.0e-9, $"Mismatch in velocity{d} field between single-core and parallel run.");
            }
            Console.WriteLine($"absolute L2 error for {solver.Pressure.Identification} field: {FieldChecker.AbsError(solver.Pressure)}");
            Assert.Less(FieldChecker.AbsError(solver.Pressure), 1.0e-9, "Mismatch in pressure field between single-core and parallel run.");
        } 

    }
}
