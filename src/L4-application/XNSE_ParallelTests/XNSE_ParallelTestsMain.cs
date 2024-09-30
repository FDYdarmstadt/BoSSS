using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;
using System.Drawing;
using BoSSS.Application.XNSE_Solver;
using XNSE_ParallelTets;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Grid;
using System.Collections;

namespace XNSE_ParallelTests {

    [TestFixture]
    public class XNSE_ParallelTests {

        [Test]
        public static void Test_ChannelFlow2D() {
            var C = Controls.Test_ChannelFlow2D(false, false);
            RunTest(C, "ChannelFlow2D_baseCase");
        }

        [Test]
        public static void Test_ChannelFlow3D() {
            var C = Controls.Test_ChannelFlow3D(false, false);
            RunTest(C, "ChannelFlow3D_baseCase");
        }


        static void Main(string[] args) {

            // to test individual setups
            var C = Controls.Test_ChannelFlow2D(false, false);
            //var C = Controls.Test_ChannelFlow3D(false, false);

            //C.PlotAgglomeration = true;
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 3;

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
                    //MomentumRes.Add(solver.CurrentResidual.Fields.Take(D).Sum(f => f.L2Norm()).MPISum());
                    //CheckLevelSetProperties(solver);
                    CompareParallelRun(solver, testcaseName);
                } catch (Exception e) {
                    Console.WriteLine($"run on {procs} procs failed");
                    Console.WriteLine(e.Message);
                    Console.WriteLine(e.StackTrace);
                    MomentumRes.Add(-1.0);
                }
            }

            //for (int i = 0; i < MomentumRes.Count; i++) {
            //    Console.WriteLine($"Momentum {i} Residual: {MomentumRes[i].Abs()}");
            //    Assert.Less(MomentumRes[i].Abs(), 1e-6, "Momentum Residual too high.");
            //}
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

            XDGField[] errorFields = new XDGField[D + 1];
            XDGField err;
            for (int d = 0; d < D; d++) {
                Console.WriteLine($"absolute L2 error for {solver.Velocity[d].Identification} field: {FieldChecker.AbsError(solver.Velocity[d])}");
                err = FieldChecker.LocalError(solver.Velocity[d]);
                errorFields[d] = err;
                //Assert.Less(FieldChecker.AbsError(solver.Velocity[d]), 1.0e-9, $"Mismatch in velocity{d} field between single-core and parallel run.");
            }
            Console.WriteLine($"absolute L2 error for {solver.Pressure.Identification} field: {FieldChecker.AbsError(solver.Pressure)}");
            err = FieldChecker.LocalError(solver.Pressure);
            errorFields[D] = err;
            //Assert.Less(FieldChecker.AbsError(solver.Pressure), 1.0e-9, "Mismatch in pressure field between single-core and parallel run.");


            BoSSS.Solution.Tecplot.Tecplot.PlotFields(errorFields, "XNSE_ParallelTests-ErrorFields", 0.0, 3);
        }

        private static void CheckLevelSetProperties(XNSE solver) {

            int D = solver.LsTrk.GridDat.SpatialDimension;

            double radius = 0.4;

            Func<double[], double> LSfunc = X => -1;
            double volumeDrop = 0.0;
            double surfaceAreaDrop = 0.0;
            if (D == 2) {
                LSfunc = (X => ((X[0] - 1.0).Pow2() + (X[1] - 1.0).Pow2()) - radius.Pow2());
                volumeDrop = Math.PI * radius.Pow(2);
                surfaceAreaDrop = 2.0 * Math.PI * radius;
            }

            if (D == 3) {
                LSfunc = (X => ((X[0] - 1.0).Pow2() + (X[1] - 1.0).Pow2() + (X[2] - 1.0).Pow2()) - radius.Pow2());
                volumeDrop = (4.0 / 3.0) * Math.PI * radius.Pow(3);
                surfaceAreaDrop = 4.0 * Math.PI * radius.Pow(2);
            }
   
            double error = ((SinglePhaseField)solver.LsTrk.LevelSets[0]).L2Error(LSfunc.Vectorize(), solver.QuadOrder());
            Console.WriteLine($"LS projection error: {error}");

            double volume = XNSEUtils.GetSpeciesArea(solver.LsTrk, solver.LsTrk.GetSpeciesId("A"), solver.QuadOrder());
            Console.WriteLine($"droplet volume: {volume} (error compared to analytic sphere: {volume - volumeDrop})");

            double surfaceArea = XNSEUtils.GetInterfaceLength(solver.LsTrk, solver.QuadOrder());
            Console.WriteLine($"droplet surface: {surfaceArea} (error compared to analytic sphere: {surfaceArea - surfaceAreaDrop})");


        }
    }
}
