using BoSSS.Application.XNSFE_Solver;
using System;
using System.Linq;
using System.Collections.Generic;
using BoSSS.Foundation.IO;
using MPI.Wrappers;
using System.IO;

namespace HangingNodesTests {
    class Program {
        static void Main(string[] args) {
            Console.WriteLine("Starting Hanging Nodes Test!");
            BoSSS.Solution.Application.InitMPI();

            //double[] sizes = new double[] { 1e0, 1e-3 };
            //byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            //int[] phases = new int[] { 1, 2, 3 };

            // to test individual setups
            double[] sizes = new double[] { 1e0 };
            byte[] setup = new byte[] { 0 };
            int[] phases = new int[] { 2 };

            bool plot = true;

            csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int procs);
            csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int rank);

            if (rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());

                Console.Write("rm");
                foreach (var pltFile in dir.GetFiles("*.txt").Concat(dir.GetFiles("*.csv"))) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
                Console.Write("rm");
                foreach (var pltFile in dir.GetFiles("*.plt")) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }



            List<double> TemperatureRes = new List<double>();
            List<double> MomentumRes = new List<double>();
            List<string> Description = new List<string>();

            foreach (double size in sizes) {
                foreach(int phase in phases) {
                    foreach(byte s in setup) {
                        string desc = String.Format("Size : {0}, Phases : {1}, Setup : {2}, Procs : {3}", size, phase, s, procs);
                        Description.Add(desc);
                        var C = HangingNodesTests.Control.TestSkeleton(size);
                        HangingNodesTests.Control.SetAMR(C, size, s);
                        HangingNodesTests.Control.SetLevelSet(C, size, phase);
                        HangingNodesTests.Control.SetParallel(C, procs);

                        if (plot) {
                            C.ImmediatePlotPeriod = 1;
                            C.SuperSampling = 3;
                        }

                        using (var solver = new XNSFE()) {
                            try {
                                solver.Init(C);
                                solver.RunSolverMode();
                                MomentumRes.Add(solver.CurrentResidual.Fields.Take(3).Sum(f => f.L2Norm()).MPISum());
                                TemperatureRes.Add(solver.CurrentResidual.Fields[3].L2Norm().MPISum());
                            } catch (Exception e) {
                                Console.WriteLine(desc + " : failed");
                                Console.WriteLine(e.Message);
                                Console.WriteLine(e.StackTrace);
                                TemperatureRes.Add(-1.0);
                                MomentumRes.Add(-1.0);
                            }
                        }
                       
                    }
                }
            }

            BoSSS.Solution.Application.FinalizeMPI();
            Console.WriteLine("Finished Hanging Nodes Test.");
            Console.WriteLine();
            Console.WriteLine("Results:");
            for(int i = 0; i < Description.Count; i++) {
                Console.WriteLine(Description[i] + " : MomRes : {0}, TempRes : {1}", MomentumRes[i], TemperatureRes[i]);
            }
        }
    }
}
