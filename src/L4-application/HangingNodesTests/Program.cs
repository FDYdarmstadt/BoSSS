﻿using BoSSS.Application.XNSFE_Solver;
using System;
using System.Linq;
using System.Collections.Generic;
using BoSSS.Foundation.IO;
using MPI.Wrappers;
using System.IO;
using NUnit.Framework;
using System.Diagnostics;

namespace HangingNodesTests {

    /// <summary>
    /// Tests correct construction of quadrules for various setups, with hanging nodes, specially handled double cut cells and grid partitioning
    /// This is a technical test specifically designed to varify the techniques used in the heated wall simulation for SFB1194 K26
    /// </summary>
    [TestFixture]    
    public class Program {
        static void Main(string[] args) {
            // mpiexec -n 2 dotnet HangingNodesTests.dll
            Console.WriteLine("Starting Hanging Nodes Test!");
            BoSSS.Solution.Application.InitMPI();

            //double[] sizes = new double[] { 1e0 };
            //byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            //int[] phases = new int[] { 1, 2, 3 };

            // to test individual setups
            double[] sizes = new double[] { 1e0 };
            byte[] setup = new byte[] { 0 };
            int[] phases = new int[] { 3 };

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
                        HangingNodesTests.Control.SetParallel(C, procs == 2 ? -procs : procs);

                        if (plot) {
                            C.ImmediatePlotPeriod = 1;
                            C.SuperSampling = 3;                            
                        }

                        using (var solver = new XNSFE()) {
                            try {
                                solver.Init(C);
                                solver.RunSolverMode();
                                if (plot) {
                                    var MultiphaseAgglomerator = solver.LsTrk.GetAgglomerator(solver.LsTrk.SpeciesIdS.ToArray(), solver.QuadOrder(), C.AgglomerationThreshold);    
                                    foreach(var spc in solver.LsTrk.SpeciesIdS) {
                                        string spcName = solver.LsTrk.GetSpeciesName(spc);
                                        var speciesAgglomerator = MultiphaseAgglomerator.GetAgglomerator(spc);
                                        speciesAgglomerator.PlotAgglomerationPairs($"agglomerationPairs-{spcName}-MPI{rank}.txt", null, true);
                                    }
           
                                }
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

        [Test]
        public static void Test1Phase() {
            double[] sizes = new double[] { 1e0, 1e-3 };
            byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            RunTest(sizes, setup, 1);
        }

        [Test]
        public static void Test2Phase() {
            double[] sizes = new double[] { 1e0, 1e-3 };
            byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            RunTest(sizes, setup, 2);
        }

        [Test]
        public static void Test3Phase() {
            double[] sizes = new double[] { 1e0, 1e-3 };
            byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            RunTest(sizes, setup, 3);
        }

        private static void RunTest(double[] sizes, byte[] setup, int phase) {

            csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int procs);
            csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int rank);

            List<double> TemperatureRes = new List<double>();
            List<double> MomentumRes = new List<double>();
            List<string> Description = new List<string>();

            foreach (double size in sizes) {
                foreach (byte s in setup) {
                    string desc = String.Format("Size : {0}, Phases : {1}, Setup : {2}, Procs : {3}", size, phase, s, procs);
                    Description.Add(desc);
                    var C = HangingNodesTests.Control.TestSkeleton(size);
                    HangingNodesTests.Control.SetAMR(C, size, s);
                    HangingNodesTests.Control.SetLevelSet(C, size, phase);
                    HangingNodesTests.Control.SetParallel(C, procs);

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

            if (procs == 2) {
                foreach (double size in sizes) {
                    foreach (byte s in setup) {
                        string desc = String.Format("Size : {0}, Phases : {1}, Setup : {2}, Procs (transpose) : {3}", size, phase, s, procs);
                        Description.Add(desc);
                        var C = HangingNodesTests.Control.TestSkeleton(size);
                        HangingNodesTests.Control.SetAMR(C, size, s);
                        HangingNodesTests.Control.SetLevelSet(C, size, phase);
                        HangingNodesTests.Control.SetParallel(C, -procs);

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

            Console.WriteLine("Finished Hanging Nodes Test with {0} procs.", procs);
            Console.WriteLine();
            Console.WriteLine("Results:");
            for (int i = 0; i < Description.Count; i++) {
                Console.WriteLine(Description[i] + " : MomRes : {0}, TempRes : {1}", MomentumRes[i], TemperatureRes[i]);
            }
            Assert.IsTrue(MomentumRes.Select(s => Math.Abs(s)).Max() < 1e-6);
            Assert.IsTrue(TemperatureRes.Select(s => Math.Abs(s)).Max() < 1e-6);
        }

    }
}
