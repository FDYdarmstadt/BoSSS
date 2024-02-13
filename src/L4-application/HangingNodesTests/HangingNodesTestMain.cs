using BoSSS.Application.XNSFE_Solver;
using System;
using System.Linq;
using System.Collections.Generic;
using BoSSS.Foundation.IO;
using MPI.Wrappers;
using System.IO;
using NUnit.Framework;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using ilPSP;

namespace HangingNodesTests {

    /// <summary>
    /// Tests correct construction of quadrules for various setups, with hanging nodes, specially handled double cut cells and grid partitioning
    /// This is a technical test specifically designed to verify the techniques used in the heated wall simulation for SFB1194 K26
    /// </summary>
    [TestFixture]    
    public class HangingNodesTestMain {
        static void Main(string[] args) {
            // mpiexec -n 2 dotnet HangingNodesTests.dll
            Console.WriteLine("Starting Hanging Nodes Test!");
            BoSSS.Solution.Application.InitMPI();

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
                                CheckLengthScales(solver, "sz" + size + "ph" + phase + "setup" + s);
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

            Console.WriteLine("Finished Hanging Nodes Test.");
            Console.WriteLine();
            Console.WriteLine("Results:");
            for(int i = 0; i < Description.Count; i++) {
                Console.WriteLine(Description[i] + " : MomRes : {0}, TempRes : {1}", MomentumRes[i], TemperatureRes[i]);
            }

            Assert.IsTrue(MomentumRes.Select(s => Math.Abs(s)).Max() < 1e-6);
            Assert.IsTrue(TemperatureRes.Select(s => Math.Abs(s)).Max() < 1e-6);

            BoSSS.Solution.Application.FinalizeMPI();
        }

        /// <summary>
        /// single phase
        /// </summary>
        [Test]
        public static void Test1Phase() {
            double[] sizes = new double[] { 1e0, 1e-3 };
            byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            RunTest(sizes, setup, 1);
        }

        /// <summary>
        /// two phases (e.g. air and water)
        /// </summary>
        [Test]
        public static void Test2Phase() {
            double[] sizes = new double[] { 1e0, 1e-3 };
            byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            RunTest(sizes, setup, 2);
        }

        /// <summary>
        /// three phases 
        /// </summary>
        [Test]
        public static void Test3Phase() {
            double[] sizes = new double[] { 1e0, 1e-3 };
            byte[] setup = new byte[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            RunTest(sizes, setup, 3);
        }

        /*
        [Test]
        public static void SerialParallelLengthScaleComp() {
            // --test=HangingNodesTests.HangingNodesTestMain.SerialParallelLengthScaleComp

             // to test individual setups
            double size = 1e0;
            byte s = 0;
            int phase = 3;

            bool plot = true;

            csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int procs);
            csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int rank);

           

            List<double> TemperatureRes = new List<double>();
            List<double> MomentumRes = new List<double>();
            List<string> Description = new List<string>();

            {
                string desc = String.Format("Size : {0}, Phases : {1}, Setup : {2}, Procs : {3}", size, phase, s, procs);
                Description.Add(desc);
                var C = HangingNodesTests.Control.TestSkeleton(size);
                HangingNodesTests.Control.SetAMR(C, size, s);
                HangingNodesTests.Control.SetLevelSet(C, size, phase);
                HangingNodesTests.Control.SetParallel(C, procs == 2 ? -procs : procs);
                if(plot) {
                    C.ImmediatePlotPeriod = 1;
                    C.SuperSampling = 3;
                }

                using(var solver = new XNSFE()) {

                    solver.Init(C);
                    solver.RunSolverMode();
                    
                    CheckLengthScales(solver, "sz" + size + "ph" + phase + "setup" + s);

                    MomentumRes.Add(solver.CurrentResidual.Fields.Take(3).Sum(f => f.L2Norm()).MPISum());
                    TemperatureRes.Add(solver.CurrentResidual.Fields[3].L2Norm().MPISum());
                }
            }
             
            Console.WriteLine("Finished Hanging Nodes Test.");
            Console.WriteLine();
            Console.WriteLine("Results:");
            for(int i = 0; i < Description.Count; i++) {
                Console.WriteLine(Description[i] + " : MomRes : {0}, TempRes : {1}", MomentumRes[i], TemperatureRes[i]);
            }

            Assert.IsTrue(MomentumRes.Select(s => Math.Abs(s)).Max() < 1e-6);
            Assert.IsTrue(TemperatureRes.Select(s => Math.Abs(s)).Max() < 1e-6);
        }
        */

        /// <summary>
        /// Length scale comparison for parallel vs. serial run using the <see cref="TestingIO"/>
        /// </summary>
        /// <remarks>
        /// Note: adaptive mesh refinement (which is used in this test) may produce different GlobalID's 
        /// for otherwise equivalent grids when running in parallel.
        /// Therefore, 
        /// - we cannot compare the GlobalID
        /// - instead, we must compare cell values based on the location  code (<see cref="BoSSS.Platform.Utils.Geom.GeomBinTreeBranchCode"/>).
        /// </remarks>
        private static void CheckLengthScales(XNSFE solver, string filename) {
            
            var species = solver.LsTrk.SpeciesIdS.ToArray();
            var Tracker = solver.LsTrk;
            var agg = Tracker.GetAgglomerator(species, solver.QuadOrder(), solver.Control.AgglomerationThreshold);
            int J = solver.GridData.iLogicalCells.NoOfLocalUpdatedCells;

            /*
            int MPIsize = solver.MPISize;
            int jCellStrange = -1;
           
            for(int j = 0; j < J; j++) {
                var Center_j = solver.GridData.iLogicalCells.GetCenter(j);
                var strangeCenter = new Vector(-0.125, 0.125);

                if((Center_j - strangeCenter).L2Norm() < 1.0e-8)
                    jCellStrange = j;
            }

            if(jCellStrange >= 0) {
                var phiDGcoord =  ((LevelSet)(Tracker.LevelSets[0])).Coordinates.GetRow(jCellStrange);
                var ot = phiDGcoord.ToConcatString("(", "|", ")");
                Console.Error.WriteLine($"Cell X={solver.GridData.iLogicalCells.GetCenter(jCellStrange)}, locIdx={jCellStrange} @ rnk{solver.MPIRank}, surf = {agg.NonAgglomeratedMetrics.CellSurface[Tracker.GetSpeciesId("B")][jCellStrange]}, bkgrCellVol = {solver.GridData.iGeomCells.GetCellVolume(jCellStrange)}, phi = {ot}");
            }
            */

            var LsChecker = new TestingIO(solver.GridData, "CellMetrics-" + filename + ".abc", false, 1);
            for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                SpeciesId spc = species[iSpc];
                string SpcName = Tracker.GetSpeciesName(spc);

                LsChecker.AddVector("LenScale-" + SpcName, agg.CellLengthScales[spc].To1DArray().Take(J).Select(a => double.IsNaN(a) ? -1.1 : a));
                LsChecker.AddVector("Vol-" + SpcName, agg.CutCellVolumes[spc].To1DArray().Take(J).Select(a => double.IsNaN(a) ? -1.2 : a));
                LsChecker.AddVector("Surf-" + SpcName, agg.CellSurface[spc].To1DArray().Take(J).Select(a => double.IsNaN(a) ? -1.3 : a));
            }

            LsChecker.DoIOnow();
            var err = LsChecker.AllRelErr();
            foreach(var kv in err) {
                if(kv.Key != "GlobalID")
                    Console.WriteLine($"    Cell Metric Comparison Error for {kv.Key} = {kv.Value}");
            }
            foreach(var kv in err) {
                if(kv.Key != "GlobalID")
                    Assert.LessOrEqual(kv.Value, 1e-10, $"Cell Metric Comparison Error for {kv.Key} = {kv.Value}, this is to high!");
            }
            
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
                    var C = Control.TestSkeleton(size);
                    Control.SetAMR(C, size, s);
                    Control.SetLevelSet(C, size, phase);
                    Control.SetParallel(C, procs);

                    using (var solver = new XNSFE()) {
                        try {
                            solver.Init(C);
                            solver.RunSolverMode();
                            MomentumRes.Add(solver.CurrentResidual.Fields.Take(3).Sum(f => f.L2Norm()).MPISum());
                            TemperatureRes.Add(solver.CurrentResidual.Fields[3].L2Norm().MPISum());
                            CheckLengthScales(solver, "sz" + size + "ph" + phase + "setup" + s);
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
                        var C = Control.TestSkeleton(size);
                        Control.SetAMR(C, size, s);
                        Control.SetLevelSet(C, size, phase);
                        Control.SetParallel(C, -procs);

                        using (var solver = new XNSFE()) {
                            try {
                                solver.Init(C);
                                solver.RunSolverMode();
                                MomentumRes.Add(solver.CurrentResidual.Fields.Take(3).Sum(f => f.L2Norm()).MPISum());
                                TemperatureRes.Add(solver.CurrentResidual.Fields[3].L2Norm().MPISum());
                                //CheckLengthScales(solver, "sz" + size + "ph" + phase + "setup" + s);
                                
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
                Console.WriteLine($"{Description[i]}: MomRes : {MomentumRes[i]}, TempRes : {TemperatureRes[i]}");
            }
            for (int i = 0; i < MomentumRes.Count; i++) {
                Assert.Less(MomentumRes[i].Abs(), 1e-6, "Momentum Residual to high.");
            }
            for (int i = 0; i < TemperatureRes.Count; i++) {
                Assert.Less(TemperatureRes[i].Abs(), 1e-6, "Temperature Residual to high.");
            }
        }

    }
}
