/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Queries;
using CNS;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.LoadBalancing;
using CNS.ShockCapturing;
using CNS.Tests;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace HilbertTest {

    /// <summary>
    /// 
    /// </summary>
    [TestFixture]
    public class HilbertTest : TestProgram<CNSControl> {

        public static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            Test();
            BoSSS.Solution.Application.FinalizeMPI();
        }

        public HilbertTest(){
            AssertWatch = new Stopwatch();
            AssertWatch.Start();
        }

        public override void Dispose() {
            CheckAssertWatch();
            AssertWatch = null;
        }

        private Stopwatch AssertWatch;

        private void CheckAssertWatch() {
            AssertWatch.Stop();
            double time = AssertWatch.Elapsed.TotalSeconds;
            //Assert.IsTrue(time < 60,"time limit of 60 seconds exceeded. There is something rotten, plz check ...");
        }

        public static void Test() {
            TestingCoordinateSamples();                 
            TestingGridDistributionEven();           
            TestingGridDistributionUneven();           
            TestingdirectHilbertEven();
            TestingdirectHilbertUneven();
            TestingGridDistributionDynamic();
            TestingReclusteringIndependency();
            ClusterHilbertTest();
        }


        /// <summary>
        /// Testing partitioning with 1 cluster, even distribution of cells among processes
        /// To ensure the result is the same as with direct Hilbert and equal costs everywhere
        /// </summary>
        [NUnitFileToCopyHack("HilbertTest/Tests.zip")]
        [Test]
        static public void TestingGridDistributionEven() {
            string dbPath = @"Tests.zip";
            //TestCase: 4x4 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_PartTest(dbPath, "7ac582f5-8913-439b-9f2b-9fbf96141d76", "b7793aee-44b6-44c7-91e7-5debd7f44c3b", 4, 4);

            using (var solver = new HilbertTest()) {

                solver.Init(control);
                solver.RunSolverMode();
                bool result = true;

                int Jloc = solver.GridData.CellPartitioning.LocalLength;
                for (int j = 0; j < Jloc; j++) {
                    double[] XC = solver.GridData.iLogicalCells.GetCenter(j);
                    double xC = XC[0];
                    double yC = XC[1];
                    switch (solver.MPIRank) {
                        case 0:
                        result &= (xC > 0) && (xC < 0.5) && (yC > 0) && (yC < 0.5);
                        break;
                        case 1:
                        result &= (xC > 0.5) && (xC < 1) && (yC > 0) && (yC < 0.5);
                        break;
                        case 2:
                        result &= (xC > 0.5) && (xC < 1) && (yC > 0.5) && (yC < 1);
                        break;
                        case 3:
                        result &= (xC > 0) && (xC < 0.5) && (yC > 0.5) && (yC < 1);
                        break;
                    }
                }
                Console.WriteLine("Test Grid Distribution even");
                Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
                Assert.IsTrue(result.MPIAnd(), "HilbertCurve or mapping (rank->Hilbertcurve) is corrupted");
            }
        }

        /// <summary>
        /// Testing bare Hilbert partitioning, even distribution of cells among processes
        /// </summary>
        [NUnitFileToCopyHack("HilbertTest/Tests.zip")]
        [Test]
        static public void TestingdirectHilbertEven() {
            //--test=HilbertTest.HilbertTest.TestingdirectHilbertEven
            //string dbPath = @"D:\Weber\BoSSS\test_db";
            string dbPath = @"Tests.zip";
            //TestCase: 4x4 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_directHilbert(dbPath, "7ac582f5-8913-439b-9f2b-9fbf96141d76", "b7793aee-44b6-44c7-91e7-5debd7f44c3b", 4, 4);

            using (var solver = new HilbertTest()) {
                solver.Init(control);
                solver.RunSolverMode();
                bool result = true;

                int Jloc = solver.GridData.CellPartitioning.LocalLength;
                for (int j = 0; j < Jloc; j++) {
                    double[] XC = solver.GridData.iLogicalCells.GetCenter(j);
                    double xC = XC[0];
                    double yC = XC[1];
                    switch (solver.MPIRank) {
                        case 0:
                        result &= (xC > 0) && (xC < 0.5) && (yC > 0) && (yC < 0.5);
                        break;
                        case 1:
                        result &= (xC > 0.5) && (xC < 1) && (yC > 0) && (yC < 0.5);
                        break;
                        case 2:
                        result &= (xC > 0.5) && (xC < 1) && (yC > 0.5) && (yC < 1);
                        break;
                        case 3:
                        result &= (xC > 0) && (xC < 0.5) && (yC > 0.5) && (yC < 1);
                        break;
                    }
                }
                Console.WriteLine("Test Hilbert: Grid Distribution even");
                Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
                Assert.IsTrue(result.MPIAnd(), "HilbertCurve or mapping of Hilbert (rank->Hilbertcurve) is corrupted");
            }
        }

        /// <summary>
        /// Testing partitioning with clusters, uneven distribution of cells among processes
        /// To ensure the result is the same as with direct Hilbert and equal costs everywhere
        /// </summary>
        [NUnitFileToCopyHack("HilbertTest/Tests.zip")]
        [Test]
        static public void TestingGridDistributionUneven() {
            string dbPath = @"Tests.zip";
            //TestCase: 3x3 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_PartTest(dbPath, "ccb23f25-04e9-467a-b667-bb3d642b6447", "9b24a2e6-2ce5-4de2-bd08-37a930f0df06", 3, 3);
            using (var solver = new HilbertTest()) {
                solver.Init(control);
                solver.RunSolverMode();
                bool result = true;

                int Jloc = solver.GridData.CellPartitioning.LocalLength;
                for (int j = 0; j < Jloc; j++) {
                    double[] XC = solver.GridData.iLogicalCells.GetCenter(j);
                    double xC = XC[0];
                    double yC = XC[1];
                    switch (solver.MPIRank) {
                        case 0:
                        result &= (xC > 0) && (xC < 0.33) && (yC > 0) && (yC < 0.67);
                        break;
                        case 1:
                        result &= (xC > 0.33) && (xC < 0.67) && (yC > 0) && (yC < 0.67);
                        break;
                        case 2:
                        result &= (xC > 0.67) && (xC < 1) && (yC > 0) && (yC < 0.67);
                        break;
                        case 3:
                        result &= (xC > 0) && (xC < 1) && (yC > 0.67) && (yC < 1);
                        break;
                    }
                }
                Console.WriteLine("Test Grid Distribution uneven");
                Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
                Assert.IsTrue(result.MPIAnd(), "Distribution pattern along HilbertCurve is corrupted");
            }
        }

        /// <summary>
        /// Testing bare Hilbert partitioning, uneven distribution of cells among processes
        /// </summary>
        [NUnitFileToCopyHack("HilbertTest/Tests.zip")]
        [Test]
        static public void TestingdirectHilbertUneven() {
            string dbPath = @"Tests.zip";
            //TestCase: 3x3 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_directHilbert(dbPath, "ccb23f25-04e9-467a-b667-bb3d642b6447", "9b24a2e6-2ce5-4de2-bd08-37a930f0df06", 3, 3);
            using (var solver = new HilbertTest()) {
                solver.Init(control);
                solver.RunSolverMode();
                bool result = true;

                int Jloc = solver.GridData.CellPartitioning.LocalLength;
                for (int j = 0; j < Jloc; j++) {
                    double[] XC = solver.GridData.iLogicalCells.GetCenter(j);
                    double xC = XC[0];
                    double yC = XC[1];
                    switch (solver.MPIRank) {
                        case 0:
                        result &= (xC > 0) && (xC < 0.33) && (yC > 0) && (yC < 0.67);
                        break;
                        case 1:
                        result &= (xC > 0.33) && (xC < 0.67) && (yC > 0) && (yC < 0.67);
                        break;
                        case 2:
                        result &= (xC > 0.67) && (xC < 1) && (yC > 0) && (yC < 0.67);
                        break;
                        case 3:
                        result &= (xC > 0) && (xC < 1) && (yC > 0.67) && (yC < 1);
                        break;
                    }
                }
                Console.WriteLine("Test Hilbert: Grid Distribution uneven");
                Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
                Assert.IsTrue(result.MPIAnd(), "Distribution pattern along HilbertCurve is corrupted");
            }
        }

        /// <summary>
        /// Testing Partition with Constraints (LTS), even distribution among processes.
        /// The test result is also valid for other Constraint-types, e.g. AV, but harder to test
        /// </summary>
        [Test]
        static public void TestingGridDistributionDynamic() {
            //string dbPath = @"D:\Weber\BoSSS\test_db";
            //TestCase: 5x4 grid, Timesteps, LTS-Cluster, PartOn,recInt,AV=false, dgdegree=0, Timestepping=LTS
            CNSControl control = ShockTube_PartTest_Dynamic(5, 4, 1, 2, true, 1);
            using (var solver = new HilbertTest()) {
                solver.Init(control);
                solver.RunSolverMode();
                bool result = false;

                long[] Gid_local = solver.GridData.CurrentGlobalIdPermutation.Values.CloneAs();
                long[] GidGl = Gid_local.MPIAllGatherv((new int[] { Gid_local.Length }).MPIAllGatherv());

                List<DGField> listOfDGFields = (List<DGField>)solver.IOFields;
                DGField field = listOfDGFields[12];
                int D = field.GridDat.SpatialDimension;
                int J = field.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                //intention:Checking if BoundaryBox of LTSCluster==1 is as expected
                //Therefore Computing BoundaryBox of LTSCluster==1
                var BB = new BoSSS.Platform.Utils.Geom.BoundingBox(D);
                var CellBB = new BoSSS.Platform.Utils.Geom.BoundingBox(D);
                for (int i = 0; i < J; i++) {

                    if (field.GetMeanValue(i) == 1) {
                        //Cell cj=solver.GridData.Cells.GetCell(i);
                        solver.GridData.iLogicalCells.GetCellBoundingBox(i, CellBB);
                        BB.AddBB(CellBB);
                    }
                }
                BB.Max = BB.Max.MPIMax();
                BB.Min = BB.Min.MPIMin();
                for (int i = 0; i < D; i++) {
                    BB.Max[i] = Math.Round(BB.Max[i] * 100) / 100;
                    BB.Min[i] = Math.Round(BB.Min[i] * 100) / 100;
                }
                double[] MaxRef = { 0.6, 1 };
                double[] MinRef = { 0, 0 };

             
                if (ItemsAreEqual(BB.Max, MaxRef) && ItemsAreEqual(BB.Min, MinRef)) {
                    //Comparing checkLTS to Distribution along HilbertCurve of Testcase
                    int[] checkLTS = { 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 2, 2, 2, 3, 3, 3 };
                    //result = ItemsAreEqual(solver.Grid.GetHilbertSortedRanks(),checkLTS);
                    long J0 = solver.GridData.CellPartitioning.i0;
                    long JE = solver.GridData.CellPartitioning.iE;

                    ulong[] discreteCenter = new ulong[D];
                    ulong[] local_HilbertIndex = new ulong[JE - J0];
                    int[] local_RankIndex = new int[JE - J0];

                    var _Grid = (GridCommons)(solver.Grid);

                    for (long j = J0; j < JE; j++) {
                        Cell Cj = _Grid.Cells[j - J0];
                        int NoOfNodes = Cj.TransformationParams.NoOfRows;
                        for (int d = 0; d < D; d++) {
                            double center = 0;
                            for (int k = 0; k < NoOfNodes; k++) {
                                center += Cj.TransformationParams[k, d];
                            }

                            center = center / ((double)NoOfNodes);
                            double centerTrf = center * Math.Pow(2, 32);

                            centerTrf = Math.Round(centerTrf);
                            if (centerTrf < 0)
                                centerTrf = 0;
                            if (centerTrf > ulong.MaxValue)
                                centerTrf = ulong.MaxValue;
                            discreteCenter[d] = (ulong)centerTrf;

                        }
                        ulong iH = ilPSP.HilbertCurve.HilbertCurve.hilbert_c2i(32, discreteCenter);
                        local_HilbertIndex[j - J0] = iH;
                        local_RankIndex[j - J0] = solver.MPIRank;
                    }
                    int[] CellsPerRank = { 5, 5, 5, 5 };
                    ulong[] HilbertIndex = local_HilbertIndex.MPIAllGatherv(CellsPerRank);
                    int[] RankIndex = local_RankIndex.MPIAllGatherv(CellsPerRank);
                    Array.Sort(HilbertIndex, RankIndex);
                    result = ItemsAreEqual(RankIndex, checkLTS);
                } else {
                    //catching error caused by changes to LTS-Clustering
                    Console.WriteLine("Unexpected result for LTS Clusters! Computation of LTS Clusters changed. Test aborted.");
                    result = false;
                }
                Console.WriteLine("Test Grid Distribution Dynamic LTS");
                Console.WriteLine("Testresult: {0}", result);
                Assert.IsTrue(result.MPIAnd(), "Dynamic Distribution along HilbertCurve is corrupted");
            }
        }

        /// <summary>
        /// Comparing simulations with&without Repartitioning while LTS-Reclustering on
        /// </summary>
        [Test]
        static public void TestingReclusteringIndependency() {

            //TestCase: 5x4 grid, AV=false, dgdegree=0, Timestepping=LTS&RK
            CNSControl ctrRepON = ShockTube_PartTest_Dynamic(5, 4, int.MaxValue, 2, true, 5);
            CNSControl ctrRepOFF = ShockTube_PartTest_Dynamic(5, 4, int.MaxValue, 2, false, 5);
            using (var solverRepON = new HilbertTest())
            using (var solverRepOFF = new HilbertTest()) {
                solverRepON.Init(ctrRepON);
                solverRepON.RunSolverMode();
                solverRepOFF.Init(ctrRepOFF);
                solverRepOFF.RunSolverMode();

                bool result = true;
                string[] varname = { "Desity", "x-Momentum", "y-Momentum", "Energy" };
                for (int i = 0; i < 4; i++) {
                    List<DGField> listOfDGFields_RepON = (List<DGField>)solverRepON.IOFields;
                    DGField variableRepON = listOfDGFields_RepON[i];
                    double L2NormRepON = variableRepON.L2Norm();
                    List<DGField> listOfDGFields_RepOFF = (List<DGField>)solverRepOFF.IOFields;
                    DGField vriableRepOFF = listOfDGFields_RepOFF[i];
                    double L2NormRepOFF = vriableRepOFF.L2Norm();

                    Console.WriteLine("{0}-L2Norm Rep ON: {1}", varname[i], L2NormRepON);
                    Console.WriteLine("{0}-L2Norm Rep OFF: {1}", varname[i], L2NormRepOFF);
                    bool normequal = (Math.Abs(L2NormRepON - L2NormRepOFF) <= 1e-14);
                    result &= normequal;
                    Console.WriteLine("{0}-L2Norm equal: {1}", varname[i], normequal);
                }
                Assert.IsTrue(result.MPIAnd(), "Repartitioning effects result!");
            }

        }

        static private bool ItemsAreEqual(double[] item1, double[] item2) {
            for (int i = 0; i < item1.Length; i++) {
                if (item1[i] != item2[i])
                    return false;
            }
            return true;
        }
        static private bool ItemsAreEqual(int[] item1, int[] item2) {
            for (int i = 0; i < item1.Length; i++) {
                if (item1[i] != item2[i])
                    return false;
            }
            return true;
        }


        static private double[] cmpBBCell(double[] old, double[] test, bool type) {
            //choose type=true for maxarg(old,test), false for minarg(old,test)
            double[] output = new double[old.Length];
            double[] cmp1;
            double[] cmp2;
            if (type) {
                cmp1 = old;
                cmp2 = test;
            } else {
                cmp1 = test;
                cmp2 = old;
            }
            for (int i = 0; i < old.Length; i++) {
                if (cmp1[i] > cmp2[i]) {
                    return old;
                }
            }
            return test;
        }

        /// <summary>
        /// Testing coordiante samples. Verification of H-Curve-Code, Coordsamples taken from "Convergence with clusterHilbert's Space Filling Curve" by ARTHUR R. BUTZ, p.133
        /// </summary>
        [Test]
        static public void TestingCoordinateSamples() {
            bool[] testresult = new bool[3];
            bool result = true;
            testresult[0] = test_both(4, 654508, new ulong[] { 4, 12, 6, 12, 12 });
            testresult[1] = test_both(4, 458751, new ulong[] { 8, 8, 0, 15, 0 });
            testresult[2] = test_both(4, 294911, new ulong[] { 7, 0, 15, 15, 0 });
            for (int i = 0; i < testresult.Length; i++) {
                Console.WriteLine("Test_i2c of Coord {0}:{1}", i, testresult[i]);
                result |= result;
            }
            Assert.IsTrue(result, "Code of HilbertCurve is corrupted");
        }

        static private bool test_both(int nBits, ulong index, ulong[] coord) {
            return test_i2c(nBits, index, coord) && test_c2i(nBits, index, coord);
        }

        static private bool test_i2c(int nBits, ulong index, ulong[] checkcoord) {
            int Dim = checkcoord.Length;
            ulong[] coord = new ulong[Dim];
            ilPSP.HilbertCurve.HilbertCurve.hilbert_i2c(nBits, index, coord);
            Console.WriteLine("coord: {0}", String.Join(",", coord));
            for (int i = 0; i < checkcoord.Length; i++) {
                if (coord[i] != checkcoord[i]) {
                    return false;
                }
            }
            return true;
        }

        static private bool test_c2i(int nBits, ulong checkindex, ulong[] coord) {
            int Dim = coord.Length;
            ulong index = ilPSP.HilbertCurve.HilbertCurve.hilbert_c2i(nBits, coord);
            Console.WriteLine("index: {0}", index);
            if (index == checkindex)
                return true;
            return false;
        }

        private static CNSControl ShockTube_PartTest(string dbPath, string SessionID, string GridID, int numOfCellsX, int numOfCellsY) {
            CNSControl c = new CNSControl();

            int dgDegree = 0;
            double sensorLimit = 1e-4;
            bool true1D = false;
            bool saveToDb = false;

            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;

            c.DynamicLoadBalancing_RedistributeAtStartup = true;
            c.GridPartType = GridPartType.clusterHilbert;

            bool AV = false;

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            if (true1D == false) {
                c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
                c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddBoundaryValue("AdiabaticSlipWall");

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = 0.25;
            c.NoOfTimesteps = 1;

            c.ProjectName = "Shock tube";
            if (true1D) {
                c.SessionName = String.Format("Shock tube, 1D, dgDegree = {0}, noOfCellsX = {1}, sensorLimit = {2:0.00E-00}", dgDegree, numOfCellsX, sensorLimit);
            } else {
                c.SessionName = String.Format("Shock tube, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}, GridPartType {7}, NoOfCores {8}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.GridPartType, ilPSP.Environment.MPIEnv.MPI_Size);
            }
            c.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(SessionID), -1);
            c.GridGuid = new Guid(GridID);

            return c;

        }

        private static CNSControl ShockTube_directHilbert(string dbPath, string SessionID, string GridID, int numOfCellsX, int numOfCellsY) {
            CNSControl c = new CNSControl();

            int dgDegree = 0;
            double sensorLimit = 1e-4;
            bool true1D = false;
            bool saveToDb = false;

            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;

            c.DynamicLoadBalancing_RedistributeAtStartup = true;
            c.GridPartType = GridPartType.Hilbert;

            bool AV = false;

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            if (true1D == false) {
                c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
                c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddBoundaryValue("AdiabaticSlipWall");

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = 0.25;
            c.NoOfTimesteps = 1;

            c.ProjectName = "Shock tube";
            if (true1D) {
                c.SessionName = String.Format("Shock tube, 1D, dgDegree = {0}, noOfCellsX = {1}, sensorLimit = {2:0.00E-00}", dgDegree, numOfCellsX, sensorLimit);
            } else {
                c.SessionName = String.Format("Shock tube, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}, GridPartType {7}, NoOfCores {8}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.GridPartType, ilPSP.Environment.MPIEnv.MPI_Size);
            }
            c.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(SessionID), -1);
            c.GridGuid = new Guid(GridID);

            return c;

        }

        private static CNSControl ShockTube_PartTest_Dynamic(int numOfCellsX, int numOfCellsY, int NoOfTimesteps, int NumberOfSubGrids, bool Repart, int RecInt) {
            CNSControl c = new CNSControl();

            int dgDegree = 0;
            double sensorLimit = 1e-4;
            bool true1D = false;
            bool saveToDb = false;

            //string dbPath = @"D:\Weber\BoSSS\test_db";
            string dbPath = null;
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;

            c.GridPartType = GridPartType.clusterHilbert;

            bool AV = false;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            c.ExplicitScheme = ExplicitSchemes.LTS;
            c.ExplicitOrder = 1;

            c.NumberOfSubGrids = NumberOfSubGrids;
            c.ReclusteringInterval = RecInt;
            c.FluxCorrection = false;

            if (Repart) {
                // Add one balance constraint for each subgrid
                c.DynamicLoadBalancing_On = true;
                c.DynamicLoadBalancing_CellCostEstimators.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
                c.DynamicLoadBalancing_ImbalanceThreshold = 0.0;
                c.DynamicLoadBalancing_Period = c.ReclusteringInterval;
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);

                if (true1D) {
                    var grid = Grid1D.LineGrid(xNodes, periodic: false);
                    // Boundary conditions
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });
                    return grid;
                } else {
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                    // Boundary conditions
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });

                    return grid;
                }
            };

            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, delegate (double[] X) {
                double x = X[0];

                if (true1D == false) {
                    double y = X[1];
                }

                if (x <= 0.5) {
                    return 1.0;
                } else {
                    return 0.125;
                }
            });
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, delegate (double[] X) {
                double x = X[0];

                if (true1D == false) {
                    double y = X[1];
                }

                if (x <= 0.5) {
                    return 1.0;
                } else {
                    return 0.1;
                }
            });
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            if (true1D == false) {
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            }


            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            if (true1D == false) {
                c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
                c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }


            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = 0.25;
            //c.dtFixed = 1.5e-3;
            c.NoOfTimesteps = NoOfTimesteps;


            c.ProjectName = String.Format("Shock tube {0} Repartitioning", (Repart ? "with" : "without"));
            if (true1D) {
                c.SessionName = String.Format("{3}, 1D, dgDegree = {0}, noOfCellsX = {1}, sensorLimit = {2:0.00E-00}", dgDegree, numOfCellsX, sensorLimit, c.ProjectName);
            } else {
                c.SessionName = String.Format("{9}, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}, GridPartType {7}, NoOfCores {8}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.GridPartType, ilPSP.Environment.MPIEnv.MPI_Size, c.ProjectName);
            }
            return c;

        }

        /// <summary>
        /// two cost cluster, dividing the domain diagonally. Using 4 processes should yield an equal partitioning.
        /// </summary>
        [Test]
        static public void ClusterHilbertTest(){
            // Arrange -- Grid
            int xMin = -1;
            int xMax = 1;
            int NoOfCores = ilPSP.Environment.MPIEnv.MPI_Size;
            int numOfCells = 4 * NoOfCores;
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCells + 1);
            var grid = Grid2D.Cartesian2DGrid(xNodes, xNodes);

            // Arrange -- CostCluster
            Func<double[], bool> Cond1 = (double[] X) => Math.Abs(X[0] - xMin) + Math.Abs(X[1] - xMin) < 2;
            Func<double[], bool> Cond2 = (double[] X) => !Cond1(X);
            var CostCluster = new List<int[]>();
            var Cluster1 = CreateCostMap(grid, Cond1);
            CostCluster.Add(Cluster1);
            var Cluster2 = CreateCostMap(grid, Cond2);
            CostCluster.Add(Cluster2);

            // Act
            var rankmap = grid.ComputePartitionHilbert(CostCluster,Functype:0);
            //PlotThisShit(costmap,CostCluster,grid);

            // Assert
            int locNoCells = numOfCells * numOfCells / NoOfCores;
            Assert.IsTrue((rankmap.Length == locNoCells).MPIAnd());
        }

        private static int[] CreateCostMap(GridCommons grid, Func<double[], bool> identifier) {
            long TotL = grid.CellPartitioning.TotalLength;
            long LocL = grid.CellPartitioning.LocalLength;
            int[] costmap = new int[LocL];
            costmap.SetAll(1);
            var theDictonary = grid.GetGlobalId2CellIndexMap();
            int D = grid.SpatialDimension;
            var cells = grid.Cells;
            foreach (var Cell in cells) {
                var centercoordinates = new double[D];
                int NoOfNodes = Cell.TransformationParams.NoOfRows;
                //Compute Barycenter of rectangular cells
                for (int d = 0; d < D; d++) {
                    double center = 0;
                    for (int k = 0; k < NoOfNodes; k++) {
                        center += Cell.TransformationParams[k, d];
                    }
                    centercoordinates[d] = center / ((double)NoOfNodes);
                }
                if (identifier(centercoordinates)) {
                    int LocalID = -1;
                    theDictonary.TryGetValue(Cell.GlobalID, out LocalID);
                    costmap[LocalID] = 10;
                }
            }
            return costmap;
        }

        /*
        /// <summary>
        /// Use this for debugging ...
        /// </summary>
        /// <param name="map"></param>
        /// <param name="costlist"></param>
        /// <param name="grid"></param>
        private static void PlotThisShit(int[] map, List<int[]> costlist, GridCommons grid) {
            long L = grid.CellPartitioning.LocalLength;
            var basis = new Basis(grid.GridData, 0);
            var MPIranks = new SinglePhaseField(basis, "MPIrank");
            var CostCluster = new SinglePhaseField(basis, "CostCluster");
            var BorderCells = new SinglePhaseField(basis, "BorderCells");
            var HilbertIdx = new SinglePhaseField(basis, "HilbertIdx");
            int[] FlatCostCluster = new int[L];
            var barray = new BitArray((int)L);
            //var theDictonary = grid.GetGlobalId2CellIndexMap();
           
            var LocHilbertIdx = grid.GetLocHilbertIdcs;

            for (int iCluster = 0; iCluster < costlist.Count(); iCluster++) {
                for (int iCell = 0; iCell < L; iCell++) {
                    if (costlist[iCluster][iCell] == 10)
                        FlatCostCluster[iCell] = iCluster;
                }
            }
            foreach (long iGlobCell in grid.GetZellsOfChangingProc) {
                int iCell = -1;
                if (iGlobCell > grid.CellPartitioning.iE || iGlobCell < grid.CellPartitioning.i0)
                    continue;
                //theDictonary.TryGetValue(iGlobCell, out iCell);
                iCell = grid.CellPartitioning.TransformIndexToLocal(iGlobCell);
                barray[iCell] = true;
            }
            for (int iCell = 0; iCell < L; iCell++) {
                MPIranks.SetMeanValue(iCell, map[iCell]);
                CostCluster.SetMeanValue(iCell, FlatCostCluster[iCell]);
                if (barray[iCell]) {
                    BorderCells.SetMeanValue(iCell, 1.0);
                }
                HilbertIdx.SetMeanValue(iCell, LocHilbertIdx[iCell]);
            }
            
            var list = new List<SinglePhaseField>();
            list.Add(MPIranks);
            list.Add(CostCluster);
            list.Add(BorderCells);
            list.Add(HilbertIdx);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(list, "BLargh.plt", 0.0, 0);
        }
        */
    }
}
