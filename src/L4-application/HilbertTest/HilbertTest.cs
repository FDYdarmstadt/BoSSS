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
using System.Collections.Generic;
using System.Diagnostics;

namespace HilbertTest {

    [TestFixture]
    public class HilbertTest : TestProgram<CNSControl> {

        public static void Main(string[] args) {
#if Debug
            System.Threading.Thread.Sleep(5000);
#endif
            SetUp();
            Test();
            Cleanup();

        }

        [TestFixtureTearDown]
        public static void Cleanup() {
            csMPI.Raw.mpiFinalize();
        }

        [Test]
        public static void Test() {
            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            //Testing coordinate samples
            bool coordresult = TestingCoordinateSamples();
            Assert.IsTrue(coordresult, "Code of HilbertCurve is corrupted");

            //Testing Partition without any Constraints, even distribution of cells among processes
            bool gridevenresult = TestingGridDistributionEven();
            Assert.IsTrue(gridevenresult, "HilbertCurve or mapping (rank->Hilbertcurve) is corrupted");

            //Testing Partition without any Constraints, uneven distribution of cells among processes
            bool gridunevenresult = TestingGridDistributionUneven();
            Assert.IsTrue(gridunevenresult, "Distribution pattern along HilbertCurve is corrupted");

            //Testing Partition yield from directHilbert without any Constraints, even distribution of cells among processes
            bool directHilbert_E = TestingdirectHilbertEven();
            Assert.IsTrue(directHilbert_E, "HilbertCurve or mapping of directHilbert (rank->Hilbertcurve) is corrupted");

            //Testing Partition yield from directHilbert without any Constraints, uneven distribution of cells among processes
            bool directHilbert_UE = TestingdirectHilbertUneven();
            Assert.IsTrue(directHilbert_UE, "Distribution pattern along HilbertCurve is corrupted");

            //Testing Partition with Constraints (LTS), even distribution among processes
            //The testresult is also valid for other Constrainttypes, e.g. AV, but harder to test
            bool gridevendynamic = TestingGridDistributionDynamic();
            Assert.IsTrue(gridevendynamic, "Dynamic Distribution along HilbertCurve is corrupted");

            //Comparing simulations with&without Repartitioning while LTS-Reclustering on
            bool RepNoEffectonResult=TestingReclusteringIndependency();
            Assert.IsTrue(RepNoEffectonResult, "Repartitioning effects result!");
        }

        static private bool TestingGridDistributionEven() {
            //string dbPath = @"D:\Weber\BoSSS\test_db";
            string dbPath = @"..\..\Tests.zip";
            //TestCase: 4x4 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_PartTest(dbPath, "7ac582f5-8913-439b-9f2b-9fbf96141d76", "b7793aee-44b6-44c7-91e7-5debd7f44c3b", 4, 4);

            var solver = new HilbertTest();
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
            return result;
        }

        static private bool TestingdirectHilbertEven() {
            //string dbPath = @"D:\Weber\BoSSS\test_db";
            string dbPath = @"..\..\Tests.zip";
            //TestCase: 4x4 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_directHilbert(dbPath, "7ac582f5-8913-439b-9f2b-9fbf96141d76", "b7793aee-44b6-44c7-91e7-5debd7f44c3b", 4, 4);

            var solver = new HilbertTest();
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
            Console.WriteLine("Test directHilbert: Grid Distribution even");
            Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
            return result;
        }

        static private bool TestingGridDistributionUneven() {
            string dbPath = @"..\..\Tests.zip";
            //TestCase: 3x3 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_PartTest(dbPath, "ccb23f25-04e9-467a-b667-bb3d642b6447", "9b24a2e6-2ce5-4de2-bd08-37a930f0df06", 3, 3);
            var solver = new HilbertTest();
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
                        result &= ((xC > 0) && (xC < 0.33) && (yC > 0) && (yC < 0.33)) ||
                        ((xC > 0) && (xC < 0.67) && (yC > 0.33) && (yC < 0.67));
                        break;
                    case 1:
                        result &= (xC > 0.33) && (xC < 1) && (yC > 0) && (yC < 0.33);
                        break;
                    case 2:
                        result &= (xC > 0.67) && (xC < 1) && (yC > 0.33) && (yC < 1);
                        break;
                    case 3:
                        result &= (xC > 0) && (xC < 0.67) && (yC > 0.67) && (yC < 1);
                        break;
                }
            }
            Console.WriteLine("Test Grid Distribution uneven");
            Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
            return result;
        }

        static private bool TestingdirectHilbertUneven() {
            string dbPath = @"..\..\Tests.zip";
            //TestCase: 3x3 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_directHilbert(dbPath, "ccb23f25-04e9-467a-b667-bb3d642b6447", "9b24a2e6-2ce5-4de2-bd08-37a930f0df06", 3, 3);
            var solver = new HilbertTest();
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
            Console.WriteLine("Test directHilbert: Grid Distribution uneven");
            Console.WriteLine("Process{0}: {1}", solver.MPIRank, result);
            return result;
        }

        static private bool TestingGridDistributionDynamic() {
            //string dbPath = @"D:\Weber\BoSSS\test_db";
            //TestCase: 5x4 grid, Timesteps, LTS-Cluster, PartOn,recInt,AV=false, dgdegree=0, Timestepping=LTS
            CNSControl control = ShockTube_PartTest_Dynamic(5,4,1,2,true,1);
            var solver = new HilbertTest();
            solver.Init(control);
            solver.RunSolverMode();
            bool result = false;

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
                int J0 = solver.GridData.CellPartitioning.i0;
                int JE = solver.GridData.CellPartitioning.iE;

                ulong[] discreteCenter = new ulong[D];
                ulong[] local_HilbertIndex = new ulong[JE - J0];
                int[] local_RankIndex = new int[JE - J0];

                var _Grid = (GridCommons)(solver.Grid);

                for (int j = J0; j < JE; j++) {
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
                result = ItemsAreEqual(RankIndex,checkLTS);
            } else {
                //catching error caused by changes to LTS-Clustering
                Console.WriteLine("Unexpected result for LTS Clusters! Computation of LTS Clusters changed. Test aborted.");
                result = false;
            }
            Console.WriteLine("Test Grid Distribution Dynamic LTS");
            Console.WriteLine("Testresult: {0}", result);
            return result;
        }

        static private bool TestingReclusteringIndependency() {
            
            //TestCase: 5x4 grid, AV=false, dgdegree=0, Timestepping=LTS&RK
            CNSControl ctrRepON = ShockTube_PartTest_Dynamic(5, 4, int.MaxValue, 2, true,5);
            CNSControl ctrRepOFF = ShockTube_PartTest_Dynamic(5, 4, int.MaxValue, 2, false,5);
            var solverRepON = new HilbertTest();
            var solverRepOFF = new HilbertTest();
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
                DGField vriableRepOFF= listOfDGFields_RepOFF[i];
                double L2NormRepOFF = vriableRepOFF.L2Norm();

                Console.WriteLine("{0}-L2Norm Rep ON: {1}", varname[i],L2NormRepON);
                Console.WriteLine("{0}-L2Norm Rep OFF: {1}", varname[i], L2NormRepOFF);
                bool normequal=(Math.Abs(L2NormRepON - L2NormRepOFF) <= 1e-14);
                result &= normequal;
                Console.WriteLine("{0}-L2Norm equal: {1}", varname[i],normequal);
            }
            return result;
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

        static private bool TestingCoordinateSamples() {
            //Validation of H-Curve-Code, Coordsamples taken from "Convergence with Hilbert's Space Filling Curve" by ARTHUR R. BUTZ, p.133
            bool[] testresult = new bool[3];
            bool result = true;
            testresult[0] = test_both(4, 654508, new ulong[] { 4, 12, 6, 12, 12 });
            testresult[1] = test_both(4, 458751, new ulong[] { 8, 8, 0, 15, 0 });
            testresult[2] = test_both(4, 294911, new ulong[] { 7, 0, 15, 15, 0 });
            for (int i = 0; i < testresult.Length; i++) {
                Console.WriteLine("Test_i2c of Coord {0}:{1}", i, testresult[i]);
                result |= result;
            }
            return result;
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
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
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
            c.GridPartType = GridPartType.directHilbert;

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
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
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

        private static CNSControl ShockTube_PartTest_Dynamic(int numOfCellsX, int numOfCellsY, int NoOfTimesteps ,int NumberOfSubGrids, bool Repart, int RecInt) {
            CNSControl c = new CNSControl();

            int dgDegree = 0;
            double sensorLimit = 1e-4;
            bool true1D = false;
            bool saveToDb = false;

            //string dbPath = @"D:\Weber\BoSSS\test_db";
            string dbPath = null;
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;

            c.GridPartType = GridPartType.Hilbert;

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
                c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
                c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
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
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
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

    }
}
