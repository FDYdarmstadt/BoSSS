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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using CNS;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.LoadBalancing;
using CNS.MaterialProperty;
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
            //Debugger.Launch();
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

            //Testing Partition without any Constraints, equal distribution of cells among processes
            bool gridevenresult = TestingGridDistributionEven();
            Assert.IsTrue(gridevenresult, "HilbertCurve or mapping: rank->Hilbertcurve is corrupted");

            //Testing Partition without any Constraints, uneven distribution of cells among processes
            bool gridunevenresult= TestingGridDistributionUneven();
            Assert.IsTrue(gridunevenresult, "Distribution pattern along 1DHilbertCurve is corrupted");
        }

        static private bool TestingGridDistributionEven() {
            //string dbPath = @"D:\Weber\BoSSS\test_db";
            string dbPath = @"..\..\Tests.zip";
            //TestCase: 4x4 grid, AV=false, dgdegree=0, Timestepping=RK1
            //CNSControl control = ShockTube_PartTest(dbPath,"4c87b1e0-b16a-4c76-ba09-47eebba66132", "fbb0b791-00b3-4da4-9768-e1fbde3c9183");
            CNSControl control = ShockTube_PartTest(dbPath, "c944043e-6bb6-4adf-86e7-2332ae09b2d0", "8651986c-15a5-46eb-b3ca-28acfb7537d1", 4,4);
            
            var solver = new HilbertTest();
            solver.Init(control);
            solver.RunSolverMode();
            bool result = true;

            int Jloc = solver.GridData.CellPartitioning.LocalLength;
            for (int j = 0; j < Jloc; j++) {
                double xC = solver.GridData.Cells.CellCenter[j, 0];
                double yC = solver.GridData.Cells.CellCenter[j, 1];
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

        static private bool TestingGridDistributionUneven() {
            string dbPath = @"..\..\Tests.zip";
            //TestCase: 3x3 grid, AV=false, dgdegree=0, Timestepping=RK1
            CNSControl control = ShockTube_PartTest(dbPath, "8fa48051-a1aa-4864-97bc-620212ac166f", "71e1c7c4-c3c8-404e-ac75-234fdba422c0", 3, 3);
            var solver = new HilbertTest();
            solver.Init(control);
            solver.RunSolverMode();
            bool result = true;

            int Jloc = solver.GridData.CellPartitioning.LocalLength;
            for (int j = 0; j < Jloc; j++) {
                double xC = solver.GridData.Cells.CellCenter[j, 0];
                double yC = solver.GridData.Cells.CellCenter[j, 1];
                switch (solver.MPIRank) {
                    case 0:
                        result &= ((xC > 0) && (xC < 0.33) && (yC > 0) && (yC < 0.33))||
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
            return test_i2c(nBits, index, coord) || test_c2i(nBits, index, coord);
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
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            if (true1D == false) {
                c.AddVariable(Variables.Momentum.yComponent, dgDegree);
                c.AddVariable(Variables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(Variables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(Variables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.AddBoundaryCondition("AdiabaticSlipWall");

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

    }
}
