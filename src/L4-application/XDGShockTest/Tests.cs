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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Queries;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using XDGShock;
using XDGShock.TimeStepping;
using XDGShock.Variables;

namespace XDGShockTest {

    /// <summary>
    /// NUnit tests for the XDGShock project
    /// </summary>
    [TestFixture]
    public static class Tests {

        #region NUnit stuff
        [TestFixtureSetUp]
        static public void Init() {
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out bool dummy);
        }

        [TestFixtureTearDown]
        static public void Cleanup() {
            //Console.Out.Dispose();
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }
        #endregion

        #region Tests
        //[Test]
        //public static void XDGShockTube_SpeciesA() {
        //    XDGShockTestMain program = null;
        //    BoSSS.Solution.Application<XDGShockControl>._Main(
        //        new string[] { "--control", "cs:XDGShockTest.XDGShockTestsControlExamples.XDGTestsShockTube(@\"A\")" },
        //        false,
        //        delegate () {
        //            program = new XDGShockTestMain();
        //            return program;
        //        });

        //    // Compare results
        //    CompareErrors(program);
        //}

        [Test]
        public static void XDGShockTube([Values("A","B")] string species) {
            var C = XDGShockTest.XDGShockTubeControlExamples.XDGTestsShockTube(species);

            using (XDGShockTestMain program = new XDGShockTestMain()) {
                program.Init(C);
                program.RunSolverMode();

                // Compare results
                CompareErrors(program);
            }
        }

        private static void CompareErrors(XDGShockTestMain program) {
            double densityL2Norm = (double)program.QueryHandler.QueryResults["L2NormDensity"];
            double densityCorrect = 0.711325556421112;

            double sensorL2Norm = (double)program.QueryHandler.QueryResults["L2NormSensor"];
            double sensorCorrect = 0.00661458018596098;

            double artificialViscosityL2Norm = (double)program.QueryHandler.QueryResults["L2NormArtificialViscosity"];
            double artificialVisocsityCorrect = 0.0149140409308147;

            double eps = 1e-14;
            Assert.IsTrue(
                 Math.Abs(densityL2Norm - densityCorrect) <= eps,
                 String.Format("Density L2-norm is {0} (threshold is {1}, difference is {2}", densityL2Norm, densityCorrect, Math.Abs(densityL2Norm - densityCorrect))
                 );

            Assert.IsTrue(
                 Math.Abs(sensorL2Norm - sensorCorrect) <= eps,
                 String.Format("Sensor L2-norm is {0} (threshold is {1}, difference is {2}", sensorL2Norm, sensorCorrect, Math.Abs(sensorL2Norm - sensorCorrect))
                 );

            Assert.IsTrue(
                 Math.Abs(artificialViscosityL2Norm - artificialVisocsityCorrect) <= eps,
                 String.Format("Artificial viscosity L2-norm is {0} (threshold is {1}, difference is {2}", artificialViscosityL2Norm, artificialVisocsityCorrect, Math.Abs(artificialViscosityL2Norm - artificialVisocsityCorrect))
                 );
        }
        #endregion
    }

    public static class XDGShockTubeControlExamples {

        public static XDGShockControl XDGTestsShockTube(string speciesName) {
            XDGShockControl c = new XDGShockControl();

            // ### Database ###
            string dbPath = null;
            //dbPath = @"c:\bosss_db";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            // ### DG degree ###
            int dgDegree = 2;

            // ### Time-Stepping ###
            c.TimeSteppingScheme = TimeSteppingSchemes.XdgRK;
            c.ExplicitOrder = 1;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            // ### Level set ###
            double shiftLevelSet;
            if (speciesName.Equals("A")) {
                shiftLevelSet = -10.0;
            } else if (speciesName.Equals("B")) {
                shiftLevelSet = 10.0;
            } else {
                throw new Exception("This should not happen!");
            }

            // Level set is shifted such that it is outside the domain
            // and only species A or B are considered
            c.LevelSetPos = X => X[0] + shiftLevelSet;
            c.AddVariable(XDGShockVariables.LevelSet, 1);

            // ### Shock-capturing ###
            if (dgDegree > 0) {
                // Sensor
                c.SensorLimit = 1e-3;
                c.SensorVariable = CompressibleVariables.Density;
                c.AddVariable(XDGShockVariables.Sensor, 0);

                // Artificial viscosity
                c.AddVariable(XDGShockVariables.ArtificialViscosity, 2);
            }

            // ### Output variables ###
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            // ### Grid ###
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 1.0;

            int numOfCellsX = 10;
            int numOfCellsY = 10;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.DefineEdgeTags((Vector X) => "AdiabaticSlipWall");
                return grid;
            };

            // ### Boundary conditions ###
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Initial conditions ###
            double densityLeft = 1.0;
            double densityRight = 0.125;
            //double xMom = 0.0;
            //double yMom = 0.0;
            double xVel = 0.0;
            double yVel = 0.0;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;
            double energyLeft = pressureLeft / (IdealGas.Air.HeatCapacityRatio - 1.0);
            double energyRight = pressureRight / (IdealGas.Air.HeatCapacityRatio - 1.0);

            double discPos = 0.5;

            // Initial conditions for conservative variables
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => X[0] <= discPos ? densityLeft : densityRight);
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.xComponent, X => xMom);
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.yComponent, X => yMom);
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Energy, X => X[0] <= discPos ? energyLeft : energyRight);

            // Initial conditions for primitive variables
            c.uAInitPrimitive = new Func<double[], double>[] {
                X => X[0] <= discPos ? densityLeft : densityRight,
                X => xVel,
                X => yVel,
                X => X[0] <= discPos ? pressureLeft : pressureRight
            };
            c.uBInitPrimitive = c.uAInitPrimitive;

            // ### Queries ###
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density));
            c.Queries.Add("L2NormSensor", QueryLibrary.L2Norm(XDGShockVariables.Sensor));
            c.Queries.Add("L2NormArtificialViscosity", QueryLibrary.L2Norm(XDGShockVariables.ArtificialViscosity));

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.dtFixed = 1e-4;

            c.Endtime = 0.25;
            c.NoOfTimesteps = 5;

            // ### Project and sessions name ###
            c.ProjectName = "XDGShockTests_ShockTube";

            string tempSessionName;

            if (c.TimeSteppingScheme == TimeSteppingSchemes.XdgRK) {
                tempSessionName = String.Format("TEST_XDGST_p{0}_xCells{1}_yCells{2}_RK{3}_s0={4}_species{5}", dgDegree, numOfCellsX, numOfCellsY, c.ExplicitOrder, c.SensorLimit, speciesName);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }
            c.SessionName = tempSessionName;

            return c;
        }
    }
}
