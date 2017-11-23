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
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using CNS.Residual;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace CNS.Tests.IBMTests {

    /// <summary>
    /// Tests for the inviscid flow around a cylinder that is given by an
    /// immersed boundary. All tests restart a converged solution on a grid a
    /// 64x32 grid on the upper half of the problem domain [-40, 40] x [0, 40].
    /// </summary>
    [TestFixture]
    public class IBMCylinderTest : TestProgram<IBMControl> {

        /// <summary>
        /// Test using zeroth order DG
        /// </summary>
        [Test]
        public static void IBMCylinder0th() {
            int dgDegree = 0;
            IBMControl control = CylinderControl(dgDegree);
            var solver = new Program();
            solver.Init(control);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1.3369E-05),
                Tuple.Create("L2ErrorXMomentum", 2.2514E-05),
                Tuple.Create("L2ErrorYMomentum", 5.4635E-06),
                Tuple.Create("L2ErrorEnergy", 4.9251E-06),
                Tuple.Create("L2ErrorEntropy", 9.7658E-02));
        }

        /// <summary>
        /// Test using first order DG
        /// </summary>
        [Test]
        public static void IBMCylinder1st() {
            int dgDegree = 1;
            IBMControl control = CylinderControl(dgDegree);
            var solver = new Program();
            solver.Init(control);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2.6335E-08),
                Tuple.Create("L2ErrorXMomentum", 9.7240E-08),
                Tuple.Create("L2ErrorYMomentum", 2.8610E-07),
                Tuple.Create("L2ErrorEnergy", 8.3690E-08),
                Tuple.Create("L2ErrorEntropy", 1.8857E-02));
        }

        /// <summary>
        /// Test using second order DG
        /// </summary>
        [Test]
        public static void IBMCylinder2nd() {
            int dgDegree = 2;
            IBMControl control = CylinderControl(dgDegree);
            var solver = new Program();
            solver.Init(control);
            solver.RunSolverMode();
            
            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2.3347E-06),
                Tuple.Create("L2ErrorXMomentum", 9.6413E-06),
                Tuple.Create("L2ErrorYMomentum", 1.0897E-05),
                Tuple.Create("L2ErrorEnergy", 3.03861E-06),
                Tuple.Create("L2ErrorEntropy", 1.1469E-03));
        }

        /// <summary>
        /// Test using third order DG
        /// </summary>
        [Test]
        public static void IBMCylinder3rd() {
            int dgDegree = 3;
            IBMControl control = CylinderControl(dgDegree);
            var solver = new Program();
            solver.Init(control);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1.4935E-06),
                Tuple.Create("L2ErrorXMomentum", 4.8652E-06),
                Tuple.Create("L2ErrorYMomentum", 6.1838E-06),
                Tuple.Create("L2ErrorEnergy", 3.7801E-06),
                Tuple.Create("L2ErrorEntropy", 3.0844E-04));
        }

        private static IBMControl CylinderControl(int dgDegree) {
            string dbPath = @"..\..\Tests\IBMTests\IBMCylinderTests.zip";
            double Mach = 0.2;
            double agglomerationThreshold = 0.3;
            int gridSize = 64;

            var restartData = new Dictionary<int, Tuple<Guid, Guid, Guid>>() {
                { 0, Tuple.Create(new Guid("ae64096b-bab4-4f63-a2cd-99f5920a11e3"), new Guid("486113d4-e700-4bdb-9f92-38a15bac5388"), new Guid("6adec616-275a-444d-86ad-1b2bc8b8cd65")) },
                { 1, Tuple.Create(new Guid("083a99ee-af0b-4948-bd0f-1f972b551b99"), new Guid("4f4e6e60-953a-43c3-b062-c4750ab0ab71"), new Guid("d0615daf-6980-4ecf-af3f-e803877fc29b")) },
                { 2, Tuple.Create(new Guid("d102c10a-0ef2-417e-af1b-4ffad5ea4fe3"), new Guid("f67d08be-3ede-4d74-9a40-829071b44e6b"), new Guid("dfe950aa-9a6b-43d4-b046-e5dd61cc67f3")) },
                { 3, Tuple.Create(new Guid("f45e8802-4e78-45b5-9aac-c155c532b6a3"), new Guid("e670a74f-efb7-41d6-959c-28289e66904e"), new Guid("02e3f737-3f03-4b5a-9387-109dd0f4ec06")) }
            };

            IBMControl c = new IBMControl();
            c.DbPath = dbPath;
            c.savetodb = false;

            int levelSetQuadratureOrder = 2 * dgDegree;

            c.ProjectName = String.Format("IBM cylinder: {0} cells, order {1}", gridSize, dgDegree);
            c.ProjectDescription = String.Format(
                "Flow around cylinder represented by a level set at Mach {0}" +
                    " with cell agglomeration threshold {1} and {2}th order" +
                    " HMF quadrature (classic variant)",
                Mach,
                agglomerationThreshold,
                levelSetQuadratureOrder);

            c.Tags.Add("Cylinder");
            c.Tags.Add("IBM");
            c.Tags.Add("Agglomeration");

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 2);

            var sessionAndGridGuid = restartData[dgDegree];
            c.RestartInfo = new Tuple<Guid, TimestepNumber>(sessionAndGridGuid.Item1, -1);
            c.GridGuid = sessionAndGridGuid.Item2;

            c.GridPartType = GridPartType.ParMETIS;
            c.GridPartOptions = "5";

            double gamma = c.EquationOfState.HeatCapacityRatio;
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, (X, t) => 1.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], (X, t) => Mach * Math.Sqrt(gamma));
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], (X, t) => 0.0);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, (X, t) => 1.0);

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.LevelSetBoundaryTag = "adiabaticSlipWall";

            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));
            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, sessionAndGridGuid.Item3));
            c.Queries.Add("L2ErrorXMomentum", QueryLibrary.L2Error(Variables.Momentum[0], sessionAndGridGuid.Item3));
            c.Queries.Add("L2ErrorYMomentum", QueryLibrary.L2Error(Variables.Momentum[1], sessionAndGridGuid.Item3));
            c.Queries.Add("L2ErrorEnergy", QueryLibrary.L2Error(Variables.Energy, sessionAndGridGuid.Item3));

            c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.SurfaceHMF_ProjectNodesToLevelSet = false;
            c.SurfaceHMF_RestrictNodes = true;
            c.SurfaceHMF_UseGaussNodes = false;
            c.VolumeHMF_NodeCountSafetyFactor = 5.0;
            c.VolumeHMF_RestrictNodes = true;
            c.VolumeHMF_UseGaussNodes = false;

            c.LevelSetQuadratureOrder = levelSetQuadratureOrder;
            c.AgglomerationThreshold = agglomerationThreshold;

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = double.MaxValue;
            c.NoOfTimesteps = 100;

            c.PrintInterval = 1;

            c.ResidualLoggerType = ResidualLoggerTypes.None;

            c.Paramstudy_CaseIdentification = new Tuple<string, object>[] {
                    new Tuple<string, object>("dgDegree", dgDegree),
                };

            return c;
        }
    }
}
