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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace CNS.Tests.IBMTests {

    [TestFixture]
    public class IBMALTSTest : TestProgram<IBMControl> {

        public static IBMControl IBMALTSTestContactDiscontinuity(double levelSetPosition, int explicitOrder, int numOfClusters) {
            IBMControl c = new IBMControl();

            c.DbPath = @"c:\bosss_db";
            c.savetodb = false;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            int numOfCellsX = 40;
            int numOfCellsY = 10;

            int dgDegree = 2;

            c.DomainType = DomainTypes.StaticImmersedBoundary;

            c.LevelSetFunction = delegate (double[] X, double t) {
                double y = X[1];
                return y - levelSetPosition;
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";

            //c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
            //c.SurfaceHMF_ProjectNodesToLevelSet = false;
            //c.SurfaceHMF_RestrictNodes = true;
            //c.SurfaceHMF_UseGaussNodes = false;
            //c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            //c.VolumeHMF_RestrictNodes = true;
            //c.VolumeHMF_UseGaussNodes = false;
            //c.LevelSetQuadratureOrder = 6;

            c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.2;
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // (A)LTS
            c.ExplicitScheme = ExplicitSchemes.LTS;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numOfClusters;
            c.ReclusteringInterval = 1;
            c.FluxCorrection = false;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);

            c.AddVariable(Variables.CFL, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "SubsonicInlet");
                grid.EdgeTagNames.Add(2, "SubsonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {   // bottom
                        return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                        return 3;
                    } else if (Math.Abs(X[0]) < 1e-14) {    // left
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                        return 2;
                    } else {
                        throw new System.Exception("Problem with definition of boundary conditions");
                    }
                });

                return grid;
            };

            Func<double[], double, double> DistanceToLine = delegate (double[] X, double t) {
                // direction vector
                Vector2D p1 = new Vector2D(0.5, 0.0);
                Vector2D p2 = new Vector2D(0.5, 1.0);
                Vector2D p = p2 - p1;

                // normal vector
                Vector2D n = new Vector2D(p.y, -p.x);
                n.Normalize();

                // angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = Math.PI / 2;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * (X[0]) + n.y * (X[1]);

                // shock speed
                double vs = 1;

                // distance to line
                double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);

                return distance;
            };

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double densityLeft = 100.0;
            double densityRight = 1.0;
            double pressure = 1.0;
            double velocityXLeft = 2.0;
            double velocityY = 0.0;

            c.AddBoundaryCondition("SubsonicInlet", Variables.Density, (X, t) => densityLeft);
            c.AddBoundaryCondition("SubsonicInlet", Variables.Velocity.xComponent, (X, t) => velocityXLeft);
            c.AddBoundaryCondition("SubsonicInlet", Variables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryCondition("SubsonicOutlet", Variables.Pressure, (X, t) => pressure);
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToLine(X, 0)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressure);
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => velocityXLeft);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => velocityY);

            // Time config 
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = 0.05;
            c.NoOfTimesteps = int.MaxValue;

            // Queries for comparison
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(Variables.Density.Name));
            c.Queries.Add("L2NormVelocityX", QueryLibrary.L2Norm(Variables.Velocity.xComponent.Name));
            c.Queries.Add("L2NormPressure", QueryLibrary.L2Norm(Variables.Pressure));

            return c;
        }

        private static Dictionary<string, object> SetupIBMALTSTest(double levelSetPosition, int explicitOrder, int numOfClusters) {
            IBMControl c = IBMALTSTestContactDiscontinuity(levelSetPosition, explicitOrder, numOfClusters);

            c.ProjectName = "IBM ALTS Tests";
            c.SessionName = String.Format("IBM ALTS test (contact discontinuity), levelSetPosition={0}, ALTS {1}/{2}", levelSetPosition, explicitOrder, numOfClusters);

            var solver = new Program();
            solver.Init(c, null);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        [Test]
        public static void IBMALTSTest1_4_pos1() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.25, explicitOrder: 1, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.790071864462E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.788854382000E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.944271909999E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest1_4_pos2() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.26, explicitOrder: 1, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.789970849048E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.788854382000E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.944271909999E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest1_4_pos3() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.281, explicitOrder: 1, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7910701220794E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719099990E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest2_4_pos1() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.25, explicitOrder: 2, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7854990740138E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719099983E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest2_4_pos2() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.26, explicitOrder: 2, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7854990029491E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719099985E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest2_4_pos3() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.281, explicitOrder: 2, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7855060368734E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719099979E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest3_4_pos1() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.25, explicitOrder: 3, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7854944137967E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819999E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719099992E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest3_4_pos2() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.26, explicitOrder: 3, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7854944287217E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719100004E-001 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest3_4_pos3() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.281, explicitOrder: 3, numOfClusters: 4),
                Tuple.Create("L2NormDensity", 6.7854937242731E+001 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819999E+000 + 1e-10),
                Tuple.Create("L2NormPressure", 8.9442719099989E-001 + 1e-10));
        }
    }
}
