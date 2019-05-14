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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using CNS.ShockCapturing;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace CNS.Tests.IBMTests {

    [TestFixture]
    public class IBMArtificialViscosityTest : TestProgram<IBMControl> {

        public static IBMControl IBMAVTestContactDiscontinuity() {
            IBMControl c = new IBMControl();

            //c.DbPath = @"c:\bosss_db";
            c.DbPath = null;
            c.savetodb = false;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            int dgDegree = 2;

            double xMin = 0.0;
            double xMax = 1.0;
            double yMin = 0.0;
            double yMax = 1.0;

            int numOfCellsX = 20;
            int numOfCellsY = 5;

            c.DomainType = DomainTypes.StaticImmersedBoundary;

            // Adjust height of cut cells such that we obtain AVCFL_cutcell = 0.5 * AVCFL
            // Here, this only depends on h_min
            double width = (xMax - xMin) / numOfCellsX;
            double height = (yMax - yMin) / numOfCellsY;
            double heightCutCell = (-2.0 * width * height) / (2.0 * height - 2.0 * Math.Sqrt(2) * width - 2.0 * Math.Sqrt(2) * height);

            double levelSetPosition = 2 * height + (height - heightCutCell);

            c.LevelSetFunction = delegate (double[] X, double t) {
                double y = X[1];
                return y - levelSetPosition;
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.2; // Using this agglomeration threshold, no cells are agglomerated (all cells are true cut cells)
            c.AddVariable(IBMVariables.LevelSet, 1);

            bool AV = true;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double sensorLimit = 1.0e-3;
            double epsilon0 = 1.0;
            double kappa = 0.5;

            Variable sensorVariable = Variables.Density;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);

            if (AV) {
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            // Runge-Kutta
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            c.AddVariable(CNSVariables.ShockSensor, 0);

            if (AV) {
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
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
                Vector p1 = new Vector(0.5, 0.0);
                Vector p2 = new Vector(0.5, 1.0);
                Vector p = p2 - p1;

                // normal vector
                Vector n = new Vector(p.y, -p.x);
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
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double densityLeft = 100.0;
            double densityRight = 1.0;
            double pressure = 1.0;
            double velocityXLeft = 2.0;
            double velocityY = 0.0;

            c.AddBoundaryValue("SubsonicInlet", Variables.Density, (X, t) => densityLeft);
            c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityXLeft);
            c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => pressure);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToLine(X, 0)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityXLeft);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);

            // Time config 
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = 5e-3;
            c.NoOfTimesteps = int.MaxValue;

            // Queries for comparison
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(Variables.Density.Name));
            c.Queries.Add("L2NormVelocityX", QueryLibrary.L2Norm(CNSVariables.Velocity.xComponent.Name));
            c.Queries.Add("L2NormPressure", QueryLibrary.L2Norm(CNSVariables.Pressure));

            return c;
        }

        private static Dictionary<string, object> SetupIBMAVTest_NoAgglomeration() {
            IBMControl c = IBMAVTestContactDiscontinuity();

            c.ProjectName = "IBM artificial viscosity tests";
            c.SessionName = String.Format("IBM artificial viscosity test (contact discontinuity), no cut cells are agglomerated");

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        private static Dictionary<string, object> SetupIBMAVTest_Agglomeration() {
            IBMControl c = IBMAVTestContactDiscontinuity();

            c.ProjectName = "IBM artificial viscosity tests";
            c.SessionName = String.Format("IBM artificial viscosity test (contact discontinuity), all cut cells are agglomerated");
            c.AgglomerationThreshold = 0.99;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        // Works only if AV projection is off (really? I don't think so)
        [Test]
        public static void IBMAVTest_NoAgglomeration() {
            CheckErrorThresholds(
                SetupIBMAVTest_NoAgglomeration(),
                //Tuple.Create("L2NormDensity", 53.8771980076417000 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.5491933384829700 + 1e-10),
                //Tuple.Create("L2NormPressure", 0.7745966692416890 + 1e-10));
                Tuple.Create("L2NormDensity", 53.8771971772865 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.54919333848297 + 1e-10),
                Tuple.Create("L2NormPressure", 0.774596669241019 + 1e-10));
        }

        [Test]
        public static void IBMAVTest_Agglomeration() {
            CheckErrorThresholds(
                SetupIBMAVTest_Agglomeration(),
                //Tuple.Create("L2NormDensity", 53.8773088665331 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.54919333848297 + 1e-10),
                //Tuple.Create("L2NormPressure", 0.774596669241494 + 1e-10));
                Tuple.Create("L2NormDensity", 53.8781974413102 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.54919333848297 + 1e-10),
                Tuple.Create("L2NormPressure", 0.774596669241489 + 1e-10));
        }
    }
}
