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
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.ShockCapturing;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace CNS.Tests.IBMTests {

    [TestFixture]
    public class IBMALTSTest : TestProgram<IBMControl> {

        public static IBMControl IBMALTSTestContactDiscontinuity(double levelSetPosition, int explicitOrder, int numOfClusters) {
            IBMControl c = new IBMControl();

            //c.DbPath = @"c:\bosss_db";
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

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
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

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            c.AddVariable(CNSVariables.CFL, 0);
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
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double densityLeft = 100.0;
            double densityRight = 1.0;
            double pressure = 1.0;
            double velocityXLeft = 2.0;
            double velocityY = 0.0;

            c.AddBoundaryValue("SubsonicInlet", CompressibleVariables.Density, (X, t) => densityLeft);
            c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityXLeft);
            c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => pressure);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceToLine(X, 0)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityXLeft);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);

            // Time config 
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            c.Endtime = 0.05;
            c.NoOfTimesteps = int.MaxValue;

            // Queries for comparison
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density.Name));
            c.Queries.Add("L2NormVelocityX", QueryLibrary.L2Norm(CNSVariables.Velocity.xComponent.Name));
            c.Queries.Add("L2NormPressure", QueryLibrary.L2Norm(CNSVariables.Pressure));

            return c;
        }

        public static IBMControl IBMShockTubeRotated(string dbPath = null, int savePeriod = 1, int dgDegree = 3, int numOfCellsX = 20, int numOfCellsY = 25, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 3, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double agg = 0.3, string restart = "False", double smoothing = 4.0) {
            IBMControl c = new IBMControl();

            // ### Database ###
            //dbPath = @"c:\bosss_db";                                          // Local

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = false;
            c.forceReclustering = true;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;

            // ### Level-set ###
            c.DomainType = DomainTypes.StaticImmersedBoundary;

            double xMin = 0.6;
            double xMax = 1.0;
            double yMin = 0.2;
            double yMax = 0.7;

            double angle = Math.PI / 6;

            double[] startOfRamp = new double[] { 0.2, 0.0 };
            double[] startOfRamp2 = new double[] { 0.0, 0.2 };

            Func<double, double, double> Ramp = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp[0]) + startOfRamp[1];
            };
            Func<double, double, double> Ramp2 = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp2[0]) + startOfRamp2[1];
            };

            c.LevelSetFunction = (X, t) => -(X[1] - Ramp(X[0], angle)) * (X[1] - Ramp2(X[0], angle));
            c.AddVariable(IBMVariables.LevelSet, 2);

            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;

            // ### Shock-Capturing ###
            bool AV = false;
            if (dgDegree >= 1) {
                AV = true;
            }
            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            double epsilon0 = 1.0;
            double kappa = 0.5;
            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.forceReclustering = false;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            // ### Output variables ###
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // ### Grid ###
            if (restart == "True") {
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("23033126-3fab-4e3e-ad55-be025358ae71"), -1);
                c.GridGuid = new Guid("f0f9dff0-8f9b-4d54-a45c-f22c1516d3e7");
            } else {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                    grid.DefineEdgeTags(X => 1);
                    return grid;
                };
            }

            // ### Boundary conditions ###
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Initial smoothing ###
            double shockAngle = angle + Math.PI / 2;
            double lengthMiddleLine = (xMax - xMin) / Math.Cos(angle);
            double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle) + xMin;
            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            double DistanceToInitialShock(double[] X, double t) {
                // direction vector
                Vector p1 = new Vector(shockPosX, Ramp(shockPosX, angle));
                Vector p2 = new Vector(p1.x - 0.1, p1.y + 0.1 / Math.Tan(angle));
                Vector p = p2 - p1;

                // normal vector
                Vector n = new Vector(p.y, -p.x);
                n.Normalize();

                // Angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = shockAngle;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * X[0] + n.y * X[1];

                // shock speed
                //double vs = 10;

                // distance to line
                //double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);
                double distance = nDotX - (0.5 * lengthMiddleLine + xMin / Math.Cos(angle));

                return distance;
            }

            // Function for smoothing the initial and top boundary conditions
            double SmoothJump(double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = smoothing * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            }

            // ### Initial conditions ###
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;
            double velocityX = 0.0;
            double velocityY = 0.0;
            double discontinuityPosition = 0.5;

            Func<double, double> Jump = x => x <= discontinuityPosition ? 0 : 1;

            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);
            }

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 0.25;
            c.NoOfTimesteps = 4;

            // ### Project and sessions name ###
            c.ProjectName = "IBMALTSTest_Rotated_Shock_Tube";

            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.LTS) {
                tempSessionName = String.Format("IBMST_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_CFLFrac{5}_ALTS{6}_{7}_re{8}_subs{9}_smooth{10}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, smoothing);
            } else if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = String.Format("IBMST_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_CFLFrac{5}_RK{6}_smooth{7}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, smoothing);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = String.Format("IBMST_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_CFLFrac{5}_AB{6}_smooth{7}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, smoothing);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }

            c.SessionName = tempSessionName;

            // Queries for comparison
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density.Name));
            c.Queries.Add("L2NormVelocityX", QueryLibrary.L2Norm(CNSVariables.Velocity.xComponent.Name));
            c.Queries.Add("L2NormPressure", QueryLibrary.L2Norm(CNSVariables.Pressure));

            return c;
        }

        private static Dictionary<string, object> SetupIBMALTSTest(double levelSetPosition, int explicitOrder, int numOfClusters) {
            IBMControl c = IBMALTSTestContactDiscontinuity(levelSetPosition, explicitOrder, numOfClusters);

            c.ProjectName = "IBM ALTS Tests";
            c.SessionName = String.Format("IBM ALTS test (contact discontinuity), levelSetPosition={0}, ALTS {1}/{2}", levelSetPosition, explicitOrder, numOfClusters);

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        private static Dictionary<string, object> SetupIBMShockTubeRotatedTest() {
            IBMControl c = IBMShockTubeRotated();

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        [Test]
        public static void IBMALTSTest1_4_pos1() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.25, explicitOrder: 1, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.790071864462E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.788854382000E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.944271909999E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8886911607585 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199984 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999888 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest1_4_pos2() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.26, explicitOrder: 1, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.789970849048E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.788854382000E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.944271909999E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8886901090685 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199984 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999926 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest1_4_pos3() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.281, explicitOrder: 1, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.7910701220794E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719099990E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.9107012207936 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199984 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999841 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest2_4_pos1() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.25, explicitOrder: 2, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.7854990740138E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719099983E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8549628200074 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199984 + 1e-10),
                Tuple.Create("L2NormPressure", 0.89442719099987 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest2_4_pos2() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.26, explicitOrder: 2, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.7854990029491E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719099985E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8549628340993 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199984 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999876 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest2_4_pos3() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.281, explicitOrder: 2, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.7855060368734E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719099979E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8550603687334 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199984 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999821 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest3_4_pos1() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.25, explicitOrder: 3, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.78549517081183E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819999E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719099992E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.854951708117 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199983 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999939 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest3_4_pos2() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.26, explicitOrder: 3, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.8549517077093E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819998E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719100004E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8549517077073 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.78885438199983 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427191000102 + 1e-10));
        }

        [Test]
        public static void IBMALTSTest3_4_pos3() {
            CheckErrorThresholds(
                SetupIBMALTSTest(levelSetPosition: 0.281, explicitOrder: 3, numOfClusters: 4),
                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2NormDensity", 6.7854937242731E+001 + 1e-10),
                //Tuple.Create("L2NormVelocityX", 1.7888543819999E+000 + 1e-10),
                //Tuple.Create("L2NormPressure", 8.9442719099989E-001 + 1e-10));

                Tuple.Create("L2NormDensity", 67.8549372427292 + 1e-10),
                Tuple.Create("L2NormVelocityX", 1.7888543819997 + 1e-10),
                Tuple.Create("L2NormPressure", 0.894427190999389 + 1e-10));
        }

        [Test]
        public static void IBM_AV_LTS_RotatedShockTubeTest() {
            CheckErrorThresholds(
                SetupIBMShockTubeRotatedTest(),
                Tuple.Create("L2NormDensity", 0.246785706630873 + 1e-13),
                Tuple.Create("L2NormVelocityX", 0.001177926751585 + 1e-13),
                Tuple.Create("L2NormPressure", 0.245646140214043 + 1e-13));
        }
    }
}
