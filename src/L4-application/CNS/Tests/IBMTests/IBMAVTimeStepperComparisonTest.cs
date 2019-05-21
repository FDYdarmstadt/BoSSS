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
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
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

    /// <summary>
    /// This test compares the L2-norm of the density fields in an rotated shock tube case
    /// after five time steps for the follwing timesteppers:
    /// - Runge-Kutta, order 1
    /// - Adams-Bashforth, order 1-3
    /// - Local time stepping, order 1-3 (reclustering forced in every time step, only 1 cluster)
    /// For an explicit order of 1, all the time steppers equal the explicit Euler scheme
    /// and should therefore produce the same results.
    /// For explicit orders greater or equal than 2, AB nud ALTS should behave the same way.
    /// </summary>

    [TestFixture]
    public class IBMAVTimeStepperComparisonTest : TestProgram<IBMControl> {

        public static IBMControl IBMAVTimeStepperTest_ShockTube(int explicitScheme, int explicitOrder) {
            IBMControl c = new IBMControl();

            // ### Database ###
            //string dbPath = @"c:\bosss_db";
            //c.DbPath = dbPath;
            //c.savetodb = dbPath != null;
            //c.saveperiod = 1;
            //c.PrintInterval = 1;
            //c.WriteLTSLog = false;
            //c.WriteLTSConsoleOutput = true;

            // ### Partitioning and load balancing ###
            int dgDegree = 2;
            int numOfCellsX = 75;
            int numOfCellsY = 55;
            double sensorLimit = 1e-3;
            double dtFixed = 0.0;
            double CFLFraction = 0.1;
            int numberOfSubGrids = 1;
            int reclusteringInterval = 1;
            int maxNumOfSubSteps = 10;
            double agg = 0.3;
            double smoothing = 4.0;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;

            // ### Level-set ###
            c.DomainType = DomainTypes.StaticImmersedBoundary;

            double xMin = 0;
            double xMax = 1.5;
            double yMin = 0;
            double yMax = 1.1;

            double[] startOfRamp = new double[] { 0.2, 0.0 };
            double[] startOfRamp2 = new double[] { 0.0, 0.2 };

            Func<double, double, double> Ramp = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp[0]) + startOfRamp[1];
            };
            Func<double, double, double> Ramp2 = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp2[0]) + startOfRamp2[1];
            };

            double angle = Math.PI / 6;
            c.LevelSetFunction = (X, t) => -(X[1] - Ramp(X[0], angle)) * (X[1] - Ramp2(X[0], angle));
            c.AddVariable(IBMVariables.LevelSet, 2);

            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;

            //c.AddVariable(IBMVariables.FluidCells, 1);
            //c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.CutCells, 1);
            //c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.SourceCells, 1);

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
            //double lambdaMax = 2.0;
            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.forceReclustering = true;
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

            //c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            //c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            //c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            //c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);

            //if (AV) {
            //    c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            //}

            //if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
            //    c.AddVariable(CNSVariables.LTSClusters, 0);
            //}
            //c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

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
                double alpha = shockAngle;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * X[0] + n.y * X[1];

                // distance to line
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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 0.25;
            c.NoOfTimesteps = 10;

            // ### Project and sessions name ###
            c.ProjectName = "IBM_AV_LTS_Test";

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

            // Queries for comparison
            c.Queries.Add("L2NormDensity", QueryLibrary.L2Norm(CompressibleVariables.Density.Name));

            return c;
        }

        private static Dictionary<string, object> Setup_IBMAVTimeStepperTest(int explicitScheme, int explicitOrder) {
            IBMControl c = IBMAVTimeStepperTest_ShockTube(explicitScheme, explicitOrder);

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        // Explicit order 1
        [Test]
        public static void IBMAVTimeStepperTest_RK1() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 1, explicitOrder: 1),
                Tuple.Create("L2NormDensity", 0.488419017555324 + 1e-14));
        }

        [Test]
        public static void IBMAVTimeStepperTest_AB1() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 2, explicitOrder: 1),
                Tuple.Create("L2NormDensity", 0.488419017555324 + 1e-14));
        }

        [Test]
        public static void IBMAVTimeStepperTest_ALTS1() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 3, explicitOrder: 1),
                Tuple.Create("L2NormDensity", 0.488419017555324 + 1e-14));
        }

        // Explicit order 2
        [Test]
        public static void IBMAVTimeStepperTest_AB2() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 2, explicitOrder: 2),
                Tuple.Create("L2NormDensity", 0.488419017763754 + 1e-14));
        }

        [Test]
        public static void IBMAVTimeStepperTest_ALTS2() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 3, explicitOrder: 2),
                Tuple.Create("L2NormDensity", 0.488419017763754 + 1e-14));
        }

        // Explicit order 3
        [Test]
        public static void IBMAVTimeStepperTest_AB3() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 2, explicitOrder: 3),
                Tuple.Create("L2NormDensity", 0.488419017764017 + 1e-14));
        }

        [Test]
        public static void IBMAVTimeStepperTest_ALTS3() {
            CheckErrorThresholds(
                Setup_IBMAVTimeStepperTest(explicitScheme: 3, explicitOrder: 3),
                Tuple.Create("L2NormDensity", 0.488419017764017 + 1e-14));
        }
    }
}