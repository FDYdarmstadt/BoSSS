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
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.Diffusion;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.LoadBalancing;
using CNS.MaterialProperty;
using CNS.Residual;
using CNS.ShockCapturing;
using CNS.Source;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CNS {

    public static class TestCases {

        public static CNSControl ShockTube(string dbPath, int savePeriod, int dgDegree, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Shock-Capturing ###
            bool AV = true;
            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            double lambdaMax = 2.0;
            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(Variables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            // ### Output variables ###
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.AddVariable(Variables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            int numOfCellsX = 200;
            int numOfCellsY = 1;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            // ### Boundary conditions ###
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // ### Initial smoothing ###
            //double crossProduct2D(double[] a, double[] b) {
            //    return a[0] * b[1] - a[1] * b[0];
            //}

            // Normal vector of initial shock
            //Vector2D normalVector = new Vector2D(1, 0);

            // Direction vector of initial shock
            //Vector2D r = new Vector2D(normalVector.y, -normalVector.x);
            //r.Normalize();

            // Distance from a point X to the initial shock
            //double[] p = new double[] { 0.5, 0.0 };

            //double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
            //    double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
            //    double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

            //    return distance;
            //}

            //double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            //Func<double, double> SmoothJump = delegate (double distance) {
            //    // smoothing should be in the range of h/p
            //    double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

            //    return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            //};

            // ### Initial conditions ###
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;
            double velocityLeft = 0.0;
            double velocityRight = 0.0;
            double discontinuityPosition = 0.5;

            Func<double, double> Jump = (x => x <= discontinuityPosition ? 0 : 1);

            //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));

            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);

            // ### Evaluation ###
            Material material = new Material(c);
            StateVector stateLeft = StateVector.FromPrimitiveQuantities(
                material, densityLeft, new Vector3D(velocityLeft, 0.0, 0.0), pressureLeft);
            StateVector stateRight = StateVector.FromPrimitiveQuantities(
                material, densityRight, new Vector3D(velocityRight, 0.0, 0.0), pressureRight);

            var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector3D(1.0, 0.0, 0.0));
            double pStar, uStar;
            riemannSolver.GetStarRegionValues(out pStar, out uStar);

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
                Variables.Density,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
                Variables.Velocity.xComponent,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
                Variables.Pressure,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            c.AddVariable(Variables.RiemannDensity, dgDegree);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = 0.25;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "Shock tube";

            if (c.DynamicLoadBalancing_On) {
                c.SessionName = String.Format("Shock tube, p={0}, {1}x{2} cells, s0={3:0.0E-00}, CFLFrac={4}, ALTS {5}/{6}/Re{7}/Sub{8}, Part={9}/Re{10}/Thresh{11}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
            } else {
                c.SessionName = String.Format("Shock tube, p={0}, {1}x{2} cells, s0={3:0.0E-00}, CFLFrac={4}, ALTS {5}/{6}/Re{7}/Sub{8}, Part={9}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString());
            }

            return c;
        }

        public static IBMControl IBMShockTube(string dbPath, int savePeriod, int dgDegree, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps, double agg) {
            IBMControl c = new IBMControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.ParMETIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Level-set ###
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.LevelSetFunction = ((X, t) => X[1] - 0.16);
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = agg;
            c.AddVariable(IBMVariables.LevelSet, 1);
            //c.AddVariable(IBMVariables.FluidCells, 1);
            //c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.CutCells, 1);
            //c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.SourceCells, 1);

            // ### Shock-Capturing ###

            bool AV = true;
            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            double epsilon0 = 1.0;
            double kappa = 0.5;
            double lambdaMax = 2.0;

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(Variables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            //c.ReclusteringInterval = c.DynamicLoadBalancing_Period;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            // ### Output variables ###
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.AddVariable(Variables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            int numOfCellsX = 50;
            int numOfCellsY = 10;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            // ### Boundary conditions ###
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // ### Initial smoothing ###
            //double crossProduct2D(double[] a, double[] b) {
            //    return a[0] * b[1] - a[1] * b[0];
            //}

            // Normal vector of initial shock
            //Vector2D normalVector = new Vector2D(1, 0);

            // Direction vector of initial shock
            //Vector2D r = new Vector2D(normalVector.y, -normalVector.x);
            //r.Normalize();

            // Distance from a point X to the initial shock
            //double[] p = new double[] { 0.5, 0.0 };

            //double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
            //    double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
            //    double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

            //    return distance;
            //}

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            // ### Initial conditions ###
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;

            Func<double, double> Jump = (x => x <= 0.5 ? 0 : 1);

            //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));

            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = 0.25;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "IBM shock tube";

            if (c.DynamicLoadBalancing_On) {
                c.SessionName = String.Format("IBM shock tube, p={0}, {1}x{2} cells, agg={3}, s0={4:0.0E-00}, CFLFrac={5}, ALTS {6}/{7}/Re{8}/Sub{9}, Part={10}/Re{11}/Thresh{12}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
            } else {
                c.SessionName = String.Format("IBM shock tube, p={0}, {1}x{2} cells, agg={3}, s0={4:0.0E-00}, CFLFrac={5}, ALTS {6}/{7}/Re{8}/Sub{9}, Part={10}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString());
            }

            return c;
        }

        public static CNSControl DoubleMachReflection(string dbPath, int savePeriod, int dgDegree, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps, bool restart = false, string sessionID = null, string gridID = null) {
            CNSControl c = new CNSControl();

            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.FluxCorrection = false;

            // Dynamic load balacing
            c.GridPartType = GridPartType.ParMETIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = c.ReclusteringInterval;

            double xMin = 0;
            double xMax = 4;
            double yMin = 0;
            double yMax = 1;
            int numOfCellsX = 400;
            int numOfCellsY = 100;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            bool AV = true;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            sensorLimit = 1e-3;
            double epsilon0 = 1.0;
            double kappa = 1.0;
            double lambdaMax = 20;

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax);
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Viscosity, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            if (dgDegree > 0) {
                c.AddVariable(Variables.Schlieren, dgDegree - 1);
            }
            if (AV) {
                c.AddVariable(Variables.ShockSensor, 0);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            // Grid
            if (!restart) {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                    //var grid = Grid2D.UnstructuredTriangleGrid(xNodes, yNodes);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[1]) < 1e-14) {   // bottom
                            if (X[0] < xWall) {         // bottom left
                                return 1;
                            } else {                    // bottom right
                                return 3;
                            }
                        } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                            return 1;
                        } else if (Math.Abs(X[0]) < 1e-14) {                    // left
                            return 1;
                        } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                            return 2;
                        } else {
                            throw new System.Exception("Boundary condition not specified");
                        }
                    });

                    return grid;
                };

            } else {
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(sessionID), -1);
                c.GridGuid = new Guid(gridID);
            }

            Func<double[], double, double> DistanceToInitialShock = delegate (double[] X, double t) {
                // direction vector
                Vector2D p1 = new Vector2D(xWall, 0.0);
                Vector2D p2 = new Vector2D(xWall + 1 / Math.Tan(Math.PI / 3), 1.0);
                Vector2D p = p2 - p1;

                // normal vector
                Vector2D n = new Vector2D(p.y, -p.x);
                n.Normalize();

                // Angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                double alpha = Math.PI / 3;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n.x * (X[0]) + n.y * (X[1]);

                // shock speed
                double vs = 10;

                // distance to line
                double distance = nDotX - (Math.Sin(alpha) * p1.x + vs * t);

                return distance;
            };

            // Function for smoothing the initial and top boundary conditions
            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            // Function for a sharp jump (no smoothing of initial and top boundary conditions)
            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            // Boundary conditions
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - Jump(X[0] - (0.1 + (X[1] + 20 * t) / 1.732)) * (8.0 - 1.4));
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 7.14471 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (7.14471 - 0.0));
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => -4.125 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (-4.125 - 0.0));
            //c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (116.5 - 1.0));

            c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(DistanceToInitialShock(X, t)) * (8.0 - 1.4));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, t)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, t)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - SmoothJump(DistanceToInitialShock(X, t)) * (116.5 - 1.0));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryCondition("SupersonicOutlet", Variables.Pressure, (X, t) => 1.0);
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            //c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (8.0 - 1.4));
            //c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => 57.157 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (57.157 - 0.0));
            //c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => -33.0 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (-33 - 0.0));
            //c.InitialValues_Evaluators.Add(Variables.Energy, X => 563.544 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (563.544 - 2.5));

            if (!restart) {
                c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(DistanceToInitialShock(X, 0)) * (8.0 - 1.4));
                c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, 0)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, 0)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(Variables.Pressure, X => 116.5 - SmoothJump(DistanceToInitialShock(X, 0)) * (116.5 - 1.0));
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.3;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "Double Mach reflection";

            if (c.DynamicLoadBalancing_On) {
                c.SessionName = String.Format("DMR, p={0}, {1}x{2} cells, s0={3:0.0E-00}, lamdaMax={4}, CFLFrac={5}, ALTS {6}/{7}/Re{8}/Sub{9}, Part={10}/Re{11}/Thresh{12}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, lambdaMax, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
            } else {
                c.SessionName = String.Format("DMR, p={0}, {1}x{2} cells, s0={3:0.0E-00}, lamdaMax={4}, CFLFrac={5}, ALTS {6}/{7}/Re{8}/Sub{9}, Part={10}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, lambdaMax, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString());
            }

            return c;
        }

        public static IBMControl IBMDoubleMachReflection(string dbPath, int savePeriod, int dgDegree, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps, double agg) {
            IBMControl c = new IBMControl();

            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            double xMin = 0.0;
            double xMax = 2.0;
            double yMin = 0.0;
            double yMax = 1.4;

            int numOfCellsX = 200;
            int numOfCellsY = 140;

            // Force cell height to be such that level set only goes through the corner of cells
            //double cellWidth = (xMax - xMin) / numOfCellsX;
            //double cellHeight = Math.Tan(60 * Math.PI / 180) * cellWidth;
            //numOfCellsY = (int)Math.Ceiling((yMax - yMin) / cellHeight);
            //yMax = yMin + numOfCellsY * cellHeight;

            //numOfCellsX = 20;
            //numOfCellsY = 5;

            //double xMin = 0.0;
            //double xMax = 0.2;

            //double yMin = 0.0;
            //double yMax = 0.05;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            double temp = numOfCellsX / xWall;
            bool resolutionOk = (temp == Math.Truncate(temp));
            if (!resolutionOk) {
                throw new Exception("Number of cells in x-direction is not applicable because of xWall!");
            }

            // Level set
            double angleInDegree = 30;
            double beta = 2 * Math.PI / 360 * angleInDegree;   // the wall has an angle of 60 degree
            double[] startOfRamp = new double[] { xWall, 0.0 };

            Func<double, double> ramp = delegate (double x) {
                return Math.Tan(beta) * (x - startOfRamp[0]) + startOfRamp[1];
            };

            // Level-set
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.LevelSetFunction = delegate (double[] X, double t) {
                return X[1] - ramp(X[0]);
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.3;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, 2);

            //c.AddVariable(IBMVariables.FluidCells, 1);
            //c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.CutCells, 1);
            //c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.SourceCells, 1);

            bool AV = true;

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Dynamic load balancing
            c.GridPartType = GridPartType.ParMETIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = c.ReclusteringInterval;
            //c.DynamicLoadBalancing_RedistributeAtStartup = true;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            sensorLimit = 1e-4;
            double epsilon0 = 1.0;
            double kappa = 1.0;
            double lambdaMax = 20;

            if (dgDegree > 0) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            }

            if (AV) {
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax);
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Viscosity, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.Rank, 0);
            if (dgDegree > 0) {
                c.AddVariable(Variables.Schlieren, dgDegree - 1);
            }
            if (AV) {
                c.AddVariable(Variables.ShockSensor, 0);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            // LTS variables
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {   // bottom
                        if (X[0] < xWall) {         // bottom left
                            return 1;
                        } else {                    // bottom right
                            return 3;
                        }
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                        return 1;
                    } else if (Math.Abs(X[0]) < 1e-14) {                    // left
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                        return 2;
                    } else {
                        throw new System.Exception("Boundary condition not specified");
                    }
                });

                return grid;
            };

            // Direction vector of initial shock (vertical)
            Vector2D r = new Vector2D(0.0, 1.0);

            // Current x-position of the shock
            double shockSpeed = 10;
            Func<double, double> getShockXPosition = delegate (double time) {
                return xWall + shockSpeed * time;
            };

            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            // Boundary conditions
            c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(X[0] - getShockXPosition(t)) * (8.0 - 1.4));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => 8.25 - SmoothJump(X[0] - getShockXPosition(t)) * (8.25 - 0.0));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => 116.5 - SmoothJump(X[0] - getShockXPosition(t)) * (116.5 - 1.0));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryCondition("SupersonicOutlet");
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(X[0] - getShockXPosition(0)) * (8.0 - 1.4));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 8.25 - SmoothJump(X[0] - getShockXPosition(0)) * (8.25 - 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 116.5 - SmoothJump(X[0] - getShockXPosition(0)) * (116.5 - 1.0));

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.1;
            c.CFLFraction = 0.1;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "Double Mach Reflection";

            if (c.DynamicLoadBalancing_On) {
                c.SessionName = String.Format("IBM DMR, p={0}, {1}x{2} cells, s0={3:0.0E-00}, lamdaMax={4}, CFLFrac={5}, ALTS {6}/{7}/Re{8}/Sub{9}, Part={10}/Re{11}/Thresh{12}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, lambdaMax, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
            } else {
                c.SessionName = String.Format("IBM DMR, p={0}, {1}x{2} cells, s0={3:0.0E-00}, lamdaMax={4}, CFLFrac={5}, ALTS {6}/{7}/Re{8}/Sub{9}, Part={10}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, lambdaMax, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString());
            }

            return c;
        }

        public static CNSControl ShockVortexInteraction(string dbPath = @"c:\bosss_db", int savePeriod = 100, int dgDegree = 2, double sensorLimit = 1e-3, double CFLFraction = 0.2, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 0, int reclusteringInterval = 0, int maxNumOfSubSteps = 0) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Shock-Capturing ###
            bool AV = true;
            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            double lambdaMax = 2.0;
            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(Variables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            // ### Output variables ###
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Temperature, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
                c.AddVariable(Variables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.AddVariable(Variables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 2;
            double yMin = 0;
            double yMax = 1;

            int numOfCellsX = 200;
            int numOfCellsY = 100;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                //var grid = Grid2D.UnstructuredTriangleGrid(xNodes, yNodes);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {   // bottom
                        return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                        return 3;
                    } else if (Math.Abs(X[0]) < 1e-14) {                    // left
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                        return 2;
                    } else {
                        throw new System.Exception("Boundary condition not specified");
                    }
                });

                return grid;
            };


            // ### Initial condtions ###

            // Parameters
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double a = 0.075;
            double b = 0.175;
            double xc = 0.25;
            double yc = 0.5;
            double p0 = 1;
            double rho0 = 1;
            double RGas = 1;
            double T0 = p0 / (RGas * rho0);
            double Mv = 0.7;
            double Ms = 1.5;
            double c0 = Math.Sqrt(gamma * p0 / rho0);
            double vm = Mv * c0;

            double Radius(double[] X) {
                return Math.Sqrt((X[0] - xc) * (X[0] - xc) + (X[1] - yc) * (X[1] - yc));
            }

            #region Check if inside vortex
            bool IsInsideVortex(double[] X) {
                bool result = false;
                if (Radius(X) <= b) {
                    result = true;
                }
                return result;
            }

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            bool IsNearVortex(double[] X) {
                bool result = false;
                if ((X[0] - xc) * (X[0] - xc) + (X[1] - yc) * (X[1] - yc) <= (b * b) + (4.0 * cellSize / Math.Max(dgDegree, 1))) {
                    result = true;
                }
                return result;
            }
            #endregion

            double vPhi(double[] X, double radius) {
                double result = 0;

                if (radius <= a) {
                    result = vm * radius / a;
                } else if (radius > a && radius <= b) {
                    result = vm * a / (a * a - b * b) * (radius - b * b / radius);
                }

                return result;
            }

            double Temperature(double radius) {
                double result = 0;
                double preFactorGlobal = (gamma - 1) / (RGas * gamma);

                if (radius <= a) {
                    result = preFactorGlobal * (vm * vm) / (a * a) * (radius * radius) / 2;
                } else if (radius > a && radius <= b) {
                    double innerCircle = (vm * vm) / 2;
                    double preFactorOuterVortex = (vm * vm) * (a * a) / ((a * a - b * b) * (a * a - b * b));
                    double partOne = (radius * radius  / 2) - (2 * b * b * Math.Log(radius)) - (b * b * b * b / (2 * radius * radius));
                    double partTwo = (a * a / 2) - (2 * b * b * Math.Log(a)) - (b * b * b * b / (2 * a * a));
                    result = preFactorGlobal * preFactorOuterVortex * (partOne - partTwo);
                }
                return result;
            }

            #region Vortex by Dumbser (2016)
            double DensityVortex(double[] X, double temp) {
                return rho0 * Math.Pow(temp / T0, 1 / (gamma - 1));
            }

            double PressureVortex(double[] X, double temp) {
                return p0 * Math.Pow(temp / T0, gamma / (gamma - 1));
            }

            double VelocityXVortex(double[] X, double radius) {
                double result = 0;
                double theta = Math.Atan2(X[1], X[0]);

                if (IsInsideVortex(X)) {
                    result = 1 - Math.Sin(theta) * vPhi(X, radius);
                }

                return result;
            }

            double VelocityYVortex(double[] X, double radius) {
                double result = 0;
                double theta = Math.Atan2(X[1], X[0]);

                if (IsInsideVortex(X)) {
                    result = Math.Cos(theta) * vPhi(X, radius);
                }

                return result;
            }
            #endregion

            // Shock
            double densityLeft = 1;
            double densityRight = ((gamma + 1) * Ms * Ms) / (2 + (gamma - 1) * Ms * Ms) * densityLeft;
            double pressureLeft = 1;
            double pressureRight = 1 + (2 * gamma) / (gamma + 1) * (Ms * Ms - 1) * pressureLeft;
            double velocityXLeft = Ms * Math.Sqrt(gamma * pressureLeft / densityLeft);
            double velocityXRight = (2 + (gamma - 1) * Ms * Ms) / ((gamma + 1) * Ms * Ms) * velocityXLeft;    // (1)
            //double velocityXRight2 = velocityXLeft * densityLeft / densityRight; // equivalent to (1)
            //double MsPostShock = Math.Sqrt((1 + ((gamma - 1) / 2) * Ms * Ms) / (gamma * Ms * Ms - (gamma - 1) / 2));
            //double velocityXRight3 = MsPostShock * Math.Sqrt(gamma * pressureRight / densityRight);     // equivalent to (1)
            double velocityYLeft = 0;
            double velocityYRight = 0;

            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double shockPosition = 0.5;

            // Current x-position of the shock
            Func<double, double> getShockXPosition = delegate (double time) {
                //return shockPosition + velocityXLeft * time;
                return shockPosition;
            };

            double DensityShock(double[] X, double t) {
                return densityLeft - SmoothJump(X[0] - getShockXPosition(t)) * (densityLeft - densityRight);
            }

            double VelocityXShock(double[] X, double t) {
                return velocityXLeft - SmoothJump(X[0] - getShockXPosition(t)) * (velocityXLeft - velocityXRight);
            }

            double VelocityYShock(double[] X, double t) {
                return velocityYLeft - SmoothJump(X[0] - getShockXPosition(t)) * (velocityYLeft - velocityYRight);
            }

            double PressureShock(double[] X, double t) {
                return pressureLeft - SmoothJump(X[0] - getShockXPosition(t)) * (pressureLeft - pressureRight);
            }

            // Stationary shock wave
            //c.InitialValues_Evaluators.Add(Variables.Density, X => DensityShock(X, 0));
            //c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => VelocityXShock(X, 0));
            //c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => VelocityYShock(X, 0));
            //c.InitialValues_Evaluators.Add(Variables.Pressure, X => PressureShock(X, 0));

            // Stationary shock wave and vortex
            c.InitialValues_Evaluators.Add(Variables.Density, X => DensityShock(X, 0) + DensityVortex(X, Temperature(Radius(X))));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => VelocityXShock(X, 0) + VelocityXVortex(X, Radius(X)));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => VelocityYShock(X, 0) + VelocityYVortex(X, Radius(X)));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => PressureShock(X, 0) + PressureVortex(X, Temperature(Radius(X))));

            // ### Boundary condtions ###
            c.AddBoundaryCondition("SupersonicInlet", Variables.Density, (X, t) => DensityShock(X, t));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.xComponent, (X, t) => VelocityXShock(X, t));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Velocity.yComponent, (X, t) => VelocityYShock(X, t));
            c.AddBoundaryCondition("SupersonicInlet", Variables.Pressure, (X, t) => PressureShock(X, t));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryCondition("SupersonicOutlet");
            c.AddBoundaryCondition("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = 0.7;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "Shock-vortex interaction";

            if (c.DynamicLoadBalancing_On) {
                c.SessionName = String.Format("Shock-vortex interaction, p={0}, {1}x{2} cells, s0={3:0.0E-00}, CFLFrac={4}, ALTS {5}/{6}/Re{7}/Sub{8}, Part={9}/Re{10}/Thresh{11}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
            } else {
                c.SessionName = String.Format("Shock-vortex interaction, p={0}, {1}x{2} cells, s0={3:0.0E-00}, CFLFrac={4}, ALTS {5}/{6}/Re{7}/Sub{8}, Part={9}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString());
            }

            return c;
        }
    }
}
