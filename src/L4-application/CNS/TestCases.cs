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
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.LoadBalancing;
using CNS.Residual;
using CNS.ShockCapturing;
using ilPSP.Utils;
using System;
using System.Diagnostics;

namespace CNS {

    public static class TestCases {

        public static CNSControl ShockTube(string dbPath = null, int savePeriod = 100, int dgDegree = 2, int numOfCellsX = 50, int numOfCellsY = 50, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, int refinementLevel = 0) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = true;
            c.forceReclustering = false;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Shock-Capturing ###
            bool AV = false;
            if (dgDegree > 0) {
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
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);    // fix lambdaMax
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

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            //int numOfCellsX = 16 * (int)Math.Pow(2, refinementLevel);
            //int numOfCellsY = 16 * (int)Math.Pow(2, refinementLevel);
            //int numOfCellsY = 4 * (int)Math.Pow(2, refinementLevel);
            //int numOfCellsY = 1;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                //var grid = Grid2D.Trapezoidal2dGrid(3, 2, 4, GenericBlas.Linspace(0, 1, 2));

                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            // ### Boundary conditions ###
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Initial smoothing ###
            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector normalVector = new Vector(1, 0);

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.Normalize();

            //Distance from a point X to the initial shock
            double[] p = new double[] { 0.5, 0.0 };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }

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
            double velocityLeft = 0.0;
            double velocityRight = 0.0;
            double discontinuityPosition = 0.5;

            Func<double, double> Jump = (x => x <= discontinuityPosition ? 0 : 1);

            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));

            //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);

            // ### Evaluation ###
            Material material = c.GetMaterial();
            StateVector stateLeft = StateVector.FromPrimitiveQuantities(
                material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
            StateVector stateRight = StateVector.FromPrimitiveQuantities(
                material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

            var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
            double pStar, uStar;
            riemannSolver.GetStarRegionValues(out pStar, out uStar);

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
                Variables.Density,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
                CNSVariables.Velocity.xComponent,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
                CNSVariables.Pressure,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            c.AddVariable(CNSVariables.RiemannDensity, dgDegree);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 25;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            // ### Project and sessions name ###
            c.ProjectName = "Shock_tube";

            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.LTS) {
                tempSessionName = String.Format("ST_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_re{7}_subs{8}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps);
            } else if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = String.Format("ST_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_RK{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = String.Format("ST_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_AB{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }

            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }


        public static CNSControl ShockTubePaper(string dbPath = null, int savePeriod = 100, int dgDegree = 2, int numOfCellsX = 20, int numOfCellsY = 20, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, int refinementLevel = 0) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            dbPath = @"e:\bosss_db_shock_tube_bug";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = true;
            c.WriteLTSConsoleOutput = false;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Shock-Capturing ###
            bool AV = false;
            if (dgDegree > 0) {
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
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;
            c.forceReclustering = false;

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

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 1;
            double yMin = 0.0;
            double yMax = 1.0;

            //int numOfCellsX = 16 * (int)Math.Pow(2, refinementLevel);
            //int numOfCellsY = 16 * (int)Math.Pow(2, refinementLevel);
            //int numOfCellsY = 4 * (int)Math.Pow(2, refinementLevel);
            //int numOfCellsY = 1;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                //var grid = Grid2D.Trapezoidal2dGrid(3, 2, 4, GenericBlas.Linspace(0, 1, 2));

                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            // ### Boundary conditions ###
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Initial smoothing ###
            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector normalVector = new Vector(1, 0);

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.Normalize();

            //Distance from a point X to the initial shock
            double[] p = new double[] { 0.5, 0.0 };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }

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
            double velocityLeft = 0.0;
            double velocityRight = 0.0;
            double discontinuityPosition = 0.5;

            Func<double, double> Jump = (x => x <= discontinuityPosition ? 0 : 1);

            //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));

            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);

            // ### Evaluation ###
            //Material material = c.GetMaterial();
            //StateVector stateLeft = StateVector.FromPrimitiveQuantities(
            //    material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
            //StateVector stateRight = StateVector.FromPrimitiveQuantities(
            //    material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

            //var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
            //double pStar, uStar;
            //riemannSolver.GetStarRegionValues(out pStar, out uStar);

            //c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
            //    Variables.Density,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            //c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
            //    CNSVariables.Velocity.xComponent,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            //c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
            //    CNSVariables.Pressure,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            //c.AddVariable(CNSVariables.RiemannDensity, dgDegree);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 0.25;
            //c.Endtime = 0.001;
            //c.NoOfTimesteps = 10;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "ST_Paper_Revision_1604";

            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.LTS) {
                tempSessionName = String.Format("ST_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_re{7}_subs{8}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps);
            } else if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = String.Format("ST_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_RK{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = String.Format("ST_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_AB{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }

            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        public static CNSControl ShockVortexInteraction(string dbPath = null, int savePeriod = 10, int dgDegree = 4, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 3, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double Mv = 0.7, double Ms = 1.5, int numOfCellsX = 200, int numOfCellsY = 100) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            //dbPath = @"\\dc1\userspace\Krueger\BoSSS\BoSSS_DBs\markus_ShockVortexInteraction";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = false;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Shock-Capturing ###
            bool AV = false;
            if (dgDegree > 0) {
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
                Variable sensorVariable = Variables.Density;
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

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.Temperature, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            c.AddVariable(CNSVariables.Schlieren, dgDegree);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 2;
            double yMin = 0;
            double yMax = 1;

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
            //double Mv = 0.7;
            //double Ms = 1.5;
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

            //bool IsNearVortex(double[] X) {
            //    bool result = false;
            //    if ((X[0] - xc) * (X[0] - xc) + (X[1] - yc) * (X[1] - yc) <= (b * b) + (4.0 * cellSize / Math.Max(dgDegree, 1))) {
            //        result = true;
            //    }
            //    return result;
            //}
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
                    double partOne = (radius * radius / 2) - (2 * b * b * Math.Log(radius)) - (b * b * b * b / (2 * radius * radius));
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
                    //result = 1 - Math.Sin(theta) * vPhi(X, radius);
                    result = -Math.Sin(theta) * vPhi(X, radius);
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
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => VelocityXShock(X, 0));
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => VelocityYShock(X, 0));
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => PressureShock(X, 0));

            // Stationary shock wave and vortex
            c.InitialValues_Evaluators.Add(Variables.Density, X => DensityShock(X, 0) + DensityVortex(X, Temperature(Radius(X))));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => VelocityXShock(X, 0) + VelocityXVortex(X, Radius(X)));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => VelocityYShock(X, 0) + VelocityYVortex(X, Radius(X)));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => PressureShock(X, 0) + PressureVortex(X, Temperature(Radius(X))));

            // ### Boundary condtions ###
            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => DensityShock(X, t));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => VelocityXShock(X, t));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => VelocityYShock(X, t));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => PressureShock(X, t));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet");
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = 0.7;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "shock_vortex_interaction";

            if (c.DynamicLoadBalancing_On) {
                c.SessionName = String.Format("SVI_p{0}_{1}x{2}cells_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_Re{7}_Sub{8}_Part={9}_Re{10}_Thresh{11}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
            } else {
                c.SessionName = String.Format("SVI_p{0}_{1}x{2}cells_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_Re{7}_Sub{8}_Part={9}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, c.GridPartType.ToString());
            }

            return c;
        }


        public static CNSControl ShockVortexInteractionHiOCFD5(string dbPath = null, int savePeriod = 1000, int dgDegree = 2, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double Mv = 0.9, double Ms = 1.5, int numOfCellsX = 200, int numOfCellsY = 100) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"e:\bosss_db_paper_revision";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = true;
            c.WriteLTSConsoleOutput = false;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_Period = 5;
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.01;
            //c.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((prog, i) => new StaticCellCostEstimator(new[] { 1, 10 }));

            // ### Shock-Capturing ###
            bool AV = false;
            if (dgDegree > 0) {
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
                Variable sensorVariable = Variables.Density;
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

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.Temperature, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            c.AddVariable(CNSVariables.Schlieren, dgDegree);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 2;
            double yMin = 0;
            double yMax = 1;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                grid.EdgeTagNames.Add(2, "SubsonicOutlet");
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

            // #########################
            // Parameters
            // #########################
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double a = 0.075;
            double b = 0.175;
            double xc = 0.25;
            double yc = 0.5;
            double p0 = 1;
            double rho0 = 1;
            double RGas = 1;
            double T0 = p0 / (RGas * rho0);
            double c0 = Math.Sqrt(gamma * p0 / rho0);
            double vm = Mv * c0;
            double Ta = T0 - (gamma - 1.0) / (RGas * gamma) * Math.Pow(vm * a / (a * a - b * b), 2) * (-0.5 * a * a + 2 * b * b * Math.Log(a) + 0.5 * b * b * b * b / (a * a) - 2 * b * b * Math.Log(b));

            // #########################
            // Shock
            // #########################
            double densityLeft = 1;
            double densityRight = (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms) * densityLeft;
            double pressureLeft = 1;
            double pressureRight = (1 + 2 * gamma / (gamma + 1) * (Ms * Ms - 1)) * pressureLeft;
            double velocityXLeft = Ms * Math.Sqrt(gamma * pressureLeft / densityLeft);
            double velocityXRight = (2 + (gamma - 1) * Ms * Ms) / ((gamma + 1) * Ms * Ms) * velocityXLeft;    // (1)
            //double velocityXRight2 = velocityXLeft * densityLeft / densityRight; // equivalent to (1)
            //double MsPostShock = Math.Sqrt((1 + ((gamma - 1) / 2) * Ms * Ms) / (gamma * Ms * Ms - (gamma - 1) / 2));
            //double velocityXRight3 = MsPostShock * Math.Sqrt(gamma * pressureRight / densityRight);     // equivalent to (1)
            double velocityYLeft = 0;
            double velocityYRight = 0;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            double shockPosition = 0.5;

            double DensityShock(double[] X) {
                return densityLeft - SmoothJump(X[0] - shockPosition) * (densityLeft - densityRight);
            }

            double VelocityXShock(double[] X) {
                return velocityXLeft - SmoothJump(X[0] - shockPosition) * (velocityXLeft - velocityXRight);
            }

            double VelocityYShock(double[] X) {
                return velocityYLeft - SmoothJump(X[0] - shockPosition) * (velocityYLeft - velocityYRight);
            }

            double PressureShock(double[] X) {
                return pressureLeft - SmoothJump(X[0] - shockPosition) * (pressureLeft - pressureRight);
            }

            // #########################
            // Vortex
            // #########################
            double Radius(double[] X) {
                return Math.Sqrt((X[0] - xc) * (X[0] - xc) + (X[1] - yc) * (X[1] - yc));
            }

            double vPhi(double r) {
                double result = 0;

                if (r <= a) {
                    result = vm * r / a;
                } else if (r > a && r <= b) {
                    result = vm * a / (a * a - b * b) * (r - b * b / r);
                } else {
                    throw new ArgumentException("Point is outside the vortex");
                }

                return result;
            }

            double TemperatureVortex(double[] X) {
                double result = 0;
                double r = Radius(X);

                if (r <= a) {
                    result = Ta - (gamma - 1.0) / (RGas * gamma) * vm * vm * 0.5 * (1.0 - r * r / (a * a));
                } else if (r > a && r <= b) {
                    result = T0 - (gamma - 1.0) / (RGas * gamma) * Math.Pow(vm * a / (a * a - b * b), 2) * (-0.5 * r * r + 2 * b * b * Math.Log(r) + 0.5 * b * b * b * b / (r * r) - 2 * b * b * Math.Log(b));
                } else {
                    throw new ArgumentException("Point is outside the vortex");
                }

                return result;
            }

            double DensityVortex(double[] X) {
                double result = rho0 * Math.Pow(TemperatureVortex(X) / T0, 1 / (gamma - 1));
                return result;
            }

            double PressureVortex(double[] X) {
                double result = p0 * Math.Pow(TemperatureVortex(X) / T0, gamma / (gamma - 1));
                return result;
            }

            double VelocityXVortex(double[] X) {
                double r = Radius(X);

                //double theta = Math.Atan2(X[1], X[0]);
                //double result = VelocityXShock(X) - vPhi(r) * Math.Sin(theta);

                double dy = X[1] - yc;
                double sin_theta = dy / r;

                double result = VelocityXShock(X) - vPhi(r) * sin_theta;
                return result;
            }

            double VelocityYVortex(double[] X) {
                double r = Radius(X);

                //double theta = Math.Atan2(X[1], X[0]);
                //double result = VelocityYShock(X) + vPhi(r) * Math.Cos(theta);

                double dx = X[0] - xc;
                double cos_theta = dx / r;

                double result = VelocityYShock(X) + vPhi(r) * cos_theta;
                return result;
            }

            // #########################
            // Initial conditions
            // #########################

            // Stationary shock wave
            //c.InitialValues_Evaluators.Add(Variables.Density, X => DensityShock(X, 0));
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => VelocityXShock(X, 0));
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => VelocityYShock(X, 0));
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => PressureShock(X, 0));

            // Stationary shock wave and vortex
            c.InitialValues_Evaluators.Add(Variables.Density, X => Radius(X) <= b ? DensityVortex(X) : DensityShock(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => Radius(X) <= b ? VelocityXVortex(X) : VelocityXShock(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => Radius(X) <= b ? VelocityYVortex(X) : VelocityYShock(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => Radius(X) <= b ? PressureVortex(X) : PressureShock(X));

            // ### Boundary condtions ###
            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => DensityShock(X));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => VelocityXShock(X));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => VelocityYShock(X));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => PressureShock(X));
            c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => PressureShock(X));
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = 0.7;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "svi";

            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.LTS) {
                tempSessionName = String.Format("SVI_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_re{7}_subs{8}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps);
            } else if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = String.Format("SVI_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_RK{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = String.Format("SVI_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_AB{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }

            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        /// <summary>
        /// Version to be submitted on the TU Darmstadt HHLR Lichtenberg cluster
        /// </summary>
        public static CNSControl ShockVortexInteractionHiOCFD5HLLR(int savePeriod = 1000, int dgDegree = 2, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double Mv = 0.9, double Ms = 1.5, int numOfCellsX = 200, int numOfCellsY = 100) {

            // Lichtenberg
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_revision_svi";
            //string restart = "False";

            CNSControl c = ShockVortexInteractionHiOCFD5(dbPath, savePeriod, dgDegree, sensorLimit, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, Mv, Ms, numOfCellsX, numOfCellsY);

            c.ProjectName = "paper_revision_svi_hllr_savePeriod1000";
            //c.Endtime = 0.3;
            //c.NoOfTimesteps = 10;

            return c;
        }


        public static CNSControl DoubleMachReflection(string dbPath = null, int savePeriod = 1, int dgDegree = 3, double xMax = 4.0, double yMax = 1.0, int numOfCellsX = 400, int numOfCellsY = 100, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 0.2, string restart = "False", int cores = int.MaxValue) {
            CNSControl c = new CNSControl();

            //dbPath = @"/work/scratch/yp19ysog/bosss_db_dmr_video";          // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = false;

            double xMin = 0;
            double yMin = 0;

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Dynamic load balacing
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = c.ReclusteringInterval;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;
            double temp = xWall / ((xMax - xMin) / numOfCellsX);
            bool resolutionOk = (temp == Math.Truncate(temp));
            if (!resolutionOk) {
                throw new Exception("Number of cells in x-direction is not applicable because of xWall!");
            }

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            bool AV;
            if (dgDegree > 0) {
                AV = true;
            } else {
                AV = false;
            }

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
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

            //c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            //c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            //c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Viscosity, dgDegree);
            //c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            //c.AddVariable(CNSVariables.Rank, 0);
            //if (dgDegree > 0) {
            //    c.AddVariable(CNSVariables.Schlieren, dgDegree - 1);
            //}
            if (AV) {
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            //c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);
            //if (AV) {
            //    c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            //}
            //if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
            //    c.AddVariable(CNSVariables.LTSClusters, 0);
            //}

            // Grid
            if (restart == "True") {
                // Restart Lichtenberg
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("a96d7c2c-fe35-4fc3-9f51-ef45185fe188"), -1);
                c.GridGuid = new Guid("c544dd46-a9d8-44c8-b5bb-10516f94f0c9");
            } else {
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
            }

            double DistanceToInitialShock(double[] X, double t) {
                // direction vector
                Vector p1 = new Vector(xWall, 0.0);
                Vector p2 = new Vector(xWall + 1 / Math.Tan(Math.PI / 3), 1.0);
                Vector p = p2 - p1;

                // normal vector
                Vector n = new Vector(p.y, -p.x);
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
            }

            // Function for smoothing the initial and top boundary conditions
            double SmoothJump(double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            }

            // Function for a sharp jump (no smoothing of initial and top boundary conditions)
            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            // Boundary conditions
            //c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => 8.0 - Jump(X[0] - (0.1 + (X[1] + 20 * t) / 1.732)) * (8.0 - 1.4));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => 7.14471 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (7.14471 - 0.0));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => -4.125 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (-4.125 - 0.0));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => 116.5 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (116.5 - 1.0));

            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(DistanceToInitialShock(X, t)) * (8.0 - 1.4));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, t)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, t)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => 116.5 - SmoothJump(DistanceToInitialShock(X, t)) * (116.5 - 1.0));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 1.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            //c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (8.0 - 1.4));
            //c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => 57.157 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (57.157 - 0.0));
            //c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => -33.0 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (-33 - 0.0));
            //c.InitialValues_Evaluators.Add(Variables.Energy, X => 563.544 - Jump(X[0] - (0.1 + (X[1] / 1.732))) * (563.544 - 2.5));

            if (restart == "False") {
                c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(DistanceToInitialShock(X, 0)) * (8.0 - 1.4));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, 0)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, 0)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 116.5 - SmoothJump(DistanceToInitialShock(X, 0)) * (116.5 - 1.0));
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = endTime;
            c.CFLFraction = CFLFraction;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "dmr";

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = string.Format("DMR_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_RK{5}_{6}cores",
                    dgDegree, numOfCellsX, numOfCellsY, sensorLimit, CFLFraction, explicitOrder, cores);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("DMR_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_AB{5}",
                    dgDegree, numOfCellsX, numOfCellsY, sensorLimit, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("DMR_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_re{7}_subs{8}",
                    dgDegree, numOfCellsX, numOfCellsY, sensorLimit, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps);
            }
            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        /// <summary>
        /// Version to be submitted on the TU Darmstadt HHLR Lichtenberg cluster
        /// </summary>
        public static CNSControl DoubleMachReflectionHHLR(int savePeriod, int dgDegree, double xMax, double yMax, int numOfCellsX, int numOfCellsY, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps, double endTime, int timeSteps, int cores) {

            // Lichtenberg
            //string dbPath = @"/home/yp19ysog/bosss_db_paper_ibmdmr2";
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_performance2";
            //string dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_ibmdmr_run3_test";
            //string dbPath = @"C:\bosss_db_paper_ibmdmr_scratch_run3_test";
            string restart = "False";

            CNSControl c = DoubleMachReflection(dbPath, savePeriod, dgDegree, xMax, yMax, numOfCellsX, numOfCellsY, sensorLimit, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, endTime, restart, cores);

            c.ProjectName = "dmr_cube_run5";
            c.NoOfTimesteps = timeSteps;

            return c;
        }

        public static IBMControl IBMShockTube(string dbPath = null, int savePeriod = 1, int dgDegree = 5, int numOfCellsX = 75, int numOfCellsY = 55, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double agg = 0.3, string restart = "False", double smoothing = 4.0) {
            IBMControl c = new IBMControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"e:\bosss_db_shock_tube_bug";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = true;
            c.WriteLTSConsoleOutput = false;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new IBMCellClassifier();
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(IBMCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 7, 7, 1 })); // HPC Cluster, 28 cores
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 10, 10, 1 })); // Lichtenberg, 64 cores
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = int.MaxValue;
            //c.DynamicLoadBalancing_RedistributeAtStartup = true;

            // ### Level-set ###
            c.DomainType = DomainTypes.StaticImmersedBoundary;

            double xMin;
            double xMax;
            double yMin;
            double yMax;

            double angle;

            double[] startOfRamp = new double[] { 0.2, 0.0 };
            double[] startOfRamp2 = new double[] { 0.0, 0.2 };

            Func<double, double, double> Ramp = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp[0]) + startOfRamp[1];
            };
            Func<double, double, double> Ramp2 = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp2[0]) + startOfRamp2[1];
            };

            bool rotatedGrid = true;

            if (rotatedGrid) {
                angle = Math.PI / 6;
                c.LevelSetFunction = (X, t) => -(X[1] - Ramp(X[0], angle)) * (X[1] - Ramp2(X[0], angle));
                c.AddVariable(IBMVariables.LevelSet, 2);

                //c.LevelSetFunction = (X, t) => X[1] - Ramp(X[0], angle);
                //c.AddVariable(IBMVariables.LevelSet, 1);
                //c.ContinuousLevelSet = true;
                xMin = 0;
                xMax = 1.5;
                yMin = 0;
                yMax = 1.1;

                //xMin = 0.5;
                //xMax = 1.0;
                //yMin = 0.2;
                //yMax = 0.7;
                //numOfCellsX = 25;
                //numOfCellsY = 25;
            } else {
                c.LevelSetFunction = (X, t) => X[1] - 0.056;
                c.AddVariable(IBMVariables.LevelSet, 1);

                xMin = 0.0;
                xMax = 1.0;
                yMin = 0.0;
                yMax = 1.0;

                numOfCellsX = 20;
                numOfCellsY = 3;

                angle = 0.0;
            }

            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;
            c.AddVariable(IBMVariables.FluidCells, 1);
            c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            c.AddVariable(IBMVariables.CutCells, 1);
            c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            c.AddVariable(IBMVariables.SourceCells, 1);

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
                Variable sensorVariable = Variables.Density;
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
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

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
            //double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle);
            double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle) + xMin;

            //double temp = shockPosX / ((xMax - xMin) / numOfCellsX);
            //bool resolutionOk = (temp == Math.Truncate(temp));
            //if (!resolutionOk) {
            //    throw new Exception("Number of cells in x-direction is not applicable because of xWall!");
            //}

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
                double vs = 10;

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
                if (rotatedGrid) {
                    c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
                    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));
                } else {
                    c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
                    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
                }
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);
            }

            // ### Evaluation ###
            //Material material = c.GetMaterial();
            //StateVector stateLeft = StateVector.FromPrimitiveQuantities(
            //    material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
            //StateVector stateRight = StateVector.FromPrimitiveQuantities(
            //    material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

            //var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
            //double pStar, uStar;
            //riemannSolver.GetStarRegionValues(out pStar, out uStar);

            //c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
            //    Variables.Density,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            //c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
            //    CNSVariables.Velocity.xComponent,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            //c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
            //    CNSVariables.Pressure,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            //c.AddVariable(CNSVariables.RiemannDensity, dgDegree);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 0.25;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "IBMST_Paper_Revision";

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

            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        public static IBMControl IBMContactDiscontinuity(string dbPath = null, int savePeriod = 10, int dgDegree = 2, int numOfCellsX = 75, int numOfCellsY = 55, double sensorLimit = int.MaxValue, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = int.MaxValue, int maxNumOfSubSteps = 0, double agg = 0.3, string restart = "False", double smoothing = 8.0) {
            IBMControl c = new IBMControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"e:\bosss_db_shock_tube_bug";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = false;


            c.SaveAgglomerationPairs = true;


            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new IBMCellClassifier();
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(IBMCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 7, 7, 1 })); // HPC Cluster, 28 cores
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 10, 10, 1 })); // Lichtenberg, 64 cores
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = int.MaxValue;
            //c.DynamicLoadBalancing_RedistributeAtStartup = true;

            // ### Level-set ###
            c.DomainType = DomainTypes.StaticImmersedBoundary;

            double xMin;
            double xMax;
            double yMin;
            double yMax;

            double angle;

            double[] startOfRamp = new double[] { 0.2, 0.0 };
            double[] startOfRamp2 = new double[] { 0.0, 0.2 };

            Func<double, double, double> Ramp = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp[0]) + startOfRamp[1];
            };
            Func<double, double, double> Ramp2 = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp2[0]) + startOfRamp2[1];
            };

            bool rotatedGrid = true;

            if (rotatedGrid) {
                angle = Math.PI / 6;
                c.LevelSetFunction = (X, t) => -(X[1] - Ramp(X[0], angle)) * (X[1] - Ramp2(X[0], angle));
                c.AddVariable(IBMVariables.LevelSet, 2);

                //c.LevelSetFunction = (X, t) => X[1] - Ramp(X[0], angle);
                //c.AddVariable(IBMVariables.LevelSet, 1);
                //c.ContinuousLevelSet = true;
                xMin = 0;
                xMax = 1.5;
                yMin = 0;
                yMax = 1.1;

                //xMin = 0.5;
                //xMax = 1.0;
                //yMin = 0.2;
                //yMax = 0.7;
                //numOfCellsX = 25;
                //numOfCellsY = 25;
            } else {
                c.LevelSetFunction = (X, t) => X[1] - 0.056;
                c.AddVariable(IBMVariables.LevelSet, 1);

                xMin = 0.0;
                xMax = 1.0;
                yMin = 0.0;
                yMax = 1.0;

                numOfCellsX = 20;
                numOfCellsY = 3;

                angle = 0.0;
            }

            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;
            c.AddVariable(IBMVariables.FluidCells, 1);
            c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            c.AddVariable(IBMVariables.CutCells, 1);
            c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            c.AddVariable(IBMVariables.SourceCells, 1);

            // ### Shock-Capturing ###
            bool AV = false;
            //if (dgDegree >= 1) {
            //    AV = true;
            //}
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
                Variable sensorVariable = Variables.Density;
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
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

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
                    grid.EdgeTagNames.Add(2, "SupersonicInlet");
                    grid.EdgeTagNames.Add(3, "SupersonicOutlet");
                    //grid.DefineEdgeTags(X => 1);

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[1]) < 1e-14) {   // bottom
                            return 2;
                        } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                            return 1;
                        } else if (Math.Abs(X[0]) < 1e-14) {                    // left
                            return 2;
                        } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                            return 3;
                        } else {
                            throw new System.Exception("Boundary condition not specified");
                        }
                    });
                    return grid;
                };
            }

            // ### Initial smoothing ###
            double shockAngle = angle + Math.PI / 2;
            double lengthMiddleLine = (xMax - xMin) / Math.Cos(angle);
            //double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle);
            double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle) + xMin;

            //double temp = shockPosX / ((xMax - xMin) / numOfCellsX);
            //bool resolutionOk = (temp == Math.Truncate(temp));
            //if (!resolutionOk) {
            //    throw new Exception("Number of cells in x-direction is not applicable because of xWall!");
            //}

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
                double vs = 10;

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
            double pressureRight = 1.0;
            double shockVelocity = 1.0;
            double velocityX = Math.Cos(angle) * shockVelocity;
            double velocityY = Math.Sin(angle) * shockVelocity;
            double discontinuityPosition = 0.5;

            // ### Boundary conditions ###
            c.AddBoundaryValue("AdiabaticSlipWall");
            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => densityLeft);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressureLeft);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicOutlet");

            Func<double, double> Jump = x => x <= discontinuityPosition ? 0 : 1;

            if (restart == "False") {
                if (rotatedGrid) {
                    c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
                    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));
                } else {
                    c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
                    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
                }
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);
            }

            // ### Evaluation ###
            //Material material = c.GetMaterial();
            //StateVector stateLeft = StateVector.FromPrimitiveQuantities(
            //    material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
            //StateVector stateRight = StateVector.FromPrimitiveQuantities(
            //    material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

            //var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
            //double pStar, uStar;
            //riemannSolver.GetStarRegionValues(out pStar, out uStar);

            //c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
            //    Variables.Density,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            //c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
            //    CNSVariables.Velocity.xComponent,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            //c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
            //    CNSVariables.Pressure,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            //c.AddVariable(CNSVariables.RiemannDensity, dgDegree);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 0.05;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "IBM_ContactDiscontinuity";

            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.LTS) {
                tempSessionName = String.Format("IBMCD_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_CFLFrac{5}_ALTS{6}_{7}_re{8}_subs{9}_smooth{10}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps, smoothing);
            } else if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = String.Format("IBMCD_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_CFLFrac{5}_RK{6}_smooth{7}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, smoothing);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = String.Format("IBMCD_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_CFLFrac{5}_AB{6}_smooth{7}", dgDegree, numOfCellsX, numOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, smoothing);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }

            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        public static IBMControl IBMShockTubePaper(string dbPath = null, int savePeriod = 100, int dgDegree = 2, int numOfCellsX = 10, int numOfCellsY = 10, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = int.MaxValue, int maxNumOfSubSteps = 0, double agg = 0.6, string restart = "False", double smoothing = 2.0) {
            IBMControl c = new IBMControl();

            // ### Database ###
            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/home/ws35kire/test_db";                               // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"e:\bosss_db_shock_tube_bug";                                          // Local
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = true;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;

            // ### Level-set ###
            c.DomainType = DomainTypes.StaticImmersedBoundary;

            double angle;

            double[] startOfRamp = new double[] { 0.2, 0.0 };
            double[] startOfRamp2 = new double[] { 0.0, 0.2 };

            Func<double, double, double> Ramp = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp[0]) + startOfRamp[1];
            };
            Func<double, double, double> Ramp2 = delegate (double x, double ang) {
                return Math.Tan(ang) * (x - startOfRamp2[0]) + startOfRamp2[1];
            };

            c.LevelSetFunction = (X, t) => X[1] - 0.0;
            c.AddVariable(IBMVariables.LevelSet, 1);

            angle = 0.0;

            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;
            c.AddVariable(IBMVariables.FluidCells, 1);
            c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            c.AddVariable(IBMVariables.CutCells, 1);
            c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            c.AddVariable(IBMVariables.SourceCells, 1);

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
                Variable sensorVariable = Variables.Density;
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
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.AddVariable(CNSVariables.Rank, 0);

            double[] yNodes = new double[numOfCellsY + 1];
            double h = 1.0 / numOfCellsY;
            yNodes[0] = -h;
            yNodes[1] = h;
            for (int i = 2; i < yNodes.Length; i++) {
                yNodes[i] = i * h;
            }

            double xMin = 0.0;
            double xMax = 1.0;
            double yMin = -h;
            double yMax = 1.0;

            // ### Grid ###
            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
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
            //double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle);
            double shockPosX = 0.5 * lengthMiddleLine * Math.Cos(angle) + xMin;

            //double temp = shockPosX / ((xMax - xMin) / numOfCellsX);
            //bool resolutionOk = (temp == Math.Truncate(temp));
            //if (!resolutionOk) {
            //    throw new Exception("Number of cells in x-direction is not applicable because of xWall!");
            //}

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
                double vs = 10;

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

                //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
                //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));

                c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);
            }

            // ### Evaluation ###
            //Material material = c.GetMaterial();
            //StateVector stateLeft = StateVector.FromPrimitiveQuantities(
            //    material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
            //StateVector stateRight = StateVector.FromPrimitiveQuantities(
            //    material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

            //var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
            //double pStar, uStar;
            //riemannSolver.GetStarRegionValues(out pStar, out uStar);

            //c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
            //    Variables.Density,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            //c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
            //    CNSVariables.Velocity.xComponent,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            //c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
            //    CNSVariables.Pressure,
            //    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            //c.AddVariable(CNSVariables.RiemannDensity, dgDegree);

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            if (dtFixed != 0.0) {
                c.dtFixed = dtFixed;
            } else {
                c.CFLFraction = CFLFraction;
            }
            c.Endtime = 0.25;
            //c.NoOfTimesteps = 10;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "IBMST_Paper_Revision_1604";

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

            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        /// <summary>
        /// Version to be submitted on the TU Darmstadt HHLR Lichtenberg cluster
        /// </summary>
        public static IBMControl IBMShockTubeHLLR(int savePeriod = 100, int dgDegree = 3, int numOfCellsX = 75, int numOfCellsY = 55, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 3, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double agg = 0.3) {

            // Lichtenberg
            //string dbPath = @"/home/yp19ysog/bosss_db_paper_ibmdmr2";
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_hllr2";
            //string dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_ibmdmr_run3_test";
            //string dbPath = @"C:\bosss_db_paper_ibmdmr_scratch_run3_test";
            string restart = "False";

            IBMControl c = IBMShockTube(dbPath, savePeriod, dgDegree, numOfCellsX, numOfCellsY, sensorLimit, dtFixed, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, agg, restart);

            c.ProjectName = "PAPER_p_study_hllr";
            //c.NoOfTimesteps = 10;

            return c;
        }

        /// <summary>
        /// Version to be submitted on the FDY HPC cluster
        /// </summary>
        public static IBMControl IBMDoubleMachReflection(string dbPath = null, int savePeriod = 1, int dgDegree = 3, int numOfCellsX = 300, int numOfCellsY = 200, double sensorLimit = 1e-4, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 3, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double agg = 0.3, double fugdeFactor = 0.5, double endTime = 0.2, double kappa = 0.5, string restart = "False") {
            //System.Threading.Thread.Sleep(10000);
            //ilPSP.Environment.StdoutOnlyOnRank0 = true;

            IBMControl c = new IBMControl();

            //c.TracingNamespaces = "BoSSS";

            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_ibmdmr";          // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            //c.saveperiod = Convert.ToInt32(50 * 0.5 / CFLFraction);
            c.PrintInterval = 1;

            c.WriteLTSLog = true;
            c.WriteLTSConsoleOutput = false;

            double xMin = 0.0;
            double xMax = 3.0;
            double yMin = 0.0;
            double yMax = 2.0;

            // Force cell height to be such that level set only goes through the corner of cells
            //double cellWidth = (xMax - xMin) / numOfCellsX;
            //double cellHeight = Math.Tan(60 * Math.PI / 180) * cellWidth;
            //numOfCellsY = (int)Math.Ceiling((yMax - yMin) / cellHeight);
            //yMax = yMin + numOfCellsY * cellHeight;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            double temp = xWall / ((xMax - xMin) / numOfCellsX);
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
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, 2);

            //c.AddVariable(IBMVariables.FluidCells, 1);
            //c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.CutCells, 1);
            //c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.SourceCells, 1);

            bool AV;
            if (dgDegree > 0) {
                AV = true;
            } else {
                AV = false;
            }

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Dynamic load balancing
            c.GridPartType = GridPartType.METIS;
            //c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new IBMCellClassifier();
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(IBMCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 7, 7, 1 })); // HPC Cluster, 28 cores
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 10, 10, 1 })); // Lichtenberg, 64 cores
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = int.MaxValue;
            //c.DynamicLoadBalancing_RedistributeAtStartup = true;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            //double kappa = 1.0;   // Only set for DMR
            //double kappa = 0.5;     // Set for all other runs

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 25);    // fix lambdaMax
                                                                                                                                                                 //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, fudgeFactor: fugdeFactor);    // dynamic lambdaMax
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Viscosity, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            if (dgDegree > 0) {
                c.AddVariable(CNSVariables.Schlieren, dgDegree - 1);
            }
            if (AV) {
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // LTS variables
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            if (restart == "True") {
                // Restart Lichtenberg "paper_ibmdmr"
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("8c6e3af9-8f53-4de6-87e5-ba1949732119"), -1);
                c.GridGuid = new Guid("0f4e4dad-7930-428f-80b1-a3ae28dc251c");
            } else {
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
            }

            // Direction vector of initial shock (vertical)
            Vector r = new Vector(0.0, 1.0);

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
            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => 8.0 - SmoothJump(X[0] - getShockXPosition(t)) * (8.0 - 1.4));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => 8.25 - SmoothJump(X[0] - getShockXPosition(t)) * (8.25 - 0.0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => 116.5 - SmoothJump(X[0] - getShockXPosition(t)) * (116.5 - 1.0));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet");
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(Variables.Density, X => 8.0 - SmoothJump(X[0] - getShockXPosition(0)) * (8.0 - 1.4));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 8.25 - SmoothJump(X[0] - getShockXPosition(0)) * (8.25 - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 116.5 - SmoothJump(X[0] - getShockXPosition(0)) * (116.5 - 1.0));
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = endTime;
            c.CFLFraction = CFLFraction;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "ibmdmr";

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = string.Format("IBMDMR_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_ff{5}_ka{6}_dtFixed{7}_CFLFrac{8}_RK{9}",
                    dgDegree, numOfCellsX, numOfCellsY, agg, sensorLimit, fugdeFactor, kappa, dtFixed, CFLFraction, explicitOrder);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("IBMDMR_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_ff{5}_ka{6}_dtFixed{7}_CFLFrac{8}_AB{9}",
                    dgDegree, numOfCellsX, numOfCellsY, agg, sensorLimit, fugdeFactor, kappa, dtFixed, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("IBMDMR_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_ff{5}_ka{6}_dtFixed{7}_CFLFrac{8}_ALTS{9}_{10}_re{11}_subs{12}",
                    dgDegree, numOfCellsX, numOfCellsY, agg, sensorLimit, fugdeFactor, kappa, dtFixed, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps);
            }
            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        /// <summary>
        /// Version to be submitted on the TU Darmstadt HHLR Lichtenberg cluster
        /// </summary>
        public static IBMControl IBMDoubleMachReflectionHHLR(int savePeriod = 1, int dgDegree = 3, int numOfCellsX = 300, int numOfCellsY = 200, double sensorLimit = 1e-4, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 3, int numberOfSubGrids = 3, int reclusteringInterval = 10, int maxNumOfSubSteps = 10, double agg = 0.3, double fugdeFactor = 0.5, double endTime = 0.2, double kappa = 0.5) {

            // Lichtenberg
            //string dbPath = @"/home/yp19ysog/bosss_db_paper_ibmdmr2";
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_ibmdmr_video";
            //string dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_ibmdmr_run3_test";
            //string dbPath = @"C:\bosss_db_paper_ibmdmr_scratch_run3_test";
            string restart = "False";

            IBMControl c = IBMDoubleMachReflection(dbPath, savePeriod, dgDegree, numOfCellsX, numOfCellsY, sensorLimit, dtFixed, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, agg, fugdeFactor, endTime, kappa, restart);

            c.ProjectName = "ibmdmr_paperReproduce_noSubStepLimit";
            //c.NoOfTimesteps = 10;

            return c;
        }

        public static IBMControl IBMGaussianBump(string dbPath = null, int savePeriod = 100, int noOfCellsY = 16, int dgDegree = 2, int lsDegree = 8, double CFL = 0.3, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 100, int maxNumOfSubSteps = 0, double epsilonX = 0.0, double epsilonY = 0.0) {
            IBMControl c = new IBMControl();

            // Solver Settings
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 1000.0;
            c.CFLFraction = CFL;
            //c.NoOfTimesteps = 200000;
            c.NoOfTimesteps = int.MaxValue;

            c.ResidualInterval = 100;
            c.ResidualLoggerType = ResidualLoggerTypes.ChangeRate | ResidualLoggerTypes.Query;
            //c.ResidualBasedTerminationCriteria.Add("changeRate_L2_abs_rhoE", 1E-8);

            // IBM Settings
            c.LevelSetBoundaryTag = "adiabaticSlipWall";
            c.LevelSetQuadratureOrder = 2 * dgDegree;
            c.AgglomerationThreshold = 0.3;

            // NEXT STEP: SET THIS BOOL TO FALSE AND JUST USE IN POSITIVE SUB_VOLUME;
            // THEN TRY BOUNDING BOX APPROACH?
            // WHY THE HELL DOES THIS CONFIGURATION FAIL!??!?!?!?
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.SurfaceHMF_ProjectNodesToLevelSet = false;
            c.SurfaceHMF_RestrictNodes = true;
            c.SurfaceHMF_UseGaussNodes = false;
            c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            c.VolumeHMF_RestrictNodes = true;
            c.VolumeHMF_UseGaussNodes = false;

            //Guid restart = new Guid(" 60688cbc-707d-4777-98e6-d237796ec14c");
            //c.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restart, -1);

            // Session Settings
            //c.DbPath = @"\\fdyprime\userspace\kraemer-eis\FDY-Cluster\dbe_bump\";
            //c.DbPath = @"/home/kraemer/GaussianBump/dbev2/";
            //c.DbPath = @"c:\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = c.DbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 100;
            c.ProjectName = "BoxHMF=" + c.CutCellQuadratureType + "_Ma=0.5_(" + 2 * noOfCellsY + "x" + noOfCellsY + ")_CFL=" + c.CFLFraction + "_lsQuadOrder=" + c.LevelSetQuadratureOrder + "_p=" + dgDegree + "_agg=" + c.AgglomerationThreshold + "_epsX=" + epsilonX + "_epsY=" + epsilonY;
            c.ProjectDescription = "GaussianBump with Ma=0.5";
            c.Tags.Add("Gaussian Bump");
            c.Tags.Add("IBM Test");

            // Solver Type
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Time-Stepping Settings
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Material Settings
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            // Primary CNSVariables
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(IBMVariables.LevelSet, lsDegree);

            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            // Grid

            //switch (noOfCells) {
            //    case 8:
            //        c.GridGuid = new Guid("7337e273-542f-4b97-b592-895ac3422621");
            //        break;
            //    case 16:
            //        c.GridGuid = new Guid("32e5a779-2aef-4ea2-bdef-b158ae785f01");
            //        break;
            //    case 32:
            //        c.GridGuid = new Guid("e96c9f83-3486-4e45-aa3b-9a436445a059");
            //        break;
            //    case 64:
            //        c.GridGuid = new Guid("a86f1b67-4fa3-48ed-b6df-dcea370eb2c0");
            //        break;
            //    default:
            //        throw new ArgumentException("Wrong Grid Input");
            //}

            c.GridFunc = delegate {
                double xBoundary = 20.0;
                double yBoundary = 20.0;
                double yBottom = 0.0;

                double[] xnodes = GenericBlas.Linspace(-xBoundary, xBoundary, 2 * noOfCellsY + 1);

                //double ySplit = 6.0;
                //int ySplitNoOfCells = (int) (0.5*noOfCells);
                //double[] ynodes1 = GenericBlas.Linspace(yBottom, ySplit, ySplitNoOfCells + 1);
                //double[] ynodes2 = GenericBlas.Linspace(ySplit, yBoundary, noOfCells-ySplitNoOfCells + 1);
                //ynodes1 = ynodes1.GetSubVector(0, ynodes1.Length - 1);
                //double[] ynodes = ArrayTools.Cat(ynodes1, ynodes2);

                double[] ynodes = GenericBlas.Linspace(yBottom, yBoundary, noOfCellsY + 1);

                GridCommons grid = Grid2D.Cartesian2DGrid(xnodes, ynodes);

                grid.EdgeTagNames.Add(1, "supersonicinlet");
                grid.EdgeTagNames.Add(2, "subsonicInlet");
                grid.EdgeTagNames.Add(3, "adiabaticSlipWall");
                grid.EdgeTagNames.Add(4, "subsonicOutlet");

                Func<double[], byte> func = delegate (double[] x) {

                    if (Math.Abs(x[0] + xBoundary) < 1e-5) { // Inflow
                        return 1;
                    } else if (Math.Abs(x[0] - xBoundary) < 1e-5) { // Outflow
                        return 1;
                    } else if (Math.Abs(x[1] - yBoundary) < 1e-5) { // Top
                        return 1;
                    } else { // Bottom
                        return 3;
                    }
                };
                grid.DefineEdgeTags(func);
                grid.Name = "IBM-[" + -xBoundary + "," + xBoundary + "]x[" + yBottom + "," + yBoundary + "]_Cells:(" + 2 * noOfCellsY + "x" + noOfCellsY + ")";
                return grid;
            };

            // Functions
            Func<double[], double, double> rho = (X, t) => 1.0;
            Func<double[], double, double> u0 = (X, t) => 1.0;
            Func<double[], double, double> u1 = (X, t) => 0.0;
            Func<double[], double, double> pressure = (X, t) => 2.8571428571428;

            // Initial Values
            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u0(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => u1(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - epsilonY - 0.01 - 0.3939 * Math.Exp(-0.5 * (X[0] - epsilonX) * (X[0] - epsilonX));

            // Supersonic boundary conditions
            c.AddBoundaryValue("adiabaticSlipWall");
            c.AddBoundaryValue("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity.xComponent, u0);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity.yComponent, u1);
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, pressure);

            // Subsonic boundary conditions
            c.AddBoundaryValue("subsonicInlet", Variables.Density, rho);
            c.AddBoundaryValue("subsonicInlet", CNSVariables.Velocity.xComponent, u0);
            c.AddBoundaryValue("subsonicInlet", CNSVariables.Velocity.yComponent, u1);
            c.AddBoundaryValue("subsonicOutlet", CNSVariables.Pressure, pressure);

            // Queries
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 2.8571428571428));

            return c;
        }

        public static IBMControl IBMForwardFacingStep(string dbPath = null, int savePeriod = 100, int dgDegree = 2, int numOfCellsX = 100, int numOfCellsY = 100, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 10, int maxNumOfSubSteps = 10, double agg = 0.3, double fugdeFactor = 0.5, double endTime = 0.5, double kappa = 0.5, string restart = "False") {
            //System.Threading.Thread.Sleep(10000);
            //ilPSP.Environment.StdoutOnlyOnRank0 = true;

            IBMControl c = new IBMControl();

            //c.TracingNamespaces = "BoSSS";

            //dbPath = @"/work/scratch/ws35kire/work_db";                       // Lichtenberg
            //dbPath = @"/work/scratch/yp19ysog/bosss_db_paper_ibmdmr";          // Lichtenberg
            dbPath = @"c:\bosss_db";                                          // Local
                                                                              //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
                                                                              //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            //c.saveperiod = Convert.ToInt32(50 * 0.5 / CFLFraction);
            c.PrintInterval = 1;

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = false;

            double xMin = 0.0;
            double xMax = 3.0 / 3;
            double yMin = 0.0;
            double yMax = 1.0;

            // Wall (corners)
            double h = (xMax - xMin) / numOfCellsX;
            //double shift = 0.8 * h;
            double shift = 0.0;
            double[] wallCorner0 = { 0.6, 0.0 };
            double[] wallCorner1 = { 0.6, 0.2 + shift };

            // Level-set
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.LevelSetFunction = delegate (double[] X, double t) {
                if (X[0] < wallCorner0[0] || (X[0] >= wallCorner1[0] && X[1] > wallCorner1[1])) {
                    return 1;
                } else {
                    //return (wallCorner0[0] - X[0]) * (wallCorner1[1] - X[1]);
                    return -1;
                }
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 3 * dgDegree;
            c.AgglomerationThreshold = agg;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, 2);

            //c.AddVariable(IBMVariables.FluidCells, 1);
            //c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.CutCells, 1);
            //c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.SourceCells, 1);

            bool AV;
            if (dgDegree > 0) {
                AV = true;
            } else {
                AV = false;
            }

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Dynamic load balancing
            c.GridPartType = GridPartType.ParMETIS;
            //c.DynamicLoadBalancing_On = true;
            //c.DynamicLoadBalancing_CellClassifier = new IBMCellClassifier();
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(IBMCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 7, 7, 1 })); // HPC Cluster, 28 cores
            ////c.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 10, 10, 1 })); // Lichtenberg, 64 cores
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = int.MaxValue;
            //c.DynamicLoadBalancing_RedistributeAtStartup = true;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            //double kappa = 1.0;   // Only set for DMR
            //double kappa = 0.5;     // Set for all other runs

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 5);    // fix lambdaMax
                                                                                                                                                                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, fudgeFactor: fugdeFactor);    // dynamic lambdaMax
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Viscosity, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            if (dgDegree > 0) {
                c.AddVariable(CNSVariables.Schlieren, dgDegree - 1);
            }
            if (AV) {
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // LTS variables
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            if (restart == "True") {
                // Restart Lichtenberg "paper_ibmdmr"
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("204ae73a-35b1-4689-8ee3-7c76353240f0"), -1);
                c.GridGuid = new Guid("7c1cfcbf-d0e3-4f29-a1f0-60ec8664ce17");
            } else {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[1]) < 1e-14) {   // bottom
                                                        //if (X[0] < wallCorner0[0]) {         // bottom left
                                                        //    return 1;
                                                        //} else {                    // bottom right
                                                        //    return 3;
                                                        //}
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
            }

            // Conditions
            double Ms = 3;
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double density = gamma;
            double pressure = 1;
            double velocityX = Math.Sqrt(gamma * pressure / density) * Ms;
            double velocityY = 0;

            // Conditions 2
            //double SmoothJump(double distance) {
            //    // smoothing should be in the range of h/p
            //    double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

            //    return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            //}

            //double densityRight = ((gamma + 1) * Ms * Ms) / (2 + (gamma - 1) * Ms * Ms) * density;
            //double pressureRight = 1 + (2 * gamma) / (gamma + 1) * (Ms * Ms - 1) * pressure;
            //double velocityXRight = (2 + (gamma - 1) * Ms * Ms) / ((gamma + 1) * Ms * Ms) * velocityX;

            //if (restart == "False") {
            //    c.InitialValues_Evaluators.Add(Variables.Density, X => density - SmoothJump(X[0] - wallCorner0[0]) * (density - densityRight));
            //    c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX - SmoothJump(X[0] - wallCorner0[0]) * (velocityX - velocityXRight));
            //    c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            //    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure - SmoothJump(X[0] - wallCorner0[0]) * (pressure - pressureRight));
            //}

            // Boundary conditions
            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet");
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(Variables.Density, X => density);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure);
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = endTime;
            c.CFLFraction = CFLFraction;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "ibmffs";

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = string.Format("IBMFFS_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_ff{5}_ka{6}_dtFixed{7}_CFLFrac{8}_RK{9}",
                    dgDegree, numOfCellsX, numOfCellsY, agg, sensorLimit, fugdeFactor, kappa, dtFixed, CFLFraction, explicitOrder);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("IBMFFS_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_ff{5}_ka{6}_dtFixed{7}_CFLFrac{8}_AB{9}",
                    dgDegree, numOfCellsX, numOfCellsY, agg, sensorLimit, fugdeFactor, kappa, dtFixed, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("IBMFFS_p{0}_xCells{1}_yCells{2}_agg{3}_s0={4:0.0E-00}_ff{5}_ka{6}_dtFixed{7}_CFLFrac{8}_ALTS{9}_{10}_re{11}_subs{12}",
                    dgDegree, numOfCellsX, numOfCellsY, agg, sensorLimit, fugdeFactor, kappa, dtFixed, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps);
            }
            if (c.DynamicLoadBalancing_On) {
                string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }

        public static CNSControl BowShock(string dbPath = null, int savePeriod = 1, int dgDegree = 0, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 10.0, string restart = "False", string gridPath = @"c:\GmshGrids\N3\grid2.msh", double lambdaMax = 15.0) {
            CNSControl c = new CNSControl();

            //System.Threading.Thread.Sleep(10000);
            //Debugger.Launch();

            //dbPath = @"/work/scratch/yp19ysog/bosss_db_dmr_video";          // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 100;
            //c.rollingSaves = 10;

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = false;

            //c.TracingNamespaces = "BoSSS.Foundation";

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Dynamic load balacing
            c.GridPartType = GridPartType.METIS;
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = c.ReclusteringInterval;

            bool AV;
            if (dgDegree > 0) {
                AV = true;
            } else {
                AV = false;
            }

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 1.0;

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                if (lambdaMax == -1.0) { // dynamic lambdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                } else { // fixed lamdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);
                }
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Viscosity, dgDegree);
            c.AddVariable(CNSVariables.Enthalpy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            c.AddVariable(CNSVariables.Rank, 0);
            //if (dgDegree > 0) {
            //    c.AddVariable(CNSVariables.Schlieren, dgDegree - 1);
            //}
            if (AV) {
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);
            //if (AV) {
            //    c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            //}
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            // Grid
            if (restart == "True") {
                // Restart Lichtenberg
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("9b3ee853-aaf2-4777-a50d-7c53c5e23ae6"), -1);
                c.GridGuid = new Guid("b1de3801-a54d-4083-8af2-e400947e626a");
            } else {
                c.GridFunc = delegate {
                    GridCommons grid = GridImporter.Import(gridPath);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    //bool IsOnCylinder(double[] X) {
                    //    bool result = false;

                    //    // Circle 1
                    //    double x0 = 0.0;
                    //    double y0 = 0.5;
                    //    double r0 = 0.5;

                    //    // Circle 2
                    //    double x1 = 0.0;
                    //    double y1 = -0.5;
                    //    double r1 = 0.5;

                    //    if ((X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0) - r0 * r0 < 1e-14) {
                    //        result = true;
                    //    }

                    //    if ((X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1) - r1 * r1 < 1e-14) {
                    //        result = true;
                    //    }

                    //    if (Math.Abs(X[0] - (-0.5)) < 1e-14) {
                    //        result = true;
                    //    }

                    //    return result;
                    //}

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[0] - 0.0) < 1e-14) {
                            return 2;
                        } else if ((X[0] - 0.0) * (X[0] - 0.0) + (X[1] - 0.0) * (X[1] - 0.0) - 1.5 * 1.5 >= 1e-14) {
                            return 1;
                        } else {
                            return 3;
                        }
                    });

                    //var gDat = new GridData(grid);
                    //var em1 = gDat.GetBoundaryEdges();
                    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                    return grid;
                };
            }

            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;
            double Mach = 4.0;
            double velocityX = Mach * Math.Sqrt(c.EquationOfState.HeatCapacityRatio * pressure / density);
            double velocityY = 0.0;

            c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 0.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(Variables.Density, X => density);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => velocityX);
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => velocityY);
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressure);
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = endTime;
            c.CFLFraction = CFLFraction;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "bowShock";

            // Extract grid name from grid path
            string gridName = gridPath.Substring(gridPath.Length - 12);
            string[] charsToRemove = new string[] { "m", "s", "h", "h", @"\", "." };
            foreach (var ch in charsToRemove) {
                gridName = gridName.Replace(ch, string.Empty);
            }

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = string.Format("BowShock_{0}_p{1}_s0={2:0.0E-00}_CFLFrac{3}_RK{4}",
                    gridName, dgDegree, sensorLimit, CFLFraction, explicitOrder);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("BowShock_p{0}_s0={1:0.0E-00}_CFLFrac{2}_AB{3}",
                    dgDegree, sensorLimit, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("BowShock_p{0}_s0={1:0.0E-00}_CFLFrac{2}_ALTS{3}_{4}_re{5}_subs{6}",
                    dgDegree, sensorLimit, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps);
            }
            if (c.DynamicLoadBalancing_On) {
                //string loadBal = String.Format("_Part={0}_Repart{1}_Thresh{2}", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                string loadBal = String.Format("_REPART", c.GridPartType.ToString(), c.DynamicLoadBalancing_Period, c.DynamicLoadBalancing_ImbalanceThreshold);
                c.SessionName = String.Concat(tempSessionName, loadBal);
            } else {
                c.SessionName = tempSessionName;
            }

            return c;
        }
    }
}