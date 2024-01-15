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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using static BoSSS.Solution.CompressibleFlowCommon.CompressibleHelperFunc;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.ShockCapturing;
using CNS.Source;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Runtime.CompilerServices;
using BoSSS.Solution;

namespace CNS {

    /// <summary>
    /// Test cases for the CNS IBM solver with the focus on
    /// shock-capturing in supersonic flows
    /// </summary>
    public static class ControlExamples_Supersonic {

        public static CNSControl ShockTube(string dbPath = null, int savePeriod = 1, int dgDegree = 2, int numOfCellsX = 50, int numOfCellsY = 50, double sensorLimit = 1e-3, double dtFixed = 1e-5, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, int refinementLevel = 0) {
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

            c.WriteLTSLog = false;
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
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
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

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
                grid.DefineEdgeTags((Vector X) => 1);
                return grid;
            };

            // ### Boundary conditions ###
            c.AddBoundaryValue("AdiabaticSlipWall");

            //// ### Initial smoothing ###
            //double crossProduct2D(double[] a, double[] b) {
            //    return a[0] * b[1] - a[1] * b[0];
            //}

            // Normal vector of initial shock
            Vector normalVector = new Vector(1, 0);

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.NormalizeInPlace();

            //Distance from a point X to the initial shock
            double[] p = new double[] { 0.5, 0.0 };

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
            double velocityLeft = 0.0;
            double velocityRight = 0.0;
            double discontinuityPosition = 0.5;

            Func<double, double> Jump = (x => x <= discontinuityPosition ? 0 : 1);

            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
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
                CompressibleVariables.Density,
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
            c.Endtime = 0.25;
            c.NoOfTimesteps = 20;

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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
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
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

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
                grid.DefineEdgeTags((Vector X) => 1);
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
            r.NormalizeInPlace();

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

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
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

        public static CNSControl ShockVortexInteractionHiOCFD5(string dbPath = null, int savePeriod = 1000, int dgDegree = 2, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double Mv = 0.9, double Ms = 1.5, int numOfCellsX = 200, int numOfCellsY = 100, double endTime = 0.7, string restart = "False") {
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
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
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.Temperature, dgDegree);
            //c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);
            c.AddVariable(CNSVariables.Schlieren, dgDegree);

            if (AV) {
                //c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            //c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 2;
            double yMin = 0;
            double yMax = 1;

            if (restart == "True") {
                // Restart Lichtenberg "paper_ibmdmr"
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("03785964-854e-4c7f-8cc3-352629732b55"), -1);
                c.GridGuid = new Guid("2169a3f6-9254-4df9-8e77-0565c3f15422");
            } else {
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
            }

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
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => Radius(X) <= b ? DensityVortex(X) : DensityShock(X));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => Radius(X) <= b ? VelocityXVortex(X) : VelocityXShock(X));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => Radius(X) <= b ? VelocityYVortex(X) : VelocityYShock(X));
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => Radius(X) <= b ? PressureVortex(X) : PressureShock(X));
            }

            // ### Boundary condtions ###
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => DensityShock(X));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => VelocityXShock(X));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => VelocityYShock(X));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => PressureShock(X));
            c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => PressureShock(X));
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = endTime;
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

        public static CNSControl StationaryShockWave_Perturbation(bool AV = true, double shockPosition = 3.14, string dbPath = null, int savePeriod = 100, int dgDegree = 2, double CFLFraction = 0.1, double Ms = 1.5, int numOfCellsX = 100, int numOfCellsY = 10, double endTime = 2.13) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"c:\bosss_db";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;

            // ### Shock-Capturing ###
            if (dgDegree > 0) {
                AV = true;
            }
            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            double sensorLimit = 1e-3;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Time-Stepping ###
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

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

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            if (AV) {
                //c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Temperature, dgDegree);
            //c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            //c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);
            //c.AddVariable(CNSVariables.Schlieren, dgDegree);

            //if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
            //    c.AddVariable(CNSVariables.LTSClusters, 0);
            //}

            //c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 6.28;
            double yMin = 0;
            double yMax = 0.03;

            numOfCellsX = 628;
            numOfCellsY = 3;

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

            // ### Shock conditions ###
            // Parameters
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double densityLeft = 1;
            double densityRight = (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms) * densityLeft;
            double pressureLeft = 1;
            double pressureRight = (1 + 2 * gamma / (gamma + 1) * (Ms * Ms - 1)) * pressureLeft;
            double velocityXLeft = Ms * Math.Sqrt(gamma * pressureLeft / densityLeft);
            double velocityXRight = (2 + (gamma - 1) * Ms * Ms) / ((gamma + 1) * Ms * Ms) * velocityXLeft;    // (1)
            //double velocityXRight2 = velocityXLeft * densityLeft / densityRight; // equivalent to (1)
            //double MsPostShock = Math.Sqrt((1 + ((gamma - 1) / 2) * Ms * Ms) / (gamma * Ms * Ms - (gamma - 1) / 2));
            //double velocityXRight3 = MsPostShock * Math.Sqrt(gamma * pressureRight / densityRight);     // equivalent to (1)

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };


            // See PDF by Yi Zhang (Meeting Geisenhofer/Kummer/Oberlack/Zhang 12/08/2020)
            double u_hat_minus = 1e-5;
            double alpha_minus = 1.0;
            double U0_minus = velocityXLeft;
            double rho0_minus = densityLeft;
            double c_minus = Math.Sqrt(gamma * pressureLeft / densityLeft);    // constant free-stream speed of sound
            double omega = (U0_minus + c_minus) * alpha_minus;
            double wave_length = 2 * Math.PI / omega;

            double u_minus(double[] X, double t) {
                return u_hat_minus * Math.Cos(alpha_minus * X[0] - (U0_minus + c_minus) * alpha_minus * t);
            }

            double p_minus(double[] X, double t) {
                return c_minus * rho0_minus * u_hat_minus * Math.Cos(alpha_minus * X[0] - (U0_minus + c_minus) * alpha_minus * t);
            }

            double rho_minus(double[] X, double t) {
                return p_minus(X, t) / c_minus / c_minus;
            }

            double Density(double[] X, double t) {
                //double tmp = densityLeft - SmoothJump(X[0] - shockPosition) * (densityLeft - densityRight);
                double tmp;
                if (X[0] <= shockPosition) {
                    tmp = densityLeft + rho_minus(X, t);
                } else {
                    tmp = densityRight;
                }
                return tmp;
            }

            double VelocityX(double[] X, double t) {
                //double tmp = velocityXLeft - SmoothJump(X[0] - shockPosition) * (velocityXLeft - velocityXRight);
                double tmp;
                if (X[0] <= shockPosition) {
                    tmp = velocityXLeft + u_minus(X, t);
                } else {
                    tmp = velocityXRight;
                }
                return tmp;
            }

            double Pressure(double[] X, double t) {
                //double tmp = pressureLeft - SmoothJump(X[0] - shockPosition) * (pressureLeft - pressureRight);
                double tmp;
                if (X[0] <= shockPosition) {
                    tmp = pressureLeft + p_minus(X, t);
                } else {
                    tmp = pressureRight;
                }
                return tmp;
            }

            // #########################
            // Initial conditions
            // #########################
            // Stationary shock wave
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => Density(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => VelocityX(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => Pressure(X, 0.0));

            // ### Boundary condtions ###
            //c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => densityLeft + rho_minus(X, t));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityXLeft + u_minus(X, t));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressureLeft + p_minus(X, t));
            //c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => pressureRight);
            //c.AddBoundaryValue("AdiabaticSlipWall");

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => Density(X, t));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => VelocityX(X, t));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => Pressure(X, t));
            c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => Pressure(X, t));
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = CFLFraction;
            c.Endtime = endTime;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "shockWavePerturbation";
            c.SessionName = String.Format("SWP_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_RK{5}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder);

            return c;
        }

        public static CNSControl AcousticWaveProjection(string dbPath = null, int savePeriod = 10000, int dgDegree = 2, double CFLFraction = 0.1, double Ms = 1.5, int numOfCellsX = 100, int numOfCellsY = 10, double endTime = 2.13) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"c:\bosss_db";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;

            // ### Time-Stepping ###
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

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

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Temperature, dgDegree);
            //c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            //c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);
            //c.AddVariable(CNSVariables.Schlieren, dgDegree);

            //if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
            //    c.AddVariable(CNSVariables.LTSClusters, 0);
            //}

            //c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            double xMin = 0;
            double xMax = 6.28;
            double yMin = 0;
            double yMax = 0.1;

            numOfCellsX = 628;
            numOfCellsY = 10;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SubsonicInlet");
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

            // ### Shock conditions###
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double densityLeft = 1;
            double densityRight = (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms) * densityLeft;
            double pressureLeft = 1;
            double pressureRight = (1 + 2 * gamma / (gamma + 1) * (Ms * Ms - 1)) * pressureLeft;
            double velocityXLeft = Ms * Math.Sqrt(gamma * pressureLeft / densityLeft);
            double velocityXRight = (2 + (gamma - 1) * Ms * Ms) / ((gamma + 1) * Ms * Ms) * velocityXLeft;    // (1)
            //double velocityXRight2 = velocityXLeft * densityLeft / densityRight; // equivalent to (1)
            //double MsPostShock = Math.Sqrt((1 + ((gamma - 1) / 2) * Ms * Ms) / (gamma * Ms * Ms - (gamma - 1) / 2));
            //double velocityXRight3 = MsPostShock * Math.Sqrt(gamma * pressureRight / densityRight);     // equivalent to (1)

            // ### Acoustic wave (perturbation) ###
            // See PDF by Yi Zhang (Meeting Geisenhofer/Kummer/Oberlack/Zhang 12/08/2020)
            double u_hat_minus = 1e-5;
            double alpha_minus = 1.0;
            double U0_minus = velocityXLeft;
            double rho0_minus = densityLeft;
            double c_minus = Math.Sqrt(gamma * pressureLeft / densityLeft);    // constant free-stream speed of sound
            double omega = (U0_minus + c_minus) * alpha_minus;
            double wave_length = 2 * Math.PI / omega;

            double u_minus(double[] X, double t) {
                return u_hat_minus * Math.Cos(alpha_minus * X[0] - (U0_minus + c_minus) * alpha_minus * t);
            }

            double p_minus(double[] X, double t) {
                return c_minus * rho0_minus * u_hat_minus * Math.Cos(alpha_minus * X[0] - (U0_minus + c_minus) * alpha_minus * t);
            }

            double rho_minus(double[] X, double t) {
                return p_minus(X, t) / c_minus / c_minus;
            }

            // ### Initial condtions ###
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rho_minus(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => u_minus(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => p_minus(X, 0.0));

            // ### Boundary condtions ###
            //c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => rho_minus(X, t));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => u_minus(X, t));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => p_minus(X, t));
            //c.AddBoundaryValue("AdiabaticSlipWall");

            c.AddBoundaryValue("SubsonicInlet", CompressibleVariables.Density, (X, t) => rho_minus(X, t));
            c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.xComponent, (X, t) => u_minus(X, t));
            c.AddBoundaryValue("SubsonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => p_minus(X, t));
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1e-4;
            c.CFLFraction = CFLFraction;
            c.Endtime = endTime;
            c.NoOfTimesteps = 0;

            // ### Project and sessions name ###
            c.ProjectName = "AcousticWave";
            c.SessionName = String.Format("AW_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}", dgDegree, numOfCellsX, numOfCellsY, c.CFLFraction, c.ExplicitOrder);

            return c;
        }
        /// <summary>
        /// First a (hopefully) sationary stationary shock wave is computed. THen after perStartTime a pertubation is added, either on the left or the rigth boundary
        /// </summary>
        /// <param name="dbPath"></param>
        /// <param name="perStartTime"></param>
        /// <param name="savePeriod"></param>
        /// <param name="s_alpha"></param>
        /// <param name="dgDegree"></param>
        /// <param name="CFLFraction"></param>
        /// <param name="MachL"></param>
        /// <param name="numOfCellsX"></param>
        /// <param name="numOfCellsY"></param>
        /// <param name="endTime"></param>
        /// <param name="shockPosition"></param>
        /// <param name="p_amp_neg"></param>
        /// <param name="p_amp_pos"></param>
        /// <param name="waveform"></param>
        /// <param name="waveLength"></param>
        /// <param name="wavePosition"></param>
        /// <param name="isRestart"></param>
        /// <param name="sessId"></param>
        /// <returns></returns>
        /// <exception cref="System.Exception"></exception>
        /// <exception cref="NotSupportedException"></exception>
        public static CNSControl AcousticWave(string dbPath = null, double perStartTime=12, int savePeriod = 10000, double s_alpha = 20, int dgDegree = 2, double CFLFraction = 0.1,
            double MachL = 1.5, int numOfCellsX = 301, int numOfCellsY = 3, double endTime = 2.13, double shockPosition = 1.5,
        double p_amp_neg = 0.0, double p_amp_pos = 1e-5, string waveform = "bump", double waveLength = 0.4, double wavePosition = -0.4, bool isRestart = false, string sessId=null) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"c:\bosss_db";
            
            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = savePeriod;
            c.ResidualInterval = savePeriod;
            c.ResidualLoggerType = BoSSS.Solution.CompressibleFlowCommon.Residual.ResidualLoggerTypes.Rigorous;
            // ### Partitioning and load balancing ###
            c.GridPartType = GridPartType.METIS;
            
            //Restart
            if (isRestart)
            {
                if (MachL == 1.5 && shockPosition == 1.5 && dgDegree == 2)
                {
                    c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("4861c701-4e55-4c3e-baf1-bb266013efb0"), new TimestepNumber(402500));
                }
                else if (MachL == 1.5 && shockPosition == 0.5 && dgDegree == 2)
                {
                    c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("cd45b796-f15a-43b4-bde2-a00854a27fad"), new TimestepNumber(397500));
                }
                else
                {
                    throw new NotSupportedException("for this parameter configuration no restart available");
                }
                perStartTime = 0.0;
            }
            // ### Time-Stepping ###
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            //if(shockPosition==1.5 && MachL = 1.5)
            //{
            //    BoSSSshell.OpenOrCreateDataBase()
            //    c.RestartInfo = new Tuple<Guid, TimestepNumber> (new Guid("d638d164 - 2857 - 4f44 - a8be - 4f5852b86267"),80);
                
            //}

            // ### Shock-Capturing ###
            bool AV = true;
            //if (dgDegree > 0) {
            //    AV = true;
            //}
            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            double sensorLimit = 1e-3;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Output variables ###
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            //c.AddVariable(CNSVariables.Entropy, dgDegree);
            //c.AddVariable(CNSVariables.Temperature, dgDegree);
            //c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            //c.AddVariable(CNSVariables.CFL, 0);
            //c.AddVariable(CNSVariables.CFLConvective, 0);
            //c.AddVariable(CNSVariables.Schlieren, dgDegree);

            if (AV) {
                //c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            //if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
            //    c.AddVariable(CNSVariables.LTSClusters, 0);
            //}

            //c.AddVariable(CNSVariables.Rank, 0);

            // ### Grid ###
            
            double xMin = 0;
            double xMax = 3.0;
            double yMin = 0;
            double yMax = 0.03;
            //if (wavePosition < shockPosition)
            //{
            //    wavePosition = xMin-waveLength;
            //}
            //else
            //{
            //    wavePosition = xMax+waveLength;
            //}

            //numOfCellsX = 600+1;
            //numOfCellsY = 6;

            c.GridFunc = delegate {
            double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
            double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
            var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

            grid.EdgeTagNames.Add(1, "SupersonicInlet");
            //grid.EdgeTagNames.Add(2, "SupersonicOutlet");
            grid.EdgeTagNames.Add(2, "AdiabaticSlipWall");

            grid.DefineEdgeTags(delegate (double[] X) {
                if (Math.Abs(X[1]) < 1e-14) {   // bottom
                    return 2;
                } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                    return 2;
                } else if (Math.Abs(X[0]) < 1e-14) {                    // left
                    return 1;
                } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                    return 1;
                } else {
                    throw new System.Exception("Boundary condition not specified");
                }
            });

            return grid;
            };

            //supersonic Base Flow - left values
            double densityL = 1; double pressureL = 1;
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double cL = Math.Sqrt(gamma * pressureL / densityL);
            double velocityL = MachL * cL;

            //get rigth values
            (double densityR, double velocityR, double pressureR, double cR, double MachR)
                = ComputeNormalShockWaveRelations(densityL, velocityL, pressureL, MachL, gamma);

            c.shockPosition = shockPosition;
            
            Func<double, double> f_waveform = x => 0; //default no pertubation
            //Waveform
            if (waveform == "1sinus")
            { //one period of sinus
                f_waveform = delegate (double X)
                {
                    if (wavePosition < X && X < wavePosition+waveLength)
                    {
                        return Math.Sin(2 * Math.PI * X / waveLength);
                    }
                    else
                    {
                        return 0;
                    }
                };
            }
            else if (waveform == "halfsinus")
            { //half of period of sinus
                f_waveform = delegate (double X)
                {
                    if (wavePosition < X && X < wavePosition+waveLength)
                    {
                        return Math.Sin(Math.PI * X / waveLength);
                    }
                    else
                    {
                        return 0;
                    }
                };
            }
            else if (waveform == "fullsinus")
            { //full sinus
                f_waveform = delegate (double X)
                {
                    return Math.Sin(2 * Math.PI * X / waveLength);
                };

            }
            else if (waveform == "bump")
            { //c ifnity bumb (differentiable everywhere)

                double L = waveLength/2;
                double bumpPos = wavePosition+L;
                Func<double, double> f_base = x => Math.Exp(-1 / (1 - x * x)) * Math.E;

                f_waveform = delegate (double X)
                {
                    if (bumpPos - L < X && X < bumpPos + L)
                    {
                        return f_base((X - bumpPos) / L);
                    }
                    else
                    {
                        return 0;
                    }
                };
            }

            // Functions for the perturbations
            var PressurePertubation = delegate (double[] X)
            {

                if (X[0] < shockPosition)
                {
                    return p_amp_neg * f_waveform(X[0] - (velocityL - cL) * X[1]) + p_amp_pos * f_waveform(X[0] - (velocityL + cL) * X[1]);
                }
                else
                {
                    return p_amp_neg * f_waveform(X[0] - (velocityR - cR) * X[1]) + p_amp_pos * f_waveform(X[0] - (velocityR + cR) * X[1]); ;
                }

            };
            var DensityPertubation = delegate (double[] X)
            {

                if (X[0] < shockPosition)
                {
                    return PressurePertubation(X) / (cL * cL);
                }
                else
                {
                    return PressurePertubation(X) / (cR * cR);
                }
            };
            var VelocityXPertubation = delegate (double[] X)
            {

                if (X[0] < shockPosition)
                {
                    return -p_amp_neg / (densityL * cL) * f_waveform(X[0] - (velocityL - cL) * X[1]) + p_amp_pos / (densityL * cL) * f_waveform(X[0] - (velocityL + cL) * X[1]);
                }
                else
                {
                    return -p_amp_neg / (densityR * cR) * f_waveform(X[0] - (velocityR - cR) * X[1]) + p_amp_pos / (densityR * cR) * f_waveform(X[0] - (velocityR + cR) * X[1]); ;
                }
            };
            //var PressurePertubation = delegate (double t) {
            //    if (wavePosition < shockPosition)
            //    {
            //        return p_amp_neg * f_waveform(t* (velocityL - cL)) + p_amp_pos * f_waveform(t* (velocityL + cL));
            //    }
            //    else
            //    {
            //        return p_amp_neg * f_waveform(t * (velocityR - cR)) + p_amp_pos * f_waveform(t * (velocityR + cR));
            //    }
            //};
            //var DensityPertubation = delegate (double t) {
            //        if (wavePosition < shockPosition)
            //        {
            //            return PressurePertubation(t) / (cL * cL);
            //        }
            //        else
            //        {
            //            return PressurePertubation(t) / (cR * cR);
            //        }
            //};
            //var VelocityXPertubation = delegate (double t) {
            //        if (wavePosition < shockPosition) 
            //        { 
            //            return p_amp_neg / (densityL * cL) * f_waveform(t * (velocityL - cL)) - p_amp_pos / (densityL * cL) * f_waveform(t * (velocityL + cL));
            //        }
            //        else
            //        {
            //            return p_amp_neg / (densityR * cR) * f_waveform(t * (velocityR - cR)) - p_amp_pos / (densityR * cR) * f_waveform(t * (velocityR + cR)); ;
            //        }
            //};

            double cellSize = (xMax - xMin) / numOfCellsX;

            // Function for smoothing the initial and top boundary conditions
            double SmoothJump(double distance)
            {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            }
            //Old
            //Func<double[], int, double> init_valsL = delegate (double[] x, int component) {
            //    var X = new double[] { x[0], 0 };

            //    if (component == 0)
            //    {
            //        return densityL + c.InitialDensityPertubation(X);
            //    }
            //    else if (component == 1)
            //    {
            //        return velocityL + c.InitialVelocityXPertubation(X);
            //    }
            //    else if (component == 2)
            //    {
            //        return pressureL + c.InitialPressurePertubation(X);
            //    }
            //    else
            //    {
            //        throw new NotSupportedException("component not supported");
            //    }
            //};
            //Func<double[], int, double> init_valsR = delegate (double[] x, int component) {
            //    var X = new double[] { x[0], 0 };
            //    if (component == 0)
            //    {
            //        return densityR + c.InitialDensityPertubation(X);
            //    }
            //    else if (component == 1)
            //    {
            //        return velocityR + c.InitialVelocityXPertubation(X);
            //    }
            //    else if (component == 2)
            //    {
            //        return pressureR + c.InitialPressurePertubation(X);
            //    }
            //    else
            //    {
            //        throw new NotSupportedException("component not supported");
            //    }
            //};
            Func<double[], int, double> init_valsL = delegate (double[] x, int component) {

                if (component == 0)
                {
                    return densityL;
                }
                else if (component == 1)
                {
                    return velocityL;
                }
                else if (component == 2)
                {
                    return pressureL;
                }
                else
                {
                    throw new NotSupportedException("component not supported");
                }
            };
            Func<double[], int, double> init_valsR = delegate (double[] x, int component) {
                if (component == 0)
                {
                    return densityR ;
                }
                else if (component == 1)
                {
                    return velocityR;
                }
                else if (component == 2)
                {
                    return pressureR ;
                }
                else
                {
                    throw new NotSupportedException("component not supported");
                }
            };

            //dirichlet boundary vals
            Func<double[], int, double> init_vals = delegate (double[] x, int component)
            {
                if (x[0] < shockPosition)
                {
                    return init_valsL(x, component);
                }
                else
                {
                    return init_valsR(x, component);
                }
            };
            Func<double[], int, double> smooth_init_vals = delegate (double[] x, int component)
            {
                return init_valsL(x, component) - SmoothJump(x[0] - shockPosition) * (init_valsL(x, component) - init_valsR(x, component));
            };
            Func<double[],double, int, double> bnd_valsL = delegate (double[] x,double t, int component)
            {
                if (t < perStartTime || wavePosition > shockPosition)
                {
                    if (component == 0)
                    {
                        return densityL;
                    }
                    else if (component == 1)
                    {
                        return velocityL;
                    }
                    else if (component == 2)
                    {
                        return pressureL;
                    }
                    else
                    {
                        throw new NotSupportedException("component not supported");
                    }
                }
                else
                {
                    if (component == 0)
                    {
                        return densityL + DensityPertubation(new double[] { x[0], t - perStartTime });
                    }
                    else if (component == 1)
                    {
                        return velocityL + VelocityXPertubation(new double[] { x[0], t - perStartTime });
                    }
                    else if (component == 2)
                    {
                        return pressureL + PressurePertubation(new double[] { x[0], t - perStartTime });
                    }
                    else
                    {
                        throw new NotSupportedException("component not supported");
                    }
                }
                };
            Func<double[],double, int, double> bnd_valsR = delegate (double[] x, double t, int component)
            {
                if (t < perStartTime || wavePosition<shockPosition)
                {
                    if (component == 0)
                    {
                        return densityR ;
                    }
                    else if (component == 1)
                    {
                        return velocityR ;
                    }
                    else if (component == 2)
                    {
                        return pressureR ;
                    }
                    else
                    {
                        throw new NotSupportedException("component not supported");
                    }
                }
                else
                {
                    if (component == 0)
                    {
                        return densityR + DensityPertubation(new double[] { x[0], t - perStartTime });
                    }
                    else if (component == 1)
                    {
                        return velocityR + VelocityXPertubation(new double[] { x[0], t - perStartTime });
                    }
                    else if (component == 2)
                    {
                        return pressureR + PressurePertubation(new double[] { x[0], t - perStartTime });
                    }
                    else
                    {
                        throw new NotSupportedException("component not supported");
                    }
                }

            };




            // ### Initial condtions ###
            //    c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => init_vals(X,0));
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => init_vals(X, 1));
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => init_vals(X, 2));

            // ### Initial condtions ###
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => smooth_init_vals(X, 0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => smooth_init_vals(X, 1));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => smooth_init_vals(X, 2));

            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 1.0);
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 1.0);

            // ### Boundary condtions ###
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => X[0] <shockPosition? bnd_valsL(X, t, 0) : bnd_valsR(X, t, 0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => X[0] < shockPosition ? bnd_valsL(X, t, 1) : bnd_valsR(X, t, 1));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => X[0] < shockPosition ? bnd_valsL(X, t, 2) : bnd_valsR(X,t,2));
            //c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => pressureR);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1e-4;
            c.CFLFraction = CFLFraction;
            c.Endtime = endTime;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "AcousticWave";
            c.SessionName = String.Format("AW_p{0}_xCells{1}_yCells{2}_sP{3}_pST{4}_wP{5}_ampneg{6}_amppos{7}_wL{8}_Mach{9}", dgDegree, numOfCellsX, numOfCellsY,c.shockPosition,perStartTime,wavePosition,p_amp_neg,p_amp_pos, waveLength,MachL);
            //c.ProjectName = "StatShockRef";
            //c.SessionName = String.Format("StatShockRef_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_sP{5}_Mach{6}", dgDegree, numOfCellsX, numOfCellsY, c.CFLFraction, c.ExplicitOrder, c.shockPosition,MachL);

            return c;
        }

        /// <summary>
        /// Version to be submitted on the TU Darmstadt HHLR Lichtenberg cluster
        /// </summary>
        public static CNSControl ShockVortexInteractionHiOCFD5HLLR(int savePeriod = 1000, int dgDegree = 4, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 3, int explicitOrder = 3, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double Mv = 0.9, double Ms = 1.5, int numOfCellsX = 600, int numOfCellsY = 300, double endTime = 0.7) {

            // Lichtenberg
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_svi_video";
            //string restart = "False";

            CNSControl c = ShockVortexInteractionHiOCFD5(dbPath, savePeriod, dgDegree, sensorLimit, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, Mv, Ms, numOfCellsX, numOfCellsY, endTime);

            c.ProjectName = "svi_video";
            //c.Endtime = endTime;
            //c.NoOfTimesteps = 10;

            return c;
        }

        public static CNSControl DoubleMachReflection(string dbPath = null, int savePeriod = int.MaxValue, int dgDegree = 2, double xMax = 4.0, double yMax = 1.0, int numOfCellsX = 800, int numOfCellsY = 200, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 0.2, string restart = "False", int cores = int.MaxValue) {
            CNSControl c = new CNSControl();

            //dbPath = @"/work/scratch/yp19ysog/bosss_db_dmr_video";          // Lichtenberg
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

            // Dynamic load balancing
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax);
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

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

                    for (int iii = 0; iii < 1; iii++) {
                        Console.WriteLine("Setting edge tags (" + iii + ")...");
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
                        MPI.Wrappers.csMPI.Raw.Barrier(MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                        Console.WriteLine("done.");
                    }
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
                n.NormalizeInPlace();

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

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => 8.0 - SmoothJump(DistanceToInitialShock(X, t)) * (8.0 - 1.4));
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
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 8.0 - SmoothJump(DistanceToInitialShock(X, 0)) * (8.0 - 1.4));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, 0)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToInitialShock(X, 0)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 116.5 - SmoothJump(DistanceToInitialShock(X, 0)) * (116.5 - 1.0));
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = endTime;
            c.CFLFraction = CFLFraction;
            c.NoOfTimesteps = 10;

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

        static float eTilde(float x) {
            float x2 = x * x;
            float x4 = x2 * x2;
            if (x < 0)
                x *= -1.0f;
            return (1.0f + x + 0.5658f * x2 + 0.143f * x4);
        }

        static float FastTanh(float x) {
            // check Appendices C.1. in paper https://arxiv.org/pdf/1702.07825.pdf
            float _eTilde = eTilde(x);
            float eTilde1 = 1 / _eTilde;
            float epsilon = 0.000001f;
            int sign = (x > epsilon) ? 1 : (x > -epsilon) ? 0 : -1;
            return sign * (_eTilde - eTilde1) / (_eTilde + eTilde1);
        }

        public static IBMControl IBMShockTube(string dbPath = null, int savePeriod = 1, int dgDegree = 2, int numOfCellsX = 75, int numOfCellsY = 55, double sensorLimit = 1e-3, double dtFixed = 0.0, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double agg = 0.3, string restart = "False", double smoothing = 4.0) {
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
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
                    grid.DefineEdgeTags((Vector X) => 1);
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
                n.NormalizeInPlace();

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
                    c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
                    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));
                } else {
                    c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
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
            c.ProjectName = "IBM_Shock_Tube";

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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
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
                n.NormalizeInPlace();

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
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => densityLeft);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressureLeft);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicOutlet");

            Func<double, double> Jump = x => x <= discontinuityPosition ? 0 : 1;

            if (restart == "False") {
                if (rotatedGrid) {
                    c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (densityLeft - densityRight));
                    c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceToInitialShock(X, t: 0.0)) * (pressureLeft - pressureRight));
                } else {
                    c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
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
                grid.DefineEdgeTags((Vector X) => 1);
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
                n.NormalizeInPlace();

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

                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 25);    // fix lambdaMax
                                                                                                                                                                    //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, fudgeFactor: fugdeFactor);    // dynamic lambdaMax
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

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
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("bc2a0355-4e40-44fb-9409-4519cc3db797"), -1);
                c.GridGuid = new Guid("3e22691f-c635-472e-bb74-dcea3729fb74");
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
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => 8.0 - SmoothJump(X[0] - getShockXPosition(t)) * (8.0 - 1.4));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => 8.25 - SmoothJump(X[0] - getShockXPosition(t)) * (8.25 - 0.0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => 0.0);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => 116.5 - SmoothJump(X[0] - getShockXPosition(t)) * (116.5 - 1.0));

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet");
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 8.0 - SmoothJump(X[0] - getShockXPosition(0)) * (8.0 - 1.4));
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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 5);    // fix lambdaMax
                                                                                                                                                                   //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, fudgeFactor: fugdeFactor);    // dynamic lambdaMax
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

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
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory, no outflow boundary condition has to be specified, as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet");
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => density);
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
            // dbg_launch();

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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                if (lambdaMax == -1.0) { // dynamic lambdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                } else { // fixed lamdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);
                }
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

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
                        } else if ((X[0] - 0.0) * (X[0] - 0.0) + (X[1] - 0.0) * (X[1] - 0.0) - 1.5 * 1.5 >= 1e-14) { // just get all other boundary edges that are walls
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

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 0.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => density);
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

        public static IBMControl IBMBowShock(string dbPath = null, int savePeriod = 100, int dgDegree = 0, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 8.0, string restart = "False", int numOfCellsX = 20, int numOfCellsY = 80, double? lambdaMax = null) {
            IBMControl c = new IBMControl();

            //double? lambdaMax = 10;

            //System.Threading.Thread.Sleep(10000);
            // dbg_launch();

            //dbPath = @"/work/scratch/yp19ysog/bosss_db_dmr_video";          // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network
            //dbPath = @"H:\geisenhofer\bosss_db_bowShock";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1000;

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

            // Grid
            double xMin = -2.0;
            double xMax = 0.0;
            double yMin = -4.0;
            double yMax = 4.0;

            // Shift grid
            double h = Math.Abs(xMax - xMin) / numOfCellsX;
            xMin = xMin - h / 2;
            xMax = xMax - h / 2;

            if (restart == "True") {
                // Restart Lichtenberg
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("d369e78c-6f34-42e1-b8b9-e3903db132c2"), -1);
                c.GridGuid = new Guid("c691d970-6e52-4dd2-9d10-e95ab99f0482");
            } else {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                            if (Math.Abs(X[1]) - 0.1 < 1e-14) { // Right boundary (part of void area)
                                return 3;
                            } else {
                                return 2;
                            }
                        } else if (Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                            return 1;
                        } else {    // Top and bottom boundary
                            return 1;
                        }
                    });

                    //var gDat = new GridData(grid);
                    //var em1 = gDat.GetBoundaryEdges();
                    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                    return grid;
                };
            }

            // Level-set
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.LevelSetFunction = delegate (double[] X, double t) {
                // Circle 1
                double x0 = 0.0;
                double y0 = 0.5;
                double r0 = 0.5;

                // Circle 2
                double x1 = 0.0;
                double y1 = -0.5;
                double r1 = 0.5;

                // Signed distance formulation
                //if (X[1] >= 0.5) {
                //    return Math.Sqrt((X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0)) - r0;
                //} else if (X[1] <= -0.5) {
                //    return Math.Sqrt((X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1)) - r1;
                //} else {
                //    return -(X[0] + 0.5);
                //}

                // Quadratic formulation
                if (X[1] >= 0.5) {
                    return (X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0) - r0 * r0;
                } else if (X[1] <= -0.5) {
                    return (X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1) - r1 * r1;
                } else {
                    return X[0] * X[0] - 0.5 * 0.5;
                }
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            int levelSetDegree = 2;
            c.LevelSetQuadratureOrder = 3 * levelSetDegree;
            c.AgglomerationThreshold = 0.3;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, levelSetDegree);

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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                if (lambdaMax == null) { // dynamic lambdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                } else { // fixed lamdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);
                }
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Enthalpy, dgDegree);

            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            if (AV) {
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            c.AddVariable(CNSVariables.CFL, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;
            double Mach = 4.0;
            double velocityX = Mach * Math.Sqrt(c.EquationOfState.HeatCapacityRatio * pressure / density);
            double velocityY = 0.0;

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 0.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => density);
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
            //c.dtFixed = 1e-3;

            c.ProjectName = "IBMBowShock";

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                if (dgDegree == 0) {
                    tempSessionName = string.Format("IBMBowShock_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_agg{5}",
                        dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, c.AgglomerationThreshold);
                } else {
                    tempSessionName = string.Format("IBMBowShock_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_s0={5:0.0E-00}_lambdaMax{6}_agg{7}_RESTART12",
                        dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, sensorLimit, lambdaMax, c.AgglomerationThreshold);
                }
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("IBMBowShock_p{0}_s0={1:0.0E-00}_CFLFrac{2}_AB{3}",
                    dgDegree, sensorLimit, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("IBMBowShock_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_ALTS{4}_{5}_re{6}_subs{7}_s0={8:0.0E-00}_lambdaMax{9}",
                    dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, sensorLimit, lambdaMax);
            }
            c.SessionName = tempSessionName;

            return c;
        }

        public static IBMControl IBMBowShockHPC(string dbPath = null, int savePeriod = 100, int dgDegree = 2, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 16.0, string restart = "False", int numOfCellsX = 80, int numOfCellsY = 320, double? lambdaMax = 15.0) {
            IBMControl c = new IBMControl();

            //double? lambdaMax = 10;

            //System.Threading.Thread.Sleep(10000);
            // dbg_launch();

            //dbPath = @"/work/scratch/yp19ysog/bosss_db_dmr_video";          // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network
            //dbPath = @"H:\geisenhofer\bosss_db_bowShock";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 100;

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

            // Grid
            double xMin = -2.0;
            double xMax = 0.0;
            double yMin = -4.0;
            double yMax = 4.0;

            // Shift grid
            double h = Math.Abs(xMax - xMin) / numOfCellsX;
            xMin = xMin - h / 2;
            xMax = xMax - h / 2;

            if (restart == "True") {
                // Restart Lichtenberg
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("a6b7d1f0-2ef5-4cca-8aca-b92228aad095"), -1);
                c.GridGuid = new Guid("2febbaeb-f611-4975-aa98-3288019507c4");
            } else {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[0] - xMax) < 1e-14) {    // Right boundary
                            if (Math.Abs(X[1]) - 0.1 < 1e-14) { // Right boundary (part of void area)
                                return 3;
                            } else {
                                return 2;
                            }
                        } else if (Math.Abs(X[0] - xMin) < 1e-14) { // Left boundary
                            return 1;
                        } else {    // Top and bottom boundary
                            return 1;
                        }
                    });

                    //var gDat = new GridData(grid);
                    //var em1 = gDat.GetBoundaryEdges();
                    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                    return grid;
                };
            }

            // Level-set
            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.LevelSetFunction = delegate (double[] X, double t) {
                // Circle 1
                double x0 = 0.0;
                double y0 = 0.5;
                double r0 = 0.5;

                // Circle 2
                double x1 = 0.0;
                double y1 = -0.5;
                double r1 = 0.5;

                // Signed distance formulation
                //if (X[1] >= 0.5) {
                //    return Math.Sqrt((X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0)) - r0;
                //} else if (X[1] <= -0.5) {
                //    return Math.Sqrt((X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1)) - r1;
                //} else {
                //    return -(X[0] + 0.5);
                //}

                // Quadratic formulation
                if (X[1] >= 0.5) {
                    return (X[0] - x0) * (X[0] - x0) + (X[1] - y0) * (X[1] - y0) - r0 * r0;
                } else if (X[1] <= -0.5) {
                    return (X[0] - x1) * (X[0] - x1) + (X[1] - y1) * (X[1] - y1) - r1 * r1;
                } else {
                    return X[0] * X[0] - 0.5 * 0.5;
                }
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            int levelSetDegree = 2;
            c.LevelSetQuadratureOrder = 3 * levelSetDegree;
            c.AgglomerationThreshold = 0.3;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, levelSetDegree);

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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                if (lambdaMax == null) { // dynamic lambdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                } else { // fixed lamdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);
                }
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Enthalpy, dgDegree);

            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            if (AV) {
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            c.AddVariable(CNSVariables.CFL, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;
            double Mach = 4.0;
            double velocityX = Mach * Math.Sqrt(c.EquationOfState.HeatCapacityRatio * pressure / density);
            double velocityY = 0.0;

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 0.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False") {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => density);
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
            //c.dtFixed = 1e-3;

            c.ProjectName = "IBMBowShock";

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                if (dgDegree == 0) {
                    tempSessionName = string.Format("IBMBowShock_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_agg{5}",
                        dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, c.AgglomerationThreshold);
                } else {
                    tempSessionName = string.Format("IBMBowShock_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_s0={5:0.0E-00}_lambdaMax{6}_agg{7}_RESTART4",
                        dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, sensorLimit, lambdaMax, c.AgglomerationThreshold);
                }
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("IBMBowShock_p{0}_s0={1:0.0E-00}_CFLFrac{2}_AB{3}",
                    dgDegree, sensorLimit, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("IBMBowShock_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_ALTS{4}_{5}_re{6}_subs{7}_s0={8:0.0E-00}_lambdaMax{9}",
                    dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, sensorLimit, lambdaMax);
            }
            c.SessionName = tempSessionName;

            return c;
        }

        /// <summary>
        /// Version to be submitted on the TU Darmstadt HHLR Lichtenberg cluster
        /// </summary>
        public static IBMControl IBMBowShockHHLR(int savePeriod, int dgDegree, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps, double endTime, int numOfCellsX, int numOfCellsY, double? lambdaMax) {

            // Lichtenberg
            string dbPath = @"/work/scratch/yp19ysog/bosss_db_bowShock";
            //string dbPath = @"/work/scratch/yp19ysog/bosss_db_ibmbowshock";
            //string dbPath = @"H:\geisenhofer\bosss_db_ibmbowshock";
            string restart = "True";

            IBMControl c = IBMBowShock(dbPath, savePeriod, dgDegree, sensorLimit, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, endTime, restart, numOfCellsX, numOfCellsY, lambdaMax);

            //c.AlternateDbPaths = new ValueTuple<string, string>[] { (@"S:\work\scratch\yp19ysog\bosss_db_bowShock", "pcmit33"), (dbPath, "") };

            //c.TracingNamespaces = "BoSSS.Solution";

            c.ProjectName = "IBMBowShock_P3";
            //c.ProjectName = "ibmbowshock_hhlr";
            //c.ProjectName = "IBMBowShock_P2_FineGrid";

            //c.NoOfTimesteps = 10;

            return c;
        }

        public static CNSControl StationaryShockWave(string dbPath = null, int savePeriod = 1000, int dgDegree = 2, int numOfCellsX = 20, int numOfCellsY = 6, double dtFixed = 1e-8) {
            CNSControl c = new CNSControl();

            // ### Database ###
            //dbPath = @"c:\bosss_db";
            dbPath = @"d:\bosss_db_ssw";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            // ### Time-Stepping ###
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // ### Physics ###
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            bool AV = false;

            // #### Activate ARTIFICIAL VISCOSITY
            if (dgDegree > 0) {
                AV = true;
            }

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }

            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            //### Shock-capturing ###
            double epsilon0 = 1.0;
            double kappa = 0.5;
            double sensorLimit = 1e-5;
            if (AV) {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);    // fix lambdaMax
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);    // dynamic lambdaMax
            }

            // ### Output variables ###
            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);

            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // ### Grid ###
            double xMin = 0;
            double xMax = 1.0;
            double yMin = 0;
            double yMax = 0.3;

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicInlet");
                //grid.EdgeTagNames.Add(2, "SubsonicOutlet");
                //grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] X) {
                    if (Math.Abs(X[1]) < 1e-14) {   // bottom
                        return 3;
                    } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {    // top
                        return 3;
                    } else if (Math.Abs(X[0]) < 1e-14) {                    // left
                        return 1;
                    } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {    // right
                        return 1;
                    } else {
                        throw new Exception("Boundary condition not specified");
                    }
                });

                return grid;
            };

            // ### Initial conditions ###

            // #########################
            // Shock
            // #########################
            //double Ms = 1.5;
            //double densityLeft = 1;
            //double pressureLeft = 1;
            //double velocityXLeft = Ms * Math.Sqrt(gamma * pressureLeft / densityLeft);
            //double velocityYLeft = 0;

            // #########################
            // Parameters
            // #########################
            double gamma = IdealGas.Air.HeatCapacityRatio;
            double Ms = 1.5;

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
            double velocityYLeft = 0.0;
            double velocityYRight = 0.0;

            double cellSize = (xMax - xMin) / numOfCellsX;

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            //Func<double, double> SmoothJump = delegate (double distance) {
            //    if (distance <= 0) {
            //        return 0;
            //    } else {
            //        return 1;
            //    }
            //};

            double shockPosition = 0.525;
            //double shockPosition = 0.5;

            //densityRight = densityRight * 1.10;
            //velocityXRight = velocityXRight * 1.10;
            //pressureRight = pressureRight * 1.10;

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

            // Initial conditions
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => X[0] <= shockPosition ? densityLeft : densityRight);
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => X[0] <= shockPosition ? velocityXLeft : velocityXRight);
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => X[0] <= shockPosition ? velocityYLeft : velocityYRight);
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => X[0] <= shockPosition ? pressureLeft : pressureRight);

            double momentumXLeft = densityLeft * velocityXLeft;
            double momentumXRight = densityRight * velocityXRight;
            double momentumYLeft = 0.0;
            double momentumYRight = 0.0;

            double innerEnergyLeft = pressureLeft / (gamma - 1);
            double innerEnergyRight = pressureRight / (gamma - 1);
            double kineticEnergyLeft = 0.5 * densityLeft * (velocityXLeft * velocityXLeft + velocityYLeft * velocityYLeft);
            double kineticEnergyRight = 0.5 * densityRight * (velocityXRight * velocityXRight + velocityYRight * velocityYRight);
            double totalEnergyLeft = innerEnergyLeft + kineticEnergyLeft;
            double totalEnergyRight = innerEnergyRight + kineticEnergyRight;

            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => X[0] <= shockPosition ? densityLeft : densityRight);
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => X[0] <= shockPosition ? velocityXLeft : velocityXRight);
            //c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => X[0] <= shockPosition ? velocityYLeft : velocityYRight);
            //c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => X[0] <= shockPosition ? pressureLeft : pressureRight);

            //c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => X[0] <= shockPosition ? densityLeft : densityRight);
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.xComponent, X => X[0] <= shockPosition ? momentumXLeft : momentumXRight);
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Momentum.yComponent, X => X[0] <= shockPosition ? momentumYLeft : momentumYRight);
            //c.InitialValues_Evaluators.Add(CompressibleVariables.Energy, X => X[0] <= shockPosition ? totalEnergyLeft : totalEnergyRight);

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => DensityShock(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => VelocityXShock(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => VelocityYShock(X));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => PressureShock(X));

            // ### Boundary conditions ###
            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => X[0] <= shockPosition ? densityLeft : densityRight);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => X[0] <= shockPosition ? velocityXLeft : velocityXRight);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => X[0] <= shockPosition ? velocityYLeft : velocityYRight);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => X[0] <= shockPosition ? pressureLeft : pressureRight);
            //c.AddBoundaryValue("SubsonicOutlet", CNSVariables.Pressure, (X, t) => X[0] <= shockPosition ? pressureLeft : pressureRight);
            //c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => X[0] <= shockPosition ? pressureLeft : pressureRight);
            c.AddBoundaryValue("AdiabaticSlipWall");


            // ### Time configuration ###
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.3;
            //c.dtFixed = 1e-6;

            c.Endtime = 2.0;
            c.NoOfTimesteps = int.MaxValue;

            // ### Project and sessions name ###
            c.ProjectName = "CNS_StationaryShockWave";

            string tempSessionName;

            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta) {
                tempSessionName = String.Format("CNS_SSW_p{0}_xCells{1}_yCells{2}_RK{3}_CFLFrac={6}_tEnd={5}_S0={4}_shockPos={7}", dgDegree, numOfCellsX, numOfCellsY, c.ExplicitOrder, sensorLimit, c.Endtime, c.CFLFraction, shockPosition);
            } else {
                throw new NotImplementedException("Session name is not available for this type of time stepper");
            }

            c.SessionName = tempSessionName;

            return c;
        }

        public static CNSControl ContactDiscontinuity(double levelSetPosition = 0.45, string dbPath = null, int dgDegree = 2, int numOfCellsX = 20, int numOfCellsY = 5, double sensorLimit = 1, bool saveToDb = false) {

            CNSControl c = new CNSControl();

            dbPath = @"c:\bosss_db";
            //dbPath = @"\\fdyprime\userspace\geisenhofer\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = 1;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            bool AV = true;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 0.5;

            Variable sensorVariable = CompressibleVariables.Density;
            c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);

            if (AV) {
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
            }

            // Runge-Kutta
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 1;
            //c.NumberOfSubGrids = 4;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

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
            c.AddVariable(CNSVariables.ShockSensor, 0);

            if (AV) {
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

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
                n.NormalizeInPlace();

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
            c.dtFixed = 1.0e-5;
            //c.CFLFraction = 0.3;
            //c.Endtime = 0.05;
            //c.NoOfTimesteps = int.MaxValue;
            c.NoOfTimesteps = 500;

            c.ProjectName = "Contact discontinuity";
            c.SessionName = String.Format("Contact discontinuity, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, CFLFraction = {3:0.00E-00}, ALTS {4}/{5}", dgDegree, numOfCellsX, numOfCellsY, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);

            return c;
        }

        public static CNSControl IBMObliqueShockTube(string dbPath = null, int dgDegree = 2, int numOfCellsX = 50, int numOfCellsY = 50, double sensorLimit = 1e-3, bool saveToDb = true) {
            IBMControl c = new IBMControl();

            dbPath = @"c:\bosss_db";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = 100;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1.3;
            double yMin = 0;
            double yMax = 0.65;
            double shockPosition = 0.8;

            bool AV = true;

            c.DomainType = DomainTypes.StaticImmersedBoundary;

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            Variable sensorVariable = CompressibleVariables.Density;
            c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            c.AddVariable(CNSVariables.ShockSensor, 0);

            if (AV) {
                double epsilon0 = 1.0;
                double kappa = 0.5;
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            }

            // Level set
            double beta = Math.PI / 16;
            double[] startOfLine1 = new double[] { 0.0, 0.05 };
            double distanceBetweenLines = 0.3;
            double[] startOfLine2 = new double[] { 0.0, startOfLine1[1] + distanceBetweenLines * Math.Cos(beta) };

            Func<double, double> line1 = delegate (double x) {
                return Math.Tan(beta) * (x - startOfLine1[0]) + startOfLine1[1];
            };
            Func<double, double> line2 = delegate (double x) {
                return Math.Tan(beta) * (x - startOfLine2[0]) + startOfLine2[1];
            };

            c.LevelSetFunction = delegate (double[] X, double t) {
                return -(X[1] - line1(X[0])) * (X[1] - line2(X[0]));
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.9;
            c.SaveAgglomerationPairs = true;
            c.AddVariable(IBMVariables.LevelSet, 2);

            // Runge-Kutta schemes
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            //c.ExplicitOrder = 3;
            //c.NumberOfSubGrids = 3;
            //c.ReclusteringInterval = 1;
            //c.FluxCorrection = false;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);

            if (AV) {
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });

                return grid;
            };

            c.AddBoundaryValue("AdiabaticSlipWall");

            double crossProduct2D(double[] a, double[] b) {
                return a[0] * b[1] - a[1] * b[0];
            }

            // Normal vector of initial shock
            Vector p1 = new Vector(0.5, line1(0.5));
            Vector p2 = new Vector(0.6, line1(0.6));
            Vector normalVector = p2 - p1;

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.NormalizeInPlace();

            // Distance from a point X to the initial shock
            double[] p = new double[] { shockPosition, line1(shockPosition) };

            double DistanceFromPointToLine(double[] X, double[] pointOnLine, double[] directionVector) {
                double[] X_minus_pointOnLine = new double[] { X[0] - pointOnLine[0], X[1] - pointOnLine[1] };
                double distance = crossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

                return distance;
            }

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            // Initial conditions
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.6;
            c.Endtime = 0.25;
            c.NoOfTimesteps = int.MaxValue;
            //c.NoOfTimesteps = 2;

            c.ProjectName = "IBM oblique shock tube";
            c.SessionName = String.Format("IBM oblique shock tube, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);

            return c;
        }


        public static CNSControl CylindricalExplosion1D(string dbPath = @"e:\bosss_db\GridOfTomorrow", int numCells = 1000) {
            int dgDegree = 0;

            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1000;
            c.PrintInterval = 10;
            c.GridPartType = GridPartType.none;

            c.ActiveOperators = Operators.Convection | Operators.CustomSource;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.CustomContinuitySources.Add(
                map => new AdHocSourceTerm(map, (X, t, state) =>
                    1.0 / X[0] * state.Momentum[0]));
            c.CustomMomentumSources[0].Add(
                map => new AdHocSourceTerm(map, (X, t, state) =>
                    1.0 / X[0] * state.Momentum[0] * state.Velocity[0]));
            c.CustomEnergySources.Add(
                map => new AdHocSourceTerm(map, (X, t, state) =>
                    1.0 / X[0] * state.Velocity[0] * (state.Energy + state.Pressure)));

            c.ProjectName = "Cylindrical explosion 1D";
            c.SessionName = String.Format(
                "Cylindrical explosion 1D, numCells={1}", dgDegree, numCells);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(0.0, 1.0, numCells + 1);
                var grid = Grid1D.LineGrid(nodes);

                grid.EdgeTagNames.Add(1, "SupersonicOutlet");
                grid.DefineEdgeTags((Vector X) => 1);

                return grid;
            };

            c.AddBoundaryValue("SupersonicOutlet");

            double R = 0.4;

            double rhoIn = 1.0;
            double rhoOut = 0.125;
            double pIn = 1.0;
            double pOut = 0.1;

            Func<double, double> Jump = (x => x < 0.0 ? 0.0 : 1.0);
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rhoIn - Jump(X[0] - R) * (rhoIn - rhoOut));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pIn - Jump(X[0] - R) * (pIn - pOut));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }

        public static CNSControl SphericalExplosion(string dbPath = @"e:\bosss_db\GridOfTomorrow", int dgDegree = 0, int numCells = 10, double sensorLimit = 1e-2) {
            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 25;
            c.PrintInterval = 1;
            c.GridPartType = GridPartType.none;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.ProjectName = "Spherical explosion";
            c.SessionName = String.Format(
                "Spherical explosion, degree={0}, numCells={1}, sensorLimit = {2:0.00E-00}", dgDegree, numCells, sensorLimit);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-1.0, 1.0, numCells + 1);
                var grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes);

                grid.EdgeTagNames.Add(1, "SupersonicOutlet");
                grid.DefineEdgeTags((Vector X) => 1);

                return grid;
            };

            double cellSize = 2.0 / numCells;
            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.5 * cellSize / Math.Max(dgDegree, 1);
                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            c.AddBoundaryValue("SupersonicOutlet");

            double radius = 0.4;
            Func<double[], double> levelSet = X => Math.Sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]) - radius;

            double rhoIn = 1.0;
            double rhoOut = 0.125;
            double pIn = 1.0;
            double pOut = 0.1;
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rhoIn - SmoothJump(levelSet(X)) * (rhoIn - rhoOut));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.zComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pIn - SmoothJump(levelSet(X)) * (pIn - pOut));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            // Shock-capturing
            {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ActiveOperators |= Operators.ArtificialViscosity;

                double epsilon0 = 1.0;
                double kappa = 1.0;
                //double lambdaMax = 2.0;    // determined manually in Visit
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(
                    c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);

                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 3);

                c.Tags.Add("Artificial viscosity");

                c.PrandtlNumber = 0.71;
                c.ReynoldsNumber = 1.0;
            }

            return c;
        }

        public static CNSControl CylindricalExplosion(string dbPath = @"e:\bosss_db\cylindricalExplosion\", int dgDegree = 0, int numCells = 20, double sensorLimit = 1e-2) {
            CNSControl c = new CNSControl();

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 25;
            c.PrintInterval = 1;
            c.GridPartType = GridPartType.none;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.ProjectName = "Cylindrical explosion";
            c.SessionName = String.Format(
                "Cylindrical explosion, degree={0}, numCells={1}, sensorLimit = {2:0.00E-00}", dgDegree, numCells, sensorLimit);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.CFL, 0);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-1.0, 1.0, numCells + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: false, periodicY: false);

                grid.EdgeTagNames.Add(1, "SupersonicOutlet");
                grid.DefineEdgeTags((Vector X) => 1);

                return grid;
            };

            double cellSize = 2.0 / numCells;
            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 1.5 * cellSize / Math.Max(dgDegree, 1);
                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            c.AddBoundaryValue("SupersonicOutlet");

            double radius = 0.4;
            Func<double[], double> levelSet = X => Math.Sqrt(X[0] * X[0] + X[1] * X[1]) - radius;

            double rhoIn = 1.0;
            double rhoOut = 0.125;
            double pIn = 1.0;
            double pOut = 0.1;
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => rhoIn - SmoothJump(levelSet(X)) * (rhoIn - rhoOut));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pIn - SmoothJump(levelSet(X)) * (pIn - pOut));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            // Shock-capturing
            {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ActiveOperators |= Operators.ArtificialViscosity;

                double epsilon0 = 1.0;
                double kappa = 1.0;

                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(
                    c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);

                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);

                c.Tags.Add("Artificial viscosity");

                c.PrandtlNumber = 0.71;
                c.ReynoldsNumber = 1.0;
            }

            return c;
        }
        public static IBMControl IBMWedgeFlow(string dbPath = null, int PrintInterval = 1000,int savePeriod = 100, int dgDegree = 0, 
            double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, 
            int numberOfSubGrids = 2, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 8.0, string restart = "False", 
            int numOfCellsX = 20, int numOfCellsY = 80, double? lambdaMax = null)
        {
            IBMControl c = new IBMControl();

            //double? lambdaMax = 10;

            //System.Threading.Thread.Sleep(10000);
            // dbg_launch();

            //dbPath = @"/work/scratch/yp19ysog/bosss_db_dmr_video";          // Lichtenberg
            //dbPath = @"c:\bosss_db";                                          // Local
            //dbPath = @"E:\geisenhofer\bosss_db_paper_ibmdmr";                   // HPC cluster
            //dbPath = @"\\dc1\userspace\geisenhofer\bosss_db_IBMShockTube";    // Network
            //dbPath = @"H:\geisenhofer\bosss_db_bowShock";

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = savePeriod;
            c.PrintInterval = PrintInterval;

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

            // Grid
            double xMin = 0;
            double xMax = 1.5;
            double yMin = 0;
            double yMax = 1.0;
            double wedge_angle = 10.0 * Math.PI / 180.0; // Angle of Wedge

            if (restart == "True")
            {
                // Restart Lichtenberg
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid("d369e78c-6f34-42e1-b8b9-e3903db132c2"), -1);
                c.GridGuid = new Guid("c691d970-6e52-4dd2-9d10-e95ab99f0482");
            }
            else
            {
                c.GridFunc = delegate {
                    double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                    if (numOfCellsX > 1 && xNodes[1] >= 0.5)
                    {
                        xNodes[1] = 0.4;
                    }
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);

                    grid.EdgeTagNames.Add(1, "SupersonicInlet");
                    grid.EdgeTagNames.Add(2, "SupersonicOutlet");
                    grid.EdgeTagNames.Add(3, "AdiabaticSlipWall");
                    //grid.EdgeTagNames.Add(4, "SubsonicOutlet");


                    grid.DefineEdgeTags(delegate (double[] X) {
                        if (Math.Abs(X[0] - xMax) < 1e-14)
                        {    // Right boundary
                             ////if((-X[0] + 0.5 + (X[1] / Math.Tan(wedge_angle)))>=0){ 
                             //if(X[1] >= Math.Tan(39.3139318 * Math.PI / 180.0)) { //-0.01 doesn't work
                             //    return 2;
                             //} else { // (part of void area)
                            return 2;
                            //return 3; // Wedge Part part of Outlet

                        }
                        else if (Math.Abs(X[0] - xMin) < 1e-14)
                        { // Left boundary
                            return 1;
                        }
                        else if (Math.Abs(X[1] - yMax) < 1e-14)
                        { // top boundary
                            return 2;
                            //if(X[0] < 0.5) {
                            //    return 1;
                            //} else {
                            //    return 2;
                            //}

                        }
                        else
                        { //if(Math.Abs(X[1] - yMin) < 1e-14) { //bottom boundary
                            return 3;
                            //if(X[0] <= 1) {
                            //    return 2;
                            //} else {
                            //    return 2;
                            //}
                        }
                    });
                    return grid;
                };
            }

            // Level-set
            c.DomainType = DomainTypes.StaticImmersedBoundary;

            //variables
            double Mach = 2.0;
            double LevelSet1Prime = Math.Tan(wedge_angle);

            c.LevelSetFunction = delegate (double[] X, double t) {
                return -X[0] + 0.5 + (X[1] / LevelSet1Prime);
            };

            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            int levelSetDegree = 1;
            c.LevelSetQuadratureOrder = 3 * levelSetDegree;
            c.AgglomerationThreshold = 0.3;
            c.SaveAgglomerationPairs = false;
            c.AddVariable(IBMVariables.LevelSet, levelSetDegree);

            bool AV;
            if (dgDegree > 0)
            {
                AV = true;
            }
            else
            {
                AV = false;
            }

            if (AV)
            {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            }
            else
            {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double epsilon0 = 1.0;
            double kappa = 1.0;

            if (AV)
            {
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                if (lambdaMax == null)
                { // dynamic lambdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                }
                else
                { // fixed lamdaMax
                    c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: lambdaMax);
                }
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Enthalpy, dgDegree);

            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);

            if (AV)
            {
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            c.AddVariable(CNSVariables.CFL, 0);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS))
            {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            // Boundary conditions
            double density = 1.0;
            double pressure = 1.0;


            double gamma = IdealGas.Air.HeatCapacityRatio;
            double speedOfSound = Math.Sqrt(gamma * pressure / density);

            double velocityX = Mach * speedOfSound;
            double velocityY = 0.0;

            double momentumX = density * velocityX;


            double innerEnergy = pressure / (gamma - 1);
            double kineticEnergy = 0.5 * density * (velocityX * velocityX + velocityY * velocityY);
            double totalEnergy = innerEnergy + kineticEnergy;

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => density);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => velocityX);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => velocityY);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => pressure);

            // In theory no outflow boundary condition has to be specified as all characteristics move downstream
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 0.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (restart == "False")
            {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => density);
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
            //c.dtFixed = 1e-3;

            c.ProjectName = "IBMWedgeFlow";

            // Session name
            string tempSessionName;
            if (c.ExplicitScheme == ExplicitSchemes.RungeKutta)
            {
                if (dgDegree == 0)
                {
                    tempSessionName = string.Format("IBMWedgeFlow_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_agg{5}",
                        dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, c.AgglomerationThreshold);
                }
                else
                {
                    tempSessionName = string.Format("IBMWedgeFlow_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_RK{4}_s0={5:0.0E-00}_lambdaMax{6}_agg{7}_RESTART12",
                        dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, sensorLimit, lambdaMax, c.AgglomerationThreshold);
                }
            }
            else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth)
            {
                tempSessionName = string.Format("IBMWedgeFlow_p{0}_s0={1:0.0E-00}_CFLFrac{2}_AB{3}",
                    dgDegree, sensorLimit, CFLFraction, explicitOrder);
            }
            else
            {
                tempSessionName = string.Format("IBMWedgeFlow_p{0}_xCells{1}_yCells{2}_CFLFrac{3}_ALTS{4}_{5}_re{6}_subs{7}_s0={8:0.0E-00}_lambdaMax{9}",
                    dgDegree, numOfCellsX, numOfCellsY, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, sensorLimit, lambdaMax);
            }
            c.SessionName = tempSessionName;

            return c;
        }
    }
}