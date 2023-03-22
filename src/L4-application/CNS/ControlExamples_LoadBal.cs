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
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.LoadBalancing;
using CNS.ShockCapturing;
using CNS.Source;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Diagnostics;

namespace CNS {

    /// <summary>
    /// Test cases for the CNS solver with the focus on
    /// load balancing
    /// </summary>
    public static class ControlExamples_LoadBal {

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

        public static CNSControl DMR_Cube(int savePeriod = int.MaxValue, int dgDegree = 2, double xMax = 4.0, double yMax = 1.0, int NoCellsX_percore = 800, int NoCellsY_percore = 200, double sensorLimit = 1e-3, double CFLFraction = 0.1, int explicitScheme = 1, int explicitOrder = 1, int numberOfSubGrids = 3, int reclusteringInterval = 1, int maxNumOfSubSteps = 0, double endTime = 0.2, string restart = "False", int cores = int.MaxValue) {
            CNSControl c = new CNSControl();

            c.WriteLTSLog = false;
            c.WriteLTSConsoleOutput = false;

            //c.TracingNamespaces = "BoSSS.Foundation.Grid.Classic";

            double xMin = 0;
            double yMin = 0;

            //Partitioning
            int[] separation = new int[] { 1, 1 };
            switch (cores) {
                case 4:
                separation = new int[] { 2, 2 };
                break;
                case 8:
                separation = new int[] { 4, 2 };
                break;
                case 16:
                separation = new int[] { 4, 4 };
                break;
                case 32:
                separation = new int[] { 8, 4 };
                break;
                case 64:
                separation = new int[] { 8, 8 };
                break;
                case 128:
                separation = new int[] { 16, 8 };
                break;
                case 256:
                separation = new int[] { 16, 16 };
                break;
                default:
                c.GridPartType = GridPartType.none;
                break;
            }
            c.GridPartType = GridPartType.Predefined;
            c.GridPartOptions = "hallo";
            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            Func<double[], int> MakeMyPartioning = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                double xspan = (xMax - xMin) / separation[0];
                double yspan = (yMax - yMin) / separation[1];
                int rank = int.MaxValue;
                int icore = 0;
                for (int i = 0; i < separation[0]; i++) {
                    for (int j = 0; j < separation[1]; j++) {
                        bool xtrue = x <= xspan * (i + 1) + xMin;
                        bool ytrue = y <= yspan * (j + 1) + yMin;
                        if (xtrue && ytrue) {
                            rank = icore;
                            return rank;
                        }
                        icore++;
                    }
                }

                return rank;
            };
            //get total number of cells for each direction of space
            int numOfCellsX = NoCellsX_percore * separation[0];
            int numOfCellsY = NoCellsY_percore * separation[1];

            // Time stepping
            c.ExplicitScheme = (ExplicitSchemes)explicitScheme;
            c.ExplicitOrder = explicitOrder;
            c.NumberOfSubGrids = numberOfSubGrids;
            c.ReclusteringInterval = reclusteringInterval;
            c.maxNumOfSubSteps = maxNumOfSubSteps;
            c.FluxCorrection = false;

            // Dynamic load balancing
            c.DynamicLoadBalancing_On = false;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = c.ReclusteringInterval;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            const double xWall = 0.2;
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
            c.AddVariable(CNSVariables.Rank, 0);

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
                    grid.AddPredefinedPartitioning("hallo", MakeMyPartioning);

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

            const double tan60 = 1.732050807568877;
            const double sin60 = 8.660254037844386e-01;// Math.Sin(Math.PI / 3);
            const double cos60 = 0.5;

            Debug.Assert((sin60 - Math.Sin(Math.PI / 3).Abs() < BLAS.MachineEps * 100));
            Debug.Assert((cos60 - Math.Cos(Math.PI / 3).Abs() < BLAS.MachineEps * 100));
            Debug.Assert((tan60 - Math.Tan(Math.PI / 3).Abs() < BLAS.MachineEps * 100));

            double DistanceToInitialShock(double[] X, double t) {
                //OptimizedHLLCFlux.DistanceToInitialShock.Start();
                // direction vector
                //Vector p1 = new Vector(xWall, 0.0);
                //Vector p2 = new Vector(xWall + 1 / tan60, 1.0);
                //Vector p = p2 - p1;
                Vector p = new Vector(1 / tan60, 1);

                // normal vector
                Vector n = new Vector(p.y, -p.x);
                n.NormalizeInPlace();

                // Angle between line and x-axis
                //double alpha = Math.Atan(Math.Abs((p2.y - p1.y)) / Math.Abs((p2.x - p1.x)));
                //double alpha = Math.PI / 3;

                // distance of a point X to the origin (normal to the line)
                double nDotX = n * X;

                // shock speed
                double vs = 10;

                // distance to line
                double distance = nDotX - (sin60 * xWall + vs * t);

                //OptimizedHLLCFlux.DistanceToInitialShock.Stop();
                return distance;
            }

            // Function for smoothing the initial and top boundary conditions
            double SmoothJump(double distance) {
                //OptimizedHLLCFlux.SmoothJump.Start();
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);

                //double retval = (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
                double retval = (FastTanh((float)(distance / maxDistance)) + 1.0) * 0.5; // ca 20%
                //OptimizedHLLCFlux.SmoothJump.Stop();
                return retval;
            }

            // Function for a sharp jump (no smoothing of initial and top boundary conditions)
            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            // Boundary conditions
            //c.AddBoundaryValue("SupersonicInlet", Variables.Density, (X, t) => 8.0 - Jump(X[0] - (0.1 + (X[1] + 20 * t) / 1.732)) * (8.0 - 1.4));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => 7.14471 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (7.14471 - 0.0));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => -4.125 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (-4.125 - 0.0));
            //c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => 116.5 - Jump(X[0] - (0.1 + (X[1] + 20.0 * t) / 1.732)) * (116.5 - 1.0));

            double DensityInlet(double[] X, double t) {
                return 8.0 - SmoothJump(DistanceToInitialShock(X, t)) * (8.0 - 1.4);
            }
            double VelocityXInlet(double[] X, double t) {
                return 8.25 * sin60 - SmoothJump(DistanceToInitialShock(X, t)) * (8.25 * sin60 - 0.0);
            }
            double VelocityYInlet(double[] X, double t) {
                return -8.25 * cos60 - SmoothJump(DistanceToInitialShock(X, t)) * (-8.25 * cos60 - 0.0);
            }
            double PressureInlet(double[] X, double t) {
                return 116.5 - SmoothJump(DistanceToInitialShock(X, t)) * (116.5 - 1.0);
            }

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, DensityInlet);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, VelocityXInlet);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, VelocityYInlet);
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, PressureInlet);

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
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 8.25 * sin60 - SmoothJump(DistanceToInitialShock(X, 0)) * (8.25 * sin60 - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => -8.25 * cos60 - SmoothJump(DistanceToInitialShock(X, 0)) * (-8.25 * cos60 - 0.0));
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
                    dgDegree, NoCellsX_percore, NoCellsY_percore, sensorLimit, CFLFraction, explicitOrder, cores);
            } else if (c.ExplicitScheme == ExplicitSchemes.AdamsBashforth) {
                tempSessionName = string.Format("DMR_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_AB{5}",
                    dgDegree, NoCellsX_percore, NoCellsY_percore, sensorLimit, CFLFraction, explicitOrder);
            } else {
                tempSessionName = string.Format("DMR_p{0}_xCells{1}_yCells{2}_s0={3:0.0E-00}_CFLFrac{4}_ALTS{5}_{6}_re{7}_subs{8}",
                    dgDegree, NoCellsX_percore, NoCellsY_percore, sensorLimit, CFLFraction, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps);
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
        public static CNSControl DoubleMachReflectionHHLR(int savePeriod, int dgDegree, double xMax, double yMax, int numOfCellsX, int numOfCellsY, double sensorLimit, double CFLFraction, int explicitScheme, int explicitOrder, int numberOfSubGrids, int reclusteringInterval, int maxNumOfSubSteps, double endTime, int timeSteps) {
            //Absturz mit 128 cores: 
            //--control "cs:CNS.TestCases.DoubleMachReflectionHHLR(2147483647, 2, 4, 1, 1280, 320, 0.001, 0.1, 1, 1, 3, 1, 0, 0.7, 100)"

            // Lichtenberg

            string restart = "False";
            int cores = ilPSP.Environment.MPIEnv.MPI_Size;

            CNSControl c = DMR_Cube(savePeriod, dgDegree, xMax, yMax, numOfCellsX, numOfCellsY, sensorLimit, CFLFraction, explicitScheme, explicitOrder, numberOfSubGrids, reclusteringInterval, maxNumOfSubSteps, endTime, restart, cores);

            string dirname = "trash";
            string linpath = @"/work/scratch/jw52xeqa/" + dirname;
            string winpath = @"W:\work\scratch\jw52xeqa\" + dirname;
            c.AlternateDbPaths = new[]{
                    new ValueTuple<string,string>(linpath, ""),
                    new ValueTuple<string,string>(winpath, "pcmit32")
            };
            c.savetodb = true;
            c.saveperiod = savePeriod;
            c.PrintInterval = 1;

            c.ProjectName = "dmr_cube_run3";
            c.NoOfTimesteps = timeSteps;

            return c;
        }

        public static CNSControl ShockTube_HilbertTest(int GPType, int dgDegree, int ExplOrder, int RecInt, int numOfCellsX, int numOfCellsY, bool LTSON, int AVratio, int Tsteps = int.MaxValue, string prjname = null, double sensorLimit = 1e-4, bool true1D = false, bool saveToDb = true, string SessionID = null, string GridID = null) {

            CNSControl c = new CNSControl();
            //c.DbPath = @"D:\Weber\BoSSS\test_db";
            //dbPath = @"e:\bosss_db\GridOfTomorrow\";
            c.DbPath = @"\\fdyprime\userspace\weber\db_Speedup_2";
            //c.DbPath = dbPath;
            //c.savetodb = dbPath != null && saveToDb;
            //c.DbPath = @"/home/yp19ysog/BoSSS_DB";

            c.savetodb = true;
            c.saveperiod = RecInt;
            c.PrintInterval = 1;

            bool AV = false;
            if (dgDegree > 0) {
                AV = true;
            }

            if (LTSON) {
                c.DynamicLoadBalancing_On = true;
                //LTS-Timestepping
                c.ExplicitScheme = ExplicitSchemes.LTS;
                c.ExplicitOrder = ExplOrder;
                c.NumberOfSubGrids = 2;//3
                c.ReclusteringInterval = RecInt;
                c.FluxCorrection = false;

                // Add one balance constraint for each subgrid
                c.DynamicLoadBalancing_CellCostEstimators.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
                c.DynamicLoadBalancing_ImbalanceThreshold = 0.0;
                c.DynamicLoadBalancing_Period = RecInt;
            } else if (!LTSON && AV && AVratio != 0) {
                c.DynamicLoadBalancing_On = true;
                //AV-Loadbalancing
                c.ExplicitScheme = ExplicitSchemes.RungeKutta;
                c.ExplicitOrder = ExplOrder;
                c.NumberOfSubGrids = 2;
                c.ReclusteringInterval = RecInt;
                c.FluxCorrection = false;

                // Add one balance constraint for each subgrid
                if (AVratio > 1) {
                    //direct cost mapping (one map with two values: 1 and AVratio)
                    c.DynamicLoadBalancing_CellCostEstimators.AddRange(ArtificialViscosityCellCostEstimator.GetStaticCostBasedEstimator(AVratio));
                } else {
                    //cluster constraint (2 maps: neutral and AV)
                    c.DynamicLoadBalancing_CellCostEstimators.AddRange(ArtificialViscosityCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
                }
                c.DynamicLoadBalancing_ImbalanceThreshold = 0.0;
                c.DynamicLoadBalancing_Period = RecInt;
            } else {
                // Explicit Timestepping
                c.ExplicitScheme = ExplicitSchemes.RungeKutta;
                c.ExplicitOrder = ExplOrder;
            }

            switch (GPType) {
                case 0:
                    c.GridPartType = GridPartType.none;
                    break;
                case 1:
                    c.GridPartType = GridPartType.ParMETIS;
                    break;
                case 2:
                    if (AVratio > 1) {
                        c.GridPartType = GridPartType.Hilbert;
                        Console.WriteLine("Hilbert is executed");
                    } else {
                        c.GridPartType = GridPartType.clusterHilbert;
                        Console.WriteLine("clusterHilbert with Clusters is executed");
                    }
                    break;
            }

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            double epsilon0 = 1.0;
            double kappa = 0.5;
            if (AV) {
                c.ActiveOperators |= Operators.ArtificialViscosity;
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            }

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            if (!true1D) {
                c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
                c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);
                if (true1D) {
                    var grid = Grid1D.LineGrid(xNodes, periodic: false);
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });
                    return grid;
                } else {
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });
                    return grid;
                }
            };

            c.AddBoundaryValue("AdiabaticSlipWall");

            #region Smoothing of initial condition
            // Normal vector of initial shock
            Vector normalVector = new Vector(1, 0);

            // Direction vector of initial shock
            Vector r = new Vector(normalVector.y, -normalVector.x);
            r.NormalizeInPlace();

            // Distance from a point X to the initial shock
            double[] point = new double[] { 0.5, 0.0 };
            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };
            #endregion

            Func<double, double> Jump = (x => x <= 0.5 ? 0 : 1);

            // Initial conditions
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;

            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            if (true1D == false) {
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1.0e-3;
            c.CFLFraction = 0.1;
            c.Endtime = 0.25;
            c.NoOfTimesteps = Tsteps;

            string LoadbalancingType = "None";
            if (LTSON) {
                LoadbalancingType = "LTS-Cluster";
            } else if (AV && AVratio == 1) {
                LoadbalancingType = "AV-Cluster";
            } else if (AV && AVratio > 1) {
                LoadbalancingType = "AV-direct";
            }
            if (prjname == null) {
                c.ProjectName = LoadbalancingType;
            } else {
                c.ProjectName = prjname;
            }
            c.SessionName = String.Format("ST, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, GridPartType {3}, LoadbalancingType {4}", dgDegree, numOfCellsX, numOfCellsY, c.GridPartType, LoadbalancingType);

            return c;
        }

        public static CNSControl DoubleMachReflection_HilbertTest(int GPType, int dgDegree, int ExplOrder, int RecInt, int numOfCellsX, int numOfCellsY, bool LTSON, int AVratio, int Tsteps = int.MaxValue, string prjname = null, double sensorLimit = 1e-3, bool restart = false, string sessionID = null, string gridID = null) {
            CNSControl c = new CNSControl();

            //c.DbPath = @"D:\Weber\BoSSS\test_db";
            //c.DbPath = @"/home/yp19ysog/BoSSS_DB";
            c.DbPath = @"\\fdyprime\userspace\weber\db_Speedup_2";
            c.savetodb = true;
            c.saveperiod = RecInt;
            c.PrintInterval = 1;

            bool AV = false;
            if (dgDegree > 0) {
                AV = true;
            }

            if (LTSON) {
                c.DynamicLoadBalancing_On = true;
                //Load is balanced according to LTS
                c.ExplicitScheme = ExplicitSchemes.LTS;
                c.ExplicitOrder = ExplOrder;
                c.NumberOfSubGrids = 3;
                c.ReclusteringInterval = RecInt;
                c.FluxCorrection = false;

                // Add one balance constraint for each subgrid
                c.DynamicLoadBalancing_CellCostEstimators.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
                c.DynamicLoadBalancing_ImbalanceThreshold = 0.0;
                c.DynamicLoadBalancing_Period = RecInt;

            } else if (!LTSON && AV && AVratio != 0) {
                c.DynamicLoadBalancing_On = true;
                //Load is balanced according to AV
                c.ExplicitScheme = ExplicitSchemes.RungeKutta;
                c.ExplicitOrder = ExplOrder;
                c.NumberOfSubGrids = 2;
                c.ReclusteringInterval = RecInt;
                c.FluxCorrection = false;

                //AV constraintmap is added
                if (AVratio > 1) {
                    //direct cost mapping (one map with two values: 1 and AVratio)
                    c.DynamicLoadBalancing_CellCostEstimators.AddRange(ArtificialViscosityCellCostEstimator.GetStaticCostBasedEstimator(AVratio));
                } else {
                    //cluster constraint (2 maps: neutral and AV)
                    c.DynamicLoadBalancing_CellCostEstimators.AddRange(ArtificialViscosityCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
                }
                c.DynamicLoadBalancing_ImbalanceThreshold = 0.0;
                c.DynamicLoadBalancing_Period = RecInt;
            } else {
                //Explicit Timestepping
                c.ExplicitScheme = ExplicitSchemes.RungeKutta;
                c.ExplicitOrder = ExplOrder;
            }

            switch (GPType) {
                case 0:
                    c.GridPartType = GridPartType.none;
                    break;
                case 1:
                    c.GridPartType = GridPartType.ParMETIS;
                    break;
                case 2:
                    if (AVratio > 1) {
                        c.GridPartType = GridPartType.Hilbert;
                        Console.WriteLine("Hilbert is executed");
                    } else {
                        c.GridPartType = GridPartType.clusterHilbert;
                        Console.WriteLine("clusterHilbert with Clusters is executed");
                    }
                    break;
            }

            double xMin = 0;
            double xMax = 4;
            double yMin = 0;
            double yMax = 1;

            // Start of the bottom wall, x = 1/6 = 0.166666, (Woodward and Colella 1984)
            // Practical choice: Should be on a cell boundary, because the boundary condition changes from
            // supersonic inflow to adiabatic wall
            double xWall = 0.16;

            double cellSize = Math.Min((xMax - xMin) / numOfCellsX, (yMax - yMin) / numOfCellsY);

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
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);

            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.Viscosity, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            if (dgDegree >= 1) {
                c.AddVariable(CNSVariables.Schlieren, dgDegree - 1);
            }
            if (AV) {
                c.AddVariable(CNSVariables.ShockSensor, 0);
                c.AddVariable(CNSVariables.ArtificialViscosity, 2);
            }

            // Time stepping variables
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }
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
                        if (Math.Abs(X[1]) < 1e-14) {// unten
                            if (X[0] < xWall) {// unten (links)
                                return 1;
                            } else {// unten (Rest)
                                return 3;
                            }
                        } else if (Math.Abs(X[1] - (yMax - yMin)) < 1e-14) {// oben
                            return 1;
                        } else if (Math.Abs(X[0]) < 1e-14) { // links
                            return 1;
                        } else if (Math.Abs(X[0] - (xMax - xMin)) < 1e-14) {// rechts
                            return 2;
                        } else {
                            throw new System.Exception("bla");
                        }
                    });

                    return grid;
                };
            } else {
                c.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(sessionID), -1);
                c.GridGuid = new Guid(gridID);
            }
            Func<double[], double, double> DistanceToLine = delegate (double[] X, double t) {
                // direction vector
                Vector p1 = new Vector(xWall, 0.0);
                Vector p2 = new Vector(xWall + 1 / Math.Tan(Math.PI / 3), 1.0);
                Vector p = p2 - p1;

                // normal vector
                Vector n = new Vector(p.y, -p.x);
                n.NormalizeInPlace();

                // angle between line and x-axis
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

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 2.0 * cellSize / Math.Max(dgDegree, 1);
                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };
            Func<double, double> Jump = (x => x < 0 ? 0 : 1);

            c.AddBoundaryValue("SupersonicInlet", CompressibleVariables.Density, (X, t) => 8.0 - SmoothJump(DistanceToLine(X, t)) * (8.0 - 1.4));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.xComponent, (X, t) => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToLine(X, t)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Velocity.yComponent, (X, t) => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToLine(X, t)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
            c.AddBoundaryValue("SupersonicInlet", CNSVariables.Pressure, (X, t) => 116.5 - SmoothJump(DistanceToLine(X, t)) * (116.5 - 1.0));
            c.AddBoundaryValue("SupersonicOutlet", CNSVariables.Pressure, (X, t) => 1.0);
            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            if (!restart) {
                c.InitialValues_Evaluators.Add(CompressibleVariables.Density, X => 8.0 - SmoothJump(DistanceToLine(X, 0)) * (8.0 - 1.4));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 8.25 * Math.Sin(Math.PI / 3) - SmoothJump(DistanceToLine(X, 0)) * (8.25 * Math.Sin(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => -8.25 * Math.Cos(Math.PI / 3) - SmoothJump(DistanceToLine(X, 0)) * (-8.25 * Math.Cos(Math.PI / 3) - 0.0));
                c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => 116.5 - SmoothJump(DistanceToLine(X, 0)) * (116.5 - 1.0));
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.Endtime = 0.25;
            //c.dtFixed = 1.0e-6;
            c.CFLFraction = 0.3;
            c.NoOfTimesteps = Tsteps;

            string LoadbalancingType = "None";
            if (LTSON) {
                LoadbalancingType = "LTS-Cluster";
            } else if (AV && AVratio == 1) {
                LoadbalancingType = "AV-Cluster";
            } else if (AV && AVratio > 1) {
                LoadbalancingType = "AV-direct";
            }
            if (prjname == null) {
                c.ProjectName = LoadbalancingType;
            } else {
                c.ProjectName = prjname;
            }
            c.SessionName = String.Format("DMR, dgDegree = {0}, noOfCellsX = {1}, noOfCellsX = {2}, GridPartType {3}, LoadbalancingType {4}", dgDegree, numOfCellsX, numOfCellsY, c.GridPartType, LoadbalancingType);

            return c;
        }

        public static CNSControl ShockTube_PredefinedGrid(string dbPath, int numOfCellsX = 3, int numOfCellsY = 3, int dgDegree = 0, double sensorLimit = 1e-4, bool true1D = false, bool saveToDb = true) {

            CNSControl c = new CNSControl();

            //dbPath = @"D:\Weber\BoSSS\test_db";
            //dbPath = @"e:\bosss_db\GridOfTomorrow\";
            //dbPath = @"\\fdyprime\userspace\geisenhofer\bosss_db\";
            c.DbPath = dbPath;
            c.savetodb = dbPath != null && saveToDb;
            c.saveperiod = int.MaxValue;
            c.PrintInterval = 1;

            // Add one balance constraint for each subgrid
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = 10;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();
            //c.DynamicLoadBalancing_CellClassifier

            // dbg_launch();

            c.GridPartType = GridPartType.Predefined;
            c.GridPartOptions = "hallo";

            bool AV = false;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            // (A)LTS
            //c.ExplicitScheme = ExplicitSchemes.LTS;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            //c.NumberOfSubGrids = 1;
            //c.ReclusteringInterval = 0;
            //c.FluxCorrection = false;

            // Add one balance constraint for each subgrid
            //c.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(LTSCellCostEstimator.Factory(c.NumberOfSubGrids));
            //c.DynamicLoadBalancing_ImbalanceThreshold = 0.1;
            //c.DynamicLoadBalancing_Period = 10;
            //c.DynamicLoadBalancing_CellClassifier = new LTSCellClassifier();

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
                Variable sensorVariable = CompressibleVariables.Density;
                c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
                //c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            }

            // Runge-Kutta schemes
            //c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            //c.ExplicitOrder = 4;

            //Adams-Bashforth
            //c.ExplicitScheme = ExplicitSchemes.AdamsBashforth;
            //c.ExplicitOrder = 3;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            c.AddVariable(CNSVariables.Rank, 0);
            if (true1D == false) {
                c.AddVariable(CompressibleVariables.Momentum.yComponent, dgDegree);
                c.AddVariable(CNSVariables.Velocity.yComponent, dgDegree);
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 2);
                }
            } else {
                if (AV) {
                    c.AddVariable(CNSVariables.ArtificialViscosity, 1);
                }
            }
            c.AddVariable(CNSVariables.CFL, 0);
            c.AddVariable(CNSVariables.CFLConvective, 0);
            if (AV) {
                c.AddVariable(CNSVariables.CFLArtificialViscosity, 0);
            }
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            Func<double[], int> MakeMyPartioning = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                int rank;
                if ((x < 0.3) && (y < 0.3)) {
                    rank = 0;
                } else if ((x < 0.3) && (y < 0.6)) {
                    rank = 1;
                } else if ((x < 0.6) && (y < 0.6)) {
                    rank = 2;
                } else {
                    rank = 3;
                }

                return rank;
            };

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, numOfCellsX + 1);

                if (true1D) {
                    var grid = Grid1D.LineGrid(xNodes, periodic: false);
                    // Boundary conditions
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });
                    return grid;
                } else {
                    double[] yNodes = GenericBlas.Linspace(yMin, yMax, numOfCellsY + 1);
                    var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                    // Boundary conditions
                    grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                    grid.DefineEdgeTags(delegate (double[] _X) {
                        return 1;
                    });

                    //grid.AddPredefinedPartitioning("hallo", new int[] { 0, 0, 1, 1, 1});
                    grid.AddPredefinedPartitioning("hallo", MakeMyPartioning);

                    return grid;
                }
            };

            c.AddBoundaryValue("AdiabaticSlipWall");

            // Initial conditions
            c.InitialValues_Evaluators.Add(CompressibleVariables.Density, delegate (double[] X) {
                double x = X[0];

                if (true1D == false) {
                    double y = X[1];
                }

                if (x <= 0.5) {
                    return 1.0;
                } else {
                    return 0.125;
                }
            });
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, delegate (double[] X) {
                double x = X[0];

                if (true1D == false) {
                    double y = X[1];
                }

                if (x <= 0.5) {
                    return 1.0;
                } else {
                    return 0.1;
                }
            });
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => 0.0);
            if (true1D == false) {
                c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => 0.0);
            }

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1.0e-3;
            c.CFLFraction = 0.3;
            c.Endtime = 0.25;
            c.NoOfTimesteps = 1;

            c.ProjectName = "Shock tube";
            if (true1D) {
                c.SessionName = String.Format("Shock tube, 1D, dgDegree = {0}, noOfCellsX = {1}, sensorLimit = {2:0.00E-00}", dgDegree, numOfCellsX, sensorLimit);
            } else {
                c.SessionName = String.Format("Shock tube, 2D, dgDegree = {0}, noOfCellsX = {1}, noOfCellsY = {2}, sensorLimit = {3:0.00E-00}, CFLFraction = {4:0.00E-00}, ALTS {5}/{6}", dgDegree, numOfCellsX, numOfCellsY, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids);
            }
            return c;
        }
    }
}