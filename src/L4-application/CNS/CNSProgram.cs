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
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.ASCIIExport;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.LoadBalancing;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using BoSSS.Solution.LoadBalancing;

namespace CNS {

    /// <summary>
    /// CNS main program
    /// </summary>
    public class CNSProgram : CNSProgram<CNSControl> {

        /// <summary>
        /// Initializes an instance of <see cref="CNSProgram"/> and starts the
        /// program execution as specified in <see cref="Application"/>
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {

            //Application.InitMPI(args);
            //CNS.Tests.ArtificialViscosity.ArtificialViscosityShockTubeTests.ToroTest1_ALTS1_3();
            //CNS.Tests.IBMTests.IBMIsentropicVortexTest.IBMVortexOneStepGaussAndStokesNoAgglomerationTest();
            //CNS.Tests.BoundaryConditions.EulerBoundaryConditionTest.TestSubsonicInletBoundaryCondition1D();
            //CNS.Tests.DiffusiveFlux.SIPGConsistency.SIPGconsistencyTest(1);
            //DeleteOldPlotFiles();
            //CNS.Tests.Ringleb.RinglebTest.RinglebIdealGasTest();

            //CNS.Tests.IBMTests.IBMALTSTest.IBMALTSTest1_4_pos1();
            //CNS.Tests.IBMTests.IBMCylinderTest.IBMCylinder0th();
            //CNS.Tests.Ringleb.RinglebTest.RinglebIdealGasTest();

            //int numOfCellsX = 8;
            //int numOfCellsY = 16;
            //var cIBMBowShock = ControlExamples_Supersonic.IBMBowShock(
            //    dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_P0_db", savePeriod: 1000, dgDegree: 0, CFLFraction: 0.1, explicitScheme: 1, explicitOrder: 1,
            //    endTime: 32, numOfCellsX: numOfCellsX, numOfCellsY: numOfCellsY);


            //var pCNS = new CNS.Program<CNS.IBM.IBMControl>();
            //pCNS.Init(cIBMBowShock);
            //pCNS.RunSolverMode();

            //Debug.Assert(false, "remove me");
            //return;

            Application<CNSControl>._Main(
                args,
                false,
                () => new CNSProgram());
        }
    }

    /// <summary>
    /// Implementation of <see cref="Application{T}"/> specific for the solution
    /// of the compressible Navier-Stokes equations.
    /// </summary>
    public class CNSProgram<T> : Application<T>, IProgram<T>
        where T : CNSControl, new() {

        /// <summary>
        /// The storage of all current variable values of the current flow
        /// field
        /// </summary>
        public CNSFieldSet WorkingSet {
            [DebuggerNonUserCode]
            get;
            [DebuggerNonUserCode]
            protected set;
        }

        /// <summary>
        /// A map that determines the active species in some point in the
        /// domain.
        /// </summary>
        public ISpeciesMap SpeciesMap {
            [DebuggerNonUserCode]
            get;
            [DebuggerNonUserCode]
            protected set;
        }

        /// <summary>
        /// The employed time stepper
        /// </summary>
        public ITimeStepper TimeStepper {
            [DebuggerNonUserCode]
            get;
            [DebuggerNonUserCode]
            protected set;
        }

        /// <summary>
        /// Creates the system of equations to be solved
        /// </summary>
        protected OperatorFactory operatorFactory;

        /// <summary>
        /// </summary>
        public Operator FullOperator {
            get;
            private set;
        }

        /// <summary>
        /// Logs residuals for all variables in <see cref="WorkingSet"/>
        /// </summary>
        protected IResidualLogger[] residualLoggers;

        /// <summary>
        /// Simulation time after restart (needed for time stepping)
        /// </summary>
        protected double StartTime = 0.0;

        /// <summary>
        /// Simulation timestep restart (needed for time stepping)
        /// </summary>
        public int TimestepNumber {
            get;
            protected set;
        }

        /// <summary>
        /// Standard constructor, <see cref="Application"/>
        /// </summary>
        public CNSProgram()
            : base() {
        }

        /// <summary>
        /// Initializes <see cref="CompressibleEnvironment"/> after the grid has been
        /// loaded regularly via <see cref="Application{T}.CreateOrLoadGrid"/>
        /// </summary>
        /// <returns></returns>
        protected override IGrid CreateOrLoadGrid() {
            using (var ht = new FuncTrace()) {
                IGrid grid = base.CreateOrLoadGrid();
                CompressibleEnvironment.Initialize(grid.SpatialDimension);
                return grid;
            }
        }

        /// <summary>
        /// Initializes the <see cref="WorkingSet"/> and adds all relevant
        /// fields to the list of IO variables
        /// </summary>
        /// <remarks>
        /// Initializes <see cref="CompressibleEnvironment"/> since it is the first
        /// method be called by <see cref="Application{T}._Main"/>.
        /// </remarks>
        protected override void CreateFields() {
            using (var ht = new FuncTrace()) {
                WorkingSet = Control.DomainType.CreateWorkingSet(GridData, Control);
                SpeciesMap = Control.DomainType.CreateSpeciesMap(WorkingSet, Control, GridData);

                CompressibleBoundaryCondMap map = GetBoundaryConditionMap();
                operatorFactory = Control.DomainType.GetOperatorFactory(
                    Control, GridData, map, WorkingSet, SpeciesMap);

                m_IOFields.AddRange(WorkingSet.AllFields);
                m_RegisteredFields.AddRange(WorkingSet.AllFields);
            }
        }

        /// <summary>
        /// Creates the correct equations depending on
        /// <see cref="CNSControl.DomainType"/>. Additionally, it creates
        /// the associated time stepper
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase gridUpdateData) {
            using (var ht = new FuncTrace()) {
                FullOperator = operatorFactory.GetJoinedOperator();

                TimeStepper = Control.ExplicitScheme.Instantiate(
                    Control,
                    operatorFactory,
                    WorkingSet,
                    ParameterMapping,
                    SpeciesMap,
                    this);

                // Resets simulation time after a restart
                TimeStepper.ResetTime(StartTime, TimestepNumber);

                // Update time information (needed for LTS runs)
                TimeStepper.UpdateTimeInfo(new TimeInformation(TimestepNumber, StartTime, -1));

                // Configure residual handling
                if (gridUpdateData == null) {
                    // Do not change these settings upon repartitioning
                    ResLogger.WriteResidualsToTextFile = true;
                    ResLogger.WriteResidualsToConsole = false;
                }
                residualLoggers = Control.ResidualLoggerType.Instantiate(
                    this,
                    Control,
                    FullOperator.ToSpatialOperator(WorkingSet),
                    WorkingSet.ConservativeVariables,
                    ParameterMapping).ToArray();

                WorkingSet.UpdateDerivedVariables(this, SpeciesMap.SubGrid.VolumeMask);
            }
        }

        /// <summary>
        /// Coordinate mapping for parameter variables for evaluations of the
        /// discrete operators. 
        /// </summary>
        public CoordinateMapping ParameterMapping {
            get {
                return (WorkingSet.ParameterFields.IsNullOrEmpty()) ? null : new CoordinateMapping(WorkingSet.ParameterFields);
            }
        }

        /// <summary>
        /// Uses <see cref="TimeStepper"/> to advance the solver on time step.
        /// </summary>
        /// <param name="TimestepNo">
        /// <see cref="Application{T}.RunSolverOneStep"/>
        /// </param>
        /// <param name="phystime">
        /// <see cref="Application{T}.RunSolverOneStep"/>
        /// </param>
        /// <param name="dt">
        /// <see cref="Application{T}.RunSolverOneStep"/>
        /// </param>
        /// <returns>
        /// <see cref="Application{T}.RunSolverOneStep"/>
        /// </returns>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            using (var ht = new FuncTrace()) {
                int printInterval = Control.PrintInterval;
                if (DatabaseDriver.MyRank == 0 && TimestepNo % printInterval == 0) {
#if DEBUG
                    Console.WriteLine();
#endif
                    Console.Write("Starting time step #" + TimestepNo + "...");
                }

                // Update shock-capturing variables before performing a time step
                // as the time step constraints (could) depend on artificial viscosity.
                // If not doing so, the artificial viscosity values from the previous
                // time step are taken (unless UpdateDerivedVariables has been called by
                // SavetoDatabase which depends on the saveperiod specified in the control file). 
                if (this.Control.ArtificialViscosityLaw != null) {
                    WorkingSet.UpdateShockCapturingVariables(this, SpeciesMap.SubGrid.VolumeMask);
                }

                // Create TimeInformation object in order to make information available for
                // the time stepper
                TimeStepper.UpdateTimeInfo(new TimeInformation(TimestepNo, phystime, dt));

                #region Evaluate operator for testing
                //// Evaluate the operator
                //if (TimestepNo % printInterval == 0)
                //{
                //    CoordinateMapping mapping = new CoordinateMapping(WorkingSet.ConservativeVariables);
                //    double[] OpAffine = new double[mapping.LocalLength];
                //    SpatialOperator spatialOperator = FullOperator.ToSpatialOperator(WorkingSet);
                //    //var ev = spatialOperator.GetEvaluatorEx(WorkingSet.ConservativeVariables, null, mapping);
                //    var ev = spatialOperator.GetEvaluatorEx(WorkingSet.ConservativeVariables, WorkingSet.ParameterFields, mapping);
                //    ev.Evaluate(1.0, 0.0, OpAffine, null);
                //    //OpAffine.SaveToTextFile(String.Format("ResidualVector_{0}.txt", TimestepNo));
                //    Console.WriteLine($"|R0|={OpAffine.MPI_L2Norm()}");
                //}

                //// Sample points
                //int noOfPoints = 1000;
                //double[] nodes = GenericBlas.Linspace(0.0, 1.0, noOfPoints);
                //MultidimensionalArray points = MultidimensionalArray.Create(noOfPoints, 2);
                //for (int i = 0; i < noOfPoints; i++) {
                //    points[i, 0] = nodes[i];
                //    points[i, 1] = 0.5;
                //}

                // FieldEvaluation
                //MultidimensionalArray results = MultidimensionalArray.Create(noOfPoints, Resi.Length);
                //for (int i = 0; i < Residuals.Length; i++) {
                //    FieldEvaluation fieldEvaluator = new FieldEvaluation((GridData)this.GridData);
                //    fieldEvaluator.Evaluate(1.0, Residuals, points, 0.0, results);
                //}

                //// StreamWriter
                //using (System.IO.StreamWriter sw = new System.IO.StreamWriter(String.Format("Residuals{0}.txt", timestepNo))) {
                //    //Console.WriteLine("x \t y \t result");
                //    sw.WriteLine("x \t y \t rho \t xMom \t yMom \t rhoE");
                //    string resultLine;
                //    for (int i = 0; i < noOfPoints; i++) {
                //        resultLine = points[i, 0] + "\t" + points[i, 1] + "\t" + results[i, 0] + "\t" + results[i, 1] + "\t" + results[i, 2] + "\t" + results[i, 3] + "\t";
                //        //Console.WriteLine(resultLine);
                //        sw.WriteLine(resultLine);
                //    }
                //    sw.Flush();
                //}
                #endregion

                using (new BlockTrace("TimeStepper.Perform", ht)) {
                    dt = TimeStepper.Perform(dt);
                    //} catch (Exception ee) {
                    //    e = ee;
                    //}
                    //e.ExceptionBcast();

                    if (DatabaseDriver.MyRank == 0 && TimestepNo % printInterval == 0) {
                        if (TimestepNo % printInterval == 0) {
                            Console.WriteLine(" done. PhysTime: {0:0.#######E-00}, dt: {1:0.#######E-00}", phystime, dt);                        }
                    }

                    IDictionary<string, double> residuals = residualLoggers.LogTimeStep(TimestepNo, dt, phystime);
                    base.TerminationKey = residualLoggers.ShouldTerminate(residuals, Control);

                    this.ResLogger.NextTimestep(
                        residualLoggers.ShouldLog(TimestepNo, Control.ResidualInterval));
                }

                if (Control.WriteLTSLog && TimeStepper is AdamsBashforthLTS) {
                    this.WriteLTSLog(dt);
                }

                //WorkingSet.ConservativeVariables[0].CoordinateVector.SaveToTextFile(String.Format("DensityCoordVec_{0}.txt", TimestepNo));
                //for (int d = 0; d < 2; d++) {
                //    WorkingSet.ConservativeVariables[d + 1].CoordinateVector.SaveToTextFile(String.Format("Momentum[{0}]CoordVec_{1}.txt", d, TimestepNo));
                //}
                //WorkingSet.ConservativeVariables[3].CoordinateVector.SaveToTextFile(String.Format("EnergyCoordVec_{0}.txt", TimestepNo));

                return dt;
            }
        }

        /// <summary>
        /// Makes sure all derived variables are updated before saving
        /// </summary>
        /// <param name="timestepno"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        protected override ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            using (var ht = new FuncTrace()) {
                WorkingSet.UpdateDerivedVariables(this, SpeciesMap.SubGrid.VolumeMask);
                return base.SaveToDatabase(timestepno, t);
            }
        }

        /// <summary>
        /// See <see cref="SaveToDatabase(TimestepNumber, double)"/>
        /// </summary>
        /// <param name="ts"></param>
        /// <param name="phystime"></param>
        void IProgram<T>.SaveToDatabase(TimestepNumber ts, double phystime) {
            this.SaveToDatabase(ts, phystime);
        }

        private PlotDriver plotDriver;

        /// <summary>
        /// Plots the current state using Tecplot (if
        /// <see cref="CompressibleEnvironment.NumberOfDimensions"/> is greater than 1)
        /// </summary>
        /// <param name="physTime">The physical (simulation time)</param>
        /// <param name="timestepNo">The time step</param>
        /// <param name="superSampling"></param>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            using (var ht = new FuncTrace()) {
                if (plotDriver == null) {
                    if (CompressibleEnvironment.NumberOfDimensions == 1) {
                        plotDriver = new CurveExportDriver(GridData, true, (uint)superSampling);
                    } else {
                        plotDriver = new Tecplot(GridData, true, false, (uint)superSampling);
                    }
                }

                WorkingSet.UpdateDerivedVariables(this, SpeciesMap.SubGrid.VolumeMask);
                plotDriver.PlotFields("CNS-" + timestepNo, physTime, m_IOFields);

                #region Write residuals to text file
                // Sample points
                //int noOfPoints = 16;
                //double[] nodes = GenericBlas.Linspace(-1.8, -0.3, noOfPoints);
                //MultidimensionalArray points = MultidimensionalArray.Create(noOfPoints, 2);
                //for (int i = 0; i < noOfPoints; i++) {
                //    points[i, 0] = nodes[i];
                //    points[i, 1] = -0.65;
                //}

                ////// FieldEvaluation
                ////MultidimensionalArray results = MultidimensionalArray.Create(noOfPoints, Residuals.Length);
                ////for (int i = 0; i < Residuals.Length; i++) {
                ////    FieldEvaluation fieldEvaluator = new FieldEvaluation((GridData)this.GridData);
                ////    fieldEvaluator.Evaluate(1.0, Residuals, points, 0.0, results);
                ////}

                ////// StreamWriter
                ////using (System.IO.StreamWriter sw = new System.IO.StreamWriter(String.Format("Residuals{0}.txt", timestepNo))) {
                ////    //Console.WriteLine("x \t y \t result");
                ////    sw.WriteLine("x \t y \t rho \t xMom \t yMom \t rhoE");
                ////    string resultLine;
                ////    for (int i = 0; i < noOfPoints; i++) {
                ////        resultLine = points[i, 0] + "\t" + points[i, 1] + "\t" + results[i, 0] + "\t" + results[i, 1] + "\t" + results[i, 2] + "\t" + results[i, 3] + "\t";
                ////        //Console.WriteLine(resultLine);
                ////        sw.WriteLine(resultLine);
                ////    }
                ////    sw.Flush();
                ////}
                #endregion

                #region Write DG fields to text file
                // Sample points
                //int noOfPoints = 16;
                //double[] nodes = GenericBlas.Linspace(-1.8, -0.3, noOfPoints);
                //MultidimensionalArray points = MultidimensionalArray.Create(noOfPoints, 2);
                //for (int i = 0; i < noOfPoints; i++) {
                //    points[i, 0] = nodes[i];
                //    points[i, 1] = -0.65;
                //}

                //// FieldEvaluation
                //MultidimensionalArray resultsFields = MultidimensionalArray.Create(noOfPoints, m_IOFields.Count());
                //for (int i = 0; i < m_IOFields.Count(); i++) {
                //    FieldEvaluation fieldEvaluator = new FieldEvaluation((GridData)this.GridData);
                //    fieldEvaluator.Evaluate(1.0, m_IOFields, points, 0.0, resultsFields);
                //}

                //// StreamWriter
                //using (System.IO.StreamWriter sw = new System.IO.StreamWriter(String.Format("DGFields{0}.txt", timestepNo))) {
                //    //Console.WriteLine("x \t y \t result");
                //    //sw.WriteLine("x \t y \t rho \t xMom \t yMom \t rhoE");
                //    string resultLine;
                //    for (int i = 0; i < noOfPoints; i++) {
                //        resultLine = points[i, 0] + "\t" + points[i, 1] + "\t" + resultsFields[i, 0] + "\t" + resultsFields[i, 1] + "\t" + resultsFields[i, 2] + "\t" + resultsFields[i, 3] + "\t";
                //        //Console.WriteLine(resultLine);
                //        sw.WriteLine(resultLine);
                //    }
                //    sw.Flush();
                //}
                //WorkingSet.ConservativeVariables[0].CoordinateVector.SaveToTextFile(String.Format("DensityCoordVec_{0}.txt", timestepNo));
                //WorkingSet.ConservativeVariables[3].CoordinateVector.SaveToTextFile(String.Format("EnergyCoordVec_{0}.txt", timestepNo));
                #endregion
            }
        }

        /// <summary>
        /// Sets the initial conditions
        /// </summary>
        protected override void SetInitial(double t) {
            using (var ht = new FuncTrace()) {
                WorkingSet.ProjectInitialValues(SpeciesMap, base.Control.InitialValues_Evaluators);
            }
        }

        private int firstTimeStepToLog;

        /// <summary>
        /// Sets the simulation time of the restart (needed for time stepping)
        /// and recomputes all derived variables
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {
            using (var ht = new FuncTrace()) {
                this.StartTime = time;
                this.TimestepNumber = timestep.MajorNumber;

                if (SpeciesMap is ImmersedSpeciesMap ibmMap) {
                    LsTrk = ibmMap.Tracker;
                }

                if (this.Control.WriteLTSLog) {
                    this.firstTimeStepToLog = timestep.MajorNumber + Control.ExplicitOrder - 1;
                }
            }
        }

        /// <summary>
        /// See <see cref="Application{T}.ComputeNewCellDistribution(int, double)"/>
        /// </summary>
        /// <param name="TimeStepNo"></param>
        /// <param name="physTime"></param>
        /// <returns></returns>
        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            using (var ht = new FuncTrace()) {
                if (SpeciesMap is ImmersedSpeciesMap ibmMap) {
                    LsTrk = ibmMap.Tracker;
                }

                return base.ComputeNewCellDistribution(TimeStepNo, physTime);
            }
        }

        /*
        /// <summary>
        /// See <see cref="ICellClassifier"/>
        /// </summary>
        /// <param name="NoOfClasses"></param>
        /// <param name="cellToPerformanceClassMap"></param>
        /// <param name="TimeStepNo"></param>
        /// <param name="physTime"></param>
        protected override void GetCellPerformanceClasses(out int NoOfClasses, out int[] cellToPerformanceClassMap, int TimeStepNo, double physTime) {
            using (var ht = new FuncTrace()) {
                if (Control.DynamicLoadBalancing_CellClassifier is ArtificialViscosityCellClassifier || TimeStepper is AdamsBashforthLTS) {
                    // Just to be sure...
                    if (this.Control.ArtificialViscosityLaw != null) {
                        WorkingSet.UpdateShockCapturingVariables(this, SpeciesMap.SubGrid.VolumeMask);
                    }
                }

                // Update clustering before cell redistribution when LTS is being used
                if (TimeStepper is AdamsBashforthLTS ABLTSTimeStepper) {
                    if (TimeStepNo % Control.DynamicLoadBalancing_Period != 0 && !(Control.DynamicLoadBalancing_RedistributeAtStartup && TimeStepNo == TimeStepNoRestart)) {
                        throw new Exception("Mismatch between time step number and dynamic load balacing period!");
                    }

                    // Just to be sure...
                    //if (this.Control.ArtificialViscosityLaw != null) {
                    //    WorkingSet.UpdateShockCapturingVariables(this, SpeciesMap.SubGrid.VolumeMask);
                    //}

                    ABLTSTimeStepper.UpdateTimeInfo(new TimeInformation(TimeStepNo, physTime, -1));
                    // LoadBal and noLoadBalRuns did not match, with this fix, it works --> probably some ABevolver were not updated correctly
                    bool reclustered = ABLTSTimeStepper.TryNewClustering(dt: -1, calledByMPIRedist: true);
                    //Debug.Assert(reclustered == true);
                    //ABLTSTimeStepper.SetReclusteredByGridRedist(true);
                }

                (NoOfClasses, cellToPerformanceClassMap) = Control.DynamicLoadBalancing_CellClassifier.ClassifyCells(this);
            }
        }
        */


        /// <summary>
        /// Constructs the boundary condition map to be used by the solver.
        /// Override this method if a specific application support specialized,
        /// custom boundary conditions.
        /// </summary>
        /// <returns></returns>
        protected virtual CompressibleBoundaryCondMap GetBoundaryConditionMap() {
            return new CompressibleBoundaryCondMap(this.GridData, Control, Control.GetMaterial());
        }

        /// <summary>
        /// saves interface points
        /// </summary>
        TextWriter LTSLogWriter;

        /// <summary>
        /// Initializes the format of the log file
        /// </summary>
        /// <param name="sessionID"></param>
        public void InitLTSLogFile(Guid sessionID) {
            // Init text writer
            if (sessionID.ToString() == "00000000-0000-0000-0000-000000000000") {
                // When not using the BoSSS database, write log to directory where the executable is stores
                LTSLogWriter = new StreamWriter("LTSLog.txt");
            } else {
                // When using the BoSSS database, write log to session directory
                LTSLogWriter = base.DatabaseDriver.FsDriver.GetNewLog("LTSLog", sessionID);
            }

            // Header
            LTSLogWriter.Write(
                "LOCAL TIME STEPPING LOG FILE" +
                "\n------------------------------------------\n" +
                "ExplicitScheme = " + Control.ExplicitScheme.ToString() +
                "\nExplicitOrder = " + Control.ExplicitOrder +
                "\nNumberOfSubGrids = " + Control.NumberOfSubGrids +
                "\nReclusteringInterval = " + Control.ReclusteringInterval +
                "\nMaxNumberOfSubSteps = " + Control.maxNumOfSubSteps +
                "\n------------------------------------------\n");
            string titleForColumns = String.Format("{0}\t{1}\t{2}", "ts", "physTime", "dt");
            AdamsBashforthLTS LTS = TimeStepper as AdamsBashforthLTS;
            for (int i = 0; i < LTS.NumberOfClustersInitial; i++) {
                titleForColumns = titleForColumns + String.Format("\t{0}\t{1}\t{2}", "c" + i + "_dt", "c" + i + "_substeps", "c" + i + "_elements");
            }
            LTSLogWriter.WriteLine(titleForColumns);
            LTSLogWriter.Flush();
        }

        /// <summary>
        /// Writes a line to the log file 
        /// initialized by <see cref="InitLTSLogFile(Guid)"/>.
        /// </summary>
        public void WriteLTSLog(double dt) {
            using (new FuncTrace()) {
                if (MPIRank != 0) {
                    return;
                }

                // Init if necessary
                if (LTSLogWriter == null) {
                    InitLTSLogFile(this.CurrentSessionInfo.ID);
                }

                // Start logging when start-up phase of LTS time stepper has been finished
                // There is also a variant when the simulation has been restarted
                bool logIt;
                if (firstTimeStepToLog > 0) {
                    logIt = TimeStepper.TimeInfo.TimeStepNumber > firstTimeStepToLog;
                } else {
                    logIt = TimeStepper.TimeInfo.TimeStepNumber > Control.ExplicitOrder - 1;
                }

                // Write a line
                if (logIt) {
                    AdamsBashforthLTS LTS = TimeStepper as AdamsBashforthLTS;
                    string line = String.Format("{0}\t{1}\t{2}", LTS.TimeInfo.TimeStepNumber, LTS.TimeInfo.PhysicalTime, dt);
                    for (int i = 0; i < LTS.NumberOfClustersInitial; i++) {
                        if (i >= LTS.CurrentClustering.NumberOfClusters) {
                            line = line + String.Format("\t0\t0\t0");   // Add zeroes if current clustering has less clusters than specified at the beginning
                        } else {
                            line = line + String.Format("\t{0}\t{1}\t{2}", LTS.log_clusterDts[i], LTS.log_clusterSubSteps[i], LTS.log_clusterElements[i]);
                        }
                    }
                    LTSLogWriter.WriteLine(line);
                    LTSLogWriter.Flush();
                }
            }
        }
    }
}
