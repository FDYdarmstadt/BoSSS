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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.ASCIIExport;
using BoSSS.Solution.Tecplot;
using CNS.Boundary;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.LoadBalancing;
using CNS.Residual;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CNS {

    /// <summary>
    /// CNS main program
    /// </summary>
    public class Program : Program<CNSControl> {

        /// <summary>
        /// Initializes an instance of <see cref="Program"/> and starts the
        /// program execution as specified in <see cref="Application"/>
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            Application<CNSControl>._Main(
                args,
                false,
                () => new Program());
        }
    }

    /// <summary>
    /// Implementation of <see cref="Application{T}"/> specific for the solution
    /// of the compressible Navier-Stokes equations.
    /// </summary>
    public class Program<T> : Application<T>, IProgram<T>
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
        public Program()
            : base() {
        }

        /// <summary>
        /// Initializes <see cref="CNSEnvironment"/> after the grid has been
        /// loaded regularly via <see cref="Application{T}.CreateOrLoadGrid"/>
        /// </summary>
        /// <returns></returns>
        protected override GridCommons CreateOrLoadGrid() {
            GridCommons grid = base.CreateOrLoadGrid();
            CNSEnvironment.Initialize(grid.SpatialDimension);
            return grid;
        }

        /// <summary>
        /// Initializes the <see cref="WorkingSet"/> and adds all relevant
        /// fields to the list of IO variables
        /// </summary>
        /// <remarks>
        /// Initializes <see cref="CNSEnvironment"/> since it is the first
        /// method be called by <see cref="Application{T}._Main"/>.
        /// </remarks>
        protected override void CreateFields() {
            WorkingSet = Control.DomainType.CreateWorkingSet(GridData, Control);
            SpeciesMap = Control.DomainType.CreateSpeciesMap(WorkingSet, Control, GridData);

            BoundaryConditionMap map = GetBoundaryConditionMap();
            operatorFactory = Control.DomainType.GetOperatorFactory(
                Control, GridData, map, WorkingSet, SpeciesMap);

            m_IOFields.AddRange(WorkingSet.AllFields);
            m_RegisteredFields.AddRange(WorkingSet.AllFields);
        }

        /// <summary>
        /// Creates the correct equations depending on
        /// <see cref="CNSControl.DomainType"/>. Additionally, it creates
        /// the associated time stepper
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase gridUpdateData) {
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

            // Configure residual handling
            if (gridUpdateData == null) {
                // Do not change these settings upon repartitioning
                ResLogger.WriteResidualsToTextFile = true;
                ResLogger.WriteResidualsToConsole = false;
            }
            residualLoggers = Control.ResidualLoggerType.Instantiate(
                this,
                Control,
                FullOperator.ToSpatialOperator(WorkingSet)).ToArray();
            
            WorkingSet.UpdateDerivedVariables(this, SpeciesMap.SubGrid.VolumeMask);
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
                    Console.Write("Starting time step #" + TimestepNo + "...");
                }

                Exception e = null;
                try {
                    dt = TimeStepper.Perform(dt);
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                if (TimestepNo % printInterval == 0) {
                    Console.WriteLine(" done. PhysTime: {0:0.#######E-00}, dt: {1:0.###E-00}", phystime, dt);
                }

                IDictionary<string, double> residuals = residualLoggers.LogTimeStep(TimestepNo, dt, phystime);
                base.TerminationKey = ShouldTerminate(residuals);

                this.ResLogger.NextTimestep(
                    residualLoggers.ShouldLog(TimestepNo, Control.ResidualInterval));

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
            WorkingSet.UpdateDerivedVariables(this, SpeciesMap.SubGrid.VolumeMask);
            return base.SaveToDatabase(timestepno, t);
        }

        /// <summary>
        /// See <see cref="SaveToDatabase(TimestepNumber, double)"/>
        /// </summary>
        /// <param name="ts"></param>
        /// <param name="phystime"></param>
        void IProgram<T>.SaveToDatabase(TimestepNumber ts, double phystime) {
            this.SaveToDatabase(ts, phystime);
        }

        private bool ShouldTerminate(IDictionary<string, double> residuals) {
            if (Control.ResidualBasedTerminationCriteria.Count > 0
                && residuals.Count > 0) {
                bool terminate = true;
                foreach (var keyThresholdPair in Control.ResidualBasedTerminationCriteria) {
                    if (!residuals.ContainsKey(keyThresholdPair.Key)) {
                        throw new Exception(String.Format(
                            "A termination criterion is based on {0} was found"
                                + " but the corresponding residual value was"
                                + " not calculated.",
                            keyThresholdPair.Key));
                    }

                    terminate &= residuals[keyThresholdPair.Key] < keyThresholdPair.Value;
                }

                if (terminate) {
                    Console.WriteLine("All residual criteria fulfilled, stopping calculation.");
                    return true;
                }
            }

            return false;
        }

        private PlotDriver plotDriver;

        /// <summary>
        /// Plots the current state using Tecplot (if
        /// <see cref="CNSEnvironment.NumberOfDimensions"/> is greater than 1)
        /// </summary>
        /// <param name="physTime">The physical (simulation time)</param>
        /// <param name="timestepNo">The time step</param>
        /// <param name="superSampling"></param>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            if (plotDriver == null) {
                if (CNSEnvironment.NumberOfDimensions == 1) {
                    plotDriver = new CurveExportDriver(GridData, true, (uint)superSampling);
                } else {
                    plotDriver = new Tecplot(GridData, true, false, (uint)superSampling);
                }
            }

            plotDriver.PlotFields("CNS-" + timestepNo, physTime, m_IOFields);
        }

        /// <summary>
        /// Sets the initial conditions
        /// </summary>
        protected override void SetInitial() {
            WorkingSet.ProjectInitialValues(SpeciesMap, base.Control.InitialValues_Evaluators);
        }

        /// <summary>
        /// Sets the simulation time of the restart (needed for time stepping)
        /// and recomputes all derived variables
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {
            this.StartTime = time;

            if (SpeciesMap is ImmersedSpeciesMap ibmMap) {
                LsTrk = ibmMap.Tracker;
            }
        }

        /// <summary>
        /// See <see cref="Application{T}.ComputeNewCellDistribution(int, double)"/>
        /// </summary>
        /// <param name="TimeStepNo"></param>
        /// <param name="physTime"></param>
        /// <returns></returns>
        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if (SpeciesMap is ImmersedSpeciesMap ibmMap) {
                LsTrk = ibmMap.Tracker;
            }

            return base.ComputeNewCellDistribution(TimeStepNo, physTime);
        }

        /// <summary>
        /// See <see cref="ICellClassifier"/>
        /// </summary>
        /// <param name="NoOfClasses"></param>
        /// <param name="cellToPerformanceClassMap"></param>
        protected override void GetCellPerformanceClasses(out int NoOfClasses, out int[] cellToPerformanceClassMap) {
            (NoOfClasses, cellToPerformanceClassMap) = Control.DynamicLoadBalancing_CellClassifier.ClassifyCells(this);
        }

        /// <summary>
        /// Constructs the boundary condition map to be used by the solver.
        /// Override this method if a specific application support specialized,
        /// custom boundary conditions.
        /// </summary>
        /// <returns></returns>
        protected virtual BoundaryConditionMap GetBoundaryConditionMap() {
            return new BoundaryConditionMap(GridData, Control);
        }
    }
}
