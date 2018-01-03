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
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using ilPSP.Tracing;
using log4net;
using NSE_SIMPLE.LowMach;
using NSE_SIMPLE.Multiphase;
using System;
using BoSSS.Solution.Utils;

namespace NSE_SIMPLE {

    /// <summary>
    /// Main class of the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) Navier-Stokes-Solver
    /// </summary>
    public partial class NSE_SIMPLEMain : BoSSS.Solution.Application<SIMPLEControl> {

        static void Main(string[] args) {
            _Main(args, false, null, delegate () {
                NSE_SIMPLEMain p = new NSE_SIMPLEMain();
                return p;
            });
        }

        /// <summary>
        /// Physical variables.
        /// </summary>
        public VariableSet WorkingSet {
            get;
            private set;
        }

        VariableMatrices m_WorkingSetMatrices = null;

        /// <summary>
        /// <see cref="VariableMatrices"/> for SIMPLE-algorithm.
        /// </summary>
        public VariableMatrices WorkingSetMatrices {
            get {
                return m_WorkingSetMatrices;
            }
        }

        /// <summary>
        /// Solver configurations read from external properties in control-file.
        /// </summary>
        public SolverConfiguration SolverConf {
            get;
            private set;
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void CreateFields() {
            using (new FuncTrace()) {

                // create fields for flow field               
                WorkingSet = new VariableSet(base.GridData, Control, base.m_IOFields, base.m_RegisteredFields);

                // create extended fields for variable density solvers
                WorkingSet.CreateExtendedVariables(Control);

                // read settings for SIMPLE-algorithm
                SolverConf = new SolverConfiguration(
                    Control, GridData, WorkingSet, base.MPIRank);

                // Plot analytic solution
                PlotAnalyticSolution(Control.AnalyticVelocityX, WorkingSet.VelBasis, "VelocityX_Ana");
                PlotAnalyticSolution(Control.AnalyticVelocityY, WorkingSet.VelBasis, "VelocityY_Ana");
                if (SolverConf.SpatialDimension == 3)
                    PlotAnalyticSolution(Control.AnalyticVelocityZ, WorkingSet.VelBasis, "VelocityZ_Ana");

                PlotAnalyticSolution(Control.AnalyticPressure, WorkingSet.PressureBasis, "Pressure_Ana");

                switch (Control.PhysicsMode) {
                    case PhysicsMode.Incompressible:
                        break;
                    case PhysicsMode.LowMach:
                        LowMachSIMPLEControl lowMachConf = Control as LowMachSIMPLEControl;
                        PlotAnalyticSolution(lowMachConf.AnalyticTemperature, WorkingSet.TemperatureBasis, "Temperature_Ana");
                        PlotAnalyticSolution(lowMachConf.AnalyticDensity, WorkingSet.TemperatureBasis, "Density_Ana");
                        break;
                    case PhysicsMode.Multiphase:
                        MultiphaseSIMPLEControl multiphaseConf = Control as MultiphaseSIMPLEControl;
                        PlotAnalyticSolution(multiphaseConf.AnalyticLevelSet, WorkingSet.LevelSetBasis, "LevelSet_Ana");
                        PlotAnalyticSolution(multiphaseConf.AnalyticDensity, WorkingSet.LevelSetBasis, "Density_Ana");
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        void PlotAnalyticSolution(Func<double[], double> AnalyticSolution, Basis PlotBasis, string Identification) {
            if (AnalyticSolution != null) {
                SinglePhaseField tmp = new SinglePhaseField(PlotBasis, Identification);
                tmp.ProjectField(AnalyticSolution);
                m_IOFields.Add(tmp);
            }
        }

        /// <summary>
        /// Initializes velocity and pressure
        /// </summary>
        protected override void SetInitial() {

            //WriteQuadNodesOrrSommerfeld();
            //InitOrrSommerfeld();

            base.SetInitial();

            //TaylorVortexHack();

            WorkingSet.Initialize(Control);

            // Create WorkingSetMatrices
            switch (Control.PhysicsMode) {
                case PhysicsMode.Incompressible:
                    break;
                case PhysicsMode.LowMach:
                    LowMachSIMPLEControl lowMachConf = Control as LowMachSIMPLEControl;
                    m_WorkingSetMatrices = new VariableMatrices(base.GridData, WorkingSet.VelBasis, WorkingSet.PressureBasis,
                        lowMachConf.EoS, WorkingSet.Temperature.Current);
                    break;
                case PhysicsMode.Multiphase:
                    MultiphaseSIMPLEControl multiphaseConf = Control as MultiphaseSIMPLEControl;
                    m_WorkingSetMatrices = new VariableMatrices(base.GridData, WorkingSet.VelBasis, WorkingSet.PressureBasis,
                        multiphaseConf.EoS, WorkingSet.Phi.Current);
                    break;
                default:
                    throw new NotImplementedException();
            }

            WorkingSet.CheckForNanOrInf(Control);
            // push start values to history
            WorkingSet.Push(Control);
        }

        /// <summary>
        /// ...
        /// </summary>
        public override void PostRestart(double time, TimestepNumber timestep) {

            //InitLogEnergyOrrSommerfeld();

            WorkingSet.Initialize(Control);

            // Create WorkingSetMatrices
            switch (Control.PhysicsMode) {
                case PhysicsMode.Incompressible:
                    break;
                case PhysicsMode.LowMach:
                    LowMachSIMPLEControl lowMachConf = Control as LowMachSIMPLEControl;
                    m_WorkingSetMatrices = new VariableMatrices(base.GridData, WorkingSet.VelBasis, WorkingSet.PressureBasis,
                        lowMachConf.EoS, WorkingSet.Temperature.Current);
                    break;
                case PhysicsMode.Multiphase:
                    MultiphaseSIMPLEControl multiphaseConf = Control as MultiphaseSIMPLEControl;
                    m_WorkingSetMatrices = new VariableMatrices(base.GridData, WorkingSet.VelBasis, WorkingSet.PressureBasis,
                        multiphaseConf.EoS, WorkingSet.Phi.Current);
                    break;
                default:
                    throw new NotImplementedException();
            }

            WorkingSet.CheckForNanOrInf(Control);
            // push start values to history
            WorkingSet.Push(Control);
        }

        /// <summary>
        /// created during <see cref="CreateEquationsAndSolvers"/>
        /// </summary>
        public SIMPLEStepStatus SIMPLEStatus;

        /// <summary>
        /// created during <see cref="CreateEquationsAndSolvers"/>
        /// </summary>
        ISIMPLEStep SIMPLEStep;

        /// <summary>
        /// 
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            using (new FuncTrace()) {

                // Create SIMPLEStatus
                SIMPLEStatus = new SIMPLEStepStatus(Control);

                // Create SIMPLEStep
                switch (Control.PhysicsMode) {
                    case PhysicsMode.Incompressible:
                        SIMPLEStep = new SIMPLEStepIncompressible(SolverConf, WorkingSet);
                        break;
                    case PhysicsMode.LowMach:
                        SIMPLEStep = new SIMPLEStepLowMach(SolverConf, WorkingSet, WorkingSetMatrices);
                        break;
                    case PhysicsMode.Multiphase:
                        SIMPLEStep = new SIMPLEStepMultiphase(SolverConf, WorkingSet, WorkingSetMatrices);
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        ILog m_Logger = LogManager.GetLogger(typeof(NSE_SIMPLEMain));

        /// <summary>
        /// 
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (var tr = new FuncTrace()) {
                tr.Info("Performing time iteration No. " + TimestepNo);

                // Set dt and SIMPLEStatus
                switch (Control.Algorithm) {
                    case SolutionAlgorithms.Steady_SIMPLE: {
                            dt = 0.0;
                            break;
                        }
                    case SolutionAlgorithms.Unsteady_SIMPLE: {
                            dt = SolverConf.dt;
                            SIMPLEStatus.NextTimestep();
                            break;
                        }
                }

                // some console-output
                if (base.MPIRank == 0) {
                    switch (Control.Algorithm) {
                        case SolutionAlgorithms.Steady_SIMPLE:
                            Console.WriteLine("Starting steady calculation...\n");
                            Console.WriteLine("Starting SIMPLE-Iterations...\n");
                            break;
                        case SolutionAlgorithms.Unsteady_SIMPLE:
                            Console.WriteLine("Starting time step #" + TimestepNo + "...\n");
                            Console.WriteLine("Starting SIMPLE-Iterations...\n");
                            break;
                        default:
                            throw new NotImplementedException();
                    }
                }

                do {
                    // do one SIMPLE iteration
                    SIMPLEStatus.NextSIMPLEIteration();
                    SIMPLEStep.OverallIteration(ref SIMPLEStatus, dt, ResLogger);

                    TerminationKey = WorkingSet.CheckForNanOrInf(Control);
                    if (TerminationKey) {
                        m_Logger.Warn("Found Nan in some field.");
                        if (base.MPIRank == 0)
                            Console.WriteLine("ERROR: Found Nan in some field.");
                        Console.ReadKey();
                    }

                    if ((Control.PhysicsMode == PhysicsMode.LowMach) && (SolverConf.Control.As<LowMachSIMPLEControl>().EdgeTagsNusselt != null))
                        CalculateNusselt(SIMPLEStatus.Timestep, base.GridData, WorkingSet.Temperature.Current, Control);

                    // save to database
                    if (SIMPLEStatus.SaveStep) {
                        SaveToDatabase(SIMPLEStatus.Timestep, phystime);
                    }

                    // calculate errors
                    int QuadDegreePressure = 20;
                    int QuadDegreeVel = 20;

                    ResLogger.ComputeL2Error(WorkingSet.Pressure, Control.AnalyticPressure, QuadDegreePressure, "p_ana");

                    ResLogger.ComputeL2Error(WorkingSet.Velocity.Current[0], Control.AnalyticVelocityX, QuadDegreeVel, "u_ana");
                    ResLogger.ComputeL2Error(WorkingSet.Velocity.Current[1], Control.AnalyticVelocityY, QuadDegreeVel, "v_ana");
                    if (SolverConf.SpatialDimension == 3)
                        ResLogger.ComputeL2Error(WorkingSet.Velocity.Current[2], Control.AnalyticVelocityZ, QuadDegreeVel, "w_ana");

                    switch (Control.PhysicsMode) {
                        case PhysicsMode.Incompressible:
                            break;
                        case PhysicsMode.LowMach:
                            LowMachSIMPLEControl lowMachConf = Control as LowMachSIMPLEControl;
                            ResLogger.ComputeL2Error(WorkingSet.Temperature.Current, lowMachConf.AnalyticTemperature, QuadDegreeVel, "T_ana");
                            ResLogger.ComputeL2Error(WorkingSet.Rho, lowMachConf.AnalyticDensity, QuadDegreeVel, "Rho_ana");
                            break;
                        case PhysicsMode.Multiphase:
                            MultiphaseSIMPLEControl multiphaseConf = Control as MultiphaseSIMPLEControl;
                            ResLogger.ComputeL2Error(WorkingSet.Phi.Current, multiphaseConf.AnalyticLevelSet, QuadDegreeVel, "Phi_ana");
                            ResLogger.ComputeL2Error(WorkingSet.Rho, multiphaseConf.AnalyticDensity, QuadDegreeVel, "Rho_ana");
                            break;
                        default:
                            throw new NotImplementedException();
                    }

                    // terminate SIMPLE in case of divergence
                    if (ResLogger.Residuals["L2Norm p'"] > 1.0E+10)
                        TerminationKey = true;

                    // push residual logger to next iteration
                    switch (Control.Algorithm) {
                        case SolutionAlgorithms.Steady_SIMPLE:
                            ResLogger.NextIteration(false);
                            break;
                        case SolutionAlgorithms.Unsteady_SIMPLE:
                            ResLogger.NextIteration(true);
                            break;
                        default:
                            throw new NotImplementedException();
                    }
                }

                while (!SIMPLEStatus.IsConverged && !SIMPLEStatus.TerminateSIMPLE && !TerminationKey);

                // determine cause for end of SIMPLE iterations
                if (SIMPLEStatus.IsConverged) {
                    tr.Info("Solution converged.");
                } else if (SIMPLEStatus.TerminateSIMPLE) {
                    if (SIMPLEStatus.SIMPLEStepNo == Control.MaxNoSIMPLEsteps) {
                        SIMPLEStatus.CntMaxNoSIMPLEsteps++;
                        m_Logger.Warn("MaxNoSIMPLEsteps are reached.");
                    } else {
                        m_Logger.Warn("Unknown reason for terminating SIMPLE iterations - should not happen.");
                    }
                } else {
                    m_Logger.Error("Solution diverged.");
                }

                // save the new timestep
                switch (Control.Algorithm) {
                    case SolutionAlgorithms.Steady_SIMPLE:
                        break;
                    case SolutionAlgorithms.Unsteady_SIMPLE:
                        WorkingSet.Push(Control);
                        ResLogger.NextTimestep(false);
                        break;
                    default:
                        throw new NotImplementedException();
                }

                // some console-output
                if (SIMPLEStatus.IsConverged) {
                    Console.WriteLine("\nINFO: Done SIMPLE-Iterations - Solution converged.\n");
                } else if (SIMPLEStatus.SIMPLEStepNo == Control.MaxNoSIMPLEsteps) {
                    Console.WriteLine("\nWARNING: Done SIMPLE-Iterations - Maximum number of SIMPLE steps was reached.\n");
                } else {
                    Console.WriteLine("\nERROR: Calculation was terminated - Solution diverged.\n");
                }

                switch (Control.Algorithm) {
                    case SolutionAlgorithms.Steady_SIMPLE:
                        Console.WriteLine("Done steady calculation.");
                        break;
                    case SolutionAlgorithms.Unsteady_SIMPLE:
                        Console.WriteLine("Done time step #" + TimestepNo + ".\n");
                        break;
                    default:
                        throw new NotImplementedException();
                }


                //LogEnergyOrrSommerfeld(TimestepNo, phystime, dt);                
                if (Control.EdgeTagsDragAndLift != null)
                    CalculateDragAndLift(phystime);

                //Log temperature history tall cavity
                //LogTemperature(phystime, this.WorkingSet.Temperature.Current);

                return dt;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override void RunSolverMode() {
            using (new FuncTrace()) {

                base.RunSolverMode();

                SIMPLEStep.Dispose();

                //Postprocessing pressure
                //=======================

                if (!double.IsNaN(Control.PressureMeanValue)) {
                    // set specified mean value of pressure
                    double NegCalcMeanValuePressure = -1.0 * WorkingSet.Pressure.GetMeanValueTotal(null);
                    WorkingSet.Pressure.AccConstant(NegCalcMeanValuePressure);
                    WorkingSet.Pressure.AccConstant(Control.PressureMeanValue);
                }

                if ((Control.PressureReferencePoint != null) && double.IsNaN(Control.PressureMeanValue)) {
                    // set pressure to zero at the specified reference point
                    double ShiftPressure = -1.0 * WorkingSet.Pressure.ProbeAt(Control.PressureReferencePoint);
                    WorkingSet.Pressure.AccConstant(ShiftPressure);
                }

                //Some more output
                //================

                if ((Control.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) && !TerminationKey) {
                    if (SIMPLEStatus.CntMaxNoSIMPLEsteps == 0) {
                        m_Logger.Info("All time steps converged.");
                        if (base.MPIRank == 0)
                            Console.WriteLine("INFO: All time steps converged.");
                    } else {
                        m_Logger.Warn("MaxNoSIMPLEsteps was reached in " + SIMPLEStatus.CntMaxNoSIMPLEsteps + " time steps!");
                        if (base.MPIRank == 0)
                            Console.WriteLine("WARNING: MaxNoSIMPLEsteps was reached in " + SIMPLEStatus.CntMaxNoSIMPLEsteps + " time steps!");
                    }
                }

                if (base.MPIRank == 0 && Log_DragAndLift != null)
                    Log_DragAndLift.Close();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="physTime"></param>
        /// <param name="timestepNo"></param>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(m_IOFields, "SIMPLE-TimeStep" + timestepNo, physTime, superSampling);
        }
    }
}