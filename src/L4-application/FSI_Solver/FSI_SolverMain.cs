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
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using FSI_Solver;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain {

        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {

            //BoSSS.Solution.Application.InitMPI();
            //BoSSS.Application.FSI_Solver.TestProgram.PeriodicTest();
            //BoSSS.Solution.Application.FinalizeMPI();
            //throw new ApplicationException( "remove me");
            
            _Main(args, false, delegate () {
                var p = new FSI_SolverMain();
                return p;
            });
            
        }

        /// <summary>
        /// Set the inital state of the simulation.
        /// </summary>
        protected override void SetInitial() {
            m_Particles = ((FSI_Control)this.Control).Particles;
            if (m_Particles.IsNullOrEmpty())
                throw new Exception("Define at least on particle");
            UpdateLevelSetParticles(phystime: 0.0);
            CreatePhysicalDataLogger();
            CreateResidualLogger();
            base.SetInitial();
        }

        private readonly int spatialDim = 2;

        /// <summary>
        /// A list for all particles
        /// </summary>
        private List<Particle> m_Particles;

        /// <summary>
        /// External access to particle list.
        /// </summary>
        public IList<Particle> Particles => m_Particles;

        /// <summary>
        /// An object with some additional methods
        /// </summary>
        readonly private FSI_Auxillary Auxillary = new FSI_Auxillary();

        /// <summary>
        /// Methods dealing with colouring and level set
        /// </summary>
        private FSI_LevelSetUpdate levelSetUpdate;

        /// <summary>
        /// Particle color field. The level set is only defined on colored cells.
        /// </summary>
        private SinglePhaseField ParticleColor;

        /// <summary>
        /// Level set distance field. 
        /// </summary>
        private SinglePhaseField LevelSetDistance;

        /// <summary>
        /// Field to make the process partition visible in exported visit files
        /// </summary>
        private SinglePhaseField ProcessColor;

        /// <summary>
        /// Sensor for shock capturing after Persson and Peraire (2006). Used for addaptive mesh refinement 
        /// </summary>
        private PerssonSensor perssonsensor;

        /// <summary>
        /// Create the colour and level set distance field. 
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();

            ParticleColor = new SinglePhaseField(new Basis(GridData, 0), "ParticleColor");
            m_RegisteredFields.Add(ParticleColor);
            m_IOFields.Add(ParticleColor);

            LevelSetDistance = new SinglePhaseField(new Basis(GridData, 0), "LevelSetDistance");
            m_RegisteredFields.Add(LevelSetDistance);
            m_IOFields.Add(LevelSetDistance);


            ProcessColor = new SinglePhaseField(new Basis(GridData, 0), "ProcessColor");
            m_RegisteredFields.Add(ProcessColor);
            m_IOFields.Add(ProcessColor);

            if (((FSI_Control)Control).UsePerssonSensor == true) {
                perssonsensor = new PerssonSensor(Pressure);
                this.IOFields.Add(perssonsensor.GetField());
            }
        }

        /// <summary>
        /// Differentiate between splitting and coupled ansatz
        /// </summary>
        private bool UseMovingMesh {
            get {
                switch (((FSI_Control)Control).Timestepper_LevelSetHandling) {
                    case LevelSetHandling.Coupled_Once:
                    case LevelSetHandling.Coupled_Iterative:
                    case LevelSetHandling.FSI_Coupled_Iterative:
                        return true;

                    case LevelSetHandling.LieSplitting:
                    case LevelSetHandling.StrangSplitting:
                    case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    case LevelSetHandling.None:
                        return false;

                    default:
                        throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
                }
            }
        }

        /// <summary>
        /// Array of all local cells with their specific color.
        /// </summary>
        [DataMember]
        private int[] cellColor = null;

        /// <summary>
        /// Array of all particles with their specific color (global).
        /// </summary>
        [DataMember]
        private int[] globalParticleColor = null;

        /// <summary>
        /// Fully coupled LieSplitting?
        /// </summary>
        [DataMember]
        private bool IsFullyCoupled => ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled || ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_Coupled_Iterative;

        [DataMember]
        private bool StaticTimestep => ((FSI_Control)Control).staticTimestep;

        /// <summary>
        /// The maximum timestep setted in the control file.
        /// </summary>
        [DataMember]
        private double DtMax => ((FSI_Control)Control).dtMax;

        /// <summary>
        /// The maximum timestep setted in the control file.
        /// </summary>
        [DataMember]
        private double StartDt => ((FSI_Control)Control).dtMax / 10;

        /// <summary>
        /// FluidDensity
        /// </summary>
        [DataMember]
        private double MaxGridLength => ((FSI_Control)Control).MaxGridLength;


        /// <summary>
        /// FluidDensity
        /// </summary>
        [DataMember]
        private double RefinementLevel => ((FSI_Control)Control).RefinementLevel;

        /// <summary>
        /// FluidDensity
        /// </summary>
        private double GetMinGridLength() {
            if (RefinementLevel == 0)
                return MaxGridLength;
            else
                return MaxGridLength / Math.Pow(2, RefinementLevel);
        }

        /// <summary>
        /// Volume of the fluid domain
        /// </summary>
        [DataMember]
        private double FluidDomainVolume => ((FSI_Control)Control).FluidDomainVolume;

        /// <summary>
        /// FluidViscosity
        /// </summary>
        [DataMember]
        private double FluidViscosity => ((FSI_Control)Control).pureDryCollisions ? 0 : ((FSI_Control)Control).PhysicalParameters.mu_A;

        /// <summary>
        /// FluidDensity
        /// </summary>
        [DataMember]
        private double FluidDensity => ((FSI_Control)Control).pureDryCollisions ? 0 : ((FSI_Control)Control).PhysicalParameters.rho_A;


        /// <summary>
        /// FluidDensity
        /// </summary>
        [DataMember]
        private double MinimalDistanceForCollision => GetMinGridLength() / 3;

        /// <summary>
        /// HydrodynConvergenceCriterion
        /// </summary>
        [DataMember]
        private double HydrodynConvergenceCriterion => ((FSI_Control)Control).hydrodynamicsConvergenceCriterion;

        /// <summary>
        /// Array with two entries (2D). [0] true: x-Periodic, [1] true: y-Periodic
        /// </summary>
        [DataMember]
        private bool[] IsPeriodic => ((FSI_Control)Control).BoundaryIsPeriodic;

        /// <summary>
        /// Array with two entries (2D). [0] true: x-Periodic, [1] true: y-Periodic
        /// </summary>
        [DataMember]
        private bool FixPosition => ((FSI_Control)Control).fixPosition;

        /// <summary>
        /// The position of the (horizontal and vertical) walls.
        /// </summary>
        [DataMember]
        private double[][] BoundaryCoordinates => ((FSI_Control)Control).BoundaryPositionPerDimension;

        /// <summary>
        /// Saves the residual of the hydrodynamic forces and torque of all particles.
        /// </summary>
        private TextWriter logHydrodynamicsResidual;

        /// <summary>
        /// Saves the physical data of all particles
        /// </summary>
        private TextWriter logPhysicalDataParticles;

        /// <summary>
        /// Creates Navier-Stokes and continuity eqution
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            if (IBM_Op != null)
                return;

            // boundary conditions
            boundaryCondMap = new IncompressibleBoundaryCondMap(GridData, Control.BoundaryValues, PhysicsMode.Incompressible);

            // choose the operators
            NSEOperatorConfiguration IBM_Op_config = new NSEOperatorConfiguration {
                convection = Control.PhysicalParameters.IncludeConvection,
                continuity = true,
                Viscous = true,
                PressureGradient = true,
                Transport = true,
                CodBlocks = new bool[] { true, true },
                DomBlocks = new bool[] { true, true },
            };

            string[] CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, spatialDim)).Cat("div");
            string[] Params = ArrayTools.Cat(VariableNames.Velocity0Vector(spatialDim), VariableNames.Velocity0MeanVector(spatialDim));
            string[] DomName = ArrayTools.Cat(VariableNames.VelocityVector(spatialDim), VariableNames.Pressure);

            string[] CodNameSelected = new string[0];
            if (IBM_Op_config.CodBlocks[0])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, spatialDim));
            if (IBM_Op_config.CodBlocks[1])
                CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(spatialDim, 1));

            string[] DomNameSelected = new string[0];
            if (IBM_Op_config.DomBlocks[0])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, spatialDim));
            if (IBM_Op_config.DomBlocks[1])
                DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(spatialDim, 1));

            IBM_Op = new XSpatialOperatorMk2(DomNameSelected, Params, CodNameSelected, (A, B, C) => HMForder, null);

            // Momentum equation
            // =============================
            // Convective part
            // =============================
            {
                if (IBM_Op_config.convection) {
                    for (int d = 0; d < spatialDim; d++) {

                        // The bulk
                        // -----------------------------
                        var comps = IBM_Op.EquationComponents[CodName[d]];
                        var convectionBulk = new LinearizedConvection(spatialDim, boundaryCondMap, d);
                        comps.Add(convectionBulk);

                        // Immersed boundary
                        // -----------------------------
                        if (((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {
                            var convectionAtIB = new Solution.NSECommon.Operator.Convection.FSI_ConvectionAtIB(d, spatialDim, LsTrk, boundaryCondMap,
                                delegate (Vector X) {
                                    throw new NotImplementedException("Currently not implemented for fixed motion");
                                },
                                UseMovingMesh);
                            comps.Add(convectionAtIB);
                        }
                        else {
                            var convectionAtIB = new Solution.NSECommon.Operator.Convection.FSI_ConvectionAtIB(d, spatialDim, LsTrk, boundaryCondMap,
                                    delegate (Vector X) {
                                        return CreateCouplingAtParticleBoundary(X);
                                    },
                                UseMovingMesh);
                            comps.Add(convectionAtIB);
                        }
                    }
                }
            }

            // Pressure part
            // =============================
            for (int d = 0; d < spatialDim; d++) {
                ICollection<IEquationComponent> comps = IBM_Op.EquationComponents[CodName[d]];

                // The bulk
                // -----------------------------
                PressureGradientLin_d pressureBulk = new PressureGradientLin_d(d, boundaryCondMap);
                comps.Add(pressureBulk);

                // Immersed boundary
                // -----------------------------
                Solution.NSECommon.Operator.Pressure.FSI_PressureAtIB pressureAtIB = new Solution.NSECommon.Operator.Pressure.FSI_PressureAtIB(d, spatialDim, LsTrk);
                comps.Add(pressureAtIB);

                // if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                if (this.Control.FixedStreamwisePeriodicBC) {
                    var presSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                    comps.Add(presSource);
                }
            }

            // Viscous part
            // =============================
            for (int d = 0; d < spatialDim; d++) {
                var comps = IBM_Op.EquationComponents[CodName[d]];
                double penalty = this.Control.AdvancedDiscretizationOptions.PenaltySafety;

                // The bulk
                // -----------------------------
                swipViscosity_Term1 viscousBulk = new swipViscosity_Term1(penalty, d, spatialDim, boundaryCondMap, ViscosityOption.ConstantViscosity, FluidViscosity, double.NaN, null);
                comps.Add(viscousBulk);

                // Immersed boundary
                // -----------------------------
                if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                    var viscousAtIB = new Solution.NSECommon.Operator.Viscosity.FSI_ViscosityAtIB(d, spatialDim, LsTrk,
                        penalty, this.ComputePenaltyIB, FluidViscosity, delegate (Vector X) {
                            throw new NotImplementedException("Currently not implemented for fixed motion");
                        });
                    comps.Add(viscousAtIB);
                }
                else {
                    var viscousAtIB = new Solution.NSECommon.Operator.Viscosity.FSI_ViscosityAtIB(d, spatialDim, LsTrk, penalty, ComputePenaltyIB, FluidViscosity,
                        delegate (Vector X) {
                            return CreateCouplingAtParticleBoundary(X);
                        }
                     );
                    comps.Add(viscousAtIB); // immersed boundary component
                }
            }

            // Continuum equation
            // =============================
            {
                for (int d = 0; d < spatialDim; d++) {
                    var src = new Divergence_DerivativeSource(d, spatialDim);
                    var flx = new Divergence_DerivativeSource_Flux(d, boundaryCondMap);
                    IBM_Op.EquationComponents["div"].Add(src);
                    IBM_Op.EquationComponents["div"].Add(flx);
                }

                if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {

                    var divPen = new Solution.NSECommon.Operator.Continuity.DivergenceAtIB(spatialDim, LsTrk, 1,
                        delegate (double[] X, double time) {
                            throw new NotImplementedException("Currently not implemented for fixed motion");
                        });
                    IBM_Op.EquationComponents["div"].Add(divPen);  // immersed boundary component
                }
                else {
                    var divPen = new Solution.NSECommon.Operator.Continuity.FSI_DivergenceAtIB(spatialDim, LsTrk,
                       delegate (Vector X) {
                           return CreateCouplingAtParticleBoundary(X);
                       });
                    IBM_Op.EquationComponents["div"].Add(divPen); // immersed boundary component 
                }
            }

            IBM_Op.Commit();

            CreateTimestepper();
        }
        int bdforder = 0;
        /// <summary>
        /// Creates the BDF-Timestepper
        /// </summary>
        private void CreateTimestepper() {
            SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;
            if (Control.PhysicalParameters.IncludeConvection) {
                SpatialOp = SpatialOperatorType.Nonlinear;
            }

            MassMatrixShapeandDependence MassMatrixShape;
            switch (((FSI_Control)Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.Coupled_Iterative:
                case LevelSetHandling.FSI_Coupled_Iterative:
                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    break;
                case LevelSetHandling.Coupled_Once:
                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.StrangSplitting:
                case LevelSetHandling.None:
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeDependent;
                    break;
                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)this.Control).Timestepper_LevelSetHandling);
            }

            int bdfOrder;
            if (Control.Timestepper_Scheme == FSI_Control.TimesteppingScheme.CrankNicolson)
                bdfOrder = -1;
            else if (Control.Timestepper_Scheme == FSI_Control.TimesteppingScheme.ImplicitEuler)
                bdfOrder = 1;
            else if (Control.Timestepper_Scheme.ToString().StartsWith("BDF"))
                bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
            else
                throw new NotImplementedException("Only Crank-Nicolson, Implicit-Euler and BDFxxx are implemented.");
            bdforder = bdfOrder;
            m_BDF_Timestepper = new XdgBDFTimestepping(
                Fields: ArrayTools.Cat(Velocity, Pressure),
                IterationResiduals: ArrayTools.Cat(ResidualMomentum, ResidualContinuity),
                LsTrk: LsTrk,
                DelayInit: true,
                _ComputeOperatorMatrix: DelComputeOperatorMatrix,
                _ComputeMassMatrix: null,
                _UpdateLevelset: DelUpdateLevelset,
                BDForder: bdfOrder,
                _LevelSetHandling: ((FSI_Control)Control).Timestepper_LevelSetHandling,
                _MassMatrixShapeandDependence: MassMatrixShape,
                _SpatialOperatorType: SpatialOp,
                _MassScale: MassScale,
                _MultigridOperatorConfig: MultigridOperatorConfig,
                _MultigridSequence: MultigridSequence,
                _SpId: FluidSpecies,
                _CutCellQuadOrder: HMForder,
                _AgglomerationThreshold: Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                _useX: true,
                nonlinconfig: Control.NonLinearSolver,
                linearconfig: Control.LinearSolver) {
                m_ResLogger = ResLogger,
                m_ResidualNames = ArrayTools.Cat(ResidualMomentum.Select(f => f.Identification), ResidualContinuity.Identification),
                IterUnderrelax = ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative ? ((FSI_Control)Control).LSunderrelax : 1.0,
                Config_LevelSetConvergenceCriterion = ((FSI_Control)Control).hydrodynamicsConvergenceCriterion,
                SessionPath = SessionPath,
                Timestepper_Init = Solution.Timestepping.TimeStepperInit.SingleInit
            };
        }

        /// <summary>
        /// Returns an array with all coupling parameters. 
        /// </summary>
        private FSI_ParameterAtIB CreateCouplingAtParticleBoundary(Vector X) {
            FSI_ParameterAtIB couplingParameters = null;
            foreach (Particle p in m_Particles) {
                bool containsParticle = m_Particles.Count == 1 ? true : p.Contains(X, 2 * MaxGridLength);
                if (containsParticle) {
                    couplingParameters = new FSI_ParameterAtIB(p, X);
                }
            }
            return couplingParameters;
        }

        /// <summary>
        /// Calls level set update depending on level set handling method.
        /// </summary>
        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {

            switch (((FSI_Control)Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.None:
                    ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                    LevSet.ProjectField(Posfunction);
                    LsTrk.UpdateTracker();
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    UpdateLevelSetParticles(phystime);
                    throw new NotImplementedException("Moving interface solver will be implemented in the near future");

                case LevelSetHandling.Coupled_Once:
                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                case LevelSetHandling.FSI_Coupled_Iterative:
                case LevelSetHandling.StrangSplitting:
                    UpdateLevelSetParticles(phystime);
                    break;

                default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
            }

            /// <summary>
            /// Computes the Residual of the forces and torque acting from to fluid to the particle. Only for coupled iterative level set handling.
            /// </summary>
            double forces_PResidual;
            if (((FSI_Control)this.Control).Timestepper_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                int iterationCounter = 0;
                ParticleHydrodynamics AllParticleHydrodynamics = new ParticleHydrodynamics(LsTrk);
                forces_PResidual = iterationCounter == 0 ? double.MaxValue : AllParticleHydrodynamics.CalculateParticleResidual(ref iterationCounter); ;
                Console.WriteLine("Current forces_PResidual:   " + forces_PResidual);
            }
            // no iterative solver, no residual
            else {
                forces_PResidual = 0;
            }
            return forces_PResidual;
        }

        /// <summary>
        /// Particle to Level-Set-Field 
        /// </summary>
        /// <param name="phystime">
        /// The current time.
        /// </param>
        private void UpdateLevelSetParticles(double phystime) {
            levelSetUpdate = new FSI_LevelSetUpdate(GridData, GetMinGridLength());
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            // Step 0
            // Check for periodic boundaries and create/delete
            // additional particles, appearing at the other side of 
            // the periodic domain
            // =======================================================
            DeleteParticlesOutsideOfDomain();
            CreateGhostParticleAtPeriodicBoundary();
            SwitchGhostAndMasterParticle();
            ProcessColoring();

            // Step 1
            // Define an array with the respective cell colors
            // =======================================================
            cellColor = cellColor == null ? InitializeColoring() : UpdateColoring();

            // Step 2
            // Delete the old level set
            // =======================================================
            DGLevSet.Current.Clear();
            // Step 3
            // Define level set per color
            // =======================================================
            CellMask allParticleMask = null;
            CellMask coloredCellMask = null;
            globalParticleColor = levelSetUpdate.DetermineGlobalParticleColor(GridData, cellColor, m_Particles);
            for (int i = 0; i < globalParticleColor.Length; i++) {
                if (globalParticleColor[i] == 0) {
                    int masterID = m_Particles[i].MasterGhostIDs[0] - 1;
                    globalParticleColor[i] = globalParticleColor[masterID];
                }
            }
            int[] _globalParticleColor = globalParticleColor.CloneAs();
            for (int p = 0; p < _globalParticleColor.Length; p++) {
                // Search for current colour on current process
                int currentColor = _globalParticleColor[p];
                bool processContainsCurrentColor = false;
                BitArray coloredCells = new BitArray(noOfLocalCells);
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (cellColor[j] == currentColor && currentColor != 0) {
                        processContainsCurrentColor = true;
                        coloredCells[j] = true;
                    }
                }
                if (processContainsCurrentColor) {
                    int[] particlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(_globalParticleColor, currentColor);
                    coloredCellMask = new CellMask(GridData, coloredCells);
                    // Save all colored cells of 
                    // any color in one cellmask
                    // -----------------------------
                    allParticleMask = allParticleMask == null ? coloredCellMask : allParticleMask.Union(coloredCellMask);
                    // Particle level set
                    // -----------------------------
                    double levelSetFunction(double[] X, double t) {
                        double levelSetFunctionOneColor = -1;
                        for (int pC = 0; pC < particlesOfCurrentColor.Length; pC++) {
                            Particle currentParticle = m_Particles[particlesOfCurrentColor[pC]];
                            if (levelSetFunctionOneColor < currentParticle.LevelSetFunction(X))
                                levelSetFunctionOneColor = currentParticle.LevelSetFunction(X);
                            _globalParticleColor[particlesOfCurrentColor[pC]] = 0;
                        }
                        return levelSetFunctionOneColor;
                    }
                    SetLevelSet(levelSetFunction, coloredCellMask, phystime);
                }
            }
            // Step 4
            // Define level set of the remaining cells ("Fluid-Cells")
            // =======================================================
            double phiFluid(double[] X, double t) {
                return -1;
            }
            CellMask fluidCells = allParticleMask != null ? allParticleMask.Complement() : CellMask.GetFullMask(GridData);
            SetLevelSet(phiFluid, fluidCells, phystime);
            // Step 5
            // Smoothing
            // =======================================================
            PerformLevelSetSmoothing(allParticleMask, fluidCells, SetFarField: true);
            // Step 6
            // Update level set tracker
            // =======================================================
            try {
                LsTrk.UpdateTracker(__NearRegionWith: 2);
            }
            catch(LevelSetCFLException e) {//hacky workaround
                if (AddedGhostParticle)
                    Console.WriteLine("Ghost particle added, exception thrown: " + e);
                else
                    throw e;
                AddedGhostParticle = false;
            }
        }

        /// <summary>
        /// Set level set based on the function phi and the current cells
        /// </summary>
        /// <param name="levelSetFunction">
        /// The level set function.
        /// </param>
        /// <param name="currentCells">
        /// The cells where the level set function is defined.
        /// </param>
        /// <param name="phystime">
        /// The current time.
        /// </param>
        private void SetLevelSet(Func<double[], double, double> levelSetFunction, CellMask currentCells, double phystime) {
            ScalarFunction Function = NonVectorizedScalarFunction.Vectorize(levelSetFunction, phystime);
            DGLevSet.Current.Clear(currentCells);
            DGLevSet.Current.ProjectField(1.0, Function, new CellQuadratureScheme(UseDefaultFactories: true, domain: currentCells));
        }

        /// <summary>
        /// Initialization of <see cref="ParticleColor"/>  based on particle geometry
        /// </summary>
        private int[] InitializeColoring() {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int JE = GridData.iLogicalCells.NoOfExternalCells + J;
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            int[] coloredCells = new int[J];
            int[] cellsExchange = new int[JE];
            for (int p = 0; p < m_Particles.Count; p++) {
                Particle currentParticle = m_Particles[p];
                for (int j = 0; j < J; j++) {
                    if (currentParticle.Contains(new Vector(CellCenters[j, 0], CellCenters[j, 1]), MaxGridLength)) {
                        ParticleColor.SetMeanValue(j, p + 1);
                        coloredCells[j] = p + 1;
                        cellsExchange[j] = coloredCells[j];
                    }
                }
            }
            cellsExchange.MPIExchange(GridData);
            GridData gridData = (GridData)GridData;
            levelSetUpdate.RecolorCellsOfNeighborParticles(coloredCells, gridData);
            SetColorDGField(coloredCells);
            return coloredCells;
        }

        /// <summary>
        /// Initialization of <see cref="ParticleColor"/>  based on particle geometry
        /// </summary>
        private void ProcessColoring() {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < J; j++) {
                ProcessColor.SetMeanValue(j, MPIRank);
            }
        }

        /// <summary>
        /// Update of <see cref="ParticleColor"/> and <see cref="LevelSetDistance"/>
        /// </summary>
        private int[] UpdateColoring() {
            int[] coloredCells = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int[] coloredCellsExchange = coloredCells.CloneAs();
            coloredCellsExchange.MPIExchange(GridData);
            levelSetUpdate.ColorNeighborCells(coloredCells, coloredCellsExchange);
            GridData gridData = (GridData)GridData;
            levelSetUpdate.RecolorCellsOfNeighborParticles(coloredCells, gridData);
            SetColorDGField(coloredCells);
            return coloredCells;
        }

        private void SetColorDGField(int[] coloredCells) {
            ushort[] regionsCode = LsTrk.Regions.RegionsCode;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < noOfLocalCells; j++) {
                ParticleColor.SetMeanValue(j, coloredCells[j]);
                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(regionsCode[j], 0));
            }
        }

        private void CreateGhostParticleAtPeriodicBoundary() {
            for (int p = 0; p < m_Particles.Count(); p++) {
                Particle currentParticle = m_Particles[p];
                List<Particle> ghostParticles = new List<Particle>();
                int[] ghostHierachy = currentParticle.MasterGhostIDs.CloneAs();
                if (!currentParticle.IsMaster)
                    continue;
                Vector particlePosition = currentParticle.Motion.GetPosition();
                int idOffset = 0;
                for (int d1 = 0; d1 < spatialDim; d1++) { // which direction?
                    if (!IsPeriodic[d1])
                        continue;
                    for(int wallID = 0; wallID < spatialDim; wallID++) { // which wall?
                        if (PeriodicOverlap(currentParticle, d1, wallID)) {
                            ghostHierachy[0] = p + 1;
                            Vector originNeighbouringDomain;
                            if (d1 == 0)
                                originNeighbouringDomain = new Vector(2 * BoundaryCoordinates[0][1 - wallID], 0);
                            else
                                originNeighbouringDomain = new Vector(0, 2 * BoundaryCoordinates[1][1 - wallID]);
                            Particle ghostParticle;
                            if (ghostHierachy[d1 + 1] == 0) {
                                ghostHierachy[d1 + 1] = m_Particles.Count() + idOffset + 1;
                                ghostParticle = currentParticle.CloneAs();
                                ghostParticle.SetGhost();
                                ghostParticle.Motion.SetGhostPosition(originNeighbouringDomain + particlePosition);
                                ghostParticles.Add(ghostParticle.CloneAs());
                            }
                            else{
                                ghostParticle = m_Particles[ghostHierachy[1] - 1];
                            }
                            if (d1 == 0) {
                                idOffset = 1;
                                if (ghostHierachy[3] != 0)
                                    continue;
                                // test for periodic boundaries in y - direction for the newly created ghost
                                for (int wallID2 = 0; wallID2 < spatialDim; wallID2++) {
                                    idOffset += 1;
                                    if (PeriodicOverlap(ghostParticle, 1, wallID2)) {
                                        originNeighbouringDomain = new Vector(0, 2 * BoundaryCoordinates[1][1 - wallID2]);
                                        ghostHierachy[3] = m_Particles.Count() + d1 + idOffset;
                                        ghostParticle = currentParticle.CloneAs();
                                        ghostParticle.SetGhost();
                                        ghostParticle.Motion.SetGhostPosition(originNeighbouringDomain + ghostParticles[0].Motion.GetPosition());
                                        ghostParticles.Add(ghostParticle.CloneAs());
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if(ghostParticles.Count >= 1) {
                    AddedGhostParticle = true;
                    currentParticle.SetGhostHierachy(ghostHierachy);
                    for (int p2 = 0; p2 < ghostParticles.Count(); p2++) {
                        ghostParticles[p2].SetGhostHierachy(ghostHierachy);
                    }
                    m_Particles.AddRange(ghostParticles);
                }
            }
        }

        bool AddedGhostParticle = false;
        private void SwitchGhostAndMasterParticle() {
            for (int p = 0; p < m_Particles.Count(); p++) {
                Particle currentParticle = m_Particles[p];
                if (!currentParticle.IsMaster)
                    continue;
                if (!IsInsideOfDomain(currentParticle)) {
                    int oldMasterID = p + 1;
                    int[] ghostHierachy = currentParticle.MasterGhostIDs;
                    for(int g = 1; g < ghostHierachy.Length; g++) {
                        if (ghostHierachy[g] <= 0)
                            continue;
                        Particle currentGhost = m_Particles[ghostHierachy[g] - 1];
                        if (IsInsideOfDomain(currentGhost)) {
                            int newMasterID = ghostHierachy[g];
                            int[] newGhostHierachy = ghostHierachy.CloneAs();
                            newGhostHierachy[0] = newMasterID;
                            newGhostHierachy[g] = oldMasterID;
                            currentGhost.SetMaster(currentParticle.Motion.CloneAs());
                            currentParticle.SetGhost();
                            for(int i = 0; i < ghostHierachy.Length; i++) {
                                if (ghostHierachy[i] <= 0)
                                    continue;
                                m_Particles[ghostHierachy[i] - 1].MasterGhostIDs = newGhostHierachy.CloneAs();
                            }
                            return;
                        }
                    }
                }
            }
        }

        private void DeleteParticlesOutsideOfDomain() {
            for (int p = 0; p < m_Particles.Count(); p++) {
                for (int d = 0; d < spatialDim; d++) {
                    for (int wallID = 0; wallID < spatialDim; wallID++) {
                        Particle currentParticle = m_Particles[p];
                        Vector particlePosition = currentParticle.Motion.GetPosition();
                        double distance = particlePosition[d] - BoundaryCoordinates[d][wallID];
                        double particleMaxLength = currentParticle.GetLengthScales().Min();
                        if (Math.Abs(distance) > particleMaxLength && !IsInsideOfDomain(currentParticle)) {
                            if (!AnyOverlap(currentParticle)) {
                                int[] ghostHierachy = currentParticle.MasterGhostIDs.CloneAs();
                                for (int g = 0; g < ghostHierachy.Length; g++) {
                                    if(ghostHierachy[g] == p + 1) {
                                        ghostHierachy[g] = 0;
                                    }
                                }
                                for (int g = 0; g < ghostHierachy.Length; g++) {
                                    if(ghostHierachy[g] > 0)
                                        m_Particles[ghostHierachy[g] - 1].MasterGhostIDs = ghostHierachy.CloneAs();
                                }
                                m_Particles.RemoveAt(p);
                                if (p >= m_Particles.Count())// already the last particle, no further action needed!
                                    return;
                                Particle lastParticle = m_Particles.Last();
                                ghostHierachy = lastParticle.MasterGhostIDs.CloneAs();
                                for (int g = 0; g < ghostHierachy.Length; g++) {
                                    if (ghostHierachy[g] == m_Particles.Count() + 1) {
                                        ghostHierachy[g] = p + 1;
                                    }
                                }
                                m_Particles.Insert(p, lastParticle);
                                m_Particles.RemoveAt(m_Particles.Count() - 1);
                                for (int g = 0; g < ghostHierachy.Length; g++) {
                                    if (ghostHierachy[g] > 0)
                                        m_Particles[ghostHierachy[g] - 1].MasterGhostIDs = ghostHierachy.CloneAs();
                                }
                            }
                        }
                    }
                }
            }
        }

        private bool PeriodicOverlap(Particle currentParticle, int d1, int d2) {
            Vector particlePosition = currentParticle.Motion.GetPosition();
            double distance = particlePosition[d1] - BoundaryCoordinates[d1][d2];
            double particleMaxLength = currentParticle.GetLengthScales().Max();
            if (Math.Abs(distance) < particleMaxLength) {
                if (d1 == 0)
                    currentParticle.ClosestPointOnOtherObjectToThis = new Vector(BoundaryCoordinates[d1][d2], particlePosition[1]);
                else
                    currentParticle.ClosestPointOnOtherObjectToThis = new Vector(particlePosition[0], BoundaryCoordinates[d1][d2]);
                FSI_Collision periodicCollision = new FSI_Collision(GetMinGridLength(), 0, 0);
                periodicCollision.CalculateMinimumDistance(currentParticle, out _, out Vector _, out Vector _, out bool Overlapping);
                return Overlapping;
            }
            return false;
        }

        private bool AnyOverlap(Particle currentParticle) {
            for (int d = 0; d < spatialDim; d++) {
                for (int wallID = 0; wallID < spatialDim; wallID++) {
                    Vector particlePosition = currentParticle.Motion.GetPosition();
                    if (d == 0)
                        currentParticle.ClosestPointOnOtherObjectToThis = new Vector(BoundaryCoordinates[d][wallID], particlePosition[1]);
                    else
                        currentParticle.ClosestPointOnOtherObjectToThis = new Vector(particlePosition[0], BoundaryCoordinates[d][wallID]);
                    FSI_Collision periodicCollision = new FSI_Collision(GetMinGridLength(), 0, 0);
                    periodicCollision.CalculateMinimumDistance(currentParticle, out _, out Vector _, out Vector _, out bool Overlapping);
                    if (Overlapping)
                        return true;    
                }   
            }
            return false;
        }

        private bool IsInsideOfDomain(Particle currentParticle) {
            Vector position = currentParticle.Motion.GetPosition();
            for (int d = 0; d < spatialDim; d++) {
                if (!IsPeriodic[d])
                    continue;
                for (int wallID = 0; wallID < spatialDim; wallID++) {
                    Vector wallNormal = new Vector(Math.Sign(BoundaryCoordinates[d][wallID]) * (1 - d), Math.Sign(BoundaryCoordinates[d][wallID]) * d);
                    Vector wallToPoint = d == 0
                        ? new Vector(position[0] - BoundaryCoordinates[0][wallID], position[1])
                        : new Vector(position[0], position[1] - BoundaryCoordinates[1][wallID]);
                    if (wallNormal * wallToPoint > 0)
                        return false;
                }
            }
            return true;
        }

        [DataMember]
        private double oldTimestep;
        [DataMember]
        bool initAddedDamping = true;

        private double CalculateTimestep(double phystime, int TimestepInt) {
            double dt = DtMax;
            if (oldTimestep == 0)
                oldTimestep = StartDt;
            if (StaticTimestep)
                return dt;
            else {
                double maxVelocityL2Norm = 1e-15;
                if (phystime == 0)
                    dt = StartDt;
                else {
                    foreach (Particle p in m_Particles) {
                        double expectedVelocity = p.Motion.GetTranslationalVelocity(0).L2Norm() + (p.Motion.GetHydrodynamicForces(0) * DtMax / p.Motion.ParticleMass).L2Norm();
                        double expectedRotation = (p.Motion.GetRotationalVelocity() + p.Motion.GetHydrodynamicTorque() * DtMax / p.MomentOfInertia) * p.GetLengthScales().Max();
                        maxVelocityL2Norm = Math.Max(maxVelocityL2Norm, expectedVelocity + Math.Abs(expectedRotation));
                    }
                    dt = Math.Min(DtMax, GetMinGridLength() / (2 * maxVelocityL2Norm));
                    if (dt / oldTimestep > 1.1)
                        dt = oldTimestep * 1.1;
                    dt = dt.MPIMin();
                    if (dt < 1e-6)
                        dt = 1e-6;
                    foreach (Particle p in m_Particles) {
                        p.Motion.AdaptParticleHistoryToNewTimeStep(dt, oldTimestep);
                    }
                }
                oldTimestep = dt;
                Console.WriteLine();
                Console.WriteLine("Starting time-step " + TimestepInt + ", size: " + dt + "...............................................................................");
                Console.WriteLine();
                return dt;
            }
        }
        
        //int postCollisionCounter = 0;
        /// <summary>
        /// runs solver one step?!
        /// </summary>
        /// <param name="TimestepInt">
        /// Timestep number
        /// </param>
        /// <param name="phystime">
        /// Physical time
        /// </param>
        /// <param name="dt">
        /// Timestep size
        /// </param>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {
                // init
                ResLogger.TimeStep = TimestepInt;
                BDFSchemeCoeffs[] m_TSCchain = BDFCommon.GetChain(bdforder);
                int S = m_TSCchain[0].S;
                dt = CalculateTimestep(phystime, TimestepInt);
                if (initAddedDamping) {
                    foreach (Particle p in m_Particles) {
                        if (p.Motion.UseAddedDamping) {
                            p.Motion.CalculateDampingTensor(p, LsTrk, FluidViscosity, FluidDensity, DtMax);
                            p.Motion.ExchangeAddedDampingTensors();
                        }
                    }
                    initAddedDamping = false;
                }
                // used later to check if there is exactly one push per timestep
                int OldPushCount = LsTrk.PushCount;

                // only particle motion & collisions, no flow solver
                // =================================================
                if (((FSI_Control)Control).pureDryCollisions) {
                    LsTrk.PushStacks(); // in other branches, called by the BDF timestepper
                    DGLevSet.Push();
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);

                    // physics
                    // -------------------------------------------------
                    ParticleHydrodynamics AllParticleHydrodynamics = new ParticleHydrodynamics(LsTrk);
                    CalculateParticleForcesAndTorque(AllParticleHydrodynamics);
                    CalculateParticleVelocity(m_Particles, dt, 0);
                    CalculateCollision(m_Particles, dt);
                    CalculateParticlePosition(dt);
                    UpdateLevelSetParticles(phystime);

                    // print
                    // -------------------------------------------------
                    Auxillary.PrintResultToConsole(m_Particles, 0, 0, phystime, TimestepInt, ((FSI_Control)Control).FluidDomainVolume);
                    LogPhysicalData(phystime, TimestepInt);

                }
                // particle motion & collisions plus flow solver
                // =================================================
                else {
                    int iterationCounter = 0;
                    double hydroDynForceTorqueResidual = double.MaxValue;
                    int minimumNumberOfIterations = 4;
                    ParticleHydrodynamics AllParticleHydrodynamics = new ParticleHydrodynamics(LsTrk);
                    while (hydroDynForceTorqueResidual > HydrodynConvergenceCriterion || iterationCounter < minimumNumberOfIterations) {
                        if (iterationCounter > ((FSI_Control)Control).maxIterationsFullyCoupled)
                            throw new ApplicationException("No convergence in coupled iterative solver, number of iterations: " + iterationCounter);
                        if (iterationCounter < 0)
                            throw new Exception("Iteration counter is negative?!");
                        if (iterationCounter > 0 && !IsFullyCoupled)
                            throw new Exception("FSI iteration in a non-coupled system, something went wrong");
                        Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                        AllParticleHydrodynamics.SaveHydrodynamicOfPreviousIteration(m_Particles);

                        if (IsFullyCoupled && iterationCounter == 0) {
                            InitializeParticlePerIteration(m_Particles, TimestepInt);
                        }
                        else {
                            m_BDF_Timestepper.Solve(phystime, dt, false);
                            CalculateParticleForcesAndTorque(AllParticleHydrodynamics);
                        }
                        if (((FSI_Control)Control).UsePerssonSensor == true) {
                            perssonsensor.Update(Pressure);
                        }
                        CalculateParticleVelocity(m_Particles, dt, iterationCounter);

                        // not a fully coupled system? -> no iteration
                        // -------------------------------------------------
                        if (!IsFullyCoupled)
                            break;

                        // residual
                        // -------------------------------------------------
                        hydroDynForceTorqueResidual = AllParticleHydrodynamics.CalculateParticleResidual(ref iterationCounter);

                        // print iteration status
                        // -------------------------------------------------
                        Auxillary.PrintResultToConsole(hydroDynForceTorqueResidual, iterationCounter);
                        LogResidual(phystime, iterationCounter, hydroDynForceTorqueResidual);
                    }


                    // collision
                    // -------------------------------------------------
                    CalculateCollision(m_Particles, dt);

                    // particle position
                    // -------------------------------------------------
                    if(!FixPosition)
                        CalculateParticlePosition(dt);

                    // print
                    // -------------------------------------------------
                    Auxillary.PrintResultToConsole(m_Particles, FluidViscosity, FluidDensity, phystime, TimestepInt, FluidDomainVolume);
                    LogPhysicalData(phystime, TimestepInt);

                    // level set tracker 
                    // -------------------------------------------------
                    if (IsFullyCoupled) {// in other branches, called by the BDF timestepper
                        LsTrk.IncreaseHistoryLength(S+1);
                        LsTrk.PushStacks();
                    }
                }

                // finalize
                // ========
                if (LsTrk.PushCount - OldPushCount != 1) {
                    throw new ApplicationException("Illegal number of level-set push actions in time-step." + (LsTrk.PushCount - OldPushCount) + " It is important that LevelSetTracker.PushStacks() is called *exactly once per time-step*, at the beginning.");
                }
                ResLogger.NextTimestep(false);
                return dt;
            }
        }

        private void CalculateParticleForcesAndTorque(ParticleHydrodynamics AllParticleHydrodynamics) {
            ParticleHydrodynamicsIntegration hydrodynamicsIntegration = new ParticleHydrodynamicsIntegration(2, Velocity, Pressure, LsTrk, FluidViscosity);
            AllParticleHydrodynamics.CalculateHydrodynamics(m_Particles, hydrodynamicsIntegration, FluidDensity, IsFullyCoupled);
        }

        /// <summary>
        /// Update of added damping tensors and prediction of hydrdynamics.
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="TimestepInt">
        /// #Timestep
        /// </param>
        internal void InitializeParticlePerIteration(List<Particle> Particles, int TimestepInt) {
            for (int p = 0; p < Particles.Count; p++) {
                Particle currentParticle = Particles[p];
                if (currentParticle.Motion.UseAddedDamping) {
                    currentParticle.Motion.UpdateDampingTensors();
                }
                currentParticle.Motion.SaveHydrodynamicsOfPreviousTimestep();
                currentParticle.Motion.PredictForceAndTorque(currentParticle.ActiveStress, currentParticle.Circumference, TimestepInt, FluidViscosity, FluidDensity,DtMax);
            }
        }

        /// <summary>
        /// Calls the calculation of the Acceleration and the velocity
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="dt">
        /// The time step
        /// </param>
        /// <param name="IterationCounter">
        /// No of iterations
        /// </param>
        internal void CalculateParticleVelocity(List<Particle> Particles, double dt, int IterationCounter) {
            int pID = 0;
            foreach (Particle p in Particles) {
                if (!p.IsMaster) {
                    continue;
                }
                if (IterationCounter == 0) {
                    p.Motion.SaveVelocityOfPreviousTimestep();
                }
                p.Motion.UpdateParticleVelocity(dt);
                for (int g = 1; g < p.MasterGhostIDs.Length; g++) {
                    if (p.MasterGhostIDs[g] < 1)
                        continue;
                    Particle ghost = Particles[p.MasterGhostIDs[g] - 1];
                    if (ghost.IsMaster)
                        throw new Exception("A ghost particle is considered to be a master, that can't be!");
                    if (IterationCounter == 0) {
                        ghost.Motion.SaveVelocityOfPreviousTimestep();
                    }
                    ghost.Motion.CopyNewVelocity(p.Motion.GetTranslationalVelocity(), p.Motion.GetRotationalVelocity());
                }
                pID += 1;
            }
        }

        /// <summary>
        /// Calls the calculation of the position.
        /// </summary>
        /// <param name="dt">
        /// The time step
        /// </param>
        private void CalculateParticlePosition(double dt) {
            for (int p = 0; p < m_Particles.Count; p++) {
                Particle particle = m_Particles[p];
                particle.Motion.UpdateParticlePositionAndAngle(dt);
            }
        }

        /// <summary>
        /// Creates a log file for the residum of the hydrodynamic forces.
        /// </summary>
        private void CreateResidualLogger() {
            if ((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                logHydrodynamicsResidual = DatabaseDriver.FsDriver.GetNewLog("HydrodynamicResidual", CurrentSessionInfo.ID);
                logHydrodynamicsResidual.WriteLine(string.Format("{0},{1},{2}", "Time", "Iteration", "Residual"));
            }
        }

        /// <summary>
        /// Creates a log file for the residum of the hydrodynamic forces.
        /// </summary>
        private void LogResidual(double phystime, int iterationCounter, double residual) {
            if ((MPIRank == 0) && (logPhysicalDataParticles != null)) {
                logHydrodynamicsResidual.WriteLine(string.Format("{0},{1},{2}", phystime, iterationCounter, residual));
                logHydrodynamicsResidual.Flush();
            }
        }

        /// <summary>
        /// Creates a log file for the physical data of the particles. Only active if a database is specified.
        /// </summary>
        private void CreatePhysicalDataLogger() {
            if ((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                logPhysicalDataParticles = DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", "timestep", "particle", "time", "posX", "posY", "angle", "velX", "velY", "rot", "fX", "fY", "T"));
            }
        }

        /// <summary>
        /// Writes the physical data of the particles to a log file.
        /// </summary>
        /// <param name = phystime>
        /// </param>
        private void LogPhysicalData(double phystime, int timestepNo) {
            if ((MPIRank == 0) && (logPhysicalDataParticles != null)) {
                for (int p = 0; p < m_Particles.Count(); p++) {
                    logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", timestepNo, p, phystime, m_Particles[p].Motion.GetPosition(0)[0], m_Particles[p].Motion.GetPosition(0)[1], m_Particles[p].Motion.GetAngle(0), m_Particles[p].Motion.GetTranslationalVelocity(0)[0], m_Particles[p].Motion.GetTranslationalVelocity(0)[1], m_Particles[p].Motion.GetRotationalVelocity(0), m_Particles[p].Motion.GetHydrodynamicForces(0)[0], m_Particles[p].Motion.GetHydrodynamicForces(0)[1], m_Particles[p].Motion.GetHydrodynamicTorque(0)));
                    logPhysicalDataParticles.Flush();
                }
            }
        }

        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {
            FSI_TimestepInfo tsi = new FSI_TimestepInfo(t, CurrentSessionInfo, timestepno, IOFields, m_Particles);
            SerialzeTester(tsi);
            return tsi;
        }

        /// <summary>
        /// Test the serialization of <see cref="FSI_TimestepInfo.Particles"/>
        /// </summary>
        private static void SerialzeTester(FSI_TimestepInfo b) {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Serialize
            };
            bool DebugSerialization = false;
            JsonReader GetJsonReader(Stream s) {
                if (DebugSerialization) {
                    return new JsonTextReader(new StreamReader(s));
                }
                else {
                    return new BsonReader(s);
                }
            }
            JsonWriter GetJsonWriter(Stream s) {
                if (DebugSerialization) {
                    return new JsonTextWriter(new StreamWriter(s));
                }
                else {
                    return new BsonWriter(s);
                }
            }
            byte[] buffer = null;
            using (var ms1 = new MemoryStream()) {
                using (var writer = GetJsonWriter(ms1)) {
                    formatter.Serialize(writer, b);
                    writer.Flush();
                    buffer = ms1.GetBuffer();
                    //writer.Close();
                }
            }
            FSI_TimestepInfo o;
            using (var ms2 = new MemoryStream(buffer)) {
                using (var reader = GetJsonReader(ms2)) {
                    o = formatter.Deserialize<FSI_TimestepInfo>(reader);
                    reader.Close();
                }
            }
            Debug.Assert(b.Particles.Length == o.Particles.Length);
            int L = b.Particles.Length;
            for (int l = 0; l < L; l++) { // loop over particles
                //Debug.Assert(GenericBlas.L2Dist((double[])b.Particles[l].Motion.GetPosition(0), (double[])o.Particles[l].Motion.GetPosition(0)) < 1e-13);
            }
        }

        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override void OnRestartTimestepInfo(TimestepInfo tsi) {
            FSI_TimestepInfo fTsi = (FSI_TimestepInfo)tsi;
            // init particles
            m_Particles = fTsi.Particles.ToList();
            UpdateLevelSetParticles(fTsi.PhysicalTime);
        }


        /// <summary>
        /// For restarting calculations, its important to reload old solutions if one uses a higher order method in time
        /// </summary>
        /// <param name="time"></param>
        /// <param name="timestep"></param>
        public override void PostRestart(double time, TimestepNumber timestep) {
            if (!((FSI_Control)Control).IsRestart)
                return;
            var fsDriver = this.DatabaseDriver.FsDriver;
            string pathToOldSessionDir = System.IO.Path.Combine(
                fsDriver.BasePath, "sessions", this.CurrentSessionInfo.RestartedFrom.ToString());
            string pathToPhysicalData = System.IO.Path.Combine(pathToOldSessionDir, "PhysicalData.txt");
            string[] records = File.ReadAllLines(pathToPhysicalData);
            int timestepIndexOffset = 0;
            for (int r = 1; r < records.Length; r++) {
                string currentLine = records[r];
                string[] currentLineFields = currentLine.Split(',');
                if (timestep.MajorNumber == Convert.ToInt32(currentLineFields[0])) {
                    timestepIndexOffset = r;
                    break;
                }
            }
            for (int p = 0; p < m_Particles.Count(); p++) {
                Particle currentParticle = m_Particles[p];
                int index = timestepIndexOffset + p;
                string currentLine = records[index];
                string[] currentLineFields = currentLine.Split(',');
                double[] position = new double[2];
                double[] translationalVelocity = new double[2];
                position[0] = Convert.ToDouble(currentLineFields[3]);
                position[1] = Convert.ToDouble(currentLineFields[4]);
                double angle = Convert.ToDouble(currentLineFields[5]) * 360 / (2 * Math.PI);
                translationalVelocity[0] = Convert.ToDouble(currentLineFields[6]);
                translationalVelocity[1] = Convert.ToDouble(currentLineFields[7]);
                double angularVelocity = Convert.ToDouble(currentLineFields[8]);
                currentParticle.Motion.InitializeParticlePositionAndAngle(position, angle);
                currentParticle.Motion.InitializeParticleVelocity(translationalVelocity, angularVelocity);
            }
            cellColor = null;
            UpdateLevelSetParticles(time);
            CreatePhysicalDataLogger();
            CreateResidualLogger();
            ProcessColoring();
        }

        /// <summary>
        /// Triggers the collision detection, which triggers the calculation of the collisions
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="CellColor">
        /// All cells on the current process with their specific colour
        /// </param>
        /// <param name="dt">
        /// Timestep
        /// </param>
        private void CalculateCollision(List<Particle> Particles, double dt) {
            if (((FSI_Control)Control).collisionModel == FSI_Control.CollisionModel.NoCollisionModel)
                return;
            if (((FSI_Control)Control).collisionModel == FSI_Control.CollisionModel.RepulsiveForce)
                throw new NotImplementedException("Repulsive force model is no longer supported, please use the momentum conservation model.");
            foreach (Particle p in Particles) {
                p.IsCollided = false;
            }
            // Only particles with the same colour a close to each other, thus, we only test for collisions within those particles.
            // Determine colour.
            // =================================================
            int[] _GlobalParticleColor = globalParticleColor.CloneAs();
            for (int i = 0; i < _GlobalParticleColor.Length; i++) {
                int CurrentColor = _GlobalParticleColor[i];
                if (CurrentColor == 0)
                    continue;
                int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(_GlobalParticleColor, CurrentColor);
                // Multiple particles with the same colour, trigger collision detection
                // =================================================
                if (ParticlesOfCurrentColor.Length >= 1 && CurrentColor != 0) {
                    Particle[] currentParticles = new Particle[ParticlesOfCurrentColor.Length];
                    for (int j = 0; j < ParticlesOfCurrentColor.Length; j++) {
                        currentParticles[j] = m_Particles[ParticlesOfCurrentColor[j]];
                    }
                    FSI_Collision _Collision = new FSI_Collision(GetMinGridLength(), ((FSI_Control)Control).CoefficientOfRestitution, dt, ((FSI_Control)Control).WallPositionPerDimension, ((FSI_Control)Control).BoundaryIsPeriodic, MinimalDistanceForCollision);
                    _Collision.CalculateCollision(currentParticles);
                }
                // Remove already investigated particles/colours from array
                // =================================================
                for (int j = 0; j < _GlobalParticleColor.Length; j++) {
                    if (_GlobalParticleColor[j] == CurrentColor)
                        _GlobalParticleColor[j] = 0;
                }
            }
            // Communicate
            // =================================================
            foreach (Particle p in m_Particles) {
                Collision_MPICommunication(p, MPISize);
            }
        }

        /// <summary>
        /// Ensures the communication between the processes after a collision. As collisions are triggered based on the (local) colouring of the cells only the owning process knows about them.
        /// Thus, it is necessary to inform the other processes about the collisions.
        /// </summary>
        /// <param name="currentParticle">
        /// The current particle.
        /// </param>
        /// <param name="MPISize">
        /// Number of mpi processes
        /// </param>
        private void Collision_MPICommunication(Particle currentParticle, int MPISize) {
            // Did a collision take place on one of the processes?
            // ===================================================
            double[] isCollidedSend = new double[1];
            isCollidedSend[0] = currentParticle.IsCollided ? 1 : 0;
            double[] isCollidedReceive = new double[MPISize];
            MPISendAndReceive(isCollidedSend, ref isCollidedReceive);

            bool noCurrentCollision = true;
            for (int i = 0; i < isCollidedReceive.Length; i++) {
                // The particle is collided, thus, copy the data from
                // the owning process.
                // ===================================================
                if (isCollidedReceive[i] != 0) {
                    Console.WriteLine("collided on " + i + " length " + isCollidedReceive.Length);
                    double[] dataSend = currentParticle.Motion.BuildSendArray();
                    int noOfVars = dataSend.Length;
                    double[] dataReceive = new double[noOfVars * MPISize];
                    MPISendAndReceive(dataSend, ref dataReceive);
                    currentParticle.Motion.WriteReceiveArray(dataReceive, offset: i * noOfVars);
                    currentParticle.IsCollided = true;
                    noCurrentCollision = false;
                }
            }
            // nothing happend
            // ===================================================
            if (noCurrentCollision) {
                currentParticle.IsCollided = false;
            }
        }

        /// <summary>
        /// MPI communication of an array.
        /// </summary>
        /// <param name="variableSend">
        /// The array to send.
        /// </param>
        /// <param name="variableReceive">
        /// An array of the data of all processes.
        /// </param>
        private void MPISendAndReceive(double[] variableSend, ref double[] variableReceive) {
            unsafe {
                fixed (double* pVariableSend = variableSend, pVariableReceive = variableReceive) {
                    csMPI.Raw.Allgather((IntPtr)pVariableSend, variableSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pVariableReceive, variableSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                }
            }
        }

        /// <summary>
        /// Adaptive mesh refinement
        /// </summary>
        /// <param name="TimestepNo">
        /// Currently unused.
        /// </param>
        /// <param name="newGrid">
        /// The adapted grid.
        /// </param>
        /// <param name="old2NewGrid">
        /// The correlation between old and new grid.
        /// </param>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            GridData gridData = (GridData)GridData;
            bool AnyChangeInGrid = false;
            List<int> CellsToRefineList = new List<int>();
            List<int[]> Coarsening = new List<int[]>();
            if (((FSI_Control)Control).AdaptiveMeshRefinement) {
                BitArray cutCells = LsTrk.Regions.GetCutCellMask().GetBitMask();
                if (TimestepNo < 1 || ((FSI_Control)Control).ConstantRefinement)
                    AnyChangeInGrid = GridRefinementController.ComputeGridChange(gridData, cutCells, GetCellMaskWithRefinementLevelsStartUpSweeps(), out CellsToRefineList, out Coarsening);
                else
                    AnyChangeInGrid = GridRefinementController.ComputeGridChange(gridData, cutCells, CHAOSGetCellMaskWithRefinementLevelsStartUpSweeps(), out CellsToRefineList, out Coarsening);
            } 
            if (AnyChangeInGrid) {
                int[] consoleRefineCoarse = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                int oldJ = GridData.CellPartitioning.TotalLength;
                Console.WriteLine("Refining " + consoleRefineCoarse[0] + " of " + oldJ + " cells");
                Console.WriteLine("Coarsening " + consoleRefineCoarse[1] + " of " + oldJ + " cells");
                Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Count().MPISum());
                newGrid = gridData.Adapt(CellsToRefineList, Coarsening, out old2NewGrid);
            }
            else {
                newGrid = null;
                old2NewGrid = null;
            }
        }

        /// <summary>
        /// Finds all cells to be refined with the perssonsensor.
        /// </summary>
        private List<Tuple<int, BitArray>> GetCellMaskWithRefinementLevelsWithPersson() {
            int refinementLevel = ((FSI_Control)Control).RefinementLevel;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            BitArray perssonCells = new BitArray(noOfLocalCells);
            SpeciesId fluidSpeciesID = LsTrk.GetSpeciesId("A");
            for (int p = 0; p < m_Particles.Count(); p++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (LsTrk.Regions.IsSpeciesPresentInCell(fluidSpeciesID, j)) {
                        double perssonValue = perssonsensor.GetValue(j);
                        double threshold = 1e-6;
                        if (!perssonCells[j] && perssonValue > threshold) {
                            Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                            perssonCells[j] = m_Particles[p].Contains(cellCenter, m_Particles[p].GetLengthScales().Max());
                        }
                        else if (!perssonCells[j]) {
                            Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                            perssonCells[j] = m_Particles[p].Contains(cellCenter, GetMinGridLength() * 4);
                        }
                    }
                    if (!perssonCells[j] && ((FSI_Control)Control).pureDryCollisions) {
                        Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                        perssonCells[j] = m_Particles[p].Contains(cellCenter, MaxGridLength / 2);
                    }
                }
            }
            
            return new List<Tuple<int, BitArray>> { new Tuple<int, BitArray>(refinementLevel, perssonCells), };
        }
        /// <summary>
        /// Finds all cells to be refined with the perssonsensor.
        /// </summary>
        private List<Tuple<int, BitArray>> CHAOSGetCellMaskWithRefinementLevelsStartUpSweeps() {
            int refinementLevel = ((FSI_Control)Control).RefinementLevel;
            int refinementLevel2 = refinementLevel >= 2 ? ((FSI_Control)Control).RefinementLevel / 2 : 1;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            BitArray particleNearRegionCells = new BitArray(noOfLocalCells);
            BitArray particleFarRegionCells = new BitArray(noOfLocalCells);
            for (int p = 0; p < m_Particles.Count(); p++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (!particleNearRegionCells[j]) {
                        Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                        particleNearRegionCells[j] = m_Particles[p].Contains(cellCenter, m_Particles[p].GetLengthScales().Min() / 2);
                    }
                    if (!particleFarRegionCells[j]) {
                        Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                        particleFarRegionCells[j] = m_Particles[p].Contains(cellCenter, m_Particles[p].GetLengthScales().Min());
                    }
                }
            }
            return new List<Tuple<int, BitArray>> { new Tuple<int, BitArray>(refinementLevel, particleNearRegionCells),
                                    new Tuple<int, BitArray>(refinementLevel2, particleNearRegionCells),
            };
        }
        /// <summary>
        /// Finds all cells to be refined with the perssonsensor.
        /// </summary>
        private List<Tuple<int, BitArray>> GetCellMaskWithRefinementLevelsStartUpSweeps() {
            int refinementLevel = ((FSI_Control)Control).RefinementLevel;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            BitArray particleNearRegionCells = new BitArray(noOfLocalCells);
            for (int p = 0; p < m_Particles.Count(); p++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (!particleNearRegionCells[j]) {
                        Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                        particleNearRegionCells[j] = m_Particles[p].Contains(cellCenter, m_Particles[p].GetLengthScales().Min());
                    }
                }
            }
            return new List<Tuple<int, BitArray>> { new Tuple<int, BitArray>(3, particleNearRegionCells), };
        }
    }
}