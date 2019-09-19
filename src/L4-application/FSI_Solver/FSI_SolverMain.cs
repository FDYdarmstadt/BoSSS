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

namespace BoSSS.Application.FSI_Solver {
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain {

        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {
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
            UpdateLevelSetParticles(phystime: 0.0);
            base.SetInitial();
        }

        /// <summary>
        /// A list for all particles
        /// </summary>
        private List<Particle> m_Particles;

        /// <summary>
        /// External access to particle list.
        /// </summary>
        public IList<Particle> GetParticles() {
            return m_Particles;
        }

        /// <summary>
        /// An object with some additional methods
        /// </summary>
        readonly private FSI_Auxillary Auxillary = new FSI_Auxillary();

        /// <summary>
        /// Particle color field. The level set is only defined on colored cells.
        /// </summary>
        private SinglePhaseField ParticleColor;

        /// <summary>
        /// Level set distance field. 
        /// </summary>
        private SinglePhaseField LevelSetDistance;

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
        }

        /// <summary>
        /// Differentiate between splitting and coupled ansatz
        /// </summary>
        private bool UseMovingMesh {
            get {
                switch (((FSI_Control)Control).Timestepper_LevelSetHandling) {
                    case LevelSetHandling.Coupled_Once:
                    case LevelSetHandling.Coupled_Iterative:
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
        private int[] cellColor = null;

        /// <summary>
        /// Array of all particles with their specific color (global).
        /// </summary>
        private int[] globalParticleColor = null;

        /// <summary>
        /// Bool to prevent restricting the level set in the stup phase.
        /// </summary>
        private bool setupLevelSet = false;

        /// <summary>
        /// Fully coupled LieSplitting?
        /// </summary>
        private bool IsFullyCoupled => ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSI_LieSplittingFullyCoupled;

        /// <summary>
        /// The maximum timestep setted in the control file.
        /// </summary>
        private double DtMax => ((FSI_Control)Control).dtMax;

        /// <summary>
        /// FluidViscosity
        /// </summary>
        private double FluidViscosity => ((FSI_Control)Control).pureDryCollisions ? 0 : ((FSI_Control)Control).PhysicalParameters.mu_A;

        /// <summary>
        /// FluidDensity
        /// </summary>
        private double FluidDensity => ((FSI_Control)Control).pureDryCollisions ? 0 : ((FSI_Control)Control).PhysicalParameters.rho_A;

        /// <summary>
        /// HydrodynConvergenceCriterion
        /// </summary>
        private double HydrodynConvergenceCriterion => ((FSI_Control)Control).hydrodynamicsConvergenceCriterion;

        /// <summary>
        /// Only for added damping. Check whether the tensors are calculated or not. (move to particle.Motion.cs in the future)
        /// </summary>
        private bool CalculatedDampingTensors = false;

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

            int spatialDim = GridData.SpatialDimension;
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
                            var convectionAtIB = new Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, spatialDim, LsTrk,
                                this.Control.AdvancedDiscretizationOptions.LFFA, boundaryCondMap,
                                delegate (double[] X, double time) {
                                    throw new NotImplementedException("Currently not implemented for fixed motion");
                                },
                                FluidDensity,
                                UseMovingMesh);
                            comps.Add(convectionAtIB);
                        }
                        else {
                            var convectionAtIB = new Solution.NSECommon.Operator.Convection.ActiveConvectionAtIB(d, spatialDim, LsTrk,
                                Control.AdvancedDiscretizationOptions.LFFA, boundaryCondMap,
                                    delegate (double[] X, double time) {
                                        return CreateCouplingAtParticleBoundary(X);
                                    },
                                FluidDensity,
                                UseMovingMesh);
                            comps.Add(convectionAtIB);
                        }
                    }
                    U0MeanRequired = true;
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
                Solution.NSECommon.Operator.Pressure.ActivePressureAtIB pressureAtIB = new Solution.NSECommon.Operator.Pressure.ActivePressureAtIB(d, spatialDim, LsTrk);
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

                    var viscousAtIB = new Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, spatialDim, LsTrk,
                        penalty, this.ComputePenaltyIB,
                        FluidViscosity / FluidDensity,
                        delegate (double[] X, double time) {
                            throw new NotImplementedException("Currently not implemented for fixed motion");
                        });
                    comps.Add(viscousAtIB);
                }
                else {
                    var viscousAtIB = new Solution.NSECommon.Operator.Viscosity.ActiveViscosityAtIB(d, spatialDim, LsTrk,
                        penalty, this.ComputePenaltyIB,
                        FluidViscosity / FluidDensity,
                        delegate (double[] X, double time) {
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
                    var divPen = new Solution.NSECommon.Operator.Continuity.ActiveDivergenceAtIB(spatialDim, LsTrk, 1,
                       delegate (double[] X, double time) {
                           return CreateCouplingAtParticleBoundary(X);
                       });
                    IBM_Op.EquationComponents["div"].Add(divPen); // immersed boundary component 
                }
            }

            IBM_Op.Commit();

            CreateTimestepper();
        }

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
                    MassMatrixShape = MassMatrixShapeandDependence.IsTimeAndSolutionDependent;
                    break;

                case LevelSetHandling.Coupled_Once:
                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
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
        private double[] CreateCouplingAtParticleBoundary(double[] X) {
            double[] couplingArray = new double[X.Length + 7];
            foreach (Particle p in m_Particles) {
                p.CalculateRadialNormalVector(X, out double[] RadialNormalVector);
                double seperateBoundaryRegions = p.activeStress != 0 ? p.SeperateBoundaryRegions(X) : 0;
                bool containsParticle = m_Particles.Count == 1 ? true : p.Contains(X, GridData.iGeomCells.h_min.Min());
                if (containsParticle) {
                    couplingArray[0] = p.Motion.translationalVelocity[0][0];
                    couplingArray[1] = p.Motion.translationalVelocity[0][1];
                    couplingArray[2] = p.Motion.rotationalVelocity[0];
                    couplingArray[3] = RadialNormalVector[0];
                    couplingArray[4] = RadialNormalVector[1];
                    couplingArray[5] = p.Motion.position[0].L2Distance(X);
                    couplingArray[6] = p.activeStress; // zero for passive particles
                    couplingArray[7] = -seperateBoundaryRegions;
                    couplingArray[8] = p.Motion.angle[0];
                }
            }
            return couplingArray;
        }
        
        /// <summary>
        /// Calls level set update depending on level set handling method.
        /// </summary>
        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {

            switch (((FSI_Control)this.Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.None:
                    ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                    LevSet.ProjectField(Posfunction);
                    LsTrk.UpdateTracker();
                    break;

                case LevelSetHandling.Coupled_Once:
                    UpdateLevelSetParticles(phystime);
                    break;

                case LevelSetHandling.Coupled_Iterative:
                    Console.WriteLine("WARNING: Coupled iterative solver is not tested!");
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                    CalculateHydrodynamicForces(m_Particles, dt);
                    UpdateLevelSetParticles(phystime);
                    foreach (Particle p in m_Particles) {
                        p.iteration_counter_P += 1;
                        p.forceAndTorque_convergence = ((FSI_Control)this.Control).hydrodynamicsConvergenceCriterion;
                    }
                    break;

                case LevelSetHandling.LieSplitting:
                    UpdateLevelSetParticles(phystime);
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    UpdateLevelSetParticles(phystime);
                    if (!CalculatedDampingTensors) {
                        foreach (Particle p in m_Particles) {
                            if (p.Motion.useAddedDamping) {
                                p.Motion.CalculateDampingTensor(p, LsTrk, FluidViscosity, FluidDensity, DtMax);
                                Auxillary.ExchangeDampingTensors(m_Particles);
                            }
                        }
                    }
                    CalculatedDampingTensors = true;
                    break;

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
                forces_PResidual = iterationCounter == 0 ? double.MaxValue : Auxillary.CalculateParticleResidual(m_Particles, ref iterationCounter);
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
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate(LsTrk);
            // Step 1
            // Define an array with the respective cell colors
            // =======================================================
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
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

                    // Get particle level set
                    // -----------------------------
                    double levelSetFunction(double[] X, double t) {
                        // Generating the correct sign
                        double levelSetFunctionOneColor = Math.Pow(-1, particlesOfCurrentColor.Length - 1);
                        // Multiplication over all particle-level-sets within the current color
                        for (int pC = 0; pC < particlesOfCurrentColor.Length; pC++) {
                            Particle currentParticle = m_Particles[particlesOfCurrentColor[pC]];
                            double tempLevelSetFunction = currentParticle.levelSetFunction(X);
                            // prevent extreme values
                            if (tempLevelSetFunction > 1)
                                tempLevelSetFunction = 1;
                            if (tempLevelSetFunction < -1 && setupLevelSet) {
                                tempLevelSetFunction = -1;
                            }
                            levelSetFunctionOneColor *= tempLevelSetFunction;
                            // Delete the particle within the current color from the particle color array
                            _globalParticleColor[particlesOfCurrentColor[pC]] = 0;
                        }
                        return levelSetFunctionOneColor;
                    }
                    // Set particle level set
                    // -----------------------------
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
            LsTrk.UpdateTracker(__NearRegionWith: 2);
            setupLevelSet = true;
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
        /// Update of <see cref="ParticleColor"/> and <see cref="LevelSetDistance"/>
        /// </summary>
        private int[] UpdateColoring() {
            // Step 1
            // Color all cells directlly related to the level set
            // ======================================================
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            ushort[] regionsCode = LsTrk.Regions.RegionsCode;
            int[] particleColor = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int[] particleColorExchange = particleColor.CloneAs();

            // No particles on current proc
            // -----------------------------
            if (m_Particles.Count == 0)
                return particleColor;
            particleColorExchange.MPIExchange(GridData);

            // Step 2
            // Color neighbour cells
            // =======================================================
            int neighbourSearchDepth = 2;
            for (int k = 0; k < neighbourSearchDepth; k++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (particleColorExchange[CellNeighbors[i]] != 0 && particleColorExchange[j] == 0) {
                            particleColor[j] = particleColorExchange[CellNeighbors[i]];
                        }
                    }
                }
                particleColorExchange = particleColor.CloneAs();
                particleColorExchange.MPIExchange(GridData);
            }

            // Step 4
            // Find neighbouring colours and recolour one of them
            // =======================================================
            // Find neighbouring cells with
            // different colours
            // -----------------------------
            int maxColor = particleColor.Max().MPIMax();
            int[,] colorToRecolorWith = new int[maxColor + 1, 2];
            for (int j = 0; j < noOfLocalCells; j++) {
                GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                for (int i = 0; i < CellNeighbors.Length; i++) {
                    if (particleColorExchange[CellNeighbors[i]] != particleColor[j] && particleColor[j] != 0 && particleColorExchange[CellNeighbors[i]] > 0) {
                        if (particleColorExchange[CellNeighbors[i]] < particleColor[j] || colorToRecolorWith[particleColor[j], 1] > particleColorExchange[CellNeighbors[i]]) {
                            colorToRecolorWith[particleColor[j], 0] = particleColor[j];
                            colorToRecolorWith[particleColor[j], 1] = particleColorExchange[CellNeighbors[i]];
                        }
                        if (particleColorExchange[CellNeighbors[i]] > particleColor[j]) {
                            if (colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 0] == 0 || colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 1] > particleColor[j]) {
                                colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 0] = particleColorExchange[CellNeighbors[i]];
                                colorToRecolorWith[particleColorExchange[CellNeighbors[i]], 1] = particleColor[j];
                            }
                        }
                    }
                }
            }

            // Communicate
            // -----------------------------
            int[][,] GlobalColorToRecolorWith = colorToRecolorWith.MPIGatherO(0);
            GlobalColorToRecolorWith = GlobalColorToRecolorWith.MPIBroadcast(0);
            for (int m = 0; m < MPISize; m++) {
                for (int i = 0; i < maxColor + 1; i++) {
                    if (GlobalColorToRecolorWith[0][i, 1] == 0 || GlobalColorToRecolorWith[0][i, 1] > GlobalColorToRecolorWith[m][i, 1] && GlobalColorToRecolorWith[m][i, 1] != 0) {
                        GlobalColorToRecolorWith[0][i, 0] = GlobalColorToRecolorWith[m][i, 0];
                        GlobalColorToRecolorWith[0][i, 1] = GlobalColorToRecolorWith[m][i, 1];
                    }
                }
            }
            colorToRecolorWith = GlobalColorToRecolorWith[0];

            // Recolor
            // -----------------------------
            for (int i = maxColor; i > 0; i--) {
                if (colorToRecolorWith[i, 0] != 0) {
                    for (int j = 0; j < noOfLocalCells; j++) {
                        if (particleColor[j] == colorToRecolorWith[i, 0]) {
                            particleColor[j] = colorToRecolorWith[i, 1];
                        }
                    }
                }
            }

            for (int j = 0; j < noOfLocalCells; j++) {
                ParticleColor.SetMeanValue(j, particleColor[j]);
                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(regionsCode[j], 0));
            }
            return particleColor;
        }

        /// <summary>
        /// Initialization of <see cref="ParticleColor"/>  based on particle geometry
        /// </summary>
        private int[] InitializeColoring() {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int JE = GridData.iLogicalCells.NoOfExternalCells + J;
            int[] cells = new int[J];
            int[] cellsExchange = new int[JE];
            List<int> coloredCells = new List<int>();
            for (int p = 0; p < m_Particles.Count; p++) {
                Particle currentParticle = m_Particles[p];
                double h_min = GridData.iGeomCells.h_min.Min() / 2;
                for (int j = 0; j < J; j++) {
                    double[] center = GridData.iLogicalCells.GetCenter(j);
                    // Check for every cell whether their center is part of a particle or not (with tolerance sqrt(h_max^2+h_min^2))
                    if (currentParticle.Contains(center, h_min, h_min)) {
                        ParticleColor.SetMeanValue(j, p + 1);
                        coloredCells.Add(j);
                        cells[j] = p + 1;
                        cellsExchange[j] = cells[j];
                    }
                }
            }
            cellsExchange.MPIExchange(GridData);
            FixNeighbourColoring(cellsExchange);
            for(int j = 0; j < J; j++) {
                cells[j] = cellsExchange[j];
            }
            return cells;
        }

        /// <summary>
        /// Checks whether there are two different colours neighbouring each other. 
        /// </summary>
        /// <param name="coloredCells">
        /// All cells with their colour, uncoloured cells are set to zero.
        /// </param>
        private void FixNeighbourColoring(int[] coloredCells) {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int i = 0; i < J; i++) {
                if (coloredCells[i] != 0) {
                    GridData.GetCellNeighbours(i, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int j = 0; j < CellNeighbors.Length; j++) {
                        if (CellNeighbors[j] < coloredCells.Max() && coloredCells[i] != coloredCells[CellNeighbors[j]] && coloredCells[CellNeighbors[j]] != 0) {
                            RecolorCells(coloredCells, coloredCells[i], coloredCells[CellNeighbors[j]]);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Recolours all cells with a specific colour.
        /// </summary>
        /// <param name="coloredCells">
        /// All cells with their colour, uncoloured cells are set to zero.
        /// </param>
        /// /// <param name="newColor">
        /// </param>
        /// /// <param name="oldColor">
        /// </param>
        private void RecolorCells(int[] coloredCells, int newColor, int oldColor) {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int i = 0; i < J; i++) {
                if (coloredCells[i] == oldColor) {
                    coloredCells[i] = newColor;
                    ParticleColor.SetMeanValue(i, newColor);
                }
            }
        }

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
                dt = GetFixedTimestep();
                Console.WriteLine("Starting time-step " + TimestepInt + "...");
                // used later to check if there is exactly one push per timestep
                int OldPushCount = LsTrk.PushCount;

                // only particle motion & collisions, no flow solver
                // =================================================
                if (((FSI_Control)Control).pureDryCollisions) {
                    if (phystime == 0) { CreatePhysicalDataLogger(); }
                    LsTrk.PushStacks(); // in other branches, called by the BDF timestepper
                    DGLevSet.Push();

                    foreach (Particle p in m_Particles) {
                        p.Motion.GetParticleDensity(p.particleDensity);
                    }
                    Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);

                    // physics
                    // -------------------------------------------------
                    CalculateHydrodynamicForces(m_Particles, dt);
                    CalculateParticleVelocity(m_Particles, dt, 0);
                    CalculateCollision(m_Particles, cellColor, dt);
                    CalculateParticlePosition(dt);
                    UpdateLevelSetParticles(phystime);

                    // print
                    // -------------------------------------------------
                    Auxillary.PrintResultToConsole(m_Particles, 0, 0, phystime, TimestepInt, out double MPIangularVelocity, out Test_Force);
                    LogPhysicalData(phystime);
                    SaveForNUnitTest(MPIangularVelocity);

                }
                // particle motion & collisions plus flow solver
                // =================================================
                else {
                    if (((FSI_Control)Control).Timestepper_LevelSetHandling != LevelSetHandling.Coupled_Iterative) {
                        if (phystime == 0) {
                            CreatePhysicalDataLogger();
                            CreateResidualLogger();
                        }
                        int iterationCounter = 0;
                        double hydroDynForceTorqueResidual = double.MaxValue;
                        foreach (Particle p in m_Particles) {
                            p.Motion.GetParticleDensity(p.particleDensity);
                        }
                        while (hydroDynForceTorqueResidual > HydrodynConvergenceCriterion) {
                            Console.WriteLine("Auxillary stuff");
                            Auxillary.CheckForMaxIterations(iterationCounter, ((FSI_Control)Control).maxIterationsFullyCoupled);
                            Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                            Auxillary.SaveOldParticleState(m_Particles, iterationCounter, ((FSI_Control)Control).hydrodynamicsConvergenceCriterion);

                            // actual physics
                            // -------------------------------------------------
                            if (IsFullyCoupled && iterationCounter == 0) {
                                Console.WriteLine("Init");
                                InitializeParticlePerIteration(m_Particles, TimestepInt);
                            }
                            else {
                                Console.WriteLine("bdf timestepper");
                                m_BDF_Timestepper.Solve(phystime, dt, false);
                                Console.WriteLine("Hydrodynamcs");
                                CalculateHydrodynamicForces(m_Particles, dt, !IsFullyCoupled);
                            }
                            Console.WriteLine("Particle Velocity");
                            CalculateParticleVelocity(m_Particles, dt, iterationCounter);

                            // not a fully coupled system? -> no iteration
                            // -------------------------------------------------
                            if (!IsFullyCoupled)
                                break;

                            // residual
                            // -------------------------------------------------
                            Console.WriteLine("Residual");
                            hydroDynForceTorqueResidual = Auxillary.CalculateParticleResidual(m_Particles, ref iterationCounter);

                            // print iteration status
                            // -------------------------------------------------
                            Console.WriteLine("Print1");
                            Auxillary.PrintResultToConsole(m_Particles, phystime, hydroDynForceTorqueResidual, iterationCounter);
                            Console.WriteLine("Log1");
                            LogResidual(phystime, iterationCounter, hydroDynForceTorqueResidual);
                        }

                        // collision
                        // -------------------------------------------------
                        Console.WriteLine("Collision");
                        CalculateCollision(m_Particles, cellColor, dt);

                        // particle position
                        // -------------------------------------------------
                        Console.WriteLine("Position");
                        CalculateParticlePosition(dt);

                        // print
                        // -------------------------------------------------
                        Console.WriteLine("Print2");
                        Auxillary.PrintResultToConsole(m_Particles, FluidViscosity, FluidDensity, phystime, TimestepInt, out double Test_RotationalVelocity, out Test_Force);
                        Console.WriteLine("Log2");
                        LogPhysicalData(phystime);

                        // Save for NUnit Test
                        // -------------------------------------------------
                        Console.WriteLine("NUnit");
                        SaveForNUnitTest(Test_RotationalVelocity);

                        // level set tracker 
                        // -------------------------------------------------
                        Console.WriteLine("LsTrk push");
                        if (IsFullyCoupled) {// in other branches, called by the BDF timestepper
                            LsTrk.IncreaseHistoryLength(1);
                            LsTrk.PushStacks();
                        }
                    }
                    else {// LevelSetHandling.Coupled_Iterative
                        foreach (Particle p in m_Particles) {
                            p.iteration_counter_P = -1;
                            p.forceAndTorque_convergence = ((FSI_Control)this.Control).hydrodynamicsConvergenceCriterion;
                        }
                        m_BDF_Timestepper.Solve(phystime, dt, false);
                    }
                }

                // finalize
                // ========
                if (LsTrk.PushCount - OldPushCount != 1) {
                    throw new ApplicationException("Illegal number of level-set push actions in time-step." + (LsTrk.PushCount - OldPushCount) + " It is important that LevelSetTracker.PushStacks() is called *exactly once per time-step*, at the beginning.");
                }
                ResLogger.NextTimestep(false);
                Console.WriteLine("Done with time-step " + TimestepInt + ".");
                return dt;
            }
        }

        /// <summary>
        /// Saves forces and rotational velocity for NUnit test.
        /// </summary>
        /// <param name="TestRotationalVelocity">
        /// </param>
        private void SaveForNUnitTest(double TestRotationalVelocity) {
            QueryHandler.ValueQuery("C_Drag", 2 * Test_Force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
            QueryHandler.ValueQuery("C_Lift", 2 * Test_Force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
            QueryHandler.ValueQuery("Angular_Velocity", TestRotationalVelocity, true); // (TestCase FlowRotationalCoupling)
        }

        /// <summary>
        /// Calls the calculation of the hydrodyn. forces and torque in the particle.cs
        /// </summary>
        /// <param name="particles">
        /// </param>
        /// /// <param name="dt">
        /// </param>
        /// /// <param name="firstIteration">
        /// </param>
        private void CalculateHydrodynamicForces(List<Particle> particles, double dt, bool firstIteration = true) {
            // Note on MPI parallelization of particle solver:
            // ===============================================
            // - hydrodynamic forces are computed for each domain and added together;
            //   e.g. in the case of particles at MPI-boundaries each processor computes his part of the integral
            // - collisions are detected globally / collision forces are thus only computed on MPI rank 0
            // - finally, the forces on all particles are summed over all MPI processors and the result is stored on all processors (MPI_Allreduce)
            // - particle motion is computed for all particles simultaneously on all processors; every processor knows every particle
            // - since the particle solver is much cheaper than the flow solver, this "not-really parallel" approach may work up to a few hundreds of particles
            // ===============================================
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            for (int p = 0; p < particles.Count(); p++) {
                Particle currentParticle = particles[p];
                if (firstIteration)
                    currentParticle.Motion.SaveHydrodynamicsOfPreviousTimestep();
                currentParticle.Motion.UpdateForcesAndTorque(Velocity, Pressure, LsTrk, currentParticle.CutCells_P(LsTrk), FluidViscosity, FluidDensity, firstIteration, dt);
            }
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
                Console.WriteLine("Predicting forces for the next timestep...");
                if (currentParticle.Motion.useAddedDamping) {
                    currentParticle.Motion.UpdateDampingTensors();
                }
                currentParticle.Motion.SaveHydrodynamicsOfPreviousTimestep();
                currentParticle.Motion.PredictForceAndTorque(currentParticle.activeStress, TimestepInt);
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
            foreach (Particle p in Particles) {
                if (IterationCounter == 0) {
                    p.Motion.SaveVelocityOfPreviousTimestep();
                }
                p.Motion.UpdateParticleVelocity(dt);
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
                logHydrodynamicsResidual.WriteLine(string.Format("{0}\t{1}\t{2}", "Time", "Iteration", "Residual"));
            }
        }

        /// <summary>
        /// Creates a log file for the residum of the hydrodynamic forces.
        /// </summary>
        private void LogResidual(double phystime, int iterationCounter, double residual) {
            if ((MPIRank == 0) && (logPhysicalDataParticles != null)) {
                logHydrodynamicsResidual.WriteLine(string.Format("{0}\t{1}\t{2}", phystime, iterationCounter, residual));
                logHydrodynamicsResidual.Flush();
            }
        }
        
        /// <summary>
        /// Creates a log file for the physical data of the particles. Only active if a database is specified.
        /// </summary>
        private void CreatePhysicalDataLogger() {
            if ((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                logPhysicalDataParticles = DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                logPhysicalDataParticles.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", "#Particle", "#Time", "Position X", "Position Y", "Angle", "Transl. Velocity X", "Transl. Velocity Y", "Rot. Velocity", "Force X", "Force Y", "Angular Momentum"));
            }
        }

        /// <summary>
        /// Writes the physical data of the particles to a log file.
        /// </summary>
        /// <param name = phystime>
        /// </param>
        private void LogPhysicalData(double phystime) {
            if ((MPIRank == 0) && (logPhysicalDataParticles != null)) {
                for (int p = 0; p < m_Particles.Count(); p++) {
                    logPhysicalDataParticles.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", p, phystime, m_Particles[p].Motion.position[0][0], m_Particles[p].Motion.position[0][1], m_Particles[p].Motion.angle[0], m_Particles[p].Motion.translationalVelocity[0][0], m_Particles[p].Motion.translationalVelocity[0][1], m_Particles[p].Motion.rotationalVelocity[0], m_Particles[p].Motion.hydrodynamicForces[0][0], m_Particles[p].Motion.hydrodynamicForces[0][1], m_Particles[p].Motion.hydrodynamicTorque[0]));
                    logPhysicalDataParticles.Flush();
                }
            }
        }

        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {
            var tsi = new FSI_TimestepInfo(t, CurrentSessionInfo, timestepno, IOFields, m_Particles);
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
                Debug.Assert(GenericBlas.L2Dist(b.Particles[l].Motion.position[0], o.Particles[l].Motion.position[0]) < 1e-13);
            }

        }

        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="m_Particles"/>) to the database
        /// </summary>
        protected override TimestepNumber RestartFromDatabase(out double time) {

            // this sux, because the database API is totally fucked up
            var db = GetDatabase();
            Guid Rst_Tsid = base.GetRestartTimestepID();
            Guid Rst_SessionId = Control.RestartInfo.Item1;
            ISessionInfo session = db.Controller.GetSessionInfo(Rst_SessionId);

            var ArschInfo = ((DatabaseDriver)(DatabaseDriver)).LoadTimestepInfo<FSI_TimestepInfo>(Rst_Tsid, session, db);

            // init particles
            m_Particles = ArschInfo.Particles.ToList();
            UpdateLevelSetParticles(ArschInfo.PhysicalTime);

            // call base shit
            var R = base.RestartFromDatabase(out time);

            // return
            return R;
        }


        /// <summary>
        /// For restarting calculations, its important to reload old solutions if one uses a higher order method in time
        /// </summary>
        /// <param name="time"></param>
        /// <param name="timestep"></param>
        public override void PostRestart(double time, TimestepNumber timestep) {
            //var fsDriver = this.DatabaseDriver.FsDriver;
            //string pathToOldSessionDir = System.IO.Path.Combine(
            //    fsDriver.BasePath, "sessions", this.CurrentSessionInfo.RestartedFrom.ToString());
            //string pathToPhysicalData = System.IO.Path.Combine(pathToOldSessionDir,"PhysicalData.txt");
            //string[] records = File.ReadAllLines(pathToPhysicalData); 

            //string line1 = File.ReadLines(pathToPhysicalData).Skip(1).Take(1).First();
            //string line2 = File.ReadLines(pathToPhysicalData).Skip(2).Take(1).First();
            //string[] fields_line1 = line1.Split('\t');
            //string[] fields_line2 = line2.Split('\t');

            //Console.WriteLine("Line 1 " + fields_line1);

            //double dt = Convert.ToDouble(fields_line2[1]) - Convert.ToDouble(fields_line1[1]);

            //int idx_restartLine = Convert.ToInt32(time/dt + 1.0);
            //string restartLine = File.ReadLines(pathToPhysicalData).Skip(idx_restartLine-1).Take(1).First();
            //double[] values = Array.ConvertAll<string, double>(restartLine.Split('\t'), double.Parse);

            //if (time == values[1]+dt)
            //{
            //    Console.WriteLine("Restarting from time " + values[1]);
            //}

            //oldPosition[0] = values[7];
            //oldPosition[1] = values[8];
            //newTransVelocity[0] = values[4];
            //newTransVelocity[1] = values[5];
            //oldTransVelocity[0] = 0;
            //oldTransVelocity[1] = 0;
            //TransVelocityN2[0] = 0;
            //TransVelocityN2[1] = 0;
            //TransVelocityN3[0] = 0;
            //TransVelocityN3[1] = 0;
            //TransVelocityN4[0] = 0;
            //TransVelocityN4[1] = 0;
            //force[0] = values[2];
            //force[1] = values[3];

            //if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty))
            //{
            //    Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
            //    string firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#Timestep", "#Time", "DragForce", "LiftForce", "VelocityX", "VelocityY", "AngularVelocity", "xPosition", "yPosition", "ParticleRe");
            //    Log_DragAndLift.WriteLine(firstline);
            //    Log_DragAndLift.WriteLine(restartLine);
            //}
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
        private void CalculateCollision(List<Particle> Particles, int[] CellColor, double dt) {
            if (((FSI_Control)Control).collisionModel == FSI_Control.CollisionModel.NoCollisionModel)
                return;
            if (((FSI_Control)Control).collisionModel == FSI_Control.CollisionModel.RepulsiveForce)
                throw new NotImplementedException("Repulsive force model is currently unsupported, please use the momentum conservation model.");

            foreach (Particle p in Particles) {
                p.isCollided = false;
            }

            // Only particles with the same colour a close to each other, thus, we only test for collisions within those particles.
            // Determine colour.
            // =================================================
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate(LsTrk);
            int[] _GlobalParticleColor = globalParticleColor.CloneAs();
            for (int i = 0; i < _GlobalParticleColor.Length; i++) {
                int CurrentColor = _GlobalParticleColor[i];
                if (CurrentColor == 0)
                    continue;
                int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesOneColor(_GlobalParticleColor, CurrentColor);

                // Multiple particles with the same colour, trigger collision detection
                // =================================================
                if (ParticlesOfCurrentColor.Length >= 1 && CurrentColor != 0) {
                    List<Particle> currentParticles = levelSetUpdate.GetParticleListOneColor(m_Particles, _GlobalParticleColor, CurrentColor);
                    FSI_Collision _Collision = new FSI_Collision(LsTrk, CurrentColor, FluidViscosity, FluidDensity, ((FSI_Control)Control).CoefficientOfRestitution, dt, LsTrk.GridDat.Cells.h_minGlobal);
                    _Collision.CalculateCollision(currentParticles, GridData, CellColor);
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
            isCollidedSend[0] = currentParticle.isCollided ? 1 : 0;
            double[] isCollidedReceive = new double[MPISize];
            MPISendAndReceive(isCollidedSend, ref isCollidedReceive);

            bool noCurrentCollision = true;
            for (int i = 0; i < isCollidedReceive.Length; i++) {
                // The particle is collided, thus, copy the data from
                // the owning process.
                // ===================================================
                if (isCollidedReceive[i] != 0) {
                    int noOfVars = 13;
                    double[] dataSend = new double[noOfVars];
                    dataSend[0] = currentParticle.Motion.rotationalVelocity[0];
                    dataSend[1] = currentParticle.Motion.translationalVelocity[0][0];
                    dataSend[2] = currentParticle.Motion.translationalVelocity[0][1];
                    dataSend[3] = currentParticle.Motion.angle[0];
                    dataSend[4] = currentParticle.Motion.position[0][0];
                    dataSend[5] = currentParticle.Motion.position[0][1];
                    dataSend[6] = currentParticle.Motion.collisionTimestep;
                    dataSend[7] = currentParticle.Motion.rotationalVelocity[1];
                    dataSend[8] = currentParticle.Motion.translationalVelocity[1][0];
                    dataSend[9] = currentParticle.Motion.translationalVelocity[1][1];
                    dataSend[10] = currentParticle.Motion.angle[1];
                    dataSend[11] = currentParticle.Motion.position[1][0];
                    dataSend[12] = currentParticle.Motion.position[1][1];

                    double[] dataReceive = new double[noOfVars * MPISize];
                    MPISendAndReceive(dataSend, ref dataReceive);

                    currentParticle.Motion.rotationalVelocity[0] = dataReceive[0 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[0][0] = dataReceive[1 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[0][1] = dataReceive[2 + i * noOfVars];
                    currentParticle.Motion.angle[0] = dataReceive[3 + i * noOfVars];
                    currentParticle.Motion.position[0][0] = dataReceive[4 + i * noOfVars];
                    currentParticle.Motion.position[0][1] = dataReceive[5 + i * noOfVars];
                    currentParticle.Motion.collisionTimestep = dataReceive[6 + i * noOfVars];
                    currentParticle.Motion.rotationalVelocity[1] = dataReceive[7 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[1][0] = dataReceive[8 + i * noOfVars];
                    currentParticle.Motion.translationalVelocity[1][1] = dataReceive[9 + i * noOfVars];
                    currentParticle.Motion.angle[1] = dataReceive[10 + i * noOfVars];
                    currentParticle.Motion.position[1][0] = dataReceive[11 + i * noOfVars];
                    currentParticle.Motion.position[1][1] = dataReceive[12 + i * noOfVars];

                    currentParticle.isCollided = true;
                    noCurrentCollision = false;
                }
            }
            // nothing happend
            // ===================================================
            if (noCurrentCollision) {
                currentParticle.isCollided = false;
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
            if (((FSI_Control)Control).AdaptiveMeshRefinement) {
                bool AnyChangeInGrid = GridRefinementController.ComputeGridChange((GridData)GridData, null, GetCellMaskWithRefinementLevels(), out List<int> CellsToRefineList, out List<int[]> Coarsening);
                if (AnyChangeInGrid) {
                    int[] consoleRefineCoarse = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                    int oldJ = this.GridData.CellPartitioning.TotalLength;
                    Console.WriteLine("       Refining " + consoleRefineCoarse[0] + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + consoleRefineCoarse[1] + " of " + oldJ + " cells");
                    newGrid = ((GridData)GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);
                }
                else {
                    newGrid = null;
                    old2NewGrid = null;
                }
            }
            else {
                newGrid = null;
                old2NewGrid = null;
            }
        }

        /// <summary>
        /// Creates the cellmask which should be refined.
        /// </summary>
        private List<Tuple<int, CellMask>> GetCellMaskWithRefinementLevels() {
            int refinementLevel = ((FSI_Control)Control).RefinementLevel;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            BitArray superFineCells = new BitArray(noOfLocalCells);
            double superFineRadius = 4 * LsTrk.GridDat.Cells.h_minGlobal;
            BitArray mediumFineCells = new BitArray(noOfLocalCells);
            double mediumFineRadius = LsTrk.GridDat.Cells.h_maxGlobal / 2;
            for (int p = 0; p < m_Particles.Count; p++) {
                Particle particle = m_Particles[p];
                for (int j = 0; j < noOfLocalCells; j++) {
                    double[] centerPoint = new double[] { CellCenters[j, 0], CellCenters[j, 1] };
                    if (!superFineCells[j])
                        superFineCells[j] = particle.Contains(centerPoint, superFineRadius);
                    if (!mediumFineCells[j])
                        mediumFineCells[j] = particle.Contains(centerPoint, mediumFineRadius);
                }
            }
            int coarseRefinementLevel = refinementLevel > 2 ? refinementLevel / 2 : 1;
            if (refinementLevel - coarseRefinementLevel > coarseRefinementLevel)
                coarseRefinementLevel += 1;
            List<Tuple<int, CellMask>> AllCellsWithMaxRefineLevel = new List<Tuple<int, CellMask>> {
                new Tuple<int, CellMask>(refinementLevel, new CellMask(GridData, superFineCells)),
                new Tuple<int, CellMask>(coarseRefinementLevel, new CellMask(GridData, mediumFineCells))
            };
            return AllCellsWithMaxRefineLevel;
        }
    }
}