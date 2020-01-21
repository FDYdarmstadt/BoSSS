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
            double[] couplingArray = new double[X.Dim + 6];
            FSI_ParameterAtIB couplingParameters;
            foreach (Particle p in m_Particles) {
                bool containsParticle = m_Particles.Count == 1 ? true : p.Contains(X, GridData.iGeomCells.h_min.Min());
                if (containsParticle) {
                    couplingParameters = new FSI_ParameterAtIB(p, X);
                    p.CalculateRadialVector(X, out Vector RadialVector, out double radialLength);
                    couplingArray[0] = p.Motion.GetTranslationalVelocity(0)[0];
                    couplingArray[1] = p.Motion.GetTranslationalVelocity(0)[1];
                    couplingArray[2] = p.Motion.GetRotationalVelocity(0);
                    couplingArray[3] = RadialVector[0];
                    couplingArray[4] = RadialVector[1];
                    couplingArray[5] = radialLength;
                    couplingArray[6] = p.ActiveStress; // zero for passive particles
                    couplingArray[7] = p.Motion.GetAngle(0);
                    return couplingParameters;
                }
            }
            return null;
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
                    break;

                case LevelSetHandling.LieSplitting:
                    UpdateLevelSetParticles(phystime);
                    break;

                case LevelSetHandling.FSI_LieSplittingFullyCoupled:
                    UpdateLevelSetParticles(phystime);
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
                Motion_AllParticles AllParticleHydrodynamics = new Motion_AllParticles(LsTrk);
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
            FSI_LevelSetUpdate levelSetUpdate = new FSI_LevelSetUpdate(LsTrk);
            // Step 1
            // Define an array with the respective cell colors
            // =======================================================
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            ExterminateExcessGhosts();
            CreateGhostsAtPeriodicBoundary();
            GhostToMaster();
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

                    // Get particle level set
                    // -----------------------------
                    double levelSetFunction(double[] X, double t) {
                        // Generating the correct sign
                        double levelSetFunctionOneColor = Math.Pow(-1, particlesOfCurrentColor.Length - 1);
                        // Multiplication over all particle-level-sets within the current color
                        for (int pC = 0; pC < particlesOfCurrentColor.Length; pC++) {
                            Particle currentParticle = m_Particles[particlesOfCurrentColor[pC]];
                            double tempLevelSetFunction = currentParticle.LevelSetFunction(X);
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
            int[] coloredCells = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int[] coloredCellsExchange = coloredCells.CloneAs();

            // No particles on current proc
            // -----------------------------
            if (m_Particles.Count == 0)
                return coloredCells;
            coloredCellsExchange.MPIExchange(GridData);

            // Step 2
            // Color neighbour cells
            // =======================================================
            int neighbourSearchDepth = 1;
            for (int k = 0; k < neighbourSearchDepth; k++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (coloredCellsExchange[CellNeighbors[i]] != 0 && coloredCellsExchange[j] == 0) {
                            coloredCells[j] = coloredCellsExchange[CellNeighbors[i]];
                        }
                    }
                }
                coloredCellsExchange = coloredCells.CloneAs();
                coloredCellsExchange.MPIExchange(GridData);
            }

            // Step 4
            // Find neighbouring colours and recolour one of them
            // =======================================================
            // Find neighbouring cells with
            // different colours
            // -----------------------------
            int maxColor = coloredCells.Max().MPIMax();
            int[,] colorToRecolorWith = new int[maxColor + 1, 2];
            for (int j = 0; j < noOfLocalCells; j++) {
                if (coloredCells[j] != 0) {
                    GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (coloredCellsExchange[CellNeighbors[i]] != coloredCells[j] && coloredCellsExchange[CellNeighbors[i]] > 0) {
                            if (coloredCellsExchange[CellNeighbors[i]] < coloredCells[j] || colorToRecolorWith[coloredCells[j], 1] > coloredCellsExchange[CellNeighbors[i]]) {
                                colorToRecolorWith[coloredCells[j], 0] = coloredCells[j];
                                colorToRecolorWith[coloredCells[j], 1] = coloredCellsExchange[CellNeighbors[i]];
                            }
                            if (coloredCellsExchange[CellNeighbors[i]] > coloredCells[j]) {
                                if (colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]], 0] == 0 || colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]], 1] > coloredCells[j]) {
                                    colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]], 0] = coloredCellsExchange[CellNeighbors[i]];
                                    colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]], 1] = coloredCells[j];
                                }
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
                        if (coloredCells[j] == colorToRecolorWith[i, 0]) {
                            coloredCells[j] = colorToRecolorWith[i, 1];
                        }
                    }
                }
            }

            // Set DG-Fields
            // -----------------------------
            for (int j = 0; j < noOfLocalCells; j++) {
                ParticleColor.SetMeanValue(j, coloredCells[j]);
                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(regionsCode[j], 0));
            }
            return coloredCells;
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
                    Vector center = new Vector(GridData.iLogicalCells.GetCenter(j));
                    if (currentParticle.Contains(center, h_min)) {
                        ParticleColor.SetMeanValue(j, p + 1);
                        coloredCells.Add(j);
                        cells[j] = p + 1;
                        cellsExchange[j] = cells[j];
                    }
                }
            }
            cellsExchange.MPIExchange(GridData);
            FixNeighbourColoring(cellsExchange);
            for (int j = 0; j < J; j++) {
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

        private void CreateGhostsAtPeriodicBoundary() {
            int spatialDim = 2;
            bool[] isPeriodic = ((FSI_Control)Control).BoundaryIsPeriodic;
            double[][] boundaryCoordinates = ((FSI_Control)Control).BoundaryPositionPerDimension;
            Vector[] ghostPositions = new Vector[3];
            for (int p = 0; p < m_Particles.Count(); p++) {
                Particle currentParticle = m_Particles[p];
                List<Particle> ghostParticles = new List<Particle>();
                int[] ghostHierachy = currentParticle.MasterGhostIDs.CloneAs();
                if (!currentParticle.IsMaster)
                    continue;
                Vector particlePosition = currentParticle.Motion.GetPosition();
                int idOffset = 0;
                for (int d1 = 0; d1 < spatialDim; d1++) { // which direction?
                    if (!isPeriodic[d1])
                        continue;
                    for(int wallID = 0; wallID < spatialDim; wallID++) { // which wall?
                        if (PeriodicOverlap(currentParticle, d1, wallID)) {
                            ghostHierachy[0] = p + 1;
                            Vector originNeighbouringDomain;
                            if (d1 == 0)
                                originNeighbouringDomain = new Vector(2 * boundaryCoordinates[0][1 - wallID], 0);
                            else
                                originNeighbouringDomain = new Vector(0, 2 * boundaryCoordinates[1][1 - wallID]);
                            Particle ghostParticle;
                            int ghostID = m_Particles.Count() + idOffset + 1;
                            if (ghostHierachy[d1 + 1] == 0) {
                                ghostPositions[d1] = originNeighbouringDomain + particlePosition;
                                ghostHierachy[d1 + 1] = ghostID;
                                ghostParticle = currentParticle.CloneAs();
                                ghostParticle.SetGhost();
                                ghostParticle.Motion.SetGhostPosition(ghostPositions[d1]);
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
                                    if (PeriodicOverlap(ghostParticle, 1, wallID2)) {
                                        originNeighbouringDomain = new Vector(0, 2 * boundaryCoordinates[1][1 - wallID2]);
                                        ghostPositions[2] = originNeighbouringDomain + ghostPositions[0];
                                        idOffset += 1;
                                        ghostID = m_Particles.Count() + d1 + idOffset;
                                        ghostHierachy[3] = ghostID;
                                        ghostParticle = currentParticle.CloneAs();
                                        ghostParticle.SetGhost();
                                        ghostParticle.Motion.SetGhostPosition(ghostPositions[2]);
                                        ghostParticles.Add(ghostParticle.CloneAs());
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if(ghostParticles.Count >= 1) {
                    currentParticle.SetGhostHierachy(ghostHierachy);
                    for (int p2 = 0; p2 < ghostParticles.Count(); p2++) {
                        ghostParticles[p2].SetGhostHierachy(ghostHierachy);
                    }
                    m_Particles.AddRange(ghostParticles);
                }
            }
        }


        private void GhostToMaster() {
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
                                m_Particles[ghostHierachy[i] - 1].MasterGhostIDs = newGhostHierachy.CloneAs();
                            }
                            return;
                        }
                    }
                }
            }
        }

        private void ExterminateExcessGhosts() {
            int spatialDim = 2;
            for (int p = 0; p < m_Particles.Count(); p++) {
                for (int d = 0; d < spatialDim; d++) {
                    for (int wallID = 0; wallID < spatialDim; wallID++) {
                        Particle currentParticle = m_Particles[p];
                        Vector particlePosition = currentParticle.Motion.GetPosition();
                        double[][] boundaryCoordinates = ((FSI_Control)Control).BoundaryPositionPerDimension;
                        double distance = particlePosition[d] - boundaryCoordinates[d][wallID];
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
            double[][] boundaryCoordinates = ((FSI_Control)Control).BoundaryPositionPerDimension;
            double distance = particlePosition[d1] - boundaryCoordinates[d1][d2];
            double particleMaxLength = currentParticle.GetLengthScales().Max();
            if (Math.Abs(distance) < particleMaxLength) {
                if (d1 == 0)
                    currentParticle.ClosestPointOnOtherObjectToThis = new Vector(boundaryCoordinates[d1][d2], particlePosition[1]);
                else
                    currentParticle.ClosestPointOnOtherObjectToThis = new Vector(particlePosition[0], boundaryCoordinates[d1][d2]);
                FSI_Collision periodicCollision = new FSI_Collision(LsTrk, 0, 0, 0);
                periodicCollision.CalculateMinimumDistance(currentParticle, out _, out Vector _, out Vector _, out bool Overlapping);
                return Overlapping;
            }
            return false;
        }

        private bool AnyOverlap(Particle currentParticle) {
            int spatialDim = 2;
            for (int d = 0; d < spatialDim; d++) {
                for (int wallID = 0; wallID < spatialDim; wallID++) {
                    Vector particlePosition = currentParticle.Motion.GetPosition();
                    double[][] boundaryCoordinates = ((FSI_Control)Control).BoundaryPositionPerDimension;
                    if (d == 0)
                        currentParticle.ClosestPointOnOtherObjectToThis = new Vector(boundaryCoordinates[d][wallID], particlePosition[1]);
                    else
                        currentParticle.ClosestPointOnOtherObjectToThis = new Vector(particlePosition[0], boundaryCoordinates[d][wallID]);
                    FSI_Collision periodicCollision = new FSI_Collision(LsTrk, 0, 0, 0);
                    periodicCollision.CalculateMinimumDistance(currentParticle, out _, out Vector _, out Vector _, out bool Overlapping);
                    if (Overlapping)
                        return true;
                }
            }
            return false;
        }

        private bool IsInsideOfDomain(Particle currentParticle) {
            int spatialDim = 2;
            bool[] isPeriodic = ((FSI_Control)Control).BoundaryIsPeriodic;
            double[][] boundaryCoordinates = ((FSI_Control)Control).BoundaryPositionPerDimension;
            Vector position = currentParticle.Motion.GetPosition();
            for (int d = 0; d < spatialDim; d++) {
                if (!isPeriodic[d])
                    continue;
                for (int wallID = 0; wallID < spatialDim; wallID++) {
                    Vector wallNormal = new Vector(Math.Sign(boundaryCoordinates[d][wallID]) * (1 - d), Math.Sign(boundaryCoordinates[d][wallID]) * d);
                    Vector wallToPoint = d == 0
                        ? new Vector(position[0] - boundaryCoordinates[0][wallID], position[1])
                        : new Vector(position[0], position[1] - boundaryCoordinates[1][wallID]);
                    if (wallNormal * wallToPoint > 0)
                        return false;
                }
            }
            return true;
        }
        
        bool initAddedDamping = true;
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
                    if (phystime == 0) { CreatePhysicalDataLogger(); }
                    LsTrk.PushStacks(); // in other branches, called by the BDF timestepper
                    DGLevSet.Push();

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
                    Auxillary.PrintResultToConsole(m_Particles, 0, 0, phystime, TimestepInt, ((FSI_Control)Control).FluidDomainVolume, out double MPIangularVelocity, out Test_Force);
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
                        int minIteration = 3;
                        Motion_AllParticles AllParticleHydrodynamics = new Motion_AllParticles(LsTrk);
                        while (hydroDynForceTorqueResidual > HydrodynConvergenceCriterion || iterationCounter < minIteration) {
                            if (iterationCounter > ((FSI_Control)Control).maxIterationsFullyCoupled)
                                throw new ApplicationException("No convergence in coupled iterative solver, number of iterations: " + iterationCounter);
                            Auxillary.ParticleState_MPICheck(m_Particles, GridData, MPISize);
                            AllParticleHydrodynamics.SaveHydrodynamicOfPreviousIteration(m_Particles);
                            if (IsFullyCoupled && iterationCounter == 0) {
                                InitializeParticlePerIteration(m_Particles, TimestepInt);
                            }
                            else {
                                VectorField<SinglePhaseField> velocityOld = Velocity.CloneAs();
                                SinglePhaseField pressureOld = Pressure.CloneAs();
                                m_BDF_Timestepper.Solve(phystime, dt, false);
                                ParticleHydrodynamicsIntegration hydrodynamicsIntegration = new ParticleHydrodynamicsIntegration(2, Velocity, Pressure, LsTrk, FluidViscosity);
                                AllParticleHydrodynamics.CalculateHydrodynamics(m_Particles, hydrodynamicsIntegration, FluidDensity, false);
                                if (iterationCounter != 1) 
                                {
                                    double underrelax = 0.1;
                                    Velocity.Scale(underrelax);
                                    Velocity.Acc((1 - underrelax), velocityOld);
                                    Pressure.Scale(underrelax);
                                    Pressure.Acc((1 - underrelax), pressureOld);
                                }
                            }
                            if (TimestepInt != 1 || iterationCounter != 0)
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
                            Auxillary.PrintResultToConsole(phystime, hydroDynForceTorqueResidual, iterationCounter);
                            LogResidual(phystime, iterationCounter, hydroDynForceTorqueResidual);
                        }
                        if (TimestepInt == 1 || IsMultiple(TimestepInt, 10)) {
                            for (int p = 0; p < m_Particles.Count(); p++) {
                                m_Particles[p].Motion.CreateStressLogger(CurrentSessionInfo, DatabaseDriver, phystime, p);
                                m_Particles[p].Motion.LogStress(phystime);
                            }
                        }
                        // collision
                        // -------------------------------------------------
                        CalculateCollision(m_Particles, cellColor, dt);

                        // particle position
                        // -------------------------------------------------
                        CalculateParticlePosition(dt);

                        // print
                        // -------------------------------------------------
                        Auxillary.PrintResultToConsole(m_Particles, FluidViscosity, FluidDensity, phystime, TimestepInt, ((FSI_Control)Control).FluidDomainVolume, out double Test_RotationalVelocity, out Test_Force);
                        LogPhysicalData(phystime);

                        // Save for NUnit Test
                        // -------------------------------------------------
                        SaveForNUnitTest(Test_RotationalVelocity);

                        // level set tracker 
                        // -------------------------------------------------
                        if (IsFullyCoupled) {// in other branches, called by the BDF timestepper
                            LsTrk.IncreaseHistoryLength(1);
                            LsTrk.PushStacks();
                        }
                    }
                    else {// LevelSetHandling.Coupled_Iterative
                        m_BDF_Timestepper.Solve(phystime, dt, false);
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

        public bool IsMultiple(int a, int b) {
            return (a % b) == 0;
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
            Motion_AllParticles AllParticleHydrodynamics = new Motion_AllParticles(LsTrk);
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
                Console.WriteLine("Predicting forces for the next timestep...");
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
            foreach (Particle p in Particles) {
                if (p.Motion.IsGhost)
                    continue;
                if (IterationCounter == 0) {
                    p.Motion.SaveVelocityOfPreviousTimestep();
                }
                p.Motion.UpdateParticleVelocity(dt);
                if (p.Motion.GetHasGhost()) {
                    Particle ghost = Particles[p.Motion.GetGhostID()-1];
                    ghost.Motion.CopyNewVelocity(p.Motion.GetTranslationalVelocity(), p.Motion.GetRotationalVelocity());
                    ghost.Motion.UpdateParticleVelocity(dt);
                }
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
            IDatabaseDriver test = DatabaseDriver;
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
                logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}", "particle", "time", "posX", "posY", "angle", "velX", "velY", "rot", "fX", "fY", "T"));
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
                    logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}", p, phystime, m_Particles[p].Motion.GetPosition(0)[0], m_Particles[p].Motion.GetPosition(0)[1], m_Particles[p].Motion.GetAngle(0), m_Particles[p].Motion.GetTranslationalVelocity(0)[0], m_Particles[p].Motion.GetTranslationalVelocity(0)[1], m_Particles[p].Motion.GetRotationalVelocity(0), m_Particles[p].Motion.GetHydrodynamicForces(0)[0], m_Particles[p].Motion.GetHydrodynamicForces(0)[1], m_Particles[p].Motion.GetHydrodynamicTorque(0)));
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
                Debug.Assert(GenericBlas.L2Dist((double[])b.Particles[l].Motion.GetPosition(0), (double[])o.Particles[l].Motion.GetPosition(0)) < 1e-13);
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
                p.IsCollided = false;
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
                    Particle[] currentParticles = new Particle[ParticlesOfCurrentColor.Length];
                    for (int j = 0; j < ParticlesOfCurrentColor.Length; j++) {
                        currentParticles[j] = m_Particles[ParticlesOfCurrentColor[j]];
                    }
                    FSI_Collision _Collision = new FSI_Collision(LsTrk, CurrentColor, ((FSI_Control)Control).CoefficientOfRestitution, dt);
                    _Collision.CalculateCollision(currentParticles, GridData, CellColor, LsTrk);
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
            if (((FSI_Control)Control).AdaptiveMeshRefinement) {
                bool AnyChangeInGrid = GridRefinementController.ComputeGridChange((GridData)GridData, null, GetCellMaskWithRefinementLevels(), out List<int> CellsToRefineList, out List<int[]> Coarsening);
                if (AnyChangeInGrid) {
                    int[] consoleRefineCoarse = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                    int oldJ = this.GridData.CellPartitioning.TotalLength;
                    Console.WriteLine("       Refining " + consoleRefineCoarse[0] + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + consoleRefineCoarse[1] + " of " + oldJ + " cells");
                    Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Count().MPISum());
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

        double h_maxStart = 0;

        /// <summary>
        /// Creates the cellmask which should be refined.
        /// </summary>
        private List<Tuple<int, BitArray>> GetCellMaskWithRefinementLevels() {
            h_maxStart = h_maxStart == 0 ? LsTrk.GridDat.Cells.h_maxGlobal : h_maxStart;
            int refinementLevel = ((FSI_Control)Control).RefinementLevel;
            double h_minStart = h_maxStart / (2 * refinementLevel);
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            BitArray coarseCells = new BitArray(noOfLocalCells);
            BitArray mediumCells = new BitArray(noOfLocalCells);
            BitArray fineCells = new BitArray(noOfLocalCells);
            double radiusCoarseCells = 2 * LsTrk.GridDat.Cells.h_maxGlobal;
            double radiusMediumCells = LsTrk.GridDat.Cells.h_maxGlobal;
            double radiusFineCells = 4 * LsTrk.GridDat.Cells.h_minGlobal;
            for (int p = 0; p < m_Particles.Count; p++) {
                Particle particle = m_Particles[p];
                for (int j = 0; j < noOfLocalCells; j++) {
                    Vector centerPoint = new Vector(CellCenters[j, 0], CellCenters[j, 1]);
                    //if (!coarseCells[j] && LsTrk.Regions.IsSpeciesPresentInCell(LsTrk.GetSpeciesId("A"), j)) {
                    //    coarseCells[j] = particle.Contains(centerPoint, radiusCoarseCells);
                    //}
                    //if (!mediumCells[j] && LsTrk.Regions.IsSpeciesPresentInCell(LsTrk.GetSpeciesId("A"), j)) {
                    //    mediumCells[j] = particle.Contains(centerPoint, radiusMediumCells);
                    //}
                    //if (!fineCells[j] && LsTrk.Regions.IsSpeciesPresentInCell(LsTrk.GetSpeciesId("A"), j)){
                    //    fineCells[j] = particle.Contains(centerPoint, radiusFineCells);
                    //}
                    if (!coarseCells[j]) {
                        coarseCells[j] = particle.Contains(centerPoint, radiusCoarseCells);
                    }
                    if (!mediumCells[j]) {
                        mediumCells[j] = particle.Contains(centerPoint, radiusMediumCells);
                    }
                    if (!fineCells[j]) {
                        fineCells[j] = particle.Contains(centerPoint, radiusFineCells);
                    }
                }
            }
            int medioumRefinementLevel = refinementLevel > 2 ? refinementLevel / 2 + 1 : 1;
            int coarseRefinementLevel = refinementLevel > 4 ? refinementLevel / 4 : 1;
            List<Tuple<int, BitArray>> AllCellsWithMaxRefineLevel = new List<Tuple<int, BitArray>> {
                new Tuple<int, BitArray>(refinementLevel, fineCells),
                new Tuple<int, BitArray>(medioumRefinementLevel, mediumCells),
                new Tuple<int, BitArray>(coarseRefinementLevel, coarseCells),
            };
            return AllCellsWithMaxRefineLevel;
        }
    }
}