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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using FSI_Solver;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    /// <summary>
    /// Fluid Structure Interaction (FSI) Solver, features:
    /// - Fluid-Particle interactions
    /// - Incompressible fluid
    /// - Rigid particles
    /// - Arbitrary geometries
    /// - Active particles
    /// </summary>
    /// <remarks>
    /// Current maintainer: B. Deußen (deussen@fdy.tu-darmstadt.de), F. Kummer
    /// </remarks>
    public class FSI_SolverMain : IBM_Solver.IBM_SolverMain {

        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {
            _Main(args, false, delegate () {
                FSI_SolverMain p = new FSI_SolverMain();
                return p;
            });
        }

        /// <summary>
        /// Set the initial state of the simulation.
        /// </summary>
        protected override void SetInitial() {
            if (((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.None) {
                throw new NotImplementedException("Currently not implemented for fixed motion");
            }

            ParticleList = ((FSI_Control)this.Control).Particles;
            if (ParticleList.IsNullOrEmpty())
                throw new Exception("Define at least on particle");

            UpdateLevelSetParticles(phystime: 0.0);
            CreatePhysicalDataLogger();
            base.SetInitial();
        }

        /// <summary>
        /// Spatial dimension
        /// </summary>
        /// <remarks>
        /// Currently hard-coded for 2D, because collision and particle solver is only implemented for two dimensions.
        /// </remarks>
        private readonly int spatialDim = 2;

        /// <summary>
        /// Iteration counter for iteration between the particle and fluid solver.
        /// </summary>
        private int iterationCounter = 0;

        /// <summary>
        /// If a dynamic time-step is used, this is the time-step size of the previous one.
        /// </summary>
        [DataMember]
        private double oldTimestep;

        /// <summary>
        /// If added damping is used, the respective tensors need to be set only in the first time-step.
        /// </summary>
        [DataMember]
        private bool initAddedDamping = true;

        /// <summary>
        /// A list for all particles
        /// </summary>
        private List<Particle> ParticleList;

        /// <summary>
        /// External access to particle list.
        /// </summary>
        public IList<Particle> Particles => ParticleList;

        /// <summary>
        /// An object with some additional methods, e.g. NaN and infinity exceptions, writing to console, MPI state check
        /// </summary>
        readonly private FSIAuxillary Auxillary = new FSIAuxillary();

        /// <summary>
        /// Methods dealing with coloring and level set
        /// </summary>
        private FSILevelSetUpdate LevelSetUpdate;

        /// <summary>
        /// Particle color field. The level set is only defined on colored cells.
        /// </summary>
        private SinglePhaseField ParticleColor;

        /// <summary>
        /// Level set distance field. 
        /// </summary>
        private SinglePhaseField LevelSetDistance;

        /// <summary>
        /// Level set distance field. 
        /// </summary>
        private SinglePhaseField CellID;

        /// <summary>
        /// Create all fields. Additional fields for the FSI-Solver are the particle color field and the level-set distance field.
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();

            ParticleColor = new SinglePhaseField(new Basis(GridData, 0), "ParticleColor");
            m_RegisteredFields.Add(ParticleColor);
            m_IOFields.Add(ParticleColor);

            LevelSetDistance = new SinglePhaseField(new Basis(GridData, 0), "LevelSetDistance");
            m_RegisteredFields.Add(LevelSetDistance);
            m_IOFields.Add(LevelSetDistance);

            CellID = new SinglePhaseField(new Basis(GridData, 0), "CellID");
            m_RegisteredFields.Add(CellID);
            m_IOFields.Add(CellID);
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
                    case LevelSetHandling.FSILieSplittingFullyCoupled:
                    return false;

                    default:
                    throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
                }
            }
        }

        /// <summary>
        /// BDF order.
        /// </summary>
        private int BDFOrder = 0;

        /// <summary>
        /// Array of all local cells with their specific color.
        /// </summary>
        [DataMember]
        private int[] CellColor = null;

        /// <summary>
        /// MPI-global array of all particles with their specific color.
        /// </summary>
        [DataMember]
        private int[] GlobalParticleColor = null;

        /// <summary>
        /// Fully coupled LieSplitting?
        /// </summary>
        [DataMember]
        private bool IsFullyCoupled => ((FSI_Control)Control).Timestepper_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled;

        /// <summary>
        /// Using static time-step or dynamic?
        /// </summary>
        [DataMember]
        private bool StaticTimestep => ((FSI_Control)Control).staticTimestep;

        /// <summary>
        /// The maximum time-step set in the control file.
        /// </summary>
        [DataMember]
        private double DtMax => ((FSI_Control)Control).dtMax;

        /// <summary>
        /// The starting time-step set in the control file. Only used in case of a dynamic time-step.
        /// </summary>
        [DataMember]
        private double StartDt => ((FSI_Control)Control).dtMax / 10;

        /// <summary>
        /// MaxGridLength
        /// </summary>
        [DataMember]
        private double MaxGridLength => ((FSI_Control)Control).MaxGridLength;

        /// <summary>
        /// RefinementLevel
        /// </summary>
        [DataMember]
        private double RefinementLevel => ((FSI_Control)Control).RefinementLevel;

        /// <summary>
        /// Calculates the minimal grid length.
        /// </summary>
        private double GetMinGridLength() {
            if (RefinementLevel == 0)
                return MaxGridLength;
            else
                return MaxGridLength / (RefinementLevel * 2);
        }

        /// <summary>
        /// Volume of the fluid domain
        /// </summary>
        [DataMember]
        private double DomainVolume => ((FSI_Control)Control).DomainVolume;

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
        /// MinimalDistanceForCollision
        /// </summary>
        [DataMember]
        private double MinimalDistanceForCollision => ((FSI_Control)Control).minDistanceThreshold;

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
        /// Fix the position of all particles in the control file.
        /// </summary>
        [DataMember]
        private bool FixPosition => ((FSI_Control)Control).fixPosition;

        /// <summary>
        /// The position of the (horizontal and vertical) boundaries.
        /// </summary>
        /// <remarks>
        /// First entry: vertical [0] or horizontal [1]
        /// Second entry: left/lower wall [0] or right/upper wall [1]
        /// </remarks>
        [DataMember]
        private double[][] BoundaryCoordinates => ((FSI_Control)Control).BoundaryPositionPerDimension;

        /// <summary>
        /// Saves the physical data of all particles
        /// </summary>
        private TextWriter logPhysicalDataParticles;

        /// <summary>
        /// Creates Navier-Stokes and continuity equation. 
        /// The nonlinear term is neglected, i.e. unsteady Stokes equation is solved, if <see cref="Control.PhysicalParameters.IncludeConvection"/> is false.
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

            IBM_Op = new XSpatialOperatorMk2(
                __DomainVar: DomNameSelected, __ParameterVar: Params, __CoDomainVar: CodNameSelected,
                QuadOrderFunc: (A, B, C) => HMForder,
                __Species: FluidSpecies.Select(id => LsTrk.GetSpeciesName(id)));
            IBM_Op.AgglomerationThreshold = this.Control.AgglomerationThreshold;

            IBM_Op.FreeMeanValue[VariableNames.Pressure] = !this.boundaryCondMap.DirichletPressureBoundary;

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
                        var convectionAtIB = new Solution.NSECommon.Operator.Convection.FSI_ConvectionAtIB(d, spatialDim, LsTrk, boundaryCondMap, ParticleList.ToArray(), UseMovingMesh, GetMinGridLength());
                        comps.Add(convectionAtIB);
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
                var viscousAtIB = new Solution.NSECommon.Operator.Viscosity.FSI_ViscosityAtIB(d, spatialDim, LsTrk, penalty, ComputePenaltyIB, FluidViscosity, ParticleList.ToArray(), GetMinGridLength());
                comps.Add(viscousAtIB);
            }

            // Continuum equation
            // =============================
            for (int d = 0; d < spatialDim; d++) {
                var src = new Divergence_DerivativeSource(d, spatialDim);
                var flx = new Divergence_DerivativeSource_Flux(d, boundaryCondMap);
                IBM_Op.EquationComponents["div"].Add(src);
                IBM_Op.EquationComponents["div"].Add(flx);
            }
            var divPen = new Solution.NSECommon.Operator.Continuity.FSI_DivergenceAtIB(spatialDim, LsTrk, ParticleList.ToArray(), GetMinGridLength());
            IBM_Op.EquationComponents["div"].Add(divPen);

            // temporal operator
            // =================
            var tempOp = new ConstantXTemporalOperator(IBM_Op, 0.0);
            foreach (var kv in this.MassScale) {
                tempOp.DiagonalScale[LsTrk.GetSpeciesName(kv.Key)].SetV(kv.Value.ToArray());
            }
            IBM_Op.TemporalOperator = tempOp;

            // Finalize
            // ========

            IBM_Op.Commit();

            CreateTimestepper();
        }

        /// <summary>
        /// Creates the BDF-Time-stepper
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
                case LevelSetHandling.StrangSplitting:
                case LevelSetHandling.FSILieSplittingFullyCoupled:
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

            BDFOrder = bdfOrder;
            m_BDF_Timestepper = new XdgBDFTimestepping(
                Fields: ArrayTools.Cat(Velocity, Pressure),
                __Parameters: IBM_Op.InvokeParameterFactory(ArrayTools.Cat(Velocity, Pressure)),
                IterationResiduals: ArrayTools.Cat(ResidualMomentum, ResidualContinuity),
                LsTrk: LsTrk,
                DelayInit: true,
                _ComputeOperatorMatrix: DelComputeOperatorMatrix,
                abstractOperator: IBM_Op,
                _UpdateLevelset: DelUpdateLevelset,
                BDForder: bdfOrder,
                _LevelSetHandling: ((FSI_Control)Control).Timestepper_LevelSetHandling,
                _MassMatrixShapeandDependence: MassMatrixShape,
                _SpatialOperatorType: SpatialOp,
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
                Timestepper_Init = TimeStepperInit.SingleInit
            };
        }

        /// <summary>
        /// Calls level set update depending on level set handling method.
        /// </summary>
        public override double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {
            double forces_PResidual = 0;
            switch (((FSI_Control)Control).Timestepper_LevelSetHandling) {
                case LevelSetHandling.None:
                ScalarFunction Posfunction = NonVectorizedScalarFunction.Vectorize(((FSI_Control)Control).MovementFunc, phystime);
                LevSet.ProjectField(Posfunction);
                LsTrk.UpdateTracker(phystime + dt);
                break;

                case LevelSetHandling.Coupled_Iterative: {
                    UpdateLevelSetParticles(phystime);
                    int iterationCounter = 0;
                    ParticleHydrodynamics AllParticleHydrodynamics = new ParticleHydrodynamics(LsTrk, spatialDim);
                    forces_PResidual = iterationCounter == 0 ? double.MaxValue : AllParticleHydrodynamics.CalculateParticleResidual(ref iterationCounter, ((FSI_Control)Control).fullyCoupledSplittingMaxIterations); ;
                    Console.WriteLine("Current forces_PResidual:   " + forces_PResidual);
                    throw new NotImplementedException("Moving interface solver will be implemented in the near future");
                }

                case LevelSetHandling.Coupled_Once:
                case LevelSetHandling.LieSplitting:
                case LevelSetHandling.FSILieSplittingFullyCoupled:
                case LevelSetHandling.StrangSplitting:
                UpdateLevelSetParticles(phystime);
                break;

                default:
                throw new ApplicationException("unknown 'LevelSetMovement': " + ((FSI_Control)Control).Timestepper_LevelSetHandling);
            }
            return forces_PResidual;
        }

        /// <summary>
        /// Level-set update.
        /// </summary>
        /// <remarks>
        /// Updates the level-set of each particle based on a point-wise max-function: levelSetFunction(X) = max(levelSetFunction(X), CurrentParticleLevelSet(X)).
        /// </remarks>
        /// <param name="phystime">
        /// The current time.
        /// </param>
        private void UpdateLevelSetParticles(double phystime) {
            if (!IsFullyCoupled || iterationCounter <= 1) {
                LevelSetUpdate = new FSILevelSetUpdate(GridData, GetMinGridLength());
                int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                CellMask allParticleMask = null;
                CellMask coloredCellMask = null;
                if (CellColor != null) {
                    DeleteParticlesOutsideOfDomain();
                    CreateDuplicateParticleAtPeriodicBoundary();
                    SwitchDuplicateAndMasterParticle();
                }

                CellColor = CellColor == null ? LevelSetUpdate.InitializeColoring(LsTrk, ParticleList.ToArray(), MaxGridLength) : LevelSetUpdate.UpdateColoring(LsTrk);
                for (int i = 0; i < noOfLocalCells; i++) {
                    CellID.SetMeanValue(i, i);
                    if (CellColor[i] != 0)
                        CellColor[i] = 1;
                }
                SetColorDGField(CellColor);

                DGLevSet.Current.Clear();
                GlobalParticleColor = LevelSetUpdate.DetermineGlobalParticleColor(GridData, CellColor, ParticleList);
                int[] globalParticleColor = GlobalParticleColor.CloneAs();

                for (int c = 0; c < globalParticleColor.Length; c++) {
                    int currentColor = globalParticleColor[c];
                    if (LevelSetUpdate.CurrentProcessContainsCurrentColor(CellColor, currentColor, noOfLocalCells)) {
                        coloredCellMask = new CellMask(GridData, LevelSetUpdate.CreateBitArrayFromColoredCells(CellColor, currentColor, noOfLocalCells));
                        allParticleMask = allParticleMask == null ? coloredCellMask : allParticleMask.Union(coloredCellMask);

                        int[] particlesOfCurrentColor = LevelSetUpdate.FindParticlesWithSameColor(globalParticleColor, currentColor);
                        double levelSetFunctionParticlesPerColor(double[] X, double t) {
                            double levelSetFunction = int.MinValue;
                            for (int p = 0; p < particlesOfCurrentColor.Length; p++) {
                                Particle currentParticle = ParticleList[particlesOfCurrentColor[p]];
                                if (IsPeriodic[0] || IsPeriodic[1]) {
                                    if (levelSetFunction < currentParticle.LevelSetFunction(X, BoundaryCoordinates))
                                        levelSetFunction = currentParticle.LevelSetFunction(X, BoundaryCoordinates);
                                } else {
                                    if (levelSetFunction < currentParticle.LevelSetFunction(X))
                                        levelSetFunction = currentParticle.LevelSetFunction(X);
                                }
                                globalParticleColor[particlesOfCurrentColor[p]] = 0;
                            }
                            return levelSetFunction;
                        }
                        SetLevelSet(levelSetFunctionParticlesPerColor, CellMask.GetFullMask(GridData), phystime);
                    }
                }

                double levelSetFunctionFluid(double[] X, double t) { return -1; }
                CellMask fluidCells = allParticleMask != null ? allParticleMask.Complement() : CellMask.GetFullMask(GridData);
                SetLevelSet(levelSetFunctionFluid, fluidCells, phystime);
                if (allParticleMask.IsNullOrEmpty())// in case of restart
                    allParticleMask = CellMask.GetFullMask(GridData);
                PerformLevelSetSmoothing(allParticleMask);
            }

            try {
                LsTrk.UpdateTracker(phystime, __NearRegionWith: 2);
            } catch (LevelSetCFLException e) {// when ghost particles are added at the opposing side of the domain, the CFL exception is thrown. However, due to periodicity this is OK. 
                if (AddedNewDuplicate)
                    Console.WriteLine("Ghost particle added! I catched the exception (which is completely OK) " + e);
                //else
                    //throw e;
                AddedNewDuplicate = false;
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
        /// Set DG field for the particle color.
        /// </summary>
        /// <param name="coloredCells"></param>
        private void SetColorDGField(int[] coloredCells) {
            ushort[] regionsCode = LsTrk.Regions.RegionsCode;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < noOfLocalCells; j++) {
                ParticleColor.SetMeanValue(j, coloredCells[j]);
                LevelSetDistance.SetMeanValue(j, LevelSetTracker.DecodeLevelSetDist(regionsCode[j], 0));
            }
        }

        /// <summary>
        /// This method creates one to three duplicate particles at the other side of a periodic boundary,
        /// which are an exact copy of the master particle. As long as the master particle exists,
        /// i.e. is inside of the domain, the duplicates copy the behavior of the master. 
        /// When the master has left the domain, one of the duplicates becomes the new master, <see cref="SwitchDuplicateAndMasterParticle()"/>.
        /// </summary>
        /// <remarks>
        /// The particles do not "know" anything about the periodicity of the domain. 
        /// Hence, a single particle would just leave the domain at one side of boundary and not appear again at the other side.
        /// As simple and straightforward solution we copy each particle, which reaches the periodic boundary. 
        /// 1 particle is created if the master is close to a single periodic boundary.
        /// 2 particles are created if the master is close to two different periodic boundaries, i.e. at an edge of the domain.
        /// 3 particles might be created in some circumstances if the particle has a certain position at the edge of the domain and appears in all four edges.
        /// </remarks>
        private void CreateDuplicateParticleAtPeriodicBoundary() {
            if (spatialDim != 2)
                throw new NotImplementedException("Periodic boundaries not implemented for more than 2D!");

            for (int p = 0; p < ParticleList.Count(); p++) {
                Particle currentParticle = ParticleList[p];
                if (!currentParticle.IsMaster)// Duplicates cannot own other duplicates.
                    continue;
                List<Particle> duplicateParticles = new List<Particle>();
                int[] duplicateHierachy = currentParticle.MasterDuplicateIDs.CloneAs();
                Vector particlePosition = currentParticle.Motion.GetPosition();
                int idOffset = 0;

                for (int d1 = 0; d1 < spatialDim; d1++) { // which direction: x or y?
                    if (!IsPeriodic[d1])
                        continue;
                    for (int wallID = 0; wallID < spatialDim; wallID++) { // which wall left or right (x) and upper or lower (y)?
                        if (ParticleHasReachedPeriodicBoundary(currentParticle, d1, wallID)) {
                            duplicateHierachy[0] = p + 1;
                            Vector originInVirtualNeighbouringDomain;
                            if (d1 == 0)
                                originInVirtualNeighbouringDomain = new Vector(2 * BoundaryCoordinates[0][1 - wallID], 0);
                            else
                                originInVirtualNeighbouringDomain = new Vector(0, 2 * BoundaryCoordinates[1][1 - wallID]);

                            Particle duplicateParticle;
                            if (DuplicateExists(duplicateHierachy, d1 + 1)) {
                                duplicateParticle = ParticleList[duplicateHierachy[d1 + 1] - 1];
                            } else {
                                duplicateHierachy[d1 + 1] = ParticleList.Count() + idOffset + 1;
                                idOffset += 1;
                                duplicateParticle = currentParticle.CloneAs();
                                duplicateParticle.SetDuplicate();
                                duplicateParticle.Motion.SetDuplicatePosition(originInVirtualNeighbouringDomain + particlePosition);
                                duplicateParticles.Add(duplicateParticle.CloneAs());
                            }

                            if (d1 == 0) {
                                // For a duplicate in x-direction, test for periodic boundaries in y-direction for the newly created duplicate. 
                                // Only necessary to check one time, hence, we use the x-direction.
                                if (!DuplicateExists(duplicateHierachy, 3)) {
                                    for (int wallID2 = 0; wallID2 < spatialDim; wallID2++) {
                                        if (ParticleHasReachedPeriodicBoundary(duplicateParticle, 1, wallID2)) {
                                            originInVirtualNeighbouringDomain = new Vector(0, 2 * BoundaryCoordinates[1][1 - wallID2]);
                                            duplicateHierachy[3] = ParticleList.Count() + idOffset + 1;
                                            idOffset += 1;
                                            Particle duplicateParticleOfDuplicateParticle = currentParticle.CloneAs();
                                            duplicateParticleOfDuplicateParticle.SetDuplicate();
                                            duplicateParticleOfDuplicateParticle.Motion.SetDuplicatePosition(originInVirtualNeighbouringDomain + duplicateParticle.Motion.GetPosition());
                                            duplicateParticles.Add(duplicateParticleOfDuplicateParticle.CloneAs());
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (duplicateParticles.Count >= 1) {
                    AddedNewDuplicate = true;
                    currentParticle.SetDuplicateHierachy(duplicateHierachy);
                    for (int p2 = 0; p2 < duplicateParticles.Count(); p2++) {
                        duplicateParticles[p2].SetDuplicateHierachy(duplicateHierachy);
                    }
                    ParticleList.AddRange(duplicateParticles);
                }
            }
        }

        /// <summary>
        /// Test whether is duplicate particle exists or not.
        /// </summary>
        /// <param name="duplicateHierachy"></param>
        /// <param name="hierachyPosition"></param>
        /// <returns></returns>
        private static bool DuplicateExists(int[] duplicateHierachy, int hierachyPosition) {
            return duplicateHierachy[hierachyPosition] > 0;
        }

        /// <summary>
        /// Needed to catch the CFL-exception in case of a new duplicate particle.
        /// </summary>
        bool AddedNewDuplicate = false;

        /// <summary>
        /// If the center of mass of a master particle has left the domain, one of its duplicates is chosen to be the new master.
        /// The center of mass of the new master has to be inside of the domain.
        /// </summary>
        private void SwitchDuplicateAndMasterParticle() {
            for (int p = 0; p < ParticleList.Count(); p++) {
                Particle currentParticle = ParticleList[p];
                if (!currentParticle.IsMaster)
                    continue;

                if (!CenterOfMassIsInsideOfDomain(currentParticle)) {
                    int oldMasterID = p + 1;
                    int[] duplicateHierachy = currentParticle.MasterDuplicateIDs;
                    for (int j = 1; j < duplicateHierachy.Length; j++) {
                        if (DuplicateExists(duplicateHierachy, j)) {
                            Particle currentDuplicate = ParticleList[duplicateHierachy[j] - 1];
                            if (CenterOfMassIsInsideOfDomain(currentDuplicate)) {
                                int newMasterID = duplicateHierachy[j];
                                int[] newDuplicateHierachy = duplicateHierachy.CloneAs();
                                newDuplicateHierachy[0] = newMasterID;
                                newDuplicateHierachy[j] = oldMasterID;
                                currentDuplicate.SetMaster(currentParticle.Motion.CloneAs());
                                currentParticle.SetDuplicate();
                                for (int i = 0; i < duplicateHierachy.Length; i++) {
                                    if (DuplicateExists(duplicateHierachy, i))
                                        ParticleList[duplicateHierachy[i] - 1].MasterDuplicateIDs = newDuplicateHierachy.CloneAs();
                                }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// If a particle has left the domain entirely it is removed from the <see cref="ParticleList"/>.
        /// </summary>
        private void DeleteParticlesOutsideOfDomain() {
            for (int p = 0; p < ParticleList.Count(); p++) {
                for (int d = 0; d < spatialDim; d++) {
                    for (int wallID = 0; wallID < spatialDim; wallID++) {
                        Particle currentParticle = ParticleList[p];
                        Vector particlePosition = currentParticle.Motion.GetPosition();
                        double distance = particlePosition[d] - BoundaryCoordinates[d][wallID];
                        double particleMaxLength = currentParticle.GetLengthScales().Min();
                        if (Math.Abs(distance) > particleMaxLength && !CenterOfMassIsInsideOfDomain(currentParticle)) {
                            if (ParticleHasLeftTheDomain(currentParticle)) {
                                int[] duplicateHierachy = currentParticle.MasterDuplicateIDs.CloneAs();
                                for (int g = 0; g < duplicateHierachy.Length; g++) {
                                    if (duplicateHierachy[g] == p + 1) {
                                        duplicateHierachy[g] = 0;
                                    }
                                }
                                for (int g = 0; g < duplicateHierachy.Length; g++) {
                                    if (duplicateHierachy[g] > 0)
                                        ParticleList[duplicateHierachy[g] - 1].MasterDuplicateIDs = duplicateHierachy.CloneAs();
                                }
                                ParticleList.RemoveAt(p);
                                if (p >= ParticleList.Count())// already the last particle, no further action needed!
                                    return;
                                // In order to keep the references in the duplicateHierachy correct, we replace the removed particle with the last particle.
                                Particle lastParticle = ParticleList.Last();
                                duplicateHierachy = lastParticle.MasterDuplicateIDs.CloneAs();
                                for (int g = 0; g < duplicateHierachy.Length; g++) {
                                    if (duplicateHierachy[g] == ParticleList.Count() + 1) {
                                        duplicateHierachy[g] = p + 1;
                                    }
                                }
                                ParticleList.Insert(p, lastParticle);
                                ParticleList.RemoveAt(ParticleList.Count() - 1);
                                for (int g = 0; g < duplicateHierachy.Length; g++) {
                                    if (duplicateHierachy[g] > 0)
                                        ParticleList[duplicateHierachy[g] - 1].MasterDuplicateIDs = duplicateHierachy.CloneAs();
                                }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Check whether the particle <paramref name="currentParticle"/> has reached a periodic boundary.
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <param name="dimension"></param>
        /// <param name="wallID"></param>
        /// <returns></returns>
        private bool ParticleHasReachedPeriodicBoundary(Particle currentParticle, int dimension, int wallID) {
            Vector particlePosition = currentParticle.Motion.GetPosition();
            double distance = particlePosition[dimension] - BoundaryCoordinates[dimension][wallID];
            double particleMaxLength = currentParticle.GetLengthScales().Max();
            if (Math.Abs(distance) < particleMaxLength) {
                if (dimension == 0)
                    currentParticle.ClosestPointOnOtherObjectToThis = new Vector(BoundaryCoordinates[dimension][wallID], particlePosition[1]);
                else
                    currentParticle.ClosestPointOnOtherObjectToThis = new Vector(particlePosition[0], BoundaryCoordinates[dimension][wallID]);
                ParticleCollision periodicCollision = new ParticleCollision(GetMinGridLength());
                periodicCollision.CalculateMinimumDistance(currentParticle, out Vector _, out Vector _, out bool Overlapping);
                return Overlapping;
            }
            return false;
        }

        /// <summary>
        /// Check whether the particle <paramref name="currentParticle"/> has left the domain.
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private bool ParticleHasLeftTheDomain(Particle currentParticle) {
            for (int d = 0; d < spatialDim; d++) {
                for (int wallID = 0; wallID < spatialDim; wallID++) {
                    Vector particlePosition = currentParticle.Motion.GetPosition();
                    if (d == 0)
                        currentParticle.ClosestPointOnOtherObjectToThis = new Vector(BoundaryCoordinates[d][wallID], particlePosition[1]);
                    else
                        currentParticle.ClosestPointOnOtherObjectToThis = new Vector(particlePosition[0], BoundaryCoordinates[d][wallID]);
                    ParticleCollision periodicCollision = new ParticleCollision(GetMinGridLength());
                    periodicCollision.CalculateMinimumDistance(currentParticle, out Vector _, out Vector _, out bool Overlapping);
                    if (Overlapping)
                        return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Check whether the center of mass of the particle <paramref name="currentParticle"/> is inside of the domain.
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private bool CenterOfMassIsInsideOfDomain(Particle currentParticle) {
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

        /// <summary>
        /// Calculates the new time-step if dynamic time-stepping is used. Depends on the movement of the particles. 
        /// </summary>
        /// <param name="phystime"></param>
        /// <param name="TimestepInt"></param>
        /// <returns></returns>
        private double CalculateTimeStepSize(double phystime, int TimestepInt) {
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
                    foreach (Particle p in ParticleList) {
                        double expectedVelocity = p.Motion.GetTranslationalVelocity(0).L2Norm() + (p.Motion.GetHydrodynamicForces(0) * DtMax / p.Motion.ParticleMass).L2Norm();
                        double expectedRotation = (p.Motion.GetRotationalVelocity() + p.Motion.GetHydrodynamicTorque() * DtMax / p.MomentOfInertia) * p.GetLengthScales().Max();
                        maxVelocityL2Norm = Math.Max(maxVelocityL2Norm, expectedVelocity + Math.Abs(expectedRotation));
                    }
                    dt = Math.Min(DtMax, GetMinGridLength() / (2 * maxVelocityL2Norm));
                    if (dt / oldTimestep > 1.1)
                        dt = oldTimestep * 1.1;
                    if (dt / oldTimestep < 0.9)
                        dt = oldTimestep * 0.9;
                    dt = dt.MPIMin();
                    if (dt < 1e-6)
                        dt = 1e-6;
                    foreach (Particle p in ParticleList) {
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

        /// <summary>
        /// runs solver one step?!
        /// </summary>
        /// <param name="TimestepInt">
        /// Time-step number
        /// </param>
        /// <param name="phystime">
        /// Physical time
        /// </param>
        /// <param name="dt">
        /// Time-step size
        /// </param>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            // initialize
            ResLogger.TimeStep = TimestepInt;
            BDFSchemeCoeffs[] m_TSCchain = BDFCommon.GetChain(BDFOrder);
            int S = m_TSCchain[0].S;
            dt = CalculateTimeStepSize(phystime, TimestepInt);
            if (initAddedDamping) {
                foreach (Particle p in ParticleList) {
                    if (p.Motion.UseAddedDamping) {
                        p.Motion.CalculateDampingTensor(p, LsTrk, FluidViscosity, FluidDensity, DtMax);
                        p.Motion.ExchangeAddedDampingTensors();
                    }
                }
                initAddedDamping = false;
            }
            // used later to check if there is exactly one push per time-step
            int OldPushCount = LsTrk.PushCount;

            // only particle motion & collisions, no flow solver
            // =================================================
            if (((FSI_Control)Control).pureDryCollisions) {
                LsTrk.PushStacks(); // in other branches, called by the BDF time-stepper
                DGLevSet.Push();
                Auxillary.ParticleState_MPICheck(ParticleList, GridData, MPISize);

                // physics
                // -------------------------------------------------
                ParticleHydrodynamics AllParticleHydrodynamics = new ParticleHydrodynamics(LsTrk, spatialDim);
                CalculateParticleForcesAndTorque(AllParticleHydrodynamics);
                CalculateParticleVelocity(ParticleList, dt, 0);
                CalculateCollision(ParticleList, dt);
                CalculateParticlePosition(dt);
                UpdateLevelSetParticles(phystime);

                // print
                // -------------------------------------------------
                Auxillary.PrintResultToConsole(ParticleList, LsTrk, 0, 0, phystime, TimestepInt, DomainVolume, ((FSI_Control)Control).FullOutputToConsole);
                LogPhysicalData(phystime, TimestepInt);

            }
            // particle motion & collisions plus flow solver
            // =================================================
            else {
                iterationCounter = 0;
                double hydroDynForceTorqueResidual = double.MaxValue;
                int minimumNumberOfIterations = 4;
                ParticleHydrodynamics AllParticleHydrodynamics = new ParticleHydrodynamics(LsTrk, spatialDim);

                Auxillary.ParticleState_MPICheck(ParticleList, GridData, MPISize);

                AllParticleHydrodynamics.SaveHydrodynamicOfPreviousIteration(ParticleList);
                if (IsFullyCoupled) 
                    InitializeParticlePerIteration(ParticleList, TimestepInt);

                while (hydroDynForceTorqueResidual > HydrodynConvergenceCriterion || iterationCounter < minimumNumberOfIterations) {
                    Stopwatch stopWatch = new Stopwatch();
                    if (((FSI_Control)Control).FullOutputToConsole) 
                        stopWatch.Start();
                    
                    AllParticleHydrodynamics.SaveHydrodynamicOfPreviousIteration(ParticleList);

                    m_BDF_Timestepper.Solve(phystime, dt, false);

                    CalculateParticleForcesAndTorque(AllParticleHydrodynamics);
                    CalculateParticleVelocity(ParticleList, dt, iterationCounter);

                    if (!IsFullyCoupled)
                        break;

                    hydroDynForceTorqueResidual = AllParticleHydrodynamics.CalculateParticleResidual(ref iterationCounter, ((FSI_Control)Control).fullyCoupledSplittingMaxIterations);

                    if (((FSI_Control)Control).FullOutputToConsole) {
                        Console.WriteLine("Iteration: {1}, Residual:  {0}", hydroDynForceTorqueResidual, iterationCounter - 1);
                        TimeSpan ts = stopWatch.Elapsed;
                        string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
                        Console.WriteLine("RunTime per Iteration" + elapsedTime);
                    }
                }

                CalculateCollision(ParticleList, dt, FixPosition);
                CheckDuplicateParticles();
                if (!FixPosition) {
                    CalculateParticlePosition(dt);
                }

                Auxillary.PrintResultToConsole(ParticleList, LsTrk, FluidViscosity, FluidDensity, phystime, TimestepInt, DomainVolume, ((FSI_Control)Control).FullOutputToConsole);
                LogPhysicalData(phystime, TimestepInt);

                if (IsFullyCoupled) {// in other branches, called by the BDF time-stepper
                    LsTrk.IncreaseHistoryLength(S + 1);
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

        private void CheckDuplicateParticles() {
            for(int p = 0; p < ParticleList.Count(); p++) {
                int[] duplicateHierachy = ParticleList[p].MasterDuplicateIDs;
                for(int p1 = 0; p1 < duplicateHierachy.Length; p1++) {
                    if(duplicateHierachy[p1] > 0) {
                        if (ParticleList[p].Motion.GetTranslationalVelocity(0).Abs() != ParticleList[duplicateHierachy[p1] - 1].Motion.GetTranslationalVelocity(0).Abs())
                            throw new Exception("Duplicate particles with unequal velocity, that cant be! Particle " + p + " and " + (duplicateHierachy[p1] - 1));
                    }
                }
            }
        }

        private void CalculateParticleForcesAndTorque(ParticleHydrodynamics AllParticleHydrodynamics) {
            ParticleHydrodynamicsIntegration hydrodynamicsIntegration = new ParticleHydrodynamicsIntegration(2, Velocity, Pressure, LsTrk, FluidViscosity);
            AllParticleHydrodynamics.CalculateHydrodynamics(ParticleList, hydrodynamicsIntegration, FluidDensity, IsFullyCoupled);
        }

        /// <summary>
        /// Update of added damping tensors and prediction of hydrodynamics.
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="TimestepInt">
        /// #Time-step
        /// </param>
        private void InitializeParticlePerIteration(List<Particle> Particles, int TimestepInt) {
            for (int p = 0; p < Particles.Count; p++) {
                Particle currentParticle = Particles[p];
                if (currentParticle.Motion.UseAddedDamping) {
                    currentParticle.Motion.UpdateDampingTensors();
                }
                currentParticle.Motion.SaveHydrodynamicsOfPreviousTimestep();
                currentParticle.Motion.PredictForceAndTorque(currentParticle.ActiveStress, currentParticle.Circumference, TimestepInt, FluidViscosity, FluidDensity, DtMax);
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
        private void CalculateParticleVelocity(List<Particle> Particles, double dt, int IterationCounter) {
            foreach (Particle p in Particles) {
                if (IterationCounter == 0) {
                    p.Motion.SaveVelocityOfPreviousTimestep();
                }

                if (p.IsMaster) {
                    p.Motion.UpdateParticleVelocity(dt);

                    for (int g = 1; g < p.MasterDuplicateIDs.Length; g++) {
                        if (p.MasterDuplicateIDs[g] >= 1) {
                            Particle ghost = Particles[p.MasterDuplicateIDs[g] - 1];
                            if (ghost.IsMaster)
                                throw new Exception("A ghost particle is considered to be a master, that can't be!");
                            ghost.Motion.CopyNewVelocity(p.Motion.GetTranslationalVelocity(), p.Motion.GetRotationalVelocity());
                        }
                    }
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
            for (int p = 0; p < ParticleList.Count; p++) {
                Particle particle = ParticleList[p];
                particle.Motion.UpdateParticlePositionAndAngle(dt);
            }
        }
                
        /// <summary>
        /// Creates a log file for the physical data of the particles. Only active if a database is specified.
        /// </summary>
        private void CreatePhysicalDataLogger() {
            if ((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                logPhysicalDataParticles = DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", "time-step", "particle", "time", "posX", "posY", "angle", "velX", "velY", "rot", "fX", "fY", "T"));
            }
        }

        /// <summary>
        /// Writes the physical data of the particles to a log file.
        /// </summary>
        /// <param name = phystime>
        /// </param>
        private void LogPhysicalData(double phystime, int timestepNo) {
            if ((MPIRank == 0) && (logPhysicalDataParticles != null)) {
                for (int p = 0; p < ParticleList.Count(); p++) {
                    logPhysicalDataParticles.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", timestepNo, p, phystime, ParticleList[p].Motion.GetPosition(0)[0], ParticleList[p].Motion.GetPosition(0)[1], ParticleList[p].Motion.GetAngle(0), ParticleList[p].Motion.GetTranslationalVelocity(0)[0], ParticleList[p].Motion.GetTranslationalVelocity(0)[1], ParticleList[p].Motion.GetRotationalVelocity(0), ParticleList[p].Motion.GetHydrodynamicForces(0)[0], ParticleList[p].Motion.GetHydrodynamicForces(0)[1], ParticleList[p].Motion.GetHydrodynamicTorque(0)));
                    logPhysicalDataParticles.Flush();
                }
            }
        }

        /// <summary>
        /// Triggers the collision detection, which triggers the calculation of the collisions. Note on parallelization: All particle operations are carried out on all processes, hence no communication is necessary.
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="dt">
        /// Time-step
        /// </param>
        /// <param name="DetermineOnlyOverlap">
        /// Set true if you are only interested in overlapping particles and not the actual distance between different particles, e.g. as check for the initialization of static particles. 
        /// </param>
        private void CalculateCollision(List<Particle> Particles, double dt, bool DetermineOnlyOverlap = false) {
            foreach (Particle p in Particles) {
                p.IsCollided = false;
            }
            ParticleCollision Collision = new ParticleCollision(GetMinGridLength(), ((FSI_Control)Control).CoefficientOfRestitution, dt, ((FSI_Control)Control).WallPositionPerDimension, ((FSI_Control)Control).BoundaryIsPeriodic, MinimalDistanceForCollision, DetermineOnlyOverlap);
            Collision.Calculate(ParticleList.ToArray());

            // The following (collision detection based on cell color might be  more efficient, however, it doesn't work with periodic boundaries.
            // Only particles with the same color a close to each other, thus, we only test for collisions within those particles.
            // Determine color.
            // =================================================
            //int[] _GlobalParticleColor = GlobalParticleColor.CloneAs();
            //for (int i = 0; i < _GlobalParticleColor.Length; i++) {
            //    int CurrentColor = _GlobalParticleColor[i];
            //    if (CurrentColor == 0)
            //        continue;
            //    int[] ParticlesOfCurrentColor = levelSetUpdate.FindParticlesWithSameColor(_GlobalParticleColor, CurrentColor);
            //    // Multiple particles with the same color, trigger collision detection
            //    // =================================================
            //    if (ParticlesOfCurrentColor.Length >= 1 && CurrentColor != 0) {
            //        Particle[] currentParticles = new Particle[ParticlesOfCurrentColor.Length];
            //        for (int j = 0; j < ParticlesOfCurrentColor.Length; j++) {
            //            currentParticles[j] = ParticleList[ParticlesOfCurrentColor[j]];
            //        }
            //        ParticleCollision Collision = new ParticleCollision(GetMinGridLength(), ((FSI_Control)Control).CoefficientOfRestitution, dt, ((FSI_Control)Control).WallPositionPerDimension, ((FSI_Control)Control).BoundaryIsPeriodic, MinimalDistanceForCollision, DetermineOnlyOverlap);
            //        Collision.Calculate(ParticleList.ToArray());
            //        //Collision.Calculate(currentParticles);
            //    }
            //    // Remove already examined particles/colors from array
            //    // =================================================
            //    for (int j = 0; j < _GlobalParticleColor.Length; j++) {
            //        if (_GlobalParticleColor[j] == CurrentColor)
            //            _GlobalParticleColor[j] = 0;
            //    }
            //}
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
                CellMask cutCells = LsTrk.Regions.GetCutCellMask();
                GridRefinementController gridRefinementController = new GridRefinementController(gridData, cutCells, null);
                AnyChangeInGrid = gridRefinementController.ComputeGridChange(GetCellsToRefine(), out CellsToRefineList, out Coarsening);
            }
            if (AnyChangeInGrid) {
                int[] consoleRefineCoarse = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                long oldJ = GridData.CellPartitioning.TotalLength;
                Console.WriteLine("Refining " + consoleRefineCoarse[0] + " of " + oldJ + " cells");
                Console.WriteLine("Coarsening " + consoleRefineCoarse[1] + " of " + oldJ + " cells");
                Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Count().MPISum());
                newGrid = gridData.Adapt(CellsToRefineList, Coarsening, out old2NewGrid);
            } else {
                newGrid = null;
                old2NewGrid = null;
            }
        }

        private int[] GetCellsToRefine() {
            int refinementLevel = ((FSI_Control)Control).RefinementLevel;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellsToRefine = new int[noOfLocalCells];
            for (int p = 0; p < ParticleList.Count(); p++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    if (cellsToRefine[j] < refinementLevel) {
                        Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                        cellsToRefine[j] = (ParticleList[p].Contains(cellCenter) && !ParticleList[p].Contains(cellCenter, -MaxGridLength)) ? refinementLevel : 0;
                    }
                }
            }
            return cellsToRefine;
        }



        /// <summary>
        /// For restarting calculations, its important to reload old solutions if one uses a higher order method in time
        /// </summary>
        /// <param name="time"></param>
        /// <param name="timestep"></param>
        public override void PostRestart(double time, TimestepNumber timestep) {
            if (!((FSI_Control)Control).IsRestart)
                return;

            IFileSystemDriver fsDriver = this.DatabaseDriver.FsDriver;
            string pathToOldSessionDir = Path.Combine(fsDriver.BasePath, "sessions", this.CurrentSessionInfo.RestartedFrom.ToString());
            string pathToPhysicalData = MPIRank == 0 ? Path.Combine(pathToOldSessionDir, "PhysicalData.txt") : "";
            pathToPhysicalData = pathToPhysicalData.MPIBroadcast(0);

            string[] records = File.ReadAllLines(pathToPhysicalData);
            int timestepIndexOffset = 0;
            for (int r = 1; r < records.Length; r++) {// 0th line does not contain data
                string currentLine = records[r];
                string[] currentLineFields = currentLine.Split(',');
                if (timestep.MajorNumber == Convert.ToInt32(currentLineFields[0])) {
                    timestepIndexOffset = r;
                    break;
                }
            }

            int historyLength = 3;
            for (int t = historyLength; t > 0; t--) {
                for (int p = 0; p < ParticleList.Count(); p++) {
                    Particle currentParticle = ParticleList[p];
                    int index = timestepIndexOffset - ParticleList.Count() * (t - 1) + p;
                    string currentLine = records[index];
                    string[] currentLineFields = currentLine.Split(',');
                    double[] position = new double[2];
                    double[] translationalVelocity = new double[2];
                    double[] force = new double[2];
                    double[] physicalData = currentLineFields.Select(eachElement => Convert.ToDouble(eachElement)).ToArray();
                    position[0] = Convert.ToDouble(currentLineFields[3]);
                    position[1] = Convert.ToDouble(currentLineFields[4]);
                    force[0] = Convert.ToDouble(currentLineFields[9]);
                    force[1] = Convert.ToDouble(currentLineFields[10]);
                    double angle = Convert.ToDouble(currentLineFields[5]) * 360 / (2 * Math.PI);
                    translationalVelocity[0] = Convert.ToDouble(currentLineFields[6]);
                    translationalVelocity[1] = Convert.ToDouble(currentLineFields[7]);
                    double angularVelocity = Convert.ToDouble(currentLineFields[8]);
                    double torque = Convert.ToDouble(currentLineFields[11]);
                    currentParticle.Motion.InitializeParticlePositionAndAngle(new double[] { physicalData[3], physicalData[4] }, physicalData[5] * 360 / (2 * Math.PI), t);
                    currentParticle.Motion.InitializeParticleVelocity(new double[] { physicalData[6], physicalData[7] }, physicalData[8], t);
                    currentParticle.Motion.InitializeParticleForceAndTorque(new double[] { physicalData[9], physicalData[10] }, physicalData[11], t, StartDt);
                }
            }
            CellColor = null;
            UpdateLevelSetParticles(time);
            CreatePhysicalDataLogger();
            ((FSI_Control)Control).IsRestart = false;
        }



        /// <summary>
        /// over-ridden in oder to save the particles (<see cref="ParticleList"/>) to the database
        /// </summary>
        protected override TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {
            FSI_TimestepInfo tsi = new FSI_TimestepInfo(t, CurrentSessionInfo, timestepno, IOFields, ParticleList);
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
                } else {
                    return new BsonReader(s);
                }
            }
            JsonWriter GetJsonWriter(Stream s) {
                if (DebugSerialization) {
                    return new JsonTextWriter(new StreamWriter(s));
                } else {
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
        /// over-ridden in oder to save the particles (<see cref="ParticleList"/>) to the database
        /// </summary>
        protected override void OnRestartTimestepInfo(TimestepInfo tsi) {
            FSI_TimestepInfo fTsi = (FSI_TimestepInfo)tsi;
            // initialize particles
            ParticleList = fTsi.Particles.ToList();
            UpdateLevelSetParticles(fTsi.PhysicalTime);
        }
    }
}