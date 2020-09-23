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

using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.IO;
using System.Diagnostics;
using System.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.IBM_Solver {

    /// <summary>
    /// A incompressible Navier-Stokes solver with the possibility of using non moving immersed boundaries.
    /// </summary>
    public class IBM_SolverMain : Application<IBM_Control> {


        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {

            _Main(args, false, delegate () {
                var p = new IBM_SolverMain();
                return p;
            });
            //Console.ReadKey();
        }

        #region instantiation
#pragma warning disable 649
        /// <summary>
        /// Pressure
        /// </summary>
        [InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField Pressure;

        /// <summary>
        /// velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            null,
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Velocity;

        /// <summary>
        /// Level-Set tracker
        /// </summary>
        [LevelSetTracker("-:A +:B", 2)]
        public LevelSetTracker LevsetTracker;

        /// <summary>
        /// the DG representation of the level set.
        /// This one is used for level-set evolution in time; it is in general discontinuous.
        /// </summary>
        [InstantiateFromControlFile("PhiDG", "PhiDG", IOListOption.ControlFileDetermined)]
        public ScalarFieldHistory<SinglePhaseField> DGLevSet;

        /// <summary>
        /// The  continuous level set field which defines the XDG space; 
        /// it is obtained from the projection of the discontinuous <see cref="DGLevSet"/> onto the 
        /// continuous element space.
        /// </summary>
        [InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        public LevelSet LevSet;

        ///// <summary>
        ///// Curvature; DG-polynomial degree should be 2 times the polynomial degree of <see cref="LevSet"/>.
        ///// </summary>
        //[InstantiateFromControlFile("Curvature", "Curvature", IOListOption.ControlFileDetermined)]
        //public SinglePhaseField Curvature;

        /// <summary>
        /// Residual of the continuity equation
        /// </summary>
        [InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ResidualContinuity;


        /// <summary>
        /// Residual in the momentum equation.
        /// </summary>
        [InstantiateFromControlFile(new string[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> ResidualMomentum;



#pragma warning restore 649
        #endregion


        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        virtual protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double rho = this.Control.PhysicalParameters.rho_A;

                int D = this.GridData.SpatialDimension;

                double[] _rho = new double[D + 1];
                if (!this.Control.IsStationary)
                    _rho.SetAll(rho);
                //No MassMatrix for the pressure
                _rho[D] = 0;

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>> {
                    { this.LsTrk.GetSpeciesId("A"), _rho }
                };

                return R;
            }
        }

        ///// <summary>
        ///// Continuous high-order finite element representation of the level-set field.
        ///// </summary>
        //SpecFemField SmoothedLevelSet;

        /// <summary>
        /// Links edge tags (<see cref="IGeometricalEdgeData.EdgeTags"/>) and
        /// boundary conditions in the control object (<see cref="BoSSS.Solution.Control.AppControl.BoundaryValues"/>).
        /// </summary>
        protected IncompressibleBoundaryCondMap boundaryCondMap;

        CoordinateVector m_CurrentSolution;

        /// <summary>
        /// Current velocity and pressure.
        /// </summary>
        protected CoordinateVector CurrentSolution {
            get {
                if (m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ArrayTools.Cat(this.Velocity, this.Pressure));
                } else {
                    for (int d = 0; d < base.GridData.SpatialDimension; d++) {
                        Debug.Assert(object.ReferenceEquals(m_CurrentSolution.Mapping.Fields[d], this.Velocity[d]));
                    }
                }

                return m_CurrentSolution;
            }
        }


        CoordinateVector m_CurrentResidual;

        /// <summary>
        /// Residual vector (residual of momentum and continuity equation) for the current solution.
        /// </summary>
        CoordinateVector CurrentResidual {
            get {
                if (m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(ArrayTools.Cat(ResidualMomentum, ResidualContinuity));
                }
                return m_CurrentResidual;
            }
        }
       
        /// <summary>
        /// equal to <see cref="PhysicalParameters.IncludeConvection"/>
        /// </summary>
        protected bool U0MeanRequired {
            get {
                return Control.PhysicalParameters.IncludeConvection;
            }
        }

        /// <summary>
        /// the spatial operator of the incompressible Navier-Stokes equation
        /// </summary>
        protected XSpatialOperatorMk2 IBM_Op;

        /// <summary>
        /// the Jacobian of <see cref="IBM_Op"/>
        /// </summary>
        protected XSpatialOperatorMk2 IBM_Op_Jacobian;

        public string SessionPath;

        /// <summary>
        /// Integration degree of HMF used throughout the application: this should ensure that
        /// only one HMF rule is created.
        /// </summary>
        public int HMForder {
            get {
                int VelDeg = this.Velocity.Max(field => field.Basis.Degree);
                int Order = (VelDeg * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2));
                Order += 1; // safety factor
                return Order;
            }
        }

        /// <summary>
        /// Ye good old time-steppa
        /// </summary>
        protected XdgBDFTimestepping m_BDF_Timestepper;

       
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

            //// Write out Multigrid Levels
            //this.MGColoring = new SinglePhaseField[base.MultigridSequence.Length];
            //for (int iLevel = 0; iLevel < base.MultigridSequence.Length; iLevel++) {
            //    this.MGColoring[iLevel] = new SinglePhaseField(new Basis(this.GridData, 0), "MGColoring_level_f" + iLevel);
            //    base.MultigridSequence[iLevel].ColorDGField(this.MGColoring[iLevel]);
            //}
            //Tecplot.PlotFields(MGColoring, "MultigridLevels", 0, 0);

            // =================================
            // create operator
            // =================================
            if (IBM_Op == null) {

                string[] CodNameSelected = new string[0];
                string[] DomNameSelected = new string[0];

                int D = this.GridData.SpatialDimension;
                boundaryCondMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

                
                var IBM_Op_config = new NSEOperatorConfiguration {
                    convection = this.Control.PhysicalParameters.IncludeConvection,
                    continuity = true,
                    Viscous = true,
                    PressureGradient = true,
                    Transport = true,
                    CodBlocks = new bool[] { true, true },
                    DomBlocks = new bool[] { true, true },
                };

                // full operator:
                var CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D)).Cat("div");
                var Params = ArrayTools.Cat(VariableNames.Velocity0Vector(D),
                                            VariableNames.Velocity0MeanVector(D));
                var DomName = ArrayTools.Cat(VariableNames.VelocityVector(D), VariableNames.Pressure);

                // selected part:
                if (IBM_Op_config.CodBlocks[0])
                    CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(0, D));
                if (IBM_Op_config.CodBlocks[1])
                    CodNameSelected = ArrayTools.Cat(CodNameSelected, CodName.GetSubVector(D, 1));

                if (IBM_Op_config.DomBlocks[0])
                    DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(0, D));
                if (IBM_Op_config.DomBlocks[1])
                    DomNameSelected = ArrayTools.Cat(DomNameSelected, DomName.GetSubVector(D, 1));

                IBM_Op = new XSpatialOperatorMk2(DomNameSelected, Params, CodNameSelected,
                    (A, B, C) => this.HMForder, 
                    FluidSpecies.Select(sId => LsTrk.GetSpeciesName(sId)));

                IBM_Op.FreeMeanValue[VariableNames.Pressure] = !this.boundaryCondMap.DirichletPressureBoundary;


                // Momentum equation
                // =================
                AddBulkEquationComponentsToIBMOp(IBM_Op_config, CodName);
                AddInterfaceEquationComponentsToIBMOp(IBM_Op_config, CodName);

                // temporal operator
                // =================

                {
                    var tempOp = new ConstantXTemporalOperator(IBM_Op, 0.0);
                    foreach (var kv in this.MassScale) {
                        tempOp.DiagonalScale[LsTrk.GetSpeciesName(kv.Key)].SetV(kv.Value.ToArray());
                    }
                    IBM_Op.TemporalOperator = tempOp;

                }


                // Finalize
                // ========
                IBM_Op.Commit();


                //IBM_Op_Jacobian = IBM_Op.GetJacobiOperator();
            }

            // ==========================
            // create timestepper
            // ==========================
            var Unknowns = ArrayTools.Cat(this.Velocity, this.Pressure);
            var Residual = ArrayTools.Cat(this.ResidualMomentum, this.ResidualContinuity);
            if (L == null)
            {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                // Creation of time-integrator (initial, no balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                if (m_BDF_Timestepper == null)
                {
                    m_BDF_Timestepper = CreateTimeStepper(Unknowns, Residual);
                }
            }
            else
            {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // restore BDF time-stepper after grid redistribution (dynamic load balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                m_BDF_Timestepper.DataRestoreAfterBalancing(L, Unknowns, Residual, this.LsTrk, base.MultigridSequence);
            }
            Debug.Assert(m_BDF_Timestepper != null);
        }

        /// <summary>
        /// Ti
        /// </summary>
        /// <param name="Unknowns"></param>
        /// <param name="Residual"></param>
        protected virtual XdgBDFTimestepping CreateTimeStepper(IEnumerable<DGField> Unknowns, IEnumerable<DGField> Residual)
        {
            LevelSetHandling lsh = LevelSetHandling.None;
            SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;

            if (this.Control.PhysicalParameters.IncludeConvection)
            {
                SpatialOp = SpatialOperatorType.Nonlinear;
            }

            int bdfOrder;
            if (this.Control.Timestepper_Scheme == IBM_Control.TimesteppingScheme.CrankNicolson)
                bdfOrder = -1;
            else if (this.Control.Timestepper_Scheme == IBM_Control.TimesteppingScheme.ImplicitEuler)
                bdfOrder = 1;
            else if (this.Control.Timestepper_Scheme.ToString().StartsWith("BDF"))
                bdfOrder = Convert.ToInt32(this.Control.Timestepper_Scheme.ToString().Substring(3));
            else
                throw new NotImplementedException("todo");

            XdgBDFTimestepping m_BDF_Timestepper = new XdgBDFTimestepping(
                Unknowns,
                this.IBM_Op.InvokeParameterFactory(Unknowns),
                Residual,
                LsTrk,
                true,
                DelComputeOperatorMatrix,
                this.IBM_Op,
                DelUpdateLevelset,
                bdfOrder,
                lsh,
                MassMatrixShapeandDependence.IsTimeDependent,
                SpatialOp,
                this.MultigridOperatorConfig,
                base.MultigridSequence,
                this.FluidSpecies,
                this.HMForder,
                this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                false, this.Control.NonLinearSolver,
                this.Control.LinearSolver
                )
            {
                m_ResLogger = base.ResLogger,
                m_ResidualNames = ArrayTools.Cat(this.ResidualMomentum.Select(
                    f => f.Identification), this.ResidualContinuity.Identification),
                Timestepper_Init = Solution.Timestepping.TimeStepperInit.MultiInit
            };
            return m_BDF_Timestepper;
        }

        void AddBulkEquationComponentsToIBMOp(NSEOperatorConfiguration IBM_Op_config, string[] CodName)
        {
            int D = this.GridData.SpatialDimension;
            // convective part:
            if (IBM_Op_config.convection)
            {
                for (int d = 0; d < D; d++)
                {

                    var comps = IBM_Op.EquationComponents[CodName[d]];

                    var ConvBulk = new Solution.NSECommon.LinearizedConvection(D, boundaryCondMap, d);
                    //var ConvBulkUp = new UpwindConvection(D, boundaryCondMap, d, Control.PhysicalParameters.rho_A);
                    comps.Add(ConvBulk); // bulk component
                }
            }

            // pressure part:
            if (IBM_Op_config.PressureGradient)
            {
                for (int d = 0; d < D; d++)
                {
                    var comps = IBM_Op.EquationComponents[CodName[d]];
                    var pres = new PressureGradientLin_d(d, boundaryCondMap);
                    comps.Add(pres); // bulk component

                    // if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                    if (this.Control.FixedStreamwisePeriodicBC)
                    {
                        var presSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                        comps.Add(presSource);
                        // Jacobian operator: not required; reason: Source-Term with no dependence on domain variables
                    }
                }
            }

            // viscous part:
            if (IBM_Op_config.Viscous)
            {
                for (int d = 0; d < D; d++)
                {
                    var comps = IBM_Op.EquationComponents[CodName[d]];

                    double penalty_bulk = this.Control.AdvancedDiscretizationOptions.PenaltySafety;


                    //var Visc = new Solution.XNSECommon.Operator.Viscosity.ViscosityInBulk_GradUTerm(penalty, 1.0, BcMap, d, D, this.Control.PhysicalParameters.mu_A, 1, ViscosityImplementation.H);
                    var Visc = new swipViscosity_Term1(penalty_bulk, d, D, boundaryCondMap,
                        ViscosityOption.ConstantViscosity,
                        this.Control.PhysicalParameters.mu_A,// / this.Control.PhysicalParameters.rho_A
                        double.NaN, null);
                    comps.Add(Visc); // bulk component GradUTerm 
                }
            }

            // Continuum equation
            // ==================
            if (IBM_Op_config.continuity)
            {
                for (int d = 0; d < D; d++)
                {

                    var src = new Divergence_DerivativeSource(d, D);
                    var flx = new Divergence_DerivativeSource_Flux(d, boundaryCondMap);
                    IBM_Op.EquationComponents["div"].Add(src);
                    IBM_Op.EquationComponents["div"].Add(flx);


                    //var presStab = new PressureStabilization(1, this.GridData.Edges.h_max_Edge, 1 / this.Control.PhysicalParameters.mu_A);
                    //IBM_Op.EquationComponents["div"].Add(presStab);
                }


                //IBM_Op.EquationComponents["div"].Add(new PressureStabilization(1, 1.0 / this.Control.PhysicalParameters.mu_A));
            }
        }

        /// <summary>
        /// Setup for equations on interface
        /// </summary>
        /// <param name="IBM_Op_config">
        /// User settings
        /// </param>
        /// <param name="CodName">
        /// Domainnames of IBM solver
        /// </param>
        protected virtual void AddInterfaceEquationComponentsToIBMOp(NSEOperatorConfiguration IBM_Op_config, string[] CodName)
        {
            int D = this.GridData.SpatialDimension;
            for (int d = 0; d < D; d++){
                var comps = IBM_Op.EquationComponents[CodName[d]];

                if (IBM_Op_config.convection){
                    var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ConvectionAtIB(
                            d, D, LsTrk, this.Control.AdvancedDiscretizationOptions.LFFA, boundaryCondMap,
                            delegate (double[] X, double time) { return new double[] { 0.0, 0.0, 0.0, 0.0 }; }, this.Control.PhysicalParameters.rho_A, false);
                    //var ConvIB = new ConvectionAtIB(LsTrk, d, D, Control.PhysicalParameters.rho_A, false);

                    comps.Add(ConvIB); // immersed boundary component
                }

                if (IBM_Op_config.PressureGradient){
                    var presLs = new BoSSS.Solution.NSECommon.Operator.Pressure.PressureFormAtIB(d, D, LsTrk);
                    comps.Add(presLs); // immersed boundary component
                }

                if (IBM_Op_config.Viscous){
                    double _D = D;
                    double penalty_mul = this.Control.AdvancedDiscretizationOptions.PenaltySafety;
                    int degU = this.Velocity[0].Basis.Degree;
                    double _p = degU;
                    double penalty_base = (_p + 1) * (_p + _D) / D;
                    double penalty = penalty_base * penalty_mul;
                    var ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ViscosityAtIB(d, D, LsTrk,
                            penalty, this.ComputePenaltyIB,
                            this.Control.PhysicalParameters.mu_A,// / this.Control.PhysicalParameters.rho_A,
                            delegate (double[] X, double time) { return new double[] { 0.0, 0.0, 0.0, 0.0 }; });
                    comps.Add(ViscLs); // immersed boundary component
                }
            }

            if (IBM_Op_config.continuity){
                var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, LsTrk, 1,
                    delegate (double[] X, double time) { return new double[] { 0.0, 0.0, 0.0, 0.0 }; });
                IBM_Op.EquationComponents["div"].Add(divPen); // immersed boundary component 
            }
        }

        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
            m_CurrentResidual = null;
            m_CurrentSolution = null;
            IBM_Op = null;
        }

        void ParameterUpdate(IEnumerable<DGField> CurrentState, IEnumerable<DGField> ParameterVar) {
            int D = this.LsTrk.GridDat.SpatialDimension;

            if (Control.PhysicalParameters.IncludeConvection) {

                var U0_CurrentState = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());

                // this method assumes that the parameter velocity is reference-equal to the current state 
                for (int d = 0; d < D; d++)
                    if (!object.ReferenceEquals(CurrentState.ElementAt(d), ParameterVar.ElementAt(d)))
                        throw new ApplicationException("internal error");

                if (this.U0MeanRequired) {
                    VectorField<SinglePhaseField> U0mean = new VectorField<SinglePhaseField>(ParameterVar.Skip(D).Take(D).Select(f => (SinglePhaseField)f).ToArray());
                    foreach (var um in U0mean)
                        Debug.Assert(um.Basis.Degree == 0);

                    ComputeAverageU(U0_CurrentState, U0mean);
                } else {
                    //U0_U0mean = new SinglePhaseField[2 * D];
                }
            }
        }

        int DelComputeOperatorMatrix_CallCounter = 0;

        /// <summary>
        /// Used by <see cref="m_BDF_Timestepper"/> to compute operator matrices (linearizations) and/or to evaluate residuals of current solution.
        /// </summary>
        protected virtual void DelComputeOperatorMatrix(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {
            DelComputeOperatorMatrix_CallCounter++;
            int D = this.LsTrk.GridDat.SpatialDimension;

            // compute operator
            //Debug.Assert(OpMatrix.InfNorm() == 0.0);
            //Debug.Assert(OpAffine.L2Norm() == 0.0);
            // Create Parameters fields
            DGField[] Params;
            {
                var U0 = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());
                SinglePhaseField[] U0_U0mean;
                if (this.U0MeanRequired) {
                    Basis U0meanBasis = new Basis(GridData, 0);
                    VectorField<SinglePhaseField> U0mean = new VectorField<SinglePhaseField>(D, U0meanBasis, "U0mean_", SinglePhaseField.Factory);
                    U0_U0mean = ArrayTools.Cat<SinglePhaseField>(U0, U0mean);
                } else {
                    U0_U0mean = new SinglePhaseField[2 * D];
                }
                Params = ArrayTools.Cat<DGField>(U0_U0mean);
            }

            m_LenScales = AgglomeratedCellLengthScales[FluidSpecies[0]];

            // create matrix and affine vector:
            if (OpMatrix != null) {
                

                // using ad-hoc linearization:
                // - - - - - - - - - - - - - - 
                ParameterUpdate(CurrentState, Params);
                var mtxBuilder = IBM_Op.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping);
                mtxBuilder.time = phystime;
                mtxBuilder.CellLengthScales[FluidSpecies[0]] = AgglomeratedCellLengthScales[FluidSpecies[0]];
                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

                // using finite difference Jacobi:
                // - - - - - - - - - - - - - - - -
                //var mtxBuilder2 = IBM_Op.GetFDJacobianBuilder(LsTrk, CurrentState, Params, Mapping,
                //    ParameterUpdate,
                //    FluidSpecies);
                //mtxBuilder2.time = phystime;
                //mtxBuilder2.SpeciesOperatorCoefficients[FluidSpecies[0]].CellLengthScales = AgglomeratedCellLengthScales[FluidSpecies[0]];
                //mtxBuilder2.ComputeMatrix(OpMatrix, OpAffine);

                // using the other kind of Jacobi:
                // - - - - - - - - - - - - - - - -
                //var mtxBuilder3 = IBM_Op_Jacobian.GetMatrixBuilder(LsTrk, Mapping, CurrentState, Mapping, FluidSpecies);
                //mtxBuilder3.time = phystime;
                //mtxBuilder3.SpeciesOperatorCoefficients[FluidSpecies[0]].CellLengthScales = AgglomeratedCellLengthScales[FluidSpecies[0]];
                //mtxBuilder3.ComputeMatrix(OpMatrix, OpAffine);

#if DEBUG
                if (DelComputeOperatorMatrix_CallCounter == 1) {
                    int[] Uidx = SaddlePointProblemMapping.GetSubvectorIndices(true, D.ForLoop(i => i));
                    int[] Pidx = SaddlePointProblemMapping.GetSubvectorIndices(true, D);
                    CoordinateMapping Umap = this.Velocity.Mapping;
                    CoordinateMapping Pmap = this.Pressure.Mapping;
                    
                    var pGrad = new BlockMsrMatrix(Umap, Pmap);
                    var divVel = new BlockMsrMatrix(Pmap, Umap);
                    OpMatrix.AccSubMatrixTo(1.0, pGrad, Uidx, default(int[]), Pidx, default(int[]));
                    OpMatrix.AccSubMatrixTo(1.0, divVel, Pidx, default(int[]), Uidx, default(int[]));

                    var pGradT = pGrad.Transpose();
                    var Err = divVel.CloneAs();
                    Err.Acc(+1.0, pGradT);
                    double ErrInfAbs = Err.InfNorm();
                    double denom = Math.Max(pGradT.InfNorm(), Math.Max(pGrad.InfNorm(), divVel.InfNorm()));
                    double ErrInfRel = ErrInfAbs / denom;
                    if (ErrInfRel >= 1e-8)
                        throw new ArithmeticException("Stokes discretization error: | div + grad^t |oo is high; absolute: " + ErrInfAbs + ", relative: " + ErrInfRel + " (denominator: " + denom + ")");
                    //Console.WriteLine("Stokes discretization error: | div - grad ^ t |oo is high; absolute: " + ErrInfAbs + ", relative: " + ErrInfRel + " (denom: " + denom + ")");
                }
#endif

            } else {
                ParameterUpdate(CurrentState, Params);
                var eval = IBM_Op.GetEvaluatorEx(LsTrk, CurrentState, Params, Mapping);
                eval.time = phystime;
                eval.CellLengthScales[FluidSpecies[0]] = AgglomeratedCellLengthScales[FluidSpecies[0]];

                eval.Evaluate(1.0, 1.0, OpAffine);

            }

            m_LenScales = null;


            if (OpMatrix != null)
                OpMatrix.CheckForNanOrInfM();
            OpAffine.CheckForNanOrInfV();


            /*
            // Set Pressure Reference Point
            if (!this.boundaryCondMap.DirichletPressureBoundary) {
                if (OpMatrix != null) {
                    IBMSolverUtils.SetPressureReferencePoint(
                        CurrentSolution.Mapping,
                        this.GridData.SpatialDimension,
                        this.LsTrk,
                        OpMatrix, OpAffine);
                    //OpMatrix.SaveToTextFileSparse("OpMatrix_3D");
                } else {
                    IBMSolverUtils.SetPressureReferencePointResidual(
                        new CoordinateVector(CurrentState),
                        this.GridData.SpatialDimension,
                        this.LsTrk,
                        OpAffine);
                }
            }
            */
        }

        public virtual double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {

            //this.LevSet.ProjectField(X => this.Control.Ph(X, phystime + dt));
            //this.LsTrk.UpdateTracker(incremental: true);

            //LevsetEvo(phystime, dt, null);

            return 0.0;
        }


        //protected TextWriter Log_DragAndLift,Log_DragAndLift_P1;
        protected double[] Test_Force = new double[3];
        protected double torque = new double();
        protected double oldtorque = new double();

        private int AnalyseCounter=1;

        //SinglePhaseField blocking = null;

        /// <summary>
        /// Depending on settings <see cref="IBM_Control.Option_Timestepper"/>, computes either one timestep or a steady-state solution.
        /// </summary>
        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using (new FuncTrace()) {

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;

                base.ResLogger.TimeStep = TimestepInt;

                dt = base.GetFixedTimestep();

                Console.WriteLine("In-stationary solve, time-step #{0}, dt = {1} ...", TimestepNo, dt);

                m_BDF_Timestepper.Solve(phystime, dt); 
                

                // Residual();
                this.ResLogger.NextTimestep(false);

                // L2 error against exact solution
                // ===============================
                this.ComputeL2Error();

                #region Get Drag and Lift Coefficiant
                /*
                if (phystime == 0 && Log_DragAndLift==null) {
                    if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                        Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                        string firstline;
                        if (this.GridData.SpatialDimension == 3) {
                            firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#Timestep", "#Time", "x-Force", "y-Force", "z-Force");
                        } else {
                            firstline = String.Format("{0}\t{1}\t{2}\t{3}", "#Timestep", "#Time", "x-Force", "y-Force");
                        }
                        Log_DragAndLift.WriteLine(firstline);
                    }
                }
                */

                Test_Force = IBMSolverUtils.GetForces(Velocity, Pressure, this.LsTrk, this.Control.PhysicalParameters.mu_A/this.Control.PhysicalParameters.rho_A);
                //oldtorque = torque;
                torque = IBMSolverUtils.GetTorque(Velocity, Pressure, this.LsTrk, this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A, this.Control.particleRadius);

                /*
                if ((base.MPIRank == 0) && (Log_DragAndLift != null)) {
                    string line;
                    if (this.GridData.SpatialDimension == 3) {
                        line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, force[0], force[1], force[2]);
                    } else {
                        line = String.Format("{0}\t{1}\t{2}\t{3}", TimestepNo, phystime, force[0], force[1]);
                    }
                    Log_DragAndLift.WriteLine(line);
                    Log_DragAndLift.Flush();
                }
                */
                
                Console.WriteLine("x-Force:   {0}", Test_Force[0]);
                Console.WriteLine("y-Force:   {0}", Test_Force[1]);
                if (this.GridData.SpatialDimension == 3)
                    Console.WriteLine("z-Force:   {0}", Test_Force[2]);
                Console.WriteLine("Torqe:   {0}", torque);
                Console.WriteLine();


                // Save for NUnit Test
                base.QueryHandler.ValueQuery("C_Drag", 2 * Test_Force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                base.QueryHandler.ValueQuery("C_Lift", 2 * Test_Force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
                #endregion

                return dt;
            }
        }

        /// <summary>
        /// Variable mapping of test and trial-space.
        /// </summary>
        protected UnsetteledCoordinateMapping SaddlePointProblemMapping {
            get {
                return this.CurrentResidual.Mapping;
            }
        }

        protected double hack_phystime;

        /// <summary>
        /// Species which represents the flow domain.
        /// </summary>
        protected SpeciesId[] FluidSpecies {
            get {
                return new SpeciesId[] { LsTrk.GetSpeciesId("A") }; // wir rechnen nur species A
            }
        }

        MultidimensionalArray m_LenScales;


        /// <summary>
        /// Custom Function to compute penalty factor for viscous terms at the immersed boundary
        /// </summary>
        /// <param name="jCellIn"></param>
        /// <param name="jCellOut"></param>
        /// <param name="cj"></param>
        /// <param name="penalty">base factor</param>
        /// <returns></returns>
        protected double ComputePenaltyIB(double penalty_base, int jCell) {
            var __gridData = (GridData)GridData;

            double hCutCellMin = m_LenScales[jCell]; // for IBM, there is no positive species!
            double hCellMin = __gridData.Cells.h_min[jCell];
            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            double µ = penalty_base / hCutCellMin;
            if (double.IsNaN(µ) || double.IsInfinity(µ))
                throw new ArgumentOutOfRangeException("Invalid penalty parameter");
            Debug.Assert(!(double.IsNaN(µ) || double.IsInfinity(µ)));
            return µ;
        }


        /// <summary>
        /// Computes average velocity in case of Navier-Stokes Equations
        /// </summary>
        /// <param name="U0"></param>
        /// <param name="U0mean"></param>
        private void ComputeAverageU(VectorField<SinglePhaseField> U0, VectorField<SinglePhaseField> U0mean) {
            using (FuncTrace ft = new FuncTrace()) {
                var CC = this.LsTrk.Regions.GetCutCellMask();
                int D = this.LsTrk.GridDat.SpatialDimension;
                double minvol = Math.Pow(this.LsTrk.GridDat.Cells.h_minGlobal, D);

                U0mean.Clear();

                int QuadDegree = this.HMForder;

                //var qh = new XQuadSchemeHelper(LsTrk, momentFittingVariant, this.FluidSpecies);
                var qh = LsTrk.GetXDGSpaceMetrics(this.FluidSpecies, QuadDegree, 1).XQuadSchemeHelper;
                foreach (var Spc in this.FluidSpecies) { // loop over species...
                    //var Spc = this.LsTrk.GetSpeciesId("B"); {
                    // shadow fields
                    var U0_Spc = U0.ToArray();
                    var U0mean_Spc = U0mean.ToArray();


                    // normal cells:
                    for (int d = 0; d < D; d++) {
                        U0mean_Spc[d].AccLaidBack(1.0, U0_Spc[d], this.LsTrk.Regions.GetSpeciesMask(Spc));
                    }


                    // cut cells
                    var scheme = qh.GetVolumeQuadScheme(Spc, IntegrationDomain: this.LsTrk.Regions.GetCutCellMask());

                    var rule = scheme.Compile(this.LsTrk.GridDat, QuadDegree);
                    CellQuadrature.GetQuadrature(new int[] { D + 1 }, // vector components: ( avg_vel[0], ... , avg_vel[D-1], cell_volume )
                        this.LsTrk.GridDat,
                        rule,
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            EvalResult.Clear();
                            for (int d = 0; d < D; d++)
                                U0_Spc[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                            var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, D);
                            Vol.SetAll(1.0);
                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            for (int i = 0; i < Length; i++) {
                                int jCell = i + i0;

                                double Volume = ResultsOfIntegration[i, D];
                                if (Math.Abs(Volume) < minvol * 1.0e-12) {
                                    // keep current value
                                    // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                                } else {
                                    for (int d = 0; d < D; d++) {
                                        double IntVal = ResultsOfIntegration[i, d];
                                        U0mean_Spc[d].SetMeanValue(jCell, IntVal / Volume);
                                    }
                                }

                            }
                        }).Execute();

                }

#if DEBUG
                {
                    var Uncut = LsTrk.Regions.GetCutCellMask().Complement();


                    VectorField<SinglePhaseField> U0mean_check = new VectorField<SinglePhaseField>(D, new Basis(LsTrk.GridDat, 0), SinglePhaseField.Factory);
                    for (int d = 0; d < D; d++) {
                        U0mean_check[d].ProjectField(1.0, U0[d].Evaluate,
                            new CellQuadratureScheme(false, Uncut).AddFixedOrderRules(LsTrk.GridDat, U0[d].Basis.Degree + 1));
                    }

                    foreach (var _Spc in this.LsTrk.SpeciesIdS) { // loop over species...
                        for (int d = 0; d < D; d++) {
                            U0mean_check[d].AccLaidBack(-1.0, U0mean[d]);
                        }
                    }

                    double checkNorm = U0mean_check.L2Norm();
                    //Debug.Assert(checkNorm < 1.0e-6);
                }
#endif


                U0mean.ForEach(F => F.CheckForNanOrInf(true, true, true));

            }
        }

        /// <summary>
        /// Tecplot output.
        /// </summary>
        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling) {
            Tecplot.PlotFields(m_RegisteredFields, "IBM_Solver" + timestepNo, physTime, superSampling);
        }

        /// <summary>
        /// DG field instantiation.
        /// </summary>
        protected override void CreateFields() {
            using (new FuncTrace()) {
                base.CreateFields();
                base.LsTrk = this.LevsetTracker;
                if (Control.CutCellQuadratureType != this.LevsetTracker.CutCellQuadratureType)
                    throw new ApplicationException();
                //if (this.Control.LevelSetSmoothing) {
                //    SmoothedLevelSet = new SpecFemField(new SpecFemBasis((GridData)LevSet.GridDat, LevSet.Basis.Degree + 1));
                //}
            }
        }

        /// <summary>
        /// Setting initial values.
        /// </summary>
        protected override void SetInitial() {

            if (true) {
                DGField mpiRank = new SinglePhaseField(new Basis(GridData, 0), "rank");
                m_IOFields.Add(mpiRank);

                for (int j = 0; j < GridData.iLogicalCells.NoOfLocalUpdatedCells; j++) {
                    mpiRank.SetMeanValue(j, DatabaseDriver.MyRank);
                }

                Console.WriteLine("Total number of cells:    {0}", Grid.NumberOfCells);
                Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Count());
                Console.WriteLine("Total number of cut cells:     {0}", LsTrk.Regions.GetCutCellMask().NoOfItemsLocally);
            }

            // Using defauls CellCostEstimateFactories          
            if (this.Control.DynamicLoadBalancing_CellCostEstimatorFactories.Count == 0) {
                Console.WriteLine("Using standard CellCostEstimatorFactories");
                Control.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication app, int noOfPerformanceClasses) {
                    Console.WriteLine("i was called");
                    int[] map = new int[] { 1, 1, 10 };
                    return new StaticCellCostEstimator(map); 
                });
                Control.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication app, int noOfPerformanceClasses) {
                    Console.WriteLine("i was called");
                    int[] map = new int[] { 1, 10, 1 };
                    return new StaticCellCostEstimator(map);
                });
                Control.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication app, int noOfPerformanceClasses) {
                    Console.WriteLine("i was called");
                    int[] map = new int[] { 10, 1, 1 };
                    return new StaticCellCostEstimator(map);
                });
            }

            // Set particle radius for exact circle integration
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.ExactCircle)
                BoSSS.Foundation.XDG.Quadrature.HMF.ExactCircleLevelSetIntegration.RADIUS = new double[] { this.Control.particleRadius };
            
            Console.WriteLine("Total number of cells:    {0}", Grid.NumberOfCells);
            Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Count().MPISum());
            base.SetInitial();

            this.LevSet.GetExtremalValues(out double LevsetMin, out double LevsetMax);
            if (LevsetMax == 0.0 && LevsetMin == 0.0) {
                // User probably does not want to use Levelset, but forgot to set it.
                LevSet.AccConstant(-1.0);
            }

            /*
            PerformLevelSetSmoothing(LsTrk.Regions.GetCutCellMask(),
                LsTrk.Regions.GetSpeciesMask("B").Except(LsTrk.Regions.GetCutCellMask()),
                false);
            LsTrk.UpdateTracker(0.0);
            */


            // =======================OUTPUT FOR GMRES=====================================
            //if(this.MPISize == 1) {
            //    Console.WriteLine("!!!GMRES solver stats are saved in .txt file!!!");
            //    if(this.Control.savetodb) {
            //        SessionPath = this.Control.DbPath + "\\sessions\\" + this.CurrentSessionInfo.ID;
            //        using(StreamWriter writer = new StreamWriter(SessionPath + "\\GMRES_Stats.txt", true)) {
            //            writer.WriteLine("#GMRESIter" + "   " + "error");
            //        }
            //    } else {
            //        SessionPath = Directory.GetCurrentDirectory();
            //        if(File.Exists("GMRES_Stats.txt")) {
            //            File.Delete("GMRES_Stats.txt");
            //        }
            //        using(StreamWriter writer = new StreamWriter("GMRES_Stats.txt", true)) {
            //            writer.WriteLine("#GMRESIter" + "   " + "error");
            //        }
            //    }
            //}

            CreateEquationsAndSolvers(null);
            After_SetInitialOrLoadRestart(0.0);
            m_BDF_Timestepper.SingleInit();
        }

        /// <summary>
        /// delegate for the initialization of previous timesteps from restart session
        /// </summary>
        /// <param name="TimestepIndex"></param>
        /// <param name="time"></param>
        /// <param name="St"></param>
        private void BDFDelayedInitLoadRestart(int TimestepIndex, double time, DGField[] St) {

            Console.WriteLine("Timestep index {0}, time {1} ", TimestepIndex, time);

            ITimestepInfo tsi_toLoad;
            if (TimestepIndex < 0) {
                throw new ArgumentOutOfRangeException("Not enough Timesteps to restart with desired Timestepper");
            } else {
                ISessionInfo reloadSession = GetDatabase().Controller.GetSessionInfo(this.CurrentSessionInfo.RestartedFrom);
                tsi_toLoad = reloadSession.Timesteps.Single(t => t.TimeStepNumber.Equals(new TimestepNumber(TimestepIndex)));
            }
            DatabaseDriver.LoadFieldData(tsi_toLoad, this.GridData, this.IOFields);

            // level-set
            // ---------
            this.DGLevSet.Current.Clear();
            this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);

            this.LsTrk.UpdateTracker(time, incremental: true);

            // solution
            // --------
            int D = this.LsTrk.GridDat.SpatialDimension;

            for (int d = 0; d < D; d++) {
                St[d] = this.Velocity[d].CloneAs();
            }
            St[D] = this.Pressure.CloneAs();
        }

        private void After_SetInitialOrLoadRestart(double time) {
            using (new FuncTrace()) {
                int D = this.GridData.SpatialDimension;
                
                // we only save 'LevSet', but not the 'DGLevSet'
                // therefore, after re-start we have to copy LevSet->DGLevSet
                this.DGLevSet.Current.Clear();
                this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);
             
                
                // we push the current state of the level-set, so we have an initial value
                this.LsTrk.UpdateTracker(time);
                this.DGLevSet.IncreaseHistoryLength(1);
                this.LsTrk.PushStacks();
                this.DGLevSet.Push();

            }
        }

        /// <summary>
        /// Ensures that the level-set field <see cref="LevSet"/> is continuous, if <see cref="IBM_Control.LevelSetSmoothing"/> is true. Note that this is not necessary if the order of the level-set function of the particles is equal to the polynomial DG order.
        /// </summary>
        protected void PerformLevelSetSmoothing(CellMask domain, CellMask NegMask, bool SetFarField) {

            if (this.Control.LevelSetSmoothing) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // smoothing on: perform some kind of C0-projection
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                var ContinuityEnforcer = new BoSSS.Solution.LevelSetTools.ContinuityProjection(
                    ContBasis: this.LevSet.Basis,
                    DGBasis: this.DGLevSet.Current.Basis,
                    gridData: GridData,
                    Option: Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG);

                //CellMask domain = this.LsTrk.Regions.GetNearFieldMask(1);

                ContinuityEnforcer.MakeContinuous(DGLevSet.Current, LevSet, domain, null, false);
                if (SetFarField)
                {
                    LevSet.Clear(NegMask);
                    LevSet.AccConstant(-1, NegMask);
                }
            } else {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // no smoothing (not recommended): copy DGLevSet -> LevSet
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                //this.LevSet.Clear(domain);
                //this.LevSet.AccLaidBack(1.0, this.DGLevSet.Current, domain);
                this.LevSet.Clear();
                this.LevSet.AccLaidBack(1.0, this.DGLevSet.Current);
            }
        }


        /// <summary>
        /// BDF timestepper init after restart
        /// </summary>
        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            base.LoadRestart(out Time, out TimestepNo);
            this.CreateEquationsAndSolvers(null);

            // =========================================
            // XDG Timestepper initialization
            // =========================================

            if (this.Control.TimeStepper_Init == Solution.Timestepping.TimeStepperInit.MultiInit) {
                // =========================================
                // XDG BDF Timestepper initialization
                // =========================================

                if (m_BDF_Timestepper != null) {
                    m_BDF_Timestepper.DelayedTimestepperInit(Time, TimestepNo.MajorNumber, this.Control.GetFixedTimestep(),
                        // delegate for the initialization of previous timesteps from restart session
                        BDFDelayedInitLoadRestart);
                }

                After_SetInitialOrLoadRestart(Time);
            } else {
                if (m_BDF_Timestepper != null) {
                    After_SetInitialOrLoadRestart(Time);
                    m_BDF_Timestepper.SingleInit();
                }
            }

        }


        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        protected MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int pVel = this.Velocity[0].Basis.Degree;
                int pPrs = this.Pressure.Basis.Degree;
                int D = this.GridData.SpatialDimension;

                if (this.Control.VelocityBlockPrecondMode != MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite
                    && this.Control.VelocityBlockPrecondMode != MultigridOperator.Mode.IdMass_DropIndefinite) {
                    throw new NotSupportedException("Invalid option for block-preconditioning of momentum equation: " + this.Control.VelocityBlockPrecondMode
                        + ". Valid options are " + MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite + " and " + MultigridOperator.Mode.IdMass_DropIndefinite + ".");

                }


                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[D + 1];

                    // configurations for velocity
                    for (int d = 0; d < D; d++) {
                        configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                            mode = this.Control.VelocityBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                        mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                        VarIndex = new int[] { D }
                    };
                }


                return configs;
            }
        }

        /// <summary>
        /// L2 norm of current solution against <see cref="IBM_Control.ExSol_Velocity"/> resp. <see cref="IBM_Control.ExSol_Pressure"/>.
        /// </summary>
        protected void ComputeL2Error() {
            if (this.Control.ExSol_Velocity_Evaluator == null && this.Control.ExSol_Pressure_Evaluator == null)
                // nothing to do
                return;


            int D = this.GridData.SpatialDimension;
            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(FluidSpecies, order, 1).XQuadSchemeHelper;

            // Velocity error
            // ==============
            if (this.Control.ExSol_Velocity_Evaluator != null) {
                double[] L2Error = new double[D];

                var spId = this.FluidSpecies.Single();


                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                for (int d = 0; d < D; d++) {
                    L2Error[d] = this.Velocity[d].L2Error(this.Control.ExSol_Velocity_Evaluator[d].Vectorize(0.0), order, scheme);
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                }
            }


            // pressure error
            // ==============
            if (this.Control.ExSol_Pressure_Evaluator != null) {

                // pass 1: mean value of pressure difference
                double DiffInt = 0;
                foreach (var spId in FluidSpecies) {

                    string spc = this.LsTrk.GetSpeciesName(spId);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    DiffInt += this.Pressure.LxError(this.Control.ExSol_Pressure_Evaluator.Vectorize(0.0), (X, a, b) => (a - b), rule);
                    //Volume +=  this.Pressure.GetSpeciesShadowField(spc).LxError(null, (a, b) => (1.0), rule);
                }
                double Volume2 = (new SubGrid(CellMask.GetFullMask(this.GridData))).Volume;
                double PressureDiffMean = DiffInt / Volume2;


                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach (var spId in this.FluidSpecies) {

                    //SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    string spc = this.LsTrk.GetSpeciesName(spId);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    double IdV = this.Pressure.LxError(this.Control.ExSol_Pressure_Evaluator.Vectorize(0.0), (X, a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, IdV.Sqrt());

                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
                }


                L2Error = L2Error.Sqrt();
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
            }
        }

        /// <summary>
        /// Cell-performance classes:
        /// - void cells are 0
        /// - non-cut fluid cells are 1
        /// - cut cells are 2
        /// </summary>
        protected override void GetCellPerformanceClasses(out int NoOfClasses, out int[] CellPerfomanceClasses, int TimeStepNo, double physTime) {
            NoOfClasses = 3;
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            CellPerfomanceClasses = new int[J];
            foreach (int j in LsTrk.Regions.GetSpeciesMask("B").ItemEnum)
                CellPerfomanceClasses[j] = 0;
            foreach (int j in LsTrk.Regions.GetSpeciesMask("A").ItemEnum)
                CellPerfomanceClasses[j] = 1;
            foreach (int j in LsTrk.Regions.GetCutCellMask().ItemEnum)
                CellPerfomanceClasses[j] = 2;
        }

        public override void PostRestart(double time, TimestepNumber timestep) {

            // Find path to PhysicalData.txt
            try {
                var fsDriver = this.DatabaseDriver.FsDriver;
                string pathToOldSessionDir = System.IO.Path.Combine(
                    fsDriver.BasePath, "sessions", this.CurrentSessionInfo.RestartedFrom.ToString());
                string pathToPhysicalData = System.IO.Path.Combine(pathToOldSessionDir, "PhysicalData.txt");
                string[] records = File.ReadAllLines(pathToPhysicalData);

                string line1 = File.ReadLines(pathToPhysicalData).Skip(1).Take(1).First();
                string line2 = File.ReadLines(pathToPhysicalData).Skip(2).Take(1).First();
                string[] fields_line1 = line1.Split('\t');
                string[] fields_line2 = line2.Split('\t');

                double dt = Convert.ToDouble(fields_line2[1]) - Convert.ToDouble(fields_line1[1]);
            } catch (FileNotFoundException) {
                Console.WriteLine("PhysicalData.txt could not be found! Assuming we start with timestep #0 ...");
            }
            //int idx_restartLine = Convert.ToInt32(time / dt + 1.0);
            //string restartLine = File.ReadLines(pathToPhysicalData).Skip(idx_restartLine - 1).Take(1).First();
            //double[] values = Array.ConvertAll<string, double>(restartLine.Split('\t'), double.Parse);

            /* string restartLine = "";
              Calculcation of dt 
             var physicalData = File.ReadLines(pathToPhysicalData);
             int count = 0;
             foreach (string line in physicalData)
             {
                 string[] fields = line.Split('\t');
                 restartLine = line;
                 if (count != 0) { 
                 if (Convert.ToDouble(fields[1]) > time)
                 {
                     break;
                 }
             }
             count++;
             }


             double dt = Convert.ToDouble(fields_line2[1]) - Convert.ToDouble(fields_line1[1]);

              Using dt to find line of restart time
             int idx_restartLine = Convert.ToInt32(time / dt + 1.0);
             string restartLine = File.ReadLines(pathToPhysicalData).Skip(idx_restartLine - 1).Take(1).First();
             double[] values = Array.ConvertAll<string, double>(restartLine.Split('\t'), double.Parse);*/

            //Adding PhysicalData.txt
            /*
            if ((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                Log_DragAndLift = base.DatabaseDriver.FsDriver.GetNewLog("PhysicalData", CurrentSessionInfo.ID);
                string firstline;
                if (this.GridData.SpatialDimension == 3) {
                    firstline = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "#Timestep", "#Time", "x-Force", "y-Force", "z-Force");
                } else {
                    firstline = String.Format("{0}\t{1}\t{2}\t{3}", "#Timestep", "#Time", "x-Force", "y-Force");
                }
                Log_DragAndLift.WriteLine(firstline);
                //Log_DragAndLift.WriteLine(restartLine);
            }
            */
        }

        protected override void Bye() {
            /*
            if (Log_DragAndLift != null) {
                try {
                    Log_DragAndLift.Flush();
                    Log_DragAndLift.Close();
                    Log_DragAndLift.Dispose();
                } catch (Exception) { }
                Log_DragAndLift = null;
            }

            if(Log_DragAndLift_P1 != null) {
                try {
                    Log_DragAndLift_P1.Flush();
                    Log_DragAndLift_P1.Close();
                    Log_DragAndLift_P1.Dispose();
                } catch (Exception) { }
                Log_DragAndLift_P1 = null;
            }
            */
        }

        /// <summary>
        /// Attention: SENSITIVE TO LEVEL INDICATOR
        /// </summary>
        readonly bool debug = true;

        /// <summary>
        /// Very primitive refinement indicator, works on a LevelSet criterion.
        /// </summary>
        /// 
        int LevelIndicator(int j, int CurrentLevel)
        {
            var LevSetCells = LsTrk.Regions.GetCutCellMask();
            var LevSetNeighbours = LsTrk.Regions.GetNearFieldMask(1);
            int DesiredLevel_j = 0;

            if (!debug) {
                if (LevSetCells.Contains(j))
                    DesiredLevel_j = 1;
            } else {
                if (LevSetCells.Contains(j)) {
                    DesiredLevel_j = 2;
                    Console.WriteLine(" ich störe");
                } else
                    if (LevSetNeighbours.Contains(j)) { DesiredLevel_j = 2; }
            }

            return DesiredLevel_j;
        }

        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {

            if (this.Control.AdaptiveMeshRefinement) {

                //if (TimestepNo > 3 && TimestepNo % 3 != 0) {
                //    newGrid = null;
                //    old2NewGrid = null;
                //    return;
                //}

                // Check grid changes
                // ==================

                //// compute curvature for levelindicator 
                //CurvatureAlgorithms.CurvatureDriver(
                //SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //CurvatureAlgorithms.FilterConfiguration.Default,
                //this.Curvature, out VectorField<SinglePhaseField> LevSetGradient, this.LsTrk,
                //this.HMForder, this.DGLevSet.Current);

                CellMask CutCells = LsTrk.Regions.GetCutCellMask();
                //CellMask CutCellNeighbors = LsTrk.Regions.GetNearFieldMask(1);
                //var CutCellArray = CutCells.ItemEnum.ToArray();
                //var CutCellNeighborsArray = CutCellNeighbors.ItemEnum.ToArray();
                //var AllCells = CutCellArray.Concat(CutCellNeighborsArray).ToArray();
                //var NoCoarseningcells = new CellMask(this.GridData, AllCells);

                // Only CutCells are NoCoarseningCells 
                GridRefinementController gridRefinementController = new GridRefinementController((GridData)(this.GridData), CutCells);
                bool AnyChange = gridRefinementController.ComputeGridChange(LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if (AnyChange) {
                    int[] glb = (new int[] {
                    CellsToRefineList.Count,
                    Coarsening.Sum(L => L.Length),
                }).MPISum();
                    NoOfCellsToRefine = glb[0];
                    NoOfCellsToCoarsen = glb[1];
                }
                int oldJ = this.GridData.CellPartitioning.TotalLength;

                // Update Grid
                // ===========

                if (AnyChange) {

                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                    newGrid = ((GridData)(this.GridData)).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                    if (this.Control.savetodb == true) {
                        //Console.WriteLine("Save adaptive Mesh...");
                        //Console.WriteLine("GridGUID:   " + newGrid.GridGuid);
                        //DatabaseDriver.SaveGrid(newGrid, base.GetDatabase());
                        //Console.WriteLine("...done");
                    }
                } else {

                    Console.WriteLine("No changes in Grid");
                    newGrid = null;
                    old2NewGrid = null;
                }

                //debug = false;

            } else {

                newGrid = null;
                old2NewGrid = null;
            }
        }

    }
}
