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
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.SpecFEM;
using MPI.Wrappers;
using BoSSS.Foundation.IO;
using System.Diagnostics;
using System.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using NUnit.Framework;
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Application.IBM_Solver {

    /// <summary>
    /// A incompressible Navier-Stokes solver with the possibility of using non moving immersed boundaries.
    /// </summary>
    public class IBM_SolverMain : Application<IBM_Control> {


        /// <summary>
        /// Application entry point.
        /// </summary>
        static void Main(string[] args) {
         
            BoSSS.Solution.Application<IBM_Control>._Main(args, false, delegate () {
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

        IDictionary<SpeciesId, IEnumerable<double>> Rho {
            get {
                double rho = this.Control.PhysicalParameters.rho_A;

                int D = this.GridData.SpatialDimension;

                double[] _rho = new double[D];
                _rho.SetAll(rho);

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho);
                return R;
            }
        }

        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        virtual protected IDictionary<SpeciesId, IEnumerable<double>> MassScale {
            get {
                double rho = this.Control.PhysicalParameters.rho_A;

                int D = this.GridData.SpatialDimension;

                double[] _rho = new double[D + 1];
                _rho.SetAll(rho);
                //No MassMatrix for the pressure
                _rho[D] = 0;

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.LsTrk.GetSpeciesId("A"), _rho);

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
        protected IncompressibleBoundaryCondMap BcMap;

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
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        protected MsrMatrix SaddlePointMatrix;

        /// <summary>
        /// output of <see cref="AssembleMatrix"/>;
        /// </summary>
        protected double[] SaddlePointRHS;

        protected bool U0MeanRequired = false;


        protected XSpatialOperator IBM_Op;

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

        protected XdgBDFTimestepping m_BDF_Timestepper;

        //SinglePhaseField[] MGColoring;

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
                BcMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

                int degU = this.Velocity[0].Basis.Degree;
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
                var Params = ArrayTools.Cat(
                     VariableNames.Velocity0Vector(D),
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

                IBM_Op = new XSpatialOperator(DomNameSelected, Params, CodNameSelected,
                    (A, B, C) => this.HMForder);
                
                // Momentum equation
                // =================

                // convective part:
                if (IBM_Op_config.convection) {
                    for (int d = 0; d < D; d++) {

                        var comps = IBM_Op.EquationComponents[CodName[d]];

                        //var ConvBulk = new Solution.XNSECommon.Operator.Convection.ConvectionInBulk_LLF(D, BcMap, d, this.Control.PhysicalParameters.rho_A, 1, IBM_Op_config.dntParams.LFFA, IBM_Op_config.dntParams.LFFB, LsTrk);
                        var ConvBulk = new Solution.NSECommon.LinearizedConvection(D, BcMap, d);
                        //IBM_Op.OnIntegratingBulk += ConvBulk.SetParameter;
                        comps.Add(ConvBulk); // bulk component

                        var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ConvectionAtIB(d, D, LsTrk, this.Control.AdvancedDiscretizationOptions.LFFA, BcMap,
                            delegate (double[] X, double time) { return new double[] { 0.0, 0.0, 0.0, 0.0 }; }, this.Control.PhysicalParameters.rho_A, false);

                        comps.Add(ConvIB); // immersed boundary component
                    }
                    this.U0MeanRequired = true;
                }

                // pressure part:
                if (IBM_Op_config.PressureGradient) {
                    for (int d = 0; d < D; d++) {
                        var comps = IBM_Op.EquationComponents[CodName[d]];
                        var pres = new PressureGradientLin_d(d, BcMap);
                        //var pres = new Solution.XNSECommon.Operator.Pressure.PressureInBulk(d, BcMap, 1, 1);
                        //IBM_Op.OnIntegratingBulk += pres.SetParameter;
                        comps.Add(pres); // bulk component

                        var presLs = new BoSSS.Solution.NSECommon.Operator.Pressure.PressureFormAtIB(d, D, LsTrk);
                        comps.Add(presLs); // immersed boundary component


                        // if periodic boundary conditions are applied a fixed pressure gradient drives the flow
                        if (this.Control.FixedStreamwisePeriodicBC) {
                            var presSource = new SrcPressureGradientLin_d(this.Control.SrcPressureGrad[d]);
                            comps.Add(presSource);
                        }

                    }
                }

                // viscous part:
                if (IBM_Op_config.Viscous) {
                    for (int d = 0; d < D; d++) {
                        var comps = IBM_Op.EquationComponents[CodName[d]];
                        double _D = D;
                        double penalty_mul = this.Control.AdvancedDiscretizationOptions.PenaltySafety;
                        double _p = degU;
                        double penalty_base = (_p + 1) * (_p + _D) / D;
                        double penalty = penalty_base * penalty_mul;
                        double penalty_bulk = this.Control.AdvancedDiscretizationOptions.PenaltySafety;


                        //var Visc = new Solution.XNSECommon.Operator.Viscosity.ViscosityInBulk_GradUTerm(penalty, 1.0, BcMap, d, D, this.Control.PhysicalParameters.mu_A, 1, ViscosityImplementation.H);
                        var Visc = new swipViscosity_Term1(penalty_bulk, d, D, BcMap,
                            ViscosityOption.ConstantViscosity,
                            this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                            double.NaN, null);
                        //delegate (double p, int i, int j, double[] cell) { return ComputePenalty(p, i, j, cell); });
                        // IBM_Op.OnIntegratingBulk += Visc.SetParameter;
                        comps.Add(Visc); // bulk component GradUTerm 
                        var ViscLs = new BoSSS.Solution.NSECommon.Operator.Viscosity.ViscosityAtIB(d, D, LsTrk,
                            penalty, this.ComputePenaltyIB,
                            this.Control.PhysicalParameters.mu_A / this.Control.PhysicalParameters.rho_A,
                            delegate (double[] X, double time) { return new double[] { 0.0, 0.0, 0.0, 0.0 }; });
                        comps.Add(ViscLs); // immersed boundary component
                    }
                }

                // Continuum equation
                // ==================
                if (IBM_Op_config.continuity) {
                    for (int d = 0; d < D; d++) {

                        //var src = new Solution.XNSECommon.Operator.Continuity.DivergenceInBulk_Volume(d, D, 1, 0, 1, false);
                        var src = new Divergence_DerivativeSource(d, D);
                        //IBM_Op.OnIntegratingBulk += src.SetParameter;
                        var flx = new Divergence_DerivativeSource_Flux(d, BcMap);
                        //IBM_Op.OnIntegratingBulk += flx.SetParameter;
                        IBM_Op.EquationComponents["div"].Add(src);
                        IBM_Op.EquationComponents["div"].Add(flx);

                        //var presStab = new PressureStabilization(1, this.GridData.Edges.h_max_Edge, 1 / this.Control.PhysicalParameters.mu_A);
                        //IBM_Op.EquationComponents["div"].Add(presStab);
                    }

                    var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB(D, LsTrk, 1,
                        delegate (double[] X, double time) { return new double[] { 0.0, 0.0, 0.0, 0.0 }; });
                    IBM_Op.EquationComponents["div"].Add(divPen); // immersed boundary component 


                    //IBM_Op.EquationComponents["div"].Add(new PressureStabilization(1, 1.0 / this.Control.PhysicalParameters.mu_A));
                }
                IBM_Op.Commit();
            }

            // ==========================
            // create timestepper
            // ==========================

            var Unknowns = ArrayTools.Cat(this.Velocity, this.Pressure);
            var Residual = ArrayTools.Cat(this.ResidualMomentum, this.ResidualContinuity);

            if (L == null) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                // Creation of time-integrator (initial, no balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++

                if (m_BDF_Timestepper == null) {
                    LevelSetHandling lsh = LevelSetHandling.None;
                    SpatialOperatorType SpatialOp = SpatialOperatorType.LinearTimeDependent;

                    if (this.Control.PhysicalParameters.IncludeConvection) {
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

                    m_BDF_Timestepper = new XdgBDFTimestepping(
                        Unknowns, Residual,
                        LsTrk, true,
                        DelComputeOperatorMatrix, null, DelUpdateLevelset,
                        bdfOrder,
                        lsh,
                        MassMatrixShapeandDependence.IsTimeDependent,
                        SpatialOp,
                        MassScale,
                        this.MultigridOperatorConfig, base.MultigridSequence,
                        this.FluidSpecies, this.HMForder,
                        this.Control.AdvancedDiscretizationOptions.CellAgglomerationThreshold,
                        false, this.Control.NonLinearSolver, this.Control.LinearSolver
                        );

                    m_BDF_Timestepper.m_ResLogger = base.ResLogger;
                    m_BDF_Timestepper.m_ResidualNames = ArrayTools.Cat(this.ResidualMomentum.Select(f => f.Identification), this.ResidualContinuity.Identification);
            
                    m_BDF_Timestepper.SessionPath = SessionPath;
                    m_BDF_Timestepper.Timestepper_Init = Solution.Timestepping.TimeStepperInit.MultiInit;
                }

            } else {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // restore BDF time-stepper after grid redistribution (dynamic load balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                m_BDF_Timestepper.DataRestoreAfterBalancing(L, Unknowns, Residual, this.LsTrk, base.MultigridSequence);
            }

            Debug.Assert(m_BDF_Timestepper != null);
        }

        public override void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            m_BDF_Timestepper.DataBackupBeforeBalancing(L);
            m_CurrentResidual = null;
            m_CurrentSolution = null;
            IBM_Op = null;
        }


        int DelComputeOperatorMatrix_CallCounter = 0;

        protected virtual void DelComputeOperatorMatrix(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime) {
            DelComputeOperatorMatrix_CallCounter++;

            // compute operator
            //Debug.Assert(OpMatrix.InfNorm() == 0.0);
            //Debug.Assert(OpAffine.L2Norm() == 0.0);

            // parameters...
            int D = this.LsTrk.GridDat.SpatialDimension;
            var U0 = new VectorField<SinglePhaseField>(CurrentState.Take(D).Select(F => (SinglePhaseField)F).ToArray());
            SinglePhaseField[] U0_U0mean;
            if (this.U0MeanRequired) {
                Basis U0meanBasis = new Basis(GridData, 0);
                VectorField<SinglePhaseField> U0mean = new VectorField<SinglePhaseField>(D, U0meanBasis, "U0mean_", SinglePhaseField.Factory);
                U0mean.Clear();

                if (this.Control.PhysicalParameters.IncludeConvection)
                    ComputeAverageU(U0, U0mean);

                U0_U0mean = ArrayTools.Cat<SinglePhaseField>(U0, U0mean);
            } else {
                U0_U0mean = new SinglePhaseField[2 * D];
            }
            var Params = ArrayTools.Cat<SinglePhaseField>(
                U0_U0mean);

            m_LenScales = AgglomeratedCellLengthScales[FluidSpecies[0]];

            // create matrix and affine vector:
            if (OpMatrix != null) {
                //IBM_Op.ComputeMatrixEx(LsTrk,
                //    Mapping, Params, Mapping,
                //    OpMatrix, OpAffine, false, phystime, true,
                //    AgglomeratedCellLengthScales,
                //    FluidSpecies);

                var mtxBuilder = IBM_Op.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping, FluidSpecies);
                mtxBuilder.time = phystime;
                mtxBuilder.SpeciesOperatorCoefficients[FluidSpecies[0]].CellLengthScales = AgglomeratedCellLengthScales[FluidSpecies[0]];

                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine);

            } else {
                var eval = IBM_Op.GetEvaluatorEx(LsTrk, CurrentState, Params, Mapping, FluidSpecies);
                eval.time = phystime;
                eval.SpeciesOperatorCoefficients[FluidSpecies[0]].CellLengthScales = AgglomeratedCellLengthScales[FluidSpecies[0]];

                eval.Evaluate(1.0, 1.0, OpAffine);

#if DEBUG
                // remark: remove this piece in a few months from now on (09may18) if no problems occur
                {
                    var check = IBM_Op.GetMatrixBuilder(LsTrk, Mapping, Params, Mapping, FluidSpecies);
                    check.time = phystime;
                    check.SpeciesOperatorCoefficients[FluidSpecies[0]].CellLengthScales = AgglomeratedCellLengthScales[FluidSpecies[0]];

                    BlockMsrMatrix checkOpMatrix = new BlockMsrMatrix(Mapping, Mapping);
                    double[] checkAffine = new double[OpAffine.Length];
                    check.ComputeMatrix(checkOpMatrix, checkAffine);

                    double[] checkResult = checkAffine.CloneAs();
                    var currentVec = new CoordinateVector(CurrentState);
                    checkOpMatrix.SpMV(1.0, new CoordinateVector(CurrentState), 1.0, checkResult);

                    double L2_dist = GenericBlas.L2DistPow2(checkResult, OpAffine).MPISum().Sqrt();
                    double RefNorm = (new double[] { checkResult.L2NormPow2(), OpAffine.L2NormPow2(), currentVec.L2NormPow2() }).MPISum().Max().Sqrt();

                    Assert.LessOrEqual(L2_dist, RefNorm * 1.0e-6);
                    Debug.Assert(L2_dist < RefNorm * 1.0e-6);
                }
#endif
            }

            m_LenScales = null;

#if DEBUG
            if (DelComputeOperatorMatrix_CallCounter == 1 && OpMatrix != null) {
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
            if (OpMatrix != null)
                OpMatrix.CheckForNanOrInfM();
            OpAffine.CheckForNanOrInfV();



            // Set Pressure Reference Point
            if (!this.BcMap.DirichletPressureBoundary) {
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


        }



        public virtual double DelUpdateLevelset(DGField[] CurrentState, double phystime, double dt, double UnderRelax, bool incremental) {

            //this.LevSet.ProjectField(X => this.Control.Ph(X, phystime + dt));
            //this.LsTrk.UpdateTracker(incremental: true);

            //LevsetEvo(phystime, dt, null);

            return 0.0;
        }

        //public virtual CutCellMetrics DelUpdateCutCellMetrics() {
        //    return new CutCellMetrics(momentFittingVariant, this.HMForder, LsTrk, this.FluidSpecies);
        //}

        //protected TextWriter Log_DragAndLift,Log_DragAndLift_P1;
        protected double[] force = new double[3];
        protected double torque = new double();
        protected double oldtorque = new double();


        //SinglePhaseField blocking = null;

        /// <summary>
        /// Depending on settings <see cref="IBM_Control.Option_Timestepper"/>, computs either one timestep or a steady-state solution.
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

                force = IBMSolverUtils.GetForces(Velocity, Pressure, this.LsTrk, this.Control.PhysicalParameters.mu_A/this.Control.PhysicalParameters.rho_A);
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
                
                Console.WriteLine("x-Force:   {0}", force[0]);
                Console.WriteLine("y-Force:   {0}", force[1]);
                if (this.GridData.SpatialDimension == 3)
                    Console.WriteLine("z-Force:   {0}", force[2]);
                Console.WriteLine("Torqe:   {0}", torque);
                Console.WriteLine();


                // Save for NUnit Test
                base.QueryHandler.ValueQuery("C_Drag", 2 * force[0], true); // Only for Diameter 1 (TestCase NSE stationary)
                base.QueryHandler.ValueQuery("C_Lift", 2 * force[1], true); // Only for Diameter 1 (TestCase NSE stationary)
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


        /*
        /// <summary>
        /// Custom Function to compute penalty factor for viscous terms, for bulk terms
        /// </summary>
        /// <param name="jCellIn"></param>
        /// <param name="jCellOut"></param>
        /// <param name="cj"></param>
        /// <param name="penalty">base factor</param>
        /// <returns></returns>
        protected double ComputePenaltyBulk(double penalty, int jCellIn, int jCellOut, MultidimensionalArray cj) {
            double muFactor; // the WTF factor
            if (jCellOut >= 0)
                muFactor = 1.0;
            else
                //muFactor = Math.Max(1, 0) / this.Control.PhysicalParameters.mu_A; // Hardcoded for single phase flows
                muFactor = Math.Max(this.Control.PhysicalParameters.mu_A, 0) / this.Control.PhysicalParameters.mu_A; // Hardcoded for single phase flows
            double penaltySizeFactor_A = 1.0 / this.m_LenScales[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / this.m_LenScales[jCellOut] : 0;
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            //if(once <= 0) {
            //    once++;
            //    Console.WriteLine("penalty: " + penalty);
            //    Console.WriteLine("penaltySizeFactor: " + penaltySizeFactor);
            //    Console.WriteLine("muFactor: " + muFactor);
            //    Console.WriteLine("total penalty: " + (penalty * penaltySizeFactor * muFactor));
            //}

            return penalty * penaltySizeFactor * muFactor;
        }
        */

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

                ilPSP.Environment.StdoutOnlyOnRank0 = false;
                Console.WriteLine("Total number of cells:    {0}", Grid.NumberOfCells);
                Console.WriteLine("Total number of DOFs:     {0}", CurrentSolution.Count());
                Console.WriteLine("Total number of cut cells:     {0}", LsTrk.Regions.GetCutCellMask().NoOfItemsLocally);

                ilPSP.Environment.StdoutOnlyOnRank0 = true;
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

            double LevsetMin, LevsetMax;
            this.LevSet.GetExtremalValues(out LevsetMin, out LevsetMax);
            if (LevsetMax == 0.0 && LevsetMin == 0.0) {
                // User probably does not want to use Levelset, but forgot to set it.
                LevSet.AccConstant(-1.0);
            }


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
            After_SetInitialOrLoadRestart();
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

            this.LsTrk.UpdateTracker(incremental: true);

            // solution
            // --------
            int D = this.LsTrk.GridDat.SpatialDimension;

            for (int d = 0; d < D; d++) {
                St[d] = this.Velocity[d].CloneAs();
            }
            St[D] = this.Pressure.CloneAs();
        }

        private void After_SetInitialOrLoadRestart() {
            using (new FuncTrace()) {
                int D = this.GridData.SpatialDimension;
                
                // we only save 'LevSet', but not the 'DGLevSet'
                // therefore, after re-start we have to copy LevSet->DGLevSet
                this.DGLevSet.Current.Clear();
                this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);
             
                
                // we push the current state of the level-set, so we have an initial value
                this.LsTrk.UpdateTracker();
                this.DGLevSet.IncreaseHistoryLength(1);
                this.LsTrk.PushStacks();
                this.DGLevSet.Push();

            }
        }

        /// <summary>
        /// Ensures that the level-set field <see cref="LevSet"/> is continuous, if <see cref="IBM_Control.LevelSetSmoothing"/> is true
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
                    Option: Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG);

                //CellMask domain = this.LsTrk.Regions.GetNearFieldMask(1);

                ContinuityEnforcer.MakeContinuous(this.DGLevSet.Current, this.LevSet, domain, null, false);
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

                After_SetInitialOrLoadRestart();
            } else {
                if (m_BDF_Timestepper != null) {
                    After_SetInitialOrLoadRestart();
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
                            Degree = Math.Max(1, pVel - iLevel),
                            mode = this.Control.VelocityBlockPrecondMode,
                            VarIndex = new int[] { d }
                        };
                    }
                    // configuration for pressure
                    configs[iLevel][D] = new MultigridOperator.ChangeOfBasisConfig() {
                        Degree = Math.Max(0, pPrs - iLevel),
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
        bool debug = true;

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
                bool AnyChange = GridRefinementController.ComputeGridChange((GridData)(this.GridData), CutCells, LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
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
