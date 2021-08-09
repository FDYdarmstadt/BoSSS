using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.Timestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Codes for various time integration schemes of <see cref="XdgBDFTimestepping"/> and <see cref="XdgRKTimestepping"/>
    /// </summary>
    public enum TimeSteppingScheme {

        /// <summary>
        /// Explicit Euler using the <see cref="XdgBDFTimestepping"/> implementation, <see cref="BDFSchemeCoeffs.ExplicitEuler"/>
        /// </summary>
        ExplicitEuler = 0,

        /// <summary>
        /// equivalent of BDF1, using the BDF-implementation, <see cref="BDFSchemeCoeffs.BDF(int)"/>
        /// </summary>
        ImplicitEuler = 1,

        /// <summary>
        /// (Implicit) Crank-Nicholson using the BDF-implementation (<see cref="XdgBDFTimestepping"/>), <see cref="BDFSchemeCoeffs.CrankNicolson"/>
        /// </summary>
        CrankNicolson = 2000,

        /// <summary>
        /// backward differentiation formula of order 2, <see cref="BDFSchemeCoeffs.BDF(int)"/>
        /// </summary>
        BDF2 = 2,

        /// <summary>
        /// backward differentiation formula of order 3, <see cref="BDFSchemeCoeffs.BDF(int)"/>
        /// </summary>
        BDF3 = 3,

        /// <summary>
        /// backward differentiation formula of order 4, <see cref="BDFSchemeCoeffs.BDF(int)"/>
        /// </summary>
        BDF4 = 4,

        /// <summary>
        /// backward differentiation formula of order 5, <see cref="BDFSchemeCoeffs.BDF(int)"/>
        /// </summary>
        BDF5 = 5,

        /// <summary>
        /// backward differentiation formula of order 6, <see cref="BDFSchemeCoeffs.BDF(int)"/>
        /// </summary>
        BDF6 = 6,

        /// <summary>
        /// explicit Runge-Kutta, <see cref="RungeKuttaScheme.RungeKutta1901"/>
        /// </summary>
        RK4 = 104,

        /// <summary>
        /// explicit Runge-Kutta, <see cref="RungeKuttaScheme.TVD3"/>
        /// </summary>
        RK3 = 103,

        /// <summary>
        /// explicit Runge-Kutta, <see cref="RungeKuttaScheme.Heun2"/> 
        /// </summary>
        RK2 = 102,

        /// <summary>
        /// explicit Runge-Kutta, <see cref="RungeKuttaScheme.ExplicitEuler"/>
        /// </summary>
        RK1 = 101,

        /// <summary>
        /// explicit Runge-Kutta, <see cref="RungeKuttaScheme.ExplicitEuler2"/>
        /// </summary>
        RK1u1 = 110,

        /// <summary>
        /// Implicit Euler using the implicit Runge-Kutta implementation, <see cref="RungeKuttaScheme.ImplicitEuler"/> 
        /// </summary>
        RK_ImplicitEuler = 201,

        /// <summary>
        /// (Implicit) Crank-Nicholson using the implicit Runge-Kutta (<see cref="XdgRKTimestepping"/>) implementation, <see cref="RungeKuttaScheme.CrankNicolson"/>
        /// </summary>
        RK_CrankNic = 202,

        /// <summary>
        /// Implicit, <see cref="RungeKuttaScheme.IMEX3"/>
        /// </summary>
        RK_IMEX3 = 203
    }


    /// <summary>
    /// Driver class which provides a simplified interface to <see cref="XdgBDFTimestepping"/> and <see cref="XdgRKTimestepping"/>
    /// </summary>
    public class XdgTimestepping {

        public TimeSteppingScheme Scheme {
            get;
            private set;
        }

        public XSpatialOperatorMk2 XdgOperator {
            get;
            private set;
        }

        XSpatialOperatorMk2 m_JacobiXdgOperator;

        XSpatialOperatorMk2 GetJacobiXdgOperator() {
            if(m_JacobiXdgOperator == null) {
                m_JacobiXdgOperator = XdgOperator.GetJacobiOperator(GridDat.SpatialDimension) as XSpatialOperatorMk2;
            }
            return m_JacobiXdgOperator;
        }

        
        public SpatialOperator DgOperator {
            get;
            private set;
        }

        SpatialOperator m_JacobiDgOperator;

        SpatialOperator GetJacobiDgOperator() {
            if(m_JacobiDgOperator == null) {
                m_JacobiDgOperator = DgOperator.GetJacobiOperator(GridDat.SpatialDimension) as SpatialOperator;
            }
            return m_JacobiDgOperator;
        }



        public ISpatialOperator Operator {
            get {
                
                if(XdgOperator != null)
                    return XdgOperator;
                if(DgOperator != null)
                    return DgOperator;

                throw new NotImplementedException();
            }
        }


        /// <summary>
        /// <see cref="XSpatialOperatorMk2.Species"/>
        /// </summary>
        public SpeciesId[] UsedSpecies {
            get {
                if(XdgOperator != null)
                    return XdgOperator.Species.Select(spcName => LsTrk.GetSpeciesId(spcName)).ToArray();
                else
                    return null;
            }
        }

        /// <summary>
        /// Initial Value
        /// </summary>
        public CoordinateMapping CurrentState {
            get {
                return TimesteppingBase.CurrentStateMapping;
            }
            
        }

        /// <summary>
        /// Residual vector during/after linear/nonlinear solver iterations
        /// </summary>
        public CoordinateMapping IterationResiduals {
            get {
                return TimesteppingBase.Residuals.Mapping;
            }
        }

        /// <summary>
        /// internal storage for operator (<see cref="Operator"/>) parameters
        /// </summary>
        public IList<DGField> Parameters {
            get;
            private set;
        }

        /// <summary>
        /// Level Set Tracker; if a standard DG operator is used, 
        /// internally a 'dummy tracker' is initiated
        /// </summary>
        public LevelSetTracker LsTrk {
            get;
            private set;
        }
        
        /// <summary>
        /// grid object
        /// </summary>
        public IGridData GridDat {
            get {
                return LsTrk.GridDat;
            }
        }

        /// <summary>
        /// Constructor for an XDG operator
        /// </summary>
        public XdgTimestepping(
            XSpatialOperatorMk2 op,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            DelUpdateLevelset _UpdateLevelset = null,
            LevelSetHandling _LevelSetHandling = LevelSetHandling.None,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig = null,
            AggregationGridData[] _MultigridSequence = null,
            double _AgglomerationThreshold = 0.1,
            LinearSolverConfig LinearSolver = null, NonLinearSolverConfig NonLinearSolver = null,
            LevelSetTracker _optTracker = null,
            IList<DGField> _Parameters = null) //
        {
            this.Scheme = __Scheme;
            this.XdgOperator = op;

            if (_Parameters.IsNullOrEmpty())
                this.Parameters = op.InvokeParameterFactory(Fields);
            else
                this.Parameters = _Parameters;

            LsTrk = _optTracker;
            foreach(var f in Fields.Cat(IterationResiduals).Cat(Parameters)) {
                if(f != null && f is XDGField xf) {
                    if(LsTrk == null) {
                        LsTrk = xf.Basis.Tracker;
                    } else {
                        if(!object.ReferenceEquals(LsTrk, xf.Basis.Tracker))
                            throw new ArgumentException();
                    }
                }
            }
            if(LsTrk == null)
                throw new ArgumentException("unable to get Level Set Tracker reference");

            if(op.AgglomerationThreshold != _AgglomerationThreshold)
                throw new ArgumentException("Mismatch between agglomeration threshold provided ");

            bool UseX = Fields.Any(f => f is XDGField) || IterationResiduals.Any(f => f is XDGField);

            SpeciesId[] spcToCompute = op.Species.Select(spcName => LsTrk.GetSpeciesId(spcName)).ToArray();

            ConstructorCommon(op, UseX,
                Fields, this.Parameters, IterationResiduals, 
                spcToCompute,
                _UpdateLevelset, 
                _LevelSetHandling,
                _MultigridOperatorConfig, 
                _MultigridSequence, 
                _AgglomerationThreshold,
                LinearSolver, NonLinearSolver);

        }       

        private void ConstructorCommon(
            ISpatialOperator op, bool UseX, 
            IEnumerable<DGField> Fields, IEnumerable<DGField> __Parameters, IEnumerable<DGField> IterationResiduals, 
            SpeciesId[] spcToCompute,
            DelUpdateLevelset _UpdateLevelset, LevelSetHandling _LevelSetHandling, 
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig, AggregationGridData[] _MultigridSequence, 
            double _AgglomerationThreshold,
            LinearSolverConfig LinearSolver, NonLinearSolverConfig NonLinearSolver) //
        {
            RungeKuttaScheme rksch;
            int bdfOrder;
            DecodeScheme(this.Scheme, out rksch, out bdfOrder);

            SpatialOperatorType _SpatialOperatorType = op.IsLinear ? SpatialOperatorType.LinearTimeDependent : SpatialOperatorType.Nonlinear;


            int quadOrder = op.QuadOrderFunction(
                Fields.Select(f => f.Basis.Degree).ToArray(),
                Parameters.Select(f => f != null ? f.Basis.Degree : 0).ToArray(),
                IterationResiduals.Select(f => f.Basis.Degree).ToArray());

            // default solvers
            // ===============
            if(LinearSolver == null) {
                LinearSolver = new LinearSolverConfig() {
                    SolverCode = LinearSolverCode.automatic
                };
            }
            if (NonLinearSolver == null) {
                NonLinearSolver = new NonLinearSolverConfig() {
                    SolverCode = NonLinearSolverCode.Newton
                };
            }

            // default Multi-Grid
            // ==================

            if (_MultigridSequence == null) {
                _MultigridSequence = new[] { CoarseningAlgorithms.ZeroAggregation(this.GridDat) };
            }

            // default level-set treatment
            // ===========================

            if (_UpdateLevelset == null) {
                _UpdateLevelset = this.UpdateLevelsetWithNothing;
                if (_LevelSetHandling != LevelSetHandling.None)
                    throw new ArgumentException($"If level-set handling is set to {_LevelSetHandling} (anything but {LevelSetHandling.None}) an updating routine must be specified.");
            }

            // default multigrid operator config
            // =================================
            if(_MultigridOperatorConfig == null) {
                int NoOfVar = Fields.Count();
                _MultigridOperatorConfig = new MultigridOperator.ChangeOfBasisConfig[1][];
                _MultigridOperatorConfig[0] = new MultigridOperator.ChangeOfBasisConfig[NoOfVar];
                for(int iVar = 0; iVar < NoOfVar; iVar++) {
                    _MultigridOperatorConfig[0][iVar] = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { Fields.ElementAt(iVar).Basis.Degree },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { iVar }
                    };
                }

            }

            // finally, create timestepper
            // ===========================

            if(bdfOrder > -1000) {
                m_BDF_Timestepper = new XdgBDFTimestepping(Fields, __Parameters, IterationResiduals,
                    LsTrk, true,
                    this.ComputeOperatorMatrix, op, _UpdateLevelset,
                    bdfOrder,
                    _LevelSetHandling,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    _SpatialOperatorType,
                    _MultigridOperatorConfig, _MultigridSequence,
                    spcToCompute, quadOrder,
                    _AgglomerationThreshold, UseX,
                    NonLinearSolver,
                    LinearSolver);

                m_BDF_Timestepper.Config_AgglomerationThreshold = _AgglomerationThreshold;
            } else {
                m_RK_Timestepper = new XdgRKTimestepping(Fields.ToArray(), __Parameters, IterationResiduals.ToArray(),
                    LsTrk,
                    this.ComputeOperatorMatrix, op, _UpdateLevelset,
                    rksch,
                    _LevelSetHandling,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    _SpatialOperatorType,
                    _MultigridOperatorConfig, _MultigridSequence,
                    spcToCompute, quadOrder,
                    _AgglomerationThreshold, UseX,
                    NonLinearSolver,
                    LinearSolver);

                m_RK_Timestepper.Config_AgglomerationThreshold = _AgglomerationThreshold;
            }
        }

        /// <summary>
        /// translates a time-stepping scheme code
        /// </summary>
        /// <param name="Scheme"></param>
        /// <param name="rksch">if <paramref name="Scheme"/> denotes a Runge-Kutta scheme, well, the Runge-Kutta scheme</param>
        /// <param name="bdfOrder">if <paramref name="Scheme"/> denotes a BDF-scheme, its order</param>
        static public void DecodeScheme(TimeSteppingScheme Scheme, out RungeKuttaScheme rksch, out int bdfOrder) {
            rksch = null;
            bdfOrder = -1000;
            if(Scheme == TimeSteppingScheme.CrankNicolson)
                bdfOrder = -1;
            else if(Scheme == TimeSteppingScheme.ExplicitEuler)
                bdfOrder = 0;
            else if(Scheme == TimeSteppingScheme.ImplicitEuler)
                bdfOrder = 1;
            else if(Scheme.ToString().StartsWith("BDF"))
                bdfOrder = Convert.ToInt32(Scheme.ToString().Substring(3));
            else if(Scheme == TimeSteppingScheme.RK1)
                rksch = RungeKuttaScheme.ExplicitEuler;
            else if(Scheme == TimeSteppingScheme.RK1u1)
                rksch = RungeKuttaScheme.ExplicitEuler2;
            else if(Scheme == TimeSteppingScheme.RK2)
                rksch = RungeKuttaScheme.Heun2;
            else if(Scheme == TimeSteppingScheme.RK3)
                rksch = RungeKuttaScheme.TVD3;
            else if(Scheme == TimeSteppingScheme.RK4)
                rksch = RungeKuttaScheme.RungeKutta1901;
            else if(Scheme == TimeSteppingScheme.RK_ImplicitEuler)
                rksch = RungeKuttaScheme.ImplicitEuler;
            else if(Scheme == TimeSteppingScheme.RK_CrankNic)
                rksch = RungeKuttaScheme.CrankNicolson;
            else if(Scheme == TimeSteppingScheme.RK_IMEX3)
                rksch = RungeKuttaScheme.IMEX3;
            else
                throw new NotImplementedException();
        }

        /// <summary>
        /// Constructor for conventional (single-phase, non-X) DG
        /// </summary>
        public XdgTimestepping(
            SpatialOperator op,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig = null,
            AggregationGridData[] _MultigridSequence = null,
            LinearSolverConfig LinearSolver = null, NonLinearSolverConfig NonLinearSolver = null,
            IList<DGField> _Parameters = null) //
        {
            this.Scheme = __Scheme;
            this.DgOperator = op;

            if (_Parameters.IsNullOrEmpty())
                this.Parameters = op.InvokeParameterFactory(Fields);
            else
                this.Parameters = _Parameters;


            var spc = CreateDummyTracker(Fields.First().GridDat);
                       
            ConstructorCommon(op, false,
                Fields, this.Parameters, IterationResiduals,
                new[] { spc },
                UpdateLevelsetWithNothing,
                LevelSetHandling.None,
                _MultigridOperatorConfig,
                _MultigridSequence,
                0.0,
                LinearSolver, NonLinearSolver);
        }

        /// <summary>
        /// some hack to help with load balancing
        /// </summary>
        internal void RecreateDummyTracker(IGridData newMesh) {
            CreateDummyTracker(newMesh);
        }

        private SpeciesId CreateDummyTracker(IGridData _gDat) {
            var gDat = (GridData) _gDat;
            var DummyLevSet = new LevelSet(new Basis(gDat, 1), "DummyPhi");
            DummyLevSet.AccConstant(-1.0);
            LsTrk = new LevelSetTracker(gDat, XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new[] { "A", "B" }, DummyLevSet);
            LsTrk.UpdateTracker(0.0, __NearRegionWith: 0);

            var spcA = LsTrk.GetSpeciesId("A");

            Debug.Assert(LsTrk.Regions.GetSpeciesMask(spcA).NoOfItemsLocally == gDat.Cells.NoOfLocalUpdatedCells);

            return spcA;
        }

        double UpdateLevelsetWithNothing(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental) {
            this.LsTrk.UpdateTracker(time + dt, incremental: true);
            return 0.0;
        }

        DGField[] JacobiParameterVars = null;

        
        /// <summary>
        /// Operator Evaluation and Linearization, <see cref="DelComputeOperatorMatrix"/>
        /// </summary>
        public void ComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] __CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time, int LsTrkHistoryIndex) {
            // compute operator
            Debug.Assert(OpAffine.L2Norm() == 0.0);

            // all kinds of checks
            if(!this.CurrentState.EqualsPartition(Mapping))
                throw new ApplicationException("something is weired");
            if(OpMtx != null) {
                if(!OpMtx.RowPartitioning.EqualsPartition(Mapping))
                    throw new ArgumentException("Codomain/Matrix Row mapping mismatch.");
                if(!OpMtx.ColPartition.EqualsPartition(Mapping))
                    throw new ArgumentException("Domain/Matrix column mapping mismatch.");
            }
                    
            if(XdgOperator != null) {
                // +++++++++++++++++++++++++++++++++++++++++++++++
                // XDG Branch: still requires length-scale-hack
                // (should be cleaned some-when in the future)
                // +++++++++++++++++++++++++++++++++++++++++++++++

                //if(XdgOperator.AgglomerationThreshold <= 0)
                //    throw new ArgumentException("Mismatch between agglomeration threshold provided ");

                if(OpMtx != null) {
                    // +++++++++++++++++++++++++++++
                    // Solver requires linearization
                    // +++++++++++++++++++++++++++++

                    Debug.Assert(OpMtx.InfNorm() == 0.0);
                    switch(XdgOperator.LinearizationHint) {

                        case LinearizationHint.AdHoc: {
                            this.XdgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

                            var mtxBuilder = XdgOperator.GetMatrixBuilder(LsTrk, Mapping, this.Parameters, Mapping, LsTrkHistoryIndex);
                            mtxBuilder.time = time;
                            mtxBuilder.MPITtransceive = true;
                            foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
                                mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                            }
                            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                            return;
                        }

                        case LinearizationHint.FDJacobi: {
                            var mtxBuilder = XdgOperator.GetFDJacobianBuilder(LsTrk, __CurrentState, this.Parameters, Mapping, LsTrkHistoryIndex);
                            mtxBuilder.time = time;
                            mtxBuilder.MPITtransceive = true;
                            if(mtxBuilder.Eval is XSpatialOperatorMk2.XEvaluatorNonlin evn) { // length-scale hack
                                foreach(var kv in AgglomeratedCellLengthScales) {
                                    evn.CellLengthScales[kv.Key] = kv.Value;
                                }
                            }
                            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                            return;
                        }

                        case LinearizationHint.GetJacobiOperator: {
                            var op = GetJacobiXdgOperator();

                            if(JacobiParameterVars == null)
                                JacobiParameterVars = op.InvokeParameterFactory(this.CurrentState);

                            op.InvokeParameterUpdate(time, __CurrentState, JacobiParameterVars);

                            var mtxBuilder = op.GetMatrixBuilder(LsTrk, Mapping, this.JacobiParameterVars, Mapping, LsTrkHistoryIndex);
                            mtxBuilder.time = time;
                            mtxBuilder.MPITtransceive = true;
                            foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
                                mtxBuilder.CellLengthScales[kv.Key] = kv.Value;
                            }
                            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                            return;
                        }
                    }
                } else {
                    // ++++++++++++++++++++++++
                    // only operator evaluation
                    // ++++++++++++++++++++++++

                   

                    this.XdgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

                    var eval = XdgOperator.GetEvaluatorEx(this.LsTrk, __CurrentState, this.Parameters, Mapping, LsTrkHistoryIndex);
                    eval.time = time;

                    eval.MPITtransceive = true;
                    foreach(var kv in AgglomeratedCellLengthScales) { // length-scale hack
                        eval.CellLengthScales[kv.Key] = kv.Value;
                    }
                    eval.Evaluate(1.0, 0.0, OpAffine);
                }



            } else if(DgOperator != null) {
                // +++++++++++++++++++++++++++++++++++++++++++++++
                // DG Branch
                // +++++++++++++++++++++++++++++++++++++++++++++++

                if(OpMtx != null) {
                    // +++++++++++++++++++++++++++++
                    // Solver requires linearization
                    // +++++++++++++++++++++++++++++

                    Debug.Assert(OpMtx.InfNorm() == 0.0);
                    switch(DgOperator.LinearizationHint) {

                        case LinearizationHint.AdHoc: {
                            this.DgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

                            var mtxBuilder = DgOperator.GetMatrixBuilder(Mapping, this.Parameters, Mapping);
                            mtxBuilder.time = time;
                            mtxBuilder.MPITtransceive = true;
                            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                            return;
                        }

                        case LinearizationHint.FDJacobi: {
                            var mtxBuilder = DgOperator.GetFDJacobianBuilder(__CurrentState, this.Parameters, Mapping);
                            mtxBuilder.time = time;
                            mtxBuilder.MPITtransceive = true;
                            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                            return;
                        }

                        case LinearizationHint.GetJacobiOperator: {
                            var op = GetJacobiDgOperator();

                            if(JacobiParameterVars == null)
                                JacobiParameterVars = op.InvokeParameterFactory(__CurrentState);
                            op.InvokeParameterUpdate(time, __CurrentState, JacobiParameterVars);

                            var mtxBuilder = op.GetMatrixBuilder(Mapping, this.JacobiParameterVars, Mapping);
                            mtxBuilder.time = time;
                            mtxBuilder.MPITtransceive = true;
                            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
                            return;
                        }
                    }
                } else {
                    // ++++++++++++++++++++++++
                    // only operator evaluation
                    // ++++++++++++++++++++++++

                    this.DgOperator.InvokeParameterUpdate(time, __CurrentState, this.Parameters.ToArray());

                    var eval = DgOperator.GetEvaluatorEx(__CurrentState, this.Parameters, Mapping);
                    eval.time = time;
                    eval.MPITtransceive = true;
                    eval.Evaluate(1.0, 0.0, OpAffine);
                }
            } else {
                throw new NotImplementedException();
            }
        }


        /// <summary>
        /// Intended For analysis purposes:
        /// Returns the Jacobi matrix at the latest/current linearization point, **not** including any temporal derivative.
        /// </summary>
        public BlockMsrMatrix GetCurrentSpatialMatrix() {
            var OpMtx = new BlockMsrMatrix(this.IterationResiduals, this.CurrentState);

            ComputeOperatorMatrix(OpMtx, new double[OpMtx.RowPartitioning.LocalLength], this.CurrentState, this.CurrentState.Fields.ToArray(),
                this.TimesteppingBase.GetAgglomeratedLengthScales(), this.TimesteppingBase.GetSimulationTime(), 1);

            return OpMtx;
        }
        
        /// <summary>
        /// Intended For analysis purposes:
        /// Returns the Jacobi matrix at the latest/current linearization point, **including** the temporal derivative
        /// </summary>
        public BlockMsrMatrix GetCurrentJacobiMatrix() {

            this.TimesteppingBase.AssembleMatrixCallback(out var OpMtx, out _, out _, CurrentState.Fields.ToArray(), true, out _);

            return OpMtx;
        }


        /// <summary>
        /// the internal object
        /// </summary>
        public XdgTimesteppingBase TimesteppingBase {
            get {
                Debug.Assert((m_BDF_Timestepper == null) != (m_RK_Timestepper == null));
                if(m_BDF_Timestepper != null)
                    return m_BDF_Timestepper;
                if(m_RK_Timestepper != null)
                    return m_RK_Timestepper;
                throw new ApplicationException("internal error");
            }
        }



        /// <summary>
        /// Returns a collection of local and global condition numbers in order to assess the operators stability
        /// </summary>
        public IDictionary<string, double> OperatorAnalysis(IEnumerable<int[]> VarGroups = null) {
            return TimesteppingBase.OperatorAnalysis(VarGroups);
        }

        public XdgBDFTimestepping m_BDF_Timestepper;

        public XdgRKTimestepping m_RK_Timestepper;
      

        /// <summary>
        /// Add logging of iteration residuals.
        /// </summary>
        public void RegisterResidualLogger(ResidualLogger _ResLogger) {
            TimesteppingBase.m_ResidualNames = this.Operator.CodomainVar.ToArray();
            TimesteppingBase.m_ResLogger = _ResLogger;
        }


        /// <summary>
        /// driver for solver calls
        /// </summary>
        /// <returns>
        /// - true: solver algorithm successfully converged
        /// - false: something went wrong
        /// </returns>
        public bool Solve(double phystime, double dt, bool SkipSolveAndEvaluateResidual = false) {
            bool success = false;
            if((m_BDF_Timestepper == null) == (m_RK_Timestepper == null))
                throw new ApplicationException();

            double[] AvailTimesBefore;
            if(TimesteppingBase.Config_LevelSetHandling != LevelSetHandling.None) {
                AvailTimesBefore = LsTrk.TimeLevelsInStack;
                Assert.IsTrue((AvailTimesBefore[0] - phystime).Abs() < dt*1e-7, "Error in Level-Set tracker time");
            }

            if(m_BDF_Timestepper != null) {
                success = m_BDF_Timestepper.Solve(phystime, dt, SkipSolveAndEvaluateResidual);
            } else {
                if (SkipSolveAndEvaluateResidual == true)
                    throw new NotSupportedException("SkipSolveAndEvaluateResidual == true is not supported for Runge-Kutta");

                success = m_RK_Timestepper.Solve(phystime, dt);
            }

            double[] AvailTimesAfter;
            if(TimesteppingBase.Config_LevelSetHandling != LevelSetHandling.None) {
                AvailTimesAfter = LsTrk.RegionsHistory.AvailabelIndices.Select((int iHist) => LsTrk.RegionsHistory[iHist].Time).ToArray();
                Assert.IsTrue((AvailTimesAfter[0] - (phystime + dt)).Abs() < dt*1e-7, "Error in Level-Set tracker time");
            }

            JacobiParameterVars = null;
            return success;
        }

        /// <summary>
        /// Step 2 of 2 for dynamic load balancing: restore this objects 
        /// status after the grid has been re-distributed.
        /// </summary>
        public void DataRestoreAfterBalancing(GridUpdateDataVaultBase L,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            AggregationGridData[] _MultigridSequence,
            ISpatialOperator abstractOperator) //
        {
            var gDat = Fields.First().GridDat;
            if(!object.ReferenceEquals(LsTrk.GridDat, gDat))
                throw new ApplicationException();
            if(LsTrk != null)
                this.LsTrk = LsTrk;

            // Reset dependent spatial operator, otherwise AMR can lead e.g. to invalid LsTrks in some equation components
            if (this.m_JacobiXdgOperator != null)
                this.m_JacobiXdgOperator = null;

            Parameters = this.Operator.InvokeParameterFactory(Fields);
            
            if(m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DataRestoreAfterBalancing(L, Fields, Parameters, IterationResiduals, LsTrk, _MultigridSequence, abstractOperator);
            } else if(m_RK_Timestepper != null) {
                throw new NotImplementedException("Load balancing and adaptive mesh refinement are not supported for Runge-Kutta XDG timestepping.");
            } else {
                throw new NotImplementedException();
            }

            if (abstractOperator.GetType() == typeof(XSpatialOperatorMk2))
                XdgOperator = (XSpatialOperatorMk2)abstractOperator;
            if (abstractOperator.GetType() == typeof(SpatialOperator))
                DgOperator = (SpatialOperator)abstractOperator;

            if (!object.ReferenceEquals(this.Operator, m_BDF_Timestepper.AbstractOperator))
                throw new ApplicationException();

        }

        /// <summary>
        /// Step 1 of 2 for dynamic load balancing: creating a backup of this objects 
        /// status in the load-balancing thing <paramref name="L"/>
        /// </summary>
        public void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {
            
            if(m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DataBackupBeforeBalancing(L);
            } else if(m_RK_Timestepper != null) {
                throw new NotImplementedException("Load balancing and adaptive mesh refinement are not supported for Runge-Kutta XDG timestepping.");
            } else {
                throw new NotImplementedException();
            }
                        
            this.LsTrk = null;
            this.Parameters = null;
        }

        
        // <summary>
        /// Number of timesteps required for restart, e.g. 1 for Runge-Kutta and implicit/explicit Euler, 2 for BDF2, etc.
        /// </summary>
        public int BurstSave {
            get {
                if(m_RK_Timestepper != null) {
                    Debug.Assert(m_BDF_Timestepper == null);
                    return 1;
                } else if(m_BDF_Timestepper != null) {
                    return m_BDF_Timestepper.GetNumberOfStages;
                } else {
                    throw new NotImplementedException();
                }

            }
        }
    }
}
