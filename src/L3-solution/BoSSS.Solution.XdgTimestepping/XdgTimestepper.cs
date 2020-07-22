using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.Timestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {


    public enum TimeSteppingScheme {
        ExplicitEuler = 0,

        /// <summary>
        /// equivalent of BDF1, using the BDF-implementation
        /// </summary>
        ImplicitEuler = 1,

        CrankNicolson = 2000,

        BDF2 = 2,

        BDF3 = 3,

        BDF4 = 4,

        BDF5 = 5,

        BDF6 = 6,

        RK4 = 104,

        RK3 = 103,

        RK2 = 102,

        RK1 = 101,

        RK1u1 = 110,

        /// <summary>
        /// Implicit Euler using the implicit Runge-Kutta implementation
        /// </summary>
        RK_ImplicitEuler = 201,

        /// <summary>
        /// Crank Nicolson using the implicit Runge-Kutta implementation
        /// </summary>
        RK_CrankNic = 202,

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
        
        public SpatialOperator DgOperator {
            get;
            private set;
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

        public IList<DGField> Parameters {
            get;
            private set;
        }


        public LevelSetTracker LsTrk {
            get;
            private set;
        }
             


        public XdgTimestepping(
            XSpatialOperatorMk2 op,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> __Parameters,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            //bool DelayInit,
            //DelComputeMassMatrix _ComputeMassMatrix,
            DelUpdateLevelset _UpdateLevelset,
            LevelSetHandling _LevelSetHandling,
            //MassMatrixShapeandDependence _MassMatrixShapeandDependence,
            IDictionary<SpeciesId, IEnumerable<double>> _MassScale,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGridData[] _MultigridSequence,
            double _AgglomerationThreshold,
            AppControlSolver control) {
            this.Scheme = __Scheme;
            this.XdgOperator = op;

            if(__Parameters == null)
                this.Parameters = new DGField[0];
            else
                this.Parameters = __Parameters.ToArray();


            foreach(var f in Fields.Cat(IterationResiduals).Cat(Parameters)) {
                if(f is XDGField xf) {
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

            bool UseX = Fields.Any(f => f is XDGField) || IterationResiduals.Any(f => f is XDGField);

            ConstructorCommon(op, UseX,
                Fields, IterationResiduals, 
                myDelComputeXOperatorMatrix, 
                _UpdateLevelset, 
                _LevelSetHandling, 
                _MassScale, 
                _MultigridOperatorConfig, 
                _MultigridSequence, 
                _AgglomerationThreshold, 
                control);

        }

        private void ConstructorCommon(ISpatialOperator op, bool UseX, IEnumerable<DGField> Fields, IEnumerable<DGField> IterationResiduals, DelComputeOperatorMatrix __delComputeOperatorMatrix, DelUpdateLevelset _UpdateLevelset, LevelSetHandling _LevelSetHandling, IDictionary<SpeciesId, IEnumerable<double>> _MassScale, MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig, AggregationGridData[] _MultigridSequence, double _AgglomerationThreshold, AppControlSolver control) {
            RungeKuttaScheme rksch = null;
            int bdfOrder = -1000;
            if(this.Scheme == TimeSteppingScheme.CrankNicolson)
                bdfOrder = -1;
            else if(this.Scheme == TimeSteppingScheme.ExplicitEuler)
                bdfOrder = 0;
            else if(this.Scheme == TimeSteppingScheme.ImplicitEuler)
                bdfOrder = 1;
            else if(this.Scheme.ToString().StartsWith("BDF"))
                bdfOrder = Convert.ToInt32(this.Scheme.ToString().Substring(3));
            else if(this.Scheme == TimeSteppingScheme.RK1)
                rksch = RungeKuttaScheme.ExplicitEuler;
            else if(this.Scheme == TimeSteppingScheme.RK1u1)
                rksch = RungeKuttaScheme.ExplicitEuler2;
            else if(this.Scheme == TimeSteppingScheme.RK2)
                rksch = RungeKuttaScheme.Heun2;
            else if(this.Scheme == TimeSteppingScheme.RK3)
                rksch = RungeKuttaScheme.TVD3;
            else if(this.Scheme == TimeSteppingScheme.RK4)
                rksch = RungeKuttaScheme.RungeKutta1901;
            else if(this.Scheme == TimeSteppingScheme.RK_ImplicitEuler)
                rksch = RungeKuttaScheme.ImplicitEuler;
            else if(this.Scheme == TimeSteppingScheme.RK_CrankNic)
                rksch = RungeKuttaScheme.CrankNicolson;
            else if(this.Scheme == TimeSteppingScheme.RK_IMEX3)
                rksch = RungeKuttaScheme.IMEX3;
            else
                throw new NotImplementedException();

            SpatialOperatorType _SpatialOperatorType = SpatialOperatorType.Nonlinear;


            int quadOrder = op.QuadOrderFunction(
                Fields.Select(f => f.Basis.Degree).ToArray(),
                Parameters.Select(f => f != null ? f.Basis.Degree : 0).ToArray(),
                IterationResiduals.Select(f => f.Basis.Degree).ToArray());

       

            if(bdfOrder > -1000) {
                m_BDF_Timestepper = new XdgBDFTimestepping(Fields, IterationResiduals,
                    LsTrk, true,
                    __delComputeOperatorMatrix, null, _UpdateLevelset,
                    bdfOrder,
                    _LevelSetHandling,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    _SpatialOperatorType,
                    _MassScale,
                    _MultigridOperatorConfig, _MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(), quadOrder,
                    control.AgglomerationThreshold, UseX,
                    control.NonLinearSolver,
                    control.LinearSolver);

                m_BDF_Timestepper.Config_AgglomerationThreshold = _AgglomerationThreshold;
            } else {
                m_RK_Timestepper = new XdgRKTimestepping(Fields.ToArray(), IterationResiduals.ToArray(),
                    LsTrk,
                    __delComputeOperatorMatrix, _UpdateLevelset,
                    rksch,
                    _LevelSetHandling,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    _SpatialOperatorType,
                    _MassScale,
                    _MultigridOperatorConfig, _MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(), quadOrder,
                    control.AgglomerationThreshold, UseX,
                    control.NonLinearSolver,
                    control.LinearSolver);

                m_RK_Timestepper.Config_AgglomerationThreshold = _AgglomerationThreshold;
            }
        }

        /// <summary>
        /// Constructor for conventional (single-phase, non-X) DG
        /// </summary>
        public XdgTimestepping(
            SpatialOperator op,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> __Parameters,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            IEnumerable<double> _MassScale,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGridData[] _MultigridSequence,
            AppControlSolver control) {
            this.Scheme = __Scheme;
            this.DgOperator = op;

            if(__Parameters == null)
                this.Parameters = new DGField[0];
            else
                this.Parameters = __Parameters.ToArray();

            CreateDummyTracker(Fields);

            var __MassScale = new Dictionary<SpeciesId, IEnumerable<double>>();
            __MassScale.Add(LsTrk.GetSpeciesId("A"), _MassScale);
            __MassScale.Add(LsTrk.GetSpeciesId("B"), _MassScale);


            ConstructorCommon(op, false,
                Fields, IterationResiduals,
                myDelComputeDgOperatorMatrix,
                UpdateLevelsetWithNothing,
                LevelSetHandling.None,
                __MassScale,
                _MultigridOperatorConfig,
                _MultigridSequence,
                0.0,
                control);
        }

        private void CreateDummyTracker(IEnumerable<DGField> Fields) {
            var gDat = Fields.First().GridDat as GridData;
            var DummyLevSet = new LevelSet(new Basis(gDat, 1), "DummyPhi");
            DummyLevSet.AccConstant(-1.0);
            LsTrk = new LevelSetTracker(gDat, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new[] { "A", "B" }, DummyLevSet);
        }

        double UpdateLevelsetWithNothing(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental) {
            this.LsTrk.UpdateTracker(incremental: true);
            return 0.0;
        }


        void myDelComputeXOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time) {
            // compute operator
            Debug.Assert(OpMtx.InfNorm() == 0.0);
            Debug.Assert(OpAffine.L2Norm() == 0.0);

            
            XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = XdgOperator.GetMatrixBuilder(this.LsTrk, Mapping, this.Parameters, Mapping, AgglomeratedCellLengthScales.Keys.ToArray());
            mtxBuilder.time = time;
            mtxBuilder.MPITtransceive = true;
            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
        }

        void myDelComputeDgOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time) {
            // compute operator
            Debug.Assert(OpMtx.InfNorm() == 0.0);
            Debug.Assert(OpAffine.L2Norm() == 0.0);


            var mtxBuilder = DgOperator.GetMatrixBuilder(Mapping, this.Parameters, Mapping);
            mtxBuilder.time = time;
            mtxBuilder.MPITtransceive = true;
            mtxBuilder.ComputeMatrix(OpMtx, OpAffine);
        }



        XdgTimesteppingBase TimesteppingBase {
            get {
                Debug.Assert((m_BDF_Timestepper == null) != (m_RK_Timestepper == null));
                if(m_BDF_Timestepper != null)
                    return m_BDF_Timestepper;
                if(m_RK_Timestepper != null)
                    return m_RK_Timestepper;
                throw new ApplicationException("internal error");
            }
        }


        public XdgBDFTimestepping m_BDF_Timestepper;



        public XdgRKTimestepping m_RK_Timestepper;


        public void Solve(double phystime, double dt) {
            if((m_BDF_Timestepper == null) == (m_RK_Timestepper == null))
                throw new ApplicationException();

            if(m_BDF_Timestepper != null) {
                m_BDF_Timestepper.Solve(phystime, dt);
            } else {
                m_RK_Timestepper.Solve(phystime, dt);
            }

        }

        /// <summary>
        /// Step 2 of 2 for dynamic load balancing: restore this objects 
        /// status after the grid has been re-distributed.
        /// </summary>
        public void DataRestoreAfterBalancing(GridUpdateDataVaultBase L,
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> Params,
            IEnumerable<DGField> IterationResiduals,
            LevelSetTracker LsTrk,
            AggregationGridData[] _MultigridSequence) //
        {
            if(LsTrk != null) {
                this.LsTrk = LsTrk;
            } else {
                CreateDummyTracker(Fields);
            }

            
            if(Params != null) {
                Parameters = Params.ToArray();
            } else {
                Parameters = new DGField[0];
            }
            
            if(m_BDF_Timestepper != null) {
                m_BDF_Timestepper.DataRestoreAfterBalancing(L, Fields, IterationResiduals, LsTrk, _MultigridSequence);
            } else if(m_RK_Timestepper != null) {
                throw new NotImplementedException("Load balancing and adaptive mesh refinement are not supported for Runge-Kutta XDG timestepping.");
            } else {
                throw new NotImplementedException();
            }


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
