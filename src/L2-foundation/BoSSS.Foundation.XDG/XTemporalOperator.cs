using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {



    /// <summary>
    /// A constant, diagonal operator
    /// </summary>
    public class ConstantXTemporalOperator : ITemporalOperator {

        static ValueTuple<string, double[]>[] Diag(XDifferentialOperatorMk2 __owner, double diagonalValue) {
            double[] diag = new double[__owner.DomainVar.Count];
            diag.SetAll(diagonalValue);
            return __owner.Species.Select(spcName => (spcName, diag.CloneAs())).ToArray();
        }


        /// <summary>
        /// Ctor, the same factor in the temporal derivative for all 
        /// </summary>
        /// <param name="__owner"></param>
        /// <param name="diagonalValue"></param>
        public ConstantXTemporalOperator(XDifferentialOperatorMk2 __owner, double diagonalValue = 1.0)
            : this(__owner, Diag(__owner, diagonalValue)) {

            owner = __owner;
            if (owner.DomainVar.Count != owner.CodomainVar.Count) {
                throw new NotSupportedException("Expecting a square operator.");
            }
        }


        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="__owner"></param>
        /// <param name="diagonal">
        /// For each species, the temporal operator/mass matrix diagonal (per variable), 
        /// i.e. if <paramref name="__owner"/> has 4 domain and 4 codomain variables,
        /// this must have 4 entries per species too, providing a single factor for the temporal derivative of each equation.
        /// </param>
        public ConstantXTemporalOperator(XDifferentialOperatorMk2 __owner, params ValueTuple<string, double[]>[] diagonal) {
            owner = __owner;
            if (owner.DomainVar.Count != owner.CodomainVar.Count) {
                throw new NotSupportedException("Expecting a square operator.");
            }

            int NoOfVar = owner.DomainVar.Count;

            m_DiagonalScale = new Dictionary<string, double[]>();
            foreach(string spc in __owner.Species) {
                m_DiagonalScale.Add(spc, new double[NoOfVar]);
            }

            foreach(var tt in diagonal) {
                string spc = tt.Item1;
                if (!m_DiagonalScale.ContainsKey(spc))
                    throw new ArgumentException($"Got species {spc} which is not contained in the spatial operator.");

                double[] vals = tt.Item2.CloneAs();
                if (vals.Length != NoOfVar)
                    throw new ArgumentException($"wrong number of entries for species {spc}");
                m_DiagonalScale[spc].SetV(vals, 1.0);
            }
        }

        XDifferentialOperatorMk2 owner;

        Dictionary<string, double[]> m_DiagonalScale;

        /// <summary>
        /// For each species, the temporal operator/mass matrix diagonal (per variable), 
        /// i.e. if the spatial operator has 4 domain and 4 codomain variables,
        /// this must have 4 entries per species too, providing a single factor for the temporal derivative of each equation.
        /// </summary>
        public IReadOnlyDictionary<string, double[]> DiagonalScale {
            get {
                return m_DiagonalScale;
            }

        }

        bool m_IsCommited;

        /// <summary>
        /// locks the configuration of the operator
        /// </summary>
        public void Commit() {
            if (m_IsCommited)
                throw new ApplicationException("'Commit' has already been called - it can be called only once in the lifetime of this object.");
            m_IsCommited = true;
        }

        /// <summary>
        /// Returns an matrix builder which always accumulates the same diagonal matrix.
        /// </summary>
        public IEvaluatorLinear GetMassMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
                     

            return new InternalBla(this, DomainVarMap, ParameterMap, CodomainVarMap);
        }


        /// <summary>
        /// Dirty hack to support e.g. the IBM solver (state sept2020) which uses only DG 
        /// fields but employs an <see cref="XDifferentialOperatorMk2"/>;
        /// may be deleted, eventually.
        /// </summary>
        public void SetTrackerHack(LevelSetTracker lstrk) {
            m_lstrk = lstrk;
        }

        LevelSetTracker m_lstrk;

        class InternalBla : IEvaluatorLinear {

            ConstantXTemporalOperator m_Owner;


            internal InternalBla(ConstantXTemporalOperator __Owner, UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
                this.DomainMapping = DomainVarMap;
                this.Parameters = ParameterMap;
                this.CodomainMapping = CodomainVarMap;
                m_Owner = __Owner;
                
                if (DomainMapping.NoOfVariables != Owner.DomainVar.Count)
                    throw new ArgumentException("Mismatch in number of domain variables.");
                if (CodomainMapping.NoOfVariables != Owner.CodomainVar.Count)
                    throw new ArgumentException("Mismatch in number of codomain variables.");

                if (ParameterMap.Count != Owner.ParameterVar.Count)
                    throw new ArgumentException("Mismatch in number of parameter variables.");

                if (!DomainVarMap.EqualsPartition(CodomainVarMap))
                    throw new ArgumentException("Only supported for square matrices - domain and codomain map are not compatible (e.g. different number of rows and columns.");

            }

            public IGridData GridData => DomainMapping.GridDat;

            public UnsetteledCoordinateMapping CodomainMapping {
                get;
                private set;
            }

            public UnsetteledCoordinateMapping DomainMapping { get; private set; }

            public IList<DGField> Parameters  { get; private set; }

            public SubGridBoundaryModes SubGridBoundaryTreatment {  get; private set; }
        

            /// <summary>
            /// No effect, operator is constant in time
            /// </summary>
            public double time {
                get;
                set;
            }
            
            /// <summary>
            /// ignored - no effect
            /// </summary>
            public bool MPITtransceive {
                get; 
                set; 
            }

            /// <summary>
            /// 
            /// </summary>
            public IDifferentialOperator Owner {
                get {
                    return m_Owner.owner;
                }
            }

            public void ActivateSubgridBoundary(CellMask mask, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                this.SubGridBoundaryTreatment = subGridBoundaryTreatment;
                if (!object.ReferenceEquals(mask.GridData, this.GridData))
                    throw new ArgumentException("grid mismatch");
                subMask = mask;
            }

            /// <summary>
            /// No effect;
            /// </summary>
            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
            }

            CellMask subMask;

            /// <summary>
            /// Computation of mass matrix, makes use of <see cref="MassMatrixFactory"/>
            /// </summary>
            public void ComputeMatrix<M, V>(M MassMatrix, V AffineOffset, double oodt = 1.0)
                where M : IMutableMatrixEx
                where V : IList<double> {

                if (!MassMatrix.RowPartitioning.EqualsPartition(this.CodomainMapping))
                    throw new ArgumentException("Mismatch between matrix columns and domain variable mapping.");
                if (!MassMatrix.ColPartition.EqualsPartition(this.DomainMapping))
                    throw new ArgumentException("Mismatch between matrix columns and domain variable mapping.");
                if (DomainMapping.NoOfVariables != CodomainMapping.NoOfVariables)
                    throw new NotSupportedException($"Mismatch between number of variables in domain ({DomainMapping.NoOfVariables}) and codomain ({CodomainMapping.NoOfVariables}).");

                int QuadratureOrder = m_Owner.owner.GetOrderFromQuadOrderFunction(DomainMapping.BasisS, XDifferentialOperatorMk2.GetBasisS(Parameters), CodomainMapping.BasisS);
                LevelSetTracker _LsTrk;
                if(m_Owner.m_lstrk != null)
                    _LsTrk = m_Owner.m_lstrk;
                else
                    _LsTrk = XDifferentialOperatorMk2.GetTracker(DomainMapping.BasisS, XDifferentialOperatorMk2.GetBasisS(Parameters), CodomainMapping.BasisS);
                SpeciesId[] _SpeciesToCompute = m_Owner.owner.Species.Select(spcName => _LsTrk.GetSpeciesId(spcName)).ToArray();


                Dictionary<SpeciesId, IEnumerable<double>> MassScale = new Dictionary<SpeciesId, IEnumerable<double>>();
                foreach(var spc in _SpeciesToCompute) {
                    double[] diag = m_Owner.m_DiagonalScale[_LsTrk.GetSpeciesName(spc)].CloneAs();
                    diag.ScaleV(oodt);
                    MassScale.Add(spc, diag);
                }

                //double MaxTime = _LsTrk.RegionsHistory.AvailabelIndices.Max((int iHist) => _LsTrk.RegionsHistory[iHist].Time.Abs());
                double[] AvailableTimes = _LsTrk.TimeLevelsInStack;
                int ii = AvailableTimes.IndexOfMin((double t) => Math.Abs(t - this.time));
                int BestTimeIdx = _LsTrk.RegionsHistory.AvailableIndices[ii];

                if(Math.Abs(Math.Abs(time - _LsTrk.RegionsHistory[BestTimeIdx].Time)) >= time * 1e-10 + 1e-10)
                    Console.WriteLine("unknown time level");


                //BestTimeIdx = 1;
                MassMatrixFactory MassFact = _LsTrk.GetXDGSpaceMetrics(_SpeciesToCompute, QuadratureOrder, HistoryIndex:BestTimeIdx).MassMatrixFactory;
                MassFact.AccMassMatrix(MassMatrix, DomainMapping, MassScale);

            }
        }

    }


    /// <summary>
    /// Temporal operator which is fully customizable, i.e. can be configured through <see cref="EquationComponents"/>
    /// in analog fashion to the spatial operator.
    /// </summary>
    public class DependentXTemporalOperator : ITemporalOperator {

        XDifferentialOperatorMk2 owner;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="__owner"></param>
        public DependentXTemporalOperator(XDifferentialOperatorMk2 __owner) {
            owner = __owner;
            InternalRepresentation = new XDifferentialOperatorMk2(__owner.DomainVar, __owner.ParameterVar, __owner.CodomainVar, __owner.QuadOrderFunction, __owner.Species.ToArray());
            InternalRepresentation.AgglomerationThreshold = __owner.AgglomerationThreshold;
            InternalRepresentation.LinearizationHint = LinearizationHint.GetJacobiOperator;
            InternalRepresentation.m_UserDefinedValues = __owner.m_UserDefinedValues;
        }

        /// <summary>
        /// Components of the Mass matrix
        /// </summary>
        public IEquationComponents EquationComponents { 
            get {
                return InternalRepresentation.EquationComponents;
            }
        
        }

        XDifferentialOperatorMk2 InternalRepresentation;

        /// <summary>
        /// locks the configuration of the operator
        /// </summary>
        public void Commit() {
            InternalRepresentation.OperatorCoefficientsProvider = owner.OperatorCoefficientsProvider;
            InternalRepresentation.LinearizationHint = LinearizationHint.GetJacobiOperator;

            InternalRepresentation.ParameterFactories.Clear();
            InternalRepresentation.ParameterFactories.AddRange(owner.ParameterFactories);

            InternalRepresentation.ParameterUpdates.Clear();
            InternalRepresentation.ParameterUpdates.AddRange(owner.ParameterUpdates);

            InternalRepresentation.EdgeQuadraturSchemeProvider = owner.EdgeQuadraturSchemeProvider;
            InternalRepresentation.VolumeQuadraturSchemeProvider = owner.VolumeQuadraturSchemeProvider;
            InternalRepresentation.SurfaceElement_EdgeQuadraturSchemeProvider = owner.SurfaceElement_EdgeQuadraturSchemeProvider;
            InternalRepresentation.SurfaceElement_VolumeQuadraturSchemeProvider = owner.SurfaceElement_VolumeQuadraturSchemeProvider;
            InternalRepresentation.GhostEdgeQuadraturSchemeProvider = owner.GhostEdgeQuadraturSchemeProvider;
            InternalRepresentation.ContactLine_VolumeQuadratureSchemeProvider = owner.ContactLine_VolumeQuadratureSchemeProvider;

            InternalRepresentation.m_UserDefinedValues = owner.m_UserDefinedValues;
            InternalRepresentation.OperatorCoefficientsProvider = owner.OperatorCoefficientsProvider;
            
            InternalRepresentation.Commit();
        }

        /// <summary>
        /// as defined by interface
        /// </summary>
        public IEvaluatorLinear GetMassMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
            //InternalRepresentation.EdgeQuadraturSchemeProvider = owner.EdgeQuadraturSchemeProvider;
            //InternalRepresentation.VolumeQuadraturSchemeProvider = owner.VolumeQuadraturSchemeProvider;
            //InternalRepresentation.SurfaceElement_EdgeQuadraturSchemeProvider = owner.SurfaceElement_EdgeQuadraturSchemeProvider;
            //InternalRepresentation.SurfaceElement_VolumeQuadraturSchemeProvider = owner.SurfaceElement_VolumeQuadraturSchemeProvider;
            //InternalRepresentation.GhostEdgeQuadraturSchemeProvider = owner.GhostEdgeQuadraturSchemeProvider;
            //InternalRepresentation.ContactLine_VolumeQuadratureSchemeProvider = owner.ContactLine_VolumeQuadratureSchemeProvider;

            InternalRepresentation.m_UserDefinedValues = owner.m_UserDefinedValues;
            InternalRepresentation.OperatorCoefficientsProvider = owner.OperatorCoefficientsProvider;

            InternalRepresentation.CurrentHomotopyValue = owner.CurrentHomotopyValue;        
            return InternalRepresentation.GetFDJacobianBuilder((CoordinateMapping)DomainVarMap, ParameterMap, (CoordinateMapping)CodomainVarMap);
        }
    }

}
