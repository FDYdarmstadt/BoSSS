using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.LinSolvers;
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

        /// <summary>
        /// Ctor, the same factor in the temporal derivative for all 
        /// </summary>
        /// <param name="__owner"></param>
        /// <param name="diagonalValue"></param>
        public ConstantXTemporalOperator(XSpatialOperatorMk2 __owner, double diagonalValue = 1.0) 
            : this(__owner, __owner.DomainVar.Select(dv => diagonalValue).ToArray()) { 
            
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
        /// Matrix diagonal per variable, i.e. if <paramref name="__owner"/> has 4 domain and 4 codomain variables,
        /// this must have 4 entries too, providing a single factor for the temporal derivative of each equation.
        /// </param>
        public ConstantXTemporalOperator(XSpatialOperatorMk2 __owner, double[] diagonal) {
            owner = __owner;
            if(owner.DomainVar.Count != owner.CodomainVar.Count) {
                throw new NotSupportedException("Expecting a square operator.");
            }

            m_diagonal = diagonal.CloneAs();
        }

        XSpatialOperatorMk2 owner;

        double[] m_diagonal;

        /// <summary>
        /// locks the configuration of the operator
        /// </summary>
        public void Commit() {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Returns an matrix builder which always accumulates the same diagonal matrix.
        /// </summary>
        public IEvaluatorLinear GetMassMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
                     

            return new InternalBla(this, DomainVarMap, ParameterMap, CodomainVarMap);
        }


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
            public ISpatialOperator Owner {
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

            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset, double oodt = 1.0)
                where M : IMutableMatrixEx
                where V : IList<double> {

                if (!Matrix.RowPartitioning.EqualsPartition(this.CodomainMapping))
                    throw new ArgumentException("Mismatch between matrix columns and domain variable mapping.");
                if (!Matrix.ColPartition.EqualsPartition(this.DomainMapping))
                    throw new ArgumentException("Mismatch between matrix columns and domain variable mapping.");
                if (DomainMapping.NoOfVariables != CodomainMapping.NoOfVariables)
                    throw new NotSupportedException($"Mismatch between number of variables in domain ({DomainMapping.NoOfVariables}) and codomain ({CodomainMapping.NoOfVariables}).");
                

                double[] VarDiag = m_Owner.m_diagonal;
                int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int NoOfVar = DomainMapping.NoOfVariables;
                if (NoOfVar != VarDiag.Length)
                    throw new NotSupportedException();

                Basis[] domBs = DomainMapping.BasisS.ToArray();
                Basis[] codBs = CodomainMapping.BasisS.ToArray();

               
                void SetCellDiag(int j) {
                    for(int iVar = 0; iVar < NoOfVar; iVar++) {
                        int Ncol = domBs[iVar].GetLength(j);
                        int Nrow = codBs[iVar].GetLength(j);
                        if(Nrow != Ncol) {
                            throw new NotSupportedException($"Diagonal block of {iVar}-th variable is not square: domain basis dimension is {Ncol}, codomain basis dimension is {Nrow}.");
                        }

                        for(int n = 0; n < Nrow; n++) {
                            int idx = DomainMapping.GlobalUniqueCoordinateIndex(iVar, j, n);
                            Matrix[idx, idx] += oodt * VarDiag[iVar];
                        }
                    }
                }


                if (subMask == null) {
                    foreach(int j in subMask.ItemEnum) {
                        SetCellDiag(j);
                    }
                } else {
                    for(int j = 0; j < J; j++) {
                        SetCellDiag(j);
                    }
                }
                
            }
        }

    }


    /// <summary>
    /// Temporal operator which is fully customizable, i.e. can be configured through <see cref="EquationComponents"/>
    /// in analog fashion to the spatial operator.
    /// </summary>
    public class DependentXTemporalOperator : ITemporalOperator {

        XSpatialOperatorMk2 owner;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="__owner"></param>
        public DependentXTemporalOperator(XSpatialOperatorMk2 __owner) {
            owner = __owner;
            InternalRepresentation = new XSpatialOperatorMk2(__owner.DomainVar, __owner.ParameterVar, __owner.CodomainVar, __owner.QuadOrderFunction, __owner.Species.ToArray());
            InternalRepresentation.LinearizationHint = LinearizationHint.AdHoc;
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

        XSpatialOperatorMk2 InternalRepresentation;

        /// <summary>
        /// locks the configuration of the operator
        /// </summary>
        public void Commit() {
            InternalRepresentation.OperatorCoefficientsProvider = owner.OperatorCoefficientsProvider;
            InternalRepresentation.LinearizationHint = LinearizationHint.AdHoc;
            InternalRepresentation.ParameterUpdate = owner.ParameterUpdate;
            
            InternalRepresentation.EdgeQuadraturSchemeProvider = owner.EdgeQuadraturSchemeProvider;
            InternalRepresentation.VolumeQuadraturSchemeProvider = owner.VolumeQuadraturSchemeProvider;
            InternalRepresentation.SurfaceElement_EdgeQuadraturSchemeProvider = owner.SurfaceElement_EdgeQuadraturSchemeProvider;
            InternalRepresentation.SurfaceElement_VolumeQuadraturSchemeProvider = owner.SurfaceElement_VolumeQuadraturSchemeProvider;
            InternalRepresentation.GhostEdgeQuadraturSchemeProvider = owner.GhostEdgeQuadraturSchemeProvider;

            InternalRepresentation.m_UserDefinedValues = owner.m_UserDefinedValues;
            InternalRepresentation.Commit();
        }

        /// <summary>
        /// as defined by interface
        /// </summary>
        public IEvaluatorLinear GetMassMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
            InternalRepresentation.EdgeQuadraturSchemeProvider = owner.EdgeQuadraturSchemeProvider;
            InternalRepresentation.VolumeQuadraturSchemeProvider = owner.VolumeQuadraturSchemeProvider;
            InternalRepresentation.SurfaceElement_EdgeQuadraturSchemeProvider = owner.SurfaceElement_EdgeQuadraturSchemeProvider;
            InternalRepresentation.SurfaceElement_VolumeQuadraturSchemeProvider = owner.SurfaceElement_VolumeQuadraturSchemeProvider;
            InternalRepresentation.GhostEdgeQuadraturSchemeProvider = owner.GhostEdgeQuadraturSchemeProvider;

            InternalRepresentation.m_UserDefinedValues = owner.m_UserDefinedValues;

            return InternalRepresentation.GetMatrixBuilder(DomainVarMap, ParameterMap, CodomainVarMap);
        }
    }

}
