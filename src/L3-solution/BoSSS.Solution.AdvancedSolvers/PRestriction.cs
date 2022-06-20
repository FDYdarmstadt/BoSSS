using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {



    /// <summary>
    /// Recursive p-Multigrid on a single mesh level;
    /// this object is not a solver on its own, but it performs the restriction of the operator matrix
    /// to a certain **subset of cells and modes**,
    /// and applies the <see cref="LowerPSolver"/> onto it.
    /// </summary>
    public class PRestriction : ISubsystemSolver {

        /// <summary>
        /// Note: at this point, we do not support everything which a <see cref="SubBlockSelector"/> can support.
        /// This is, because not all selections, which <see cref="SubBlockSelector"/> can do are supported by <see cref="ICoordinateMapping"/> at this point.
        /// </summary>
        class MgOperatorRestriction : IOperatorMappingPair {


            public MgOperatorRestriction(IOperatorMappingPair Unrestricted, int[] RestrictedDeg, Func<int[]> __GetCellRestriction, Func<BlockMsrMatrix> __GetExtRows, bool __MPIself) {
                
                var subBlock = new SubBlockSelector(Unrestricted.DgMapping);

                if(__GetCellRestriction != null) {
                    jLoc2Glob = __GetCellRestriction();
                    subBlock.CellSelector(jLoc2Glob);
                }
                subBlock.SetModeSelector((int jCell, int iVar, int iSpc, int Deg) => Deg <= RestrictedDeg[iVar]);

                m_MPIself = __MPIself;
                m_RestrictedDeg = RestrictedDeg;
                m_Unrestricted = Unrestricted;
                m_Mask = new BlockMask(subBlock, __GetExtRows?.Invoke());

                var idx = m_Mask.GlobalIndices_Internal;
                for(int i = 1; i < idx.Count; i++) {
                    if (idx[i - 1] >= idx[i])
                        throw new Exception("ghdshj: " + Unrestricted.GetType());
                }
            }

            int[] jLoc2Glob = null;

            IOperatorMappingPair m_Unrestricted;
            BlockMask m_Mask;
            int[] m_RestrictedDeg;
            bool m_MPIself;

            BlockMsrMatrix m_OperatorMatrix;

            //Func<int[]> GetCellRestriction = null;
            //Func<BlockMsrMatrix> GetExtRows = null;

            public BlockMsrMatrix OperatorMatrix {
                get {
                    if(m_OperatorMatrix == null) {
                        if(!m_MPIself)
                            m_OperatorMatrix = m_Mask.GetSubBlockMatrix(m_Unrestricted.OperatorMatrix, m_Unrestricted.OperatorMatrix.MPI_Comm);
                        else
                            m_OperatorMatrix = m_Mask.GetSubBlockMatrix_MpiSelf(m_Unrestricted.OperatorMatrix);

                        //m_OperatorMatrix.SaveToTextFileSparse("Mtx-k" + m_RestrictedDeg[0] + ".txt");
                    }
                    return m_OperatorMatrix;
                }
            }

            class RestrictionMapping : BlockPartitioning, ICoordinateMapping {

                public RestrictionMapping(IBlockPartitioning basePart, MgOperatorRestriction ownerHack) : base(basePart) {
                    m_ownerHack = ownerHack;
                }

                private int GetNp(int p) {
                    int SpacDim = m_ownerHack.DgMapping.SpatialDimension;
                    switch(SpacDim) {
                        case 1: return p + 1;
                        case 2: return (p * p + 3 * p + 2) / 2;
                        case 3: return (p * p * p + 6 * p * p + 11 * p + 6) / 6;
                        default: throw new ArgumentOutOfRangeException("wtf?Spacialdim=1,2,3 expected");
                    }
                }


                MgOperatorRestriction m_ownerHack;

                public int NoOfVariables => m_ownerHack.m_Unrestricted.DgMapping.NoOfVariables;

                public int[] DgDegree => m_ownerHack.m_RestrictedDeg;

                public int NoOfExternalCells {
                    get {
                        if (m_ownerHack.m_MPIself)
                            return 0;
                        else 
                            return m_ownerHack.m_Unrestricted.DgMapping.NoOfExternalCells;
                    }
                }

                public int SpatialDimension => m_ownerHack.m_Unrestricted.DgMapping.SpatialDimension;

                public int NoOfLocalUpdatedCells => m_ownerHack.m_Unrestricted.DgMapping.NoOfLocalUpdatedCells;

                public int LocalCellCount => m_ownerHack.m_Unrestricted.DgMapping.LocalCellCount;

                public SpeciesId[] UsedSpecies => m_ownerHack.m_Unrestricted.DgMapping.UsedSpecies;

                public int GetLength(int jLoc) {
                    int NoOfVar = this.NoOfVariables;
                    int NoOfSpc = this.GetNoOfSpecies(jLoc);
                    int len = 0;
                    for(int i = 0; i < NoOfVar; i++) {
                        len += GetNp(this.m_ownerHack.m_RestrictedDeg[i]);
                    }
                    len *= NoOfSpc;
                    return len;
                }

                public int GetNoOfSpecies(int jCell) {
                    int _jCell = CellIdxTrafo(jCell);
                    return m_ownerHack.m_Unrestricted.DgMapping.GetNoOfSpecies(_jCell);
                }

                private int CellIdxTrafo(int jCell) {
                    int _jCell;
                    if(m_ownerHack.jLoc2Glob != null)
                        _jCell = m_ownerHack.jLoc2Glob[jCell];
                    else
                        _jCell = jCell;
                    return _jCell;
                }

                public int GetSpeciesIndex(int jCell, SpeciesId SId) {
                    int _jCell = CellIdxTrafo(jCell);
                    return m_ownerHack.m_Unrestricted.DgMapping.GetSpeciesIndex(_jCell, SId);
                }

                public long GlobalUniqueIndex(int ifld, int jCell, int jSpec, int n) {
                    return this.i0 + LocalUniqueIndex(ifld, jCell, jSpec, n);
                }

                public long GlobalUniqueIndex(int ifld, int jCell, int n) {
                    return this.i0 + LocalUniqueIndex(ifld, jCell, n);
                }

                public bool IsXDGvariable(int iVar) {
                    return m_ownerHack.m_Unrestricted.DgMapping.IsXDGvariable(iVar);
                }

                public int LocalUniqueIndex(int ifld, int jCell, int iSpec, int n) {
                    int ret = (int)(base.GetBlockI0(jCell) - base.i0);
                    int NoOfSpc = this.GetNoOfSpecies(jCell);
                    for(int i = 0; i < ifld; i++)
                        ret += GetNp(this.m_ownerHack.m_RestrictedDeg[i]) * NoOfSpc;
                    if(iSpec > 0)
                        ret += GetNp(this.m_ownerHack.m_RestrictedDeg[ifld]) * iSpec;
                    ret += n;
                    return ret;
                }

                public int LocalUniqueIndex(int ifld, int jCell, int n) {
                    int ret = (int)(base.GetBlockI0(jCell) - base.i0);
                    int NoOfSpc = this.GetNoOfSpecies(jCell);
                    for(int i = 0; i < ifld; i++)
                        ret += GetNp(this.m_ownerHack.m_RestrictedDeg[i]) * NoOfSpc;
                    ret += n;
                    return ret;
                }
            }

            RestrictionMapping m_rm;

            public ICoordinateMapping DgMapping {
                get {
                    if(m_rm == null) {
                        m_rm = new RestrictionMapping(this.OperatorMatrix._RowPartitioning, this);
                    }
                    return m_rm;
                }
            }

           
            public BlockMask BlkMask {
                get {
                    return m_Mask;
                }
            }
        }

        /// <summary>
        /// Optional restriction to a subset of cells/matrix blocks
        /// </summary>
        public Func<int[]> GetCellRestriction = null;

        /// <summary>
        /// Matrix containing external rows, if <see cref="GetCellRestriction"/> provides external cells
        /// </summary>
        public Func<BlockMsrMatrix> GetExtRows = null;


        /// <summary>
        /// Memory for external rows, if <see cref="GetCellRestriction"/> provides external cells;
        /// (see <see cref="MPIexchange{T}.Vector_Ext"/>, <see cref="MPIexchangeInverse{T}.Vector_Ext"/>)
        /// </summary>
        public Func<(double[] Xext, double[] RHSext)> GetExtMem = null;

        /// <summary>
        /// if true, the nested solver is working serially (only on the respective core), i.e. on the SELF communicator;
        /// </summary>
        public bool RestrictToMPIself = false;


        public int IterationsInNested {
            get {
                return LowerPSolver?.IterationsInNested ?? 0;
            }
        }

        public int ThisLevelIterations {
            get {
                return LowerPSolver?.ThisLevelIterations ?? 0;
            }
        }

        public bool Converged {
            get {
                return LowerPSolver?.Converged ?? false;
            }
        }

        public object Clone() {
            return new PRestriction() {
                LowerPSolver = this.LowerPSolver?.CloneAs()
            };
        }

        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        public void Init(MultigridOperator op) {
            InitImpl(op);
        }

        MgOperatorRestriction m_MgOperatorRestriction;

        void InitImpl(IOperatorMappingPair op) {
            this.Dispose();
            this.ResetStat();

            //int[] RestrictedDeg = op.DgMapping.DgDegree.Select(p => p - 1).ToArray();
            if (RestrictedDeg == null)
                throw new NotSupportedException("'RestrictedDeg' must be initialized before this object can be used.");
            if (RestrictedDeg.Length != op.DgMapping.NoOfVariables)
                throw new NotSupportedException("Configuration Error: length of 'RestrictedDeg' must match number of variables.");
            if(RestrictedDeg.Min() < 0)
                throw new NotSupportedException("Cannot reduce DG Degree below 0.");
            if ((this.GetExtRows != null) != (this.GetExtMem != null))
                throw new ArgumentException("Illegal configuration; if external matrix is provided, also external memory must be provided, and vise-versa.");

            for (int iVar = 0; iVar < op.DgMapping.NoOfVariables; iVar++) {
                if(RestrictedDeg[iVar] > op.DgMapping.DgDegree[iVar])
                    throw new NotSupportedException($"Configuration Error: 'RestrictedDeg[{iVar}] == {RestrictedDeg[iVar]}' but the maximum supported degree for the respective variable is {op.DgMapping.DgDegree[iVar]}..");
            }

            Len = op.DgMapping.LocalLength;
            m_MgOperatorRestriction = new MgOperatorRestriction(op, RestrictedDeg, this.GetCellRestriction, this.GetExtRows, this.RestrictToMPIself);
            LenSub = m_MgOperatorRestriction.DgMapping.TotalLength;

            if(LowerPSolver != null && LenSub > 0)
                LowerPSolver.Init(m_MgOperatorRestriction);

        }

        int Len;

        long LenSub;

        /// <summary>
        /// DG polynomial degrees in the restricted system
        /// </summary>
        public int[] RestrictedDeg;


        public void ResetStat() {
            if(LowerPSolver != null)
                LowerPSolver.ResetStat();
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            if (X.Count != Len)
                throw new ArgumentException("mismatch in solution vector length");
            if (B.Count != Len)
                throw new ArgumentException("mismatch in right-hand-side vector length");
      
            if (LenSub <= 0)
                return; // empty solver: might happen in certain IBM-calculations

            double[] xLo, bLo;
            double[] xExt, bExt;
            if (this.GetExtMem != null) {
                (xExt, bExt) = GetExtMem();
                xLo = m_MgOperatorRestriction.BlkMask.GetSubVec(xExt, X);
                bLo = m_MgOperatorRestriction.BlkMask.GetSubVec(bExt, B);
            } else {
                xExt = null;
                bExt = null;
                xLo = m_MgOperatorRestriction.BlkMask.GetSubVec(X);
                bLo = m_MgOperatorRestriction.BlkMask.GetSubVec(B);
            }


            if(LowerPSolver != null)
                LowerPSolver.Solve(xLo, bLo);

            if(xExt != null)
                m_MgOperatorRestriction.BlkMask.AccSubVec(xLo, xExt, X);
            else 
                m_MgOperatorRestriction.BlkMask.AccSubVec(xLo, X);
        }

        public long UsedMemory() {
            return LowerPSolver.UsedMemory();
        }

        public void Dispose() {
            if(LowerPSolver != null)
                LowerPSolver.Dispose();
        }

        /// <summary>
        /// solver for the restricted system
        /// </summary>
        public ISubsystemSolver LowerPSolver;


        public BlockMask RestrictionMask {
            get {
                return m_MgOperatorRestriction.BlkMask;
            }
        }

        /// <summary>
        /// matrix of sub-system
        /// </summary>
        public BlockMsrMatrix SubsysMatrix {
            get {
                return m_MgOperatorRestriction.OperatorMatrix;
            }
        }

    }
}
