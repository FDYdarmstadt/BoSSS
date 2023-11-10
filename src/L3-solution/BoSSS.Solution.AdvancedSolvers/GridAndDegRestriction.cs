using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {



    /// <summary>
    /// Restriction of a <see cref="IOperatorMappingPair"/> onto a sub-mesh with lower polynomial degree, 
    /// i.e. this object is not a solver on its own, but it performs the restriction of the operator matrix
    /// to a certain **subset of cells and modes**.
    /// - sub-mes restriction and degree restriction can be used independent,
    ///   i.e. it is also possible to use only one of them
    /// - the solution of the sub-system is performed by <see cref="LowerPSolver"/>
    /// - the polynomial degree hierarchy is defined by setting <see cref="RestrictedDeg"/>
    /// - mesh restriction is performed by <see cref="GetCellRestriction"/>
    /// </summary>
    public class GridAndDegRestriction : ISubsystemSolver {

        /*
        static int counter_objectId = 1;

        public int objectId;

        /// <summary>
        /// ctor
        /// </summary>
        public GridAndDegRestriction() {
            objectId = counter_objectId;
            counter_objectId++;
        }
        */

        /// <summary>
        /// Note: at this point, we do not support everything which a <see cref="SubBlockSelector"/> can support.
        /// This is, because not all selections, which <see cref="SubBlockSelector"/> can do are supported by <see cref="ICoordinateMapping"/> at this point.
        /// </summary>
        class MgOperatorRestriction : IOperatorMappingPair {


            public MgOperatorRestriction(IOperatorMappingPair Unrestricted, int[] RestrictedDeg, Func<int[]> __GetCellRestriction, Func<BlockMsrMatrix> __GetExtRows, bool __MPIself) {
                using (new FuncTrace()) {
                    var subBlock = new SubBlockSelector(Unrestricted.DgMapping);

                    if (__GetCellRestriction != null) {
                        jLoc2Glob = __GetCellRestriction();
                        subBlock.CellSelector(jLoc2Glob);
                    }
                    subBlock.SetModeSelector((int jCell, int iVar, int iSpc, int Deg) => Deg <= RestrictedDeg[iVar]);

                    m_MPIself = __MPIself;
                    m_RestrictedDeg = RestrictedDeg;
                    m_Unrestricted = Unrestricted;
                    m_Mask = new BlockMask(subBlock, __GetExtRows?.Invoke());

                    var idx = m_Mask.GlobalIndices_Internal;
                    for (int i = 1; i < idx.Count; i++) {
                        if (idx[i - 1] >= idx[i])
                            throw new Exception("ghdshj: " + Unrestricted.GetType());
                    }
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
                        default: throw new ArgumentOutOfRangeException("Unknown spatial dimension: " + SpacDim);
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

                public int NoOfLocalUpdatedCells {
                    get {
                        if (m_ownerHack.jLoc2Glob != null)
                            return m_ownerHack.jLoc2Glob.Length;
                        else
                            return m_ownerHack.m_Unrestricted.DgMapping.NoOfLocalUpdatedCells;
                    }
                }
                public int LocalCellCount {
                    get {
                        if (m_ownerHack.jLoc2Glob != null) {
                            if (this.MpiSize == 1)
                                return NoOfLocalUpdatedCells;
                            else
                                throw new NotImplementedException();
                        } else {
                            return m_ownerHack.m_Unrestricted.DgMapping.LocalCellCount;
                        }
                    }
                }

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
        /// Optional restriction to a subset of cells/matrix blocks:
        /// returns local indices of cells which should be in the selection;
        /// if null, all cells are selected
        /// </summary>
        public Func<int[]> GetCellRestriction = null;

        /// <summary>
        /// Matrix containing external rows, if <see cref="GetCellRestriction"/> provides external cells;
        /// these rows would then be copied to this processor
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

        bool? m_converged;

        public bool Converged {
            get {
                return m_converged ?? false;
            }
        }

        public object Clone() {
            return new GridAndDegRestriction() {
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
            using (var tr = new FuncTrace()) {
                this.Dispose();
                this.ResetStat();

                //int[] RestrictedDeg = op.DgMapping.DgDegree.Select(p => p - 1).ToArray();
                if (RestrictedDeg == null)
                    throw new NotSupportedException("'RestrictedDeg' must be initialized before this object can be used.");
                if (RestrictedDeg.Length != op.DgMapping.NoOfVariables)
                    throw new NotSupportedException("Configuration Error: length of 'RestrictedDeg' must match number of variables.");
                if (RestrictedDeg.Min() < 0)
                    throw new NotSupportedException("Cannot reduce DG Degree below 0.");
                if ((this.GetExtRows != null) != (this.GetExtMem != null))
                    throw new ArgumentException("Illegal configuration; if external matrix is provided, also external memory must be provided, and vise-versa.");

                for (int iVar = 0; iVar < op.DgMapping.NoOfVariables; iVar++) {
                    if (RestrictedDeg[iVar] > op.DgMapping.DgDegree[iVar])
                        throw new NotSupportedException($"Configuration Error: 'RestrictedDeg[{iVar}] == {RestrictedDeg[iVar]}' but the maximum supported degree for the respective variable is {op.DgMapping.DgDegree[iVar]}..");
                }
                if (op.DgMapping.MPI_Comm != op.OperatorMatrix.MPI_Comm)
                    throw new ApplicationException("mismatch of MPI communicators");


                Len = op.DgMapping.LocalLength;
                m_MgOperatorRestriction = new MgOperatorRestriction(op, RestrictedDeg, this.GetCellRestriction, this.GetExtRows, this.RestrictToMPIself);
                LenSub = m_MgOperatorRestriction.DgMapping.TotalLength;

                tr.Info("MPI rank of DG mapping " + m_MgOperatorRestriction.DgMapping.MpiRank + " of " + m_MgOperatorRestriction.DgMapping.MpiSize);
                if (m_MgOperatorRestriction.DgMapping.MPI_Comm != m_MgOperatorRestriction.OperatorMatrix.MPI_Comm)
                    throw new ApplicationException("mismatch of MPI communicators");

                if (m_MgOperatorRestriction.DgMapping.MPI_Comm != m_MgOperatorRestriction.OperatorMatrix.MPI_Comm)
                    throw new ApplicationException("mismatch of MPI communicators");
                

                if (this.GetCellRestriction != null) {
                    // Note, fk, 08jul22: the following exception would wrongfully occur when there are completely empty cells;
                    // then, empty blocks are eliminated from the matrix and the number of matrix blocks might be less than the number of cells.

                    //if (m_MgOperatorRestriction.OperatorMatrix._RowPartitioning.LocalNoOfBlocks != this.GetCellRestriction().Length)
                    //    throw new Exception($"Blocks in Operator Matrix: {m_MgOperatorRestriction.OperatorMatrix._RowPartitioning.LocalNoOfBlocks}," +
                    //        $" number of cell restriction: {this.GetCellRestriction().Length}");

                    //Console.Error.WriteLine(
                    //    $"Rank{op.DgMapping.MpiRank}of{op.DgMapping.MpiSize}: " +
                    //    $"Blocks in Operator Matrix: {m_MgOperatorRestriction.OperatorMatrix._RowPartitioning.LocalNoOfBlocks}," +
                    //        $" number of cell restriction: {this.GetCellRestriction().Length}");
                }

                if (LowerPSolver != null && LenSub > 0)
                    LowerPSolver.Init(m_MgOperatorRestriction);
            }
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


            if (LowerPSolver != null) {
                double normRHS = bLo.MPI_L2Norm(m_MgOperatorRestriction.OperatorMatrix.MPI_Comm);
                //csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIrnk);
                //Console.Error.WriteLine($"  ++++ R{MPIrnk} Restriction: RHS is {normRHS}; input RHS norm is {B.MPI_L2Norm(m_MgOperatorRestriction.OperatorMatrix.MPI_Comm)}");
                if (normRHS <= 0) {
                    //B.SaveToTextFile("Br" + MPIrnk + "-" + Guid.NewGuid().ToString() + ".txt", m_MgOperatorRestriction.OperatorMatrix.MPI_Comm);
                    //X.ClearEntries();

                    m_converged = true;
                    return;
                } else {
                    LowerPSolver.Solve(xLo, bLo);
                    m_converged = LowerPSolver.Converged;
                }
            }

            if(xExt != null)
                m_MgOperatorRestriction.BlkMask.AccSubVec(xLo, xExt, X);
            else 
                m_MgOperatorRestriction.BlkMask.AccSubVec(xLo, X);
        }

        public long UsedMemory() {
            return LowerPSolver?.UsedMemory() ?? 0;
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

        /// <summary>
        /// Access to sub-system
        /// </summary>
        public IOperatorMappingPair OperatorRestriction {
            get {
                return m_MgOperatorRestriction;
            }
        }

    }
}
