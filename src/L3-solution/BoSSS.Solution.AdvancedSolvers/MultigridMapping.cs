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
using System.Text;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using MPI.Wrappers;
using System.Numerics;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Comm;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// For each aggregation grid level, this mapping defines a bijection between variable index,
    /// DG mode and aggregation cell index and a unique index.
    /// </summary>
    public class MultigridMapping : IBlockPartitioning, ICoordinateMapping {

        /// <summary>
        /// Base grid on which the problem is defined.
        /// </summary>
        public GridData GridData {
            get {
                return (GridData)(this.ProblemMapping.GridDat);
            }
        }

        /// <summary>
        /// Coordinate mapping on the original grid.
        /// </summary>
        public UnsetteledCoordinateMapping ProblemMapping {
            get;
            //private set;
        }


        /// <summary>
        /// Aggregation basis on this level, for each variable.
        /// </summary>
        public AggregationGridBasis[] AggBasis {
            get;
            //private set;
        }

        /// <summary>
        /// Number of variables
        /// </summary>
        public int NoOfVariables {
            get {
                return AggBasis.Length;
            }
        }


        /// <summary>
        /// aggregation grid on this level
        /// </summary>
        public AggregationGridData AggGrid {
            get {
                return AggBasis[0].AggGrid;
            }
        }
        
        /// <summary>
        /// Number of DOF's stored on this  MPI process
        /// </summary>
        public int LocalLength {
            get {
                int Jup = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                if (this.m_i0 != null) {
                    return checked((int)(this.m_i0[Jup]));
                } else {
                    return Jup * this.MaximalLength;
                }
            }
        }
        

        /// <summary>
        /// Partitioning of the vector among MPI processes.
        /// </summary>
        public IPartitioning Partitioning {
            get;
        }

        /// <summary>
        /// Total number of DOF's over all MPI processes.
        /// </summary>
        public long TotalLength {
            get {
                return this.Partitioning.TotalLength;
            }
        }

        int[] m_DgDegree;

        /// <summary>
        /// For each DG basis in the <see cref="ProblemMapping"/>, the DG polynomial degree
        /// which is used on this level.
        /// </summary>
        public int[] DgDegree {
            get {
                return m_DgDegree.CloneAs();
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="__ProblemMapping">
        /// Mapping of original problem, equal to <see cref="ProblemMapping"/>.
        /// </param>
        /// <param name="__aggGrdB">
        /// Sequence of aggregation grid DG basis objects, correlates to original mapping
        /// </param>
        /// <param name="DgDegrees"></param>
        public MultigridMapping(UnsetteledCoordinateMapping __ProblemMapping, AggregationGridBasis[] __aggGrdB, int[] DgDegrees) {
            using (var tr =new FuncTrace()) {
                // check args
                // ===========
                tr.Info("Entering multigrid mapping " + __aggGrdB[0] + " Basis length " + __aggGrdB.Length);

                if (__aggGrdB.Length != __ProblemMapping.BasisS.Count)
                    throw new ArgumentException("Mismatch between number of multigrid basis objects and number of variables in problem mapping.");

                var mgGrid = __aggGrdB[0].AggGrid;
                for (int iVar = 0; iVar < __aggGrdB.Length; iVar++) {
                    if (__ProblemMapping.BasisS[iVar] is BoSSS.Foundation.XDG.XDGBasis) {
                        if (!(__aggGrdB[iVar] is XdgAggregationBasis))
                            throw new ArgumentException();

                    } else if (__ProblemMapping.BasisS[iVar] is BoSSS.Foundation.Basis) {
                        //if((__aggGrdB[iVar] is XdgAggregationBasis))
                        //    throw new ArgumentException();
                    }

                    if (!object.ReferenceEquals(mgGrid, __aggGrdB[iVar].AggGrid))
                        throw new ArgumentException("Basis object must all be defined on the same grid.");
                }

                // find basis with maximum degree
                // ==============================
                this.ProblemMapping = __ProblemMapping;
                //this.MaxBasis = ProblemMapping.BasisS.ElementAtMax(basis => basis.Degree);
                this.m_DgDegree = DgDegrees.CloneAs();
                //if(!this.MaxBasis.IsSubBasis(__aggGrdB.DGBasis)) {
                //    throw new ArgumentException("Basis on aggregation grid is insufficient;");
                //}

                if (m_DgDegree.Length != __ProblemMapping.BasisS.Count()) {
                    throw new ArgumentException("Wrong number of DG degrees.");
                }
                for (int i = 0; i < m_DgDegree.Length; i++) {
                    if (m_DgDegree[i] < 0)
                        throw new ArgumentException("DG degree must be greater of equal to 0.");
                    if (m_DgDegree[i] > __ProblemMapping.BasisS[i].Degree)
                        throw new ArgumentException("DG degree on sub-level can not exceed DG degree of the original problem.");
                }
                
                // create basis for this level
                // ===========================
                this.AggBasis = __aggGrdB;

                // min/max length
                // ==============
                {
                    
                    int Smin = 0;
                    int Smax = 0;
                    int Nofields = this.m_DgDegree.Length;
                    
                    for (int ifld = 0; ifld < Nofields; ifld++) {
                        Smin += this.AggBasis[ifld].GetMinimalLength(this.m_DgDegree[ifld]);
                        Smax += this.AggBasis[ifld].GetMaximalLength(this.m_DgDegree[ifld]);
                    }
                    this.MinimalLength = Smin.MPIMin();
                    this.MaximalLength = Smax.MPIMax();
                }

                // offsets
                // =======
                if (this.MinimalLength != this.MaximalLength) {
                    int JAGGloc = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    int JAGGtot = this.AggGrid.iLogicalCells.Count;                
                    long[] __i0Tmp = new long[JAGGtot];

                    HashSet<int> BlockLen = new HashSet<int>();

                    int LL = 0;
                    for (int jag = 0; jag < JAGGloc; jag++) {
                        int S = 0;
                        for (int i = 0; i < m_DgDegree.Length; i++) {
                            S += this.AggBasis[i].GetLength(jag, m_DgDegree[i]);
                        }
                        //if (S == 0) Console.WriteLine("zero at j="+jag);
                        if (jag < JAGGloc - 1)
                            __i0Tmp[jag + 1] = __i0Tmp[jag] + S;
                        else
                            LL = checked((int)(__i0Tmp[jag] + S));

                        BlockLen.Add(S);
                    }
                    this.Partitioning = new Partitioning(LL);
                    long i0Part = Partitioning.i0;
                    m_i0 = new int[JAGGtot + 1];
                    for (int jag = 0; jag < JAGGloc; jag++) { // loop over local cells
                        m_i0[jag] = (int)__i0Tmp[jag]; // store local i0 index
                        __i0Tmp[jag] += i0Part; // convert to global index
                    }
                    m_i0[JAGGloc] = LL;
                    __i0Tmp.MPIExchange(this.AggGrid);
                    
                    // compute global cell i0's in the external range
                    m_i0_ExtGlob = new long[JAGGtot - JAGGloc];
                    Array.Copy(__i0Tmp, JAGGloc, m_i0_ExtGlob, 0, JAGGtot - JAGGloc);
                    
                    // compute local cell i0's in the external range (very confusing)
                    for (int jag = JAGGloc; jag < JAGGtot; jag++) {
                        int S = 0;
                        for (int i = 0; i < m_DgDegree.Length; i++) {
                            S += this.AggBasis[i].GetLength(jag, m_DgDegree[i]);
                        }
                        m_i0[jag + 1] = this.m_i0[jag] + S;

                        BlockLen.Add(S);
                    }

                    // build look-up for block types
                    m_Len2SublockType = new Dictionary<int, int>();
                    m_Subblk_i0 = new int[BlockLen.Count][];
                    m_SubblkLen = new int[BlockLen.Count][];
                    int type = 0;
                    foreach (int S in BlockLen) {
                        m_Subblk_i0[type] = new int[] { 0 };
                        m_SubblkLen[type] = new int[] { S };
                        m_Len2SublockType.Add(S, type);
                        type++;
                    }

//#if DEBUG
                    for (int jag = JAGGloc; jag < JAGGtot; jag++) { // loop over external cells
                        int S = 0;
                        for (int i = 0; i < m_DgDegree.Length; i++) {
                            S += this.AggBasis[i].GetLength(jag, m_DgDegree[i]);
                        }

                        if (S > 0 && Partitioning.IsInLocalRange(m_i0_ExtGlob[jag - JAGGloc]) == true)
                            throw new ApplicationException($"Multigrid level {this.AggGrid.MgLevel}: 'm_i0_ExtGlob[{jag - JAGGloc}]' seems to be in local range: got {m_i0_ExtGlob[jag - JAGGloc]}, must be outside of {Partitioning.i0}--{Partitioning.iE}; number of locally updated cells is {JAGGloc}, total number is {JAGGtot}");
                    }
//#endif

                } else {
                    m_Subblk_i0 = new int[][] { new int[] { 0 } };
                    m_SubblkLen = new int[][] { new int[] { this.MaximalLength } };

                    this.Partitioning = new Partitioning(this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells * this.MaximalLength);
                }


                Debug.Assert(Partitioning != null);
                Debug.Assert(Partitioning.LocalLength == this.LocalLength);
                
            }
        }

        Dictionary<int, int> m_Len2SublockType;

        int[][] m_Subblk_i0;
        int[][] m_SubblkLen;

        /// <summary>
        /// For each aggregation cell, the local vector index of the first DG coordinate in this cell.
        /// - index: local aggregation cell index.
        /// </summary>
        int[] m_i0;

        /// <summary>
        /// For each external aggregation cell, the global vector index of the first DG coordinate in this cell.
        /// - index: external aggregation cell index, minus local number of aggregation cells
        /// </summary>
        long[] m_i0_ExtGlob;
        
        /// <summary>
        /// Number of degrees-of-freedom in cell <paramref name="jCell"/> , for the <paramref name="iVar"/>-th variable.
        /// </summary>
        public int GetLengthForVar(int jCell, int iVar) {
            int Nofields = this.m_DgDegree.Length;
            int S = 0;
            for (int ifld = 0; ifld < Nofields; ifld++) {
                int pField = this.m_DgDegree[ifld];
                S += this.AggBasis[ifld].GetLength(jCell, pField);
            }
            return S;
        }


        /// <summary>
        /// For this mapping, the maximum DOF used per cell over all cells.
        /// </summary>
        public int MaximalLength {
            get;
            private set;
        }

        /// <summary>
        /// For this mapping, the minimum DOF used per cell over all cells.
        /// </summary>
        public int MinimalLength {
            get;
            private set;
        }


        public long TotalNoOfBlocks {
            get {
                return AggGrid.CellPartitioning.TotalLength;
            }
        }

        public int LocalNoOfBlocks {
            get {
                return AggGrid.CellPartitioning.LocalLength;
            }
        }

        /// <summary>
        /// Gets the first block on this process
        /// </summary>
        public long FirstBlock {
            get {
                return AggGrid.CellPartitioning.i0;
            }
        }
        /// <summary>
        /// Tests, if Max and Minblocksize are equal
        /// </summary>
        public bool AllBlockSizesEqual {
            get {
                return (MaximalLength == MinimalLength);
            }
        }

        public MPI_Comm MPI_Comm {
            get {
                return Partitioning.MPI_Comm;
            }
        }

        public int MpiSize {
            get {
                return Partitioning.MpiSize;
            }
        }

        public int MpiRank {
            get {
                return Partitioning.MpiRank;
            }
        }

        public long i0 {
            get {
                return Partitioning.i0;
            }
        }

        public long iE {
            get {
                return Partitioning.iE;
            }
        }

        public bool IsMutable {
            get {
                return false;
            }
        }

        public int NoOfExternalCells {
            get {
                return AggGrid.iLogicalCells.NoOfExternalCells;
            }
        }

        public int SpatialDimension {
            get {
                return AggGrid.SpatialDimension;
            }
        }

        public int NoOfLocalUpdatedCells {
            get {
                return AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            }
        }

        public int LocalCellCount {
            get {
                return NoOfLocalUpdatedCells + NoOfExternalCells;
            }
        }

        public int LocalUniqueIndex(int ifld, int jCell, int n) {
            Debug.Assert(ifld >= 0 && ifld < this.m_DgDegree.Length);
            Debug.Assert(jCell >= 0 && jCell < (this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells + this.AggGrid.iLogicalCells.NoOfExternalCells));
            Debug.Assert((n >= 0 && (n < this.AggBasis[ifld].GetLength(jCell, this.m_DgDegree[ifld]))) || this.AggBasis[ifld].GetLength(jCell, this.m_DgDegree[ifld])==0); // 0<= n < n_max oder empty cell in case of IBM

            int S;
            if(this.m_i0 != null) {

                S = this.m_i0[jCell];
            } else {
                Debug.Assert(this.MaximalLength == this.MinimalLength);
                S = jCell * this.MaximalLength;
            }
            //if (ilPSP.Environment.MPIEnv.MPI_Rank == 1) {
                //Console.WriteLine("proc:{0}, iVar:{1}, iCell:{2}, {3}", ilPSP.Environment.MPIEnv.MPI_Rank, ifld, jCell, S);
            //}

            
            for (int iF = 0; iF < ifld; iF++)
                S += this.AggBasis[iF].GetLength(jCell, this.m_DgDegree[iF]);
            S += n;
            return S;
        }

        public int LocalUniqueIndex(int ifld, int jCell, int jSpec, int n) {
            int Np_tot = this.AggBasis[ifld].GetLength(jCell, this.m_DgDegree[ifld]);
            int NoOfSpec = AggBasis[ifld].GetNoOfSpecies(jCell);
            int Np_Spec = Np_tot / NoOfSpec;
            int n_agg=jSpec*Np_Spec + n;

            return LocalUniqueIndex(ifld, jCell, n_agg);
        }

        public long GlobalUniqueIndex(int ifld, int jCell, int jSpec, int n) {
            int Np_tot = this.AggBasis[ifld].GetLength(jCell, this.m_DgDegree[ifld]);
            int NoOfSpec = AggBasis[ifld].GetNoOfSpecies(jCell);
            int Np_Spec = Np_tot / NoOfSpec;
            int n_agg = jSpec * Np_Spec + n;

            return GlobalUniqueIndex(ifld, jCell, n_agg);
        }

        public long GlobalUniqueIndex(int ifld, int jCell, int n) {
            Debug.Assert(ifld >= 0 && ifld < this.m_DgDegree.Length);
            Debug.Assert(jCell >= 0 && jCell < (this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells + this.AggGrid.iLogicalCells.NoOfExternalCells));        
            Debug.Assert(n >= 0 && n < this.AggBasis[ifld].GetLength(jCell, this.m_DgDegree[ifld]));

            int Jup = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            if(jCell < Jup) {
                return this.i0 + LocalUniqueIndex(ifld, jCell, n);
            } else {
                if(this.m_i0 == null) {
                    long jGlb = this.AggGrid.iParallel.GlobalIndicesExternalCells[jCell - Jup];
                    long S = jGlb * this.MaximalLength;
                    for(int iF = 0; iF < ifld; iF++)
                        S += this.AggBasis[iF].GetLength(jCell, this.m_DgDegree[iF]);
                    S += n;
                    return S;
                } else {
                    long S = m_i0_ExtGlob[jCell - Jup];
                    for (int iF = 0; iF < ifld; iF++)
                        S += this.AggBasis[iF].GetLength(jCell, this.m_DgDegree[iF]);
                    S += n;
                    return S;
                }

            }

        }

        /// <summary>
        /// Prolongation/Injection operator to finer grid level.
        /// </summary>
        public BlockMsrMatrix GetProlongationOperator(MultigridMapping finerLevel) {
            using (new FuncTrace()) {

                // Argument checking
                // =================

                if (!object.ReferenceEquals(finerLevel.AggGrid, this.AggGrid.ParentGrid))
                    throw new ArgumentException("Only prolongation/injection to next level is supported.");
                if (finerLevel.AggBasis.Length != this.AggBasis.Length) {
                    throw new ArgumentException("");
                }
                int NoOfVar = this.AggBasis.Length;

                MultidimensionalArray[][] InjOp = new MultidimensionalArray[NoOfVar][];
                AggregationGridBasis[] B = new AggregationGridBasis[NoOfVar];
                bool[] useX = new bool[NoOfVar];
                int[] DegreeS = new int[NoOfVar];
                int[] DegreeSfine = new int[NoOfVar];
                
                for (int iVar = 0; iVar < NoOfVar; iVar++) {
                    InjOp[iVar] = this.AggBasis[iVar].InjectionOperator;
                    B[iVar] = AggBasis[iVar];
                    DegreeS[iVar] = this.DgDegree[iVar];
                    DegreeSfine[iVar] = finerLevel.DgDegree[iVar];
                    if (DegreeSfine[iVar] < DegreeS[iVar])
                        throw new ArgumentException("Lower DG degree on finer grid is not supported by this method ");
                    useX[iVar] = this.AggBasis[iVar] is XdgAggregationBasis;
                    if (useX[iVar] != (finerLevel.AggBasis[iVar] is XdgAggregationBasis))
                        throw new ArgumentException("XDG / DG mismatch between this and finer level for " + iVar + "-th variable.");
                }

                XdgAggregationBasis XB = null;
                XdgAggregationBasis XBf = null;
                int[][,] spcIdxMap = null;
                SpeciesId[][] spc = null;
                //SpeciesId[][] spcf = null;
                for(int iVar = 0; iVar < NoOfVar; iVar++) {
                    if(useX[iVar]) {
                        XB = (XdgAggregationBasis)(B[iVar]);
                        XBf = (XdgAggregationBasis)(finerLevel.AggBasis[iVar]);
                        spcIdxMap = XB.SpeciesIndexMapping;
                        spc = XB.AggCellsSpecies;
                        //spcf = XBf.AggCellsSpecies;
                        break;
                    }
                }

                int[] Np = this.AggBasis[0].GetNp();
                int[] Np_fine = finerLevel.AggBasis[0].GetNp();
                


                // create matrix
                // =============

                // init retval
                var PrlgMtx = new BlockMsrMatrix(finerLevel, this);

                int[][] C2F = this.AggGrid.jCellCoarse2jCellFine;
                int JCoarse = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                //Debug.Assert((JCoarse == C2F.Length) || ());
                for(int jc = 0; jc < JCoarse; jc++) { // loop over coarse cells...
                    int[] AggCell = C2F[jc];
                    int I = AggCell.Length;
                    
                    for(int iVar = 0; iVar < NoOfVar; iVar++) {
                        int DgDeg = DegreeS[iVar];
                        int DgDegF = DegreeSfine[iVar];
                        MultidimensionalArray Inj_iVar_jc = InjOp[iVar][jc];
                        Debug.Assert(Inj_iVar_jc.GetLength(0) == I);

                        bool useX_iVar = false;
                        if(useX[iVar]) {
                            if (spcIdxMap[jc] != null)
                                useX_iVar = true;
                        }


                        if(useX_iVar) {
                            //throw new NotImplementedException("todo");
                            
                            int NoOfSpc = XB.GetNoOfSpecies(jc);
                            int Np_col = Np[DgDeg];
                            Debug.Assert(Np_col*NoOfSpc == B[iVar].GetLength(jc, DgDeg));

                            for(int iSpc = 0; iSpc < NoOfSpc; iSpc++) { // loop over species
                                SpeciesId spc_jc_i = spc[jc][iSpc];

                                long Col0 = this.GlobalUniqueIndex(iVar, jc, Np_col*iSpc);


                                for (int i = 0; i < I; i++) { // loop over finer cells
                                    int jf = AggCell[i];

                                    int iSpc_Row = XBf.GetSpeciesIndex(jf, spc_jc_i);
                                    if(iSpc_Row < 0) {
                                        // nothing to do
                                        continue;
                                    }

                                    int Np_row = Np_fine[DgDegF];
                                    Debug.Assert(Np_row*XBf.GetNoOfSpecies(jf)  == finerLevel.AggBasis[iVar].GetLength(jf, DgDegF));

                                    long Row0 = finerLevel.GlobalUniqueIndex(iVar, jf, Np_row*iSpc_Row);

                                    PrlgMtx.AccBlock(Row0, Col0, 1.0, Inj_iVar_jc.ExtractSubArrayShallow(new[] { i, 0, 0 }, new[] { i - 1, Np_row - 1, Np_col - 1 }));
                                }
                            }
                            
                        } else {
                            // ++++++++++++++++++
                            // standard DG branch
                            // ++++++++++++++++++

                            int Np_col = Np[DgDeg];
                            Debug.Assert(Np_col == B[iVar].GetLength(jc, DgDeg));
                            long Col0 = this.GlobalUniqueIndex(iVar, jc, 0);

                            for(int i = 0; i < I; i++) { // loop over finer cells
                                int jf = AggCell[i];
                                int Np_row = Np_fine[DgDegF];
                                Debug.Assert(Np_row == finerLevel.AggBasis[iVar].GetLength(jf, DgDegF));

                                long Row0 = finerLevel.GlobalUniqueIndex(iVar, jf, 0);

                                PrlgMtx.AccBlock(Row0, Col0, 1.0, Inj_iVar_jc.ExtractSubArrayShallow(new[] { i, 0, 0 }, new[] { i - 1, Np_row - 1, Np_col - 1 }));
                                
                            }
                        }
                    }
                }


                // return
                // ======

                return PrlgMtx;
            }
        }


        /// <summary>
        /// Prolongation or Restriction from *any* other level to this level.
        /// </summary>
        /// <param name="otherLevel">
        /// Other level, from which one wants to prolongate/restrict.
        /// </param>
        /// <returns>
        /// A matrix, where 
        /// - row correspond to this mapping
        /// - columns correspond to <paramref name="otherLevel"/>
        /// If the other level is coarser, this is a prolongation; if the other level is finer, it is the restriction in the L2-sense,
        /// for standard DG; for XDG, in some other norm  determined by the cut-cell shape.
        /// </returns>
        public BlockMsrMatrix FromOtherLevelMatrix(MultigridMapping otherLevel) {
            using (new FuncTrace()) {
                BlockMsrMatrix PrlgMtx;
                {
                    PrlgMtx = new BlockMsrMatrix(otherLevel, otherLevel.ProblemMapping);
                    for (int ifld = 0; ifld < otherLevel.AggBasis.Length; ifld++)
                        otherLevel.AggBasis[ifld].GetRestrictionMatrix(PrlgMtx, otherLevel, ifld);  // prolongate from the other level to the full grid
                    PrlgMtx = PrlgMtx.Transpose();
                }

                BlockMsrMatrix RestMtx;
                {
                    RestMtx = new BlockMsrMatrix(this, this.ProblemMapping);
                    for (int ifld = 0; ifld < this.AggBasis.Length; ifld++)
                        this.AggBasis[ifld].GetRestrictionMatrix(RestMtx, this, ifld);  //     ... and restrict to this level          
                }
                var result = BlockMsrMatrix.Multiply(RestMtx, PrlgMtx);
#if DEBUG
            {
                var resultT = result.Transpose();
                BlockMsrMatrix ShoudBeId;
                if(result.RowPartitioning.TotalLength < result.ColPartition.TotalLength)
                    ShoudBeId = BlockMsrMatrix.Multiply(result, resultT);
                else
                    ShoudBeId = BlockMsrMatrix.Multiply(resultT, result);

                ShoudBeId.AccEyeSp(-1.0);

                double ShouldBeID_Norm = ShoudBeId.InfNorm();
                Debug.Assert(ShouldBeID_Norm < 1.0e-8);
                //Console.WriteLine("Id norm {0} \t (level {1})", ShouldBeID_Norm, this.AggGrid.MgLevel);
            }
#endif
                return result;
            }
        }

        /// <summary>
        /// Returns global unique indices which correlate to a certain sub-set of this mapping's 
        /// basises (<see cref="UnsetteledCoordinateMapping.BasisS"/>).
        /// </summary>
        /// <param name="Fields">
        /// Indices into <see cref="AggBasis"/>
        /// </param>
        public long[] GetSubvectorIndices(params int[] Fields) {
            ilPSP.MPICollectiveWatchDog.Watch();
            var map = this.ProblemMapping;

            var _BasisS = map.BasisS.ToArray();
            for(int iVar = 0; iVar < _BasisS.Length; iVar++) {
                Basis b = _BasisS[iVar];

                if(b is BoSSS.Foundation.XDG.XDGBasis) {
                    if(!b.IsSubBasis(((XdgAggregationBasis)(this.AggBasis[iVar])).XDGBasis))
                        throw new ArgumentException();
                } else {
                    if(!b.IsSubBasis(this.AggBasis[iVar].DGBasis))
                        throw new ArgumentException();
                }
            }
            var ag = this.AggGrid;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            List<long> R = new List<long>();
            Partitioning p = new Partitioning(JAGG);
            long j0_aggCell = p.i0;

            if(this.MaximalLength == this.MinimalLength) {
                // ++++++++++++++++++++++++++++++
                // case: constant length per cell
                // ++++++++++++++++++++++++++++++


                Debug.Assert(this.m_DgDegree.Length == this.AggBasis.Length);
                int[] DofVar = new int[this.m_DgDegree.Length];
                for (int i = 0; i < DofVar.Length; i++)
                    DofVar[i] = AggBasis[i].GetLength(0, m_DgDegree[i]);

                int L = DofVar.Sum();
                int[] Offset = new int[Fields.Length];
                for(int i = 0; i < Fields.Length; i++) {
                    int iField = Fields[i];
                    if(iField < 0 || iField > _BasisS.Length)
                        throw new IndexOutOfRangeException();
                    Offset[i] = iField.ForLoop(k => DofVar[k]).Sum();
                }
                

                for(int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over aggregate cells

                    for(int i = 0; i < Fields.Length; i++) {
                        int iField = Fields[i];
                        int N = DofVar[iField];

                        N = ag.iLogicalCells.AggregateCellToParts[jAgg].Max(j => DofVar[iField]);

                        long i0 = L * (j0_aggCell + jAgg) + Offset[i];

                        for(int n = 0; n < N; n++) {
                            long iX = i0 + n;
                            R.Add(iX);
                            Debug.Assert(this.Partitioning.IsInLocalRange(iX));
                            Debug.Assert(iX - j0_aggCell*L == this.LocalUniqueIndex(iField, jAgg, n));
                        }
                    }
                }
            } else {
                // ++++++++++++++++++++++++++++++
                // case: variable length per cell
                // ++++++++++++++++++++++++++++++

                int Nofields = this.m_DgDegree.Length;

                //int[] Nfld = new int[Nofields];
                //int[] Ofst = new int[Nofields];

                long i0Proc = this.Partitioning.i0;

                for(int jAgg = 0; jAgg < JAGG; jAgg++) {

                    //for(int ifld = 0; ifld < Nofields; ifld++) {
                    //    int pField = this.m_DgDegree[ifld];
                    //    Nfld[ifld] = this.AggBasis[ifld].GetLength(jAgg, pField);
                    //    if(ifld > 0)
                    //        Ofst[ifld] = Ofst[ifld - 1] + Nfld[ifld - 1];
                    //}
                    //int Ncell = Ofst[Nofields - 1] + Nfld[Nofields - 1];
                        

                    for(int i = 0; i < Fields.Length; i++) {
                        int iField = Fields[i];
                        int pField = this.m_DgDegree[iField];
                        int Nfld = this.AggBasis[iField].GetLength(jAgg, pField);


                        

                        //int N_iField = Nfld[iField];
                        //int N0 = Ofst[iField];

                        int i0Loc = this.LocalUniqueIndex(iField, jAgg, 0);

                        for(int n = 0; n < Nfld; n++) {
                            long iX = i0Proc + i0Loc + n;
                            R.Add(iX);
                            Debug.Assert(this.Partitioning.IsInLocalRange(iX));
                            Debug.Assert(iX - i0Proc == this.LocalUniqueIndex(iField, jAgg, n));
                        }
                    }
                    
                }
            }
            
            return R.ToArray();
        }

        

        /// <summary>
        /// Gets index vecor of all ghost cells available on this proc
        /// </summary>
        /// <param name="Fields"></param>
        /// <returns></returns>
        public long[] GetSubvectorIndices_Ext(params int[] Fields) {
            int Locoffset = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int[] LocCellIdxExt = this.AggGrid.iLogicalCells.NoOfExternalCells.ForLoop(i=>i + Locoffset);
            List<long> R = new List<long>();
            foreach(int jCell in LocCellIdxExt) {
                foreach (int fld in Fields) {
                    long i0_Block = GlobalUniqueIndex(fld, jCell, 0);
                    int N = this.AggBasis[fld].GetLength(jCell, this.DgDegree[fld]);
                    for(int i = 0; i < N; i++) {
                        R.Add(i0_Block + i);
                    }
                }
            }
            return R.ToArray();
        }

        /// <summary>
        /// Returns global unique indices which correlate to a certain species and basises.
        /// </summary>
        /// <param name="Fields">
        /// Indices into <see cref="AggBasis"/>
        /// </param>
        /// <param name="Species">
        /// </param>
        /// <returns>a list of global (over all MPI processes) unique indices.</returns>
        public long[] GetSubvectorIndices(Foundation.XDG.SpeciesId Species, params int[] Fields) {
            var map = this.ProblemMapping;

            var _BasisS = map.BasisS.ToArray();
            //XdgAggregationBasis XaggBasis;
            //if(this.AggBasis is XdgAggregationBasis) {
            //    // super-quick and super-dirty fix
            //    XaggBasis = (XdgAggregationBasis)(this.AggBasis);

            //    foreach(var b in _BasisS) {
            //        if(!b.IsSubBasis(XaggBasis.XDGBasis))
            //            throw new ArgumentException();
            //    }
            //} else {
            //    throw new NotSupportedException();
            //}
            XdgAggregationBasis[] XaggBasis = new XdgAggregationBasis[_BasisS.Length];
            for(int iVar = 0; iVar < _BasisS.Length; iVar++) {
                Basis b = _BasisS[iVar];

                if(b is BoSSS.Foundation.XDG.XDGBasis) {
                    if(!b.IsSubBasis(((XdgAggregationBasis)(this.AggBasis[iVar])).XDGBasis))
                        throw new ArgumentException();
                    XaggBasis[iVar] = ((XdgAggregationBasis)(this.AggBasis[iVar]));
                } else {
                    //if(!b.IsSubBasis(this.AggBasis[iVar].DGBasis))
                    //    throw new ArgumentException();
                    throw new NotSupportedException();
                }
            }

            var ag = this.AggGrid;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            Partitioning p = new Partitioning(JAGG);
            List<long> R = new List<long>();
            long j0_aggCell = p.i0;

            {
                int Nofields = this.m_DgDegree.Length;

                int[] Nfld = new int[Nofields];
                int[] Ofst = new int[Nofields];

                long i0 = this.Partitioning.i0;

                for(int jAgg = 0; jAgg < JAGG; jAgg++) {
                    
                    
                    for(int ifld = 0; ifld < Nofields; ifld++) {
                        int pField = this.m_DgDegree[ifld];
                        Nfld[ifld] = this.AggBasis[ifld].GetLength(jAgg, pField);
                        if(ifld > 0)
                            Ofst[ifld] = Ofst[ifld - 1] + Nfld[ifld - 1];
                    }
                    int Ncell = Ofst[Nofields - 1] + Nfld[Nofields - 1];

                    for(int i = 0; i < Fields.Length; i++) {
                        int iSpc = XaggBasis[i].GetSpeciesIndex(jAgg, Species);
                        if(iSpc >= 0) {

                            int NoSpc = XaggBasis[i].GetNoOfSpecies(jAgg);

                            int iField = Fields[i];
                            int pField = this.m_DgDegree[iField];

                            int N_iField = Nfld[iField];
                            int N0 = Ofst[iField];

                            int N_Spc = N_iField / NoSpc;
                            Debug.Assert(N_iField % NoSpc == 0);


                            for(int n = 0; n < N_Spc; n++) {
                                long iX = i0 + n + N0 + N_Spc * iSpc;
                                R.Add(iX);
                                Debug.Assert(this.Partitioning.IsInLocalRange(iX));
                                Debug.Assert(iX - j0_aggCell == this.LocalUniqueIndex(iField, jAgg, n + N_Spc * iSpc));
                            }
                        }
                    }
                    i0 += Ncell;
                }
            }

            return R.ToArray();
        }

        public int[] GetSubblk_i0(int blockType) {
            return m_Subblk_i0[blockType];
        }

        public int[] GetSubblkLen(int blockType) {
            return m_SubblkLen[blockType];
        }

        public int GetBlockType(long iBlock) {
            this.AggGrid.CellPartitioning.TestIfInLocalRange(iBlock);
            if (this.m_i0 == null) {
                Debug.Assert(MaximalLength == MinimalLength);
                return 0;
            } else {
                int iBlockLoc = this.AggGrid.CellPartitioning.TransformIndexToLocal(iBlock);
                int S = m_i0[iBlockLoc + 1] - m_i0[iBlockLoc];
                return m_Len2SublockType[S];
            }
        }

        public long GetBlockI0(long iBlock) {
            this.AggGrid.CellPartitioning.TestIfInLocalRange(iBlock);
            int iBlockLoc = this.AggGrid.CellPartitioning.TransformIndexToLocal(iBlock);
            if (m_i0 == null) {
                Debug.Assert(MaximalLength == MinimalLength);
                Debug.Assert(MaximalLength * iBlockLoc + this.i0 == iBlock * MaximalLength);
                return iBlock * MaximalLength;
            } else {
                return m_i0[iBlockLoc] + this.i0;
            }
        }

        /// <summary>
        /// Number of degrees-of-freedom, for all variables, in cell <paramref name="jCell"/>.
        /// </summary>
        /// <param name="jCell">Local cell index.</param>
        public int GetLength(int jCell) {
            int Nofields = this.m_DgDegree.Length;
            int S = 0;
            for (int ifld = 0; ifld < Nofields; ifld++) {
                int pField = this.m_DgDegree[ifld];
                S += this.AggBasis[ifld].GetLength(jCell, pField);
            }
            return S;
        }

        public int GetBlockLen(long iBlock) {
            this.AggGrid.CellPartitioning.TestIfInLocalRange(iBlock);
            int iBlockLoc = this.AggGrid.CellPartitioning.TransformIndexToLocal(iBlock);
            return this.GetLength(iBlockLoc);
        }

        public long GetBlockIndex(long i) {
            this.TestIfInLocalRange(i);
            if (MaximalLength == MinimalLength) {
                Debug.Assert(m_i0 == null);
                return i / MaximalLength;
            } else {
                int iLoc = this.TransformIndexToLocal(i);

                int iBlockLoc = Array.BinarySearch<int>(m_i0, iLoc);
                if (iBlockLoc < 0) {
                    iBlockLoc = (~iBlockLoc) - 1;
                }

                // Play it safe in case of potentially empty blocks
                Debug.Assert(this.m_i0.Length - 1 == this.LocalNoOfBlocks + this.AggGrid.iLogicalCells.NoOfExternalCells);
                while (m_i0[iBlockLoc] == m_i0[iBlockLoc + 1] && iBlockLoc < this.m_i0.Length - 2) {
                    iBlockLoc++;
                }

                Debug.Assert(iBlockLoc < this.LocalNoOfBlocks);
                Debug.Assert(iLoc >= m_i0[iBlockLoc]);
                Debug.Assert(iLoc < m_i0[iBlockLoc + 1]);
                return iBlockLoc + FirstBlock;
            }
        }

        public long GetFirstBlock(int proc) {
            return this.AggGrid.CellPartitioning.GetI0Offest(proc);
        }

        public int GetLocalNoOfBlocks(int proc) {
            return this.AggGrid.CellPartitioning.GetLocalLength(proc);
        }

        public int FindProcessForBlock(long iBlk) {
            return this.AggGrid.CellPartitioning.FindProcess(iBlk);
        }

        public IBlockPartitioning GetImmutableBlockPartitioning() {
            return this;
        }

        public long GetI0Offest(int proc) {
            return Partitioning.GetI0Offest(proc);
        }

        public IPartitioning GetImmutablePartition() {
            return this.Partitioning;
        }

        public int FindProcess(int index) {
            return this.Partitioning.FindProcess(index);
        }

        public int FindProcess(long index) {
            return this.Partitioning.FindProcess(index);
        }

        public bool IsInLocalRange(long i) {
            return this.Partitioning.IsInLocalRange(i);
        }

        public int GetLocalLength(int proc) {
            return this.Partitioning.GetLocalLength(proc);
        }

        public int Global2Local(long i) {
            return checked((int)(i - this.i0));
        }

        public bool IsXDGvariable(int iVar) {
            return AggBasis[iVar] is XdgAggregationBasis;
        }

        public int GetSpeciesIndex(int jCell, SpeciesId SId) {
            for(int iVar = 0; iVar < NoOfVariables; iVar++) {
                if(AggBasis[iVar] is XdgAggregationBasis xb)
                    return xb.GetSpeciesIndex(jCell, SId);
            }
            throw new NotSupportedException("Only DG variables; no species defined.");
        }

        public int GetNoOfSpecies(int jCell) {
            for(int iVar = 0; iVar < NoOfVariables; iVar++) {
                if(AggBasis[iVar] is XdgAggregationBasis xb)
                    return xb.GetNoOfSpecies(jCell);
            }
            return 1;
            //if(AggBasis[0] is XdgAggregationBasis xb)
            //    return xb.GetNoOfSpecies(jCell);
            //else
            //    return 1;
        }

        /// <summary>
        /// All used XDG species.
        /// Index: enumeration over species.
        /// </summary>
        public SpeciesId[] UsedSpecies {
            get {
                SpeciesId[] ret = null;
                for(int iVar = 0; iVar < NoOfVariables; iVar++) {
                    if(IsXDGvariable(iVar)) {
                        if(ret == null) {
                            ret = ((XdgAggregationBasis)AggBasis[iVar]).UsedSpecies.CloneAs();
                        } else {
                            if(!ret.SetEquals(((XdgAggregationBasis)AggBasis[iVar]).UsedSpecies))
                                throw new NotSupportedException("Different Variables seem to be defined for different species; not supported yet.");
                        }

                    }

                }
                return ret;
            }
        }
    }
}
