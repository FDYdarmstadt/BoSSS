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
using MathNet.Numerics.Algorithms.LinearAlgebra;
using System.Numerics;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Comm;
using System.Collections;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// This contains a block Selection, which can be specified by the user with the hierachical Selectors: every Matrix is subdivided into cell, variable, species and mode blocks.
    /// </summary>
    public abstract class SubBlockSelectorBase {

        /// <summary>
        /// Specifies, which blocks in a matrix shall be selected. Blocksubdivision Default: Selects all blocks.  
        /// </summary>
        /// <param name="map"></param>
        public SubBlockSelectorBase(MultigridMapping map) {
            if (map == null)
                throw new Exception("empty mapping! This will not end well ...");
            m_map = map;
            this.CellSelector();
            this.VariableSelector();
            this.SpeciesSelector();
            this.ModeSelector();
        }

        //internal get
        protected MultigridMapping m_map;

        //internal set
        protected Func<int, bool> m_CellFilter = null;
        protected Func<int, int, bool> m_VariableFilter = null;
        protected Func<int, int, int, bool> m_SpeciesFilter = null;
        protected Func<int, int, int, int, bool> m_ModeFilter = null;

        #region CellSelector
        /// <summary>
        /// Selects all aggregation cell blocks
        /// </summary>
        /// <returns></returns>
        public SubBlockSelectorBase CellSelector() {
            int NoOfCells = m_NoLocalCells;
            
            this.m_CellFilter = GetAllInstruction(NoOfCells);
            return this;
        }

        /// <summary>
        /// Selects cell according to global/local cell index
        /// </summary>
        /// <param name="CellIdx"></param>
        /// <param name="global"></param>
        /// <returns></returns>
        public SubBlockSelectorBase CellSelector(int CellIdx, bool global = true) {

            int LocNoOfBlocks = m_NoLocalCells;
            int GlobNoOfBlocks = m_NoTotalCells;
            if (global) {
                if (CellIdx >= GlobNoOfBlocks || CellIdx < 0)
                    throw new ArgumentOutOfRangeException(CellIdx + " is greater then global No of Blocks: " + GlobNoOfBlocks + " or smaller than 0");
            } else {
                if (CellIdx >= LocNoOfBlocks || CellIdx < 0)
                    throw new ArgumentOutOfRangeException(CellIdx + " is greater then Local No of Blocks: " + LocNoOfBlocks + " or smaller than 0");
            }

            if (global) {
                //translate into local index ...
                CellIdx -= m_i0;

                int idxfound = 0;
                if (m_i0 <= CellIdx && m_iE >= CellIdx)
                    idxfound = 1;
                Debug.Assert(idxfound.MPISum() == 1);

                if (m_i0 > CellIdx || m_iE < CellIdx) {
                    this.m_CellFilter = GetDoNothingInstruction();
                    return this;
                }
            }

            this.m_CellFilter = GetIntInstruction(CellIdx);
            return this;
        }

        /// <summary>
        /// Selects List of cells according to global/local index
        /// </summary>
        /// <param name="ListOfCellIdx"></param>
        /// <param name="global"></param>
        /// <returns></returns>
        public SubBlockSelectorBase CellSelector<V>(V ListOfCellIdx, bool global = true)
            where V : IList<int>{
            int LocNoOfBlocks = m_NoLocalCells;
            int GlobNoOfBlocks = m_NoTotalCells;

            foreach (int CellIdx in ListOfCellIdx) {
                if (global) {
                    if (CellIdx >= GlobNoOfBlocks || CellIdx < 0)
                        throw new ArgumentOutOfRangeException(CellIdx + " is greater then global No of Blocks: " + GlobNoOfBlocks + " or smaller than 0");
                } else {
                    if (CellIdx >= LocNoOfBlocks || CellIdx < 0)
                        throw new ArgumentOutOfRangeException(CellIdx + " is greater then Local No of Blocks: " + LocNoOfBlocks + " or smaller than 0");
                }
            }


            List<int> tmpList = new List<int>();
            if (global) {
                foreach (int CellIdx in ListOfCellIdx) {
                    int tmpIdx = CellIdx;
                    tmpIdx -= m_i0;

                    int idxfound = 0;
                    if (m_i0 <= tmpIdx && m_iE >= tmpIdx)
                        idxfound = 1;
                    Debug.Assert(idxfound.MPISum() == 1);

                    if (idxfound == 1)
                        tmpList.Add(tmpIdx);
                }
            } else {
                tmpList = ListOfCellIdx.ToList();
            }
            if (tmpList.Count <= 0)
                this.m_CellFilter = GetDoNothingInstruction();
            if(tmpList.Count>0)
                this.m_CellFilter = GetListInstruction(tmpList);
            return this;
        }

        /// <summary>
        /// Selects all cells in a cell mask.
        /// </summary>
        public SubBlockSelectorBase CellSelector(CellMask CM) {
            if (CM == null) {
                m_CellFilter = jCell => true;
            } else {
                var BitMask = CM.GetBitMask();
                this.m_CellFilter = (int jCell) => BitMask[jCell];
            }
            return this;
        }
        #endregion

        #region VariableSelector
        
        /// <summary>
        /// Selects all Variable blocks
        /// </summary>
        /// <returns></returns>
        public SubBlockSelectorBase VariableSelector() {
            this.m_VariableFilter = delegate (int iCell, int iVar) {
#if DEBUG
                int NoOfVar = m_AggBS.Length;
                Debug.Assert(iVar >= 0);
                Debug.Assert(iVar < NoOfVar);
#endif
                return true;
            };
            return this;
        }

        /// <summary>
        /// 
        /// </summary>
        public SubBlockSelectorBase VariableSelector(params int[] SetOfVariables) {
            return VariableSelector((ICollection<int>)SetOfVariables);
        }

        /// <summary>
        /// 
        /// </summary>
        public SubBlockSelectorBase VariableSelector(int iVariable) {
            if (iVariable < 0)
                throw new ArgumentOutOfRangeException("Variable index cannot be negative.");
            if (iVariable >= m_NoOfVar)
                throw new ArgumentOutOfRangeException("Variable index is larger than number of variables..");

            return VariableSelector(new int[] { iVariable });
        }


        /// <summary>
        /// 
        /// </summary>
        public SubBlockSelectorBase VariableSelector(IEnumerable<int> SetOfVariables) {
            int[] Variables = SetOfVariables.ToArray();
            if (!Variables.IsSet())
                throw new ArgumentOutOfRangeException("Input is not a set - some variable index is listed twice.");
            if (Variables.Min() < 0)
                throw new ArgumentOutOfRangeException("Variable index cannot be negative.");
            if (Variables.Max() >= m_NoOfVar)
                throw new ArgumentOutOfRangeException("Some variable index is larger than number of variables..");

            int NoOfVar = m_AggBS.Length;

            var VarSelector = GetListInstruction(Variables);

            this.m_VariableFilter = delegate (int iCell, int iVar) {
                // Kommentar an Jens: deine orginale Impl. war nicht ideal: bei jedem Aufruf wird 'GetInstruction neu aufgerufen'
                //return GetListInstruction(ListOfVariables)(iVar); // nach ausserhalb verschoben
                return VarSelector(iVar);
            };
            return this;
        }
#endregion

#region SpeciesSelector
        /// <summary>
        /// Selects all species blocks
        /// </summary>
        /// <returns></returns>
        public SubBlockSelectorBase SpeciesSelector() {
            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                int NoOfSpec = m_AggBS[iVar].GetNoOfSpecies(iCell);
                return GetAllInstruction(NoOfSpec)(iSpec);
            };
            return this;
        }

        /// <summary>
        /// Selects Species by <see cref="SpeciesId"/>.
        /// </summary>
        /// <param name="Species"></param>
        /// <returns></returns>
        public SubBlockSelectorBase SpeciesSelector(SpeciesId SId) {
            //for (int v = 0; v < m_map.NoOfVariables; v++) {
            //    if (this.m_map.AggBasis[v].GetType() == typeof(XdgAggregationBasis))
            //        Console.WriteLine("WARNING: Variable {0} has no XdgBasis and thus may be not selected",v);
            //}
            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                int SpcIdx = ((XdgAggregationBasis)this.m_map.AggBasis[iVar]).GetSpeciesIndex(iCell, SId);
                return GetIntInstruction(SpcIdx)(iSpec);
            };
            return this;
        }

        /// <summary>
        /// Selects multiple species by <see cref="IEnumerable<SpeciesId>"/>.
        /// </summary>
        /// <param name="SetOfSpecies"></param>
        /// <returns></returns>
        public SubBlockSelectorBase SpeciesSelector(IEnumerable<SpeciesId> SetOfSpecies) {
            //for (int v = 0; v < m_map.NoOfVariables; v++) {
            //    if (this.m_map.AggBasis[v].GetType() == typeof(XdgAggregationBasis))
            //        Console.WriteLine("WARNING: Variable {0} has no XdgBasis and thus may be not selected",v);
            //}
            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                List<int> SpcInt = new List<int>();
                foreach(SpeciesId spi in SetOfSpecies) {
                    int SpcIdx = ((XdgAggregationBasis)this.m_map.AggBasis[iVar]).GetSpeciesIndex(iCell, spi);
                    SpcInt.Add(SpcIdx);
                }
                return GetListInstruction(SpcInt)(iSpec);
            };
            return this;
        }

#endregion

#region ModeSelector

        /// <summary>
        /// Selects all Modes
        /// </summary>
        /// <returns></returns>
        public SubBlockSelectorBase ModeSelector() {
            this.m_ModeFilter = delegate (int iCell, int iVar, int iSpec, int iMode) {
#if DEBUG
                int maxDG = m_DGdegree[iVar];
                Debug.Assert(iMode <= maxDG);
#endif
                return true;
            };
            return this;
        }

        /// <summary>
        /// Selects Modes according to instruction.
        /// </summary>
        public SubBlockSelectorBase ModeSelector(Func<int, bool> boolinstruct) {
            this.m_ModeFilter = delegate (int iCell, int iVar, int iSpec, int pDeg) {
                return boolinstruct(pDeg);
            };
            return this;
        }

        /// <summary>
        /// Selects Modes according to instruction.
        /// </summary>
        public SubBlockSelectorBase ModeSelector(Func<int, int, int, int, bool> boolinstruct) {
            this.m_ModeFilter = boolinstruct;
            return this;
        }


        #endregion

        #region BasisInstructions

        private Func<int, bool> GetListInstruction(IEnumerable<int> intList) {
            if(!intList.IsSet())
                throw new ArgumentException("List contains duplicates!");
            bool[] Selector = new bool[intList.Max() + 1];
            foreach(int i in intList) {
                Selector[i] = true;
            }

            return delegate (int Index) {
                if(Index >= Selector.Length) {
                    return false;
                } else {
                    return Selector[Index];
                }
            };
        }

        private Func<int, bool> GetIntInstruction(int thisint) {
            return delegate (int Index) {
                return thisint == Index;
            };
        }

        
        private Func<int, bool> GetAllInstruction(int MaxIdx) {
            //delegate, which always return true is also possible
            //List<int> tmp = new List<int>();
            //for (int iLoc = 0; iLoc < MaxIdx; iLoc++) {
            //    tmp.Add(iLoc);
            //}
            //return GetListInstruction(tmp);
            return delegate (int i) {
                return i < MaxIdx;
            };
        }

        private Func<int, bool> GetDoNothingInstruction() {
            return delegate (int idx) {
                return false;
            };
        }

        #endregion

        #region Getta

        /// <summary>
        /// returns a delegate \f$ (j, \gamma, \mathfrak{s}, k ) \mapsto \{ \textrm{true},\textrm{false} \} \f$,
        /// which determines whether a specific cell-variable-species-mode tuple is in a mapping or not.
        /// Delegate arguments are:
        /// - 1st argument: cell index \f$ j \f$
        /// - 2nd argument: variable index \f$ \gamma \f$
        /// - 3rd argument: species index \f$ \mathfrak{s} \f$
        /// - 4th argument: DG polynomial degree \f$ k \f$
        /// </summary>
        public Func<int, int, int, int, bool> CombinedFilter {
            get {
                return delegate (int iCell, int iVar, int iSpec, int iMod) {
                    bool result = false;
                    result &= m_CellFilter(iCell);
                    result &= m_VariableFilter(iCell, iVar);
                    result &= m_SpeciesFilter(iCell, iVar, iSpec);
                    result &= m_ModeFilter(iCell, iVar, iSpec, iMod);
                    return result;
                };
            }
        }

        public Func<int, bool> CellFilter {
            get { return m_CellFilter; }
        }

        public Func<int, int, bool> VariableFilter {
            get { return m_VariableFilter; }
        }

        public Func<int, int, int, bool> SpeciesFilter {
            get { return m_SpeciesFilter; }
        }

        public Func<int, int, int, int, bool> ModeFilter {
            get { return m_ModeFilter; }
        }


        /// <summary>
        /// gets the multigrid operator on which this selector shall work on 
        /// </summary>
        public MultigridMapping GetMapping {
            get {
                return m_map;
            }
        }

        private AggregationGridBasis[] m_AggBS {
            get {
                return m_map.AggBasis;
            }
        }

        private int[] m_DGdegree {
            get {
                return m_map.DgDegree;
            }
        }

        private int m_NoLocalCells {
            get {
                return m_map.LocalNoOfBlocks + m_map.AggBasis[0].AggGrid.iLogicalCells.NoOfExternalCells;
            }
        }

        private int m_NoTotalCells {
            get {
                return m_map.TotalNoOfBlocks;
            }
        }
        private int m_NoOfVar {
            get {
                return m_map.NoOfVariables;
            }
        }

        private int m_i0 {
            get {
                return m_map.i0;
            }
        }

        private int m_iE {
            get {
                return m_map.iE;
            }
        }

        

        #endregion
    }


    public struct extNi0
    {
        public extNi0(int Li0, int Gi0, int Si0, int N) {
            m_Li0 = Li0;
            m_Gi0 = Gi0;
            m_Si0 = Si0;
            m_N = N;
        }
        private int m_Li0, m_Gi0, m_Si0, m_N;
        public int N {
            get {
                return m_N;
            }
        }
        public int Li0 {
            get {
                return m_Li0;
            }
        }
        public int Gi0 {
            get {
                return m_Gi0;
            }
        }
        public int Si0 {
            get {
                return m_Si0;
            }
        }
    }

    /// <summary>
    /// This abstract class is the unification of internal and external cell masking.
    /// Therefore <see cref="m_NoOfCells"/> and <see cref="m_CellOffset"/> have to be overriden by inheriting classes.
    /// There are two inheriting classes: <see cref="BlockMaskExt"/> and <see cref="BlockMaskLoc"/>,
    /// which handle the masking of external cells and internal cells respectively.
    /// The purpose of this class is to provide the index lists and structs,
    /// which correspond to underlying selection through the <see cref="SubBlockSelector"/>.
    /// </summary>
    public abstract class BlockMaskBase {

        /// <summary>
        /// auxiliary structure. Stores Offsets and Length of something.
        /// </summary>
        
        private struct Ni0 {
            public Ni0(int i0, int N) {
                m_i0 = i0;
                m_N = N;
            }
            private int m_i0, m_N;
            public int N {
                get {
                    return m_N;
                }
            }
            public int i0 {
                get {
                    return m_i0;
                }
            }
        }

        /// <summary>
        /// Generates Block Mask from Sub block selection. 
        /// Sub matrix can be generated of the multigrid operator overgiven to Sub block selection.
        /// </summary>
        /// <param name="SBS"></param>
        public BlockMaskBase(SubBlockSelector SBS, MPI_Comm MPIcomm) {
            m_map = SBS.GetMapping;
            m_sbs = SBS;
            m_AggBS = m_map.AggBasis;
            m_DGdegree = m_map.DgDegree;
            m_Ni0 = Ni0Gen();
            m_NoOfVariables = m_AggBS.Length;
            //Testen ob es cells gibt wo Var<>NoOfVar, das würde dementsprechend auch m_DG beeinflussen
            m_NoOfSpecies = new int[m_NoOfCells][];
            for (int iCell = 0; iCell < m_NoOfCells; iCell++) {
                m_NoOfSpecies[iCell] = new int[m_NoOfVariables];
                for (int iVar = 0; iVar < m_NoOfVariables; iVar++) {
                    m_NoOfSpecies[iCell][iVar] = m_AggBS[iVar].GetNoOfSpecies(iCell);
                }
            }
            //TestForQuadraticMatrix();   
            GenerateAllMasks();
        }

        protected abstract int m_NoOfCells {
            get;
        }

        protected abstract int m_CellOffset {
            get;
        }

        protected abstract int m_LocalLength {
            get;
        }

        //internal get
        private SubBlockSelectorBase m_sbs;
        protected MultigridMapping m_map;
        private AggregationGridBasis[] m_AggBS;
        private int[] m_DGdegree;
        private int m_NoOfVariables;
        private int[][] m_NoOfSpecies;
        private MPI_Comm m_MPIcomm;
        //internal set
        public List<int> m_GlobalMask = null;
        public List<int> m_LocalMask = null;
        public List<int> m_SubBlockMask = null;
        public extNi0[][][][] m_StructuredNi0 = null;
        public int m_Ni0Len;
        private int m_MaskLen;
        private Ni0[] m_Ni0;

        private bool TestForQuadraticMatrix() {
            throw new NotSupportedException("Supports only quadratic matrices");
        }

        private void GenerateAllMasks() {
            int NoOfCells = m_NoOfCells;
            int NoOfVariables = m_NoOfVariables;
            int[][] NoOfSpecies = m_NoOfSpecies;
            int[] DGdegreeP1 = m_DGdegree.CloneAs();
            for (int iDG = 0; iDG < DGdegreeP1.Length; iDG++) {
                DGdegreeP1[iDG] += 1;
            }

            List<extNi0> ListNi0 = new List<extNi0>();
            List<int> Globalint = new List<int>();
            List<int> Localint = new List<int>();
            List<int> SubBlockIdx = new List<int>();

            int SubOffset = 0;
            int Ni0Length = 0;
            int MaskLen = 0;
            int prevLocie = m_map.LocalNoOfBlocks; 
            var tmpCell = new List<extNi0[][][]>();

            // local caching of filter functions 
            // ensures that functions are not re-allocated during the loops
            var CellInstruction = m_sbs.CellFilter;
            var VarInstruction = m_sbs.VariableFilter;
            var SpecInstruction = m_sbs.SpeciesFilter;
            var ModeInstruction = m_sbs.ModeFilter;

            // loop over cells...
            for (int iLoc=0; iLoc < NoOfCells; iLoc++) {
                int jLoc = m_CellOffset + iLoc; //sic:to address correctly, external cells offset has to be concidered, you know ...
                if (!m_sbs.CellFilter(jLoc))
                    continue;
                var tmpVar = new List<extNi0[][]>();

                // loop over variables...
                for (int iVar = 0; iVar < NoOfVariables; iVar++) {
                    if (!VarInstruction(jLoc, iVar))
                        continue;
                    var tmpSpc = new List<extNi0[]>();

                    // loop over species...
                    for (int iSpc = 0; iSpc < NoOfSpecies[iLoc][iVar]; iSpc++) {
                        if (!SpecInstruction(jLoc, iVar, iSpc))
                            continue;
                        int GlobalOffset = m_map.GlobalUniqueIndex(iVar, jLoc, iSpc, 0);
                        int LocalOffset = m_map.LocalUniqueIndex(iVar, jLoc, iSpc, 0);
                        var tmpMod = new List<extNi0>();

                        // loop over polynomial degrees...
                        for (int degree = 0; degree < DGdegreeP1[iVar]; degree++) {
                            if (ModeInstruction(jLoc, iVar, iSpc, degree)) {
                                int GlobalModeOffset = m_Ni0[degree].i0 + GlobalOffset;
                                int LocalModeOffset = m_Ni0[degree].i0 + LocalOffset;
                                int ModeLength = m_Ni0[degree].N;
                                var newNi0 = new extNi0(LocalModeOffset, GlobalModeOffset, SubOffset, ModeLength);
                                    SubOffset += ModeLength;
                                    // Fill int lists
                                    for (int i = 0; i < newNi0.N; i++) {
                                        Globalint.Add(newNi0.Gi0 + i);
                                        Localint.Add(newNi0.Li0 + i);
                                        SubBlockIdx.Add(newNi0.Si0 + i);
                                        MaskLen++;
                                    }
                                    // Fill Ni0 Lists
                                    tmpMod.Add(newNi0);
                                    Ni0Length++;
                                    ListNi0.Add(newNi0);
                                Debug.Assert(m_map.LocalUniqueIndex(iVar, jLoc, iSpc, GetNp(degree) - 1) == LocalModeOffset + ModeLength - 1);
                            }
                        }
                        tmpSpc.Add(tmpMod.ToArray());
                    }
                    tmpVar.Add(tmpSpc.ToArray());
                }
                    tmpCell.Add(tmpVar.ToArray());
            }
            var tmpStructNi0 = tmpCell.ToArray();

            int NumOfNi0 = 0;
#if DEBUG
            for (int iCell=0; iCell < tmpStructNi0.Length; iCell++){
                for(int iVar=0; iVar < tmpStructNi0[iCell].Length; iVar++){
                    for(int iSpc=0; iSpc < tmpStructNi0[iCell][iVar].Length; iSpc++){
                        NumOfNi0+=tmpStructNi0[iCell][iVar][iSpc].Length;
                    }
                }
            }

            bool emptysel = true;
            for(int iLoc = 0; iLoc < NoOfCells; iLoc++) {
                emptysel &= !m_sbs.CellFilter(iLoc);
            }
#endif
            if (emptysel)
                Console.WriteLine("WARNING: no cells selceted with {0}", this.ToString());
            if (!emptysel) {
                Debug.Assert(ListNi0.GroupBy(x => x.Li0).Any(g => g.Count() == 1));
                Debug.Assert(ListNi0.GroupBy(x => x.Gi0).Any(g => g.Count() == 1));
                Debug.Assert(ListNi0.Count() == NumOfNi0);
                Debug.Assert(MaskLen <= m_LocalLength);
                Debug.Assert(Localint.GroupBy(x => x).Any(g => g.Count() == 1));
                Debug.Assert(Globalint.GroupBy(x => x).Any(g => g.Count() == 1));
                Debug.Assert(Localint.Count() == MaskLen);
                Debug.Assert(SubBlockIdx.Count() == MaskLen);
                Debug.Assert(SubBlockIdx.GroupBy(x => x).Any(g => g.Count() == 1));
            }

            m_GlobalMask = Globalint;
            m_LocalMask = Localint;
            m_StructuredNi0 = tmpStructNi0;
            m_MaskLen = MaskLen;
            m_Ni0Len = Ni0Length;
            m_SubBlockMask = SubBlockIdx;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public int[] GetCellwiseLocalidx(int iCell) {
            List<int> cellidx = new List<int>();
            for (int i = 0; i < m_StructuredNi0[iCell].Length; i++) {
                for (int j = 0; j < m_StructuredNi0[iCell][i].Length; j++) {
                    for (int k = 0; k < m_StructuredNi0[iCell][i][j].Length; k++) {
                        extNi0 block = m_StructuredNi0[iCell][i][j][k];
                        for (int m = 0; m < block.N; m++) {
                            cellidx.Add(block.Li0 + m);
                        }
                    }
                }
            }
            int[] array = cellidx.ToArray();
            Debug.Assert(array.GroupBy(x => x).Any(g => g.Count() == 1));
            return array;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public int[] GetCellwiseGlobalidx(int iCell) {
            List<int> cellidx = new List<int>();
            for (int i = 0; i < m_StructuredNi0[iCell].Length; i++) {
                for (int j = 0; j < m_StructuredNi0[iCell][i].Length; j++) {
                    for (int k = 0; k < m_StructuredNi0[iCell][i][j].Length; k++) {
                        extNi0 block = m_StructuredNi0[iCell][i][j][k];
                        for (int m = 0; m < block.N; m++) {
                            cellidx.Add(block.Gi0 + m);
                        }
                    }
                }
            }
            int[] array = cellidx.ToArray();
            Debug.Assert(array.GroupBy(x => x).Any(g => g.Count() == 1));
            return array;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public int GetCellwiseLength(int iCell) {
            int len = 0;
            for (int i = 0; i < m_StructuredNi0[iCell].Length; i++) {
                for (int j = 0; j < m_StructuredNi0[iCell][i].Length; j++) {
                    for (int k = 0; k < m_StructuredNi0[iCell][i][j].Length; k++) {
                        extNi0 block = m_StructuredNi0[iCell][i][j][k];
                        len += block.N;
                    }
                }
            }
            return len;
        }

        public List<int> GetAllSubMatrixCellLength() {
            var intList = new List<int>();
            for (int iCell = 0; iCell < m_StructuredNi0.Length; iCell++) {
                intList.Add(GetCellwiseLength(iCell));
            }
            return intList;
        }

        private int GetCellwiseSubIdx(int iCell) {
            return m_StructuredNi0[iCell][0][0][0].Si0;
        }

        public List<int> GetAllSubMatrixCellOffsets() {
            var intList = new List<int>();
            for (int iCell = 0; iCell < m_StructuredNi0.Length; iCell++)
                intList.Add(GetCellwiseSubIdx(iCell));
            return intList;
        }

        /// <summary>
        /// gets subblock offset relative to the parent cellblock
        /// </summary>
        /// <param name="iCell"></param>
        /// <param name="iVar"></param>
        /// <param name="iSpc"></param>
        /// <param name="iMod"></param>
        /// <returns></returns>
        public int GetRelativeSubBlockOffset(int iCell, int iVar, int iSpc, int iMod) {
            return m_StructuredNi0[iCell][iVar][iSpc][iMod].Si0 - m_StructuredNi0[iCell][0][0][0].Si0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        private extNi0[][][] GetSubBlockNi0Cellwise(int iCell) {
            return m_StructuredNi0[iCell];
        }

        /// <summary>
        /// Pre generates Offsets and Length for a default Variable block with the maximal DGdegree, to save computation time.
        /// </summary>
        /// <returns></returns>
        private Ni0[] Ni0Gen() {
            int maxDG = m_DGdegree.Max();
            var ModeOffsetNLengt = new List<Ni0>();
            for (int p = 0; p < maxDG + 1; p++) {
                int aux = p == 0 ? 0 : GetNp(p - 1);
                int pOffset = aux;
                int pLength = GetNp(p) - aux;
                Debug.Assert(p == 0 || ModeOffsetNLengt.Last().i0 < pOffset);
                Debug.Assert(p == 0 || ModeOffsetNLengt.Last().N < pLength);
                ModeOffsetNLengt.Add(new Ni0(pOffset, pLength));
            }
            Ni0[] Ni0array = ModeOffsetNLengt.ToArray();
            Debug.Assert(Ni0array.GroupBy(x => x.i0).Any(g => g.Count() == 1));
            return Ni0array;
        }

        private int GetNp(int p) {
            int Np = -1;
            int SpacDim = m_map.AggGrid.SpatialDimension;
            Debug.Assert(p >= 0);
            switch (SpacDim) {
                case 1:
                    Np = p + 1 + p + 1;
                    break;
                case 2:
                    Np = (p * p + 3 * p + 2) / 2;
                    break;
                case 3:
                    Np = (p * p * p + 6 * p * p + 11 * p + 6) / 6;
                    break;
                default:
                    throw new Exception("wtf?Spacialdim=1,2,3 expected");
            }
            return Np;
        }

        /// <summary>
        /// Returns the local number of subblocks rowwise.
        /// </summary>
        public int LocalLength {
            get {
                return m_Ni0Len;
            }
        }

        /// <summary>
        /// Returns local number of DOF of Sub Matrix
        /// </summary>
        public int LocalDOF {
            get {
                return m_MaskLen;
            }
        }

    }

}
