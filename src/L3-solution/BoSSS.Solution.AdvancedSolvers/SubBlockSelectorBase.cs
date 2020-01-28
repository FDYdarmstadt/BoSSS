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
        public SubBlockSelectorBase CellSelector(List<int> ListOfCellIdx, bool global = true) {
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
                tmpList = ListOfCellIdx;
            }
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
        /// 
        /// </summary>
        /// <param name="Species"></param>
        /// <returns></returns>
        public SubBlockSelectorBase SpeciesSelector(int Species) {
            // ist position von species in matrix consistent? Species A ist immer an erster Stelle, etc.?
            // Florian: NEIN
            // Diese Methode ist eher sinnlos - Selektion sollte mittels SpeciesId gemacht weren
            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                return GetIntInstruction(Species)(iSpec);
            };
            return this;
        }

        /// <summary>
        /// not implemented yet
        /// </summary>
        /// <param name="Species"></param>
        /// <returns></returns>
        public SubBlockSelectorBase SpeciesSelector(string Species) {
            int intSpecies = -1;
            throw new NotImplementedException("to do. Use other Selector instead");
            // how can be decided, which is species A or B?
            // Translate char to intSpecies and ...
            return SpeciesSelector(intSpecies); ;
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
    /// The block mask, which are generated through this class, presumes that the target Multigridoperator(defined through SubBlockSelector) is a quadratic Matrix.
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
        /// It also contains methods to translate Vectors from full matrix to sub matrix and vice versa.
        /// </summary>
        /// <param name="SBS"></param>
        public BlockMaskBase(SubBlockSelectorBase SBS, MPI_Comm MPIcomm) {
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
            var tmpCell = new List<extNi0[][][]>();

            // local caching of filter functions 
            // ensures that functions are not re-allocated during the loops
            var CellInstruction = m_sbs.CellFilter;
            var VarInstruction = m_sbs.VariableFilter;
            var SpecInstruction = m_sbs.SpeciesFilter;
            var ModeInstruction = m_sbs.ModeFilter;

            // loop over cells...
            for (int iLoc=0; iLoc < NoOfCells; iLoc++) {
                int jLoc = m_CellOffset + iLoc; //sic:to address external cells offset has to be concidered, you know ...
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
#endif
            Debug.Assert(ListNi0.GroupBy(x => x.Li0).Any(g => g.Count() == 1));
            Debug.Assert(ListNi0.GroupBy(x => x.Gi0).Any(g => g.Count() == 1));
            Debug.Assert(ListNi0.Count() == NumOfNi0);
            Debug.Assert(MaskLen <= m_map.LocalLength);
            Debug.Assert(Localint.GroupBy(x => x).Any(g => g.Count() == 1));
            Debug.Assert(Globalint.GroupBy(x => x).Any(g => g.Count() == 1));
            Debug.Assert(Localint.Count() == MaskLen);

            m_GlobalMask = Globalint;
            m_LocalMask = Localint;
            m_StructuredNi0 = tmpStructNi0;
            m_MaskLen = MaskLen;
            m_Ni0Len = Ni0Length;
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

        ///// <summary>
        ///// Get Array of Cellblocks
        ///// </summary>
        ///// <param name="target"></param>
        ///// <param name="ignoreCellCoupling"></param>
        ///// <param name="ignoreVarCoupling"></param>
        ///// <param name="ignoreSpecCoupling"></param>
        ///// <returns></returns>
        //public MultidimensionalArray[] GetSubBlocks(BlockMsrMatrix target, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {
        //    int NoOfCells = m_StructuredNi0.Length;
        //    int size = ignoreCellCoupling ? NoOfCells : NoOfCells * NoOfCells;

        //    MultidimensionalArray[] Sblocks = new MultidimensionalArray[size];

        //    int auxIdx = 0;
        //    for (int iLoc = 0; iLoc < m_StructuredNi0.Length; iLoc++) {
        //        for (int jLoc = 0; jLoc < m_StructuredNi0.Length; jLoc++) {
        //            if (ignoreCellCoupling && jLoc != iLoc) {
        //                continue;
        //            }
        //            int CellBlockLen = GetCellwiseLength(jLoc);
        //            Sblocks[auxIdx] = MultidimensionalArray.Create(CellBlockLen, CellBlockLen);
        //            for (int iVar = 0; iVar < m_StructuredNi0[iLoc].Length; iVar++) {
        //                for (int jVar = 0; jVar < m_StructuredNi0[jLoc].Length; jVar++) {
        //                    if (ignoreVarCoupling && jVar != iVar) {
        //                        continue;
        //                    }
        //                    for (int iSpc = 0; iSpc < m_StructuredNi0[iLoc][iVar].Length; iSpc++) {
        //                        for (int jSpc = 0; jSpc < m_StructuredNi0[jLoc][jVar].Length; jSpc++) {
        //                            if (ignoreSpecCoupling && jSpc != iSpc) {
        //                                continue;
        //                            }
        //                            for (int iMode = 0; iMode < m_StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
        //                                for (int jMode = 0; jMode < m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {
        //                                    extNi0 RowNi0 = m_StructuredNi0[iLoc][iVar][iSpc][iMode];
        //                                    extNi0 ColNi0 = m_StructuredNi0[jLoc][jVar][jSpc][jMode];
        //                                    int Targeti0 = RowNi0.Gi0;
        //                                    int Targetj0 = ColNi0.Gi0;
        //                                    int Subi0 = GetRelativeSubBlockOffset(iLoc, iVar, iSpc, iMode);
        //                                    int Subj0 = GetRelativeSubBlockOffset(jLoc, jVar, jSpc, jMode);
        //                                    int Subie = Subi0 + RowNi0.N - 1;
        //                                    int Subje = Subj0 + ColNi0.N - 1;

        //                                    target.ReadBlock(Targeti0, Targetj0,
        //                                        Sblocks[auxIdx].ExtractSubArrayShallow(new int[] { Subi0, Subj0 }, new int[] { Subie, Subje }));
        //                                }
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //            auxIdx++;
        //        }
        //    }
        //    return Sblocks;
        //}

        ///// <summary>
        ///// If you want nothing special. Take this one. If you want only diagonal block matrix choose one of the other methods
        ///// </summary>
        ///// <returns></returns>
        //public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix target) {
        //    int Loclength = m_MaskLen;

        //    var tmpN = GetAllSubMatrixCellLength();
        //    var tmpi0 = GetAllSubMatrixCellOffsets();

        //    BlockPartitioning localBlocking = new BlockPartitioning(Loclength, tmpi0.ToArray(), tmpN.ToArray(), m_MPIcomm, i0isLocal: true);
        //    var SubMSR = new BlockMsrMatrix(localBlocking);
        //    target.AccSubMatrixTo(1.0, SubMSR, m_GlobalMask, default(int[]), m_GlobalMask, default(int[]));
        //    return SubMSR;
        //}

        ///// <summary>
        ///// Coupling can be ignored. If you want full output take the basic GetSubBlockMatrix() method, which will be faster.
        ///// </summary>
        ///// <param name="target"></param>
        ///// <param name="ignoreCellCoupling"></param>
        ///// <param name="ignoreVarCoupling"></param>
        ///// <param name="ignoreSpecCoupling"></param>
        ///// <returns></returns>
        //public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix target, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {
        //    int Loclength = m_MaskLen;

        //    var tmpN = GetAllSubMatrixCellLength();
        //    var tmpi0 = GetAllSubMatrixCellOffsets();

        //    BlockPartitioning localBlocking = new BlockPartitioning(Loclength, tmpi0.ToArray(), tmpN.ToArray(), m_MPIcomm, i0isLocal: true);
        //    var SubMSR = new BlockMsrMatrix(localBlocking);

        //    int SubRowIdx = 0;
        //    int SubColIdx = 0;

        //    int auxIdx = 0;
        //    for (int iLoc = 0; iLoc < m_StructuredNi0.Length; iLoc++) {
        //        for (int jLoc = 0; jLoc < m_StructuredNi0.Length; jLoc++) {
        //            if (ignoreCellCoupling && jLoc != iLoc) {
        //                continue;
        //            }
        //            for (int iVar = 0; iVar < m_StructuredNi0[iLoc].Length; iVar++) {
        //                for (int jVar = 0; jVar < m_StructuredNi0[jLoc].Length; jVar++) {
        //                    if (ignoreVarCoupling && jVar != iVar) {
        //                        continue;
        //                    }
        //                    for (int iSpc = 0; iSpc < m_StructuredNi0[iLoc][iVar].Length; iSpc++) {
        //                        for (int jSpc = 0; jSpc < m_StructuredNi0[jLoc][jVar].Length; jSpc++) {
        //                            if (ignoreSpecCoupling && jSpc != iSpc) {
        //                                continue;
        //                            }
        //                            for (int iMode = 0; iMode < m_StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
        //                                for (int jMode = 0; jMode < m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {

        //                                    int[] RowIdx = new int[] { iLoc, iVar, iSpc, iMode };
        //                                    int[] ColIdx = new int[] { jLoc, jVar, jSpc, jMode };
        //                                    extNi0 RowNi0 = m_StructuredNi0[iLoc][iVar][iSpc][iMode];
        //                                    extNi0 ColNi0 = m_StructuredNi0[jLoc][jVar][jSpc][jMode];
        //                                    int Targeti0 = RowNi0.Gi0;
        //                                    int Targetj0 = ColNi0.Gi0;

        //                                    var tmpBlock = MultidimensionalArray.Create(RowNi0.N, ColNi0.N);
        //                                    target.ReadBlock(Targeti0, Targetj0,
        //                                        tmpBlock);
        //                                    SubMSR.AccBlock(SubRowIdx, SubColIdx, 1, tmpBlock);
        //                                    SubColIdx += ColNi0.N;
        //                                }
        //                                SubRowIdx += m_StructuredNi0[iLoc][iVar][iSpc][iMode].N;
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //            auxIdx++;
        //        }
        //        auxIdx++;
        //    }
        //    return SubMSR;
        //}

        ///// <summary>
        ///// Translates back your masked vector to the full one
        ///// </summary>
        ///// <typeparam name="V"></typeparam>
        ///// <typeparam name="W"></typeparam>
        ///// <param name="targetVector">output</param>
        ///// <param name="accVector"></param>
        //public void AccVecToFull<V, W>(W accVector, V targetVector)
        //    where V : IList<double>
        //    where W : IList<double> {
        //    if (targetVector.Count() != m_map.LocalLength)
        //        throw new ArgumentException("Length of targetVector != Length of original");
        //    if (accVector.Count() != m_MaskLen)
        //        throw new ArgumentException("accVector length is not equal to length of mask");
        //    targetVector.AccV(1.0, accVector, m_LocalMask, default(int[]));
        //}

        ///// <summary>
        ///// Translates back your masked vector to the full one cellwise. i.e. Useful, if you solved blockwise
        ///// </summary>
        ///// <typeparam name="V"></typeparam>
        ///// <typeparam name="W"></typeparam>
        ///// <param name="targetVector">output</param>
        ///// <param name="accVector"></param>
        ///// <param name="iCell">starts with 0</param>
        //public void AccVecCellwiseToFull<V, W>(W accVector, int iCell, V targetVector)
        //    where V : IList<double>
        //    where W : IList<double> {
        //    if (iCell > m_StructuredNi0.Length - 1)
        //        throw new ArgumentOutOfRangeException("iCell is greater than Cells in mask");
        //    if (accVector.Count() != GetCellwiseLength(iCell))
        //        throw new ArgumentException("accVector length is not equal to length of mask");
        //    if (targetVector.Count() != m_map.LocalLength)
        //        throw new ArgumentException("Length of targetVector not equal Length of original");
        //    var Cidx = GetCellwiseLocalidx(iCell);
        //    Debug.Assert(accVector.Count() == Cidx.Count());
        //    targetVector.AccV(1.0, accVector, Cidx, default(int[]));
        //}

        ///// <summary>
        ///// Translates a full vector to a vector corresponding to matrix mask cellwise
        ///// </summary>
        ///// <param name="fullVector"></param>
        ///// <param name="iCell"></param>
        ///// <returns></returns>
        //public double[] GetVectorCellwise(IList<double> fullVector, int iCell) {
        //    if (iCell > m_StructuredNi0.Length - 1)
        //        throw new ArgumentOutOfRangeException("iCell is greater than number of cells in mask");
        //    if (fullVector.Count() != m_map.LocalLength)
        //        throw new ArgumentException("Length of targetVector not equal Length of original");
        //    double[] subVector = new double[GetCellwiseLength(iCell)];
        //    var Cidx = GetCellwiseLocalidx(iCell);
        //    Debug.Assert(subVector.Length == Cidx.Length);
        //    ArrayTools.GetSubVector<int[], int[], double>(fullVector, subVector, Cidx);
        //    return subVector;
        //}

        ///// <summary>
        ///// Translates a full vector to a vector corresponding to matrix mask
        ///// </summary>
        ///// <param name="fullVector"></param>
        //public double[] GetSubBlockVec(IList<double> fullVector) {
        //    if (fullVector.Count() != m_map.LocalLength)
        //        throw new ArgumentException("Length of targetVector not equal Length of original");
        //    int Len = m_MaskLen;
        //    double[] subVector = new double[m_MaskLen];
        //    subVector.AccV(1.0, fullVector, default(int[]), m_LocalMask);
        //    return subVector;
        //}

    }

}
