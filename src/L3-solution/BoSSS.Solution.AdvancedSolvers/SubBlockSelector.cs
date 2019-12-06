﻿/* =======================================================================
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
    public class SubBlockSelector {

        /// <summary>
        /// Specifies, which blocks in a matrix shall be selected. Blocksubdivision Default: Selects all blocks.  
        /// </summary>
        /// <param name="map"></param>
        public SubBlockSelector() {
            this.CellSelector();
            this.VariableSelector();
            this.SpeciesSelector();
            this.ModeSelector();
        }

        //internal get
        private MultigridMapping m_map;
        //internal set
        private Func<int, bool> m_CellSelInstruction = null;
        private Func<int, int, bool> m_VariableSelInstruction = null;
        private Func<int, int, int, bool> m_SpeciesSelInstruction = null;
        private Func<int, int, int, int, bool> m_ModeSelInstruction = null;

        #region CellSelector
        /// <summary>
        /// Selects all aggregation cell blocks
        /// </summary>
        /// <returns></returns>
        public SubBlockSelector CellSelector() {
            int NoOfCells = m_map.LocalNoOfBlocks;
            this.m_CellSelInstruction = GetAllInstruction(NoOfCells);
            return this;
        }

        /// <summary>
        /// Selects cells according to global/local cell index
        /// </summary>
        /// <param name="CellIdx"></param>
        /// <param name="global"></param>
        /// <returns></returns>
        public SubBlockSelector CellSelector(int CellIdx, bool global = true) {
            int LocNoOfBlocks = m_map.LocalNoOfBlocks;
            int GlobNoOfBlocks = m_map.TotalNoOfBlocks;
            if (global) {
                if (CellIdx >= GlobNoOfBlocks || CellIdx < 0)
                    throw new ArgumentOutOfRangeException(CellIdx + " is greater then global No of Blocks: " + GlobNoOfBlocks + " or smaller than 0");
            } else {
                if (CellIdx >= LocNoOfBlocks || CellIdx < 0)
                throw new ArgumentOutOfRangeException(CellIdx + " is greater then Local No of Blocks: " + LocNoOfBlocks + " or smaller than 0");
            }

            if (global) {
                //translate into local index ...
                CellIdx -= m_map.i0;

                int idxfound = 0;
                if (m_map.i0 <= CellIdx && m_map.iE >= CellIdx)
                    idxfound = 1;
                Debug.Assert(idxfound.MPISum()==1);

                if (m_map.i0 > CellIdx || m_map.iE < CellIdx) {
                    this.m_CellSelInstruction = GetDoNothingInstruction();
                    return this;
                }
            }
            this.m_CellSelInstruction = GetIntInstruction(CellIdx);
            return this;
        }

        /// <summary>
        /// Selects List of cells according to global/local index
        /// </summary>
        /// <param name="ListOfCellIdx"></param>
        /// <param name="global"></param>
        /// <returns></returns>
        public SubBlockSelector CellSelector(List<int> ListOfCellIdx, bool global=true) {
            int LocNoOfBlocks = m_map.LocalNoOfBlocks;
            int GlobNoOfBlocks = m_map.TotalNoOfBlocks;

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
                    tmpIdx -= m_map.i0;

                    int idxfound = 0;
                    if (m_map.i0 <= tmpIdx && m_map.iE >= tmpIdx)
                        idxfound = 1;
                    Debug.Assert(idxfound.MPISum()==1);

                    if(idxfound==1)
                    tmpList.Add(tmpIdx);
                }
            } else {
                tmpList = ListOfCellIdx;
            }
            this.m_CellSelInstruction = GetListInstruction(tmpList);
            return this;
        }

        /// <summary>
        /// TBD ...
        /// </summary>
        /// <param name="CM"></param>
        /// <returns></returns>
        public SubBlockSelector CellSelector(CellMask CM) {
            throw new NotImplementedException();
            this.m_CellSelInstruction = delegate (int iCell) {

                return false;
            };
            return this;
        }
        #endregion

        #region VariableSelector
        /// <summary>
        /// Selects all Variable blocks
        /// </summary>
        /// <returns></returns>
        public SubBlockSelector VariableSelector() {
            this.m_VariableSelInstruction = delegate (int iCell, int iVar) {
                int NoOfVar = m_AggBS.Length;
                return GetAllInstruction(NoOfVar)(iVar);
            };
            return this;
        }

        public SubBlockSelector VariableSelector(List<int> ListOfVariables) {
            int NoOfVar = m_AggBS.Length;
            int maxValue = ListOfVariables.ToArray().Max();
            int minValue = ListOfVariables.ToArray().Min();
            if ((maxValue >= NoOfVar) || (minValue < 0))
                throw new ArgumentOutOfRangeException();

            this.m_VariableSelInstruction = delegate (int iCell, int iVar) {

                return GetListInstruction(ListOfVariables)(iVar);
            };
            return this;
        }
        #endregion

        #region SpeciesSelector
        /// <summary>
        /// Selects all species blocks
        /// </summary>
        /// <returns></returns>
        public SubBlockSelector SpeciesSelector() {
            this.m_SpeciesSelInstruction = delegate (int iCell, int iVar, int iSpec) {
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
        public SubBlockSelector SpeciesSelector(int Species) {
            // ist position von species in matrix consistent? Species A ist immer an erster Stelle, etc.?
            this.m_SpeciesSelInstruction = delegate (int iCell, int iVar, int iSpec) {
                return GetIntInstruction(Species)(iSpec);
            };
            return this;
        }

        /// <summary>
        /// not implemented yet
        /// </summary>
        /// <param name="Species"></param>
        /// <returns></returns>
        public SubBlockSelector SpeciesSelector(char Species) {
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
        public SubBlockSelector ModeSelector() {
            this.m_ModeSelInstruction = delegate (int iCell, int iVar, int iSpec, int iMode) {
                int maxDG = m_DG[iVar];
                return GetAllInstruction(maxDG+1)(iMode);
            };
            return this;
        }

        /// <summary>
        /// Selects Modes according to instruction.
        /// </summary>
        /// <param name="boolinstruct"></param>
        /// <returns></returns>
        public SubBlockSelector ModeSelector(Func<int,bool> boolinstruct) {
            this.m_ModeSelInstruction = delegate (int iCell, int iVar, int iSpec, int iMode) {
                return boolinstruct(iMode);
            };
            return this;
        }

        #endregion

        #region BasisInstructions
        private Func<int, bool> GetListInstruction(List<int> intList) {
            //test for dublicates
            int[] tmp = intList.ToArray();
            if (tmp.GroupBy(x => x).Any(g => g.Count() > 1))
                throw new ArgumentException("List contains dublicates!");

            return delegate (int Index) {
                bool selectit = false;
                foreach (int cell in intList) {
                    selectit = (cell == Index);
                    if (selectit)
                        break;
                }
                return selectit;
            };
        }

        private Func<int, bool> GetIntInstruction(int thisint) {
            return delegate (int Index) {
                return thisint == Index;
            };
        }

        private Func<int, bool> GetAllInstruction(int MaxIdx) {
            //delegate, which always return true is also possible
            List<int> tmp = new List<int>();
            for (int iLoc = 0; iLoc < MaxIdx; iLoc++) {
                tmp.Add(iLoc);
            }
            return GetListInstruction(tmp);
        }

        private Func<int, bool> GetDoNothingInstruction() {
            return delegate (int idx) {
                return false;
            };
        }

        #endregion

        #region Getta
        public Func<int, int, int, int, bool> GetAllInstructions {
            get {
                return delegate (int iCell, int iVar, int iSpec, int iMod) {
                    bool result=false;
                    result &= m_CellSelInstruction(iCell);
                    result &= m_VariableSelInstruction(iCell,iVar);
                    result &= m_SpeciesSelInstruction(iCell, iVar, iSpec);
                    result &= m_ModeSelInstruction(iCell, iVar, iSpec, iMod);
                    return result;
                };
            }
        }

        public Func<int, bool> GetCellInstruction {
            get { return m_CellSelInstruction; }
        }

        public Func<int, int, bool> GetVarInstruction {
            get { return m_VariableSelInstruction; }
        }

        public Func<int, int ,int ,bool> GetSpecInstruction {
            get { return m_SpeciesSelInstruction; }
        }

        public Func<int, int, int, int,bool> GetModeInstruction {
            get { return m_ModeSelInstruction; }
        }


        /// <summary>
        /// gets the multigrid operator on which this selector shall work on 
        /// </summary>
        public MultigridMapping SetMapping {
            set {
                m_map=value;
            }
        }

        private AggregationGridBasis[] m_AggBS {
            get {
                if (m_map != null) {
                    return m_map.AggBasis;
                } else {
                    throw new Exception("No Mapping set.");
                }
            }
        }

        private int[] m_DG {
            get {
                if (m_map != null) {
                    return m_map.DgDegree;
                } else {
                    throw new Exception("No Mapping set.");
                }
            }
        }

        #endregion
    }

    /// <summary>
    /// The block mask, which are generated through this class, presumes that the target Multigridoperator(defined through SubBlockSelector) is a quadratic Matrix.
    /// </summary>
    public class BlockMask {

        /// <summary>
        /// auxiliary structure. Stores Offsets and Length of something.
        /// </summary>
        private struct extNi0 {
            public extNi0(int Li0, int Gi0, int Si0, int N) {
                m_Li0 = Li0;
                m_Gi0 = Gi0;
                m_Si0 = Si0;
                m_N = N;
            }
            private int m_Li0, m_Gi0 ,m_Si0, m_N;
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
        /// Generates Block Mask from Sub block selection. Sub matrix can be generated of the mutligrid operator overgiven to Sub block selection. It also contains methods to translate Vectors from full matrix to sub matrix and vice versa.
        /// </summary>
        /// <param name="SBS"></param>
        public BlockMask(SubBlockSelector SBS, MultigridMapping map) {
            SBS.SetMapping=map;
            m_sbs = SBS;
            m_AggBS = m_map.AggBasis;
            m_DG = m_map.DgDegree;
            m_Ni0 = Ni0Gen();
            m_NoOfCells = m_map.LocalNoOfBlocks;
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

        //internal get
        private SubBlockSelector m_sbs;
        private MultigridMapping m_map;
        private AggregationGridBasis[] m_AggBS;
        private int[] m_DG;
        private int m_NoOfCells;
        private int m_NoOfVariables;
        private int[][] m_NoOfSpecies;
        //internal set
        private List<int> m_GlobalMask = null;
        private List<int> m_LocalMask = null;
        private extNi0[][][][] m_StructuredNi0 = null;
        private int m_Ni0Len;
        private int m_MaskLen;
        private Ni0[] m_Ni0;

        private bool TestForQuadraticMatrix() {
            throw new NotSupportedException("Supports only quadratic matrices");
        }

        private void GenerateAllMasks() {
            int NoOfCells = m_NoOfCells;
            int NoOfVariables = m_NoOfVariables;
            int[][] NoOfSpecies = m_NoOfSpecies;
            int[] NoOfDG = m_DG.CloneAs();
            for(int iDG=0;iDG < NoOfDG.Length; iDG++) {
                NoOfDG[iDG] += 1;
            }

            List<extNi0> ListNi0 = new List<extNi0>();
            List<int> Globalint = new List<int>();
            List<int> Localint = new List<int>();
            List<int> SubBlockIdx = new List<int>();

            int SubOffset = 0;
            int Ni0Length = 0;
            int MaskLen = 0;
            var tmpCell = new List<extNi0[][][]>();
            for (int jLoc = 0; jLoc < NoOfCells; jLoc++) {
                if (!m_sbs.GetCellInstruction(jLoc))
                    continue;
                var tmpVar = new List<extNi0[][]>();
                for (int iVar = 0; iVar < NoOfVariables; iVar++) {
                    if (!m_sbs.GetVarInstruction(jLoc, iVar))
                        continue;
                    var tmpSpc = new List<extNi0[]>();
                    for (int iSpc = 0; iSpc < NoOfSpecies[jLoc][iVar]; iSpc++) {
                        if (!m_sbs.GetSpecInstruction(jLoc, iVar, iSpc))
                            continue;
                        int GlobalOffset = m_map.GlobalUniqueIndex(iVar, jLoc, iSpc, 0);
                        int LocalOffset = m_map.LocalUniqueIndex(iVar, jLoc, iSpc, 0);
                        var tmpMod = new List<extNi0>();
                        for (int iMode = 0; iMode < NoOfDG[iVar]; iMode++) {
                            if (m_sbs.GetModeInstruction(jLoc, iVar, iSpc, iMode)) {
                                int GlobalModeOffset = m_Ni0[iMode].i0 + GlobalOffset;
                                int LocalModeOffset = m_Ni0[iMode].i0 + LocalOffset;
                                int ModeLength = m_Ni0[iMode].N;
                                var newNi0 = new extNi0(LocalModeOffset, GlobalModeOffset,SubOffset, ModeLength);
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
                                Debug.Assert(m_map.LocalUniqueIndex(iVar, jLoc, iSpc, GetNp(iMode)-1) == LocalModeOffset + ModeLength-1);
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
        private int[] GetCellwiseLocalidx(int iCell) {
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
        private int GetCellwiseLength(int iCell) {
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

        private List<int> GetAllSubMatrixCellLength() {
            var intList = new List<int>();
            for (int iCell = 0; iCell < m_StructuredNi0.Length; iCell++) {
                intList.Add(GetCellwiseLength(iCell));
            }
            return intList;
        }

        private int GetCellwiseSubIdx(int iCell) {
            return m_StructuredNi0[iCell][0][0][0].Si0;
        }

        private List<int> GetAllSubMatrixCellOffsets() {
            var intList = new List<int>();
            for (int iCell = 0; iCell < m_StructuredNi0.Length;iCell++)
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
        private int GetRelativeSubBlockOffset(int iCell, int iVar, int iSpc, int iMod) {
            return m_StructuredNi0[iCell][iVar][iSpc][iMod].Si0 - m_StructuredNi0[iCell][0][0][0].Si0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        private extNi0[][][] GetSubBlockNi0Cellwise(int iCell) {

            //int subOffset = 0;
            //List<Ni0[][]> tmpVar = new List<Ni0[][]>();
            //for (int iVar = 0; iVar < m_StructuredNi0[iCell].Length; iVar++) {
            //    List<Ni0[]> tmpSpc = new List<Ni0[]>();
            //    for (int iSpc = 0; iSpc < m_StructuredNi0[iCell][iVar].Length; iSpc++) {
            //        List<Ni0> tmpMod = new List<Ni0>();
            //        for (int iMod = 0; iMod < m_StructuredNi0[iCell][iVar][iSpc].Length; iMod++) {
            //            extNi0 block = m_StructuredNi0[iCell][iVar][iSpc][iMod];
            //            int subN = block.N;
            //            var subNi0 = new Ni0(subOffset, subN);
            //            tmpMod.Add(subNi0);
            //            subOffset += subN;
            //        }
            //        tmpSpc.Add(tmpMod.ToArray());
            //    }
            //    tmpVar.Add(tmpSpc.ToArray());
            //}

            return m_StructuredNi0[iCell];
        }


        /// <summary>
        /// Pre generates Offsets and Length for a default Variable block with the maximal DGdegree, to save computation time.
        /// </summary>
        /// <returns></returns>
        private Ni0[] Ni0Gen() {
            int maxDG = m_DG.Max();
            var ModeOffsetNLengt = new List<Ni0>();
            for (int p = 0; p < maxDG+1; p++) {
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
            Debug.Assert(p>=0);
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
        /// Returns the total number of subblocks row or columnwise.
        /// </summary>
        public int TotalLength {
            get {
                return m_Ni0Len.MPISum();
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

        /// <summary>
        /// returns total number of DOF of Sub Matrix
        /// </summary>
        public int TotalDOF {
            get {
                return m_MaskLen.MPISum();
            }
        }

        /// <summary>
        /// Get Array of Cellblocks
        /// </summary>
        /// <param name="target"></param>
        /// <param name="ignoreCellCoupling"></param>
        /// <param name="ignoreVarCoupling"></param>
        /// <param name="ignoreSpecCoupling"></param>
        /// <returns></returns>
        public MultidimensionalArray[] GetSubBlocks(BlockMsrMatrix target, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {
            int NoOfCells = m_StructuredNi0.Length;
            int size = ignoreCellCoupling ? NoOfCells : NoOfCells * NoOfCells;

            MultidimensionalArray[] Sblocks = new MultidimensionalArray[size];

            int auxIdx = 0;
            for (int iLoc = 0; iLoc < m_StructuredNi0.Length; iLoc++) {
                for (int jLoc = 0; jLoc < m_StructuredNi0.Length; jLoc++) {
                    if (ignoreCellCoupling && jLoc != iLoc) {
                        continue;
                    }
                    int CellBlockLen = GetCellwiseLength(jLoc);
                    Sblocks[auxIdx] = MultidimensionalArray.Create(CellBlockLen, CellBlockLen);
                    for (int iVar = 0; iVar < m_StructuredNi0[iLoc].Length; iVar++) {
                        for (int jVar = 0; jVar < m_StructuredNi0[jLoc].Length; jVar++) {
                            if (ignoreVarCoupling && jVar != iVar) {
                                continue;
                            }
                            for (int iSpc = 0; iSpc < m_StructuredNi0[iLoc][iVar].Length; iSpc++) {
                                for (int jSpc = 0; jSpc < m_StructuredNi0[jLoc][jVar].Length; jSpc++) {
                                    if (ignoreSpecCoupling && jSpc != iSpc) {
                                        continue;
                                    }
                                    for (int iMode = 0; iMode < m_StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
                                        for (int jMode = 0; jMode < m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {
                                            extNi0 RowNi0 = m_StructuredNi0[iLoc][iVar][iSpc][iMode];
                                            extNi0 ColNi0 = m_StructuredNi0[jLoc][jVar][jSpc][jMode];
                                            int Targeti0 = RowNi0.Gi0;
                                            int Targetj0 = ColNi0.Gi0;
                                            int Subi0 = GetRelativeSubBlockOffset(iLoc, iVar, iSpc, iMode);
                                            int Subj0 = GetRelativeSubBlockOffset(jLoc, jVar, jSpc, jMode);
                                            int Subie = Subi0 + RowNi0.N-1;
                                            int Subje = Subj0 + ColNi0.N-1;

                                            target.ReadBlock(Targeti0, Targetj0,
                                                Sblocks[auxIdx].ExtractSubArrayShallow(new int[] { Subi0, Subj0 }, new int[] { Subie, Subje }));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    auxIdx++;
                }
            }
            return Sblocks;
        }

        /// <summary>
        /// If you want nothing special. Take this one. If you want only diagonal block matrix choose one of the other methods
        /// </summary>
        /// <returns></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix target) {
            int Loclength = m_MaskLen;

            var tmpN = GetAllSubMatrixCellLength();
            var tmpi0 = GetAllSubMatrixCellOffsets();

            BlockPartitioning localBlocking = new BlockPartitioning(Loclength, tmpi0.ToArray(), tmpN.ToArray(), m_map.MPI_Comm, i0isLocal: true);
            var SubMSR = new BlockMsrMatrix(localBlocking);
            target.AccSubMatrixTo(1.0, SubMSR, m_GlobalMask, default(int[]), m_GlobalMask, default(int[]));
            return SubMSR;
        }

        /// <summary>
        /// Coupling can be ignored. If you want full output take the basic GetSubBlockMatrix() method, which will be faster.
        /// </summary>
        /// <param name="target"></param>
        /// <param name="ignoreCellCoupling"></param>
        /// <param name="ignoreVarCoupling"></param>
        /// <param name="ignoreSpecCoupling"></param>
        /// <returns></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix target, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {
            int Loclength = m_MaskLen;

            var tmpN = GetAllSubMatrixCellLength();
            var tmpi0 = GetAllSubMatrixCellOffsets();

            BlockPartitioning localBlocking = new BlockPartitioning(Loclength, tmpi0.ToArray(), tmpN.ToArray(), m_map.MPI_Comm, i0isLocal: true);
            var SubMSR = new BlockMsrMatrix(localBlocking);

            int SubRowIdx = 0;
            int SubColIdx = 0;

            int auxIdx = 0;
            for (int iLoc = 0; iLoc < m_StructuredNi0.Length; iLoc++) {
                for (int jLoc = 0; jLoc < m_StructuredNi0.Length; jLoc++) {
                    if (ignoreCellCoupling && jLoc != iLoc) {
                        continue;
                    }
                    for (int iVar = 0; iVar < m_StructuredNi0[iLoc].Length; iVar++) {
                        for (int jVar = 0; jVar < m_StructuredNi0[jLoc].Length; jVar++) {
                            if (ignoreVarCoupling && jVar != iVar) {
                                continue;
                            }
                            for (int iSpc = 0; iSpc < m_StructuredNi0[iLoc][iVar].Length; iSpc++) {
                                for (int jSpc = 0; jSpc < m_StructuredNi0[jLoc][jVar].Length; jSpc++) {
                                    if (ignoreSpecCoupling && jSpc != iSpc) {
                                        continue;
                                    }
                                    for (int iMode = 0; iMode < m_StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
                                        for (int jMode = 0; jMode < m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {

                                            int[] RowIdx = new int[] { iLoc, iVar, iSpc, iMode };
                                            int[] ColIdx = new int[] { jLoc, jVar, jSpc, jMode };
                                            extNi0 RowNi0 = m_StructuredNi0[iLoc][iVar][iSpc][iMode];
                                            extNi0 ColNi0 = m_StructuredNi0[jLoc][jVar][jSpc][jMode];
                                            int Targeti0 = RowNi0.Gi0;
                                            int Targetj0 = ColNi0.Gi0;

                                            var tmpBlock = MultidimensionalArray.Create(RowNi0.N, ColNi0.N);
                                            target.ReadBlock(Targeti0, Targetj0,
                                                tmpBlock);
                                            SubMSR.AccBlock(SubRowIdx, SubColIdx, 1, tmpBlock);
                                            SubColIdx += ColNi0.N;
                                        }
                                        SubRowIdx += m_StructuredNi0[iLoc][iVar][iSpc][iMode].N;
                                    }
                                }
                            }
                        }
                    }
                    auxIdx++;
                }
                auxIdx++;
            }
            return SubMSR;
        }

        /// <summary>
        /// Translates back your masked vector to the full one
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="targetVector">output</param>
        /// <param name="accVector"></param>
        public void AccVecToFull<V, W>(W accVector, V targetVector)
            where V : IList<double>
            where W : IList<double> {
            if (targetVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector != Length of original");
            if (accVector.Count() != m_MaskLen)
                throw new ArgumentException("accVector length is not equal to length of mask");
            targetVector.AccV(1.0, accVector, m_LocalMask, default(int[]));
        }

        /// <summary>
        /// Translates back your masked vector to the full one cellwise. i.e. Useful, if you solved blockwise
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="targetVector">output</param>
        /// <param name="accVector"></param>
        /// <param name="iCell">starts with 0</param>
        public void AccVecCellwiseToFull<V, W>(W accVector, int iCell, V targetVector)
            where V : IList<double>
            where W : IList<double> {
            if (iCell > m_StructuredNi0.Length - 1)
                throw new ArgumentOutOfRangeException("iCell is greater than Cells in mask");
            if (accVector.Count() != GetCellwiseLength(iCell))
                throw new ArgumentException("accVector length is not equal to length of mask");
            if(targetVector.Count()!= m_map.LocalLength)
                throw new ArgumentException("Length of targetVector not equal Length of original");
            var Cidx=GetCellwiseLocalidx(iCell);
            Debug.Assert(accVector.Count()==Cidx.Count());
            targetVector.AccV(1.0, accVector, Cidx, default(int[]));
        }

        /// <summary>
        /// Translates a full vector to a vector corresponding to matrix mask cellwise
        /// </summary>
        /// <param name="fullVector"></param>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public double[] GetVectorCellwise(IList<double> fullVector, int iCell) {
            if (iCell > m_StructuredNi0.Length - 1)
                throw new ArgumentOutOfRangeException("iCell is greater than number of cells in mask");
            if (fullVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector not equal Length of original");
            double[] subVector = new double[GetCellwiseLength(iCell)];
            var Cidx=GetCellwiseLocalidx(iCell);
            Debug.Assert(subVector.Length == Cidx.Length);
            ArrayTools.GetSubVector<int[], int[], double>(fullVector, subVector, Cidx);
            return subVector;
        }

        /// <summary>
        /// Translates a full vector to a vector corresponding to matrix mask
        /// </summary>
        /// <param name="fullVector"></param>
        public double[] GetSubBlockVec(IList<double> fullVector) {
            if (fullVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector not equal Length of original");
            int Len = m_MaskLen;
            double[] subVector = new double[m_MaskLen];
            subVector.AccV(1.0, fullVector, default(int[]), m_LocalMask);
            return subVector;
        }

    }
}
