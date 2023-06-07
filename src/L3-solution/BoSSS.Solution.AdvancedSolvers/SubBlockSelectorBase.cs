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
using System.Collections;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Common interface for **coordinate mappings**
    /// </summary>
    public interface ICoordinateMapping : IBlockPartitioning {

        /// <summary>
        /// 
        /// </summary>
        int NoOfVariables { get; }

        /// <summary>
        /// true, if the <paramref name="iVar"/>-th variable is XDG; otherwise, it must be DG
        /// </summary>
        bool IsXDGvariable(int iVar);

        /// <summary>
        /// index of species <paramref name="SId"/> in cell <paramref name="jCell"/>
        /// </summary>
        int GetSpeciesIndex(int jCell, SpeciesId SId);

        /// <summary>
        /// 
        /// </summary>
        int GetNoOfSpecies(int jCell);

        /// <summary>
        /// DG polynomial degree (aka order) for each variable
        /// </summary>
        int[] DgDegree {
            get;
        }

        /// <summary>
        /// Number of locally stored external cells - no computations are carried out for
        /// that cells, but their values are needed.
        /// 
        /// see also <see cref="ILogicalCellData.NoOfExternalCells"/>
        /// </summary>
        int NoOfExternalCells {
            get;
        }

        /// <summary>
        /// Number of degrees-of-freedom per cell
        /// </summary>
        /// <param name="jLoc">local cell index</param>
        /// <returns></returns>
        int GetLength(int jLoc);


        /// <summary>
        /// Mapping from a quadruple (<paramref name="ifld"/>,<paramref name="jCell"/>,<paramref name="iSpec"/>,<paramref name="n"/>) to a 
        /// MPI-local linear index range 
        /// </summary>
        /// <param name="ifld">field/variable index</param>
        /// <param name="jCell">cell index</param>
        /// <param name="n">DG/XDG mode index</param>
        /// <param name="iSpec">XDG species index</param>
        /// <returns>
        /// local index, i.e. starting at 0 on all MPI processes 
        /// </returns>
        int LocalUniqueIndex(int ifld, int jCell, int iSpec, int n);

        /// <summary>
        /// MPI-global version of <see cref="LocalUniqueIndex(int, int, int, int)"/>
        /// </summary>
        /// <returns>
        /// A global index, i.e. it a different range, 
        /// from <see cref="IPartitioning.i0"/> (including) to <see cref="IPartitioning.iE"/> (excluding),
        /// on each MPI process.
        /// </returns>
        long GlobalUniqueIndex(int ifld, int jCell, int jSpec, int n);

        /// <summary>
        /// Mapping from a triple (<paramref name="ifld"/>,<paramref name="jCell"/>,<paramref name="n"/>) to a linear index range 
        /// </summary>
        /// <param name="ifld">field/variable index</param>
        /// <param name="jCell">cell index</param>
        /// <param name="n">DG/XDG mode index</param>
        /// <returns>
        /// local index, i.e. starting at 0 on all MPI processes 
        /// </returns>
        int LocalUniqueIndex(int ifld, int jCell, int n);

        /// <summary>
        /// MPI-global version of <see cref="LocalUniqueIndex(int, int, int)"/>
        /// </summary>
        /// <returns>
        /// A global index, i.e. it a different range, 
        /// from <see cref="IPartitioning.i0"/> (including) to <see cref="IPartitioning.iE"/> (excluding),
        /// on each MPI process.
        /// </returns>
        long GlobalUniqueIndex(int ifld, int jCell, int n);

        /// <summary>
        /// 1D, 2D or 3D;
        /// </summary>
        int SpatialDimension { get;  }


        /// <summary>
        /// Alias for <see cref="IBlockPartitioning.LocalNoOfBlocks"/>:
        /// Number of locally updated cells - the cells which are computed on
        /// this processor (in contrast, see <see cref="NoOfExternalCells"/>); 
        /// see also <see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>
        /// </summary>
        int NoOfLocalUpdatedCells {
            get;
        }
                

        /// <summary>
        /// <see cref="NoOfExternalCells"/> plus <see cref="NoOfLocalUpdatedCells"/>; see also <see cref="ILogicalCellData.Count"/>
        /// </summary>
        int LocalCellCount {
            get;
        }

        /// <summary>
        /// All used XDG species.
        /// Index: enumeration over species.
        /// </summary>
        SpeciesId[] UsedSpecies {
            get;
        }
    }


    /// <summary>
    /// This contains a block Selection, which can be specified by the user with the hierachical Selectors: every Matrix is subdivided into cell, variable, species and mode blocks.
    /// </summary>
    public abstract class SubBlockSelectorBase {

        /// <summary>
        /// Specifies, which blocks in a matrix shall be selected. Blocksubdivision Default: Selects all blocks.  
        /// </summary>
        /// <param name="map"></param>
        public SubBlockSelectorBase(ICoordinateMapping map) {
            if (map == null)
                throw new ArgumentNullException("empty mapping! This will not end well ...");
            m_map = map;
            this.CellSelector();
            this.SetVariableSelector();
            this.SetSpeciesSelector();
            this.SetModeSelector();
        }

        //internal get
        protected ICoordinateMapping m_map;

        /// <summary>
        /// Selector for cells.
        /// - 1st argument: local cell index
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        protected Func<int, bool> m_CellFilter = null;

        /// <summary>
        /// Selector for cells and variables.
        /// - 1st argument: local cell index
        /// - 2nd argument: variable index
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        protected Func<int, int, bool> m_VariableFilter = null;

        /// <summary>
        /// Selector for cells and variables and species.
        /// - 1st argument: cell index
        /// - 2nd argument: variable index
        /// - 3rd argument: species index
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        protected Func<int, int, int, bool> m_SpeciesFilter = null;

        /// <summary>
        /// Selector for cells and variables and species.
        /// - 1st argument: cell index
        /// - 2nd argument: variable index
        /// - 3rd argument: species index
        /// - 4th argument: DG polynomial degree
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
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
        /// <returns></returns>
        public SubBlockSelectorBase CellSelector(int CellIdx) {

            int LocNoOfBlocks = m_NoLocalCells;
            long GlobNoOfBlocks = m_NoTotalCells;
//            if (global) {
//               if (CellIdx >= GlobNoOfBlocks || CellIdx < 0)
//                    throw new ArgumentOutOfRangeException(CellIdx + " is greater then global No of Blocks: " + GlobNoOfBlocks + " or smaller than 0");
//            } else {
                if (CellIdx >= LocNoOfBlocks || CellIdx < 0)
                    throw new ArgumentOutOfRangeException(CellIdx + " is greater then Local No of Blocks: " + LocNoOfBlocks + " or smaller than 0");
 //           }

            //if (global) {
            //    //translate into local index ...
            //    CellIdx -= m_i0;

            //    int idxfound = 0;
            //    if (m_i0 <= CellIdx && m_iE >= CellIdx)
            //        idxfound = 1;
            //    Debug.Assert(idxfound.MPISum() == 1);

            //    if (m_i0 > CellIdx || m_iE < CellIdx) {
            //        this.m_CellFilter = GetDoNothingInstruction();
            //        return this;
            //    }
            //}

            this.m_CellFilter = GetIntInstruction(CellIdx);
            return this;
        }

        /// <summary>
        /// Selects List of cells according to global/local index
        /// </summary>
        /// <param name="ListOfCellIdx"></param>
        /// <param name="global"></param>
        /// <returns></returns>
        public SubBlockSelectorBase CellSelector<V>(V ListOfCellIdx, bool global = false)
            where V : IList<int> //
        {
            int LocNoOfBlocks = m_NoLocalCells;
            long GlobNoOfBlocks = m_NoTotalCells;

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
                    int tmpIdx = (int)(CellIdx - Cell_j0);

                    int idxfound = 0;
                    if (Cell_j0 <= tmpIdx && Cell_jE >= tmpIdx)
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
        public SubBlockSelectorBase SetVariableSelector() {
            this.m_VariableFilter = delegate (int iCell, int iVar) {
#if DEBUG
                int NoOfVar = m_map.NoOfVariables;
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
        public SubBlockSelectorBase SetVariableSelector(params int[] SetOfVariables) {
            return SetVariableSelector((ICollection<int>)SetOfVariables);
        }

        /// <summary>
        /// 
        /// </summary>
        public SubBlockSelectorBase SetVariableSelector(int iVariable) {
            if (iVariable < 0)
                throw new ArgumentOutOfRangeException("Variable index cannot be negative.");
            if (iVariable >= m_NoOfVar)
                throw new ArgumentOutOfRangeException("Variable index is larger than number of variables..");

            return SetVariableSelector(new int[] { iVariable });
        }


        /// <summary>
        /// 
        /// </summary>
        public SubBlockSelectorBase SetVariableSelector(IEnumerable<int> SetOfVariables) {
            int[] Variables = SetOfVariables.ToArray();
            if (!Variables.IsSet())
                throw new ArgumentOutOfRangeException("Input is not a set - some variable index is listed twice.");
            if (Variables.Min() < 0)
                throw new ArgumentOutOfRangeException("Variable index cannot be negative.");
            if (Variables.Max() >= m_NoOfVar)
                throw new ArgumentOutOfRangeException("Some variable index is larger than number of variables..");

            int NoOfVar = m_map.NoOfVariables;

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
        public SubBlockSelectorBase SetSpeciesSelector() {
            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                int NoOfSpec = m_map.GetNoOfSpecies(iCell);
                return GetAllInstruction(NoOfSpec)(iSpec);
            };
            return this;
        }

        /// <summary>
        /// Selects Species by <see cref="SpeciesId"/>.
        /// </summary>
        public SubBlockSelectorBase SetSpeciesSelector(SpeciesId SId) {

            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                if (!this.m_map.IsXDGvariable(iVar))
                    throw new NotSupportedException("You tried to select a species within a non-xdg field!");
                int SpcIdx = this.m_map.GetSpeciesIndex(iCell, SId);
                if (SpcIdx < 0)
                    return GetDoNothingInstruction()(iSpec);
                else
                    return GetIntInstruction(SpcIdx)(iSpec);
            };

            return this;
        }

        /*
        /// <summary>
        /// Selects Species by <see cref="SpeciesId"/>.
        /// </summary>
        public SubBlockSelectorBase SpeciesSelector(string Species) {

            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                if (!this.m_map.IsXDGvariable(iVar))
                    throw new NotSupportedException("You tried to select a species within a non-xdg field!");
                SpeciesId SpecId = this.m_map.GetSpeciesId(Species);
                int SpcIdx = this.m_map.GetSpeciesIndex(iCell, SpecId);
                if (SpcIdx < 0)
                    return GetDoNothingInstruction()(iSpec);
                else
                    return GetIntInstruction(SpcIdx)(iSpec);
            };

            return this;
        }
        */

        /// <summary>
        /// Selects multiple species by <see cref="IEnumerable{SpeciesId}"/>.
        /// </summary>
        public SubBlockSelectorBase SetSpeciesSelector(IEnumerable<SpeciesId> SetOfSpecies) {
            //for (int v = 0; v < m_map.NoOfVariables; v++) {
            //    if (this.m_map.AggBasis[v].GetType() == typeof(XdgAggregationBasis))
            //        Console.WriteLine("WARNING: Variable {0} has no XdgBasis and thus may be not selected",v);
            //}
            this.m_SpeciesFilter = delegate (int iCell, int iVar, int iSpec) {
                List<int> SpcInt = new List<int>();
                foreach(SpeciesId spi in SetOfSpecies) {
                    int SpcIdx = m_map.GetSpeciesIndex(iCell, spi);
                    if (SpcIdx <0)
                        continue;
                    SpcInt.Add(SpcIdx);
                }
                if (SpcInt.Count <= 0)
                    return GetDoNothingInstruction()(iSpec);
                else
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
        public SubBlockSelectorBase SetModeSelector() {
            this.m_ModeFilter = delegate (int iCell, int iVar, int iSpec, int pDeg) {
#if DEBUG
                int maxDG = m_DGdegree[iVar];
                Debug.Assert(pDeg <= maxDG);
#endif
                return true;
            };
            return this;
        }

        /// <summary>
        /// Selects Modes according to instruction.
        /// </summary>
        public SubBlockSelectorBase SetModeSelector(Func<int, bool> boolinstruct) {
            this.m_ModeFilter = delegate (int iCell, int iVar, int iSpec, int pDeg) {
                return boolinstruct(pDeg);
            };
            return this;
        }

        /// <summary>
        /// Selects Modes according to instruction.
        /// </summary>
        public SubBlockSelectorBase SetModeSelector(Func<int, int, int, int, bool> boolinstruct) {
            this.m_ModeFilter = boolinstruct;
            return this;
        }


        #endregion

        #region BasisInstructions

        private Func<int, bool> GetListInstruction(IEnumerable<int> intList) {
#if DEBUG
            if(!intList.IsSet())
                throw new ArgumentException("List contains duplicates!");
#endif
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

        /// <summary>
        /// Selector for cells.
        /// - 1st argument: local cell index
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        public Func<int, bool> CellFilter {
            get { return m_CellFilter; }
        }

        /// <summary>
        /// Selector for cells and variables.
        /// - 1st argument: local cell index
        /// - 2nd argument: variable index
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        public Func<int, int, bool> VariableFilter {
            get { return m_VariableFilter; }
        }

        /// <summary>
        /// Selector for cells and variables and species.
        /// - 1st argument: cell index
        /// - 2nd argument: variable index
        /// - 3rd argument: species index
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        public Func<int, int, int, bool> SpeciesFilter {
            get { return m_SpeciesFilter; }
        }


        /// <summary>
        /// Selector for cells and variables and species.
        /// - 1st argument: cell index
        /// - 2nd argument: variable index
        /// - 3rd argument: species index
        /// - 4th argument: DG polynomial degree
        /// - return value: respective DOFs should be included (true) or excluded (false)
        /// </summary>
        public Func<int, int, int, int, bool> ModeFilter {
            get { return m_ModeFilter; }
        }


        /// <summary>
        /// gets the multigrid operator on which this selector shall work on 
        /// </summary>
        public ICoordinateMapping Mapping {
            get {
                return m_map;
            }
        }

        /*
        private AggregationGridBasis[] m_AggBS {
            get {
                return m_map.AggBasis;
            }
        }
        */

#if DEBUG
        private int[] m_DGdegree {
            get {
                return m_map.DgDegree;
            }
        }
#endif
        private int m_NoLocalCells {
            get {
                return m_map.LocalNoOfBlocks + m_map.NoOfExternalCells;
            }
        }

        private long m_NoTotalCells {
            get {
                return m_map.TotalNoOfBlocks;
            }
        }
        private int m_NoOfVar {
            get {
                return m_map.NoOfVariables;
            }
        }

        private long Cell_j0 {
            get {
                return m_map.FirstBlock;
            }
        }

        private long Cell_jE {
            get {
                return m_map.LocalNoOfBlocks + Cell_j0;
            }
        }

        

        #endregion
    }


    public struct extNi0
    {
        /// <summary>
        /// stores offsets and length of a sub block.
        /// <paramref name="Li0"/> refers to local index of <see cref="IBlockPartitioning"/>.
        /// <paramref name="Gi0"/> refers to global index of <see cref="IBlockPartitioning"/>.
        /// <paramref name="Si0"/> is a local index,
        /// which numbers the entries of this mask consecutively.
        /// </summary>
        /// <param name="Li0">local offset</param>
        /// <param name="Gi0">global offset</param>
        /// <param name="Si0">offset relative to mask</param>
        /// <param name="N">length of block</param>
        public extNi0(int Li0, long Gi0, int Si0, int N) {
            m_Li0 = Li0;
            m_Gi0 = Gi0;
            m_Si0 = Si0;
            m_N = N;
        }
        private int m_Li0, m_Si0, m_N;
        private long m_Gi0;
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
        public long Gi0 {
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
    /// Therefore <see cref="m_NoOfCells"/> and <see cref="m_CellOffset"/> have to be overridden by inheriting classes.
    /// There are two inheriting classes: <see cref="BlockMask.BlockMaskExt"/> and <see cref="BlockMask.BlockMaskLoc"/>,
    /// which handle the masking of external cells and internal cells respectively.
    /// The purpose of this class is to provide the index lists and structs,
    /// which correspond to underlying selection through the <see cref="SubBlockSelector"/>.
    /// </summary>
    public abstract class BlockMaskBase {

        /// <summary>
        /// auxiliary structure. Stores Offsets and Length of subblocks.
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
        /// Generates Block Mask (index lists) from Sub block selection based on a multigrid mapping. 
        /// abstract parts are individualized by child classes: <see cref="BlockMask.BlockMaskLoc"/> and <see cref="BlockMask.BlockMaskExt"/>
        /// </summary>
        public BlockMaskBase(SubBlockSelector SBS, MPI_Comm MPIcomm) {
            m_sbs = SBS;
            m_DGdegree = m_map.DgDegree;
            m_Ni0 = Ni0Gen();
            m_MPIcomm = MPIcomm;
            //Testen ob es cells gibt wo Var<>NoOfVar, das würde dementsprechend auch m_DG beeinflussen
            m_NoOfSpecies = new int[m_NoOfCells][];
            for (int iCell = 0; iCell < m_NoOfCells; iCell++) {
                m_NoOfSpecies[iCell] = new int[m_NoOfVariables];
                for (int iVar = 0; iVar < m_NoOfVariables; iVar++) {
                    m_NoOfSpecies[iCell][iVar] = m_map.GetNoOfSpecies(iCell + m_CellOffset);
                }
            }
        }

        // the subsequent members are different
        // dependent on mask type: local or external

        protected abstract int m_NoOfCells {
            get;
        }

        protected abstract int m_CellOffset {
            get;
        }

        protected abstract int m_LocalLength {
            get;
        }

        protected abstract int m_SubBlockOffset {
            get;
        }

        // ============
        // internal get
        // ============
        private SubBlockSelectorBase m_sbs;

        protected ICoordinateMapping m_map => m_sbs.Mapping;

        //private AggregationGridBasis[] m_AggBS;
        private int[] m_DGdegree;
        private int m_NoOfVariables => m_map.NoOfVariables;
        
        private int[][] m_NoOfSpecies;

        private MPI_Comm m_MPIcomm;

        // ============
        // internal set
        // ============
        #region internal set

        /// <summary>
        /// global indices in this mask
        /// </summary>
        public List<long> m_GlobalMask = null;

        /// <summary>
        /// local indices in this mask
        /// </summary>
        public List<int> m_LocalMask = null;

        /// <summary>
        /// subblock indices in this mask
        /// </summary>
        public List<int> m_SubBlockMask = null;

        /// <summary>
        /// stores offsets (of local, global and subblock numbering) and lengths of dg blocks
        /// structure mimics subblock hierarchy:
        /// - 1st idx : cells
        /// - 2nd idx : variables
        /// - 3rd idx : species
        /// - 4th idx : dg blocks
        /// content : offset (of local, global and subblock numbering) and length of dg block
        /// </summary>
        public extNi0[][][][] m_StructuredNi0 = null;

        /// <summary>
        /// number of DG subblocks in this mask
        /// </summary>
        public int m_Ni0Len;

        /// <summary>
        /// DOF within this mask
        /// </summary>
        private int m_MaskLen;

        /// <summary>
        /// same content as <see cref="m_StructuredNi0"/> but stored in a plain array
        /// </summary>
        private Ni0[] m_Ni0;
        #endregion


        /// <summary>
        /// Pre generate offsets and lengths for the DG blocks.
        /// Used in <see cref="GenerateAllMasks"/>.
        /// </summary>
        /// <returns>Ni0 for a default Variable up to the maximal occurring DGdegree</returns>
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
            int SpacDim = m_map.SpatialDimension;
            Debug.Assert(p >= 0);
            switch (SpacDim) {
                case 1: return p + 1;
                case 2: return (p * p + 3 * p + 2) / 2;
                case 3: return (p * p * p + 6 * p * p + 11 * p + 6) / 6;
                default: throw new Exception("wtf?Spacialdim=1,2,3 expected");
            }
        }

        /// <summary>
        /// The core of the masking
        /// Generates index lists and the Ni0-struct-list corresponding to mask
        /// and is called by the child classes: local and external mask
        /// Note: the smallest sub blocks are DG blocks!
        /// </summary>
        protected void GenerateAllMasks() {
            int NoOfCells = m_NoOfCells;
            int NoOfVariables = m_NoOfVariables;
            int[][] NoOfSpecies = m_NoOfSpecies;
            int[] DGdegreeP1 = m_DGdegree.CloneAs();
            for (int iDG = 0; iDG < DGdegreeP1.Length; iDG++) {
                DGdegreeP1[iDG] += 1;
            }

            List<extNi0> ListNi0 = new List<extNi0>();
            List<long> Globalint = new List<long>();
            List<int> Localint = new List<int>();
            List<int> SubBlockIdx = new List<int>();

            int SubOffset = m_SubBlockOffset; // 0 for local mask and BMLoc.LocalDof for external mask
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
            bool emptysel = true; // for debugging, do not look at selections containing no cells, possible if external cells shall be considered but block has none
            int totNoOfSpeciesInSelection = 0; // for debugging, in case of IBM, there may be empty selection, due to 0 species in cells 

            //ilPSP.Environment.StdoutOnlyOnRank0 = false;

            // loop over cells...
            for (int iLoc = 0; iLoc < NoOfCells; iLoc++) {
                int jLoc = m_CellOffset + iLoc; //to address correctly, external cells offset has to be considered, you know ...
                emptysel &= !CellInstruction(jLoc); //for testing if the entire selection is empty, which hopefully only can happen at the level of cells
                if (!CellInstruction(jLoc))
                    continue;
                var tmpVar = new List<extNi0[][]>();

                // loop over variables...
                for (int iVar = 0; iVar < NoOfVariables; iVar++) {
                    if (!VarInstruction(jLoc, iVar))
                        continue;
                    var tmpSpc = new List<extNi0[]>();

                    // loop over species...
                    for (int iSpc = 0; iSpc < NoOfSpecies[iLoc][iVar]; iSpc++) {
                        totNoOfSpeciesInSelection += NoOfSpecies[iLoc][iVar];
                        if (!SpecInstruction(jLoc, iVar, iSpc))
                            continue;
                        int LocalOffset = m_map.LocalUniqueIndex(iVar, jLoc, iSpc, 0);
                        long GlobalOffset = m_map.GlobalUniqueIndex(iVar, jLoc, iSpc, 0);
                        var tmpMod = new List<extNi0>();

                        // loop over polynomial degrees...
                        for (int degree = 0; degree < DGdegreeP1[iVar]; degree++) {
                            if (ModeInstruction(jLoc, iVar, iSpc, degree)) {
                                long GlobalModeOffset = m_Ni0[degree].i0 + GlobalOffset;
                                int LocalModeOffset = m_Ni0[degree].i0 + LocalOffset;
                                int ModeLength = m_Ni0[degree].N;
                                var newNi0 = new extNi0(LocalModeOffset, GlobalModeOffset, SubOffset, ModeLength);
                                SubOffset += ModeLength;

                                // Fill int lists
                                for (int m = 0; m < newNi0.N; m++) {
                                    Globalint.Add(newNi0.Gi0 + m);
                                    Localint.Add(newNi0.Li0 + m);
                                    SubBlockIdx.Add(newNi0.Si0 + m);
                                    MaskLen++;
                                }

                                // Fill Ni0 Lists
                                tmpMod.Add(newNi0);
                                Ni0Length++;
                                ListNi0.Add(newNi0);
                                Debug.Assert(m_map.LocalUniqueIndex(iVar, jLoc, iSpc, GetNp(degree) - 1) == LocalModeOffset + ModeLength - 1);
                            }
                        }
                        if (tmpMod.Count > 0)
                            tmpSpc.Add(tmpMod.ToArray());
                    }
                    if (tmpSpc.Count > 0)
                        tmpVar.Add(tmpSpc.ToArray());
                }
                if (tmpVar.Count > 0)
                    tmpCell.Add(tmpVar.ToArray());
            }
            var tmpStructNi0 = tmpCell.ToArray();

            int NumOfNi0 = 0;
#if DEBUG
            for(int iCell = 0; iCell < tmpStructNi0.Length; iCell++) {
                for(int iVar = 0; iVar < tmpStructNi0[iCell].Length; iVar++) {
                    for(int iSpc = 0; iSpc < tmpStructNi0[iCell][iVar].Length; iSpc++) {
                        NumOfNi0 += tmpStructNi0[iCell][iVar][iSpc].Length;
                    }
                }
            }
#endif
            // an empty selection is allowed,
            // e.g. consider a combination of empty external and non empty local mask
            // In case of IBM, selections with no species are also allowed
            if (!emptysel && totNoOfSpeciesInSelection > 0) {

                // `ListNi0.Count` filters some phatological use case, e.g. very coarse meshes

                Debug.Assert(ListNi0.Count <= 0 || ListNi0.GroupBy(x => x.Li0).Any(g => g.Count() == 1)); // test for uniqueness of local index
                Debug.Assert(ListNi0.Count <= 0 || ListNi0.GroupBy(x => x.Gi0).Any(g => g.Count() == 1)); // test for uniqueness of global index
                Debug.Assert(ListNi0.Count() == NumOfNi0);
                Debug.Assert(MaskLen <= m_LocalLength);
                Debug.Assert(Localint.Count <= 0 || Localint.GroupBy(x => x).Any(g => g.Count() == 1));
                Debug.Assert(Globalint.Count <= 0 || Globalint.GroupBy(x => x).Any(g => g.Count() == 1));
                Debug.Assert(Localint.Count() == MaskLen);
                Debug.Assert(SubBlockIdx.Count() == MaskLen);
                Debug.Assert(SubBlockIdx.Count <= 0 || SubBlockIdx.GroupBy(x => x).Any(g => g.Count() == 1));
            }
            
            m_GlobalMask = Globalint;
            m_LocalMask = Localint;
            m_StructuredNi0 = tmpStructNi0;
            m_MaskLen = MaskLen;
            m_Ni0Len = Ni0Length;
            m_SubBlockMask = SubBlockIdx;
        }

        /// <summary>
        /// Get the local index list of a cell within this mask
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public int[] GetLocalidcOfCell(int iCell) {
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
        /// Get the Global index list of a cell within this mask
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public long[] GetGlobalidcOfCell(int iCell) {
            List<long> cellidx = new List<long>();
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
            long[] array = cellidx.ToArray();
            Debug.Assert(array.GroupBy(x => x).Any(g => g.Count() == 1));
            return array;
        }

        /// <summary>
        /// Get the length of a cell within this mask
        /// </summary>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public int GetLengthOfCell(int iCell) {
            int len = 0;
            Debug.Assert(iCell< m_StructuredNi0.Length);
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

        /// <summary>
        /// Get Block length for cells in mask
        /// </summary>
        /// <returns></returns>
        public List<int> GetAllSubMatrixCellLength() {
            var SubM_N = new List<int>();
            for (int iCell = 0; iCell < m_StructuredNi0.Length; iCell++) {
                SubM_N.Add(GetLengthOfCell(iCell));
            }
            return SubM_N;
        }

        /// <summary>
        /// Gets all offsets of subblock index for cellblocks in this mask.
        /// </summary>
        /// <returns></returns>
        public List<long> GetAllSubMatrixCellOffsets() {
            var intList = new List<long>();
            for (int iCell = 0; iCell < m_StructuredNi0.Length; iCell++)
                intList.Add(m_StructuredNi0[iCell][0][0][0].Si0);
            return intList;
        }

        /// <summary>
        /// gets sub-block offset relative to the parent cell-block
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
        /// Number of cell blocks in mask
        /// </summary>
        public int NoOfCells {
            get {
                return m_StructuredNi0.Length;
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
