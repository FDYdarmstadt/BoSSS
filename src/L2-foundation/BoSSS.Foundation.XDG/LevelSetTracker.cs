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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;
using System.Collections.ObjectModel;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// The level set - tracker manages the tracking of narrow - bands of up to 4 level set.
    /// It is a perquisite for
    /// the memory management of cut-cell fields, i.e. <see cref="XDGField"/> and <see cref="XDGBasis"/>
    /// </summary>
    /// <remarks>
    /// After changing one ore more level sets, the <see cref="UpdateTracker()"/>-method must be
    /// called.
    /// </remarks>
    public partial class LevelSetTracker : IObservable<LevelSetTracker.LevelSetRegions>, IDisposable {

        /// <summary>
        /// Region code (see <see cref="m_LevSetRegions"/>) indicating that all level set values in a cell are in
        /// the positive far region
        /// </summary>
        public const ushort AllFARplus = 0xffff;

        /// <summary>
        /// Region code (see <see cref="m_LevSetRegions"/>) indicating that all level set values in a cell are in
        /// the negative far region
        /// </summary>
        public const ushort AllFARminus = 0x1111;

        /// <summary>
        /// Region code (see <see cref="m_LevSetRegions"/>) indicating that all level sets cut the actual cell 
        /// </summary>
        public const ushort AllCut = 0x8888;

        /// <summary>
        /// Creates a level set tracker for just one level set.
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>
        public LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, string[] _SpeciesTable, ILevelSet levSet1) {
            ConstructorCommon(BackgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, levSet1);
        }

        /// <summary>
        /// Creates a level set tracker for just one level set.
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        /// <param name="BruteForceDivisions"></param>
        /// <param name="BruteForceOrder"></param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>
        public LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, int BruteForceDivisions, int BruteForceOrder, string[] _SpeciesTable, ILevelSet levSet1) {
            ConstructorCommon(BackgroundGrid, __NearRegionWidth, BruteForceOrder, BruteForceDivisions, _SpeciesTable, cutCellquadType, levSet1);
        }

        /// <summary>
        /// Creates a level set tracker for two level sets.
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="levSet2"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>
        public LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, string[,] _SpeciesTable, ILevelSet levSet1, ILevelSet levSet2) {
            ConstructorCommon(BackgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, levSet1, levSet2);
        }

        /// <summary>
        /// Creates a level set tracker for three level sets.
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="levSet2"></param>
        /// <param name="levSet3"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>
        public LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, string[,,] _SpeciesTable, ILevelSet levSet1, ILevelSet levSet2, ILevelSet levSet3) {
            ConstructorCommon(BackgroundGrid, cutCellquadType, __NearRegionWidth, SpeciesTable, levSet1, levSet2, levSet3);
        }

        /// <summary>
        /// Creates a level set tracker for four level sets.
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="levSet2"></param>
        /// <param name="levSet3"></param>
        /// <param name="levSet4"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>
        public LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, string[,,,] _SpeciesTable, ILevelSet levSet1, ILevelSet levSet2, ILevelSet levSet3, ILevelSet levSet4) {
            ConstructorCommon(BackgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, levSet1, levSet2, levSet3, levSet4);
        }

        /// <summary>
        /// Implementation of the constructor with default values for the brute
        /// force quadrature (see <see cref="SetBruteForceQuadratureRules"/>).
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="__NearRegionWidth"></param>
        /// <param name="SpeciesTable"></param>
        /// <param name="levSets"></param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>       
        private LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, Array SpeciesTable, params ILevelSet[] levSets) {
            ConstructorCommon(BackgroundGrid, __NearRegionWidth, 5, 2, SpeciesTable, cutCellquadType, levSets);
        }

        /// <summary>
        /// Implementation of the constructor with default values for the brute
        /// force quadrature (see <see cref="SetBruteForceQuadratureRules"/>).
        /// </summary>
        private void ConstructorCommon(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, Array SpeciesTable, params ILevelSet[] levSets) {
            ConstructorCommon(BackgroundGrid, __NearRegionWidth, 5, 2, SpeciesTable, cutCellquadType, levSets);
        }

        /// <summary>
        /// Increases <see cref="HistoryLength"/> to <paramref name="ReqLeng"/>, if the latter is smaler.
        /// </summary>
        /// <param name="ReqLeng">The requested length</param>
        /// <returns>
        /// the actual value of <see cref="HistoryLength"/>
        /// </returns>
        public int IncreaseHistoryLength(int ReqLeng) {
            if(ReqLeng < 0)
                throw new ArgumentException();
            HistoryLength = Math.Max(ReqLeng, HistoryLength);
            return HistoryLength;
        }


        /// <summary>
        /// Number of previous states in the various history stacks (<see cref="DataHistories"/>, <see cref="LevelSetHistories"/>, <see cref="RegionsHistory"/>, etc.).
        /// </summary>
        public int HistoryLength {
            get {
                int L = RegionsHistory.HistoryLength;

                for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                    Debug.Assert(m_DataHistories[iLs].HistoryLength == L);
                    Debug.Assert(m_LevelSetHistories[iLs].HistoryLength == L);
                }
                Debug.Assert(m_QuadFactoryHelpersHistory.HistoryLength == L);
                Debug.Assert(m_XDGSpaceMetricsHistory.HistoryLength == L);
                return L;
            }
            set {
                int currentLen = this.HistoryLength;
                if(value != currentLen) {
                    RegionsHistory.HistoryLength = value;
                    for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                        m_DataHistories[iLs].HistoryLength = value;
                        m_LevelSetHistories[iLs].HistoryLength = value;
                    }
                    m_QuadFactoryHelpersHistory.HistoryLength = value;
                    m_XDGSpaceMetricsHistory.HistoryLength = value;
                }
                Debug.Assert(HistoryLength == value);
            }
        }

        /// <summary>
        /// Number of times <see cref="PushStacks"/> was called;
        /// </summary>
        public int PushCount {
            get {
                int L = RegionsHistory.PushCount;

                for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                    Debug.Assert(m_DataHistories[iLs].PushCount == L);
                    Debug.Assert(m_LevelSetHistories[iLs].PushCount == L);
                }
                Debug.Assert(m_QuadFactoryHelpersHistory.PushCount == L);
                Debug.Assert(m_XDGSpaceMetricsHistory.PushCount == L);

                return L;
            }
        }

        /// <summary>
        /// The number of previous time-steps available in the various stacks, see also <see cref="HistoryStack{T}.GetPopulatedLength"/>.
        /// </summary>
        public int PopulatedHistoryLength {
            get {
                int L = RegionsHistory.GetPopulatedLength();

                for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                    Debug.Assert(m_DataHistories[iLs].GetPopulatedLength() == L);
                    Debug.Assert(m_LevelSetHistories[iLs].GetPopulatedLength() == L);
                }
                Debug.Assert(m_QuadFactoryHelpersHistory.GetPopulatedLength() == L);
                Debug.Assert(m_XDGSpaceMetricsHistory.GetPopulatedLength() == L);
                
                return L;
            }
        }

        /// <summary>
        /// Implementation of the constructor;
        /// </summary>
        private void ConstructorCommon(GridData griData, int __NearRegionWidth, int BruteForceDivisions, int BruteForceOrder, Array SpeciesTable, XQuadFactoryHelper.MomentFittingVariants cutCellQuadratureType, params ILevelSet[] levSets) {
            // check args, init members
            // ========================
            m_gDat = griData;
            if (__NearRegionWidth < 0 || __NearRegionWidth > 6)
                throw new ArgumentException("near region width must be between 0 (including) and 7 (excluding)", "__NearRegionWidth");
            m_NearRegionWidth = __NearRegionWidth;
            
            if (SpeciesTable.Rank != levSets.Length)
                throw new ArgumentException("rank of species table must match number of level sets.", "SpeciesTable");
            m_SpeciesTable = SpeciesTable;
            m_LevSetAllowedMovement = new int[levSets.Length];
            ArrayTools.SetAll(m_LevSetAllowedMovement, 1);
            this.CutCellQuadratureType = cutCellQuadratureType;

            // init stacks
            // ===========

            var _LevelSetHistories = new HistoryStack<ILevelSet>[levSets.Length];
            var _DataHistories = new HistoryStack<LevelSetData>[levSets.Length];
            for(int iLs =0; iLs < levSets.Length; iLs++) {
                _LevelSetHistories[iLs] = new HistoryStack<ILevelSet>(levSets[iLs]);
                _DataHistories[iLs] = new HistoryStack<LevelSetData>(new LevelSetData(this, iLs));
            }
            m_DataHistories = _DataHistories.ToList().AsReadOnly();
            m_LevelSetHistories = _LevelSetHistories.ToList().AsReadOnly();
            m_RegionsHistory = new HistoryStack<LevelSetRegions>(new LevelSetRegions(this));
            m_QuadFactoryHelpersHistory = new HistoryStack<Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper>>(
                new Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper>());
            m_XDGSpaceMetricsHistory = new HistoryStack<Dictionary<Tuple<SpeciesId[], XQuadFactoryHelper.MomentFittingVariants, int>, XDGSpaceMetrics>>(NewXDGSpaceMetricsCache());

            this.IncreaseHistoryLength(1); // at least one previous time-step is required to support update of XDG fields

            // 1st tracker update
            // ==================

            SpeciesId invalid;
            invalid.cntnt = int.MinValue;
            ArrayTools.SetAll(m_SpeciesIndex2Id, invalid);
            ArrayTools.SetAll(m_SpeciesId2Index, int.MinValue);
            CollectSpeciesNamesRecursive(0, new int[m_SpeciesTable.Rank]);
            ComputeNoOfSpeciesRecursive(0, new int[4]);
            CollectLevelSetSignCodes();
            ComputeGhostTable();


            UpdateTracker();
            PushStacks();
        }

        /// <summary>
        /// The type of integration used in cut cells.
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType {
            get;
            private set;
        }

        /// <summary>
        /// used by <see cref="ContainesSpecies"/>;
        ///  - 1st index: <em>ID</em>: Species Id minus <see cref="___SpeciesIDOffest"/>;
        ///  - 2nd index: <em>SIGN</em>: <see cref="LevelsetCellSignCode.lsSig"/>;
        /// Content: true, if species with given <em>ID</em> is contained in a cell with level set sign code <em>SIGN</em>;
        /// </summary>
        bool[][] m_GhostTable;

        /// <summary>
        /// initializes <see cref="m_GhostTable"/>
        /// </summary>
        void ComputeGhostTable() {
            m_GhostTable = new bool[this.m_SpeciesNames.Count][];
            int NoOfLevSets = this.NoOfLevelSets;

            foreach (String specNmn in this.m_SpeciesNames) {
                SpeciesId specId = GetSpeciesId(specNmn);
                LevelSetSignCode[] speciesCodes = GetLevelSetSignCodes(specNmn);

                bool[] Table = new bool[81];
                m_GhostTable[specId.cntnt - ___SpeciesIDOffest] = Table;

                LevelsetCellSignCode cd;
                for (cd.lsSig = 0; cd.lsSig < 81; cd.lsSig++) {

                    foreach (var scd in speciesCodes) {
                        if (cd.IsContained(scd, NoOfLevSets))
                            Table[cd.lsSig] = true;
                    }
                }
            }
        }

        /// <summary>
        /// Detects whether a species in some cell 
        /// </summary>
        /// <param name="id">
        /// species id;
        /// </param>
        /// <param name="cd">
        /// acquire by <see cref="LevelsetCellSignCode.Extract"/>;
        /// </param>
        /// <returns>
        /// true, if species with given <paramref name="id"/> is contained in a cell with level set sign code 
        /// <paramref name="cd"/>; Not that, if the return value is false, the cell still may have DOF allocated
        /// for the species;
        /// </returns>
        public bool ContainesSpecies(SpeciesId id, LevelsetCellSignCode cd) {
            return m_GhostTable[id.cntnt - ___SpeciesIDOffest][cd.lsSig];
        }



        /// <summary>
        /// initializes <see cref="m_SpeciesNames"/>
        /// </summary>
        /// <param name="recdeph"></param>
        /// <param name="indices"></param>
        void CollectSpeciesNamesRecursive(int recdeph, int[] indices) {
            if (m_SpeciesTable.GetLength(recdeph) != 2
                || m_SpeciesTable.GetLowerBound(recdeph) != 0)
                throw new ArgumentException("length of each dimension in species table must be 2, lower bound must be 0");

            for (indices[recdeph] = 0; indices[recdeph] < 2; indices[recdeph]++) {
                if (recdeph < (m_SpeciesTable.Rank - 1)) {
                    CollectSpeciesNamesRecursive(recdeph + 1, indices);
                } else {
                    string species = (string)m_SpeciesTable.GetValue(indices);
                    if (!m_SpeciesNames.Contains(species))
                        m_SpeciesNames.Add(species);
                }
            }
        }

        /// <summary>
        /// returns a collection of all species (identified by their names) which can 
        /// possibly occur in a cell with <paramref name="ReducedRegionCode"/>
        /// </summary>
        /// <param name="ReducedRegionCode">
        /// the reduced region code (<see cref="ReducedRegionCode.Extract"/>);
        /// </param>
        /// <returns></returns>
        public ICollection<string> CollectPossibleSpecies(ReducedRegionCode ReducedRegionCode) {
            if (ReducedRegionCode.rrc < 0 || ReducedRegionCode.rrc >= 81)
                throw new ArgumentOutOfRangeException("ReducedRegionCode must satisfy: 0 <= ReducedRegionCode < 81;");

            int[] signs = new int[4];
            int b = 3 * 3 * 3;
            for (int i = 3; i >= 0; i--) {
                int sign_i = ReducedRegionCode.rrc / b;
                switch (sign_i) {
                    case 2:
                    signs[i] = 1;
                    break;
                    case 1:
                    signs[i] = -1;
                    break;
                    case 0:
                    signs[i] = 0;
                    break;
                    default:
                    throw new ApplicationException("algorithm error");
                }

                ReducedRegionCode.rrc -= sign_i * b;
                b /= 3;
            }

            List<string> ret = new List<string>();
            CollectPossibleSpeciesRecursive(signs, ret, 0, new int[m_SpeciesTable.Rank]);
            return ret;
        }

        /// <summary>
        /// collects the names of species that can occur in a cell
        /// </summary>
        /// <param name="LevSetSigns"></param>
        /// <returns>
        /// a collection containing the names of all species that can occur in a cell
        /// with level sets with signs <paramref name="LevSetSigns"/>;
        /// </returns>
        public ICollection<string> CollectPossibleSpecies(params LevelsetSign[] LevSetSigns) {
            if (LevSetSigns.Length != NoOfLevelSets)
                throw new ArgumentException("length must match number of level sets;", "LevSetSigns");

            int[] signs = new int[LevSetSigns.Length];
            for (int i = 0; i < LevSetSigns.Length; i++)
                signs[i] = (int)LevSetSigns[i];

            List<string> ret = new List<string>();
            CollectPossibleSpeciesRecursive(signs, ret, 0, new int[m_SpeciesTable.Rank]);
            return ret;
        }


        void CollectPossibleSpeciesRecursive(int[] LevSetSigns, ICollection<string> outp, int r, int[] indices) {
            for (indices[r] = 0; indices[r] < 2; indices[r]++) {
                if (r < (m_SpeciesTable.Rank - 1)) {
                    CollectPossibleSpeciesRecursive(LevSetSigns, outp, r + 1, indices);
                } else {
                    bool maskfits = true;
                    for (int i = 0; i < m_SpeciesTable.Rank; i++) {
                        switch (LevSetSigns[i]) {
                            case 0:
                            break;

                            case 1:
                            // pos. sign;
                            if (indices[i] == 0)
                                maskfits = false;
                            break;

                            case -1:
                            // neg. sign:
                            if (indices[i] == 1)
                                maskfits = false;
                            break;

                            default:
                            throw new ArgumentOutOfRangeException();
                        }
                    }

                    string species = (string)m_SpeciesTable.GetValue(indices);
                    if (maskfits && !outp.Contains(species))
                        outp.Add(species);
                }
            }
        }

        /// <summary>
        /// for each (reduced) region code, this method
        /// computes the possible number of species in that cell (<see cref="m_NoOfSpecies"/>)
        /// and assigned
        /// </summary>
        /// <param name="recdeph"></param>
        /// <param name="ReducedRegionSign"></param>
        void ComputeNoOfSpeciesRecursive(int recdeph, int[] ReducedRegionSign) {
            List<string> speciesList = new List<string>();


            for (ReducedRegionSign[recdeph] = -1; ReducedRegionSign[recdeph] <= 1; ReducedRegionSign[recdeph]++) {
                if (recdeph >= 3) {
                    speciesList.Clear();
                    CollectPossibleSpeciesRecursive(ReducedRegionSign, speciesList, 0, new int[m_SpeciesTable.Rank]);
                    ReducedRegionCode ReducedRegionCode;
                    SetNumSpecies(speciesList.Count, ReducedRegionSign, out ReducedRegionCode);

                    LevelSetSignCode signcode;
                    for (signcode.val = 0; signcode.val < 16; signcode.val++) {
                        string species = GetSpecies(signcode);
                        SpeciesId SpeciesId = GetSpeciesId(species);

                        int i0 = speciesList.IndexOf(species);
                        SetSpeciesIndex(ReducedRegionCode.rrc, signcode, i0, SpeciesId, speciesList.Count);
                    }
                } else {
                    ComputeNoOfSpeciesRecursive(recdeph + 1, ReducedRegionSign);
                }
            }
        }

        List<string> m_SpeciesNames = new List<string>();

        internal const int ___SpeciesIDOffest = 11111;

        /// <summary>
        /// The set of all species names tracked by this object.
        /// </summary>
        public IList<string> SpeciesNames {
            get {
                return m_SpeciesNames.ToArray();
            }
        }

        /// <summary>
        /// The set of all species id's tracked by this object.
        /// </summary>
        public IList<SpeciesId> SpeciesIdS {
            get {
                var ret = new List<SpeciesId>();
                foreach (var sp in m_SpeciesNames) {
                    ret.Add(GetSpeciesId(sp));
                }
                return ret;
            }
        }

        /// <summary>
        /// the number of all different species;
        /// </summary>
        public int TotalNoOfSpecies {
            get {
                return m_SpeciesNames.Count;
            }
        }


        /// <summary>
        /// returns the name of the species with ID-number <paramref name="SpeciesId"/>;
        /// </summary>
        /// <param name="SpeciesId"></param>
        /// <returns></returns>
        public string GetSpeciesName(SpeciesId SpeciesId) {
            return m_SpeciesNames[SpeciesId.cntnt - ___SpeciesIDOffest];
        }

        /// <summary>
        /// returns the ID-number which the level set tracker (this object) 
        /// has assigned to the species with name <paramref name="species"/>.
        /// </summary>
        /// <param name="species"></param>
        /// <returns>species ID-number</returns>
        public SpeciesId GetSpeciesId(string species) {
            int r = m_SpeciesNames.IndexOf(species);
            if (r < 0)
                throw new ArgumentException("unknown species name");
            SpeciesId ret;
            ret.cntnt = r + ___SpeciesIDOffest;
            return ret;
        }

        /// <summary>
        /// see <see cref="SpeciesTable"/>
        /// </summary>
        Array m_SpeciesTable;

        /// <summary>
        /// Gets a copy of the species table, a multidimensional string-array.<br/>
        /// The species table defines, which species corresponds with 
        /// which combination of level set - signs.
        /// </summary>
        /// <remarks>
        /// The rank (number of dimensions) of this array is equal to the number of 
        /// level-sets (<see cref="LevelSets"/>).
        /// Valid indices are only 0 and 1, 0 corresponds to a negative sign an 1 to the positive sign.
        /// E.g. if there are two level sets <em>G</em><sub>1</sub> and <em>G</em><sub>1</sub>
        /// and for some point in space, their sign is -1 and 1, than the species 
        /// at this point is identified by entry [0,1].
        /// </remarks>
        public Array SpeciesTable {
            get {
                return (Array)m_SpeciesTable.Clone();
            }
        }

        HistoryStack<LevelSetRegions> m_RegionsHistory;

        /// <summary>
        /// Information about the level set sign for current and previous states (see <see cref="PushStacks"/>).
        /// </summary>
        public HistoryStack<LevelSetRegions> RegionsHistory {
            get {
                return m_RegionsHistory;
            }
        }

        /// <summary>
        /// Information about the current level set sign in each cell.
        /// </summary>
        public LevelSetRegions Regions {
            get {
                return RegionsHistory.Current;
            }
        }

        

        /// <summary>
        /// 
        /// </summary>
        /// <param name="LevSetSignCode"></param>
        /// <returns></returns>
        public string GetSpecies(LevelSetSignCode LevSetSignCode) {
            int[] indices = new int[m_SpeciesTable.Rank];
            if (indices.Length >= 4 && ((LevSetSignCode.val & 0x8) != 0))
                indices[3] = 1;
            if (indices.Length >= 3 && ((LevSetSignCode.val & 0x4) != 0))
                indices[2] = 1;
            if (indices.Length >= 2 && ((LevSetSignCode.val & 0x2) != 0))
                indices[1] = 1;
            if (indices.Length >= 1 && ((LevSetSignCode.val & 0x1) != 0))
                indices[0] = 1;

            return (string)m_SpeciesTable.GetValue(indices);
        }

        /// <summary>
        /// computes the id of the species that is set for a specific sign
        /// </summary>
        /// <param name="levelSetValues">signs of all level set functions</param>
        /// <returns></returns>
        public SpeciesId GetSpeciesIdFromSign(params double[] levelSetValues) {
            if (levelSetValues.Length != NoOfLevelSets)
                throw new ArgumentException();
            var lssc = LevelSetSignCode.ComputeLevelSetBytecode(levelSetValues);
            return this.GetSpeciesId(this.GetSpecies(lssc));
        }

        /// <summary>
        /// Determines all sign codes (i.e. level set sign combinations, see
        /// <see cref="GetSpecies(LevelSetSignCode)"/>) in <see cref="SpeciesTable"/> that
        /// identify a given species with name <paramref name="species"/>
        /// <example>
        /// Consider two level sets "Phi_0" and "Phi_1" (with level set indices 0 and
        /// 1, respectively) and three species "a", "b" and "c". Consider the
        /// following <see cref="SpeciesTable"/>:
        /// [0, 0] = "a"
        /// [0, 1] = "b"
        /// [1, 0] = "c"
        /// [1, 1] = "c"
        /// That is, if "Phi_0" and "Phi_1" are both negative, we have species "a". If
        /// "Phi_0" is negative and "Phi_1" is positive, we have "b". If "Phi_0" is
        /// positive, we have "c", irrespective of the value of "Phi_1". Thus,
        /// there are two sign combinations for "c" (namely [+,-] and [+,+]).
        /// This translates to the sign codes 1 and 3 and the return value
        /// would be [1, 3] in such a case.
        /// </example>
        /// </summary>
        /// <param name="species">
        /// The considered species
        /// </param>
        /// <returns>
        /// A list of sign codes identifying <paramref name="species"/>
        /// </returns>
        public LevelSetSignCode[] GetLevelSetSignCodes(string species) {
            return GetLevelSetSignCodes(GetSpeciesId(species));
        }

        /// <summary>
        /// see <see cref="GetLevelSetSignCodes(string)"/>
        /// </summary>
        public LevelSetSignCode[] GetLevelSetSignCodes(SpeciesId species) {
            int idx = species.cntnt - ___SpeciesIDOffest;
            return m_LevelSetSignCodes[idx];
        }


        void CollectLevelSetSignCodes() {
            m_LevelSetSignCodes = new LevelSetSignCode[this.SpeciesNames.Count][];

            foreach (var species in this.SpeciesNames) {
                var SpeciesId = GetSpeciesId(species);

                int idx = SpeciesId.cntnt - ___SpeciesIDOffest;

                List<LevelSetSignCode> result = new List<LevelSetSignCode>();

                int[] indices = new int[m_SpeciesTable.Rank];
                int maxNoOfSpecies = 0x1 << m_SpeciesTable.Rank;
                for (int i = 0; i < maxNoOfSpecies; i++) {
                    for (int j = 0; j < m_SpeciesTable.Rank; j++) {
                        indices[j] = ((0x1 << j) & i) >> j;
                    }

                    if (species.Equals(m_SpeciesTable.GetValue(indices))) {
                        LevelSetSignCode cd;
                        cd.val = i;
                        result.Add(cd);
                    }
                }

                m_LevelSetSignCodes[idx] = result.ToArray();

            }

        }

        LevelSetSignCode[][] m_LevelSetSignCodes;


        int m_NearRegionWidth = 1;

        /// <summary>
        /// Width, in Number of cells, of the near field (set as an argument of <see cref="UpdateTracker()"/>);
        /// </summary>
        public int NearRegionWidth {
            get {
                return m_NearRegionWidth;
            }
            //set {
            //    if (value < 0 || value > 6)
            //        throw new ArgumentException("maximum width for near region is 4 cells.");
            //    m_NearRegionWidth = value;
            //}
        }



        /// <summary>
        /// The level sets stack;
        /// - 1st index: index of level-set
        /// - 2nd index: index into stack
        /// </summary>
        public IList<HistoryStack<ILevelSet>> LevelSetHistories {
            get {
                return m_LevelSetHistories.ToList();
            }
        }

        /// <summary>
        /// Advances the history stacks which store previous states of the level-set tracker
        /// (see <see cref="DataHistories"/>, <see cref="LevelSetHistories"/>, <see cref="RegionsHistory"/>, etc.).
        /// </summary>
        public void PushStacks() {
            int NoOfLs = LevelSets.Count;
            int PHL = PopulatedHistoryLength;
            int HL = this.HistoryLength;

            Debug.Assert(NoOfLs == m_LevelSets.Count);
            for(int iLs = 0; iLs < NoOfLs; iLs++) {
                m_LevelSetHistories[iLs].Push((ls1) => ls1, (ls1, ls0) => ls1.CloneAs());
            }

            Debug.Assert(NoOfLs == m_DataHistories.Count);
            for(int iLs = 0; iLs < NoOfLs; iLs++) {
                m_DataHistories[iLs].Push((data1) => new LevelSetData(this, iLs), (data1, data0) => data1);

                // fix the history index...
                for(int iStack = 1; iStack > Math.Max(-PHL - 1, -HL); iStack--) {
                    m_DataHistories[iLs][iStack].m_HistoryIndex = iStack;
                }
            }

            m_RegionsHistory.Push((r1) => r1.CloneAs(), (r1, r0) => r1);

            m_QuadFactoryHelpersHistory.Push((r1) => new Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper>(), (r1, r0) => r1);

            m_XDGSpaceMetricsHistory.Push((r1) => NewXDGSpaceMetricsCache(), (r1,r0) => r1);

#if DEBUG
            for(int iLs = 0; iLs < NoOfLs; iLs++) {
                Debug.Assert(object.ReferenceEquals(LevelSets[iLs], LevelSetHistories[iLs].Current));
            }
            PHL = PopulatedHistoryLength;
            for(int iH = 1; iH > -PHL + 1; iH--) {
                for(int iLs = 0; iLs < NoOfLs; iLs++) {
                    Debug.Assert(!object.ReferenceEquals(LevelSetHistories[iLs][iH], LevelSetHistories[iLs][iH - 1]));
                }
                Debug.Assert(!object.ReferenceEquals(RegionsHistory[iH], RegionsHistory[iH - 1]));
            }

            for(int iH = 1; iH > -PHL; iH--) {
                for(int iLs = 0; iLs < NoOfLs; iLs++) {
                    Debug.Assert(DataHistories[iLs][iH].HistoryIndex == iH);
                }
            }
#endif
        }

        /// <summary>
        /// Number of used level-sets fields, between 1 and 4.
        /// </summary>
        public int NoOfLevelSets {
            get {
                int L = m_LevelSetHistories.Count;
                Debug.Assert(L == LevelSets.Count);
                Debug.Assert(L == LevelSetHistories.Count);
                Debug.Assert(L == DataHistories.Count);
                return L;
            }
        }


        ReadOnlyCollection<ILevelSet> m_LevelSets = null;

        /// <summary>
        /// The level sets, identic to the top of the level-set stack 
        /// - 1st index: index of level-set
        /// - 2nd index: index into stack
        /// </summary>
        public IList<ILevelSet> LevelSets {
            get {
                if(m_LevelSets == null) { // accelerate access if a user accesses this very often
                    var _LevelSets = new ILevelSet[m_LevelSetHistories.Count];
                    for(int iLs = 0; iLs < m_LevelSetHistories.Count; iLs++) {
                        _LevelSets[iLs] = m_LevelSetHistories[iLs].Current;
                    }
                    m_LevelSets = _LevelSets.ToList().AsReadOnly();
                }
#if DEBUG
                Debug.Assert(m_LevelSets.Count == m_LevelSetHistories.Count);
                for(int iLs = 0; iLs < m_LevelSetHistories.Count; iLs++) {
                    Debug.Assert(object.ReferenceEquals(m_LevelSets[iLs], m_LevelSetHistories[iLs].Current));
                }
#endif
                return m_LevelSets;
            }
        }


        GridData m_gDat;

        /// <summary>
        /// Reference to the background grid of the XDG space.
        /// </summary>
        public GridData GridDat {
            get {
                return m_gDat;
            }
        }

        ReadOnlyCollection<HistoryStack<ILevelSet>> m_LevelSetHistories;
       

        /// <summary>
        /// returns the (possible) number of species in 
        /// a cell with region code <paramref name="RegionCode"/>.
        /// </summary>
        /// <param name="RegionCode">
        /// region code for the cell;
        /// </param>
        /// <param name="_ReducedRegionCode">
        /// on exit, the reduced region code for <paramref name="RegionCode"/>
        /// in an 3-adic representation; This number 
        /// is later on required as an input for
        /// <see cref="GetSpeciesIndex(SpeciesId,int)"/>;
        /// The reduced region code in 3-adic representation;
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// Here, three states (the reduced region) for each of the for level sets are considered:
        /// <list type="bullet">
        ///   <item>positive far (FAR+)</item>
        ///   <item>negative far (FAR-)</item>
        ///   <item>positive near, negative near or cut</item>
        /// </list>
        /// This implies, that also for cells in the near region, memory is allocated for more than one
        /// species.
        /// </remarks>
        public int GetNoOfSpeciesByRegionCode(ushort RegionCode, out ReducedRegionCode _ReducedRegionCode) {
            _ReducedRegionCode = ReducedRegionCode.Extract(RegionCode);
            return m_NoOfSpecies[_ReducedRegionCode.rrc];
        }

        /// <summary>
        /// Number of species, for each reduced region code;
        /// </summary>
        /// <remarks>
        /// This array, which contains exactly 3<sup>4</sup> = 81 entries, stores the number of
        /// species for each possible combination of level sets; <br/>
        /// For each level set, in each cell, we distinct between three different states:
        /// <list type="bullet">
        ///   <item>the cell if in the positive FAR region: <i>distance</i> == 7, i.e. the <i>code</i> is 0xf</item>
        ///   <item>the cell if in the negative FAR region: <i>distance</i> == -1, i.e the <i>code</i> is 0x1</item>
        ///   <item>the cell is cutted or in the near reagion: -6 &lt; <i>distance</i> &lt; 6</item>
        /// </list>
        /// So, there are in total 3*4 states for some cell (for 4 level sets);
        /// </remarks>
        /// <see cref="GetNoOfSpecies"/>
        int[] m_NoOfSpecies = new int[81];

        /// <summary>
        /// - There are 3 states for each level-set in each cell: completely positive, completely negative or cut;
        /// - Considering at max. 4 level-sets, this leads to 3<sup>4</sup> = 81 different states.
        /// - For the sign of all 4 level-sets, there are 2<sup>4</sup> = 16 different states.
        /// </summary>
        int[] m_SpeciesIndex = new int[81 * 16];


        SpeciesId[] m_SpeciesIndex2Id = new SpeciesId[81 * 16];


        int[] m_SpeciesId2Index = new int[81 * 16];


        /// <summary>
        /// 
        /// </summary>
        /// <param name="NumSpecies"></param>
        /// <param name="ReducedRegionSign">
        /// Signs of the of reduced region, for each level set, i.e.
        /// <list type="bullet">
        ///   <item><b>positive</b>: the cell is in the positive FAR region</item>
        ///   <item><b>negative</b>: the cell is in the negative FAR region</item>
        ///   <item><b>0</b>: the cell is cut or in the positive or negative near region;</item>
        /// </list>
        /// </param>
        /// <param name="_ReducedRegionCode">
        /// on exit, the reduced region code for <paramref name="ReducedRegionSign"/>
        /// </param>
        private void SetNumSpecies(int NumSpecies, int[] ReducedRegionSign, out ReducedRegionCode _ReducedRegionCode) {
            int ReducedRegionCode = 0;

            if (ReducedRegionSign[3] > 0)
                ReducedRegionCode += 2;  // reduced region sign for level set 4 is +FAR => reduced region is 2
            if (ReducedRegionSign[3] < 0)
                ReducedRegionCode += 1;  // reduced region sign for level set 4 is -FAR => reduced region is 1
            ReducedRegionCode *= 3;

            if (ReducedRegionSign[2] > 0)
                ReducedRegionCode += 2;  // reduced region sign for level set 3 is +FAR => reduced region is 2
            if (ReducedRegionSign[2] < 0)
                ReducedRegionCode += 1;  // reduced region sign for level set 3 is -FAR => reduced region is 1
            ReducedRegionCode *= 3;

            if (ReducedRegionSign[1] > 0)
                ReducedRegionCode += 2;  // reduced region sign for level set 2 is +FAR => reduced region is 2
            if (ReducedRegionSign[1] < 0)
                ReducedRegionCode += 1;  // reduced region sign for level set 2 is -FAR => reduced region is 1
            ReducedRegionCode *= 3;

            if (ReducedRegionSign[0] > 0)
                ReducedRegionCode += 2;  // reduced region sign for level set 1 is +FAR => reduced region is 2
            if (ReducedRegionSign[0] < 0)
                ReducedRegionCode += 1;  // reduced region sign for level set 1 is -FAR => reduced region is 1

            _ReducedRegionCode.rrc = ReducedRegionCode;
            m_NoOfSpecies[ReducedRegionCode] = NumSpecies;
        }

        

        /// <summary>
        /// the index of some species, with level set signs <paramref name="LevelSetSignBytecode"/>
        /// and the region <paramref name="rrc"/>.
        /// </summary>
        /// <param name="rrc">
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpecies"/>;
        /// </param>
        /// <param name="LevelSetSignBytecode"></param>
        /// <returns></returns>
        public int GetSpeciesIndex(ReducedRegionCode rrc, LevelSetSignCode LevelSetSignBytecode) {
            int ReducedRegionCode = rrc.rrc << 4;    // equal to int *= 16
            ReducedRegionCode = ReducedRegionCode | LevelSetSignBytecode.val; // equal to ind += LevelSetSignBytecode, if  0 <= LevelSetSignBytecode < 16
            return m_SpeciesIndex[ReducedRegionCode];
        }

        /// <summary>
        /// this function is the inverse to <see cref="GetSpeciesIndex(SpeciesId,int)"/>;
        /// </summary>
        /// <param name="_ReducedRegionCode">
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpecies"/>;
        /// </param>
        /// <param name="SpeciesIndex"></param>
        /// <returns></returns>
        /// <remarks>
        /// this function is the inverse to <see cref="GetSpeciesIdFromIndex"/>
        /// </remarks>
        public SpeciesId GetSpeciesIdFromIndex(ReducedRegionCode _ReducedRegionCode, int SpeciesIndex) {
            int ReducedRegionCode = _ReducedRegionCode.rrc << 4;    // equal to int *= 16
            ReducedRegionCode = ReducedRegionCode | SpeciesIndex; // equal to ind += LevelSetSignBytecode, if  0 <= LevelSetSignBytecode < 16
            return m_SpeciesIndex2Id[ReducedRegionCode];
        }
        
        /// <summary>
        /// converts a species id (<paramref name="_SpeciesId"/>) into a species index;
        /// Note that the sorting of species is not constant (over the grid), so the return value
        /// depends on the reduced region code <paramref name="_ReducedRegionCode"/>;
        /// </summary>
        /// <param name="_ReducedRegionCode">
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpecies"/>;
        /// </param>
        /// <param name="_SpeciesId"></param>
        /// <returns>
        /// a negative value indicates that the species is not present within the given <paramref name="_ReducedRegionCode"/>;
        /// </returns>
        /// <remarks>
        /// this function is the inverse to <see cref="GetSpeciesIdFromIndex"/>
        /// </remarks>
        public int GetSpeciesIndex(ReducedRegionCode _ReducedRegionCode, SpeciesId _SpeciesId) {
            int SpeciesId = _SpeciesId.cntnt - ___SpeciesIDOffest;
            if (SpeciesId < 0 || SpeciesId >= 16)
                throw new ArgumentOutOfRangeException("illegal ID");
            int __ReducedRegionCode = _ReducedRegionCode.rrc << 4;
            __ReducedRegionCode |= SpeciesId;
            return m_SpeciesId2Index[__ReducedRegionCode];
        }

        

        private void SetSpeciesIndex(int ind, LevelSetSignCode LevelSetSignBytecode, int SpeciesIndex, SpeciesId speciesId, int NoOfSpecies) {
            int ind2 = ind << 4;
            ind2 = ind2 | LevelSetSignBytecode.val;
            m_SpeciesIndex[ind2] = SpeciesIndex;

            if (SpeciesIndex >= 0 && SpeciesIndex < NoOfSpecies)
                m_SpeciesIndex2Id[ind * 16 + SpeciesIndex] = speciesId;
            m_SpeciesId2Index[ind * 16 + speciesId.cntnt - ___SpeciesIDOffest] = SpeciesIndex;
        }

        ///// <summary>
        ///// 
        ///// </summary>
        //public int GetLevelSetIndex(ILevelSet levSet) {
        //    return Array.IndexOf<ILevelSet>(m_LevelSets, levSet);
        //}

        /// <summary>
        /// 
        /// </summary>
        /// <param name="levSetInd"></param>
        /// <returns></returns>
        public ushort LevelSetBitmask(int levSetInd) {
            if (levSetInd < 0 || levSetInd >= NoOfLevelSets)
                throw new ArgumentException("unknown level set.", "levset");

            return (ushort)(0xF << (4 * levSetInd));
        }


        /// <summary>
        /// Extracts the distance layer index for level-set <paramref name="levSetInd"/>
        /// from the code <paramref name="code"/> (see <see cref="LevelSetRegionsInfo.LevelSetRegions"/>).
        /// </summary>
        /// <param name="code"></param>
        /// <param name="levSetInd"></param>
        /// <returns></returns>
        /// <remarks>
        /// For some cell, the level set distance is 0 for all cells
        /// which are cut by the level set.
        /// </remarks>
        public static int DecodeLevelSetDist(ushort code, int levSetInd) {
            int R = (((int)code >> (4 * levSetInd)) & 0xf) - 8;
            Debug.Assert(-7 <= R);
            Debug.Assert(R <= +7);
            return R;
        }

        /// <summary>
        /// Encodes the distance layer index <paramref name="dist"/> for level-set <paramref name="levSetInd"/>
        /// int the code <paramref name="code"/> (<see cref="LevelSetRegionsInfo.LevelSetRegions"/>).
        /// </summary>
        static public void EncodeLevelSetDist(ref ushort code, int dist, int levSetInd) {
            Debug.Assert(-7 <= dist);
            Debug.Assert(dist <= +7);
            int c = ((dist + 8) << (4 * levSetInd));
            int msk = (0xf << (4 * levSetInd));
            code = (ushort)((code & ~msk) | c);
        }
        
             



        /// <summary>
        /// index: reference element index
        /// </summary>
        NodeSet[] TestNodes;

        /// <summary>
        /// 1st index: reference element index, correlates with 1st index of <see cref="TestNodes"/><br/>
        /// 2nd index: face index
        /// </summary>
        int[][] TestNodesPerFace;

        /// <summary>
        /// a counter that is increased every time when <see cref="UpdateTracker()"/>
        /// is called;
        /// </summary>
        public int VersionCnt {
            get {
                return m_VersionCnt;
            }
        }

        int m_VersionCnt = 0;

        /// <summary>
        /// debug/test code
        /// </summary>
        public void TestRegions() {


            for (int i = 0; i < NoOfLevelSets; i++) { // loop over level sets
                var data = this.m_DataHistories[i].Current;
                for (int j = 0; j < m_gDat.Cells.NoOfCells; j++) {
                    int iKref = this.GridDat.Cells.GetRefElementIndex(j);

                    MultidimensionalArray levSetVals = data.GetLevSetValues(TestNodes[iKref], j, 1);

                    bool containsPos = false;
                    bool containsNeg = false;

                    for (int k = 0; k < levSetVals.GetLength(1); k++) {
                        double v = levSetVals[0, k];

                        if (v <= 0)
                            containsNeg = true;
                        if (v >= 0)
                            containsPos = true;
                    }


                    int dist = LevelSetTracker.DecodeLevelSetDist(Regions.m_LevSetRegions[j], i);

                    if (dist > 0 && containsNeg)
                        throw new ApplicationException();
                    if (dist < 0 && containsPos)
                        throw new ApplicationException();
                    if (dist == 0 && !(containsNeg && containsPos))
                        throw new ApplicationException();
                }
            }

        }

        /// <summary>
        /// checks if the level set may has moved more than one cell
        /// </summary>
        int CheckLevelSetCFL(int LevSetIdx) {
            ushort[] oldCode = this.RegionsHistory[0].m_LevSetRegions;

            CellMask newCut = this.RegionsHistory[1].GetCutCellSubgrid4LevSet(LevSetIdx).VolumeMask;

            int fail_count = 0;

            // for all cells that are cut by the levelset,
            // check whether they are in Near - region of the previous state;
            foreach (var chunk in newCut) {
                for (int i = 0; i < chunk.Len; i++) {
                    int j = i + chunk.i0;

                    int old_dist = LevelSetTracker.DecodeLevelSetDist(oldCode[j], LevSetIdx);
                    if (Math.Abs(old_dist) > 1)
                        fail_count++;
                }
            }

            int failCountGlobal = int.MaxValue;
            unsafe
            {
                MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&fail_count), (IntPtr)(&failCountGlobal), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.INT, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
            }

            return failCountGlobal;
        }

        

        /// <summary>
        /// Early stage of dynamic load balancing, backup of data before the grid cells are re-distributed.
        /// </summary>
        public int[] BackupBeforeLoadBalance() {
            throw new NotImplementedException("todo"); 
            /*
            using (new FuncTrace()) {
                int Jup = this.GridDat.CellPartitioning.LocalLength;
                int[] tempRegionData = new int[Jup];

                var oldRegion = this.RegionsHistory[0].m_LevSetRegions;
                var newRegion = this.RegionsHistory[1].m_LevSetRegions;

                for (int j = 0; j < Jup; j++) {
                    ushort a = oldRegion[j];
                    ushort b = newRegion[j];
                    int C = a | ((int)b << 16);
                    tempRegionData[j] = C;
                }
                
#if DEBUG
                RestoreAfterLoadBalance(this.m_VersionCnt, tempRegionData, false);
                int JA = m_gDat.Cells.NoOfCells;
                var restored1 = RegionsHistory[0].m_LevSetRegions;
                var restored2 = Regions.m_LevSetRegions;
                for (int j = 0; j < JA; j++) {
                    Debug.Assert(restored1[j] == oldRegion[j]);
                    Debug.Assert(restored2[j] == newRegion[j]);
                }
#endif 

                return tempRegionData;
            }
            */
        }

        /// <summary>
        /// Late stage of dynamic load balancing, restoring data after the grid cells were re-distributed.
        /// </summary>
        public void RestoreAfterLoadBalance(int versionNo, int[] ExchangeData, bool ObUp = true) {
            throw new NotImplementedException("todo"); 
            /*
            using(new FuncTrace()) {
                int J = m_gDat.Cells.NoOfLocalUpdatedCells;
                int JA = m_gDat.Cells.NoOfCells;
                Debug.Assert(ExchangeData.Length == J);

                PreviousRegions = new LevelSetRegions(this) {
                    m_LevSetRegions = new ushort[JA],
                    Version = versionNo - 1
                };
                _Regions = new LevelSetRegions(this) {
                    m_LevSetRegions = new ushort[JA],
                    Version = versionNo
                };
                ushort[] oldR = PreviousRegions.m_LevSetRegions;
                ushort[] newR = _Regions.m_LevSetRegions;

                for(int j = 0; j < J; j++) {
                    int C = ExchangeData[j];
                    ushort ra = (ushort)(C & 0x0000FFFF);
                    ushort rb = (ushort)((C & 0xFFFF0000) >> 16);
                    oldR[j] = ra;
                    newR[j] = rb;
                }

                MPIUpdate(oldR, this.GridDat);
                MPIUpdate(newR, this.GridDat);

                PreviousRegions.Recalc_LenToNextchange();
                _Regions.Recalc_LenToNextchange();

                if(ObUp)
                    this.ObserverUpdate();
            }
            */
        }

        /// <summary>
        /// Late stage of mesh adaptation (which may includes dynamic load balancing), 
        /// restoring data after the grid cells were re-distributed.
        /// </summary>
        public void RestoreAfterMeshAdaptation(int versionNo, int[][] ExchangeData, bool ObUp = true) {
            throw new NotImplementedException("todo"); 
            /*
            using(new FuncTrace()) {
                int J = m_gDat.Cells.NoOfLocalUpdatedCells;
                int JA = m_gDat.Cells.NoOfCells;
                Debug.Assert(ExchangeData.Length == J);

                int NoOfLevSets = this.LevelSets.Count;

                PreviousRegions = new LevelSetRegions(this) {
                    m_LevSetRegions = new ushort[JA],
                    Version = versionNo - 1
                };
                Regions = new LevelSetRegions(this) {
                    m_LevSetRegions = new ushort[JA],
                    Version = versionNo
                };
                ushort[] oldR = PreviousRegions.m_LevSetRegions;
                ushort[] newR = _Regions.m_LevSetRegions;

                for(int j = 0; j < J; j++) {
                    int[] CL = ExchangeData[j];
                    int L = CL.Length;
                    ushort ra = 0, rb = 0;
                    if(L > 1) {
                        // +++++++++++++++++++++++
                        // cell coarsening
                        // +++++++++++++++++++++++

                        // combine cell codes (by taking the minimum)
                        
                        for(int iLs = 0; iLs < NoOfLevSets; iLs++) {
                            int NewDistLs = int.MaxValue;
                            int OldDistLs = int.MaxValue;

                            for(int l = 0; l < L; l++) { // loop over all cells in coarsening cell cluster
                                int C = CL[l];
                                ushort ra_l = (ushort)(C & 0x0000FFFF);
                                ushort rb_l = (ushort)((C & 0xFFFF0000) >> 16);

                                int NewLsDist_l = DecodeLevelSetDist(ra_l, iLs);
                                if(Math.Abs(NewDistLs) > Math.Abs(NewLsDist_l)) {
                                    NewDistLs = NewLsDist_l;
                                }

                                int OldLsDist_l = DecodeLevelSetDist(rb_l, iLs);
                                if(Math.Abs(OldDistLs) > Math.Abs(OldLsDist_l)) {
                                    OldDistLs = OldLsDist_l;
                                }
                            }
                            Debug.Assert(Math.Abs(NewDistLs) <= 7);
                            Debug.Assert(Math.Abs(OldDistLs) <= 7);
                            
                            EncodeLevelSetDist(ref ra, NewDistLs, iLs);
                            EncodeLevelSetDist(ref rb, OldDistLs, iLs);
                        }

                    } else {
                        // ++++++++++++++++++++++++++++++++++++++++
                        // cell is either refined or conserved 
                        // ++++++++++++++++++++++++++++++++++++++++

                        int C = CL[0];
                        ra = (ushort)(C & 0x0000FFFF);
                        rb = (ushort)((C & 0xFFFF0000) >> 16);
                    }

                    oldR[j] = ra;
                    newR[j] = rb;
                }

                MPIUpdate(oldR, this.GridDat);
                MPIUpdate(newR, this.GridDat);

                PreviousRegions.Recalc_LenToNextchange();
                Regions.Recalc_LenToNextchange();

                if(ObUp)
                    this.ObserverUpdate();
            }
            */
        }



        /// <summary>
        /// Must be called after changing the level-set;
        /// Invoking this method 
        /// updates the members <see cref="m_LevSetRegions"/>, <see cref="m_LevSetRegionsUnsigned"/>
        /// and <see cref="m_LenToNextChange"/>; <br/>
        /// After this update, for every <see cref="XDGField"/> the method
        /// <see cref="XDGField.UpdateMemory"/> must be invoked;
        /// </summary>
        /// <param name="__NearRegionWith">
        /// new width of near region;
        /// </param>
        /// <param name="__LevSetAllowedMovement">
        /// new values for the allowed level-set movement, see <see cref="GetLevSetAllowedMovement"/>;
        /// If this value is set to a number higher than the <see cref="NearRegionWidth"/>, the CFL test for the 
        /// corresponding 
        /// level set is omitted. 
        /// </param>
        /// <param name="incremental">
        /// If true, the distance of near-cells can only increase by one.
        /// E.g., if the level-set leaves the domain, setting this parameter to true ensures that the respective 
        /// boundary cells are treated as near-cells for one 
        /// <see cref="UpdateTracker(int, bool, int[])"/>-cycle.
        /// This maintains a 'correct' near field when level-set leaves domain, which is necessary e.g. for operator matrix update.
        /// 
        /// Furthermore, it speeds up the detection of cut cells: the detection of cut cells is limited to the near-region,
        /// which drastically reduces the amount of level-set evaluations.
        /// </param>
        public void UpdateTracker(int __NearRegionWith = -1, bool incremental = false, params int[] __LevSetAllowedMovement) {
            using (var tr = new FuncTrace()) {
                ilPSP.MPICollectiveWatchDog.Watch();

                if (this.NearRegionWidth <= 0 && incremental == true) {
                    throw new NotSupportedException("Incremental update requires a near-region width of at least 1.");
                }
                CellMask oldNearMask = null;
                if (incremental)
                    oldNearMask = this.RegionsHistory[1].GetNearFieldMask(m_NearRegionWidth); // index '[1]' is NO mistake!
                else
                    oldNearMask = null;


                // set default values
                if (__NearRegionWith < 0)
                    __NearRegionWith = this.m_NearRegionWidth;
                if (__LevSetAllowedMovement == null || __LevSetAllowedMovement.Length <= 0)
                    __LevSetAllowedMovement = this.m_LevSetAllowedMovement;

                // check args
                if (__LevSetAllowedMovement.Length != m_LevSetAllowedMovement.Length)
                    throw new ArgumentException("length must be equal to number of level sets.", "__LevSetAllowedMovement");
                for (int i = 0; i < __LevSetAllowedMovement.Length; i++)
                    if (__LevSetAllowedMovement[i] < 0)
                        throw new ArgumentOutOfRangeException("__LevSetAllowedMovement", " each entry must be grater or equal to zero.");
                if (__NearRegionWith < 0 || __NearRegionWith > 6)
                    throw new ArgumentException("maximum width for near region is 4 cells.", "NearRegionWith");
                m_NearRegionWidth = __NearRegionWith;

                m_VersionCnt++;
                int J = m_gDat.Cells.NoOfLocalUpdatedCells;
                int JA = m_gDat.Cells.NoOfCells;
                int D = m_gDat.Grid.SpatialDimension;
                var smplx = m_gDat.Grid.RefElements;
                //int NoOfSmplxVertice = smplx.NoOfVertices;
                var Krefs = m_gDat.Grid.RefElements;
                int[] NoOfSmplxVertice = Krefs.Select(Kref => Kref.NoOfVertices).ToArray();
                int NoOfLevSets = this.NoOfLevelSets;
                
                Regions.Version = m_VersionCnt;

                ushort[] VertexMarker, LevSetRegions, LevSetRegionsUnsigned;
                BitArray[] LevSetNeg;
                ushort[,] VertexMarkerExternal;
                #region UpdateTracker_INIT
                using (new BlockTrace("INIT", tr)) {

                    // init memory
                    // ===========
                    Regions.m_LevSetRegions = new ushort[JA];
                    LevSetRegions = Regions.m_LevSetRegions;
                    Regions.InvalidateCaches();
                    LevSetRegionsUnsigned = new ushort[JA];


                    // necessary to avoid a 'LevelSetCFLException' on the first
                    // call to this method
                    for (int i = 0; i < JA; i++) {
                        LevSetRegions[i] = AllCut;
                    }

                    // initialize everything to FAR = +7
                    for (int i = 0; i < JA; i++)
                        LevSetRegionsUnsigned[i] = AllFARplus;

                    VertexMarker = new ushort[m_gDat.Vertices.Count];
                    for (int i = VertexMarker.Length - 1; i >= 0; i--)
                        VertexMarker[i] = AllFARplus;
                    VertexMarkerExternal = new ushort[JA - J, Krefs.Max(Kref => Kref.NoOfVertices)];

                    LevSetNeg = new BitArray[NoOfLevSets]; // true marks a cell in which the level set field is completely negative
                    for (int i = 0; i < LevSetNeg.Length; i++)
                        LevSetNeg[i] = new BitArray(J);

                    // define test vertices
                    // ====================


                    if (TestNodes == null) {
                        TestNodes = new NodeSet[Krefs.Length];
                        TestNodesPerFace = new int[Krefs.Length][];

                        for (int iKref = 0; iKref < Krefs.Length; iKref++) {
                            var Kref = Krefs[iKref];

                            QuadRule rule = Kref.FaceRefElement.GetBruteForceQuadRule(4, 2);

                            int myD = Math.Max(GridDat.Grid.GetRefElement(0).FaceRefElement.SpatialDimension, 1);
                            int noOfNodesPerEdge = Kref.FaceRefElement.NoOfVertices + rule.NoOfNodes;

                            var FaceNodes = MultidimensionalArray.Create(noOfNodesPerEdge, myD);
                            for (int i = 0; i < Kref.FaceRefElement.NoOfVertices; i++) {
                                for (int d = 0; d < myD; d++) {
                                    FaceNodes[i, d] = Kref.FaceRefElement.Vertices[i, d];
                                }
                            }

                            int offset = Kref.FaceRefElement.NoOfVertices;
                            for (int i = 0; i < Kref.NoOfFaces; i++) {
                                for (int j = 0; j < rule.NoOfNodes; j++) {
                                    for (int d = 0; d < myD; d++) {
                                        FaceNodes[offset + j, d] = rule.Nodes[j, d];
                                    }
                                }
                            }

                            TestNodesPerFace[iKref] = new int[Kref.NoOfFaces];
                            TestNodesPerFace[iKref].SetAll(noOfNodesPerEdge);

                            MultidimensionalArray TstVtx = MultidimensionalArray.Create(noOfNodesPerEdge * Kref.NoOfFaces, D);
                            for (int iFace = 0; iFace < Kref.NoOfFaces; iFace++) {
                                Kref.GetFaceTrafo(iFace).Transform(
                                    FaceNodes, TstVtx.ExtractSubArrayShallow(
                                        new int[] { iFace * noOfNodesPerEdge, 0 },
                                        new int[] { (iFace + 1) * noOfNodesPerEdge - 1, D - 1 }));
                            }


                            // Guarantees that both cells are recognized if the level
                            // set directly passes through an edge.
                            // WARNING: This factor is completely arbitrary. This check
                            // should be replaced by something more sophisticated soon
                            // since it does not cover certain corner cases.
                            TstVtx.Scale(1.005);
                            TestNodes[iKref] = new NodeSet(Kref, TstVtx);
                        }
                    }



                    // clear cached level set values (level set may has changed)
                    foreach(var lsdh in DataHistories) {
                        lsdh.Current.ClearCaches();
                    }
                    
                }
                #endregion

                // evaluate level sets / find cut cells
                // ====================================
                #region UpdateTracker_FIND_CUT_CELLS
                using (new BlockTrace("FIND_CUT_CELLS", tr)) {
                    
                    int[][] VerticeInd = m_gDat.Cells.CellVertices;
                    {
                        // 1st sweep: find cut cells
                        // =========================

                        CellMask SearchMask;
                        if (incremental) {
                            

                            ushort[] __PrevLevSetRegions = RegionsHistory[0].m_LevSetRegions;

                            for (int levSetind = 0; levSetind < NoOfLevSets; levSetind++) {
                                BitArray _LevSetNeg = LevSetNeg[levSetind];
                                for (int j = 0; j < J; j++) {
                                    int dist = DecodeLevelSetDist(__PrevLevSetRegions[j], levSetind);
                                    _LevSetNeg[j] = dist < 0;
                                }
                            }

                            SearchMask = oldNearMask;
                        } else {
                            SearchMask = CellMask.GetFullMask(m_gDat);
                        }

                        int MaxVecLen = (int)Math.Ceiling(16000.0 / ((double)TestNodes.Max(ns => ns.NoOfNodes)));
                        MaxVecLen = Math.Max(1, MaxVecLen);
                        

                        foreach(var t_j0_Len in SearchMask.GetGeometricCellChunks(MaxVecLen, CellInfo.RefElementIndex_Mask | CellInfo.CellType_Mask)) { // loop over all cells in the search mask...
                            int j = t_j0_Len.Item1;
                            int VecLen = t_j0_Len.Item2;


                            int iKref = m_gDat.Cells.GetRefElementIndex(j);
                            var Kref = Krefs[iKref];
                            int noOfEdges = Kref.NoOfFaces;
                            int[] _TestNodesPerFace = this.TestNodesPerFace[iKref];

                            // loop over level sets ...
                            for (int levSetind = NoOfLevSets - 1; levSetind >= 0; levSetind--) {
                                var data = this.m_DataHistories[levSetind].Current;
                                MultidimensionalArray levSetVal = data.GetLevSetValues(this.TestNodes[iKref], j, VecLen);

                                for (int jj = 0; jj < VecLen; jj++) {
                                    bool Pos = false;
                                    bool Neg = false;

                                    int nodeIndex = 0;
                                    for (int e = 0; e < noOfEdges; e++) {
                                        bool PosEdge = false;
                                        bool NegEdge = false;

                                        for (int k = 0; k < _TestNodesPerFace[e]; k++) {
                                            double v = levSetVal[jj, nodeIndex];

                                            if (v < 0) {
                                                NegEdge = true;
                                            } else if (v > 0) {
                                                PosEdge = true;
                                            }

                                            nodeIndex++;
                                        }

                                        Pos |= PosEdge;
                                        Neg |= NegEdge;

                                        // Save sign of the edge
                                        {
                                            int edge = -1;
                                            for (int ee = 0; ee < GridDat.Cells.Cells2Edges[j + jj].Length; ee++) {
                                                int signedEdgeIndex = GridDat.Cells.Cells2Edges[j + jj][ee];
                                                int edgeIndex = Math.Abs(signedEdgeIndex) - 1;
                                                int inOut = signedEdgeIndex > 0 ? 0 : 1;
                                                if (GridDat.Edges.FaceIndices[edgeIndex, inOut] == e) {
                                                    edge = edgeIndex;
                                                }
                                            }

                                            if (edge < 0) {
                                                throw new Exception(
                                                    "Could not determine edge index; This should not have happened");
                                            }
                                        }
                                    }

                                    if ((Pos && Neg) || (!Pos && !Neg)) {
                                        // cell j+jj is cut by level set

                                        // code cell:
                                        EncodeLevelSetDist(ref LevSetRegionsUnsigned[j + jj], 0, levSetind);

                                    }

                                    if (Neg == true && Pos == false) {
                                        LevSetNeg[levSetind][j + jj] = true;
                                    } else {
                                        LevSetNeg[levSetind][j + jj] = false;
                                    }
                                }
                            }
                            j += VecLen;
                        }
                    }
                    #endregion

                    using (new BlockTrace("CELL_SWEEPS", tr)) {

                        // MPI update (cells)
                        // ==================

                        MPIUpdate(LevSetRegionsUnsigned, m_gDat);

                        // 2nd sweep: mark vertices
                        // ========================

                        // loop over level sets ...
                        for (int levSetind = NoOfLevSets - 1; levSetind >= 0; levSetind--) {
                            for (int _j = 0; _j < JA; _j++) {
                                int[] _VerticeInd = VerticeInd[_j];
                                int _NoOfSmplxVertice = _VerticeInd.Length;
                                Debug.Assert(_NoOfSmplxVertice == m_gDat.Cells.GetRefElement(_j).NoOfVertices);

                                if (DecodeLevelSetDist(LevSetRegionsUnsigned[_j], levSetind) == 0) {

                                    // code vertices:
                                    for (int k = 0; k < _NoOfSmplxVertice; k++) {
                                        int iVtx = _VerticeInd[k];
                                        EncodeLevelSetDist(ref VertexMarker[iVtx], 0, levSetind);
                                    }
                                }
                            }
                        }

                        // MPI Update (vertex markers)
                        // ===========================
                        MPIUpdateVertex(VertexMarker, VertexMarkerExternal, m_gDat, NoOfLevSets);
                    }


                    // find near and far regions
                    // =========================
                    {
                        ushort[] vtxMarkers = new ushort[NoOfSmplxVertice.Max()];
                        int[] vtsInd = new int[NoOfSmplxVertice.Max()];

                        for (int dist = 1; dist <= m_NearRegionWidth; dist++) {

                            // For every cell
                            for (int j = 0; j < J; j++) {
                                int[] _VerticeInd = VerticeInd[j];
                                int _NoOfSmplxVertice = _VerticeInd.Length;
                                Debug.Assert(_NoOfSmplxVertice == m_gDat.Cells.GetRefElement(j).NoOfVertices);

                                //for (int k = 0; k < NoOfSmplxVertice; k++)
                                //    vtxMarkers[k] = VertexMarker[VerticeInd[j, k]];
                                ushort[] __PrevLevSetRegions = null;
                                if (incremental)
                                    __PrevLevSetRegions = RegionsHistory[0].m_LevSetRegions;

                                for (int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) {

                                    int dJ = DecodeLevelSetDist(LevSetRegionsUnsigned[j], levSetInd);
                                    if (dJ >= dist) {

                                        int mindv = int.MaxValue;
                                        for (int k = 0; k < _NoOfSmplxVertice; k++)
                                            mindv = Math.Min(mindv, DecodeLevelSetDist(VertexMarker[_VerticeInd[k]], levSetInd));
                                        if (incremental) {
                                            // The following line ensures that when the level-set is leaving the computational domain,
                                            // we still have a near-cell on the boundary.
                                            int prevDist = Math.Abs(DecodeLevelSetDist(__PrevLevSetRegions[j], levSetInd));
                                            mindv = Math.Min(mindv, prevDist);
                                        }
                                        mindv++;

                                        if (mindv == dist) {
                                            EncodeLevelSetDist(ref LevSetRegionsUnsigned[j], mindv, levSetInd);

                                            for (int k = 0; k < _NoOfSmplxVertice; k++) {
                                                int dVk = DecodeLevelSetDist(VertexMarker[_VerticeInd[k]], levSetInd);
                                                dVk = Math.Min(dVk, mindv);
                                                EncodeLevelSetDist(ref VertexMarker[_VerticeInd[k]], dVk, levSetInd);
                                            }
                                        }
                                    }
                                }
                            }

                            // MPI update
                            // ----------
                            MPIUpdate(LevSetRegionsUnsigned, m_gDat);
                            MPIUpdateVertex(VertexMarker, VertexMarkerExternal, m_gDat, NoOfLevSets);
                        }
                    }

                    // set the sign
                    // ============
                    {
                        Array.Copy(LevSetRegionsUnsigned, LevSetRegions, J);

                        for (int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) {
                            for (int j = 0; j < J; j++) {

                                if (LevSetNeg[levSetInd][j]) {
                                    int dist = DecodeLevelSetDist(LevSetRegions[j], levSetInd);
                                    dist *= -1;
                                    EncodeLevelSetDist(ref LevSetRegions[j], dist, levSetInd);
                                }
                            }
                        }

                        MPIUpdate(LevSetRegions, m_gDat);
                    }

                    // recalculate m_LenToNextChange
                    // =============================
                    RegionsHistory.Current.Recalc_LenToNextchange();
                }

               
                // forget values that are not correct anymore
                // ==========================================
                this.m_QuadFactoryHelpersHistory.Current.Clear();
                this.m_XDGSpaceMetricsHistory.Current.Clear();
                
                // check the LevelSet CFL
                // ======================

                bool throwCFL;
                int[] fail;
                using (new BlockTrace("CFL_CHECK", tr)) {

                    fail = new int[NoOfLevSets];
                    if(this.PopulatedHistoryLength > 0) {
                        m_LevSetAllowedMovement = __LevSetAllowedMovement;

                        // cannot be moved down because we need the OLD subgrid
                        throwCFL = false;
                        for(int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) {

                            if(m_LevSetAllowedMovement[levSetInd] <= m_NearRegionWidth)
                                fail[levSetInd] = CheckLevelSetCFL(levSetInd);
                            if(fail[levSetInd] > 0)
                                throwCFL = true;

                        }
                    } else {
                        throwCFL = false;
                    }
                }

                // update memory of XDG fields, etc.
                // =================================

                ObserverUpdate();

                // throw exception, if levelset CFL violated
                // =========================================
                if (throwCFL) {
                    LevelSetCFLException exception = new LevelSetCFLException(fail);
                    foreach (var reference in m_Observers) {
                        IObserver<LevelSetRegions> observer = reference.Target;
                        if (observer != null) {
                            observer.OnError(exception);
                        }
                    }

                    throw exception;
                }

            }
        }


        /// <summary>
        /// Clears all internal references for this object, to make sure that any attempt to use it leads to an exception.
        /// </summary>
        public void Invalidate() {
            this.m_RegionsHistory = null;
            this.m_SpeciesNames = null;
            this.m_SpeciesIndex2Id = null;
            this.m_SpeciesIndex = null;
            this.m_gDat = null;
            this.m_GhostTable = null;
            this.m_DataHistories = null;
            this.m_LevelSetHistories = null;
            this.m_LevelSetSignCodes = null;
            this.m_LevSetAllowedMovement = null;
            this.m_NoOfSpecies = null;
            this.m_Observers = null;
            this.m_XDGSpaceMetricsHistory = null;
            this.m_QuadFactoryHelpersHistory = null;
            this.m_SpeciesId2Index = null;
            this.TestNodes = null;
            this.TestNodesPerFace = null;
        }


        private void ObserverUpdate() {
            using (new FuncTrace()) {

                // Remove obsolete observers from list...
                // ======================================
                for (int i = 0; i < m_Observers.Count; i++) {
                    if (!m_Observers[i].IsAlive) {
                        m_Observers.RemoveAt(i);
                        i--;
                    }
                }

                // update memory of all registered fields
                // =====================================
                // a disadvantage of this notification-by-weak-ref -- approach
                // is that the 'UpdateMemory' may be called also for
                // objects that are already unused but not yet collected...
                // A solution would be to call GC.Collect(), but it is not known
                // whether a GC run or update of unused memory is more expensive.

                
                // call the update method of all active fields
                foreach (var reference in m_Observers) {
                    IObserver<LevelSetRegions> observer = reference.Target;
                    if (observer != null) {
                        reference.Target.OnNext(Regions);
                    }
                }
            }
        }

        static void MPIUpdateVertex(ushort[] VertexMarker, ushort[,] VertexMarkerExternal, GridData grdDat, int NoOfLevelSets) {
            int Size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out Size);
            if (Size <= 1)
                return;

            int[] sndProc = grdDat.Parallel.ProcessesToSendTo;
            int[] rvcProc = grdDat.Parallel.ProcessesToReceiveFrom;

            MPI_Request[] rqst = new MPI_Request[sndProc.Length + rvcProc.Length];
            MPI_Status[] staTussies = new MPI_Status[rqst.Length];

            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);

            int N = VertexMarkerExternal.GetLength(1);

            var VerticeInd = grdDat.Cells.CellVertices;

            int J = grdDat.Cells.NoOfLocalUpdatedCells;
            int JA = grdDat.Cells.NoOfCells;


            unsafe
            {
                fixed (ushort* pCA = &VertexMarkerExternal[0, 0])
                {
                    ushort[][,] SendBuffers = new ushort[sndProc.Length][,];
                    GCHandle[] SendBufferPin = new GCHandle[sndProc.Length];
                    int cnt = 0;

                    try {

                        // Sending ...
                        // -----------

                        // over all processes to which we have to send data to ...
                        for (int i = 0; i < sndProc.Length; i++) {

                            // destination processor and comm list
                            int pDest = sndProc[i];
                            var commList = grdDat.Parallel.SendCommLists[pDest];
                            int Len = commList.Length;

                            // build&fill send buffer
                            SendBuffers[i] = new ushort[Len, N];
                            var SendBuffer = SendBuffers[i];
                            for (int l = 0; l < Len; l++) {
                                int[] cellVtx = VerticeInd[commList[l]];
                                int _N = cellVtx.Length;
                                for (int n = 0; n < _N; n++) {
                                    int iVtx = cellVtx[n];
                                    SendBuffer[l, n] = VertexMarker[iVtx];
                                }
                            }

                            // MPI send
                            SendBufferPin[i] = GCHandle.Alloc(SendBuffers[i], GCHandleType.Pinned);
                            cnt++;
                            csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(SendBuffers[i], 0),
                                Len * sizeof(ushort) * N, csMPI.Raw._DATATYPE.BYTE, pDest,
                                44442 + MyRank,
                                csMPI.Raw._COMM.WORLD,
                                out rqst[i]);
                        }


                        // Receiving ...
                        // -------------

                        // over all processes from which we receive data...
                        for (int i = 0; i < rvcProc.Length; i++) {

                            // Source processor and insert index and no of elements to receive ...
                            int pOrigin = rvcProc[i];
                            int iInsert = grdDat.Parallel.RcvCommListsInsertIndex[pOrigin] - J;
                            int Len = grdDat.Parallel.RcvCommListsNoOfItems[pOrigin];

                            // MPI receive
                            csMPI.Raw.Irecv((IntPtr)(pCA + iInsert * N),
                                Len * sizeof(ushort) * N, csMPI.Raw._DATATYPE.BYTE, pOrigin,
                                44442 + pOrigin,
                                csMPI.Raw._COMM.WORLD,
                                out rqst[i + sndProc.Length]);
                        }

                        // Wait for comm to finish
                        // -----------------------

                        csMPI.Raw.Waitall(rqst.Length, rqst, staTussies);

                    } finally {

                        // release GC handles
                        // ==================

                        for (int i = 0; i < cnt; i++)
                            SendBufferPin[i].Free();
                    }
                }
            }

            for (int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                for (int j = J; j < JA; j++) {
                    int _j = j - J;

                    int[] _VerticeInd = VerticeInd[j];
                    int _N = _VerticeInd.Length;

                    for (int n = 0; n < _N; n++) {
                        var ext = DecodeLevelSetDist(VertexMarkerExternal[_j, n], iLs);
                        var in_ = DecodeLevelSetDist(VertexMarker[_VerticeInd[n]], iLs);

                        int iNeu = Math.Min(ext, in_);
                        EncodeLevelSetDist(ref VertexMarker[_VerticeInd[n]], iNeu, iLs);
                    }

                }
            }

        }


        static void MPIUpdate(ushort[] codesArray, GridData grdDat) {
            Comm.VectorTransceiver_Ext.MPIExchange<ushort[], ushort>(codesArray, grdDat);
        }

        int[] m_LevSetAllowedMovement;


        #region IObservable<LevelSetInfo> Members

        /// <summary>
        /// See <see cref="Subscribe"/>;
        /// </summary>
        private List<BoSSS.Platform.WeakReference<IObserver<LevelSetRegions>>> m_Observers =
            new List<BoSSS.Platform.WeakReference<IObserver<LevelSetRegions>>>();

        /// <summary>
        /// At construction time, objects may register themselves at the level
        /// set tracker in order to be updated when the level set information
        /// changes. 
        /// </summary>
        /// <param name="observer">
        /// The observer that wants to register.
        /// </param>
        /// <remarks>
        /// The reference to the subscribed will be "weak" (see
        /// <see cref="WeakReference{T}"/>) so that garbage collection will not
        /// be affected by the subscription.
        /// </remarks>
        public IDisposable Subscribe(IObserver<LevelSetRegions> observer) {
            BoSSS.Platform.WeakReference<IObserver<LevelSetRegions>> reference =
                new BoSSS.Platform.WeakReference<IObserver<LevelSetRegions>>(observer);
            m_Observers.Add(reference);
            return new Unsubscriber(this, reference);
        }

        /// <summary>
        /// Unsubscriber for objects that have subscribed to objects this
        /// object through <see cref="Subscribe"/>.
        /// </summary>
        private class Unsubscriber : IDisposable {

            /// <summary>
            /// The creator of this object.
            /// </summary>
            private LevelSetTracker owner;

            /// <summary>
            /// The observer that may want to unregister itself.
            /// </summary>
            private BoSSS.Platform.WeakReference<IObserver<LevelSetRegions>> observer;

            /// <summary>
            /// Creates a new unsubscriber.
            /// </summary>
            /// <param name="owner">
            /// The creator of this object.
            /// </param>
            /// <param name="observer">
            /// The observer that may want to unregister itself.
            /// </param>
            public Unsubscriber(LevelSetTracker owner, BoSSS.Platform.WeakReference<IObserver<LevelSetRegions>> observer) {
                this.owner = owner;
                this.observer = observer;
            }

            #region IDisposable Members

            /// <summary>
            /// Removes the given observer from the list of registered
            /// observers of the observed object.
            /// </summary>
            public void Dispose() {
                if (observer.IsAlive && owner.m_Observers.Contains(observer)) {
                    owner.m_Observers.Remove(observer);
                }
            }

            #endregion
        }

        #endregion

        #region IDisposable Members

        /// <summary>
        /// Disposes the references to all registered observers (see
        /// <see cref="Subscribe"/>)
        /// </summary>
        public void Dispose() {
            foreach (var reference in m_Observers) {
                IObserver<LevelSetRegions> observer = reference.Target;
                if (observer != null) {
                    observer.OnCompleted();
                }

                reference.Dispose();
            }

            m_Observers = null;
        }

        #endregion
    }



}
