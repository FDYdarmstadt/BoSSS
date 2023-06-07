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
using System.IO;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// The level set - tracker manages the tracking of narrow - bands of up to 4 level set.
    /// It is a perquisite for
    /// the memory management of cut-cell fields, i.e. <see cref="XDGField"/> and <see cref="XDGBasis"/>
    /// </summary>
    /// <remarks>
    /// After changing one ore more level sets, the <see cref="UpdateTracker"/>-method must be
    /// called.
    /// </remarks>
    public partial class LevelSetTracker : IObservable<LevelSetTracker.LevelSetRegions>, IDisposable {

        /// <summary>
        /// Region code (see <see cref="LevelSetRegions.m_LevSetRegions"/>) indicating that all level set values in a cell are in
        /// the positive far region
        /// </summary>
        public const ushort AllFARplus = 0xffff;

        /// <summary>
        /// Region code (see <see cref="LevelSetRegions.m_LevSetRegions"/>) indicating that all level set values in a cell are in
        /// the negative far region
        /// </summary>
        public const ushort AllFARminus = 0x1111;

        /// <summary>
        /// Region code (see <see cref="LevelSetRegions.m_LevSetRegions"/>) indicating that all level sets cut the actual cell 
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
        /// 
        /// </summary>
        /// <param name="BackgroundGrid"></param>
        /// <param name="__NearRegionWidth"></param>
        /// <param name="SpeciesTable"></param>
        /// <param name="levSets"></param>
        /// <param name="cutCellquadType">
        /// the type of integration in cut-cells; if more than one type is required within a single application, two <see cref="LevelSetTracker"/>'s should be used.
        /// </param>       
        public LevelSetTracker(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, Array SpeciesTable, params ILevelSet[] levSets) {
            ConstructorCommon(BackgroundGrid, __NearRegionWidth, SpeciesTable, cutCellquadType, levSets);
        }

        
        private void ConstructorCommon(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, Array SpeciesTable, params ILevelSet[] levSets) {
            ConstructorCommon(BackgroundGrid, __NearRegionWidth, SpeciesTable, cutCellquadType, levSets);
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
        /// indexes of all available times, i.e. all valid indexes to access any <see cref="HistoryStack{T}"/>
        /// </summary>
        public int[] PopulatedHistoryIndices {
            get {
                int L = HistoryLength;
                return RegionsHistory.PopulatedIndices;
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
        private void ConstructorCommon(GridData griData, int __NearRegionWidth, Array SpeciesTable, XQuadFactoryHelper.MomentFittingVariants cutCellQuadratureType, params ILevelSet[] levSets) {
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


            UpdateTracker(0.0);
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
        bool[,] m_GhostTable;

        /// <summary>
        /// initializes <see cref="m_GhostTable"/>
        /// </summary>
        void ComputeGhostTable() {
            m_GhostTable = new bool[this.m_SpeciesNames.Count,81];
            int NoOfLevSets = this.NoOfLevelSets;

            foreach (String specNmn in this.m_SpeciesNames) {
                SpeciesId specId = GetSpeciesId(specNmn);
                LevelSetSignCode[] speciesCodes = GetLevelSetSignCodes(specNmn);

                bool[] Table = new bool[81];

                LevelsetCellSignCode cd;
                for (cd.lsSig = 0; cd.lsSig < 81; cd.lsSig++) {

                    foreach (var scd in speciesCodes) {
                        if (cd.IsContained(scd, NoOfLevSets))
                            Table[cd.lsSig] = true;
                    }
                }
                m_GhostTable.SetRow(specId.cntnt - ___SpeciesIDOffest, Table);
            }
        }

        /// <summary>
        /// Detects whether a species is present in some cell 
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
            return m_GhostTable[id.cntnt - ___SpeciesIDOffest, cd.lsSig];
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
        /// Gets a copy of the species table, a multidimensional string-array.
        /// 
        /// The species table defines, which species corresponds with 
        /// which combination of level set - signs.
        /// </summary>
        /// <remarks>
        /// The rank (number of dimensions) of this array is equal to the number of level-sets (<see cref="LevelSets"/>).
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
        /// Returns all <see cref="LevelSetRegions.Time"/> entries in the <see cref="RegionsHistory"/> stack.
        /// </summary>
        public double[] TimeLevelsInStack {
            get {
                return RegionsHistory.AvailableIndices.Select((int iHist) => RegionsHistory[iHist].Time).ToArray();
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
        /// Returns the species that are speratated by Level Set <paramref name="levSetIdx"/>, as defined in <see cref="m_SpeciesTable"/>. 
        /// </summary>
        /// <param name="levSetIdx"></param>
        /// <returns></returns>
        public IEnumerable<string> GetSpeciesSeparatedByLevSet(int levSetIdx) {
            LinkedList<SpeciesPair> speciesPairsOfLevelSet = FindSeparatedSpeciesPairs(levSetIdx);
            LinkedList<string> speciesOfLevelSet = new LinkedList<string>();
            foreach (SpeciesPair pair in speciesPairsOfLevelSet) {
                if (!speciesOfLevelSet.Contains(pair.NegativeSpecies)) {
                    speciesOfLevelSet.AddLast(pair.NegativeSpecies);
                }
                if (!speciesOfLevelSet.Contains(pair.PositiveSpecies)) {
                    speciesOfLevelSet.AddLast(pair.PositiveSpecies);
                }
            }
            return speciesOfLevelSet;
        }

        /// <summary>
        /// Returns the species that are speratated by Level Set <paramref name="levSetIdx"/>, as defined in <see cref="m_SpeciesTable"/>. 
        /// </summary>
        /// <param name="levSetIdx"></param>
        /// <returns>
        /// Species Pair. First species is negative species, second species is positive species.
        /// </returns>
        public (string, string)[] GetSpeciesPairsSeparatedByLevSet(int levSetIdx) {
            LinkedList<SpeciesPair> speciesOfLevelSet = FindSeparatedSpeciesPairs(levSetIdx);
            (string, string)[] pairs = new (string, string)[speciesOfLevelSet.Count];
            int i = 0;
            foreach (SpeciesPair pair in speciesOfLevelSet) {
                pairs[i].Item1 = pair.NegativeSpecies;
                pairs[i].Item2 = pair.PositiveSpecies;
                ++i;
            }
            return pairs;
        }


        /// <summary>
        /// Returns species that are separated by levelSet with No. <paramref name="levSetIdx"/>
        /// </summary>
        LinkedList<SpeciesPair> FindSeparatedSpeciesPairs(int levSetIdx) {
            int[] levelSetSigns = new int[SpeciesTable.Rank];
            LinkedList<SpeciesPair> speciesOfLevelSet = new LinkedList<SpeciesPair>();
            FindSeparatedSpecies(levelSetSigns, 0, levSetIdx, speciesOfLevelSet);
            return speciesOfLevelSet;
        }

        struct SpeciesPair : IEquatable<SpeciesPair> {
            public string PositiveSpecies;
            public string NegativeSpecies;

            public SpeciesPair(string negativeSpecies, string positiveSpecies) {
                this.PositiveSpecies = positiveSpecies;
                this.NegativeSpecies = negativeSpecies;
            }

            public bool Equals(SpeciesPair other) {
                return other.PositiveSpecies == PositiveSpecies && other.NegativeSpecies == NegativeSpecies;
            }
        }

        void FindSeparatedSpecies(int[] levelSetSigns, int level, int levSetIdx, LinkedList<SpeciesPair> species) {
            if (level == levelSetSigns.Length) {
                levelSetSigns[levSetIdx] = 0;
                string negativeSpecies = (string)SpeciesTable.GetValue(levelSetSigns);
                levelSetSigns[levSetIdx] = 1;
                string positiveSpecies = (string)SpeciesTable.GetValue(levelSetSigns);
                if (negativeSpecies != positiveSpecies) {
                    SpeciesPair pair = new SpeciesPair(negativeSpecies, positiveSpecies);
                    if (!species.Contains(pair))
                        species.AddLast(pair);
                }
            } else if (level == levSetIdx) {
                FindSeparatedSpecies(levelSetSigns, level + 1, levSetIdx, species);
            } else {
                levelSetSigns[level] = 0;
                FindSeparatedSpecies(levelSetSigns, level + 1, levSetIdx, species);
                levelSetSigns[level] = 1;
                FindSeparatedSpecies(levelSetSigns, level + 1, levSetIdx, species);
            }
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
        /// Width, in Number of cells, of the near field (set as an argument of <see cref="UpdateTracker"/>);
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
        public IReadOnlyList<HistoryStack<ILevelSet>> LevelSetHistories {
            get {
                return m_LevelSetHistories;
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
        /// throws away the top of the history, reverts to the latest pushed state (<see cref="PushStacks"/>,
        /// see also <see cref="DataHistories"/>, <see cref="LevelSetHistories"/>, <see cref="RegionsHistory"/>, etc.).
        /// 
        /// The top state of the stack can also be back-upped/resored by
        /// - <see cref="BackupTimeLevel(int)"/>
        /// - <see cref="ReplaceCurrentTimeLevel(TrackerBackup)"/>
        /// </summary>
        public void PopStacks() {
            int NoOfLs = LevelSets.Count;
            int PHL = PopulatedHistoryLength;
            int HL = this.HistoryLength;


            Debug.Assert(NoOfLs == m_LevelSets.Count);
            for(int iLs = 0; iLs < NoOfLs; iLs++) {
                LevelSetHistories[iLs].Pop((ls1, ls0) => {
                    ls1.CopyFrom(ls0);
                    return ls1;
                });//.Push((ls1) => ls1, (ls1, ls0) => ls1.CloneAs());
            }

            Debug.Assert(NoOfLs == m_DataHistories.Count);
            for(int iLs = 0; iLs < NoOfLs; iLs++) {
                DataHistories[iLs].Pop((data1, data0) => data0);//.Push((data1) => new LevelSetData(this, iLs), (data1, data0) => data1);

                // fix the history index...
                for(int iStack = 1; iStack > -m_DataHistories[iLs].GetPopulatedLength(); iStack--) {
                    m_DataHistories[iLs][iStack].m_HistoryIndex = iStack;
                }
            }

            m_RegionsHistory.Pop((r1, r0) => r0); // .Push((r1) => r1.CloneAs(), (r1, r0) => r1);

            m_QuadFactoryHelpersHistory.Pop((r1, r0) => r0);

            m_XDGSpaceMetricsHistory.Pop((r1, r0) => r0);

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
        /// The level sets, identical to the top of the level-set stack (<see cref="LevelSetHistories"/>)
        /// - 1st index: index of level-set
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
        /// <see cref="GetSpeciesIndex(ReducedRegionCode, SpeciesId)"/>;
        /// The reduced region code in 3-adic representation;
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// Here, three states (the reduced region) for each of the for level sets are considered:
        /// - positive far (FAR+)
        /// - negative far (FAR-)
        /// - positive near, negative near or cut
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
        /// - the cell if in the positive FAR region: <i>distance</i> == 7, i.e. the <i>code</i> is 0xf
        /// - the cell if in the negative FAR region: <i>distance</i> == -1, i.e the <i>code</i> is 0x1
        /// - the cell is cut or in the near region: -6 &lt; <i>distance</i> &lt; 6
        /// So, there are in total 3*4 states for some cell (for 4 level sets);
        /// </remarks>
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
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpeciesByRegionCode(ushort, out ReducedRegionCode)"/>;
        /// </param>
        /// <param name="LevelSetSignBytecode"></param>
        /// <returns></returns>
        public int GetSpeciesIndex(ReducedRegionCode rrc, LevelSetSignCode LevelSetSignBytecode) {
            int ReducedRegionCode = rrc.rrc << 4;    // equal to int *= 16
            ReducedRegionCode = ReducedRegionCode | LevelSetSignBytecode.val; // equal to ind += LevelSetSignBytecode, if  0 <= LevelSetSignBytecode < 16
            return m_SpeciesIndex[ReducedRegionCode];
        }

        /// <summary>
        /// this function is the inverse to <see cref="GetSpeciesIndex(ReducedRegionCode, LevelSetSignCode)"/>;
        /// </summary>
        /// <param name="_ReducedRegionCode">
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpeciesByRegionCode(ushort, out ReducedRegionCode)"/>;
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
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpeciesByRegionCode(ushort, out ReducedRegionCode)"/>;
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

        ///// <summary>
        ///// 
        ///// </summary>
        ///// <param name="levSetInd"></param>
        ///// <returns></returns>
        //public ushort LevelSetBitmask(int levSetInd) {
        //    if (levSetInd < 0 || levSetInd >= NoOfLevelSets)
        //        throw new ArgumentException("unknown level set.", "levset");

        //    return (ushort)(0xF << (4 * levSetInd));
        //}


        /// <summary>
        /// Extracts the distance layer index for level-set <paramref name="levSetInd"/>
        /// from the code <paramref name="code"/> (see <see cref="LevelSetRegions.m_LevSetRegions"/>).
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
        /// int the code <paramref name="code"/> (<see cref="LevelSetRegions.m_LevSetRegions"/>).
        /// </summary>
        static public void EncodeLevelSetDist(ref ushort code, int dist, int levSetInd) {
            Debug.Assert(-7 <= dist);
            Debug.Assert(dist <= +7);
            int c = ((dist + 8) << (4 * levSetInd));
            int msk = (0xf << (4 * levSetInd));
            code = (ushort)((code & ~msk) | c);
        }
        
             



        /// <summary>
        /// Node-set for determination of cut-cells.
        /// - index: reference element index
        /// </summary>
        /// <remarks>
        /// Structured as follows:
        /// - a set of quadrature nodes on each face
        /// - nodes slightly exterior to the cell
        /// </remarks>
        NodeSet[] TestNodes;

        /// <summary>
        /// Quadrature weights 
        /// - index: reference element index, correlates with 1st index of <see cref="TestNodes"/>
        /// </summary>
        MultidimensionalArray[] TestNodes_QuadWeights;

        /// <summary>
        /// Fore each node set in <see cref="TestNodes"/>, their correlation to the element faces
        /// - 1st index: reference element index, correlates with 1st index of <see cref="TestNodes"/>
        /// - 2nd index: face index
        /// </summary>
        int[][] TestNodesPerFace;

        BernsteinTransformator[] m_TestTransformer;

        /// <summary>
        /// Transformator for each level set to transform the legendre coefficients to bernstein coefficients
        /// we then use there superior geometric properties to search for cut cells.
        /// here the bernstein coefficients lie on a slight offset of the reference element
        /// </summary>
        BernsteinTransformator[] TestTransformer { 
            get{
                if (m_TestTransformer == null) {
                    m_TestTransformer = new BernsteinTransformator[this.LevelSets.Count];
                    for(int i = 0; i< this.LevelSets.Count; i++ ) {
                        if (this.LevelSets[i] is LevelSet ls) {
                            m_TestTransformer[i] = new BernsteinTransformator(ls.Basis, 0.005);
                        }
                    }
                }
                return m_TestTransformer;
            } 
        }

        BernsteinTransformator[] m_TestTransformerEdges;

        /// <summary>
        /// Transformator for each level set to transform the legendre coefficients to bernstein coefficients
        /// we then use there superior geometric properties to search for cut cells.
        /// here the bernstein coefficients lie exactly on the edges of the reference element
        /// </summary>
        BernsteinTransformator[] TestTransformerEdges {
            get {
                if (m_TestTransformerEdges == null) {
                    m_TestTransformerEdges = new BernsteinTransformator[this.LevelSets.Count];
                    for (int i = 0; i < this.LevelSets.Count; i++) {
                        if (this.LevelSets[i] is LevelSet ls)
                            m_TestTransformerEdges[i] = new BernsteinTransformator(ls.Basis);
                    }
                }
                return m_TestTransformerEdges;
            }
        }

        private void DefineTestNodes() {
            int D = this.GridDat.SpatialDimension;
            var Krefs = m_gDat.Grid.RefElements;
            int MaxLsDegree = this.LevelSets.Select(LevSet => (LevSet is DGField dgLS ? dgLS.Basis.Degree : 2)).Max();

            if(TestNodes == null) { // only the first time in object lifetime ...
                TestNodes = new NodeSet[Krefs.Length];
                TestNodesPerFace = new int[Krefs.Length][];
                TestNodes_QuadWeights = new MultidimensionalArray[Krefs.Length];

                for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                    var Kref = Krefs[iKref];
                    int myD = Math.Max(GridDat.Grid.GetRefElement(0).FaceRefElement.SpatialDimension, 1);
                    int NoOfFaces = Kref.NoOfFaces;
                    TestNodesPerFace[iKref] = new int[NoOfFaces];
                    
                    // transformation from face to reference element
                    int TransformFaceNodes(double scaling, int iFace, MultidimensionalArray FaceNodesIn, MultidimensionalArray TstVtxOut, int OffsetIntoTestVtx) {
                        int NN = FaceNodesIn.GetLength(0);
                        var Res = TstVtxOut.ExtractSubArrayShallow(
                                        new int[] { OffsetIntoTestVtx, 0 },
                                        new int[] { OffsetIntoTestVtx + NN - 1, D - 1 });

                        Kref.GetFaceTrafo(iFace).Transform(FaceNodesIn, Res);

                        Res.Scale(scaling);
                        return NN;
                    }

     

                    // various sets on faces
                    QuadRule BruteRule = Kref.FaceRefElement.GetBruteForceQuadRule(4, 2); // brute-force rule on face
                    //QuadRule BruteRule = Kref.FaceRefElement.GetBruteForceQuadRule(2, 1); // brute-force rule on face
                    NodeSet corners = Kref.FaceRefElement.Vertices;
                    QuadRule GaussRule = Kref.FaceRefElement.GetQuadratureRule(MaxLsDegree * 2); // Gauss rule on face
                    //QuadRule GaussRule = Kref.FaceRefElement.GetQuadratureRule(MaxLsDegree);
                    TestNodes_QuadWeights[iKref] = GaussRule.Weights;

                    // allocate memory form test node set 
                    int TotNumberOfNodes = NoOfFaces * (corners.NoOfNodes + GaussRule.NoOfNodes + BruteRule.NoOfNodes);
                    NodeSet TstVtx = new NodeSet(Kref, TotNumberOfNodes, D, true);
                    int offset = 0;


                    // step 1: place gauss rules on edges
                    // ----------------------------------
                    for(int iFace = 0; iFace < NoOfFaces; iFace++) {
                        int NN = TransformFaceNodes(1.0, iFace, GaussRule.Nodes, TstVtx, offset);
                        TestNodesPerFace[iKref][iFace] = NN;
                        offset += NN;
                    }
                                       
                    
                    // step 2: Place knots slightly outside the perimeter of the element
                    // ----------------------------------
                    for(int iFace = 0; iFace < NoOfFaces; iFace++) {
                        offset += TransformFaceNodes(1.005, iFace, corners, TstVtx, offset);
                        offset += TransformFaceNodes(1.005, iFace, BruteRule.Nodes, TstVtx, offset);
                    }

                    // record the newly created node set
                    // ---------------------------------
                    TstVtx.LockForever();
                    TestNodes[iKref] = TstVtx;
                }
            }
        }




        /// <summary>
        /// a counter that is increased every time when <see cref="UpdateTracker"/> is called;
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
                for (int j = 0; j < m_gDat.Cells.Count; j++) {
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
            using(new FuncTrace()) {
                ushort[] oldCode = this.RegionsHistory[0].m_LevSetRegions;
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                var msk = new BitArray(J);
           
                CellMask newCut = this.RegionsHistory[1].GetCutCellSubgrid4LevSet(LevSetIdx).VolumeMask;
                int fail_count = 0;

                // for all cells that are cut by the levelset,
                // check whether they are in Near - region of the previous state;
                foreach(int j in newCut.ItemEnum) {
                    int old_dist = LevelSetTracker.DecodeLevelSetDist(oldCode[j], LevSetIdx);
                if (Math.Abs(old_dist) > 1) {
                        fail_count++;
                        msk[j] = true;
                    }
                }

                int failCountGlobal = fail_count.MPISum();
            if (failCountGlobal > 0 && GridDat.MpiRank == 0)
                    (new CellMask(this.GridDat, msk)).SaveToTextFile("fail.csv", WriteHeader: false);

                return failCountGlobal;
            }
        }


        /// <summary>
        /// The minimum amount of data required to restore the level-set-tracker state
        /// after mesh adaptation
        /// </summary>
        public class EssentialTrackerBackup {
            
            /// <summary>
            /// clone of the level-set fields, <see cref="m_LevelSetHistories"/>
            /// </summary>
            public SinglePhaseField[] LevelSets;

            /// <summary>
            /// level-set version index, <see cref="LevelSetRegions.Version"/>
            /// </summary>
            public int Version;
            
            /// <summary>
            /// associated physical time, <see cref="LevelSetRegions.Time"/>
            /// </summary>
            public double time;

        }


        /// <summary>
        /// Full backup of the tracker state at a respective time-step.
        /// </summary>
        public class TrackerBackup : EssentialTrackerBackup, ICloneable {
            /// <summary>
            /// backup of region code for each cell, <see cref="LevelSetRegions.RegionsCode"/>
            /// </summary>
            public ushort[] Regions;


            /// <summary>
            /// backup of <see cref="LevelSetRegions.m_LevSetCoincidingFaces"/>
            /// </summary>
            public (int iLevSet, int iFace)[][] LevSetCoincidingFaces;

            /// <summary>
            /// non-shallow cloning
            /// </summary>
            public object Clone() {
                var clone_LevSetCoincidingFaces = new (int iLevSet, int iFace)[LevSetCoincidingFaces.Length][];
                for(int j = 0; j < clone_LevSetCoincidingFaces.Length; j++) {
                    var lscf_j = LevSetCoincidingFaces[j];
                    if(lscf_j != null) {
                        clone_LevSetCoincidingFaces[j] = lscf_j.CloneAs();
                        Debug.Assert(lscf_j.Length <= 0 || !object.ReferenceEquals(lscf_j[0], clone_LevSetCoincidingFaces[j][0]));
                    }
                }


                return new TrackerBackup() {
                    LevelSets = this.LevelSets.Select(ls => ls.CloneAs()).ToArray(),
                    Version = this.Version,
                    time = this.time,
                    Regions = this.Regions.CloneAs(),
                    LevSetCoincidingFaces = clone_LevSetCoincidingFaces
                };
            }
        }

        
        /// <summary>
        /// Backup of the internal state of the level-set tracker for a certain history stack index (<paramref name="iHistory"/>).
        /// Not the entire state is backed up (cell masks, cut-cell quadrature rules, mass matrices, ..., are lost),
        /// but everything essential to restore a certain state.
        /// </summary>
        /// <param name="iHistory">History stack index.</param>
        /// <returns>
        /// Can be used as input for <see cref="ReplaceCurrentTimeLevel(TrackerBackup)"/> or <see cref="ReplaceCurrentTimeLevel(EssentialTrackerBackup)"/>.
        /// </returns>
        public TrackerBackup BackupTimeLevel(int iHistory) {
            int Jup = this.GridDat.Cells.NoOfLocalUpdatedCells;
            
            ushort[] RegionClone = new ushort[Jup];
            Array.Copy(this.RegionsHistory[iHistory].RegionsCode, 0, RegionClone, 0, Jup);

            var LSCF = this.RegionsHistory[iHistory].m_LevSetCoincidingFaces;
            (int iLevSet, int iFace)[][] _LevSetCoincidingFaces;
            _LevSetCoincidingFaces = new (int iLevSet, int iFace)[Jup][];
            if(LSCF != null) {
                for(int j = 0; j < Jup; j++) {
                    if (_LevSetCoincidingFaces[j] != null)
                        _LevSetCoincidingFaces[j] = LSCF[j].CloneAs();
                    else
                        _LevSetCoincidingFaces[j] = null;
                }
            }

            int NoOfLevelSets = this.NoOfLevelSets;
            LevelSet[] LevSetClones = new LevelSet[NoOfLevelSets];
            for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                LevelSet Ls = (LevelSet)(this.LevelSetHistories[iLs][iHistory]);
                LevSetClones[iLs] = Ls;
            }

            return new TrackerBackup() {
                LevelSets = LevSetClones,
                LevSetCoincidingFaces = _LevSetCoincidingFaces,
                Regions = RegionClone,
                Version = this.RegionsHistory[iHistory].Version,
                time = this.RegionsHistory[iHistory].Time
            };
        }

        /// <summary>
        /// Counterpart of <see cref="BackupTimeLevel(int)"/>, 
        /// but the region codes (<see cref="LevelSetRegions.RegionsCode"/>) are
        /// re-computed by calling <see cref="UpdateTracker"/>.
        /// </summary>
        /// <remarks>
        /// Used for **mesh adaptation**,
        /// i.e. when the mesh changes and region codes must be re-computed.
        /// </remarks>
        public void ReplaceCurrentTimeLevel(EssentialTrackerBackup backup) {

            SinglePhaseField[] LevSet = backup.LevelSets;
            int VersionCounter = backup.Version;
            double time = backup.time;

            if(LevSet.Length != this.NoOfLevelSets)
                throw new ArgumentOutOfRangeException();
            int NoOfLevelSet = this.NoOfLevelSets;

            // invalidate everything we got so far
            // ===================================

            this.Regions.InvalidateCaches();
            for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                this.DataHistories[iLs].Current.ClearCaches();
            }

            m_QuadFactoryHelpersHistory.Current.Clear();

            m_XDGSpaceMetricsHistory.Current.Clear();

            // set level-set data
            // ==================
                        
            for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                LevelSet Ls = (LevelSet)(this.LevelSets[iLs]);
                Ls.Clear();
                Ls.Acc(1.0, LevSet[iLs]);
            }

            // update tracker
            // ==============
            UpdateTracker(time);
            this.Regions.Version = VersionCounter;
            this.Regions.Time = time;
            this.m_VersionCnt = VersionCounter;

        }

        /// <summary>
        /// Counterpart of <see cref="BackupTimeLevel(int)"/>, using the full backup.
        /// </summary>
        /// <param name="fullBackup"></param>
        /// <remarks>
        /// Used for **mesh redistribution**  (MPI load balancing)
        /// i.e. when the mesh remains essentially constant, but the (global and local) indexing of cells changes.
        /// </remarks>
        public void ReplaceCurrentTimeLevel(TrackerBackup fullBackup) {
            SinglePhaseField[] LevSet = fullBackup.LevelSets;
            ushort[] RegionCode = fullBackup.Regions;
            (int iLevSet, int iFace)[][] LevSetCoincidingFaces = fullBackup.LevSetCoincidingFaces;
            int VersionCounter = fullBackup.Version;
            double time = fullBackup.time;

            if(LevSet.Length != this.NoOfLevelSets)
                throw new ArgumentOutOfRangeException();
            int NoOfLevelSet = this.NoOfLevelSets;
            if(RegionCode.Length != this.GridDat.Cells.NoOfLocalUpdatedCells)
                throw new ArgumentOutOfRangeException();

            // invalidate everything we got so far
            // ===================================

            this.Regions.InvalidateCaches();
            for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                this.DataHistories[iLs].Current.ClearCaches();
            }

            m_QuadFactoryHelpersHistory.Current.Clear();

            m_XDGSpaceMetricsHistory.Current.Clear();

            // set level-set data
            // ==================
                        
            for(int iLs = 0; iLs < NoOfLevelSets; iLs++) {
                LevelSet Ls = (LevelSet)(this.LevelSets[iLs]);
                Ls.Clear();
                Ls.Acc(1.0, LevSet[iLs]);
            }

            // update region code
            // ==================

            Array.Copy(RegionCode, 0, this.Regions.m_LevSetRegions, 0, RegionCode.Length);
            if(LevSetCoincidingFaces.Any(e => e != null && e.Length > 0)) {
                this.Regions.m_LevSetCoincidingFaces = new (int iLevSet, int iFace)[RegionCode.Length][];
                Array.Copy(LevSetCoincidingFaces, 0, this.Regions.m_LevSetCoincidingFaces, 0, RegionCode.Length);
            } else {
                this.Regions.m_LevSetCoincidingFaces = null;
            }
            this.Regions.Version = VersionCounter;
            this.m_VersionCnt = VersionCounter;
            this.Regions.Time = time;
            MPIUpdate(this.Regions.m_LevSetRegions, this.GridDat);
            this.Regions.Recalc_LenToNextchange();
        }
        
        /// <summary>
        /// Must be called after changing the level-set;
        /// Invoking this method updates the state of cells (e.g. cut, -near, +near, etc.), <see cref="Regions"/>.
        /// After this update, for every <see cref="XDGField"/> the method <see cref="XDGField.OnNext"/> must be invoked;
        /// </summary>
        /// <param name="__NearRegionWith">
        /// new width of near region; if smaller than 0, unchanged
        /// </param>
        /// <param name="__LevSetAllowedMovement">
        /// new values for the allowed level-set movement.
        /// If this value is set to a number higher than the <see cref="NearRegionWidth"/>, the CFL test for the 
        /// corresponding 
        /// level set is omitted. 
        /// </param>
        /// <param name="incremental">
        /// If true, the distance of near-cells can only increase by one.
        /// E.g., if the level-set leaves the domain, setting this parameter to true ensures that the respective 
        /// boundary cells are treated as near-cells for one 
        /// <see cref="UpdateTracker(double, int, bool, int[])"/>-cycle.
        /// This maintains a 'correct' near field when level-set leaves domain, which is necessary e.g. for operator matrix update.
        /// 
        /// Furthermore, it speeds up the detection of cut cells: the detection of cut cells is limited to the near-region,
        /// which drastically reduces the amount of level-set evaluations.
        /// 
        /// Also, detection of topology changes/collisions *require* incremental update set to true.
        /// </param>
        /// <param name="PhysTime">
        /// physical time associated with current Level Set state, see <see cref="LevelSetRegions.Time"/>
        /// </param>
        public void UpdateTracker(double PhysTime, int __NearRegionWith = -1, bool incremental = false, params int[] __LevSetAllowedMovement) {
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
                int JA = m_gDat.Cells.Count;
                int D = m_gDat.Grid.SpatialDimension;
                var smplx = m_gDat.Grid.RefElements;
                //int NoOfSmplxVertice = smplx.NoOfVertices;
                var Krefs = m_gDat.Grid.RefElements;
                int[] NoOfSmplxVertice = Krefs.Select(Kref => Kref.NoOfVertices).ToArray();
                int NoOfLevSets = this.NoOfLevelSets;
                int[][] VerticeInd = m_gDat.Cells.CellVertices;

                Regions.Version = m_VersionCnt;
                Regions.Time = PhysTime;

                ushort[] VertexMarker, LevSetRegions, LevSetRegionsUnsigned;
                BitArray[] LevSetNeg;
                ushort[,] VertexMarkerExternal;
                (int iLevSet, int iFace)[][] _LevSetCoincidingFaces = null;

                


                // init & first time calls
                // =======================
                #region UpdateTracker_INIT
                using (new BlockTrace("INIT", tr)) {

                    // init memory
                    // -----------
                    Regions.m_LevSetRegions_b4Update = Regions.m_LevSetRegions;
                    Regions.m_LevSetRegions = new ushort[JA];
                    LevSetRegions = Regions.m_LevSetRegions;
                    Regions.InvalidateCaches();
                    LevSetRegionsUnsigned = new ushort[JA];
                    Regions.m_LevSetCoincidingFaces = null;

                    // necessary to avoid a 'LevelSetCFLException' on the first
                    // call to this method
                    for(int i = 0; i < JA; i++) {
                        LevSetRegions[i] = AllCut;
                    }

                    // initialize everything to FAR = +7
                    for(int i = 0; i < JA; i++)
                        LevSetRegionsUnsigned[i] = AllFARplus;

                    VertexMarker = new ushort[m_gDat.Vertices.Count];
                    for(int i = VertexMarker.Length - 1; i >= 0; i--)
                        VertexMarker[i] = AllFARplus;
                    VertexMarkerExternal = new ushort[JA - J, Krefs.Max(Kref => Kref.NoOfVertices)];

                    LevSetNeg = new BitArray[NoOfLevSets]; // true marks a cell in which the level set field is completely negative
                    for(int i = 0; i < LevSetNeg.Length; i++)
                        LevSetNeg[i] = new BitArray(J);

                    // define test vertices
                    DefineTestNodes();

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
                    // cell sweep: find cut cells
                    // ==========================

                    CellMask SearchMask;
                    var TempCutCellsBitmaskS = NoOfLevSets.ForLoop(iLs => new BitArray(J));
                    if (incremental) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // incremental update: use the old Near-Band as a search-mask
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
                        // ++++++++++++++++++++++++++++++++++++
                        // Full update: search in entire domain
                        // ++++++++++++++++++++++++++++++++++++
                        SearchMask = CellMask.GetFullMask(m_gDat);
                    }

                    int MaxVecLen = (int)Math.Ceiling(16000.0 / ((double)TestNodes.Max(ns => ns.NoOfNodes)));
                    MaxVecLen = Math.Max(1, MaxVecLen);
                    var eps = BLAS.MachineEps*10;


                    var Gchnks = SearchMask.GetGeometricCellChunks(MaxVecLen, CellInfo.RefElementIndex_Mask | CellInfo.CellType_Mask);
                    foreach (var t_j0_Len in Gchnks) { // loop over all cells in the search mask...
                        int j = t_j0_Len.Item1;
                        int VecLen = t_j0_Len.Item2;


                        int iKref = m_gDat.Cells.GetRefElementIndex(j);
                        var Kref = Krefs[iKref];
                        int noOfFaces = Kref.NoOfFaces;
                        int NoOfNodes = TestNodes[iKref].NoOfNodes;
                        int[] _TestNodesPerFace = this.TestNodesPerFace[iKref];
                        var quadWeights = TestNodes_QuadWeights[iKref];

                        // loop over level sets ...
                        for (int levSetind = NoOfLevSets - 1; levSetind >= 0; levSetind--) {
                            var TempCutCellsBitmask = TempCutCellsBitmaskS[levSetind];


                            if (this.m_DataHistories[levSetind].Current.LevelSet is LevelSet ls) {
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // Use the accelerated bernstein cut cell finding technique for dg levelsets
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                var data = this.m_DataHistories[levSetind].Current;
                                NodeSet EdgeNodes = new NodeSet(Kref, this.TestNodes[iKref].ExtractSubArrayShallow(new int[] {0 , 0}, new int[] { _TestNodesPerFace.Sum() - 1, D - 1}), false); // only use edge nodes
                                MultidimensionalArray levSetVal = data.GetLevSetValues(EdgeNodes, j, VecLen);

                                // loop over all cells in this chunk
                                for (int jj = j; jj < j + VecLen; jj++) {
                                    double[] modalVals = ls.Coordinates.GetRow(jj);
                                    var TransformerEdges = this.TestTransformerEdges[levSetind];
                                    bool Pos = false;
                                    bool Neg = false;

                                    
                                    #region edges  
                                    // loop over nodes on edges...
                                    int nodeIndex = 0;
                                    for (int e = 0; e < noOfFaces; e++) {
                                        bool PosEdge = false;
                                        bool NegEdge = false;

                                        double quadResult = 0.0;
                                        for (int k = 0; k < _TestNodesPerFace[e]; k++) {
                                            double v = levSetVal[jj - j, nodeIndex];

                                            if (v < 0) {
                                                NegEdge = true;
                                            } else if (v > 0) {
                                                PosEdge = true;
                                            }

                                            quadResult += v * v * quadWeights[k]; // weight might not even be necessary to test only for positivity

                                            nodeIndex++;
                                        }

                                        Pos |= PosEdge;
                                        Neg |= NegEdge;

                                        // detect an edge which coincides with the zero-level-set
                                        if (quadResult < eps) {
                                            if (_LevSetCoincidingFaces == null)
                                                _LevSetCoincidingFaces = new (int iLevSet, int iFace)[J][];
                                            (levSetind, e).AddToArray(ref _LevSetCoincidingFaces[jj]);
                                        }

                                    } // end of edges loop
                                    #endregion
                                    

                                    /* Seems to be not robust enough... using the "old" procedure for now
                                    // loop over nodes on edges...
                                    double[] bernsteinValsEdges = new double[TransformerEdges.Destination.Polynomials[iKref].Count];
                                    TransformerEdges.Origin2Dest[iKref].MatVecMul(1.0, modalVals, 0.0, bernsteinValsEdges);

                                    for (int e = 0; e < Kref.NoOfFaces; e++) {
                                        bool PosEdge = false;
                                        bool NegEdge = false;

                                        double quadResult = 0.0;
                                        foreach (int k in TransformerEdges.FaceCoefficients[iKref][e]) {
                                            double v = bernsteinValsEdges[k];

                                            if (v < 0) {
                                                NegEdge = true;
                                            } else if (v > 0) {
                                                PosEdge = true;
                                            }

                                            // another option scale v by some appropriate scaling factor, the problem here can be,
                                            // that for very small physical cells a coinciding edge is detected, even though this is not really the case.
                                            quadResult += v * v; //Math.Abs(v); // if all edge control points are zero that edge is necessarily zero as well  
                                        }

                                        Pos |= PosEdge;
                                        Neg |= NegEdge;

                                        // detect an edge which coincides with the zero-level-set
                                        if (quadResult < eps) {
                                            if (_LevSetCoincidingFaces == null)
                                                _LevSetCoincidingFaces = new (int iLevSet, int iFace)[J][];
                                            (levSetind, e).AddToArray(ref _LevSetCoincidingFaces[jj]);
                                        }
                                    } // end of edges loop 
                                    */

                                    // check also with slight offset and inside the cell
                                    var Transformer = this.TestTransformer[levSetind];
                                    double[] bernsteinVals = new double[Transformer.Destination.Polynomials[iKref].Count];
                                    Transformer.Origin2Dest[iKref].MatVecMul(1.0, modalVals, 0.0, bernsteinVals);

                                    double Min = bernsteinVals.Min();
                                    double Max = bernsteinVals.Max();

                                    Pos |= Max > 0;
                                    Neg |= Min < 0;

                                    // either clear sign change or all zeros
                                    if (Pos && Neg || (!Pos && !Neg)) {
                                        // cell jj is cut by level set
                                        // code cell:

                                        
                                        EncodeLevelSetDist(ref LevSetRegionsUnsigned[jj], 0, levSetind);
                                        TempCutCellsBitmask[jj] = true;
                                    }

                                    // detect purely negative cells
                                    if (Max < 0.0) {
                                            LevSetNeg[levSetind][jj] = true;
                                    } else {
                                        LevSetNeg[levSetind][jj] = false;
                                    }
                                    
                                }
                            } else {
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // Use the old cut-cell detection,
                                // which just evaluates the level-set at certain points in the cell.
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                                var data = this.m_DataHistories[levSetind].Current;
                                MultidimensionalArray levSetVal = data.GetLevSetValues(this.TestNodes[iKref], j, VecLen);

                                for (int jj = 0; jj < VecLen; jj++) {
                                    bool Pos = false;
                                    bool Neg = false;

                                    // loop over nodes on edges...
                                    int nodeIndex = 0;
                                    for (int e = 0; e < noOfFaces; e++) {
                                        bool PosEdge = false;
                                        bool NegEdge = false;

                                        double quadResult = 0.0;
                                        for (int k = 0; k < _TestNodesPerFace[e]; k++) {
                                            double v = levSetVal[jj, nodeIndex];

                                            if (v < 0) {
                                                NegEdge = true;
                                            } else if (v > 0) {
                                                PosEdge = true;
                                            }

                                            quadResult += v * v * quadWeights[k]; // weight might not even be necessary to test only for positivity

                                            nodeIndex++;
                                        }

                                        Pos |= PosEdge;
                                        Neg |= NegEdge;

                                        // detect an edge which coincides with the zero-level-set
                                        if (quadResult < eps) {
                                            if (_LevSetCoincidingFaces == null)
                                                _LevSetCoincidingFaces = new (int iLevSet, int iFace)[J][];
                                            (levSetind, e).AddToArray(ref _LevSetCoincidingFaces[j + jj]);
                                        }

                                    } // end of edges loop

                                    // loop over remaining Nodes...
                                    for (; nodeIndex < NoOfNodes; nodeIndex++) {
                                        double v = levSetVal[jj, nodeIndex];

                                        bool PosNode = false;
                                        bool NegNode = false;
                                        if (v < 0) {
                                            NegNode = true;
                                        } else if (v > 0) {
                                            PosNode = true;
                                        }

                                        Pos |= PosNode;
                                        Neg |= NegNode;
                                    }


                                    if ((Pos && Neg) || (!Pos && !Neg)) {
                                        // cell j+jj is cut by level set

                                        // code cell:
                                        EncodeLevelSetDist(ref LevSetRegionsUnsigned[j + jj], 0, levSetind);
                                        TempCutCellsBitmask[j + jj] = true;
                                    }

                                    if (Neg == true && Pos == false) {
                                        LevSetNeg[levSetind][j + jj] = true;
                                    } else {
                                        LevSetNeg[levSetind][j + jj] = false;
                                    }
                                }
                            }
                        }
                        j += VecLen;

                    }
                    

                    if (_LevSetCoincidingFaces != null)
                        Regions.m_LevSetCoincidingFaces = _LevSetCoincidingFaces;
                }

                //test.AddVector()


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
                #endregion

                // find near and far regions
                // =========================
                using (new BlockTrace("NEAR_SWEEPS", tr)) {
                    
                    int[][] Neighbours = m_gDat.iLogicalCells.CellNeighbours;

                    for(int dist = 1; dist <= m_NearRegionWidth; dist++) { // loop over near region bands...

                        long NumberOfChanges = 1;
                        while(NumberOfChanges > 0) { // sweep until we have no more updates
                            NumberOfChanges = 0;

                            // Pre-sweep: ensure correct distance in the vertex markers
                            // --------------------------------------------------------

                            // rem: this pre-sweep is only to handle periodic & MPI-parallel cases
                            // there is certainly a smarter way to integrate the loops

                            // ensures that the distance of each vertex is lower or equal to the cell distance
                            for(int j = 0; j < J; j++) { // sweep over cells;
                                for(int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) { // loop over level sets...
                                    int dJ = DecodeLevelSetDist(LevSetRegionsUnsigned[j], levSetInd);
                                    if(dJ <= dist) {
                                        foreach(int kk in VerticeInd[j]) {
                                            int dVk = DecodeLevelSetDist(VertexMarker[kk], levSetInd);
                                            dVk = Math.Min(dVk, dJ);
                                            EncodeLevelSetDist(ref VertexMarker[kk], dVk, levSetInd);
                                        }
                                    }
                                }
                            }

                            //using(var stw = new StreamWriter("Vertice-" + dist + ".csv")) {
                            //    var vCoords = m_gDat.Vertices.Coordinates;
                            //    for(int k = 0; k < m_gDat.Vertices.Count; k++) {
                            //        stw.Write(vCoords[k, 0]);
                            //        stw.Write(" ");
                            //        stw.Write(vCoords[k, 1]);
                            //        stw.Write(" ");
                            //        double dVk = DecodeLevelSetDist(VertexMarker[k], 0);
                            //        stw.Write(dVk);
                            //        stw.WriteLine();
                            //    }
                            //}


                            // MPI update
                            // ----------
                            MPIUpdateVertex(VertexMarker, VertexMarkerExternal, m_gDat, NoOfLevSets);

                            // Main sweep: set cells and neighbors
                            // -----------------------------------
                            for(int j = 0; j < J; j++) { // sweep over cells; in this sweep, we are setting band 'dist'
                                int[] _VerticeInd = VerticeInd[j];
                                int _NoOfSmplxVertice = _VerticeInd.Length;
                                Debug.Assert(_NoOfSmplxVertice == m_gDat.Cells.GetRefElement(j).NoOfVertices);

                                //for (int k = 0; k < NoOfSmplxVertice; k++)
                                //    vtxMarkers[k] = VertexMarker[VerticeInd[j, k]];
                                ushort[] __PrevLevSetRegions = null;
                                if(incremental)
                                    __PrevLevSetRegions = RegionsHistory[0].m_LevSetRegions;

                                for(int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) { // loop over level sets...

                                    int dJ = DecodeLevelSetDist(LevSetRegionsUnsigned[j], levSetInd); // distance of 'cell j'
                                    if(dJ >= dist) {
                                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        // cell has not been assigned to band 'dist' or closer yet
                                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                        int[] jNeighs = Neighbours[j];


                                        // vertex-to-cell propagation; find minimum of: { distance of all vertices, cell distance so far }
                                        int mindv = int.MaxValue;
                                        for(int k = 0; k < _NoOfSmplxVertice; k++)
                                            mindv = Math.Min(mindv, DecodeLevelSetDist(VertexMarker[_VerticeInd[k]], levSetInd));
                                        // cell-to-cell propagation; find minimum of: { distance of neighbor cells, cell distance so far }
                                        foreach(int jN in jNeighs) {
                                            int distN = DecodeLevelSetDist(LevSetRegionsUnsigned[jN], levSetInd);
                                            mindv = Math.Min(mindv, distN);
                                        }
                                        if(incremental) {
                                            // The following line ensures that when the level-set is leaving the computational domain,
                                            // we still have a near-cell on the boundary.
                                            int prevDist = Math.Abs(DecodeLevelSetDist(__PrevLevSetRegions[j], levSetInd));
                                            mindv = Math.Min(mindv, prevDist);
                                        }
                                        mindv++; // increase layer

                                        if(mindv == dist) {
                                            EncodeLevelSetDist(ref LevSetRegionsUnsigned[j], mindv, levSetInd);
                                            if(mindv != dJ)
                                                NumberOfChanges++;

                                            // propagate back form cells to vertices
                                            for(int k = 0; k < _NoOfSmplxVertice; k++) {
                                                int dVk = DecodeLevelSetDist(VertexMarker[_VerticeInd[k]], levSetInd);
                                                dVk = Math.Min(dVk, mindv);
                                                EncodeLevelSetDist(ref VertexMarker[_VerticeInd[k]], dVk, levSetInd);
                                            }

                                            // treatment of neighbor cells: robustness against 
                                            //  - hanging nodes (for e.g. 3:1 cell neighborship at an edge) and 
                                            //  - periodic boundary conditions
                                            foreach(int jN in jNeighs) {
                                                int distN = DecodeLevelSetDist(LevSetRegionsUnsigned[jN], levSetInd);
                                                if(distN > dist + 1) {
                                                    EncodeLevelSetDist(ref LevSetRegionsUnsigned[jN], dist + 1, levSetInd);

                                                    if(jN < J) {
                                                        foreach(int kk in VerticeInd[jN]) {
                                                            int dVk = DecodeLevelSetDist(VertexMarker[kk], levSetInd);
                                                            dVk = Math.Min(dVk, dist + 1);
                                                            EncodeLevelSetDist(ref VertexMarker[kk], dVk, levSetInd);
                                                        }
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }
                            } // end of cell loop

                            // MPI update
                            // ----------
                            MPIUpdate(LevSetRegionsUnsigned, m_gDat);
                            MPIUpdateVertex(VertexMarker, VertexMarkerExternal, m_gDat, NoOfLevSets);



                            NumberOfChanges = NumberOfChanges.MPISum();
                        }
                        
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



                // forget values that are not correct anymore
                // ==========================================
                this.m_QuadFactoryHelpersHistory.Current.Clear();
                this.m_XDGSpaceMetricsHistory.Current.Clear();

                // Check Level-Set Topology
                // ========================
                var TopologyProblems = new List<(int iLevSet, int j, int Neigh, int dist_j, int dist_neigh)>();
                using (new BlockTrace("TOPOLOGY_CHECK", tr)) {
                    int[][] Neighbours = m_gDat.iLogicalCells.CellNeighbours;


                    for (int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) {
                        for (int j = 0; j < J; j++) {
                            int dist = DecodeLevelSetDist(LevSetRegions[j], levSetInd);
                            //Console.WriteLine($"dist[{j}] = {dist}");
                            int distSign = Math.Sign(dist);

                            if (Math.Abs(distSign) > 0) {
                                int[] Neighs = Neighbours[j];
                                foreach (var jNeigh in Neighs) {
                                    int distNeigh = DecodeLevelSetDist(LevSetRegions[jNeigh], levSetInd); //lsTrk.Regions.GetLevelSetDistance(0, jNeigh);

                                    int distNeighSign = Math.Sign(distNeigh);

                                    if (Math.Abs(distNeighSign) != 0 && distNeighSign != distSign) {
                                        TopologyProblems.Add((levSetInd, j, jNeigh, dist, distNeigh));
                                        //Console.WriteLine($"Topology error in Cell {j}; contact of purly positive/negative domain across an edge without a cut cell in between. distance of cell {j} is {dist}, distance of neighbour cell {jNeigh} is {distNeigh}.");
                                    }
                                }
                            }
                        }
                    }

                   

                }

                // check the LevelSet CFL
                // ======================

                bool throwCFL;
                int[] fail;
                using (new BlockTrace("CFL_CHECK", tr)) {

                    fail = new int[NoOfLevSets];
                    if (this.PopulatedHistoryLength > 0) {
                        m_LevSetAllowedMovement = __LevSetAllowedMovement;

                        // cannot be moved down because we need the OLD subgrid
                        throwCFL = false;
                        for (int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) {

                            if (m_LevSetAllowedMovement[levSetInd] <= m_NearRegionWidth)
                                fail[levSetInd] = CheckLevelSetCFL(levSetInd);
                            if (fail[levSetInd] > 0)
                                throwCFL = true;

                        }
                    } else {
                        throwCFL = false;
                    }
                }


                // update memory of XDG fields, etc.
                // =================================
                using (new BlockTrace("ObserverUpdate", tr)) {
                    ObserverUpdate(null);
                }

                // throw exceptions, if levelset topology/CFL violated
                // ===================================================

                //
                // Ein schöner Gruß von Florian
                // an jeden, der nochmal daran denken sollte die `LevelSetCFLException` auszukommentieren
                // und das Ganze dann in den Hauptzweig pusht: 
                // Du hast einen Freiflug gewonnen, vom Dach des Maschinenbau-Gebäudes!
                // Herzlichen Glückwunsch!
                //
                // A big Salute
                // from Florian to everyone who should dare to comment out the `LevelSetCFLException` 
                // and push it into the main branch: 
                // You have won a free flight, from the roof of the mechanical engineering building! 
                // Congratulations!
                //


                if (TopologyProblems.Count > 0) {
                    LevelSetTopologyException exception = new LevelSetTopologyException(TopologyProblems);
                    foreach (var reference in m_Observers) {
                        IObserver<LevelSetRegions> observer = reference.Target;
                        if (observer != null) {
                            observer.OnError(exception);
                        }
                    }
                    throw exception;
                }

                if (throwCFL) {
                    LevelSetCFLException exception = new LevelSetCFLException(fail);
                    foreach(var reference in m_Observers) {
                        IObserver<LevelSetRegions> observer = reference.Target;
                        if(observer != null) {
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
            this.TestNodes_QuadWeights = null;
            this.TestNodesPerFace = null;
            this.m_LevelSets = null;
        }

        /// <summary>
        /// Calls the <see cref="IObserver{LevelSetRegions}.OnNext(LevelSetRegions)"/> for all observers.
        /// </summary>
        private void ObserverUpdate(BehaveUnder_LevSetMoovement? updateBehaveOverride) {
            int rnk = ilPSP.Environment.MPIEnv.MPI_Rank;
            int sz = this.GridDat.MpiSize;

            // Remove obsolete observers from list...
            // ======================================

            // MPI synchronization of observers...

            var ObserversRefs = new List<IObserver<LevelSetRegions>>();
            {
                int NoObservers = m_Observers.Count;

                int[] checkNoObservers = (new[] { NoObservers, -NoObservers }).MPIMax();
                if (NoObservers != checkNoObservers[0] || NoObservers != -checkNoObservers[1])
                    throw new ApplicationException("MPI parallelization bug: number of observers of LevelSetTracker is not equal among MPI processors.");

                int[] AliveObservers = new int[NoObservers];
                for (int i = 0; i < NoObservers; i++) {
                    ObserversRefs.Add(m_Observers[i].Target); // having a local ref. prevents the GC of collecting the object while we do MPI sync.
                    AliveObservers[i] = ObserversRefs[i] != null ? 0xff : 0;
                }

                // if observer was killed on any rank, we must also kill it locally!
                int[] GlobalAliveObservers = AliveObservers.MPIMin();

                // observed fields must be synchronized over MPI processors, otherwise deadlocks may occur

                int ii = 0;
                for (int i = 0; i < NoObservers; i++) {
                    if (GlobalAliveObservers[ii] <= 0) {
                        m_Observers.RemoveAt(i);
                        ObserversRefs.RemoveAt(i);
                        i--;
                        NoObservers--;
                    }
                    ii++;
                }
            } //*/


            // update memory of all registered fields
            // =====================================
            // a disadvantage of this notification-by-weak-ref -- approach
            // is that the 'UpdateMemory' may be called also for
            // objects that are already unused but not yet collected...
            // A solution would be to call GC.Collect(), but it is not known
            // whether a GC run or update of unused memory is more expensive.


            // call the update method of all active fields
            foreach (var t in ObserversRefs) {

                BehaveUnder_LevSetMoovement bkup = BehaveUnder_LevSetMoovement.JustReallocate;
                XDGField xDG = t as XDGField;
                if(xDG != null && updateBehaveOverride != null) {
                    bkup = xDG.UpdateBehaviour;
                    xDG.UpdateBehaviour = updateBehaveOverride.Value;
                }

                t.OnNext(Regions);

                if(xDG != null && updateBehaveOverride != null) {
                    xDG.UpdateBehaviour = bkup;
                }
            }
            
            
        }

        /// <summary>
        /// Dirty Hack to support dynamic load balancing,
        /// Calls the <see cref="IObserver{LevelSetRegions}.OnNext"/> for all observers
        /// </summary>
        public void ObserverHack() {
            using (new FuncTrace()) {
                this.ObserverUpdate(BehaveUnder_LevSetMoovement.JustReallocate);
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
            int JA = grdDat.Cells.Count;


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
                            for (int l = 0; l < Len; l++) { // loop over cells in send-buffer for 'pDest'...
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
        /// <see cref="BoSSS.Platform.WeakReference{T}"/>) so that garbage collection will not
        /// be affected by the subscription.
        /// </remarks>
        public IDisposable Subscribe(IObserver<LevelSetRegions> observer) {
            MPICollectiveWatchDog.Watch();
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
