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
    public partial class LevelSetTracker : IObservable<LevelSetTracker.LevelSetRegionsInfo>, IDisposable {

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
        /// <param name="Context"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        public LevelSetTracker(GridData Context, int __NearRegionWidth, string[] _SpeciesTable, ILevelSet levSet1) {
            ConstructorCommon(Context, __NearRegionWidth, _SpeciesTable, levSet1);
        }

        /// <summary>
        /// Creates a level set tracker for just one level set.
        /// </summary>
        /// <param name="Context"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        /// <param name="BruteForceDivisions"></param>
        /// <param name="BruteForceOrder"></param>
        public LevelSetTracker(GridData Context, int __NearRegionWidth, int BruteForceDivisions, int BruteForceOrder, string[] _SpeciesTable, ILevelSet levSet1) {
            ConstructorCommon(Context, __NearRegionWidth, BruteForceOrder, BruteForceDivisions, _SpeciesTable, levSet1);
        }

        /// <summary>
        /// Creates a level set tracker for two level sets.
        /// </summary>
        /// <param name="Context"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="levSet2"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        public LevelSetTracker(GridData Context, int __NearRegionWidth, string[,] _SpeciesTable, ILevelSet levSet1, ILevelSet levSet2) {
            ConstructorCommon(Context, __NearRegionWidth, _SpeciesTable, levSet1, levSet2);
        }

        /// <summary>
        /// Creates a level set tracker for three level sets.
        /// </summary>
        /// <param name="Context"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="levSet2"></param>
        /// <param name="levSet3"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        public LevelSetTracker(GridData Context, int __NearRegionWidth, string[,,] _SpeciesTable, ILevelSet levSet1, ILevelSet levSet2, ILevelSet levSet3) {
            ConstructorCommon(Context, __NearRegionWidth, SpeciesTable, levSet1, levSet2, levSet3);
        }

        /// <summary>
        /// Creates a level set tracker for four level sets.
        /// </summary>
        /// <param name="Context"></param>
        /// <param name="_SpeciesTable">The species table, see <see cref="SpeciesTable"/>;</param>
        /// <param name="levSet1"></param>
        /// <param name="levSet2"></param>
        /// <param name="levSet3"></param>
        /// <param name="levSet4"></param>
        /// <param name="__NearRegionWidth">
        /// width of near region, in number of cells
        /// </param>
        public LevelSetTracker(GridData Context, int __NearRegionWidth, string[,,,] _SpeciesTable, ILevelSet levSet1, ILevelSet levSet2, ILevelSet levSet3, ILevelSet levSet4) {
            ConstructorCommon(Context, __NearRegionWidth, _SpeciesTable, levSet1, levSet2, levSet3, levSet4);
        }

        /// <summary>
        /// Implementation of the constructor with default values for the brute
        /// force quadrature (see <see cref="SetBruteForceQuadratureRules"/>).
        /// </summary>
        /// <param name="Context"></param>
        /// <param name="__NearRegionWidth"></param>
        /// <param name="SpeciesTable"></param>
        /// <param name="levSets"></param>
        private LevelSetTracker(GridData Context, int __NearRegionWidth, Array SpeciesTable, params ILevelSet[] levSets) {
            ConstructorCommon(Context, __NearRegionWidth, 5, 2, SpeciesTable, levSets);
        }

        /// <summary>
        /// Implementation of the constructor with default values for the brute
        /// force quadrature (see <see cref="SetBruteForceQuadratureRules"/>).
        /// </summary>
        /// <param name="Context"></param>
        /// <param name="__NearRegionWidth"></param>
        /// <param name="SpeciesTable"></param>
        /// <param name="levSets"></param>
        private void ConstructorCommon(GridData Context, int __NearRegionWidth, Array SpeciesTable, params ILevelSet[] levSets) {
            ConstructorCommon(Context, __NearRegionWidth, 5, 2, SpeciesTable, levSets);
        }

        /// <summary>
        /// Implementation of the constructor;
        /// </summary>
        private void ConstructorCommon(GridData griData, int __NearRegionWidth, int BruteForceDivisions, int BruteForceOrder, Array SpeciesTable, params ILevelSet[] levSets) {
            // check args, init members
            // ========================
            m_LevelSets = levSets;
            m_gDat = griData;
            if (__NearRegionWidth < 0 || __NearRegionWidth > 6)
                throw new ArgumentException("near region width must be between 0 (including) and 7 (excluding)", "__NearRegionWidth");
            m_NearRegionWidth = __NearRegionWidth;
            
            m_LevelSetValc = new LevSetValueCache[m_LevelSets.Length];
            m_LevelSetGradientsCache = new LevelSetGradientCache[m_LevelSets.Length];
            m_LevelSetNormalsCache = new LevelSetNormalsCache[m_LevelSets.Length];
            m_LevelSetReferenceGradientsCache = new LevelSetReferenceGradientCache[m_LevelSets.Length];
            m_LevelSetReferenceNormalsCache = new LevelSetReferenceNormalsCache[m_LevelSets.Length];
            m_LevelSetReferenceCurvatureCache = new LevelSetReferenceCurvatureCache[m_LevelSets.Length];
            m_LevelSetReferenceHessianCache = new LevelSetReferenceHessianCache[m_LevelSets.Length];
            for (int i = 0; i < m_LevelSets.Length; i++) {
                m_LevelSetValc[i] = new LevSetValueCache(m_LevelSets[i], griData);
                m_LevelSetGradientsCache[i] = new LevelSetGradientCache(m_LevelSets[i], griData);
                m_LevelSetNormalsCache[i] = new LevelSetNormalsCache(this, i);
                m_LevelSetReferenceGradientsCache[i] = new LevelSetReferenceGradientCache(m_LevelSets[i], this);
                m_LevelSetReferenceNormalsCache[i] = new LevelSetReferenceNormalsCache(this, i);
                m_LevelSetReferenceCurvatureCache[i] = new LevelSetReferenceCurvatureCache(this, i);
                m_LevelSetReferenceHessianCache[i] = new LevelSetReferenceHessianCache(this, i);
            }

            // 
            //foreach (SpeciesConfig sc in speciesconf) {
            //    if (m_SpeciesNames.Contains(sc.SpeciesName))
            //        throw new ArgumentException("species named \"" + sc.SpeciesName + "\" occures twice.", "speciesconf");
            //    m_SpeciesNames.Add(sc.SpeciesName);

            //    foreach (bool[] state in sc.LevelSetState) {
            //        if (state.Length != levSets.Length)
            //            throw new ArgumentException("illegal configuration.", "speciesconf");
            //    }
            //}

            if (SpeciesTable.Rank != levSets.Length)
                throw new ArgumentException("rank of species table must match number of level sets.", "SpeciesTable");
            m_SpeciesTable = SpeciesTable;
            SpeciesId invalid;
            invalid.cntnt = int.MinValue;
            ArrayTools.SetAll(m_SpeciesIndex2Id, invalid);
            ArrayTools.SetAll(m_SpeciesId2Index, int.MinValue);
            CollectSpeciesNamesRecursive(0, new int[m_SpeciesTable.Rank]);
            ComputeNoOfSpeciesRecursive(0, new int[4]);
            CollectLevelSetSignCodes();
            ComputeGhostTable();


            // init 
            // ====
            m_LevSetAllowedMovement = new int[levSets.Length];
            ArrayTools.SetAll(m_LevSetAllowedMovement, 1);
            UpdateTracker();

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
            int NoOfLevSets = m_LevelSets.Length;

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
            if (LevSetSigns.Length != m_LevelSets.Length)
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

        /// <summary>
        /// Information about the current level set sign in each cell.
        /// </summary>
        public LevelSetRegionsInfo _Regions {
            get;
            private set;
        }

        /// <summary>
        /// Information about the previous level set sign in each cell.
        /// </summary>
        public LevelSetRegionsInfo PreviousRegions {
            get;
            private set;
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
            if (levelSetValues.Length != this.m_LevelSets.Length)
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
        /// the level sets
        /// </summary>
        public IList<ILevelSet> LevelSets {
            get {
                return (IList<ILevelSet>)m_LevelSets.Clone();
            }
        }

        GridData m_gDat;

        /// <summary>
        /// sdjk
        /// </summary>
        public GridData GridDat {
            get {
                return m_gDat;
            }
        }

        ILevelSet[] m_LevelSets;

        LevSetValueCache[] m_LevelSetValc;

        LevelSetGradientCache[] m_LevelSetGradientsCache;

        LevelSetNormalsCache[] m_LevelSetNormalsCache;

        LevelSetReferenceGradientCache[] m_LevelSetReferenceGradientsCache;

        LevelSetReferenceNormalsCache[] m_LevelSetReferenceNormalsCache;

        LevelSetReferenceCurvatureCache[] m_LevelSetReferenceCurvatureCache;

        LevelSetReferenceHessianCache[] m_LevelSetReferenceHessianCache;

        /// <summary>
        /// returns the (possible) number of species in cell <paramref name="j"/>;
        /// </summary>
        /// <param name="j">
        /// local cell index;
        /// </param>
        /// <param name="ReducedRegionCode">
        /// on exit, the reduced region code for cell <paramref name="j"/>
        /// in an 3-adic representation; This number 
        /// is later on required as an input for <see cref="GetSpeciesIndex(SpeciesId,int)"/>;
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// Here, three states for each of the for level sets are considered:
        /// <list type="bullet">
        ///   <item>positive far (FAR+)</item>
        ///   <item>negative far (FAR-)</item>
        ///   <item>positive near, negative near or cut</item>
        /// </list>
        /// This implies, that also for cells in the near region, memory is allocated for more than one
        /// species.
        /// </remarks>
        public int GetNoOfSpecies(int j, out ReducedRegionCode ReducedRegionCode) {
            ushort celJ = _Regions.m_LevSetRegions[j];
            return GetNoOfSpeciesByRegionCode(celJ, out ReducedRegionCode);
        }

        /// <summary>
        /// returns the (possible) number of species in cell <paramref name="j"/>;
        /// </summary>
        /// <param name="j">
        /// local cell index;
        /// </param>
        /// <returns></returns>
        public int GetNoOfSpecies(int j) {
            ReducedRegionCode dummy;
            int NoOf = GetNoOfSpecies(j, out dummy);
            return NoOf;
        }

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
        ///   <item>positive near, negative near or cutted</item>
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
        /// speciec for each possible combination of level sets; <br/>
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
        /// There are 3 states for each levelset in each cell:
        /// Completly positive, Completly negative or cutted;
        /// Considering at max. 4 Levelsets, this leads to 3<sup>4</sup> = 81 different states.<br/>
        /// For the sign of all 4 level sets, there are 2<sup>4</sup> = 16 different states.<br/>
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
        ///   <item><b>0</b>: the cell is cutted or in the positive or negative near region;</item>
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

        ///// <summary>
        ///// Encodes the signs of the level set values into a byte code
        ///// </summary>
        ///// <param name="levelSetValues">
        ///// The values of the level sets at a given point
        ///// </param>
        ///// <returns>
        ///// A byte code where the n-th bit (starting from the least significant
        ///// bit) is 0 if <paramref name="levelSetValues"/>[n] is less than zero
        ///// and 1 otherwise.
        ///// </returns>
        //public static int ComputeLevelSetBytecode(params double[] levelSetValues) {
        //    if (levelSetValues.Length > 4) {
        //        throw new ArgumentException("Currently, a maximum of 4 level sets is supported", "levelSetValues");
        //    }

        //    int result = 0;
        //    if (levelSetValues.Length > 0 && levelSetValues[0] >= 0) {
        //        result += 1;
        //    }
        //    if (levelSetValues.Length > 1 && levelSetValues[1] >= 0) {
        //        result += 2;
        //    }
        //    if (levelSetValues.Length > 2 && levelSetValues[2] >= 0) {
        //        result += 4;
        //    }
        //    if (levelSetValues.Length > 3 && levelSetValues[3] >= 0) {
        //        result += 8;
        //    }
        //    return result;
        //}

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
        /// this function is the inverse to <see cref="GetSpeciesIndex(SpeciesId,int)"/>;
        /// </summary>
        /// <param name="_ReducedRegionCode">
        /// the value returned by the 2nd parameter of <see cref="GetNoOfSpecies"/>;
        /// </param>
        /// <param name="SpeciesIndex"></param>
        /// <returns></returns>
        /// <remarks>
        /// this function is the inverse to <see cref="GetSpeciesIndex(SpeciesId, int)"/>
        /// </remarks>
        /// <param name="jCell">
        /// local cell index.
        /// </param>
        public SpeciesId GetSpeciesIdFromIndex(int jCell, int SpeciesIndex) {
            ReducedRegionCode rrc;
            int NoOfSpc = GetNoOfSpecies(jCell, out rrc);
            if (SpeciesIndex >= NoOfSpc) {
                SpeciesId invalid;
                invalid.cntnt = int.MinValue;
                return invalid;
            } else {
                return GetSpeciesIdFromIndex(rrc, SpeciesIndex);
            }
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

        /// <summary>
        /// Returns the index of species '<paramref name="_SpeciesId"/>' in the cell '<paramref name="jCell"/>'.
        /// </summary>
        /// <param name="jCell">
        /// a local cell index
        /// </param>
        /// <param name="_SpeciesId">
        /// species identification
        /// </param>
        public int GetSpeciesIndex(SpeciesId _SpeciesId, int jCell) {
            ReducedRegionCode rrc;
            int NoOfSpc = this.GetNoOfSpecies(jCell, out rrc);
            return GetSpeciesIndex(rrc, _SpeciesId);
        }


        private void SetSpeciesIndex(int ind, LevelSetSignCode LevelSetSignBytecode, int SpeciesIndex, SpeciesId speciesId, int NoOfSpecies) {
            int ind2 = ind << 4;
            ind2 = ind2 | LevelSetSignBytecode.val;
            m_SpeciesIndex[ind2] = SpeciesIndex;

            if (SpeciesIndex >= 0 && SpeciesIndex < NoOfSpecies)
                m_SpeciesIndex2Id[ind * 16 + SpeciesIndex] = speciesId;
            m_SpeciesId2Index[ind * 16 + speciesId.cntnt - ___SpeciesIDOffest] = SpeciesIndex;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="levSet"></param>
        /// <returns></returns>
        public int GetLevelSetIndex(ILevelSet levSet) {
            return Array.IndexOf<ILevelSet>(m_LevelSets, levSet);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="levSetInd"></param>
        /// <returns></returns>
        public ushort LevelSetBitmask(int levSetInd) {
            if (levSetInd < 0 || levSetInd >= m_LevelSets.Length)
                throw new ArgumentException("unknown level set.", "levset");

            return (ushort)(0xF << (4 * levSetInd));
        }


        /// <summary>
        /// extracts the distance layer index for level set <paramref name="levSetInd"/>
        /// from the code <paramref name="code"/> (see <see cref="m_LevSetRegions"/>).
        /// </summary>
        /// <param name="code"></param>
        /// <param name="levSetInd"></param>
        /// <returns></returns>
        /// <remarks>
        /// For some cell, the level set distance is 0 for all cells
        /// which are cut by the level set.
        /// </remarks>
        public static int DecodeLevelSetDist(ushort code, int levSetInd) {
            return (((int)code >> (4 * levSetInd)) & 0xf) - 8;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dist"></param>
        /// <param name="levSetInd"></param>
        /// <param name="code"></param>
        /// <returns></returns>
        static public void EncodeLevelSetDist(ref ushort code, int dist, int levSetInd) {
            int c = ((dist + 8) << (4 * levSetInd));
            int msk = (0xf << (4 * levSetInd));
            code = (ushort)((code & ~msk) | c);
        }







        /// <summary>
        /// 
        /// </summary>
        class LevSetValueCache : Caching.CacheLogic_CNs {

            internal LevSetValueCache(ILevelSet levSet, GridData gridData)
                : base(gridData) {
                m_LevSet = levSet;
            }

            ILevelSet m_LevSet;

            /// <summary>
            /// 
            /// </summary>
            protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                m_LevSet.Evaluate(j0, Len, N, output);
            }

            /// <summary>
            /// 
            /// </summary>
            protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                return MultidimensionalArray.Create(Len, NS.NoOfNodes);
            }
        }

        /// <summary>
        /// cached evaluation of the level set fields.
        /// </summary>
        /// <param name="levSetInd"></param>
        /// <param name="NS"></param>
        /// <param name="j0">local index of first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns>
        /// values of the level set field at the nodes of the specified node set.
        /// </returns>
        public MultidimensionalArray GetLevSetValues(int levSetInd, NodeSet NS, int j0, int Len) {
            return m_LevelSetValc[levSetInd].GetValue_Cell(NS, j0, Len);
        }

        /// <summary>
        /// Caches the gradient of a level set
        /// </summary>
        private class LevelSetGradientCache : Caching.CacheLogic_CNs {

            /// <summary>
            /// The level set in question
            /// </summary>
            private ILevelSet levelSet;

            /// <summary>
            /// Constructs a value cache for the evaluation of the gradient of
            /// the given level set.
            /// </summary>
            /// <param name="levelSet">
            /// The level set in question
            /// </param>
            /// <param name="nsc">
            /// The node set controller holding the evaluation nodes
            /// </param>
            internal LevelSetGradientCache(ILevelSet levelSet, GridData gridData)
                : base(gridData) //
            {
                this.levelSet = levelSet;
            }

            /// <summary>
            /// <see cref="ILevelSet.EvaluateGradient"/>
            /// </summary>
            /// <param name="NodeSetIndex">
            /// <see cref="ILevelSet.EvaluateGradient"/>
            /// </param>
            /// <param name="j0">
            /// <see cref="ILevelSet.EvaluateGradient"/>
            /// </param>
            /// <param name="Len">
            /// <see cref="ILevelSet.EvaluateGradient"/>
            /// </param>
            /// <param name="output">
            /// <see cref="ILevelSet.EvaluateGradient"/>
            /// </param>
            protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                levelSet.EvaluateGradient(j0, Len, N, output);
            }

            /// <summary>
            /// <see cref="NodeSetController.ByCellValueCache{MultidimensionalArray}.Allocate"/>
            /// </summary>
            protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet N) {
                return MultidimensionalArray.Create(Len, N.NoOfNodes, N.SpatialDimension);
            }


        }

        /// <summary>
        /// Calculates the gradients of the level set in every affected cell
        /// and every point contained the node set identified by
        /// <paramref name="NodeSetIndex"/>. The layout of the resulting array
        /// is equivalent to <see cref="ILevelSet.EvaluateGradient"/>.
        /// </summary>
        /// <param name="levSetInd">The index of the level set</param>
        /// <param name="NS">
        /// Nodes to evaluate at.
        /// </param>
        /// <param name="j0">
        /// The first cell to be evaluated
        /// </param>
        /// <param name="Len">
        /// The number of cell to be evaluated
        /// </param>
        /// <returns></returns>
        public MultidimensionalArray GetLevelSetGradients(int levSetInd, NodeSet NS, int j0, int Len) {
            return m_LevelSetGradientsCache[levSetInd].GetValue_Cell(NS, j0, Len);
        }

        private class LevelSetReferenceGradientCache : Caching.CacheLogic_CNs {

            private ILevelSet levelSet;

            private LevelSetTracker owner;

            internal LevelSetReferenceGradientCache(ILevelSet levelSet, LevelSetTracker owner)
                : base(owner.GridDat) //
            {
                this.levelSet = levelSet;
                this.owner = owner;
            }

            protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                LevelSet levelSetField = levelSet as LevelSet;
                if (levelSetField == null) {
                    ComputeValuesNonField(levelSet, N, j0, Len, output);
                } else {
                    ComputeValuesField(levelSetField, N, j0, Len, output);
                }
            }

            protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                int D = base.GridData.SpatialDimension;
                return MultidimensionalArray.Create(Len, NS.NoOfNodes, D);
            }

            private void ComputeValuesField(LevelSet levelSet, NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                int D = levelSet.Basis.GridDat.SpatialDimension;
                MultidimensionalArray grad = levelSet.Basis.EvaluateGradient(NS);
                int noOfNodes = grad.GetLength(0);

                unsafe
                {
                    fixed (double* pGrad = &grad.Storage[0], pRes = &output.Storage[0])
                    {
                        double* pResCur = pRes;
                        double* pGradCur = pGrad;
                        for (int i = 0; i < Len; i++) {
                            for (int j = 0; j < noOfNodes; j++) {
                                //double norm = 0.0;
                                for (int d = 0; d < D; d++) {
                                    for (int k = 0; k < levelSet.Basis.MinimalLength; k++) {
                                        *(pResCur + d) += levelSet.Coordinates[i + j0, k] * *(pGradCur + grad.Index(j, k, d));
                                    }
                                }

                                pResCur += D;
                            }
                        }
                    }
                }

                //// Reference implementation
                //for (int i = 0; i < Len; i++) {
                //    for (int j = 0; j < noOfNodes; j++) {
                //        for (int d = 0; d < D; d++) {
                //            for (int k = 0; k < levelSetField.Basis.MinimalLength; k++) {
                //                output[i, j, d] += levelSetField.Coordinates[i + j0, k] * grad[j, k, d];
                //            }
                //        }
                //    }
                //}
            }

            private void ComputeValuesNonField(ILevelSet levelSet, NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                int noOfNodes = NS.NoOfNodes;
                int D = NS.SpatialDimension;
                var R = owner.GridDat.Cells.Transformation;
                var JacDet = owner.GridDat.ChefBasis.Scaling;
                var Cells = owner.GridDat.Cells;

                MultidimensionalArray physGradient = owner.GetLevelSetGradients(0, NS, j0, Len);
                for (int i = 0; i < Len; i++) {
                    int jCell = j0 + i;
                    if (Cells.IsCellAffineLinear(jCell)) {
                        double det = JacDet[j0 + i];

                        for (int j = 0; j < noOfNodes; j++) {
                            for (int d = 0; d < D; d++) {
                                double r = 0.0;
                                for (int dd = 0; dd < D; dd++) {
                                    r += R[i + j0, dd, d] * physGradient[i, j, dd];
                                }

                                output[i, j, d] = r / det;
                            }
                        }
                    } else {
                        throw new NotImplementedException("todo: nonlinear cell");
                    }
                }
            }



        }

        public MultidimensionalArray GetLevelSetReferenceGradients(int levSetInd, NodeSet NS, int j0, int Len) {
            return m_LevelSetReferenceGradientsCache[levSetInd].GetValue_Cell(NS, j0, Len);
        }

        private class LevelSetNormalsCache : Caching.CacheLogic_CNs {

            private LevelSetTracker owner;

            private int levelSetIndex;

            public LevelSetNormalsCache(LevelSetTracker owner, int levelSetIndex)
                : base(owner.GridDat) //
            {
                this.owner = owner;
                this.levelSetIndex = levelSetIndex;
            }

            protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                MultidimensionalArray gradient =
                    owner.m_LevelSetGradientsCache[levelSetIndex].GetValue_Cell(NS, j0, Len);
                //gradient.Storage.CopyTo(output.Storage, 0);
                Debug.Assert(gradient.Dimension == 3 && output.Dimension == 3);
                Debug.Assert(gradient.GetLength(1) == output.GetLength(1));
                Debug.Assert(gradient.GetLength(2) == output.GetLength(2));
                Debug.Assert(output.GetLength(0) == Len);
                output.Set(gradient.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, gradient.GetLength(1) - 1, gradient.GetLength(2) - 1 }));

                int N = gradient.GetLength(1);
                int D = gradient.GetLength(2);

                for (int i = 0; i < Len; i++) {
                    for (int j = 0; j < gradient.GetLength(1); j++) {
                        double normOfGradient = 0.0;
                        for (int k = 0; k < D; k++) {
                            normOfGradient += gradient[i, j, k] * gradient[i, j, k];
                        }

                        // Avoid NaN. Normal is zero in this case
                        if (normOfGradient == 0.0) {
                            continue;
                        }

                        double OOnormOfGradient = 1.0 / Math.Sqrt(normOfGradient);

                        for (int k = 0; k < D; k++) {
                            output[i, j, k] *= OOnormOfGradient;
                        }
                    }
                }
            }

            protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                return MultidimensionalArray.Create(Len, NS.NoOfNodes, NS.SpatialDimension);
            }
        }

        /// <summary>
        /// Variant of <see cref="GetLevelSetGradients"/> which normalizes the
        /// gradients before returning them
        /// </summary>
        public MultidimensionalArray GetLevelSetNormals(int levSetInd, NodeSet NS, int j0, int Len) {
            return m_LevelSetNormalsCache[levSetInd].GetValue_Cell(NS, j0, Len);
        }

        private class LevelSetReferenceNormalsCache : Caching.CacheLogic_CNs {

            private LevelSetTracker owner;

            private int levelSetIndex;

            public LevelSetReferenceNormalsCache(LevelSetTracker owner, int levelSetIndex)
                : base(owner.GridDat) //
            {
                this.owner = owner;
                this.levelSetIndex = levelSetIndex;
            }

            protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                MultidimensionalArray gradient =
                    owner.m_LevelSetReferenceGradientsCache[levelSetIndex].GetValue_Cell(NS, j0, Len);
                gradient.Storage.CopyTo(output.Storage, 0);

                for (int i = 0; i < gradient.GetLength(0); i++) {
                    for (int j = 0; j < gradient.GetLength(1); j++) {
                        double normOfGradient = 0.0;
                        for (int d = 0; d < gradient.GetLength(2); d++) {
                            normOfGradient += gradient[i, j, d] * gradient[i, j, d];
                        }

                        // Avoid NaN. Normal is zero in this case
                        if (normOfGradient == 0.0) {
                            continue;
                        }

                        normOfGradient = Math.Sqrt(normOfGradient);

                        for (int d = 0; d < gradient.GetLength(2); d++) {
                            output[i, j, d] = gradient[i, j, d] / normOfGradient;
                        }
                    }
                }
            }

            protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                return MultidimensionalArray.Create(Len, NS.NoOfNodes, NS.SpatialDimension);
            }
        }

        public MultidimensionalArray GetLevelSetReferenceNormals(int levSetInd, NodeSet NS, int j0, int Len) {
            return m_LevelSetReferenceNormalsCache[levSetInd].GetValue_Cell(NS, j0, Len);
        }

        public MultidimensionalArray GetLevelSetNormalReferenceToPhysicalMetrics(int levSetInd, NodeSet NS, int j0, int Len) {
            MultidimensionalArray physGradients = GetLevelSetGradients(levSetInd, NS, j0, Len);
            MultidimensionalArray refGradients = GetLevelSetReferenceGradients(levSetInd, NS, j0, Len);
            int noOfNodes = physGradients.GetLength(1);
            int D = GridDat.Grid.SpatialDimension;
            var OneOverSqrt_AbsJacobiDet = GridDat.ChefBasis.Scaling;
            MultidimensionalArray result = MultidimensionalArray.Create(Len, noOfNodes);

            for (int i = 0; i < Len; i++) {
                int jCell = j0 + i;

                if (GridDat.Cells.IsCellAffineLinear(jCell)) {
                    double sc = OneOverSqrt_AbsJacobiDet[jCell];
                    for (int j = 0; j < noOfNodes; j++) {
                        double normPhys = 0.0;
                        double normRef = 0.0;

                        for (int d = 0; d < D; d++) {
                            normPhys += physGradients[i, j, d] * physGradients[i, j, d];
                            normRef += refGradients[i, j, d] * refGradients[i, j, d];
                        }

                        result[i, j] = Math.Sqrt(normRef / normPhys) * sc;
                    }
                } else {

                    throw new NotImplementedException("nonlinear cell: todo");
                }
            }


            //{
            //    double erracc = 0;

            //    var NormalsRef = this.GetLevelSetReferenceNormals(levSetInd, NodeSetIndex, j0, Len);
            //    var NoramlsPhys = this.GetLevelSetNormals(levSetInd, NodeSetIndex, j0, Len);

            //    MultidimensionalArray Jacobi = MultidimensionalArray.Create(D, D);
            //    double[] v1 = new double[D];


            //    for (int i = 0; i < Len; i++) {
            //        int jCell = i + j0;
            //        for (int d1 = 0; d1 < D; d1++) {
            //            for (int d2 = 0; d2 < D; d2++) {
            //                Jacobi[d1, d2] = Ctx.GridDat.Transformation[jCell, d1 + D*d2];
            //            }
            //        }

            //        for (int iNode = 0; iNode < noOfNodes; iNode++) {

            //            for (int d1 = 0; d1 < D; d1++) {
            //                double acc = 0;
            //                for (int d2 = 0; d2 < D; d2++) {
            //                    acc += Jacobi[d1, d2]*NormalsRef[i, iNode, d2];
            //                }
            //                v1[d1] = acc;
            //            }

            //            double metrix = 0;
            //            for (int d = 0; d < D; d++)
            //                metrix += v1[d]*NoramlsPhys[i, iNode, d];

            //            double anderemetrix = result[i, iNode];

            //            erracc += (metrix - anderemetrix).Pow2();
            //        }
            //    }

            //    Console.WriteLine("metrix test: " + erracc);
            //}




            return result;
        }

        ///// <summary>
        ///// 
        ///// </summary>
        //class LevelSetBitmaskCache : NodeSetController.ByCellValueCache<byte[]> {


        //    internal LevelSetBitmaskCache(LevelSetTracker owner)
        //        : base(owner.m_context.NSC) {
        //        m_owner = owner;
        //    }

        //    LevelSetTracker m_owner;


        //    /// <summary>
        //    /// 
        //    /// </summary>
        //    /// <param name="NodeSetIndex"></param>
        //    /// <param name="j0"></param>
        //    /// <param name="Len"></param>
        //    /// <param name="output"></param>
        //    protected override void ComputeValues(int NodeSetIndex, int j0, int Len, byte[] output) {
        //        Array.Clear(output, 0, output.Length);

        //        for (int i = 0; i < m_owner.m_LevelSets.Length; i++) {
        //            MultidimensionalArray LevSetVal = m_owner.m_LevelSetValc[i].GetValue(NodeSetIndex, j0, Len);

        //            int NoOfNodes = LevSetVal.GetLength(1);
        //            byte PosMask = m_owner.LevelSetPositiveMask(i);
        //            byte NegMask = m_owner.LevelSetPositiveMask(i);

        //            for (int j = 0; j < Len; j++) {
        //                for (int k = 0; k < NoOfNodes; k++) {
        //                    double v = LevSetVal[j, k];

        //                    if (v <= 0) output[j] |= NegMask;
        //                    if (v >= 0) output[j] |= PosMask;
        //                }
        //            }

        //        }
        //    }

        //    /// <summary>
        //    /// 
        //    /// </summary>
        //    /// <param name="Len"></param>
        //    /// <param name="NoOfNodes"></param>
        //    /// <param name="D"></param>
        //    /// <returns></returns>
        //    protected override byte[] Allocate(int Len, int NoOfNodes, int D) {
        //        return new byte[Len];
        //    }
        //}

        ///// <summary>
        ///// used by <see cref="GetLevSetBitMask"/>;
        ///// </summary>
        //LevelSetBitmaskCache m_LevelSetBitmaskCache;

        ///// <summary>
        ///// 
        ///// </summary>
        ///// <param name="NodeSetIndex"></param>
        ///// <param name="j0"></param>
        ///// <param name="Len"></param>
        ///// <returns>
        ///// A bitmask that indicates wether some cell is in the positive domain, 
        ///// in the negative domain, or cut by a level set (see remaks).
        ///// </returns>
        ///// <remarks>
        ///// With respect to some level set, there are three possible states for some cell:
        ///// <list type="bullet">
        /////   <item>(A): the cell lies completely in the positive domain of the level set, or...</item>
        /////   <item>(B): the cell lies completely within the negative domain of the level set, or...</item>
        /////   <item>(C): the cell is cut by the level set.</item>
        ///// </list>
        ///// These states, for each cell <i>j</i> and each  
        ///// level set <i>G</i> in <see cref="LevelSets"/>, are identified by bitmasks: Therefore, let be...
        ///// <list type="bullet">
        /////   <item>
        /////     <i>Gp</i> the positive set bitmask, i.e. <i>Gp</i>=<see cref="LevelSetPositiveMask"/>(<i>G</i>) and 
        /////   </item>
        /////   <item>
        /////     <i>Gn</i> the negative set bitmask, i.e. <i>Gn</i>=<see cref="LevelSetNegativeMask"/>(<i>G</i>).
        /////   </item>
        ///// </list>
        ///// Then, cases (A), (B) and (C) are encoded in the following way:
        ///// <list type="bullet">
        /////   <item>
        /////   (A) is the case if, and only if, (<i>Gp</i> BAND <code>LevelSetBitMasks</code>[<i>j</i>]) is true
        /////   and (<i>Gn</i> BAND <code>LevelSetBitMasks</code>[<i>j</i>]) is false, 
        /////   </item>
        /////   <item>
        /////   (B) is the case if, and only if, (<i>Gn</i> BAND <code>LevelSetBitMasks</code>[<i>j</i>]) is true
        /////   and (<i>Gp</i> BAND <code>LevelSetBitMasks</code>[<i>j</i>]) is false and
        /////   </item>
        /////   <item>
        /////   (C) is the case if, and only if, (<i>Gp</i> BAND <code>LevelSetBitMasks</code>[<i>j</i>]) is true
        /////   and (<i>Gn</i> BAND <code>LevelSetBitMasks</code>[<i>j</i>]) is true. 
        /////   </item>
        ///// </list>
        ///// Here, BAND denotes bitwise AND, and <i>j</i> a local cell index.</remarks>
        //public byte[] GetLevSetBitMask(int NodeSetIndex, int j0, int Len) {
        //    return m_LevelSetBitmaskCache.GetValue(NodeSetIndex, j0, Len);
        //}

        private class LevelSetReferenceHessianCache : Caching.CacheLogic_CNs {

            private LevelSetTracker owner;

            private int levelSetIndex;

            public LevelSetReferenceHessianCache(LevelSetTracker owner, int levelSetIndex)
                : base(owner.GridDat) //
            {
                this.owner = owner;
                this.levelSetIndex = levelSetIndex;
            }



            protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                LevelSet LevSet = (LevelSet)this.owner.LevelSets[this.levelSetIndex];

                var BasisHessian = LevSet.Basis.Evaluate2ndDeriv(NS);

                Debug.Assert(output.GetLength(0) == Len);
                int N = LevSet.Basis.Length;
                Debug.Assert(BasisHessian.GetLength(1) == N);
                int D = this.owner.m_gDat.SpatialDimension;
                Debug.Assert(D == BasisHessian.GetLength(2));
                Debug.Assert(D == BasisHessian.GetLength(3));
                Debug.Assert(D == output.GetLength(2));
                Debug.Assert(D == output.GetLength(3));
                int K = output.GetLength(1); // No of nodes
                Debug.Assert(K == BasisHessian.GetLength(0));

                var Coordinates = ((MultidimensionalArray)LevSet.Coordinates).ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, N - 1 });
                output.Multiply(1.0, BasisHessian, Coordinates, 0.0, "jkdr", "kndr", "jn");
            }

            protected override MultidimensionalArray Allocate(int j0, int Len, NodeSet N) {
                return MultidimensionalArray.Create(Len, N.NoOfNodes, N.SpatialDimension, N.SpatialDimension);
            }
        }

        /// <summary>
        /// the Hessian with respect to reference coordinates, of the level set function 
        /// </summary>
        /// <param name="levSetInd"></param>
        /// <param name="NodeSet"></param>
        /// <param name="j0">first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns></returns>
        public MultidimensionalArray GetLevelSetReferenceHessian(int levSetInd, NodeSet nodes, int j0, int Len) {
            return m_LevelSetReferenceHessianCache[levSetInd].GetValue_Cell(nodes, j0, Len);
        }

        private class LevelSetReferenceCurvatureCache : Caching.CacheLogic_CNs {

            private LevelSetTracker owner;

            private int levelSetIndex;

            public LevelSetReferenceCurvatureCache(LevelSetTracker owner, int levelSetIndex)
                : base(owner.GridDat) //
            {
                this.owner = owner;
                this.levelSetIndex = levelSetIndex;
            }

            protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {

                MultidimensionalArray Phi = owner.GetLevSetValues(this.levelSetIndex, NS, j0, Len);
                MultidimensionalArray GradPhi = owner.GetLevelSetReferenceGradients(this.levelSetIndex, NS, j0, Len);
                MultidimensionalArray HessPhi = owner.GetLevelSetReferenceHessian(this.levelSetIndex, NS, j0, Len);

                MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
                MultidimensionalArray Laplace = new MultidimensionalArray(2);
                MultidimensionalArray Q = new MultidimensionalArray(3);

                int K = output.GetLength(1);
                int D = GradPhi.GetLength(2);
                Debug.Assert(D == this.owner.GridDat.SpatialDimension);

                ooNormGrad.Allocate(Len, K);
                Laplace.Allocate(Len, K);
                Q.Allocate(Len, K, D);


                // compute the monstrous formula
                // -----------------------------

                // norm of Gradient:
                for (int d = 0; d < D; d++) {
                    var GradPhi_d = GradPhi.ExtractSubArrayShallow(-1, -1, d);
                    ooNormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                }
                ooNormGrad.ApplyAll(x => 1.0 / Math.Sqrt(x));

                // laplacian of phi:
                for (int d = 0; d < D; d++) {
                    var HessPhi_d_d = HessPhi.ExtractSubArrayShallow(-1, -1, d, d);
                    Laplace.Acc(1.0, HessPhi_d_d);
                }

                // result = Laplacian(phi)/|Grad phi|
                output.Multiply(1.0, Laplace, ooNormGrad, 0.0, "ik", "ik", "ik");


                // result = Grad(1/|Grad(phi)|)
                for (int d1 = 0; d1 < D; d1++) {
                    var Qd = Q.ExtractSubArrayShallow(-1, -1, d1);

                    for (int d2 = 0; d2 < D; d2++) {
                        var Grad_d2 = GradPhi.ExtractSubArrayShallow(-1, -1, d2);
                        var Hess_d2_d1 = HessPhi.ExtractSubArrayShallow(-1, -1, d2, d1);

                        Qd.Multiply(-1.0, Grad_d2, Hess_d2_d1, 1.0, "ik", "ik", "ik");
                    }
                }

                ooNormGrad.ApplyAll(x => x * x * x);

                output.Multiply(1.0, GradPhi, Q, ooNormGrad, 1.0, "ik", "ikd", "ikd", "ik");
            }

            protected override MultidimensionalArray Allocate(int j0, int Len, NodeSet N) {
                return MultidimensionalArray.Create(Len, N.NoOfNodes);
            }
        }

        public MultidimensionalArray GetLevelSetReferenceCurvature(int levSetInd, NodeSet NS, int j0, int Len) {
            return m_LevelSetReferenceCurvatureCache[levSetInd].GetValue_Cell(NS, j0, Len);
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


            for (int i = 0; i < m_LevelSets.Length; i++) {
                for (int j = 0; j < m_gDat.Cells.NoOfCells; j++) {
                    int iKref = this.GridDat.Cells.GetRefElementIndex(j);

                    MultidimensionalArray levSetVals = this.GetLevSetValues(i, TestNodes[iKref], j, 1);

                    bool containsPos = false;
                    bool containsNeg = false;

                    for (int k = 0; k < levSetVals.GetLength(1); k++) {
                        double v = levSetVals[0, k];

                        if (v <= 0)
                            containsNeg = true;
                        if (v >= 0)
                            containsPos = true;
                    }


                    int dist = LevelSetTracker.DecodeLevelSetDist(_Regions.m_LevSetRegions[j], i);

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
            ushort[] oldCode = this.PreviousRegions.m_LevSetRegions;

            CellMask newCut = this._Regions.GetCutCellSubgrid4LevSet(LevSetIdx).VolumeMask;

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
            using (new FuncTrace()) {
                int Jup = this.GridDat.CellPartitioning.LocalLength;
                int[] tempRegionData = new int[Jup];

                var oldRegion = PreviousRegions.m_LevSetRegions;
                var newRegion = _Regions.m_LevSetRegions;

                for (int j = 0; j < Jup; j++) {
                    ushort a = oldRegion[j];
                    ushort b = newRegion[j];
                    int C = a | ((int)b << 16);
                    tempRegionData[j] = C;
                }
                
#if DEBUG
                RestoreAfterLoadBalance(this.m_VersionCnt, tempRegionData, false);
                int JA = m_gDat.Cells.NoOfCells;
                var restored1 = PreviousRegions.m_LevSetRegions;
                var restored2 = _Regions.m_LevSetRegions;
                for (int j = 0; j < JA; j++) {
                    Debug.Assert(restored1[j] == oldRegion[j]);
                    Debug.Assert(restored2[j] == newRegion[j]);
                }
#endif 

                return tempRegionData;
            }
        }

        /// <summary>
        /// Late stage of dynamic load balancing, restoring data after the grid cells were re-distributed.
        /// </summary>
        public void RestoreAfterLoadBalance(int versionNo, int[] ExchangeData, bool ObUp = true) {
            using (new FuncTrace()) {
                int J = m_gDat.Cells.NoOfLocalUpdatedCells;
                int JA = m_gDat.Cells.NoOfCells;
                Debug.Assert(ExchangeData.Length == J);

                PreviousRegions = new LevelSetRegionsInfo(this) {
                    m_LevSetRegions = new ushort[JA],
                    Version = versionNo - 1
                };
                _Regions = new LevelSetRegionsInfo(this) {
                    m_LevSetRegions = new ushort[JA],
                    Version = versionNo
                };
                ushort[] oldR = PreviousRegions.m_LevSetRegions;
                ushort[] newR = _Regions.m_LevSetRegions;

                for (int j = 0; j < J; j++) {
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
                    oldNearMask = this._Regions.GetNearFieldMask(m_NearRegionWidth);
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
                int NoOfLevSets = m_LevelSets.Length;

                PreviousRegions = _Regions;
                _Regions = new LevelSetRegionsInfo(this);
                _Regions.Version = m_VersionCnt;

                ushort[] VertexMarker, LevSetRegions, LevSetRegionsUnsigned;
                BitArray[] LevSetNeg;
                ushort[,] VertexMarkerExternal;
                #region UpdateTracker_INIT
                using (new BlockTrace("INIT", tr)) {

                    // init memory
                    // ===========
                    _Regions.m_LevSetRegions = new ushort[JA];
                    LevSetRegions = _Regions.m_LevSetRegions;
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

                    LevSetNeg = new BitArray[m_LevelSets.Length]; // true marks a cell in which the level set field is completely 
                                                                  // negative
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
                    foreach (var ch in m_LevelSetValc) {
                        ch.Clear();
                    }
                    foreach (var ch in m_LevelSetGradientsCache) {
                        ch.Clear();
                    }
                    foreach (var ch in m_LevelSetReferenceGradientsCache) {
                        ch.Clear();
                    }
                    foreach (var ch in m_LevelSetNormalsCache) {
                        ch.Clear();
                    }
                    foreach (var ch in m_LevelSetReferenceNormalsCache) {
                        ch.Clear();
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
                            int NoofLevSets = m_LevelSets.Length;

                            ushort[] __PrevLevSetRegions = PreviousRegions.m_LevSetRegions;

                            for (int levSetind = 0; levSetind < NoofLevSets; levSetind++) {
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
                            for (int levSetind = m_LevelSets.Length - 1; levSetind >= 0; levSetind--) {

                                MultidimensionalArray levSetVal = GetLevSetValues(levSetind, this.TestNodes[iKref], j, VecLen);

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
                        for (int levSetind = m_LevelSets.Length - 1; levSetind >= 0; levSetind--) {
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
                        MPIUpdateVertex(VertexMarker, VertexMarkerExternal, m_gDat, m_LevelSets.Length);
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
                                    __PrevLevSetRegions = PreviousRegions.m_LevSetRegions;

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
                            MPIUpdateVertex(VertexMarker, VertexMarkerExternal, m_gDat, m_LevelSets.Length);
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
                    _Regions.Recalc_LenToNextchange();
                }

                if (PreviousRegions == null)
                    PreviousRegions = _Regions;

                // forget values that are not correct anymore
                // ==========================================
                //this.VolFraction = null;
                //this.SpecisArea = null;
                //this.m_TotalVolume = null;
                this.m_QuadFactoryHelpers = new Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper>();

                // check the LevelSet CFL
                // ======================

                bool throwCFL;
                int[] fail;
                using (new BlockTrace("CFL_CHECK", tr)) {

                    {
                        m_LevSetAllowedMovement = __LevSetAllowedMovement;

                        // cannot be moved down because we need the OLD subgrid
                        fail = new int[NoOfLevSets];
                        throwCFL = false;
                        for (int levSetInd = 0; levSetInd < NoOfLevSets; levSetInd++) {

                            if (m_LevSetAllowedMovement[levSetInd] <= m_NearRegionWidth)
                                fail[levSetInd] = CheckLevelSetCFL(levSetInd);
                            if (fail[levSetInd] > 0)
                                throwCFL = true;

                        }
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
                        IObserver<LevelSetRegionsInfo> observer = reference.Target;
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
            this._Regions = null;
            this.PreviousRegions = null;
            this.m_SpeciesNames = null;
            this.m_SpeciesIndex2Id = null;
            this.m_SpeciesIndex = null;
            this.m_gDat = null;
            this.m_GhostTable = null;
            this.m_LevelSetGradientsCache = null;
            this.m_LevelSetNormalsCache = null;
            this.m_LevelSetReferenceCurvatureCache = null;
            this.m_LevelSetReferenceGradientsCache = null;
            this.m_LevelSetReferenceHessianCache = null;
            this.m_LevelSetReferenceNormalsCache = null;
            this.m_LevelSets = null;
            this.m_LevelSetSignCodes = null;
            this.m_LevelSetValc = null;
            this.m_LevSetAllowedMovement = null;
            this.m_NoOfSpecies = null;
            this.m_Observers = null;
            this.m_QuadFactoryHelpers = null;
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
                    IObserver<LevelSetRegionsInfo> observer = reference.Target;
                    if (observer != null) {
                        reference.Target.OnNext(_Regions);
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
        private List<BoSSS.Platform.WeakReference<IObserver<LevelSetRegionsInfo>>> m_Observers =
            new List<BoSSS.Platform.WeakReference<IObserver<LevelSetRegionsInfo>>>();

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
        public IDisposable Subscribe(IObserver<LevelSetRegionsInfo> observer) {
            BoSSS.Platform.WeakReference<IObserver<LevelSetRegionsInfo>> reference =
                new BoSSS.Platform.WeakReference<IObserver<LevelSetRegionsInfo>>(observer);
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
            private BoSSS.Platform.WeakReference<IObserver<LevelSetRegionsInfo>> observer;

            /// <summary>
            /// Creates a new unsubscriber.
            /// </summary>
            /// <param name="owner">
            /// The creator of this object.
            /// </param>
            /// <param name="observer">
            /// The observer that may want to unregister itself.
            /// </param>
            public Unsubscriber(LevelSetTracker owner, BoSSS.Platform.WeakReference<IObserver<LevelSetRegionsInfo>> observer) {
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
                IObserver<LevelSetRegionsInfo> observer = reference.Target;
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
