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
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using MPI.Wrappers;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Foundation.XDG {
    
    
    /// <summary>
    /// common functionality for the framing classes <see cref="BoSSS.Foundation.XDG.XSpatialOperator.SpeciesFrameMatrix{M}"/> and <see cref="BoSSS.Foundation.XDG.XSpatialOperator.SpeciesFrameVector{V}"/>.
    /// </summary>
    public sealed class FrameBase {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="lsTrk">level set tracker</param>
        /// <param name="spcId">species which should be framed</param>
        /// <param name="__FullMap"></param>
        /// <param name="SupportExternal">
        /// true: support also indices related to external/ghost cells
        /// </param>
        public FrameBase(LevelSetTracker.LevelSetRegions regions, SpeciesId spcId, UnsetteledCoordinateMapping __FullMap, bool SupportExternal) {
            m_Regions = regions;
            m_spcId = spcId;
            FrameMap = CreateMap(__FullMap);
            FullMap = __FullMap;
            m_supportExternal = SupportExternal;
                        
            Setup_Frame2Full(); // Frame2Full is almost used every time, so we set it up immediately
            //                     (in contrast, Full2Frame is set up on demand)
        }



        /// <summary>
        /// Frame mapping, i.e. only the DG-basis of the related <see cref="Species"/>.
        /// </summary>
        public UnsetteledCoordinateMapping FrameMap {
            get;
            private set;
        }
                
        /// <summary>
        /// Full mapping, i.e. including all species.
        /// </summary>
        public UnsetteledCoordinateMapping FullMap {
            get;
            private set;
        }

        /// <summary>
        /// Replaces all <see cref="XDGBasis"/>-objects in the mapping <paramref name="full"/>
        /// with their non-cut counterparts (<see cref="XDGBasis.NonX_Basis"/>.
        /// </summary>
        static UnsetteledCoordinateMapping CreateMap(UnsetteledCoordinateMapping full) {
            var basisS = from b in full.BasisS
                         select ((b is XDGBasis) ? ((XDGBasis)b).NonX_Basis : b);
            return new UnsetteledCoordinateMapping(basisS.ToArray<Basis>());
        }

        LevelSetTracker.LevelSetRegions m_Regions;
        SpeciesId m_spcId;

        /// <summary>
        /// regions halt
        /// </summary>
        public LevelSetTracker.LevelSetRegions Regions {
            get {
                return m_Regions;
            }
        }

        /// <summary>
        /// the 'framed' species.
        /// </summary>
        public SpeciesId Species {
            get {
                return m_spcId;
            }
        }

        bool m_supportExternal;
        int[] Frame2Full_Lookup;

        void Setup_Frame2Full() {
            //if (FrameMap.Rank == 0)
            //    Debugger.Launch();
            //csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            //FrameMap.GridDat.Cells.NoOfLocalUpdatedCells


            int L = FrameMap.LocalLength;
            Frame2Full_Lookup = new int[L];


            Basis[] BasisS = FullMap.BasisS.ToArray();
            XDGBasis[] XBasisS = new XDGBasis[BasisS.Length];
            for (int i = 0; i < BasisS.Length; i++)
                XBasisS[i] = BasisS[i] as XDGBasis;




            var GridDat = this.FullMap.GridDat;

            Frame2Full_Lookup = new int[L];
            Frame2Full_Lookup.SetAll(int.MinValue);
            int Jup = FrameMap.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int Jtot = FrameMap.GridDat.iLogicalCells.Count;
            int G = BasisS.Length;
            for (int j = 0; j < Jup; j++) { // loop over all cells ...
                int NoOfSpc = m_Regions.GetNoOfSpecies(j);
                int iSpc = m_Regions.GetSpeciesIndex(m_spcId, j);

                for (int g = 0; g < G; g++) { // loop over all basises
                    Basis b = BasisS[g];
                    XDGBasis xb = XBasisS[g];
                    bool b_is_XDGbasis = (xb != null);

                    int n0, N;
                    if (b_is_XDGbasis) {
                        if (iSpc < 0) {
                            n0 = int.MaxValue;
                            N = int.MinValue;
                        } else {
                            N = xb.DOFperSpeciesPerCell;
                            n0 = N*iSpc;
                        }
                    } else {
                        n0 = 0;
                        N = b.GetLength(j);
                    }

                    if (N >= 0) {
                        int iFull = (int)FullMap.GlobalUniqueCoordinateIndex(g, j, n0);
                        int iLoc = FrameMap.LocalUniqueCoordinateIndex(g, j, 0);

                        for (int n = 0; n < N; n++) {
                            Frame2Full_Lookup[iLoc + n] = iFull + n;
                            Debug.Assert(Frame2Full_Lookup[iLoc + n] < 0 || FullMap.IsInLocalRange(Frame2Full_Lookup[iLoc + n]));
                        }
                    }
                }
            }

#if DEBUG
            var _Frame2Full_Lookup = new int[L];
            for (int iLoc = 0; iLoc < L; iLoc++) { // loop over all indices of the frame...
                int j, g, n;
                FrameMap.LocalFieldCoordinateIndex(iLoc, out g, out j, out n);

                int NoOfSpc = m_Regions.GetNoOfSpecies(j);
                Basis b = BasisS[g];
                XDGBasis xb = XBasisS[g];
                bool b_is_XDGbasis = (xb != null);
                int iSpc = m_Regions.GetSpeciesIndex(m_spcId, j);

                if (b_is_XDGbasis) {
                    if (iSpc < 0) {
                        _Frame2Full_Lookup[iLoc] = int.MinValue;
                    } else {
                        int Nsep = xb.DOFperSpeciesPerCell;
                        _Frame2Full_Lookup[iLoc] =  (int)(FullMap.GlobalUniqueCoordinateIndex(g, j, iSpc * Nsep + n));
                    }
                } else {
                    _Frame2Full_Lookup[iLoc] = (int)(FullMap.GlobalUniqueCoordinateIndex(g, j, n));

                }
                Debug.Assert(_Frame2Full_Lookup[iLoc] == Frame2Full_Lookup[iLoc]);
            }

            for (int l = 0; l < L; l++) {
                Debug.Assert(Frame2Full_Lookup[l] < 0 || FullMap.IsInLocalRange(Frame2Full_Lookup[l]));
            }
#endif


            /*
            if (m_supportExternal) {
                int[] PeerProcess = GridDat.Parallel.ProcessesToReceiveFrom;
                int TT = PeerProcess.Length;
                int DOFperCell = FrameMap.MaxTotalNoOfCoordinatesPerCell;

                ExtRangeStart = new int[TT];
                ExtRangeLen = new int[TT];
                ExtRangeTrafo = new int[TT][];

                for(int tt = 0; tt < TT; tt++) {
                    int otherRank = PeerProcess[tt];
                    int jFristCell = GridDat.Parallel.m_RcvCommListsInsertIndex[otherRank];
                    int Jrank = GridDat.Parallel.m_RcvCommListsNoOfItems[otherRank];

                    
                    int offset = (int)(FrameMap.GlobalUniqueCoordinateIndex(0, jFristCell, 0));
                    ExtRangeStart[tt] = offset;
                    ExtRangeLen[tt] = Jrank*DOFperCell;

                    
                    ExtRangeTrafo[tt] = new int[ExtRangeLen[tt]];
                    ExtRangeTrafo[tt].SetAll(int.MinValue);

                    for (int j = 0; j < Jrank; j++) { // loop over all external cells that this process receives from 'otherRank' ...
                        int jCell = j + jFristCell;
                        Debug.Assert(jCell >= Jup);
                        Debug.Assert(jCell < Jtot);
                        ReducedRegionCode rrc;
                        int NoOfSpc = m_lsTrk.GetNoOfSpecies(jCell, out rrc);
                        int iSpc = m_lsTrk.GetSpeciesIndexFromId(rrc, m_spcId);

                        for (int g = 0; g < G; g++) { // loop over all basises
                            Basis b = BasisS[g];
                            XDGBasis xb = XBasisS[g];
                            bool b_is_XDGbasis = (xb != null);

                            int n0, N;
                            if (b_is_XDGbasis) {
                                if (iSpc < 0) {
                                    n0 = int.MaxValue;
                                    N = int.MinValue;
                                } else {
                                    N = xb.DOFperSpeciesPerCell;
                                    n0 = N*iSpc;
                                }
                            } else {
                                n0 = 0;
                                N = b.GetLength(j);
                            }

                            if (N >= 0) {
                                int iFull = (int)(FullMap.GlobalUniqueCoordinateIndex(g, jCell, n0));
                                int iFrame = (int)(FrameMap.GlobalUniqueCoordinateIndex(g, jCell, 0));

                                for (int n = 0; n < N; n++) {
                                    ExtRangeTrafo[tt][iFrame + n - offset] = iFull + n;
                                }
                            }
                        }
                    }
                }
            }
             */

            {
                int NoOfSpc = this.m_Regions.SpeciesIdS.Count;
                int N_frame = this.FrameMap.MaxTotalNoOfCoordinatesPerCell;
                var FramBasisS = this.FrameMap.BasisS.ToArray();
                var FullBasisS = this.FullMap.BasisS.ToArray();


                int[] NBs = FramBasisS.Select(b => b.Length).ToArray();
                int[] NBf = FullBasisS.Select(b => b.MaximalLength).ToArray();

                Debug.Assert(NBf.Length == NBs.Length);
                for (int i = 0; i < NBf.Length; i++) {
                    Debug.Assert(NBf[i] % NBs[i] == 0);
                }
               
                bool[] XBasis = FullBasisS.Select(b => b is XDGBasis).ToArray();
                
                this.IndexTrafoWithinCell = new int[NoOfSpc, this.FrameMap.MaxTotalNoOfCoordinatesPerCell];

                for (int iSpc = 0; iSpc < NoOfSpc; iSpc++) { // loop over all possible species indices
                    int n_frame = 0;
                    int n_full = 0;

                    for (int ifld = 0; ifld < NBs.Length; ifld++) { // loop over all basises in the mapping

                        if (XBasis[ifld]) {
                            for (int n = 0; n < NBs[ifld]; n++) {
                                this.IndexTrafoWithinCell[iSpc, n_frame + n] = n_full + iSpc*NBs[ifld] + n;
                            }
                        } else {
                            for (int n = 0; n < NBs[ifld]; n++) {
                                this.IndexTrafoWithinCell[iSpc, n_frame + n] = n_full + n;
                            }
                        }

                        n_frame += NBs[ifld];
                        n_full += NBf[ifld];
                        
                    }
                }
            }
        }


        int[,] IndexTrafoWithinCell;


        /*
        int[] ExtRangeStart;
        int[] ExtRangeLen;
        int[][] ExtRangeTrafo;
        */

        /// <summary>
        /// conversion of framed index (index into this object) to full index (object which is framed),
        /// for local indices
        /// </summary>
        /// <param name="iFrame">
        /// a local index
        /// </param>
        public int Frame2Full_Loc(int iFrame) {
            if (iFrame < 0 || iFrame >= FrameMap.LocalLength)
                throw new IndexOutOfRangeException();
            
            int iGlob = Frame2Full_Lookup[iFrame];
            iGlob -= FullMap.i0;
            return iGlob;
        }

        /// <summary>
        /// conversion of framed index (index into this object) to full index (object which is framed). 
        /// </summary>
        /// <param name="iFrame">
        /// a global index
        /// </param>
        public int Frame2Full(int iFrame) {
            if (FrameMap.IsInLocalRange(iFrame)) {
                int iLoc = FrameMap.TransformIndexToLocal(iFrame);
                return Frame2Full_Lookup[iLoc];
            } else {
                if (!m_supportExternal) {
                    throw new IndexOutOfRangeException();
                } else {
                    /*
                    int TT = ExtRangeStart.Length;
                    Debug.Assert(ExtRangeLen.Length == TT);
                    Debug.Assert(ExtRangeTrafo.Length == TT);

                    for (int tt = 0; tt < TT; tt++) { // since the number of 'peer processes' is usually low, this loop should not harm to much
                        int iX = iFrame - ExtRangeStart[tt];
                        if (iX >= 0 && iX < ExtRangeLen[tt]) {
                            return ExtRangeTrafo[tt][iX];
                        }
                    }
                    throw new IndexOutOfRangeException();
                     */
                    var Parallel = this.FrameMap.GridDat.iParallel;

                    int N_frame = FrameMap.MaxTotalNoOfCoordinatesPerCell;
                    int N_full = FullMap.MaxTotalNoOfCoordinatesPerCell;


                    int jCellGlob = iFrame/N_frame;
                    int jCellLoc, iSpc;
                    if (m_Last_jCellGlob != jCellGlob) {
                        jCellLoc = Parallel.Global2LocalIdx[jCellGlob];
                        iSpc = m_Regions.GetSpeciesIndex(this.m_spcId, jCellLoc);
                        m_Last_jCellGlob = jCellGlob;
                        m_Last_jCellLoc = jCellLoc;
                        m_Last_iSpc = iSpc;
                    } else {
                        jCellLoc = m_Last_jCellLoc;
                        iSpc = m_Last_iSpc;
                    }

                    int n_frame = iFrame - N_frame*(jCellGlob);
                    int n_full = this.IndexTrafoWithinCell[iSpc, n_frame];

                    return jCellGlob*N_full + n_full;
                }
            }
        }

        int m_Last_jCellGlob = int.MinValue;
        int m_Last_jCellLoc = int.MinValue;
        int m_Last_iSpc = int.MinValue;



        /*

        int[] Full2Frame_Lookup;

        void Setup_Full2Frame() {
            
            int L = FullMap.LocalLength;
            Full2Frame_Lookup = new int[L];
            Basis[] BasisS = FullMap.BasisS.ToArray();
            XDGBasis[] XBasisS = new XDGBasis[BasisS.Length];
            for (int i = 0; i < BasisS.Length; i++)
                XBasisS[i] = BasisS[i] as XDGBasis;


            for (int iLoc = 0; iLoc < L; iLoc++) {
                int j, g, n;
                FullMap.LocalFieldCoordinateIndex(iLoc, out g, out j, out n);
                ReducedRegionCode rrc;
                int NoOfSpc = m_lsTrk.GetNoOfSpecies(j, out rrc);
                int iSpc = m_lsTrk.GetSpeciesIndexFromId(rrc, m_spcId);
                Basis b = BasisS[g];
                XDGBasis xb = XBasisS[g];
                bool b_is_XDGbasis = (xb != null);

                if (b_is_XDGbasis) {
                    Full2Frame_Lookup[iLoc] = (int)FrameMap.GlobalUniqueCoordinateIndex(g, j, n);
                } else {
                    if(iSpc < 0) {
                        Full2Frame_Lookup[iLoc] = int.MinValue;
                    } else {
                        int Nsep = xb.DOFperSpeciesPerCell;
                        Full2Frame_Lookup[iLoc] = (int)FrameMap.GlobalUniqueCoordinateIndex(g, j, n - iSpc * Nsep);
                    }
                }
            }
        }


        /// <summary>
        /// conversion of full index (object which is framed) to framed index (index into this object). 
        /// </summary>
        public int Full2Frame(int iFull) {
            FullMap.TestIfInLocalRange(iFull);
            int iLoc = FullMap.TransformToLocal(iFull);
            if (Full2Frame_Lookup == null)
                Setup_Full2Frame();
            return Full2Frame_Lookup[iLoc];
        }
         */
    }


    /// <summary>
    /// XDG-related extension methods for <see cref="UnsetteledCoordinateMapping"/>
    /// </summary>
    static public class UnsetteledCoordinateMapping_Extensions {


        /// <summary>
        /// returns global unique indices which correlate to a certain sub-set of ...
        /// <list type="bullet">
        /// <item>the mappings's basis (<paramref name="mapping"/>, <see cref="UnsetteledCoordinateMapping.BasisS"/>).</item>
        /// <item>species (<paramref name="_SpcIds"/>).</item>
        /// <item>cells (<paramref name="cm"/>).</item>
        /// </list>
        /// </summary>
        /// <param name="mapping">
        /// the mapping
        /// </param>
        /// <param name="Fields">
        /// determines which fields should be in the sub-vector-indices-list
        /// </param>
        /// <param name="_SpcIds">
        /// determines which species should be in the sub-vector-indices-list; null indicates all species.
        /// </param>
        /// <param name="cm">
        /// determines which cells should be in the sub-vector-indices-list; null indicates all cells.
        /// </param>
        /// <param name="regions">
        /// Determines which XDG-coordinate belongs to which species.
        /// </param>
        /// <param name="presenceTest">
        /// if true, cells where some species has zero measure are excluded from the sub-vector-indices-list.
        /// </param>
        /// <param name="drk">
        /// a cell-agglomeration object: if null, ignored; if provided,
        /// cells which are eliminated by the agglomeration are excluded from the sub-vector-indices-list
        /// </param>
        /// <returns>a list of global (over all MPI processes) unique indices.</returns>
        static public int[] GetSubvectorIndices(this UnsetteledCoordinateMapping mapping, LevelSetTracker.LevelSetRegions regions, int[] Fields,
            ICollection<SpeciesId> _SpcIds = null, CellMask cm = null, bool presenceTest = true, MultiphaseCellAgglomerator drk = null) {

            #region Check Arguments
            // =========
            SpeciesId[] SpcIds = null;
            if (_SpcIds == null) {
                // take all species, if arg is null
                SpcIds = regions.SpeciesIdS.ToArray();
            }
            else {
                SpcIds = _SpcIds.ToArray();
            }

            if (SpcIds.Length <= 0)
                // return will be empty for no species
                return new int[0]; 
            #endregion

            #region collect cells which are excluded by agglomeration
            // =================================================

            int[][] ExcludeLists;
            if (drk != null) {
                ExcludeLists = new int[regions.SpeciesIdS.Count][]; //  new Dictionary<SpeciesId, int[]>();
                foreach (var id in SpcIds) {
                    int[] ExcludeList = drk.GetAgglomerator(id).AggInfo.AgglomerationPairs.Select(pair => pair.jCellSource).ToArray();
                    Array.Sort(ExcludeList);
                    ExcludeLists[id.cntnt - LevelSetTracker.___SpeciesIDOffest] = ExcludeList;
                }
            }
            else {
                ExcludeLists = null;
            } 
            #endregion
            
            #region determine over which cells we want to loop
            // ==========================================
            CellMask loopCells;
            if ((new HashSet<SpeciesId>(regions.SpeciesIdS)).SetEquals(new HashSet<SpeciesId>(SpcIds))) {
                if (cm == null) {
                    loopCells = CellMask.GetFullMask(regions.GridDat);
                }
                else {
                    loopCells = cm;
                }
            }
            else {
                loopCells = regions.GetSpeciesMask(SpcIds[0]);
                for (int i = 1; i < SpcIds.Length; i++)
                    loopCells = loopCells.Intersect(regions.GetSpeciesMask(SpcIds[i]));
                if (cm != null)
                    loopCells = loopCells.Intersect(cm);
            } 
            #endregion

            #region build list
            // ==========

            var R = new List<int>();

            // which basis is XDG?
            var _BasisS = mapping.BasisS.ToArray();
            XDGBasis[] _XBasisS = new XDGBasis[_BasisS.Length];
            for (int i = 0; i < _BasisS.Length; i++)
                _XBasisS[i] = _BasisS[i] as XDGBasis;


            Tuple<int, SpeciesId>[] iSpcS = new Tuple<int, SpeciesId>[SpcIds.Length];
            int KK;
            int Len_iSpcS = iSpcS.Length;


            int[] ExcludeListsPointer = new int[regions.SpeciesIdS.Count];


            // loop over cells?
            foreach (int jCell in loopCells.ItemEnum) {

                int NoOfSpc = regions.GetNoOfSpecies(jCell);

                #region  
                //
                // in cell 'jCell', for each species in 'SpcIds' ...
                // * is the species present in 'jCell'?
                // * what is the species index in 'jCell'?
                // * we also have to consider that due to agglomeration, the respective cut-cell for the species might have been agglomerated.
                // so we determine
                // * 'KK' is the number of species (from 'SpcIds') which is actually present and not agglomerated in 'jCell'
                // * for i in (0 ... KK-1), iSpcS[i] is a tuple of (Species-index-in-cell, SpeciesID)
                //
                // yes, the following code sucks:
                KK = 0;
                for (int k = 0; k < Len_iSpcS; k++) { // loop over all species (provided as argument)...
                    SpeciesId id = SpcIds[k];
                    int __iSpcS = regions.GetSpeciesIndex(id, jCell);
                    if (__iSpcS >= 0) {
                        // species index is non-negative, so indices for the species are allocated ...
                        if (!presenceTest || regions.IsSpeciesPresentInCell(id, jCell)) {
                            // species is also present (i.e. has a greater-than-zero measure ...

                            if(drk == null) {
                                // no agglomeration exclution
                                iSpcS[KK] = new Tuple<int, SpeciesId>(__iSpcS, id);
                                KK++;
                            } else {
                                int q = id.cntnt - LevelSetTracker.___SpeciesIDOffest;
                                var el = ExcludeLists[q];

                                while(ExcludeListsPointer[q] < el.Length
                                    && el[ExcludeListsPointer[q]] < jCell) {
                                    ExcludeListsPointer[q]++;
                                }
                                Debug.Assert(ExcludeListsPointer[q] >= el.Length || el[ExcludeListsPointer[q]] >= jCell);

                                if(ExcludeListsPointer[q] >= el.Length || el[ExcludeListsPointer[q]] != jCell) {
                                    // cell is also not agglomerated ...
                                    iSpcS[KK] = new Tuple<int, SpeciesId>(__iSpcS, id);
                                    KK++;
                                }
                            }
                        }
                    }
                }
                if (KK > 1)
                    Array.Sort<Tuple<int, SpeciesId>>(iSpcS, 0, KK, new ilPSP.FuncComparer<Tuple<int, SpeciesId>>((a, b) => a.Item1 - b.Item1));
                // --------- esc (end of sucking code) 
                #endregion

                for (int k = 0; k < KK; k++) { // loop over species in cell
                    int iSpc = iSpcS[k].Item1; 
                    
                    for (int i = 0; i < Fields.Length; i++) {
                        int iField = Fields[i];

                        if (_XBasisS[iField] == null) {
                            int N = _BasisS[iField].GetLength(jCell);
                            int i0 = mapping.LocalUniqueCoordinateIndex(iField, jCell, 0);
                            for (int n = 0; n < N; n++) {
                                R.Add(i0 + n);
                            }
                        }
                        else {
                            int N = _XBasisS[iField].DOFperSpeciesPerCell;
                            int i0 = mapping.LocalUniqueCoordinateIndex(iField, jCell, 0) + iSpc * N;
                            for (int n = 0; n < N; n++) {
                                R.Add(i0 + n);
                            }
                        }
                    }
                }

            } 
            #endregion

            return R.ToArray();
        }

        /*

        /// <summary>
        /// 
        /// </summary>
        /// <param name="map"></param>
        /// <param name="vec"></param>
        static public void ClearIrrelevantNearfieldDOFs<T>(this UnsetteledCoordinateMapping map, T vec)
            where T : IList<double> 
        {
            //
            if(vec.Count != map.LocalLength)
                throw new ArgumentException("Mismatch between vector length and NUpdate in mapping.");

            foreach(var B in map.BasisS) {
                if(!(B is XDGBasis))
                    throw new NotSupportedException("All basis objects in the map are expected to be an XDG basis.");
            }

            int[] NSpc = map.BasisS.Select(b => ((XDGBasis)b).NonX_Basis.Length).ToArray();
            int NF = NSpc.Length;
            //int Ntot_NonX = NSpc.Sum();
            //int Ntot_X = map.BasisS.Sum(b => b.MaximalLength);
            int[] NXdg = map.BasisS.Select(b => b.MaximalLength).ToArray();

            var LsTrk = ((XDGBasis)map.BasisS[0]).Tracker;
            var gDat = LsTrk.GridDat;
            int J = gDat.Cells.NoOfLocalUpdatedCells;

            int NoLevSet = LsTrk.LevelSets.Count;
            double[] _levSetSign = new double[4];
            double[] levSetSign = new double[4];

            int NoPresentSpc = 0;
            int[] PresentSpc = new int[LsTrk.SpeciesIdS.Count];
            
            for(int j = 0; j < J; j++) {
                ReducedRegionCode rrc;
                int NoOfSpecies;
                NoOfSpecies = LsTrk.GetNoOfSpecies(j, out rrc);

                if(NoOfSpecies > 1) {
                    
                    for(int iLs = 0; iLs < NoLevSet; iLs++) {
                        _levSetSign[iLs] = LevelSetTracker.DecodeLevelSetDist(LsTrk.m_LevSetRegions[j], iLs);
                    }
                    Array.Copy(_levSetSign, levSetSign, 4);
                    for(int iLs = 0; iLs < NoLevSet; iLs++) {
                        if(_levSetSign[iLs] != 0.0) {
                            
                        }

                    }


                    for(int iSpc = 0; iSpc < NoOfSpecies; iSpc++) {
                        

                    }
                    
                } else {
                    for(int ifld = 0; ifld < NF; ifld++) {
                        for(int n = NSpc[ifld]; n < NXdg[ifld]; n++) {
                            int idx = map.LocalUniqueCoordinateIndex(ifld, j, n);
                            Debug.Assert(vec[idx] == 0); // wenn es nur eine species in Zelle 'j' gibt,
                            //                              sollte da gar kein Eintrag ungleich 0.0 existieren dürfen.
                        }
                    }

                }
            }
         

        }*/
    }
}
