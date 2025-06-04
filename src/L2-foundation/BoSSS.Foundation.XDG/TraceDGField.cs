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
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP.Utils;
using ilPSP;
using MPI.Wrappers;
using BoSSS.Foundation.Voronoi;
using static System.Net.Mime.MediaTypeNames;

namespace BoSSS.Foundation.XDG {


    /// <summary>
    /// class for trace-fields, i.e., fields on the level-set surface  (<see cref="TraceDGBasis"/>);
    /// </summary>
    public partial class TraceDGField : SinglePhaseField, IObserver<LevelSetTracker.LevelSetRegions> {

        /// <summary>
        /// a factory that creates <see cref="Trace"/>-DG-fields.
        /// </summary>
        /// <param name="__Basis">
        /// The basis that is used for this field;
        /// Must be a <see cref="TraceDGBasis"/>-object.
        /// </param>
        /// <param name="__Identification">
        /// identification string for this field;
        /// This can be null or empty, 
        /// however, if IO should be performed for this object, the identification must be unique 
        /// within a given context
        /// </param>
        /// <returns>a <see cref="TraceDGField"/>-instance</returns>
        public static new TraceDGField Factory(Basis __Basis, String __Identification) {
            return new TraceDGField((TraceDGBasis)__Basis, __Identification);
        }

        /// <summary> constructor </summary>
        /// <param name="basis"></param>
        public TraceDGField(TraceDGBasis basis)
            : this(basis, null) {
        }

        /// <summary> constructor </summary>
        /// <param name="basis"></param>
        /// <param name="Identification">
        /// identification string for this field;
        /// 
        /// This can be null or empty, 
        /// however, if IO should be performed for this object, the identification must be unique 
        /// within a given context
        /// </param>
        public TraceDGField(TraceDGBasis basis, string Identification)
            : base(basis, Identification) {
            m_TraceBasis = basis;

            // allocate memory // For the moment, just use a MultidimensionalArray.
            // ---------------

            //int J = this.GridDat.iLogicalCells.Count;
            //m_Coordinates = new FieldStorage(J, m_TraceBasis.MinimalLength, m_TraceBasis.MaximalLength);
            //m_Coordinates.BeginResize(m_TraceBasis.MaximalLength);
            //for (int j = 0; j < J; j++) {
            //    m_Coordinates.Resize(j, m_TraceBasis.GetLength(j));
            //}
            //m_Coordinates.FinishResize();
            //m_TrackerVersionCnt = m_CCBasis.Tracker.VersionCnt;


            // register field with level set tracker
            // -------------------------------------
            //this.m_TrackerVersionCnt = m_CCBasis.Tracker.Regions.Version;
            m_TraceBasis.Tracker.Subscribe(this);
            this.m_TrackerVersion = m_TraceBasis.Tracker.VersionCnt;
            this.OnNext(m_TraceBasis.Tracker.Regions); // initialize data structures.
            Debug.Assert(this.m_TrackerVersion == m_TraceBasis.Tracker.VersionCnt);
        }

        TraceDGBasis m_TraceBasis;

        /// <summary>
        /// Identical to <see cref="DGField.Basis"/> but return an instance of
        /// <see cref="XDGBasis"/> instead of <see cref="Basis"/>
        /// </summary>
        /// <remarks>
        /// A more elegant (i.e. polymorphic) way to accomplish this would be
        /// to make the super class generic but since the functionality is
        /// identical and so polymorphism is not an issue
        /// </remarks>
        public new TraceDGBasis Basis {
            get {
                Debug.Assert(Object.ReferenceEquals(m_TraceBasis, base.Basis));
                return m_TraceBasis;
            }
        }

        


        /// <summary>
        /// see <see cref="DGField.Acc(double,DGField,CellMask)"/>;
        /// </summary>
        /// <remarks>
        /// the other field <paramref name="a"/> can be either a <see cref="XDGField"/>
        /// or a <see cref="XDGField"/>.
        /// </remarks>
        public override void AccLaidBack(double mult, DGField a, CellMask cm = null) {
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            
            IEnumerator<Chunk> cmEnu = (cm == null) ? null : cm.GetEnumerator();

            {
                Basis a_Basis = a.Basis;
                IMatrix _aCoordinates = a.Coordinates;

               
                while (true) {

                    Chunk cnk;
                    if (cmEnu == null) {
                        cnk.i0 = 0;
                        cnk.Len = J;
                    } else {
                        if (cmEnu.MoveNext())
                            cnk = cmEnu.Current;
                        else
                            break;
                    }

                    int JE = cnk.i0 + cnk.Len;
                    for (int j = cnk.i0; j < JE; j++) {
                        int Nj = m_TraceBasis.GetLength(j);
                        int Na = a_Basis.GetLength(j);
                        int No = Math.Min(Nj, Na);

                        for (int n = 0; n < No; n++) {
                            this.Coordinates[j, n] += mult * _aCoordinates[j, n];
                        }
                    }

                    if (cmEnu == null)
                        // full mask - done with one outer loop
                        break;
                }
            } 
        }

        /// <summary>
        /// see <see cref="DGField.Acc(double,DGField,CellMask)"/>;
        /// </summary>
        /// <remarks>
        /// the other field <paramref name="a"/> can be either a <see cref="XDGField"/>
        /// or a <see cref="XDGField"/>.
        /// </remarks>
        public override void Acc(double mult, DGField a, CellMask cm) {
            if (!a.Basis.Equals(this.Basis)) {
                throw new ArgumentException("DG basis of other field must be a equal to basis of this field.", "a");
            }

            AccLaidBack(mult, a, cm); // optimization possible
        }





        /// <summary>
        /// guess what?
        /// </summary>
        public override object Clone() {
            TraceDGField r = new TraceDGField(m_TraceBasis, this.Identification);
            r.Acc(1.0, this);
            r.m_UpdateBehaviour = this.m_UpdateBehaviour;
            return r;
        }

        /// <summary>
        /// guess what?
        /// </summary>
        new public TraceDGField CloneAs() {
            return (TraceDGField)Clone();
        }

        BehaveUnder_LevSetMoovement m_UpdateBehaviour = BehaveUnder_LevSetMoovement.PreserveMemory;

        /// <summary>
        /// defines the Behavior of the DG coordinates during a <see cref="LevelSetTracker.UpdateTracker"/>-call
        /// </summary>
        public BehaveUnder_LevSetMoovement UpdateBehaviour {
            get {
                return m_UpdateBehaviour;
            }
            set {
                m_UpdateBehaviour = value;
            }
        }


        // I am not sure if we still need it in future moving interface case
        //void AutoExtrapolate(SubGrid oldSpeciesSubGrid)
        //{
        //    LevelSetTracker LsTrk = m_TraceBasis.Tracker;
        //    CellMask allNearMask = m_TraceBasis.Tracker.Regions.GetNearFieldMask(m_TraceBasis.Tracker.NearRegionWidth);


        //    CellMask ExtrapolateTo = allNearMask.Intersect(LsTrk.Regions.GetSpeciesMask(Id));
        //    CellMask ExtrapolateFrom = oldSpeciesSubGrid.VolumeMask;

        //    this.CellExtrapolation(ExtrapolateTo, ExtrapolateFrom);
        //}


        /*

        int[] m_MPISendBufSize;

        int[] m_MPIRecvBufSize;

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.GetMPISendBufferSize"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <returns></returns>
        public override int GetMPISendBufferSize(int proc) {
            return m_MPISendBufSize[proc];
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.GetMPIRecvBufferSize"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <returns></returns>
        public override int GetMPIRecvBufferSize(int proc) {
            return m_MPIRecvBufSize[proc];
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.FillMPISendBuffer"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <param name="Buffer"></param>
        /// <param name="st"></param>
        /// <returns></returns>
        public override int FillMPISendBuffer(int proc, double[] Buffer, int st) {
            int[] CellIndexList = this.GridDat.iParallel.SendCommLists[proc];

            int I = CellIndexList.Length;
            int N = m_Coordinates.NMin;
            double[] bStor = m_Coordinates.m_BaseStorage;
            double[][] extStor = m_Coordinates.m_ExtendedStorage;

            int l = 0;
            for (int i = 0; i < I; i++) {
                int j = CellIndexList[i];
                int i0 = m_Coordinates.IndBase(j, 0);
                Array.Copy(bStor, i0, Buffer, st + l, N);
                l += N;

                double[] ext = extStor[j];
                if (ext != null) {
                    Array.Copy(ext, 0, Buffer, st + l, ext.Length);
                    l += ext.Length;
                }
            }

            // test code
            if (m_MPISendBufSize[proc] != l)
                throw new ApplicationException("internal error.");

            return l;
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.CopyFromMPIrecvBuffer"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <param name="Buffer"></param>
        /// <param name="st"></param>
        /// <returns></returns>
        public override int CopyFromMPIrecvBuffer(int proc, double[] Buffer, int st) {
            int N = m_Coordinates.NMin;
            double[] bStor = m_Coordinates.m_BaseStorage;
            double[][] extStor = m_Coordinates.m_ExtendedStorage;
            int j_insert = this.GridDat.iParallel.RcvCommListsInsertIndex[proc];
            int Len = this.GridDat.iParallel.RcvCommListsNoOfItems[proc];
            Len += j_insert;

            int l = 0;
            for (int j = j_insert; j < Len; j++) {

                int i0 = m_Coordinates.IndBase(j, 0);
                Array.Copy(Buffer, st + l, bStor, i0, N);
                l += N;

                double[] ext = extStor[j];
                if (ext != null) {
                    Array.Copy(Buffer, st + l, ext, 0, ext.Length);
                    l += ext.Length;
                }
            }

            // test code
            if (m_MPIRecvBufSize[proc] != l)
                throw new ApplicationException("internal error.");

            return l;
        }

        */

        /// <summary>
        /// initializes this field to be a copy of another field
        /// </summary>
        /// <param name="other">
        /// must be a <see cref="XDGField"/>-object
        /// </param>
        public override void CopyFrom(DGField other) {
            if (!(other is TraceDGField))
                throw new ApplicationException("unable to copy, because the other field is no TraceDG field");

            if (!other.Basis.Equals(this.Basis))
                throw new ApplicationException("unable to copy, because the DG polynomial basis of other field is different.");

            this.Clear();
            this.Acc(1.0, other);
        }

      
        #region IObserver<LevelSetInfo> Members

        /// <summary>
        /// Do nothing.
        /// </summary>
        public void OnCompleted() {
        }

        /// <summary>
        /// Do nothing.
        /// </summary>
        /// <param name="error"></param>
        public void OnError(Exception error) {
        }

        /// <summary>
        /// most recent entry for <see cref="LevelSetTracker.VersionCnt"/>
        /// </summary>
        public int m_TrackerVersion;

        /// <summary>
        /// Updated the data structure of this cut-cell DG field to reflect the
        /// latest status of the level set.
        /// </summary>
        /// <param name="levelSetStatus"></param>
        public void OnNext(LevelSetTracker.LevelSetRegions levelSetStatus) {
            int J = this.GridDat.iLogicalCells.Count;
            LevelSetTracker trk = m_TraceBasis.Tracker;

            int oldTrackerVersion = m_TrackerVersion;
            m_TrackerVersion = trk.VersionCnt;

            /*
            m_Coordinates.BeginResize(m_TraceBasis.MaximalLength);

            if (m_UpdateBehaviour == BehaveUnder_LevSetMoovement.PreserveMemory || m_UpdateBehaviour == BehaveUnder_LevSetMoovement.AutoExtrapolate) {
                if (trk.HistoryLength < 1)
                    throw new NotSupportedException("LevelSettracker must have a history length >= 1 in order to support 'PreserveMemory' or 'AutoExtrapolate'.");
            }

            // rearrange DG coordinates if regions have changed
            // ================================================
            if ((m_UpdateBehaviour == BehaveUnder_LevSetMoovement.PreserveMemory || m_UpdateBehaviour == BehaveUnder_LevSetMoovement.AutoExtrapolate)
                && (m_TrackerVersion - oldTrackerVersion) > 0 && levelSetStatus.m_LevSetRegions_b4Update != null) {

                if((m_TrackerVersion - oldTrackerVersion) > 1) {
                    string message = $"The update behavior '{BehaveUnder_LevSetMoovement.PreserveMemory}' and '{BehaveUnder_LevSetMoovement.AutoExtrapolate}' do not work if every tracker update is paired with a 'PushStacks()' call (DGfield '" + (this.Identification ?? "no name set") + "').";
                    //Console.WriteLine(message);
                    throw new NotSupportedException(message);
                }

                // rearrange DG coordinates, preserve State of each species 
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

                double[] coordsFull = new double[m_TraceBasis.MaximalLength];
                ushort[] OldRegionCode = levelSetStatus.m_LevSetRegions_b4Update;
                ushort[] NewRegionCode = levelSetStatus.RegionsCode;
                Debug.Assert(!object.ReferenceEquals(OldRegionCode, NewRegionCode));

                int Nsep = m_TraceBasis.DOFperSpeciesPerCell;
                int i0CmnFull = trk.TotalNoOfSpecies * Nsep;
                //int[] i0SepFull = new int[trk.TotalNoOfSpecies];
                //for( int iSepc = 0; iSepc < i0SepFull.Length; iSepc++) i0SepFull[iSepc] = iSepc*Nsep;

                for (int j = 0; j < J; j++) {

                    ushort oldCd = OldRegionCode[j];
                    ushort newCd = NewRegionCode[j];
                                       
                    if (oldCd != newCd  // quick pre-test
                         && ReducedRegionCode.Extract(oldCd) != ReducedRegionCode.Extract(newCd)) {
                        // something changed

                        Array.Clear(coordsFull, 0, coordsFull.Length);

                        // save coordinates
                        // ----------------
                        {
                            ReducedRegionCode OldInd;
                            int OldNo = trk.GetNoOfSpeciesByRegionCode(oldCd, out OldInd);

                            // separate coordinates
                            for (int iSpec = 0; iSpec < OldNo; iSpec++) {
                                SpeciesId SpecId = trk.GetSpeciesIdFromIndex(OldInd, iSpec);
                                int iSpecGlob = SpecId.cntnt - LevelSetTracker.___SpeciesIDOffest;
                                int i0SepOld = iSpec * Nsep;
                                int i0SepFull = iSpecGlob * Nsep;

                                for (int _n = 0; _n < Nsep; _n++)
                                    coordsFull[_n + i0SepFull] = m_Coordinates[j, _n + i0SepOld];
                            }
                        }

                        // resize array
                        // ------------
                        m_Coordinates.Resize(j, m_TraceBasis.GetLength(j));

                        // write back coordinates in new order
                        // -----------------------------------
                        {
                            ReducedRegionCode NewInd;
                            int NewNo = trk.GetNoOfSpeciesByRegionCode(newCd, out NewInd);

                            // separate coordinates
                            for (int iSpec = 0; iSpec < NewNo; iSpec++) {
                                SpeciesId SpecId = trk.GetSpeciesIdFromIndex(NewInd, iSpec);
                                int iSpecGlob = SpecId.cntnt - LevelSetTracker.___SpeciesIDOffest;
                                int i0SepNew = iSpec * Nsep;
                                int i0SepFull = iSpecGlob * Nsep;

                                for (int _n = 0; _n < Nsep; _n++)
                                    m_Coordinates[j, _n + i0SepNew] = coordsFull[_n + i0SepFull];
                            }
                        }
                    } else {
                        m_Coordinates.Resize(j, m_TraceBasis.GetLength(j)); // need to call resize in every case
                    }
                }
            } else {
                // just allocate/free memory
                // +++++++++++++++++++++++++

                for (int j = 0; j < J; j++) {
                    int l = m_TraceBasis.GetLength(j);
                    m_Coordinates.Resize(j, l);
                }

            }

            m_Coordinates.FinishResize();
            */


            // update MPI buffer size
            // ======================
            //int size = this.GridDat.CellPartitioning.MpiSize;
            //if (size > 1) {
            //    if (m_MPIRecvBufSize == null)
            //        m_MPIRecvBufSize = new int[size];
            //    if (m_MPISendBufSize == null)
            //        m_MPISendBufSize = new int[size];

            //    for (int p = 0; p < size; p++) {
            //        // send list
            //        {
            //            int[] senditems = this.GridDat.iParallel.SendCommLists[p];
            //            if (senditems != null) {
            //                int L = senditems.Length;

            //                int sz = 0;
            //                for (int l = 0; l < L; l++)
            //                    sz += m_TraceBasis.GetLength(senditems[l]);
            //                m_MPISendBufSize[p] = sz;
            //            } else {
            //                m_MPISendBufSize[p] = int.MinValue;
            //            }
            //        }

            //        // receive list
            //        {
            //            int L = this.GridDat.iParallel.RcvCommListsNoOfItems[p];
            //            if (L > 0) {
            //                int j0 = this.GridDat.iParallel.RcvCommListsInsertIndex[p];
            //                L += j0;
            //                int sz = 0;
            //                for (int j = j0; j < L; j++) {
            //                    sz += m_TraceBasis.GetLength(j);
            //                }
            //                m_MPIRecvBufSize[p] = sz;
            //            } else {
            //                m_MPIRecvBufSize[p] = int.MinValue;
            //            }
            //        }
            //    }
            //}

            // do Extrapolation, if necessary
            //if (m_UpdateBehaviour == BehaveUnder_LevSetMoovement.AutoExtrapolate && trk.PopulatedHistoryLength >= 1) {
            //    AutoExtrapolate(trk.RegionsHistory[0].GetNearFieldSubgrid(trk.NearRegionWidth));
            //}
        }

        #endregion
    }
}
