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
    /// see <see cref="XDGField.UpdateBehaviour"/>;
    /// </summary>
    public enum BehaveUnder_LevSetMoovement {

        /// <summary>
        /// just reallocates the memory structure of the cut-cell field;
        /// The behavior in the near band is undefined;
        /// </summary>
        JustReallocate,

        /// <summary>
        /// the values of the DG coordinates are preserved; For new cells (i.e. 
        /// cells where a new species just entered during the last level set movement),
        /// the DG coordinates of the new species are unassigned.
        /// </summary>
        PreserveMemory,

        /// <summary>
        /// most expensive option: DG coordinates of new cells (i.e. 
        /// cells where a new species just entered during the last level set movement),
        /// are extrapolated from cells where the DG coordinates for the corresponding species are known.
        /// </summary>
        AutoExtrapolate
    }


    /// <summary>
    /// a DG field for a cut-cell -basis (<see cref="XDGBasis"/>);
    /// </summary>
    public partial class XDGField : DGField, IObserver<LevelSetTracker.LevelSetRegions> {

        /// <summary>
        /// a factory that creates <see cref="XDGField"/>-DG-fields.
        /// </summary>
        /// <param name="__Basis">
        /// The basis that is used for this field;
        /// Must be a <see cref="XDGBasis"/>-object.
        /// </param>
        /// <param name="__Identification">
        /// identification string for this field;
        /// This can be null or empty, 
        /// however, if IO should be performed for this object, the identification must be unique 
        /// within a given context
        /// </param>
        /// <returns>a <see cref="SinglePhaseField"/>-instance</returns>
        public static XDGField Factory(Basis __Basis, String __Identification) {
            return new XDGField((XDGBasis)__Basis, __Identification);
        }

        /// <summary> constructor </summary>
        /// <param name="basis"></param>
        public XDGField(XDGBasis basis)
            : this(basis, null) {
        }

        /// <summary> constructor </summary>
        /// <param name="basis"></param>
        /// <param name="Identification">
        /// identification string for this field;<br/>
        /// This can be null or empty, 
        /// however, if IO should be performed for this object, the identification must be unique 
        /// within a given context
        /// </param>
        public XDGField(XDGBasis basis, string Identification)
            : base(basis, Identification) {
            m_CCBasis = basis;

            // allocate memory
            // ---------------

            int J = this.GridDat.iLogicalCells.Count;
            m_Coordinates = new FieldStorage(J, m_CCBasis.MinimalLength, m_CCBasis.MaximalLength);
            m_Coordinates.BeginResize(m_CCBasis.MaximalLength);
            for (int j = 0; j < J; j++) {
                m_Coordinates.Resize(j, m_CCBasis.GetLength(j));
            }
            m_Coordinates.FinishResize();
            //m_TrackerVersionCnt = m_CCBasis.Tracker.VersionCnt;


            // register field with level set tracker
            // -------------------------------------
            //this.m_TrackerVersionCnt = m_CCBasis.Tracker.Regions.Version;
            m_CCBasis.Tracker.Subscribe(this);
            this.m_TrackerVersion = m_CCBasis.Tracker.VersionCnt;
            this.OnNext(m_CCBasis.Tracker.Regions); // initialize data structures.
            Debug.Assert(this.m_TrackerVersion == m_CCBasis.Tracker.VersionCnt);
        }

        XDGBasis m_CCBasis;

        /// <summary>
        /// Identical to <see cref="DGField.Basis"/> but return an instance of
        /// <see cref="XDGBasis"/> instead of <see cref="Basis"/>
        /// </summary>
        /// <remarks>
        /// A more elegant (i.e. polymorphic) way to accomplish this would be
        /// to make the super class generic but since the functionality is
        /// identical and so polymorphism is not an issue
        /// </remarks>
        public new XDGBasis Basis {
            get {
                return (XDGBasis)base.Basis;
            }
        }

        //int m_TrackerVersionCnt = int.MinValue;

        ///// <summary>
        ///// dirty hack used by dynamic load balancing.
        ///// </summary>
        //internal void Override_TrackerVersionCnt(int i) {
        //    m_TrackerVersionCnt = i;
        //}

        /// <summary>
        /// caches the return value of <see cref="GetSpeciesShadowField(SpeciesId)"/>;
        /// </summary>
        SortedDictionary<int, SpeciesShadowField> m_SpeciesShadowFields = new SortedDictionary<int, SpeciesShadowField>();

        /// <summary>
        /// creates the shadow field for the given Species
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        /// <remarks>
        /// return value is cached
        /// </remarks>
        public SpeciesShadowField GetSpeciesShadowField(SpeciesId id) {
            MPICollectiveWatchDog.Watch();
            if (!m_SpeciesShadowFields.ContainsKey(id.cntnt)) {
                SpeciesShadowField ssf = new SpeciesShadowField(this, id);
                m_SpeciesShadowFields.Add(id.cntnt, ssf);
            }
            return m_SpeciesShadowFields[id.cntnt];
        }

        /// <summary>
        /// see <see cref="GetSpeciesShadowField(SpeciesId)"/>;
        /// </summary>
        public SpeciesShadowField GetSpeciesShadowField(string SpecName) {
            return GetSpeciesShadowField(m_CCBasis.Tracker.GetSpeciesId(SpecName));
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
            int Nsep = m_CCBasis.DOFperSpeciesPerCell;

            IEnumerator<Chunk> cmEnu = (cm == null) ? null : cm.GetEnumerator();

            if (a is ConventionalDGField) {
                // version which works with SinglePhaseField's
                // + + + + + + + + + + + + + + + + + + + + + + 

                IMatrix _aCoordinates = a.Coordinates;

                int _NSepMin = Math.Min(a.Basis.MaximalLength, Nsep);
                //int Nmin = Math.Min(Nsep + Ncom, a.Basis.MaximalLength);


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
                        ReducedRegionCode rrc;
                        int NoOfSpec = m_CCBasis.Tracker.Regions.GetNoOfSpecies(j, out rrc);

                        int Nmin = Math.Min(Nsep, a.Basis.GetLength(j));

                        int n0;
                        for (int iSpec = 0; iSpec < NoOfSpec; iSpec++) {

                            n0 = iSpec * Nsep;
                            for (int n = 0; n < _NSepMin; n++) {
                                m_Coordinates[j, n + n0] += mult * _aCoordinates[j, n];
                            }
                        }
                        n0 = NoOfSpec * Nsep;
                        for (int n = Nsep; n < Nmin; n++) {
                            m_Coordinates[j, n + n0 - Nsep] += mult * _aCoordinates[j, n];
                        }
                    }

                    if (cmEnu == null)
                        // full mask - done with one outer loop
                        break;
                }
            } else if (a is XDGField) {
                // version which works with XdgField's
                // + + + + + + + + + + + + + + + + + + + +

                XDGField _a = (XDGField)a;
                FieldStorage _aCoordinates = _a.m_Coordinates;

                if (!object.ReferenceEquals(_a.Basis.Tracker, this.Basis.Tracker))
                    throw new ArgumentException("XDG basis of other object may not be assigned to a different tracker.");


                int _Nsep = _a.m_CCBasis.DOFperSpeciesPerCell;

                int _NSepMin = Math.Min(_Nsep, Nsep);

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
                        ReducedRegionCode rrc;
                        int NoOfSpec = m_CCBasis.Tracker.Regions.GetNoOfSpecies(j, out rrc); // both fields own the same tracker,
                        // so the number of species in cell j is equal for
                        // both fields, '_a' and 'this'.
                        int n0, _n0;
                        for (int iSpec = 0; iSpec < NoOfSpec; iSpec++) {

                            n0 = iSpec * Nsep;
                            _n0 = iSpec * _Nsep;
                            for (int n = 0; n < _NSepMin; n++) {
                                m_Coordinates[j, n + n0] += mult * _aCoordinates[j, n + _n0];
                            }
                        }
                        n0 = NoOfSpec * Nsep;
                    }

                    if (cmEnu == null)
                        break;
                }

            } else {
                throw new NotSupportedException("XdgField.Acc(...) does not support fields of type " + a.GetType().FullName);
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
        /// Evaluates
        /// \f$ 
        /// f_A(x_i) - f_B(x_i)
        /// \f$ 
        /// in all nodes with spatial coordinates
        /// \f$ x_i\f$  defined in
        /// <paramref name="NodeSet"/> in all cells in the range
        /// [<paramref name="j0"/>; <paramref name="j0"/> + <paramref name="Len"/>]
        /// where \f$ f_A(x)\f$  and
        /// \f$ f_B(x)\f$  are the values of the
        /// polynomial defined by the coefficients stored for the species
        /// identified by <paramref name="SpeciesIdA"/> and
        /// <paramref name="SpeciesIdB"/>, respectively. The result will then
        /// be stored in <paramref name="result"/>
        /// </summary>
        /// <param name="j0">
        /// The first cell in the evaluation range.
        /// </param>
        /// <param name="Len">
        /// The length of the evaluation range.
        /// </param>
        /// <param name="NodeSet">
        /// The node set.
        /// </param>
        /// <param name="result">
        /// On exit, the result of the calculation will have been accumulated
        /// here according to the formula
        /// <paramref name="result"/> := <paramref name="ResultPreScale"/>*<paramref name="result"/> + jump
        /// <list type="bullet">
        ///     <item>1st index: Cell index minus <paramref name="j0"/></item>
        ///     <item>2nd index: Node index</item>
        /// </list>
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// An offset for the first index of <paramref name="result"/>
        /// </param>
        /// <param name="ResultPreScale">
        /// Scaling of the original content of <paramref name="result"/>. Use
        /// 0.0 to overwrite existing value
        /// </param>
        /// <param name="SpeciesIdA">
        /// Id of the first species to be evaluated
        /// </param>
        /// <param name="SpeciesIdB">
        /// Id of the second species to be evaluated
        /// </param>
        public void EvaluateJumpHeight(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale, SpeciesId SpeciesIdA, SpeciesId SpeciesIdB) {
            //MultidimensionalArray basisValues = m_Basis.Evaluate(NodeSet);
            //int NoOfNodesPerCell = basisValues.GetLength(0);
            int NoOfNodes = NodeSet.NoOfNodes;
            if (result.GetLength(1) != NoOfNodes)
                throw new ArgumentException();

            var _result = result.ExtractSubArrayShallow(new int[] { ResultCellindexOffset, 0, }, new int[] { ResultCellindexOffset + Len - 1, NoOfNodes - 1 });
            if (ResultPreScale == 0.0)
                _result.Clear();
            else if (ResultPreScale != 1.0)
                _result.Scale(ResultPreScale);

            this.GetSpeciesShadowField(SpeciesIdB).Evaluate(j0, Len, NodeSet, _result, 1.0);
            this.GetSpeciesShadowField(SpeciesIdA).Evaluate(j0, Len, NodeSet, _result, -1.0);
            _result.Scale(-1.0);
        }

        

        /// <summary>
        /// Evaluates the cut-cell DG - field;
        /// </summary>
        /// <param name="_j0">local index of the first cell to evaluate</param>
        /// <param name="_Len">Number of cells to evaluate</param>
        /// <param name="_NodeSet">
        /// as usual, the node set;
        /// </param>
        /// <param name="_result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="_ResultPreScale"/>;<br/>
        /// 1st index: cell index minus <paramref name="_j0"/>;<br/>
        /// 2nd index: node index;
        /// </param>
        /// <param name="_ResultCellindexOffset">
        /// an offset for the first index of <paramref name="_result"/>;
        /// </param>
        /// <param name="_ResultPreScale">
        /// see <paramref name="_result"/>
        /// </param>
        public override void Evaluate(int _j0, int _Len, NodeSet _NodeSet, MultidimensionalArray _result, int _ResultCellindexOffset, double _ResultPreScale) {
            int _M = _NodeSet.NoOfNodes; // number of nodes per cell

            if (_result.Dimension != 2)
                throw new ArgumentOutOfRangeException("result", "dimension of result array must be 2");
            if (_result.GetLength(1) != _M)
                throw new ArgumentOutOfRangeException();

            MultidimensionalArray[] m_Evaluate_SpeciesEvalBuffer = new MultidimensionalArray[0];

            GenericEval(_j0, _Len, _NodeSet, _result, _ResultCellindexOffset, _ResultPreScale,
                DGField.EvaluateInternal,
                ref m_Evaluate_SpeciesEvalBuffer,
                __M => new int[] { 1, __M },
                delegate (MultidimensionalArray R, int offset, int m, int SpcInd, MultidimensionalArray[] SR) {
                    double r = R[offset, m] * _ResultPreScale;
                    r += SR[SpcInd][0, m];
                    R[offset, m] = r;
                });
        }

        /// <summary>
        /// <see cref="DGField.EvaluateEdge(int, int, NodeSet, MultidimensionalArray, MultidimensionalArray, MultidimensionalArray, MultidimensionalArray, MultidimensionalArray, MultidimensionalArray, int, double)"/>
        /// </summary>
        public override void EvaluateEdge(int e0, int Len, NodeSet NS, MultidimensionalArray ValueIN, MultidimensionalArray ValueOT, MultidimensionalArray MeanValueIN, MultidimensionalArray MeanValueOT, MultidimensionalArray GradientIN, MultidimensionalArray GradientOT, int ResultIndexOffset, double ResultPreScale) {
            throw new NotImplementedException();
        }



        /// <summary>
        /// Evaluates the gradient of the field;
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="NodeSet">
        /// as usual, the node set;
        /// </param>
        /// <param name="result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="ResultPreScale"/>;
        /// 1st index: cell index minus <paramref name="j0"/>;
        /// 2nd index: node index;
        /// 3rd index: spatial coordinate;
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// an offset for the first index of <paramref name="result"/>;
        /// </param>
        /// <param name="ResultPreScale">
        /// see <paramref name="result"/>
        /// </param>
        public override void EvaluateGradient(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale) {

            int D = this.GridDat.SpatialDimension; // spatial dimension
            int M = NodeSet.NoOfNodes; // number of nodes per cell

            if (result.Dimension != 3)
                throw new ArgumentOutOfRangeException("result", "dimension of result array must be 3");
            if (result.GetLength(1) != M)
                throw new ArgumentOutOfRangeException();
            if (result.GetLength(2) != D)
                throw new ArgumentOutOfRangeException();

            GenericEval(j0, Len, NodeSet, result, ResultCellindexOffset, ResultPreScale,
                DGField.EvaluateGradientInternal,
                ref m_EvaluateGradient_SpeciesEvalBuffer,
                _M => new int[] { 1, _M, D },
                delegate (MultidimensionalArray R, int offset, int m, int SpcInd, MultidimensionalArray[] SR) {
                    for (int d = 0; d < D; d++) {
                        double r = R[offset, m, d] * ResultPreScale;
                        r += SR[SpcInd][0, m, d];
                        R[offset, m, d] = r;
                    }
                });
        }

        MultidimensionalArray[] m_EvaluateGradient_SpeciesEvalBuffer = new MultidimensionalArray[0];


        /// <summary>
        /// evaluates the mean value over a cell;
        /// of course, the mean value doesn't depend on node set or anything like that,
        /// so no information about that has to be provided.
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="ResultPreScale"/>;
        /// 1st index: cell index minus <paramref name="j0"/>;
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// an offset for the first index of <paramref name="result"/>;
        /// </param>
        /// <param name="ResultPreScale">
        /// see <paramref name="result"/>
        /// </param> 
        public override void EvaluateMean(int j0, int Len, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale) {
            throw new NotImplementedException("will come soon");
        }

        public void EvaluateMeanAB(int j0, int Len, MultidimensionalArray result, LevelSetTracker levelSetTracker, int ResultCellindexOffset = 0, double ResultPreScale = 0.0) {
            CellMask cutCells = levelSetTracker.Regions.GetCutCellMask();
            CellMask speciesA = levelSetTracker.Regions.GetSpeciesMask("A");

            for (int j = 0; j < Len; j++) {
                int cell = j0 + j;
                double temp;

                // Take average of species A and B in cut cells
                if (cutCells.Contains(cell)) {
                    double valueA = this.GetSpeciesShadowField("A").GetMeanValue(cell);
                    double valueB = this.GetSpeciesShadowField("B").GetMeanValue(cell);
                    temp = 0.5 * (valueA + valueB);
                } else {
                    if (speciesA.Contains(cell)) {
                        temp = this.GetSpeciesShadowField("A").GetMeanValue(cell);
                    } else {
                        temp = this.GetSpeciesShadowField("B").GetMeanValue(cell);
                    }
                }

                if (ResultPreScale != 0.0) {
                    result[j + ResultCellindexOffset] *= ResultPreScale;
                }
                //result[j + ResultCellindexOffset] = GetMeanValue(j + j0);
                result[j + ResultCellindexOffset] = temp;
            }
        }

        /// <summary>
        /// the DG coordinated of the cut field;
        /// Note that this is a sparse structure, where not all elements are assigned;
        /// </summary>
        /// <remarks>
        /// 
        /// </remarks>
        public override IMatrix Coordinates {
            get {
                return m_Coordinates;
            }
        }

        /// <summary>
        /// guess what?
        /// </summary>
        public override object Clone() {
            XDGField r = new XDGField(m_CCBasis, this.Identification);
            r.m_Coordinates = (FieldStorage)m_Coordinates.Clone();
            //r.m_TrackerVersionCnt = this.m_TrackerVersionCnt;
            r.m_UpdateBehaviour = this.m_UpdateBehaviour;
            return r;
        }

        /// <summary>
        /// guess what?
        /// </summary>
        new public XDGField CloneAs() {
            return (XDGField)Clone();
        }

        FieldStorage m_Coordinates;

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

        /// <summary>
        /// Gets the DG coordinates for a specific species in a specific cell;
        /// </summary>
        /// <param name="jCell">
        /// local cell index
        /// </param>
        /// <param name="specId">
        /// species ID
        /// </param>
        /// <param name="ouput_DGcords4spec">
        /// on exit, the DG coordinates for species <paramref name="specId"/>
        /// </param>
        /// <returns>
        /// true, if the cell <paramref name="jCell"/> containes DOF for species <paramref name="specId"/>,
        /// otherwise false;
        /// </returns>
        public bool GetCoordinates4Species(int jCell, SpeciesId specId, double[] ouput_DGcords4spec) {
            LevelSetTracker trk = m_CCBasis.Tracker;
            ReducedRegionCode rrc;
            int NoOfSpec = trk.Regions.GetNoOfSpecies(jCell, out rrc);
            int SpeciedIdx = m_CCBasis.Tracker.GetSpeciesIndex(rrc, specId);


            if (SpeciedIdx < 0) {
                // species not allocated
                Array.Clear(ouput_DGcords4spec, 0, ouput_DGcords4spec.Length);
                return false;
            } else {
                int NSep = m_CCBasis.DOFperSpeciesPerCell;


                int n0 = NSep * SpeciedIdx;
                for (int n = 0; n < NSep; n++) {
                    ouput_DGcords4spec[n] = this.m_Coordinates[jCell, n + n0];
                }

                return true;
            }
        }

        /*
        /// <summary>
        /// Sets the DG coordinates for a specific species in a specific cell;
        /// </summary>
        /// <param name="jCell">
        /// local cell index
        /// </param>
        /// <param name="specId">
        /// species ID
        /// </param>
        /// <param name="input_DGcords4spec">
        /// the DG coordinates that should be assigned for species <paramref name="specId"/>
        /// </param>
        /// <returns>
        /// true, if the operation was successful, i.e. if the cell <paramref name="jCell"/> containes DOF for species <paramref name="specId"/>,
        /// otherwise false;
        /// </returns>
        public bool SetCoordinates4Species(int jCell, SpeciesId specId, double[] input_DGcords4spec) {
            LevelSetTracker trk = m_CCBasis.Tracker;
            ReducedRegionCode rrc;
            int NoOfSpec = trk.Regions.GetNoOfSpecies(jCell, out rrc);
            int SpeciedIdx = m_CCBasis.Tracker.GetSpeciesIndex(rrc, specId);

            if (SpeciedIdx < 0) {
                // species not allocated
                return false;
            } else {
                int NSep = m_CCBasis.DOFperSpeciesPerCell;

                int n0 = NSep * SpeciedIdx;
                for (int n = 0; n < NSep; n++) {
                    this.m_Coordinates[jCell, n + n0] = input_DGcords4spec[n];
                }

                return true;
            }
        }
        */

        void AutoExtrapolateSpecies(SpeciesId Id, SubGrid oldSpeciesSubGrid) {
            LevelSetTracker LsTrk = m_CCBasis.Tracker;
            CellMask allNearMask = m_CCBasis.Tracker.Regions.GetNearFieldMask(m_CCBasis.Tracker.NearRegionWidth);
            
            var SpeciesField = this.GetSpeciesShadowField(Id);

            CellMask ExtrapolateTo = allNearMask.Intersect(LsTrk.Regions.GetSpeciesMask(Id));
            CellMask ExtrapolateFrom = oldSpeciesSubGrid.VolumeMask;

            SpeciesField.CellExtrapolation(ExtrapolateTo, ExtrapolateFrom);
        }

        /// <summary>
        /// sets all DG coordinates that are associated with the species with
        /// ID - number <paramref name="SpeciesId"/> (see <see cref="LevelSetTracker.GetSpeciesId"/>) 
        /// to 0.0. Note that this method does not clear the common DG coordinates (DG - coordinates
        /// common to all species);
        /// </summary>
        /// <param name="SpeciesId"></param>
        public void ClearSpecies(SpeciesId SpeciesId) {
            {
                // check argument
                try {
                    m_CCBasis.Tracker.GetSpeciesName(SpeciesId);
                } catch (Exception) {
                    throw new ArgumentException("Invalid SpeciesId id.");
                }
            }

            LevelSetTracker trk = m_CCBasis.Tracker;
            int Nsep = m_CCBasis.DOFperSpeciesPerCell;

            int J = this.GridDat.iLogicalCells.Count;
            for (int j = 0; j < J; j++) {
                ReducedRegionCode redRegionCode = ReducedRegionCode.Extract(trk.Regions.m_LevSetRegions[j]);
                int spec_idx = trk.GetSpeciesIndex(redRegionCode, SpeciesId);

                if (spec_idx >= 0) {
                    int n0 = spec_idx * Nsep;

                    for (int n = 0; n < Nsep; n++) {
                        m_Coordinates[j, n0 + n] = 0;
                    }
                }
            }
        }


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

        /// <summary>
        /// Not accurate for cut-cells.
        /// </summary>
        public override double GetMeanValue(int j) {
             if (j < 0 || j >= Basis.GridDat.iLogicalCells.Count)
                throw new ArgumentException($"cell index out of range: got {j}, expecting between 0 and {Basis.GridDat.iLogicalCells.Count}, MPI rank {Basis.GridDat.MpiRank}.", "j");


            //throw new NotImplementedException();

            //jhjkhjhk


            var tracker = m_CCBasis.Tracker;
            int NoOfSpc = tracker.Regions.GetNoOfSpecies(j);
            if(NoOfSpc <= 1) {
                return base.GetMeanValue(j);
            } else {
                double acc = 0.0;
                double denom = 0.0;

                for(int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                    SpeciesId sid = tracker.Regions.GetSpeciesIdFromIndex(j, iSpc);
                    if(tracker.Regions.IsSpeciesPresentInCell(sid, j)) {
                        acc += GetSpeciesShadowField(sid).GetMeanValue(j);
                        denom += 1.0;
                    }
                }

                if(denom > 0)
                    return acc / denom;
                else
                    return 0.0;
            }



            //SubGrid subgrd = tracker.GetCutCellSubGrid();
            //int j_subgrd = subgrd.LocalCellIndex2SubgridIndex[j];

            //if (j_subgrd < 0) {
            //    return base.GetMeanValue(j);
            //} else {
            //    //if (this.GridDat.Cells.IsCellAffineLinear(j)) {
            //    //    double bv = m_Basis.Polynomials[0].Coeff[0];
            //    //    double sc = this.GridDat.Cells.OneOverSqrt_AbsDetTransformation[j];
            //    //    //return (bv * sc * this.Coordinates[j, 0]);
            //        double Vol = this.GridDat.Cells.GetCellVolume(j);

            //        ReducedRegionCode rrc;
            //        int NoOfSpec = tracker.GetNoOfSpecies(j, out rrc);
            //        double Acc = 0;



            //        double facSum = 0;
            //        for (int iSpec = 0; iSpec < NoOfSpec; iSpec++) {
            //            double factor = tracker.GetSpeciesArea(j, iSpec) / Vol;
            //            Debug.Assert(factor < 1.001);
            //            CoordinatesAcc += m_Coordinates[j, NSep * iSpec] * factor;
            //            facSum += factor;
            //        }
            //        Debug.Assert(facSum < 1.001);


            //    //    return (bv * sc * CoordinatesAcc);
            //    //} else {
            //    //    throw new NotImplementedException
            //    //}
            //}
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="j"></param>
        /// <param name="v"></param>
        public override void SetMeanValue(int j, double v) {
            //base.SetMeanValue(j, v);
            throw new NotImplementedException();
        }

        public void SetMeanValueAB(int j, double v) {
            this.GetSpeciesShadowField("A").SetMeanValue(j, v);
            this.GetSpeciesShadowField("B").SetMeanValue(j, v);
        }

        /// <summary>
        /// performs the projection of <paramref name="func"/> for each species.
        /// </summary>
        override public void ProjectField(double alpha, ScalarFunction func, CellQuadratureScheme scheme = null) {
            using (new FuncTrace()) {
                var lsTrk = this.Basis.Tracker;

                foreach (var spc in lsTrk.SpeciesIdS) {
                    var domain = lsTrk.Regions.GetSpeciesSubGrid(spc).VolumeMask;
                    CellQuadratureScheme _scheme;
                    if (scheme == null) {
                        _scheme = new CellQuadratureScheme(UseDefaultFactories: true, domain: domain);
                    } else {
                        if (scheme.Domain != null)
                            domain = domain.Intersect(scheme.Domain);
                        _scheme = new CellQuadratureScheme(UseDefaultFactories: false, domain: domain);
                        scheme.FactoryChain.ForEach(fact_dom => _scheme.AddFactory(fact_dom.RuleFactory, fact_dom.Domain));
                    }

                    this.GetSpeciesShadowField(spc).ProjectField(alpha, func, _scheme);
                }
            }
        }

        /// <summary>
        /// performs the projection of <paramref name="f"/> for each species.
        /// </summary>
        public override void ProjectFunction(double alpha, Func<Vector, double[], int, double> f, CellQuadratureScheme scheme, params DGField[] U) {
            
            using (new FuncTrace()) {
                var lsTrk = this.Basis.Tracker;

                foreach (var spc in lsTrk.SpeciesIdS) {
                    var domain = lsTrk.Regions.GetSpeciesSubGrid(spc).VolumeMask;
                    CellQuadratureScheme _scheme;
                    if (scheme == null) {
                        _scheme = new CellQuadratureScheme(UseDefaultFactories: true, domain: domain);
                    } else {
                        if (scheme.Domain != null)
                            domain = domain.Intersect(scheme.Domain);
                        _scheme = new CellQuadratureScheme(UseDefaultFactories: false, domain: domain);
                        scheme.FactoryChain.ForEach(fact_dom => _scheme.AddFactory(fact_dom.RuleFactory, fact_dom.Domain));
                    }

                    this.GetSpeciesShadowField(spc).ProjectFunction(alpha, f, _scheme, U);
                }
            }
        }



        /// <summary>
        /// Add a constant to this field
        /// </summary>
        public override void AccConstant(double a, CellMask cm) {
            foreach(var spc in this.Basis.Tracker.SpeciesIdS) {
                var sf = this.GetSpeciesShadowField(spc);
                sf.AccConstant(a, cm);
            }
        }


        /// <summary>
        /// symbolic derivation, cell by cell and species by species;
        /// accumulates the derivative of DG field <paramref name="f"/> 
        /// (along the <paramref name="d"/>-th axis) times <paramref name="alpha"/>
        /// to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* \f$ \frac{\partial}{\partial x_d}\f$ <paramref name="f"/>;
        /// </summary>
        /// <param name="f"></param>
        /// <param name="d">
        /// 0 for the x-derivative, 1 for the y-derivative, 2 for the z-derivative
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="f"/>;
        /// </param>
        /// <param name="em">
        /// An optional restriction to the domain in which the derivative is computed (it may, e.g.
        /// be only required in boundary cells, so a computation over the whole domain 
        /// would be a waste of computation power. A proper execution mask for this case would be e.g. 
        /// <see cref="BoSSS.Foundation.Grid.Classic.GridData.BoundaryCells"/>.)
        /// 
        /// if null, the computation is carried out in the whole domain;
        /// </param>
        /// <remarks>
        /// The derivative is calculated by a cell-by-cell (symbolic) derivation of the DG polynomials, therefore the
        /// (effective) DG polynomial degree is one lower than the degree of <paramref name="f"/>;<br/>
        /// In comparison to <see cref="DerivativeByFlux"/>, this method should be much faster,
        /// because no quadrature is involved;<br/>
        /// The field <paramref name="f"/> can be either a <see cref="XDGField"/> or a <see cref="XDGField"/>;
        /// </remarks>
        public override void Derivative(double alpha, DGField f, int d, CellMask em) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                if (this.Basis.Degree < f.Basis.Degree - 1)
                    throw new ArgumentException("cannot compute derivative because of incompatible basis functions.", "f");
                if (f.Basis.GetType() != this.Basis.GetType())
                    throw new ArgumentException("cannot compute derivative because of incompatible basis functions.", "f");

                SpeciesId[] SpcIds = this.Basis.Tracker.SpeciesIdS.ToArray();
                DGField[] f_spc = new DGField[SpcIds.Length];
                SpeciesShadowField[] this_sh = new SpeciesShadowField[SpcIds.Length];
                CellMask[] emS = new CellMask[SpcIds.Length];

                for (int iSpc = 0; iSpc < SpcIds.Length; iSpc++) {
#if DEBUG
                    // avoid watchdog errors in DEBUG caused by runtime difference
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
#endif
                    SpeciesId SpcId = SpcIds[iSpc];
                    if (f is XDGField)
                        f_spc[iSpc] = ((XDGField)f).GetSpeciesShadowField(SpcId);
                    else
                        f_spc[iSpc] = f;
                    this_sh[iSpc] = this.GetSpeciesShadowField(SpcId);

#if DEBUG
                    // avoid watchdog errors in DEBUG caused by runtime difference
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
#endif
                    CellMask _em = this.Basis.Tracker.Regions.GetSpeciesSubGrid(SpcId).VolumeMask;
                    if (em != null)
                        emS[iSpc] = _em.Intersect(em);
                    else
                        emS[iSpc] = _em;


                }

                for (int iSpc = 0; iSpc < SpcIds.Length; iSpc++) {
#if DEBUG
                    // avoid watchdog errors in DEBUG caused by runtime difference
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
#endif
                    this_sh[iSpc].Derivative(alpha, f_spc[iSpc], d, emS[iSpc]);
                }
            }
        }


        /// <summary>
        /// equal to implementation in base class, see
        /// <see cref="DGField.DerivativeByFlux(double,DGField,int,SubGrid,SubGridBoundaryModes)"/>;
        /// The DG coordinates of all species are computed equally;
        /// </summary>
        public override void DerivativeByFlux(double alpha, DGField f, int d,
            SubGrid optionalSubGrid = null,
            SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary) {

            foreach (var SpcId in this.Basis.Tracker.SpeciesIdS) {
                DGField _f;
                if (f is XDGField)
                    _f = ((XDGField)f).GetSpeciesShadowField(SpcId);
                else
                    _f = f;

                SubGrid sgrd = this.Basis.Tracker.Regions.GetSpeciesSubGrid(SpcId);
                if (optionalSubGrid != null)
                    sgrd = new SubGrid(sgrd.VolumeMask.Intersect(optionalSubGrid.VolumeMask));

                this.GetSpeciesShadowField(SpcId).DerivativeByFlux(alpha, _f, d, sgrd, bndMode);
            }
        }


        /// <summary>
        /// initializes this field to be a copy of another field
        /// </summary>
        /// <param name="other">
        /// must be a <see cref="XDGField"/>-object
        /// </param>
        public override void CopyFrom(DGField other) {
            if (!(other is XDGField))
                throw new ApplicationException("unable to copy, because the other field is no XDG field");

            if (!other.Basis.Equals(this.Basis))
                throw new ApplicationException("unable to copy, because the DG polynomial basis of other field is different.");

            XDGField o = other as XDGField;
            FieldStorage oc = o.m_Coordinates;
            FieldStorage tc = this.m_Coordinates;

            tc.CopyFrom(oc);
        }

        /// <summary>
        /// executes <see cref="DGField.ProjectProduct(double,DGField,DGField,CellMask,bool)"/>, for each species.
        /// </summary>
        public override void ProjectProduct(double alpha, DGField a, DGField b, CellMask em, bool accumulateResult) {


            var tracker = this.m_CCBasis.Tracker;
            if (a is XDGField) {
                if (!object.ReferenceEquals(((XDGField)a).Basis.Tracker, tracker))
                    throw new ArgumentException("XDG field 'a' may not be assigned to a different tracker.", "a");
            }
            if (b is XDGField) {
                if (!object.ReferenceEquals(((XDGField)b).Basis.Tracker, tracker))
                    throw new ArgumentException("XDG field 'b' may not be assigned to a different tracker.", "b");
            }

            foreach (SpeciesId sp in tracker.SpeciesIdS) {
                var Grd = tracker.Regions.GetSpeciesSubGrid(sp);

                DGField a_sp;
                if (a is XDGField) {
                    a_sp = ((XDGField)a).GetSpeciesShadowField(sp);
                } else {
                    a_sp = a;
                }

                DGField b_sp;
                if (b is XDGField) {
                    b_sp = ((XDGField)b).GetSpeciesShadowField(sp);
                } else {
                    b_sp = b;
                }

                //CellMask ExeMask;
                //if (em == null)
                //    ExeMask = Grd.VolumeMask;
                //else
                //    ExeMask = CellMask.Intersect(Grd.VolumeMask, em);

                this.GetSpeciesShadowField(sp).ProjectProduct(alpha, a_sp, b_sp, em, accumulateResult);
            }
        }

        /// <summary>
        /// executes <see cref="DGField.ProjectQuotient(double,DGField,DGField,CellMask,bool)"/>, for each species.
        /// </summary>
        public override void ProjectQuotient(double alpha, DGField a, DGField b, CellMask em, bool accumulateResult) {
            var tracker = this.m_CCBasis.Tracker;
            if (a is XDGField) {
                if (!object.ReferenceEquals(((XDGField)a).Basis.Tracker, tracker))
                    throw new ArgumentException("XDG field 'a' may not be assigned to a different tracker.", "a");
            }
            if (b is XDGField) {
                if (!object.ReferenceEquals(((XDGField)b).Basis.Tracker, tracker))
                    throw new ArgumentException("XDG field 'b' may not be assigned to a different tracker.", "b");
            }

            foreach (SpeciesId sp in tracker.SpeciesIdS) {
                var Grd = tracker.Regions.GetSpeciesSubGrid(sp);

                DGField a_sp;
                if (a is XDGField) {
                    a_sp = ((XDGField)a).GetSpeciesShadowField(sp);
                } else {
                    a_sp = a;
                }

                DGField b_sp;
                if (b is XDGField) {
                    b_sp = ((XDGField)b).GetSpeciesShadowField(sp);
                } else {
                    b_sp = b;
                }

                CellMask ExeMask;
                if (em == null)
                    ExeMask = Grd.VolumeMask;
                else
                    ExeMask = CellMask.Intersect(Grd.VolumeMask, em);

                this.GetSpeciesShadowField(sp).ProjectQuotient(alpha, a_sp, b_sp, em, accumulateResult);
            }
        }


        /// <summary>
        /// executes <see cref="DGField.ProjectAbs(double, CellMask, DGField[])"/>, for each species.
        /// </summary>
        public override void ProjectAbs(double alpha, CellMask em, params DGField[] vec) {

            var tracker = this.m_CCBasis.Tracker;
            int D = vec.Length;
            for (int d = 0; d < D; d++) {
                DGField vec_comp = vec[d];
                if (vec_comp is XDGField) {
                    if (!object.ReferenceEquals(((XDGField)vec_comp).Basis.Tracker, tracker))
                        throw new ArgumentException("XDG field component may not be assigned to a different tracker.");
                }
            }

            foreach (SpeciesId sp in tracker.SpeciesIdS) {
                var Grd = tracker.Regions.GetSpeciesSubGrid(sp);

                DGField[] vec_sp = new DGField[D];
                for (int d = 0; d < D; d++) {
                    DGField vec_comp = vec[d];
                    if (vec_comp is XDGField) {
                        vec_sp[d] = ((XDGField)vec_comp).GetSpeciesShadowField(sp);
                    } else {
                        vec_sp[d] = vec_comp;
                    }
                }

                CellMask ExeMask;
                if (em == null)
                    ExeMask = Grd.VolumeMask;
                else
                    ExeMask = CellMask.Intersect(Grd.VolumeMask, em);

                this.GetSpeciesShadowField(sp).ProjectAbs(alpha, em, vec_sp);
            }
        }

        /// <summary>
        /// todo
        /// </summary>
        public override void ProjectPow(double alpha, DGField f, double pow, CellMask em) {
            throw new NotImplementedException("todo");
        }


        /// <summary>
        /// Converts a XDG-Field to a SinglePhaseField, while ignoring any negative Side-Effects like Gibbs' phenomena.
        /// Can be executed on a subgrid only.
        /// <param name="BasisDegreeMultiplicator">
        /// Factor for the Basis Degree of the returned Field.
        /// A higher Polynomial Degree may be chosen to better represent kinks and jumps.
        /// </param>
        /// <param name="cellMask">
        /// Returns Field only on a subgrid
        /// </param>
        /// </summary>
        public SinglePhaseField ProjectToSinglePhaseField(int BasisDegreeMultiplicator = 2, CellMask cellMask = null) {
            var QRs = GridDat.iGeomCells.RefElements.Select(Kref => Kref.GetQuadratureRule(this.Basis.Degree * BasisDegreeMultiplicator));
            SinglePhaseField FieldReturn = new SinglePhaseField(this.Basis.NonX_Basis, this.Identification);
            CellQuadratureScheme CQS = new CellQuadratureScheme(false, cellMask).AddFixedRuleS(QRs);
            FieldReturn.ProjectField(1.0,
                this.Evaluate,
                CQS);
            return FieldReturn;
        }


        public double L2NormAllSpecies(CellMask cm = null) {
            

            double acc = 0;
            foreach (SpeciesId spc in this.Basis.Tracker.SpeciesIdS) {
                acc += L2NormSpecies(spc, cm).Pow2();
            }

            return acc.Sqrt();
        }


        /// <summary>
        /// <see cref="L2NormSpecies(SpeciesId, CellMask)"/>
        /// </summary>
        public double L2NormSpecies(string spcNmn, CellMask cm = null) {
            SpeciesId spc = this.Basis.Tracker.GetSpeciesId(spcNmn);
            return L2NormSpecies(spc,cm);
        }

        /// <summary>
        /// Canonical L2 Norm in the XDG space, i.e.
        /// ```math
        ///   \left( \int_{\Omega \cap \mathfrak{ s} \cap \text{ CM} }
        ///        u(\vec{ x})
        ///   \text{dV} \right)^{1/2},
        /// ```
        /// where 
        /// $`\mathfrak{s}`$ denotes the domain of species <paramref name="spc"/> and
        /// $`\text{ CM} }`$ denotes the optional cell mask <paramref name="cm"/>. 
        /// </summary>
        /// <remarks>
        /// The foundation for computing the L2-norm in cut cells is the following relation:
        /// For any arbitrary (i.e. non-orthonormal) DG or XDG basis,
        /// the norm of a Field 
        /// $`u = sum_{j} \phi_{j} \tilde{u}_{j}`$ is given as:
        /// ```math
        ///         \left\| u \right\|_{L^2}^2 =
        ///         (u, u) = 
        ///     \int_\Omega 
        ///         \left( \sum_{j} \phi_{j} \tilde{u}_{j} \right) 
        ///         \left( \sum_{l} \phi_{l} \tilde{u}_{l} \right) 
        ///     dV = 
        ///     \sum_{j l} \tilde{u}_{j} \tilde{u}_{l} ( \phi_{j}, \phi_{l} )
        ///     =
        ///       \tilde{u}^T M \tilde{u},
        /// ```
        /// where $`M `$ denotes the mass matrix ($` M_{j l} = ( \phi_ { j}, \phi_ { l} )  `$).
        /// </remarks>
        public double L2NormSpecies(SpeciesId spc, CellMask cm = null) {
            int deg = this.Basis.Degree;
            var trk = this.Basis.Tracker;
            int N = this.Basis.NonX_Basis.Length; // DOFs per cell per species.
            double acc = 0; // accumulator for the norm

            if (cm == null)
                cm = trk.Regions.GetSpeciesMask(spc);
            
            if (cm.MaskType != MaskType.Logical)
                throw new ArgumentException();


            double[] GetCoords(int j) {
                int iSpc = trk.Regions.GetSpeciesIndex(spc, j);
                if (iSpc < 0)
                    return null;
                else
                    return this.Coordinates.GetRowPart(j, iSpc * N, N);
            }


            // obtain the mass matrix required for the non-orthonormal cells
            // =============================================================

            var AvailOrders = trk.GetCachedOrders().Where(order => order >= 2 * deg);
            int order2Pick = AvailOrders.Any() ? AvailOrders.Min() : 2 * deg;

            //Console.WriteLine(this.Identification + ": available XDG quad orders: " + AvailOrders.ToConcatString("", ",", ";") + " using order: " + order2Pick);

            var MMF = trk.GetXDGSpaceMetrics(spc, order2Pick).MassMatrixFactory;
            var MMblox = MMF.GetMassMatrixBlocks(this.Basis.NonX_Basis, spc);



            // 1st, do all the cut cells with non-id mass matrix
            // =================================================

            BitArray cmMask = cm.GetBitMask().CloneAs(); // we need to modify the mask
            double[] tmp = new double[N];

            for (int iSub = 0; iSub < MMblox.jSub2jCell.Length; iSub++) {
                int jCell = MMblox.jSub2jCell[iSub];
                if (cmMask[jCell] == false)
                    continue; // not in mask
                
                double[] Coords = GetCoords(jCell);
                if (Coords == null)
                    continue; // species not present in cell; no contribution.

                var MM_j = MMblox.MassMatrixBlocks.ExtractSubArrayShallow(new int[] { iSub, 0, 0 }, new int[] { iSub - 1, N - 1, N - 1 });

                MM_j.GEMV(1.0, Coords, 0.0, tmp);
                double res_j = Coords.InnerProd(tmp);
                acc += res_j;

                cmMask[jCell] = false; // mark that we have taken care of this cell
            }

            // 2nd pass, do all normal cells
            // =============================

            foreach(int jCell in cm.ItemEnum) {
                if (cmMask[jCell] == false)
                    continue; // already done.

                double[] Coords = GetCoords(jCell);

                if (Coords == null)
                    continue; // species not present in cell; no contribution.
                // in this cell, we have an orthonormal basis, i.e. the mass matrix is the identity

                acc += Coords.L2NormPow2();
            }


            // return
            // ======

            acc = acc.MPISum().Sqrt();
            return acc;
        }


        ///// <summary>
        ///// L2 Error; for the cut cells, some precise quadrature is used
        ///// </summary>
        //public override double L2Error(ScalarFunction function, CellQuadratureScheme qr = null) {
        //    using (new FuncTrace()) {
        //        if (qr == null)
        //            qr = new CellQuadratureScheme();
        //        if (qr.FactoryChain.Count() <= 0) {
        //            var splx = this.m_context.Grid.GridSimplex;
        //            var trk = this.Basis.Tracker;

        //            qr.AddStandardRule(splx, (CellMask)null)
        //              .AddFixedRule(splx, trk.RecommendedVolQR, trk.GetCutCellSubGrid().VolumeMask);

        //        }

        //        return base.L2Error(function, qr);
        //    }
        //}

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
        int m_TrackerVersion;

        /// <summary>
        /// Updated the data structure of this cut-cell DG field to reflect the
        /// latest status of the level set.
        /// </summary>
        /// <param name="levelSetStatus"></param>
        public void OnNext(LevelSetTracker.LevelSetRegions levelSetStatus) {
            int J = this.GridDat.iLogicalCells.Count;
            LevelSetTracker trk = m_CCBasis.Tracker;

            int oldTrackerVersion = m_TrackerVersion;
            m_TrackerVersion = trk.VersionCnt;

            m_Coordinates.BeginResize(m_CCBasis.MaximalLength);

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

                double[] coordsFull = new double[m_CCBasis.MaximalLength];
                ushort[] OldRegionCode = levelSetStatus.m_LevSetRegions_b4Update;
                ushort[] NewRegionCode = levelSetStatus.RegionsCode;
                Debug.Assert(!object.ReferenceEquals(OldRegionCode, NewRegionCode));

                int Nsep = m_CCBasis.DOFperSpeciesPerCell;
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
                        m_Coordinates.Resize(j, m_CCBasis.GetLength(j));

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
                        m_Coordinates.Resize(j, m_CCBasis.GetLength(j)); // need to call resize in every case
                    }
                }
            } else {
                // just allocate/free memory
                // +++++++++++++++++++++++++

                for (int j = 0; j < J; j++) {
                    int l = m_CCBasis.GetLength(j);
                    m_Coordinates.Resize(j, l);
                }

            }

            m_Coordinates.FinishResize();

            if (m_UpdateBehaviour == BehaveUnder_LevSetMoovement.JustReallocate)
                // the DOF in the near band will be crap anyway...
                m_Coordinates.Clear();

            //m_TrackerVersionCnt = levelSetStatus.Version;

            // update MPI buffer size
            // ======================
            int size = this.GridDat.CellPartitioning.MpiSize;
            if (size > 1) {
                if (m_MPIRecvBufSize == null)
                    m_MPIRecvBufSize = new int[size];
                if (m_MPISendBufSize == null)
                    m_MPISendBufSize = new int[size];

                for (int p = 0; p < size; p++) {
                    // send list
                    {
                        int[] senditems = this.GridDat.iParallel.SendCommLists[p];
                        if (senditems != null) {
                            int L = senditems.Length;

                            int sz = 0;
                            for (int l = 0; l < L; l++)
                                sz += m_CCBasis.GetLength(senditems[l]);
                            m_MPISendBufSize[p] = sz;
                        } else {
                            m_MPISendBufSize[p] = int.MinValue;
                        }
                    }

                    // receive list
                    {
                        int L = this.GridDat.iParallel.RcvCommListsNoOfItems[p];
                        if (L > 0) {
                            int j0 = this.GridDat.iParallel.RcvCommListsInsertIndex[p];
                            L += j0;
                            int sz = 0;
                            for (int j = j0; j < L; j++) {
                                sz += m_CCBasis.GetLength(j);
                            }
                            m_MPIRecvBufSize[p] = sz;
                        } else {
                            m_MPIRecvBufSize[p] = int.MinValue;
                        }
                    }
                }
            }

            // do Extrapolation, if necessary
            if (m_UpdateBehaviour == BehaveUnder_LevSetMoovement.AutoExtrapolate && trk.PopulatedHistoryLength >= 1) {
                foreach (var species in m_CCBasis.Tracker.SpeciesIdS) {
                    AutoExtrapolateSpecies(species, trk.RegionsHistory[0].GetSpeciesSubGrid(species));
                }
            }
        }

        #endregion
    }
}
