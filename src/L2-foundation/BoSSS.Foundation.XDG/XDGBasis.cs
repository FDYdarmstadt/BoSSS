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
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using System.Runtime.Serialization;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// A Basis which depends on one ore more level sets (see <see cref="Tracker"/>).
    /// </summary>
    public partial class XDGBasis : Basis {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="Degree"></param>
        /// <param name="levSetTracker"></param>
        public XDGBasis(LevelSetTracker levSetTracker, int Degree) :
            base(levSetTracker.GridDat, Degree) {

            m_Tracker = levSetTracker;
            base.IsOrthonormal = false;

            // minimal Length
            // ==============
            this.m_MinimalLength = base.Polynomials[0].Count;
            if(base.MaximalLength != base.MinimalLength)
                throw new NotSupportedException();

            // maximal Length
            // ==============

            this.m_MaxNoOfSpecies = m_Tracker.TotalNoOfSpecies;
            this.m_MaximalLength = this.m_MinimalLength*m_MaxNoOfSpecies;


            // caches
            // ======
            m_ByCellEvalCache = new ByCellEvalCache(this);
            m_ByCellEvalGradCache = new ByCellEvalGradCache(this);
        }

        [DataMember]
        int m_MaxNoOfSpecies;

        /// <summary>
        /// returns the DG coordinate index for the first separate DG coordinate
        /// (DG coordinates which are separate for each species) 
        /// for the species with index <paramref name="SpeciesInd"/>.
        /// </summary>
        /// <param name="SpeciesInd"></param>
        /// <returns></returns>
        /// <remarks>
        /// Note that the index of some species may vary from cell to cell: E.g.
        /// consider a case with two species, "A" and "B". In cells which contain
        /// only species "A", the species index for "A" is 0 and the species index for "B" 
        /// is not assigned. <br/>
        /// In cell which contain both, "A" and "B", the index for "A" may be 1 and for "B" it may be 0.
        /// </remarks>
        public int GetSpeciesI0(int SpeciesInd) {
            if (SpeciesInd < 0 || SpeciesInd >= m_MaxNoOfSpecies)
                throw new IndexOutOfRangeException();

            return SpeciesInd * this.m_MinimalLength;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="jCell"></param>
        /// <returns></returns>
        public override int GetLength(int jCell) {
            ushort code = m_Tracker.Regions.m_LevSetRegions[jCell];
            if (code == LevelSetTracker.AllFARplus || code == LevelSetTracker.AllFARminus) {
                // default: uncut cell
                return base.GetLength(jCell);
            } else {
                // cut cell
                ReducedRegionCode dummy;
                return (m_Tracker.Regions.GetNoOfSpecies(jCell, out dummy) * base.GetLength(jCell));
            }
        }



        /// <summary>
        /// the number of degrees-of-freedom for each species, for each cell.
        /// </summary>
        public int DOFperSpeciesPerCell {
            get {
                return this.m_MinimalLength;
            }
        }

        /// <summary>
        /// determines whether this DG basis is a sub-basis (in a vector-space sense)
        /// of another basis <paramref name="other"/>.
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public override bool IsSubBasis(Basis other) {
            XDGBasis o = other as XDGBasis;
            if (o != null) {
                if (!object.ReferenceEquals(o.m_Tracker, this.m_Tracker))
                    // different LevelSet-tracker -> no Subbasis
                    return false;


            } else {
                // a cut-cell - basis cannot be a subbasis of a 'normal' DG basis.
                return false;
            }

            // a (normal) DG basis is also a sub-basis of a cut-cell - basis!
            if (this.Polynomials.Count > other.Polynomials.Count)
                return false;

            // check polynomials
            for(int iKref = this.Polynomials.Count - 1; iKref >= 0; iKref--) {
                if(other.Polynomials[iKref].Count < this.Polynomials[iKref].Count)
                    return false;

                for(int j = this.Polynomials[iKref].Count - 1; j >= 0; j--)
                    if(!this.Polynomials[iKref][j].Equals(other.Polynomials[iKref][j]))
                        return false;
            }

            return true;
        }


        /// <summary>
        /// see <see cref="MinimalLength"/>;
        /// </summary>
        [DataMember]
        int m_MinimalLength;

        /// <summary>
        /// see <see cref="MaximalLength"/>;
        /// </summary>
        [DataMember]
        int m_MaximalLength;

        /// <summary>
        /// This is the maximum number of degrees-of-freedom that can occur
        /// in one cell, i.e. in a cell which contains all species
        /// tracked by <see cref="Tracker"/> at once;
        /// </summary>
        public override int MaximalLength {
            get {
                return m_MaximalLength;
            }
        }

        /// <summary>
        /// The number of basis functions in uncut cells, 
        /// which is equal to the number of basis polynomials
        /// <see cref="BoSSS.Foundation.Basis.Polynomials"/>
        /// </summary>
        public override int MinimalLength {
            get {
                return m_MinimalLength;
            }
        }

        [NonSerialized]
        LevelSetTracker m_Tracker;

        /// <summary>
        /// the used tracker for the level set
        /// </summary>
        public LevelSetTracker Tracker {
            get {
                return m_Tracker;
            }
        }

        /// <summary>
        /// two <see cref="XDGBasis"/>-objects are equal, if their polynomial list
        /// (<see cref="Basis.Polynomials"/>) is equal (equal Guid for each entry),
        /// and if <see cref="DOFperSpeciesPerCell"/> coincide.
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            if (obj.GetType() != typeof(XDGBasis))
                return false;

            XDGBasis other = (XDGBasis)obj;

            if (!object.ReferenceEquals(other.Tracker, this.Tracker))
                return false;

            if(other.Polynomials.Count != this.Polynomials.Count)
                return false;

            if (other.DOFperSpeciesPerCell != this.DOFperSpeciesPerCell)
                return false;

            for(int iKref = other.Polynomials.Count - 1; iKref >= 0; iKref--) {
                if(other.Polynomials[iKref].Count != this.Polynomials[iKref].Count)
                    return false;

                for(int j = other.Polynomials[iKref].Count - 1; j >= 0; j--)
                    if(!this.Polynomials[iKref][j].Equals(other.Polynomials[iKref][j]))
                        return false;
            }

            return true;
        }

        /// <summary>
        /// default implementation;
        /// </summary>
        public override int GetHashCode() {
            return this.MaximalLength + 333;
        }

        [NonSerialized]
        Basis m_NonX_Basis;

        /// <summary>
        /// creates the standard (not extended) DG basis of the same polynomial degree as this one.
        /// </summary>
        public Basis NonX_Basis {
            get {
                if(m_NonX_Basis== null) {
                    m_NonX_Basis = new Basis(this.GridDat, this.Degree);
                }
                return m_NonX_Basis;
            }
        }


        /// <summary>
        /// evaluation of basis function within cells
        /// </summary>
        public override MultidimensionalArray CellEval(NodeSet nodes, int j0, int Len) {
            return m_ByCellEvalCache.GetValue_Cell(nodes, j0, Len);
        }

        [NonSerialized]
        ByCellEvalCache m_ByCellEvalCache;

        abstract class ByCellCache : Caching.CacheLogic_CNs {

            internal ByCellCache(XDGBasis owner)
                : base(owner.GridDat) //
            {
                m_Owner = owner;
            }

            protected XDGBasis m_Owner;

            bool CheckAll(int j0, int Len) {
                if( Len <= 0)
                    return true;

                int l0 = m_Owner.GetLength(j0);
                for (int i = 1; i < Len; i++) {
                    if (m_Owner.GetLength(j0 + i) != l0)
                        return false;
                }
                return true;
            }


            protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {

                Debug.Assert(CheckAll(j0, Len), "basis length must be equal for all cells from j0 to j0+Len");
                
                int l0 = m_Owner.GetLength(j0);
                var trk = m_Owner.Tracker;

                // evaluation of un-cut polynomials
                // ================================
                var BasisValues = NonXEval(NS, j0, Len);
                
                
                // evaluation of XDG 
                // =================

                if (l0 == m_Owner.NonX_Basis.Length) {
                    // uncut region -> evaluation is equal to uncut basis
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++

                    output.Set(BasisValues);

                } else {
                    // modulated basis polynomials
                    // +++++++++++++++++++++++++++

                    int NoOfLevSets = trk.LevelSets.Count;
                    int M = output.GetLength(1); // number of nodes
                    output.Clear();

                    ReducedRegionCode rrc;
                    int NoOfSpecies = trk.Regions.GetNoOfSpecies(j0, out rrc);
                    int Nsep = m_Owner.DOFperSpeciesPerCell;
                    


                    for (int j = 0; j < Len; j++) { // loop over cells


                        int jCell = j + j0;

                        ushort RegionCode = trk.Regions.m_LevSetRegions[jCell];


                        int dist = 0;
                        for (int i = 0; i < NoOfLevSets; i++) {
                            dist = LevelSetTracker.DecodeLevelSetDist(RegionCode, i);
                            if (dist == 0) {
                                m_levSetVals[i] = trk.DataHistories[i].Current.GetLevSetValues(NS, jCell, 1);
                            } else {
                                //_levSetVals[i] = m_CCBasis.Tracker.GetLevSetValues(i, NodeSet, jCell, 1);
                                m_levSetVals[i] = null;
                                m_levSetSign[i] = dist;
                            }
                        }


                        LevelSetSignCode levset_bytecode;
                        int SpecInd;
                        if (dist != 0) {
                            // near - field: 
                            // Basis contains zero-entries
                            // ++++++++++++++++++++++++++++

                            levset_bytecode = LevelSetSignCode.ComputeLevelSetBytecode(m_levSetSign);
                            SpecInd = trk.GetSpeciesIndex(rrc, levset_bytecode);
                        } else {
                            levset_bytecode = default(LevelSetSignCode);
                            SpecInd = 0;
                        }


                        for (int m = 0; m < M; m++) { // loop over nodes

                            if (dist == 0) {
                                // cut cells
                                // Modulation necessary
                                // +++++++++++++++++++++

                                for (int i = 0; i < NoOfLevSets; i++) {
                                    if (m_levSetVals[i] != null)
                                        m_levSetSign[i] = m_levSetVals[i][0, m];
                                }

                                // re-compute species index
                                levset_bytecode = LevelSetSignCode.ComputeLevelSetBytecode(m_levSetSign);
                                SpecInd = trk.GetSpeciesIndex(rrc, levset_bytecode);
                            }


                            // separate coordinates
                            {
                                int n0 = Nsep * SpecInd;
                                for (int n = 0; n < Nsep; n++) { // loop over basis polynomials
                                    //output[j, m, n0 + n] = BasisValues[m, n];
                                    Operation(output, BasisValues, j, m, n0 + n, n);
                                }
                            }
                        }
                    }
                }
            }

            abstract protected void Operation(MultidimensionalArray output, MultidimensionalArray BasisValues, int j, int m, int nA, int nB);

            abstract protected MultidimensionalArray NonXEval(NodeSet NS, int j0, int Len);


            MultidimensionalArray[] m_levSetVals = new MultidimensionalArray[4];
            double[] m_levSetSign = new double[4];

        }

        class ByCellEvalCache : ByCellCache {

            internal ByCellEvalCache(XDGBasis o) : base(o) {
            }


            protected override void Operation(MultidimensionalArray output, MultidimensionalArray BasisValues, int j, int m, int nA, int nB) {
                output[j, m, nA] = BasisValues[j, m, nB];
            }

            protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet N) {
                return MultidimensionalArray.Create(Len, N.NoOfNodes, m_Owner.GetLength(i0));
            }

            protected override MultidimensionalArray NonXEval(NodeSet NodeSet, int j0, int Len) {
                return m_Owner.NonX_Basis.CellEval(NodeSet, j0, Len);
            }
        }


        /// <summary>
        /// evaluation of basis function gradient within cells
        /// </summary>
        public override MultidimensionalArray CellEvalGradient(NodeSet NodeSet, int j0, int Len) {
            return m_ByCellEvalGradCache.GetValue_Cell(NodeSet, j0, Len);
        }

        [NonSerialized]
        ByCellEvalGradCache m_ByCellEvalGradCache;

        class ByCellEvalGradCache : ByCellCache {
            internal ByCellEvalGradCache(XDGBasis o) : base(o) { }

            protected override void Operation(MultidimensionalArray output, MultidimensionalArray BasisValues, int j, int m, int nA, int nB) {
                int D = BasisValues.GetLength(3);
                Debug.Assert(D == m_Owner.GridDat.SpatialDimension);
                Debug.Assert(D == output.GetLength(3));
                for (int d = 0; d < D; d++)
                    output[j, m, nA, d] = BasisValues[j, m, nB, d];
            }

            protected override MultidimensionalArray Allocate(int j0, int Len, NodeSet N) {
                return MultidimensionalArray.Create(Len, N.NoOfNodes, m_Owner.GetLength(j0), N.SpatialDimension);
            }

            protected override MultidimensionalArray NonXEval(NodeSet NodeSet, int j0, int Len) {
                return m_Owner.NonX_Basis.CellEvalGradient(NodeSet, j0, Len);
            }
        }
    }
}
