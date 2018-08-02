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
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation {

    /// <summary>
    /// An orthonormal, modal basis \f$ \phi_{j,n} \f$ for a DG approximation of absolute degree <see cref="Degree"/>.
    /// </summary>
    public partial class Basis {

        /// <summary>
        /// See <see cref="GridDat"/>
        /// </summary>
        protected IGridData m_GridDat;

        /// <summary>
        /// reference to grid of this basis
        /// </summary>
        public IGridData GridDat {
            get {
                return m_GridDat;
            }
        }

        /// <summary>
        /// common basis
        /// </summary>
        public GridData.BasisData Data {
            get {
                return m_GridDat.ChefBasis;
            }
        }


        /// <summary>
        /// True if the basis forms an orthonormal system; otherwise false. Eh' klar...
        /// </summary>
        public bool IsOrthonormal {
            get;
            protected set;
        }


        /// <summary>
        /// creates a new basis
        /// </summary>
        /// <param name="grd">the <see cref="GridData"/> that stores grid information</param>
        /// <param name="degree">highest polynomial degree of basis polynomials</param>
        public Basis(IGridData grd, int degree) {
            this.m_GridDat = grd;
            var PolyList = grd.ChefBasis.GetOrthonormalPolynomials(degree);
            this.Polynomials = PolyList.ToList<PolynomialList>().AsReadOnly();
            this.IsOrthonormal = true;

            this.MinimalLength = Polynomials.Min(pl => pl.Count);
            this.MaximalLength = Polynomials.Max(pl => pl.Count);
            this.Degree = this.Polynomials.Max(pl => pl.MaxAbsoluteDegree);

        }

        /// <summary>
        /// Checks if this <see cref="Basis"/> is a sub-basis of another basis (<paramref name="other"/>).
        /// This is the case if, and only if, the <see cref="Polynomials"/>-array (lets say, of length N) of
        /// this object is equal to the "0 to N-1" subarray of <paramref name="other"/>.
        /// The comparison of the array entries is done by using (see 
        /// <see cref="BoSSS.Foundation.Grid.Polynomial.Equals"/>).
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        virtual public bool IsSubBasis(Basis other) {

            if (!object.ReferenceEquals(other.m_GridDat, this.m_GridDat))
                return false;

            int L = m_GridDat.iGeomCells.RefElements.Length;
            Debug.Assert(this.Polynomials.Count == L);
            Debug.Assert(other.Polynomials.Count == L);

            for (int l = 0; l < L; l++) {
                int N = this.Polynomials[l].Count;

                if (other.Polynomials[l].Count < N)
                    return false;

                for (int i = N - 1; i >= 0; i--)
                    if (!this.Polynomials[l][i].Equals(other.Polynomials[l][i]))
                        return false;
            }

            return true;
        }

        /// <summary>
        /// two <see cref="Basis"/>-objects are equal, if their polynomial list
        /// (<see cref="Polynomials"/>) is equal (equal Guid for each entry).
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            if (obj.GetType() != typeof(Basis))
                return false;

            Basis other = (Basis)obj;

            if (!object.ReferenceEquals(other.m_GridDat, this.m_GridDat))
                return false;

            int L = m_GridDat.iGeomCells.RefElements.Length;
            Debug.Assert(this.Polynomials.Count == L);
            Debug.Assert(other.Polynomials.Count == L);

            for (int l = 0; l < L; l++) {
                int N = this.Polynomials[l].Count;
                if (other.Polynomials[l].Count != N)
                    return false;

                for (int i = N - 1; i >= 0; i--)
                    if (!this.Polynomials[l][i].Equals(other.Polynomials[l][i]))
                        return false;
            }


            return true;
        }

        /// <summary>
        /// default implementation;
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        /// <summary>
        /// the DG basis functions for the reference elements
        /// <list type="bullet">
        ///     <item>1st index: reference element index (see <see cref="BoSSS.Foundation.Grid.GridCommons.RefElements"/>)</item>
        ///     <item>2nd index: Polynomial index</item>
        /// </list>
        /// </summary>
        public IList<PolynomialList> Polynomials {
            get;
            private set;
        }


        /// <summary>
        /// returns the number of basis functions in the cell <paramref name="jCell"/>;
        /// </summary>
        /// <param name="jCell">a local cell index;</param>
        virtual public int GetLength(int jCell) {
            int jGeom = this.GridDat.GetGeometricCellIndices(jCell).First();
            int iKref = this.GridDat.iGeomCells.GetRefElementIndex(jGeom);
            int N = this.Polynomials[iKref].Count;
#if DEBUG
            foreach (int jG in this.GridDat.GetGeometricCellIndices(jCell)) {
                int _iKref = this.GridDat.iGeomCells.GetRefElementIndex(jGeom);
                int _N = this.Polynomials[_iKref].Count;
                Debug.Assert(N == _N);
            }
#endif
            return N;
        }

        /// <summary>
        /// Returns at list of indices into <see cref="Polynomials"/>[$iKref], where
        /// $iKref ist the reference element index for cell <paramref name="jCell"/>.
        /// In particular, returns the indices of all polynomials with the given
        /// <paramref name="degree"/>.
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="degree"></param>
        /// <returns></returns>
        public IEnumerable<int> GetPolynomialIndicesForDegree(int jCell, int degree) {
            int jGeom = this.GridDat.GetGeometricCellIndices(jCell).First();
            int iKref = this.GridDat.iGeomCells.GetRefElementIndex(jGeom);

            for (int index = 0; index < Polynomials[iKref].Count; index++) {
                if (Polynomials[iKref][index].AbsoluteDegree == degree) {
                    yield return index;
                }
            }
        }

        /// <summary>
        /// highest polynomial degree of all basis polynomials
        /// </summary>
        public int Degree {
            get;
            private set;
        }

        /// <summary>
        /// the minimal number of basis polynomials (per cell) over all cells;
        /// </summary>
        public virtual int MinimalLength {
            get;
            private set;
        }

        /// <summary>
        /// the maximum number of basis polynomials (per cell) over all cells;
        /// </summary>
        public virtual int MaximalLength {
            get;
            private set;
        }

        /// <summary>
        /// if <see cref="MinimalLength"/>==<see cref="MaximalLength"/>,
        /// this value; otherwise an exception.
        /// </summary>
        public int Length {
            get {
                if (MinimalLength != MaximalLength)
                    throw new NotSupportedException("not supported on variable-length basis objects");

                return MaximalLength;
            }
        }


        /// <summary>
        /// Evaluates all polynomials in reference space
        /// in this basis (see <see cref="Polynomials"/>);
        /// </summary>
        /// <param name="Ns"></param>
        /// <returns>
        /// <list type="bullet">
        ///   <item>1st index: node index</item>
        ///   <item>2nd index: polynomial index</item>
        /// </list>
        /// </returns>
        public virtual MultidimensionalArray Evaluate(NodeSet Ns) {
            int iKref = Ns.GetVolumeRefElementIndex(this.GridDat);
            int N = this.Polynomials[iKref].Count;
            MultidimensionalArray Values = this.GridDat.ChefBasis.BasisValues.GetValues(Ns, this.Degree);

            if (Values.GetLength(1) == N)
                return Values;
            else
                return Values.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Ns.NoOfNodes - 1, N - 1 });
        }


        /// <summary>
        /// Evaluates the gradient of all polynomials in reference space in this basis (see <see cref="Polynomials"/>);
        /// </summary>
        /// <param name="Ns"></param>
        /// <returns>
        /// <list type="bullet">
        ///   <item>1st index: node index</item>
        ///   <item>2nd index: polynomial index</item>
        ///   <item>3rd index: spatial dimension</item>
        /// </list>
        /// </returns>
        public MultidimensionalArray EvaluateGradient(NodeSet Ns) {
            int iKref = Ns.GetVolumeRefElementIndex(this.GridDat);
            int N = this.Polynomials[iKref].Count;
            int D = this.GridDat.SpatialDimension;
            MultidimensionalArray Values = this.GridDat.ChefBasis.BasisGradientValues.GetValues(Ns, this.Degree);

            if (Values.GetLength(1) == N)
                return Values;
            else
                return Values.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Ns.NoOfNodes - 1, N - 1, D - 1 });
        }


        /// <summary>
        /// A vectorized evaluation for the 2nd derivative, i.e. the Hessian, in reference coordinates: \f$ \partial^2_{\vec{\xi}} \phi_n \f$.
        /// </summary>
        /// <param name="Ns"></param>
        /// <returns>
        /// An array containing the value of the second derivatives of \f$ \phi_n \f$ at the nodes <paramref name="Ns"/>.
        /// <list type="bullet">
        ///   <item>1st index: node index \f$ m \f$</item>
        ///   <item>2nd index: polynomial index \f$ n \f$</item>
        ///   <item>3rd index: spatial direction of 1st derivation, \f$ k \f$ </item>
        ///   <item>4th index: spatial direction of 2nd derivation, \f$ l \f$ </item>
        /// </list>
        /// The \f$ (m,n,k,l) \f$-th entry is equal to 
        /// \f$ \frac{\partial}{\partial \xi_k} \frac{\partial}{\partial \xi_l} \phi_n (\vec{\xi}_m) \f$,
        /// where \f$ \vec{\xi}_m \f$ is the \f$ m \f$-th vector in the nodeset <paramref name="Ns"/>.
        /// </returns>        
        public MultidimensionalArray Evaluate2ndDeriv(NodeSet Ns) {

            int iKref = Ns.GetVolumeRefElementIndex(this.GridDat);
            int N = this.Polynomials[iKref].Count;
            int D = this.GridDat.SpatialDimension;
            MultidimensionalArray Values = this.GridDat.ChefBasis.BasisHessianValues.GetValues(Ns, this.Degree);

            if (Values.GetLength(1) == N)
                return Values;
            else
                return Values.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Ns.NoOfNodes - 1, N - 1, D - 1, D - 1 });

        }

        /// <summary>
        /// A vectorized evaluation for the 2nd derivative, i.e. the Hessian, in physical coordinates: \f$ \partial^2_{\vec{x}} \phi_{j n} \f$.
        /// </summary>
        /// <param name="Nodes"></param>
        /// <param name="j0">Index of first cell in the vector</param>
        /// <param name="Len">Length of the vector</param>
        /// <returns>
        ///  - 1st index: cell
        ///  - 2nd index: node index <em>m</em>
        ///  - 3rd index: polynomial index <em>n</em>
        ///  - 4th index: spatial direction of 1st derivation, <em>k</em>
        ///  - 5th index: spatial direction of 2nd derivation, <em>l</em>
        /// So, the entry [m,n,k,l] =
        /// \f$ 
        /// \frac{\partial}{\partial x_k} \frac{\partial}{\partial x_l} \phi_n (\vec{x}_m)
        ///  \f$, where \f$ \vec{x}_m \f$ is the
        /// <em>m</em>-th vector in the nodeset #<paramref name="Nodes"/>.
        /// </returns>
        public MultidimensionalArray CellEval2ndDeriv(NodeSet Nodes, int j0, int Len) {

            int iKref = Nodes.GetVolumeRefElementIndex(this.GridDat);
            int N = this.Polynomials[iKref].Count;
            int D = this.GridDat.SpatialDimension;
            MultidimensionalArray Values = this.GridDat.ChefBasis.CellBasisHessianValues.GetValue_Cell(Nodes, j0, Len, this.Degree);

            if (Values.GetLength(2) == N)
                return Values;
            else
                return Values.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0, 0 }, new int[] { Len - 1, Nodes.NoOfNodes - 1, N - 1, D - 1, D - 1 });

        }

        /// <summary>
        /// Evaluates all polynomials in this basis (see <see cref="Polynomials"/>) within specified cells.
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="j0">index of first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns>
        /// <list type="bullet">
        ///   <item>1st index: cell index</item>
        ///   <item>2nd index: node index</item>
        ///   <item>3rd index: polynomial index</item>
        /// </list>
        /// </returns>
        public virtual MultidimensionalArray CellEval(NodeSet nodes, int j0, int Len) {
            var ret = this.GridDat.ChefBasis.CellBasisValues.GetValue_Cell(nodes, j0, Len, this.Degree);
            if (ret.GetLength(0) != Len) {
                int M = ret.GetLength(1);
                int N = ret.GetLength(2);

                ret = ret.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, M - 1, N - 1 });
            }
            return ret;
        }

        /// <summary>
        /// Evaluates all polynomials in this basis (see <see cref="Polynomials"/>) at specified edges.
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="e0">index of first edge to evaluate</param>
        /// <param name="Len">number of edges to evaluate</param>
        /// <returns>
        /// <list type="bullet">
        ///   <item>1st index: cell index</item>
        ///   <item>2nd index: node index</item>
        ///   <item>3rd index: polynomial index</item>
        /// </list>
        /// </returns>
        public virtual Tuple<MultidimensionalArray, MultidimensionalArray> EdgeEval(NodeSet nodes, int e0, int Len) {
            var ret = this.GridDat.ChefBasis.CellBasisValues.GetValue_EdgeDV(nodes, e0, Len, this.Degree);
            Debug.Assert(Len == ret.Item1.GetLength(0));
            Debug.Assert(nodes.NoOfNodes == ret.Item1.GetLength(1));
            Debug.Assert(Len == ret.Item2.GetLength(0));
            Debug.Assert(nodes.NoOfNodes == ret.Item2.GetLength(1));

            if (ret.Item1.GetLength(0) != Len) {
                int M = ret.Item1.GetLength(1);
                int N = ret.Item1.GetLength(2);


                ret = new Tuple<MultidimensionalArray, MultidimensionalArray>(
                    ret.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, M - 1, N - 1 }),
                    ret.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, M - 1, N - 1 }));
            }
            return ret;
        }



        /// <summary>
        /// Evaluates the gradient of all polynomials in this basis (see <see cref="Polynomials"/>) within specified cells.
        /// </summary>
        /// <param name="NS"></param>
        /// <param name="j0">index of first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns>
        /// <list type="bullet">
        ///   <item>1st index: cell index</item>
        ///   <item>2nd index: node index</item>
        ///   <item>3rd index: polynomial index</item>
        ///   <item>4th index: spatial dimension</item>
        /// </list>
        /// </returns>
        public virtual MultidimensionalArray CellEvalGradient(NodeSet NS, int j0, int Len) {
            var ret = this.GridDat.ChefBasis.CellBasisGradientValues.GetValue_Cell(NS, j0, Len, this.Degree);
            if (ret.GetLength(0) != Len) {
                int M = ret.GetLength(1);
                int N = ret.GetLength(2);
                int D = ret.GetLength(3);

                ret = ret.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Len - 1, M - 1, N - 1, D - 1 });
            }
            return ret;
        }

        /// <summary>
        /// Evaluates the gradient of all polynomials in this basis (see <see cref="Polynomials"/>) at specified edges.
        /// </summary>
        /// <param name="NS"></param>
        /// <param name="j0">index of first edge to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns>
        /// For each 'in'- and 'out'-cell, arrays;<br/>
        /// <list type="bullet">
        ///   <item>1st index: cell index</item>
        ///   <item>2nd index: node index</item>
        ///   <item>3rd index: polynomial index</item>
        ///   <item>4th index: spatial dimension</item>
        /// </list>
        /// </returns>
        public virtual Tuple<MultidimensionalArray, MultidimensionalArray> EdgeEvalGradient(NodeSet NS, int j0, int Len) {
            var ret = this.GridDat.ChefBasis.CellBasisGradientValues.GetValue_EdgeDV(NS, j0, Len, this.Degree);
            if (ret.Item1.GetLength(0) != Len) {
                int M = ret.Item1.GetLength(1);
                int N = ret.Item1.GetLength(2);
                int D = ret.Item1.GetLength(3);

                ret = new Tuple<MultidimensionalArray, MultidimensionalArray>(
                    ret.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Len - 1, M - 1, N - 1, D - 1 }),
                    ret.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Len - 1, M - 1, N - 1, D - 1 }));
            }
            return ret;
        }


        /// <summary>
        /// This class contains all necessary information to recreate a Basis.
        /// </summary>
        [Serializable]
        [DataContract]
        public class BasisInitializer : Initializer<Basis> {

            /// <summary>
            /// The unique id of the grid.
            /// </summary>
            [DataMember]
            public Guid GridGuid;

            /// <summary>
            /// The degree of the basis
            /// </summary>
            [DataMember]
            public int Degree;

            /// <summary>
            /// Instantiates the basis described by this object.
            /// </summary>
            /// <param name="c"></param>
            /// <returns></returns>
            public override Basis Initialize(IInitializationContext c) {
                Basis basis;
                if (c.TryGetValue(this, out basis))
                    return basis;

                if (c.GridData == null)
                    throw new ArgumentException();
                if (!c.GridData.GridID.Equals(this.GridGuid))
                    throw new ArgumentException("Wrong grid.");

                Basis bb = new Basis(c.GridData, this.Degree);
                myInstance = bb;
                c.Add(this, bb);
                return bb;
            }

            /// <summary>
            /// Compares the given object <paramref name="other"/> with respect
            /// to the <see cref="GridGuid"/> and the <see cref="Degree"/>.
            /// </summary>
            /// <param name="other"></param>
            /// <returns></returns>
            public override bool Equals(Initializer<Basis> other) {
                BasisInitializer initializer = other as BasisInitializer;
                return (initializer != null)
                    && (initializer.GridGuid.Equals(this.GridGuid))
                    && (initializer.Degree == this.Degree);
            }

            /// <summary>
            /// Computes a hash code based on <see cref="GridGuid"/> and
            /// <see cref="Degree"/>.
            /// </summary>
            public override int GetHashCode() {
                // http://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
                int hash = 9311; // a prime number
                hash += 319993 * GridGuid.GetHashCode();
                hash += 319993 * this.Degree;
                return hash;
            }

            /// <summary>
            /// An instance of the basis of described by this object; prevents
            /// garbage collection.
            /// </summary>
            [NonSerialized]
            private Basis myInstance;
        }

        /// <summary>
        /// See <see cref="Initializer"/>
        /// </summary>
        private BasisInitializer m_Initializer;

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. Note
        /// that it is essential that this member always returns the SAME
        /// object (reference-equals)!
        /// </summary>
        public virtual BasisInitializer Initializer {
            get {
                if (m_Initializer == null) {
                    m_Initializer = new BasisInitializer() {
                        GridGuid = this.GridDat.GridID,
                        Degree = this.Degree
                    };
                }
                return m_Initializer;
            }
        }


    }
}
