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
using System.Linq;
using BoSSS.Foundation.Grid;
using System.Diagnostics;

namespace BoSSS.Foundation {

    ///// <summary>
    ///// delegate used by the constructor of <see cref="VectorField{T}"/>;
    ///// One implementation is <see cref="SinglePhaseField.Factory"/>
    ///// </summary>
    ///// <param name="b"></param>
    ///// <param name="id"></param>
    ///// <returns></returns>
    //public delegate T FieldFactory<T>(Basis b, string id) where T : DGField;

    /// <summary>
    /// A vector field (stage-1 tensor field) as composition
    /// of scalar fields (<see cref="DGField"/>);
    /// </summary>
    public class VectorField<T> : ICloneable, IList<T> where T : DGField {

        /// <summary>
        /// performs <see cref="DGField.CheckForNanOrInf"/> for each component;
        /// </summary>
        public long CheckForNanOrInf(bool CheckForInf, bool CheckForNan, bool ExceptionIfFound) {
            long ret = -1;
            foreach (DGField f in this.m_Components) {
                ret = f.CheckForNanOrInf(CheckForInf, CheckForNan, ExceptionIfFound);
                if (ret >= 0)
                    return ret;
            }

            return ret;
        }

        /// <summary>
        /// sets all components of the vector field to 0.0
        /// </summary>
        public void Clear() {
            foreach (DGField f in m_Components)
                f.Clear();
        }

        /// <summary>
        /// sets all coordinates of cells in <paramref name="cellMask"/> of this field to 0;
        /// </summary>
        public void Clear(CellMask cellMask) {
            foreach (DGField f in m_Components)
                f.Clear(cellMask);
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="D">
        /// number of components, i.d. dimension (see <see cref="Dim"/>)
        /// </param>
        /// <param name="b"></param>
        /// <param name="id">
        /// identification string ('name') for the vector 
        /// field; the name of the components is constructed from this name
        /// </param>
        /// <param name="fac">
        /// factory for instantiation of <see cref="DGField"/>-objects; By the
        /// factory model, it becomes possible to use this container for
        /// different subclasses of <see cref="DGField"/>.
        /// </param>
        public VectorField(int D, Basis b, string id, Func<Basis, string, T> fac) {
            m_Components = new T[D];
            for (int d = 0; d < D; d++) {
                m_Components[d] = fac(b, id + "[" + d + "]");
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="D">
        /// number of components, i.d. dimension (see <see cref="Dim"/>)
        /// </param>
        /// <param name="b"></param>
        /// <param name="fac">
        /// factory for the instantiation of <see cref="DGField"/>-objects; By
        /// the factory model, it becomes possible to use this container for
        /// different subclasses of <see cref="DGField"/>.
        /// </param>
        public VectorField(int D, Basis b, Func<Basis, string, T> fac)
            : this(D, b, "", fac) {
        }

        /// <summary>
        /// ctor;
        /// </summary>
        /// <param name="components">components of the vector field</param>
        public VectorField(params T[] components) {
            var g = components[0].GridDat;
            for(int d = 1; d < components.Length; d++)
                if(!object.ReferenceEquals(g, components[d].GridDat))
                    throw new ArgumentException("all components must be defined on the same grid.");

            m_Components = (T[]) components.Clone();
            
        }

        /// <summary>
        /// dimension (number of vector components) of the vector field.
        /// </summary>
        public int Dim {
            get {
                return m_Components.Length;
            }
        }

        /// <summary>
        /// grid data object for this vector field.
        /// </summary>
        public IGridData GridDat {
            get {
                var g = this.m_Components[0].GridDat;
#if DEBUG
                for(int d = 1; d < this.m_Components.Length; d++)
                    Debug.Assert(object.ReferenceEquals(g, m_Components[d].GridDat));
#endif
                return g;
            }
        }

        /// <summary>
        /// Components
        /// </summary>
        T[] m_Components;

        /// <summary>
        /// returns the components of this vector field as some array of <see cref="DGField"/>-objects
        /// </summary>
        /// <returns></returns>
        public T[] ToArray() {
            return (T[])m_Components.Clone();
        }

        /// <summary>
        /// identifications of all components (<see cref="DGField.Identification"/>);
        /// </summary>
        public string[] Identifications {
            get {
                return m_Components.Select(x => x.Identification).ToArray();
            }
        }

        CoordinateMapping m_CoordMap = null; 

        /// <summary>
        /// returns a <see cref="CoordinateMapping"/> for all DG fields in this vector
        /// </summary>
        public CoordinateMapping Mapping {
            get {
                if (m_CoordMap == null) {
                    m_CoordMap = new CoordinateMapping(this.m_Components);
                }
                return m_CoordMap;
            }
        }

        CoordinateVector m_CoordVec = null;

        /// <summary>
        /// returns a <see cref="CoordinateVector"/> of the DG coordinates of
        /// the fields in this vector field
        /// </summary>
        public CoordinateVector CoordinateVector {
            get {
                if (m_CoordVec == null) {
                    m_CoordVec = new CoordinateVector(Mapping);
                }
                return m_CoordVec;
            }
        }

        /// <summary>
        /// access the <paramref name="d"/>-th component of this 
        /// vector filed
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        public T this[int d] {
            get {
                return m_Components[d];
            }
            set {
                throw new NotSupportedException("read-only list");
            }
        }

        /// <summary>
        /// adds an other vector field
        /// <paramref name="a"/>*<paramref name="mult"/> to this field;
        /// </summary>
        public void Acc(double mult, VectorField<T> a, CellMask m = null) {
            if (a.Dim != this.Dim)
                throw new ArgumentException("a is of another dimension", "a");

            for (int d = 0; d < m_Components.Length; d++)
                m_Components[d].Acc(mult, a.m_Components[d], m);
        }

        /// <summary>
        /// 'lax' accumulation (see
        /// <see cref="DGField.AccLaidBack(double, DGField, CellMask)"/>) of an
        /// other vector field
        /// <paramref name="a"/>*<paramref name="mult"/> to this field;
        /// </summary>
        public void AccLaidBack(double mult, VectorField<T> a, CellMask m = null) {
            if (a.Dim != this.Dim)
                throw new ArgumentException("a is of another dimension", "a");

            for (int d = 0; d < m_Components.Length; d++)
                m_Components[d].AccLaidBack(mult, a.m_Components[d], m);
        }

        /// <summary>
        /// multiplies this vector field by the number <paramref name="alpha"/>
        /// </summary>
        /// <param name="alpha">scaling factor</param>
        /// <param name="m">optional cell mask</param>
        public void Scale(double alpha, CellMask m = null) {
            foreach (DGField f in m_Components)
                f.Scale(alpha, m);
        }

        /// <summary>
        /// multiplies this object by <paramref name="f"/>*<paramref name="alpha"/>;
        /// Note that this is a relatively heavy operation.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        public void ScalePointwise(double alpha, DGField f) {
            for (int d = 0; d < m_Components.Length; d++) {
                // Workaround due to XDGFields
                DGField dummy = m_Components[d].CloneAs();
                m_Components[d].Clear();

                m_Components[d].ProjectProduct(alpha, dummy, f, null, true);
            }
        }

        /// <summary>
        /// adding two vector fields
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        static public VectorField<T> operator +(VectorField<T> a, VectorField<T> b) {
            VectorField<T> ret = (VectorField<T>)b.Clone();
            ret.Acc(1.0, a);
            return ret;
        }

        /// <summary>
        /// subtracting two vector fields
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        static public VectorField<T> operator -(VectorField<T> a, VectorField<T> b) {
            VectorField<T> ret = (VectorField<T>)b.Clone();
            ret.Scale(-1);
            ret.Acc(1.0, a);
            return ret;
        }


        /// <summary>
        /// returns <paramref name="a"/>*<paramref name="alpha"/>
        /// </summary>
        /// <param name="a">a vector field</param>
        /// <param name="alpha">scalar</param>
        /// <returns></returns>
        static public VectorField<T> operator *(VectorField<T> a, double alpha) {
            VectorField<T> ret = (VectorField<T>)a.Clone();
            ret.Scale(alpha);
            return ret;
        }

        /// <summary>
        /// returns <paramref name="a"/>*<paramref name="alpha"/>
        /// </summary>
        /// <param name="a">a vector field</param>
        /// <param name="alpha">scalar</param>
        /// <returns></returns>
        static public VectorField<T> operator *(double alpha, VectorField<T> a) {
            return a * alpha;
        }

        /// <summary>
        /// returns <paramref name="a"/>*<paramref name="alpha"/>
        /// </summary>
        /// <param name="a">a vector field</param>
        /// <param name="alpha">a scalar field</param>
        /// <returns></returns>
        static public VectorField<T> operator *(VectorField<T> a, DGField alpha) {
            T[] ret_comps = new T[a.Dim];
            for (int d = 0; d < ret_comps.Length; d++) {
                ret_comps[d] = (T)a.m_Components[d].Clone();
                ret_comps[d].Clear();
                ret_comps[d].ProjectProduct(1.0, a.m_Components[d], alpha);
            }

            return new VectorField<T>(ret_comps);
        }

        /// <summary>
        /// returns <paramref name="a"/>*<paramref name="alpha"/>
        /// </summary>
        /// <param name="a">a vector field</param>
        /// <param name="alpha">a scalar field</param>
        /// <returns></returns>
        static public VectorField<T> operator *(DGField alpha, VectorField<T> a) {
            return a * alpha;
        }

        /// <summary>
        /// non-shallow copy 
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            return this.CloneAs();
        }

        /// <summary>
        /// non-shallow copy 
        /// </summary>
        /// <returns></returns>
        public VectorField<T> CloneAs() {
            T[] c = new T[m_Components.Length];
            for (int d = 0; d < m_Components.Length; d++) {
                c[d] = (T)m_Components[d].Clone();
            }

            return new VectorField<T>(c);
        }

        /// <summary>
        /// copies the values of <paramref name="other"/> to this object
        /// (non-shallow) if possible
        /// </summary>
        /// <param name="other"></param>
        public void CopyFrom(VectorField<T> other) {
            if (this.Dim != other.Dim)
                throw new ApplicationException("unable to copy from other vector field - dimension mismatch.");

            for (int d = 0; d < m_Components.Length; d++) {
                this[d].CopyFrom(other[d]);
            }
        }

        /// <summary>
        /// accumulates the curl of a 3D DG vector field <paramref name="vec"/> 
        /// times <paramref name="alpha"/> to this vector field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* Curl(<paramref name="vec"/>)
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="vec"></param>
        /// <param name="em">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain  would be a waste of
        /// computational power. A proper execution mask for this case would be
        /// e.g. <see cref="BoSSS.Foundation.Grid.GridData.BoundaryCells"/>.)
        /// <br/>if null, the computation is carried out in the whole domain
        /// </param>
        /// <remarks>
        /// This method is based on
        /// <see cref="DGField.Derivative(double,DGField,int)"/>, i.e. it
        /// calculates derivatives by analytic cell-by-cell derivation of the
        /// DG polynomials;
        /// </remarks>
        /* <seealso cref="DGField.Curl2D{T}(double, VectorField{T}, CellMask)"/> */
        public void Curl3D(double alpha, VectorField<T> vec, CellMask em = null) {
            if (vec.Dim != 3)
                throw new ArgumentException("this method works only for 3-dimensional vector fields.", "vec");
            if (this.Dim != 3)
                throw new ApplicationException("this vector field must be 3-dimensional.");

            this.m_Components[0].Derivative(alpha, vec.m_Components[2], 1, em);
            this.m_Components[0].Derivative(-alpha, vec.m_Components[1], 2, em);

            this.m_Components[1].Derivative(alpha, vec.m_Components[0], 2, em);
            this.m_Components[1].Derivative(-alpha, vec.m_Components[2], 0, em);

            this.m_Components[2].Derivative(alpha, vec.m_Components[1], 0, em);
            this.m_Components[2].Derivative(-alpha, vec.m_Components[0], 1, em);
        }

        /// <summary>
        /// accumulates the curl of 3D DG vector field <paramref name="vec"/> 
        /// times <paramref name="alpha"/>
        /// to this vector field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* Curl(<paramref name="vec"/>)
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="vec"></param>
        /// <param name="optionalSubGrid">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="Grid.GridData.BoundaryCells"/>.)
        /// <br/> if null, the computation is carried out in the whole domain.
        /// </param>
        /// <param name="bndMode"></param>
        /// <remarks>
        /// This method is based on
        /// <see cref="DGField.DerivativeByFlux(double,DGField,int,SubGrid,SpatialOperator.SubGridBoundaryModes)"/>,
        /// i.e. it calculates derivatives by central-difference fluxes;
        /// </remarks>
        /* <seealso cref="DGField.Curl2DByFlux{T}"/> */
        public void Curl3DByFlux(double alpha, VectorField<T> vec,
            SubGrid optionalSubGrid = null, SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {
            if (vec.Dim != 3)
                throw new ArgumentException("this method works only for 3-dimensional vector fields.", "vec");
            if (this.Dim != 3)
                throw new ApplicationException("this vector field must be 3-dimensional.");

            this.m_Components[0].DerivativeByFlux(alpha, vec.m_Components[2], 1, optionalSubGrid, bndMode);
            this.m_Components[0].DerivativeByFlux(-alpha, vec.m_Components[1], 2, optionalSubGrid, bndMode);

            this.m_Components[1].DerivativeByFlux(alpha, vec.m_Components[0], 2, optionalSubGrid, bndMode);
            this.m_Components[1].DerivativeByFlux(-alpha, vec.m_Components[2], 0, optionalSubGrid, bndMode);

            this.m_Components[2].DerivativeByFlux(alpha, vec.m_Components[1], 0, optionalSubGrid, bndMode);
            this.m_Components[2].DerivativeByFlux(-alpha, vec.m_Components[0], 1, optionalSubGrid, bndMode);
        }

        /// <summary>
        /// this = this + <paramref name="alpha"/>*Gradient of
        /// <paramref name="f"/>
        /// </summary>
        /// <param name="alpha">scaling</param>
        /// <param name="f"></param>
        /// <param name="em">
        /// optional mask, null is the full domain
        /// </param>
        /// <remarks>
        /// this function is based on the
        /// <see cref="DGField.Derivative(double,DGField,int,CellMask)"/>-method
        /// </remarks>
        public void Gradient(double alpha, DGField f, CellMask em = null) {
            int D = f.Basis.GridDat.SpatialDimension;

            if (this.Dim != D)
                throw new NotSupportedException("this field has not the right spatial dimension");

            for (int d = 0; d < D; d++) {
                this.m_Components[d].Derivative(alpha, f, d, em);
            }
        }


        /// <summary>
        /// Accumulates the derivative of <paramref name="f"/> with respect to
        /// x_i to the i-th component of this vector field by making use of
        /// <see cref="DGField.DerivativeByFlux"/>
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="f"></param>
        /// <param name="optionalSubGrid">
        /// An optional restriction to the domain in which the derivative is
        /// computed (it may, e.g. be only required in boundary cells, so a
        /// computation over the whole domain would be a waste of computational
        /// power. A proper execution mask would be see e.g. 
        /// <see cref="Grid.GridData.BoundaryCells"/>.)
        /// <br/>
        /// if null, the computation is carried out in the whole domain.
        /// </param>
        /// <param name="bndMode">
        /// </param>
        public void GradientByFlux(
            double alpha,
            DGField f,
            SubGrid optionalSubGrid = null,
            SpatialOperator.SubGridBoundaryModes bndMode = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {

            int D = f.Basis.GridDat.SpatialDimension;

            if (this.Dim != D)
                throw new NotSupportedException("this field has not the right spatial dimension");

            for (int d = 0; d < D; d++) {
                this.m_Components[d].DerivativeByFlux(alpha, f, d, optionalSubGrid, bndMode);
            }
        }

        #region IEnumerable<T> Members

        /// <summary>
        /// as defined by interface;
        /// </summary>
        public IEnumerator<T> GetEnumerator() {
            return ((IEnumerable<T>)m_Components).GetEnumerator();
        }

        #endregion

        #region IEnumerable Members

        /// <summary>
        /// as defined by interface;
        /// </summary>
        IEnumerator IEnumerable.GetEnumerator() {
            return m_Components.GetEnumerator();
        }

        #endregion

        /// <summary>
        /// The inner product of vector fields <paramref name="a"/> and <paramref name="b"/>;
        /// </summary>
        static public double InnerProduct(VectorField<T> a, VectorField<T> b) {
            if (a.Dim != b.Dim)
                throw new ArgumentException("mismatch in dimension");

            double acc = 0;
            for (int d = 0; d < a.Dim; d++) {
                acc += DGField.InnerProduct(a[d], b[d]);
            }
            return acc;
        }

        #region IList<T> Members

        /// <summary>
        /// index of <paramref name="item"/>;
        /// </summary>
        public int IndexOf(T item) {
            return Array.IndexOf(this.m_Components, item);
        }

        /// <summary>
        /// not supported (read - only);
        /// </summary>
        public void Insert(int index, T item) {
            throw new NotSupportedException("read-only list.");
        }

        /// <summary>
        /// not supported (read - only);
        /// </summary>
        public void RemoveAt(int index) {
            throw new NotSupportedException("read-only list.");
        }

        #endregion

        #region ICollection<T> Members

        /// <summary>
        /// not supported (read - only);
        /// </summary>
        public void Add(T item) {
            throw new NotSupportedException("read-only collection.");
        }

        /// <summary>
        /// %
        /// </summary>
        public bool Contains(T item) {
            return (IndexOf(item) >= 0);
        }

        /// <summary>
        /// %
        /// </summary>
        public void CopyTo(T[] array, int arrayIndex) {
            for (int i = 0; i < this.Dim; i++)
                array[arrayIndex + i] = this.m_Components[i];
        }

        /// <summary>
        /// alias for <see cref="Dim"/>;
        /// </summary>
        public int Count {
            get {
                return this.Dim;
            }
        }

        /// <summary>
        /// always true
        /// </summary>
        public bool IsReadOnly {
            get {
                return true;
            }
        }

        /// <summary>
        /// not supported (read - only);
        /// </summary>
        public bool Remove(T item) {
            throw new NotSupportedException("read-only collection.");
        }

        #endregion

        /// <summary>
        /// 2-Norm of the vector field.
        /// </summary>
        /// <param name="CM">
        /// optional restriction of the computational domain
        /// </param>
        /// <returns></returns>
        public double L2Norm(CellMask CM = null) {
            double acc = 0;

            for (int d = 0; d < Dim; d++) {
                double nd = m_Components[d].L2Norm(CM);
                acc += nd * nd;
            }

            return Math.Sqrt(acc);
        }

    }
}
