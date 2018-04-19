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
using System.Text;
using ilPSP;

namespace BoSSS.Foundation.Grid {
    
    
    
    /// <summary>
    /// A list of polynomials. The purpose of this class is to optimize the 
    /// evaluation of polynomials and their derivatives. Since this is usually done
    /// for multiple polynomials at once, optimizations can be made.
    /// </summary>
    public class PolynomialList : IList<Polynomial> {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="polys"></param>
        public PolynomialList(IEnumerable<Polynomial> polys) {
            this.m_polys = polys.ToArray();
            if(polys.Count() < 1)
                throw new ArgumentException("Empty lists are not supported.");

            this.SpatialDimension = m_polys[0].SpatialDimension;
            this.MaxAbsoluteDegree = -1;
            for(int i = 0; i < m_polys.Length; i++) {
                if(m_polys[i].SpatialDimension != this.SpatialDimension) {
                    throw new ArgumentException("All polynomials must have the same spatial dimension.");
                }

                this.MaxAbsoluteDegree = Math.Max(this.MaxAbsoluteDegree, m_polys[i].AbsoluteDegree);
            }

            this.Values = new Caching.CacheLogicImpl_Ns(this.GetValues);
        }

        Polynomial[] m_polys;

        /// <summary>
        /// Spatial dimension in which all polynomials are defined.
        /// </summary>
        public int SpatialDimension {
            get;
            private set;
        }

        /// <summary>
        /// The maximum <see cref="Polynomial.AbsoluteDegree"/> over all polynomials.
        /// </summary>
        public int MaxAbsoluteDegree {
            get;
            private set;
        }

#region ILIST_STUFF
        
        /// <summary>
        /// As defined by interface <see cref="IList{T}"/>.
        /// </summary>
        public int IndexOf(Polynomial item) {
            return  Array.IndexOf(this.m_polys, item);
        }

        /// <summary>
        /// Not supported - this is read-only.
        /// </summary>
        public void Insert(int index, Polynomial item) {
            throw new NotSupportedException("This list is read-only.");
        }

        /// <summary>
        /// Not supported - this is read-only.
        /// </summary>
        public void RemoveAt(int index) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Returns the <paramref name="index"/>-th polynomial.
        /// Setting is not supported, since this list is read-only.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public Polynomial this[int index] {
            get {
                return m_polys[index];
            }
            set {
                throw new NotSupportedException("This list is read-only.");
            }
        }

        /// <summary>
        /// Not supported - this is read-only.
        /// </summary>
        public void Add(Polynomial item) {
            throw new NotSupportedException("This list is read-only.");
        }

        /// <summary>
        /// Not supported - this is read-only.
        /// </summary>
        public void Clear() {
            throw new NotSupportedException("This list is read-only.");
        }


        /// <summary>
        /// As defined by interface <see cref="IList{T}"/>.
        /// </summary>
        public bool Contains(Polynomial item) {
            return this.IndexOf(item) >= 0;
        }

        /// <summary>
        /// As defined by interface <see cref="IList{T}"/>.
        /// </summary>
        public void CopyTo(Polynomial[] array, int arrayIndex) {
            Array.Copy(this.m_polys, 0, array, arrayIndex, this.m_polys.Length);
        }

        /// <summary>
        /// Number of polynomials in the list.
        /// </summary>
        public int Count {
            get {
                return this.m_polys.Length;
            }
        }

        /// <summary>
        /// yes.
        /// </summary>
        public bool IsReadOnly {
            get {
                return true;
            }
        }

        /// <summary>
        /// Not supported - this is read-only.
        /// </summary>
        public bool Remove(Polynomial item) {
            throw new NotSupportedException("This list is read-only.");
        }

        /// <summary>
        /// As defined by interface <see cref="IList{T}"/>.
        /// </summary>
        public IEnumerator<Polynomial> GetEnumerator() {
            return this.m_polys.AsEnumerable().GetEnumerator();
        }

        /// <summary>
        /// As defined by interface <see cref="IList{T}"/>.
        /// </summary>
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return this.m_polys.GetEnumerator();
        }

#endregion ILIST_STUFF


        /// <summary>
        /// returns the number of basis polynomials
        /// which are of degree <paramref name="deg"/>;
        /// (see also <see cref="FirstPolynomialIndexforDegree"/>);
        /// </summary>
        /// <param name="deg"></param>
        /// <returns></returns>
        public int NoOfPolynomialsPerDegree(int deg) {
            if(deg < 0 || deg > this.MaxAbsoluteDegree)
                throw new ArgumentException("out of range", "deg");

            int n = 0;
            foreach(Polynomial p in this.m_polys) {
                if(p.AbsoluteDegree == deg)
                    n++;
            }
            return n;
        }

        /// <summary>
        /// returns the index of the first polynomial (see <see cref="Polynomials"/>)
        /// which is of degree <paramref name="deg"/>;
        /// (see also <see cref="NoOfPolynomialsPerDegree"/>);
        /// </summary>
        /// <param name="deg"></param>
        /// <returns></returns>
        public int FirstPolynomialIndexforDegree(int deg) {
            if(deg < 0 || deg > this.MaxAbsoluteDegree)
                throw new ArgumentException("out of range", "deg");

            for(int i = 0; i < this.m_polys.Length; i++) {
                if(this.m_polys[i].AbsoluteDegree == deg)
                    return i;
            }
            throw new ApplicationException("fatal error: should never happen;");
        }

        /// <summary>
        /// Non-cached evaluation.
        /// </summary>
        /// <param name="NS">Node set.</param>
        /// <param name="R">
        /// Output
        ///  - 1st index: node 
        ///  - 2nd index: polynomial
        /// </param>
        public void Evaluate(NodeSet NS, MultidimensionalArray R) {
            int D = NS.SpatialDimension;
            if(D != this.SpatialDimension)
                throw new ArgumentException("Wrong spatial dimension of node set.");
            int K = NS.NoOfNodes;
            int N = this.Count;

            if(R.Dimension != 2)
                throw new ArgumentException("Wrong dimension of output array.");
            if(R.GetLength(0) != K)
                throw new ArgumentException("Wrong 1st length of output array.");
            if(R.GetLength(1) != N)
                throw new ArgumentException("Wrong 2nd length of output array.");

            // use at least cached monomials: this ensures that the monomials are 
            // available for the highest degree that will be needed
            Caching.MonomialCache.Instance.GetMonomials(NS, this.MaxAbsoluteDegree); 

            for(int n = 0; n < N; n++) {
                this[n].Evaluate(R.ExtractSubArrayShallow(-1, n), NS);
            }
        }

        MultidimensionalArray GetValues(NodeSet NS) {
            MultidimensionalArray R = MultidimensionalArray.Create(NS.NoOfNodes, this.Count);
            this.Evaluate(NS, R);
            return R;
        }

        /// <summary>
        /// Cached evaluation.
        /// </summary>
        public Caching.CacheLogicImpl_Ns Values {
            get;
            private set;
        }

    }
     
}
