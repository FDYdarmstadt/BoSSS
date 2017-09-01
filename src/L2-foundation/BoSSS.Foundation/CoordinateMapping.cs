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

namespace BoSSS.Foundation {
    


    /// <summary>
    /// This class provides bijective mappings:
    /// <list type="bullet">
    ///   <item>
    ///    First, between <em>local field coordinate indices</em> and <em>local unique indices</em>,
    ///    by the methods <see cref="LocalUniqueCoordinateIndex"/> and <see cref="LocalFieldCoordinateIndex"/>, 
    ///    and ...
    ///   </item>
    ///   <item>
    ///    second, between <em>local unique indices</em> and <em>global unique indices</em>,
    ///    by the methods <see cref="UnsetteledCoordinateMapping.Global2LocalIndex"/> and <see cref="UnsetteledCoordinateMapping.Local2GlobalIndex"/>.
    ///   </item>
    /// </list>
    /// In easy words, that means it mapps the DG coordinate from a list of fields into one long, one-dimensional vector
    /// and it converts between indices that are valid only on the actual MPI process (local indices)
    /// and indices that are valid amoung all MPI processes in the current MPI communicator (global indices).<br/>
    /// In contrast to <see cref="UnsetteledCoordinateMapping"/>, this is
    /// a coordinate mapping that is bound to a specific list of fields
    /// (instead a list of <see cref="Basis"/>-objects);
    /// </summary>
    /// <remarks>
    /// A "<em>local index</em>" is an index that is valid only on the actual MPI process,
    /// while a "<em>gloabl index</em>" is an index that is valid among all MPI processes 
    /// (in the corresponding MPI communicator).<br/>
    /// A "<em>local field coordinate index</em>" is a tuple (<i>f</i>,<i>j</i>,<i>n</i>), where...<br/>
    /// <list type="bullet">
    ///   <item><i>f</i> is a <see cref="DGField"/>,</item>
    ///   <item><i>j</i> is a local cell index and </item>
    ///   <item><i>n</i> is a polynomial index;</item>
    /// </list>
    /// A "<em>local unique index</em>" is a single integer that is unique for every field in this mapping, 
    /// for every local cell and for every polynomial;
    /// </remarks>
    public class CoordinateMapping : UnsetteledCoordinateMapping, IList<DGField> {

        /// <summary>
        /// clears all DG fields in the mapping
        /// </summary>
        public void Clear() {
            foreach (DGField f in this.m_Fields)
                f.Clear();
        }



        static Basis[] GetBasisList(IList<DGField> _fields) {
            Basis[] ret = new Basis[_fields.Count];
            for (int i = 0; i < _fields.Count; i++)
                ret[i] = _fields[i].Basis;
            return ret;
        }


        /// <summary>
        /// Constructs a new mapping from an ordered list of fields;
        /// </summary>
        /// <param name="_fields">the list of <see cref="DGField"/>'s for this mapping.</param>
        public CoordinateMapping(params DGField[] _fields)
            : this((IList<DGField>) _fields) { }


        /// <summary>
        /// Constructs a new mapping from an ordered list of fields;
        /// </summary>
        /// <param name="_fields">the list of <see cref="DGField"/>'s for this mapping.</param>
        public CoordinateMapping(IList<DGField> _fields) :
            base(GetBasisList(_fields)) {
            m_Fields = new DGField[_fields.Count];
            for (int i = 0; i < m_Fields.Length; i++)
                m_Fields[i] = _fields[i];
        }


        /// <summary>
        /// the fields in this mapping;
        /// </summary>
        internal DGField[] m_Fields;

        //IReadOnlyList<DGField> m_RdonlyFields;


        /// <summary>
        /// the fields in this mapping; 
        /// </summary>
        /// <remarks>
        /// By this list, an <i>field index</i> for each field in this mapping is defined;
        /// </remarks>
        public IList<DGField> Fields {
            get {
                //if(m_RdonlyFields == null) {
                //    m_RdonlyFields = m_Fields.ToList().AsReadOnly();
                //}
                //return m_RdonlyFields;
                return m_Fields;
            }
        }


        /// <summary>
        /// computes a local unique coordinate index ("local" means local on this processor);
        /// this index is unique over all fields (in this mapping), over all cells, over all basis functions, 
        /// but it's only locally (on this processor) valid.
        /// A local index in the update range (smaller than <see cref="UnsetteledCoordinateMapping.LocalLength"/>) can be converted into 
        /// a global index by adding <see cref="ilPSP.Partitioning.i0"/>.
        /// </summary>
        /// <param name="f">the field</param>
        /// <param name="j">local cell index</param>
        /// <param name="n">basis index.</param>
        /// <returns>
        /// the local coordinate index for the given parameters.
        /// </returns>
        /// <remarks>
        /// This method is not supported and will throw an exception 
        /// if <see cref="Fields"/>==null.
        /// </remarks>
        public int LocalUniqueCoordinateIndex(DGField f, int j, int n) {
            if (m_Fields == null)
                throw new ApplicationException("calling this method is not supported, because the Fields - property is null; Use other constructor;");
            int find = Array.IndexOf<DGField>(m_Fields, f);
            return LocalUniqueCoordinateIndex(find, j, n);
        }

        /*
        /// <summary>
        /// computes a global unique coordinate index ("global" means over all MPI processors);
        /// this index is unique over all fields (in this mapping), over all cells, over all basis functions;
        /// </summary>
        /// <param name="f">the field</param>
        /// <param name="j">local cell index</param>
        /// <param name="n">basis index.</param>
        /// <returns>
        /// the global DG coordinate index for the given parameters.
        /// </returns>
        public long GlobalUniqueCoordinateIndex(DGField f, int j, int n) {
            int iloc = LocalUniqueCoordinateIndex(f, j, n);
            return Local2GlobalIndex(iloc);
        }
        */

        /// <summary>
        /// inverse mapping of <see cref="LocalUniqueCoordinateIndex"/>;
        /// </summary>
        /// <param name="Index">local unique coordinate index</param>
        /// <param name="f">on exit, the field which <paramref name="Index"/> belongs to.</param>
        /// <param name="j">on exit, the local cell index</param>
        /// <param name="n">on exit, the basis function index</param>
        /// <returns></returns>
        public void LocalFieldCoordinateIndex(int Index, out DGField f, out int j, out int n) {
            int ifld;
            LocalFieldCoordinateIndex(Index, out ifld, out j, out n);
            f = m_Fields[ifld];
        }


        /// <summary>
        /// Two <see cref="CoordinateMapping"/> are equal if all entries of their <see cref="Fields"/>-properties are
        /// equal in reference. Note that this is a very strict condition; if the mappings should be compared just on
        /// <see cref="UnsetteledCoordinateMapping"/>-level, the <see cref="UnsetteledCoordinateMapping.EqualsUnsetteled"/>-comparison
        /// should be used.
        /// </summary>
        public override bool Equals(object obj) {
            if (!base.Equals(obj))
                return false;
            
            CoordinateMapping othr = obj as CoordinateMapping;
            if (othr == null)
                return false;


            if (othr.m_Fields.Length != this.m_Fields.Length)
                return false;
            for (int i = 0; i < othr.m_Fields.Length; i++)
                if (!object.ReferenceEquals(this.m_Fields[i], othr.m_Fields[i]))
                    return false;

            return true;
        }
        
        /// <summary>
        /// calls base implementation
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }
        
        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>.
        /// </summary>
        public int IndexOf(DGField item) {
            return Array.IndexOf(this.m_Fields,item);
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>, but not supported
        /// </summary>
        public void Insert(int index, DGField item) {
            throw new NotSupportedException("items are fixed");
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>, but not supported
        /// </summary>
        public void RemoveAt(int index) {
            throw new NotSupportedException("items are fixed");
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>, setting is not supported
        /// </summary>
        public DGField this[int index] {
            get {
                return m_Fields[index];
            }
            set {
                throw new NotSupportedException("items are fixed");
            }
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>, but not supported
        /// </summary>
        public void Add(DGField item) {
            throw new NotSupportedException("items are fixed");
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>.
        /// </summary>
        public bool Contains(DGField item) {
            return this.Fields.Contains(item);
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>.
        /// </summary>
        public void CopyTo(DGField[] array, int arrayIndex) {
            this.Fields.CopyTo(array, arrayIndex);
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>.
        /// </summary>
        public int Count {
            get { return m_Fields.Length; }
        }

        /// <summary>
        /// true!
        /// </summary>
        public bool IsReadOnly {
            get { return true; }
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>, but not supported
        /// </summary>
        public bool Remove(DGField item) {
            throw new NotSupportedException("items are fixed.");
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>.
        /// </summary>
        public IEnumerator<DGField> GetEnumerator() {
            return this.Fields.GetEnumerator();
        }

        /// <summary>
        /// as defined by interface <see cref="IList{Field}"/>.
        /// </summary>
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return this.Fields.GetEnumerator();
        }
    }
}
