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
using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// Represents as set of vertices (i.e., each vertex is unique) and
    /// associates them with unique ids, e.g. in order to allow for an
    /// efficient checking for duplicates. Additionally, sets can be nested
    /// which enables the usage of "layers of vertices". This can be useful if
    /// vertices should be kept in distinct sets while, at the same time, a
    /// vertex must not be contained in both sets (i.e. a set and its parent
    /// set are distinct). 
    /// </summary>
    /// <remarks>
    /// Throughout this class, the word "local" refers to properties associated
    /// with this set alone (i.e., not considering potential ancestors) while
    /// "global" refers to properties with this set <b>and</b> all of its
    /// ancestors.
    /// </remarks>
    public class NestedVertexSet {

        /// <summary>
        /// The minimum distance below which two vertices are considered equal.
        /// </summary>
        private const double EPS = 1e-14;

        /// <summary>
        /// The vertex storage. The implementation ensures that this is
        /// actually a set (and not a list) of vertices.
        /// </summary>
        /// <remarks>
        /// We cannot use one of the built-in set implementations since this
        /// would either be rather inefficient or, if using an ordered set,
        /// require a complete ordering of the vertices which is not possible
        /// (without loosing the ability to look for vertices that are the
        /// same up to <see cref="EPS"/>).
        /// </remarks>
        private List<Vertex> listOfVertices = new List<Vertex>();

        /// <summary>
        /// A mapping between the norm of a vertex (see
        /// <see cref="Vertex.Norm"/>) and the indices to all vertices in
        /// <see cref="listOfVertices"/> that have the same norm.
        /// </summary>
        private OrderedMultiDictionary<double, int> normToVertexIndexMap;

        /// <summary>
        /// Cache for <see cref="LocalVertices"/>.
        /// </summary>
        private double[] localVerticesAsArray;

        /// <summary>
        /// Constructs an empty set without any ancestors.
        /// </summary>
        /// <param name="spatialDimension">
        /// The spatial dimension of the vertices.
        /// </param>
        public NestedVertexSet(int spatialDimension) {
            this.SpatialDimension = spatialDimension;
            normToVertexIndexMap = new OrderedMultiDictionary<double, int>();
        }

        /// <summary>
        /// Constructs an "empty" set based on <paramref name="parrentSet"/>.
        /// </summary>
        /// <param name="parrentSet">
        /// A parrent containing vertices that cannot be added to this set
        /// because they're considered duplicates.
        /// </param>
        public NestedVertexSet(NestedVertexSet parrentSet)
            : this(parrentSet.SpatialDimension) {
            this.Parrent = parrentSet;
        }

        /// <summary>
        /// The number of vertices contained in this set (i.e., not counting
        /// the vertices in <see cref="Parrent"/>).
        /// </summary>
        public int NumberOfLocalVertices {
            get {
                return listOfVertices.Count;
            }
        }

        /// <summary>
        /// The spatial dimension of the vertices contained in this set.
        /// </summary>
        public int SpatialDimension {
            get;
            private set;
        }

        /// <summary>
        /// An (optional) parrent set. Any vertices in <see cref="Parrent"/>
        /// cannot be added to this set because they're considered duplicates.
        /// </summary>
        public NestedVertexSet Parrent {
            get;
            private set;
        }

        /// <summary>
        /// Returns an array of all vertices contained in this set (but not the
        /// vertices contained in <see cref="Parrent"/>). This operation is
        /// cached (that is, if no vertices have been added since the last
        /// call).
        /// </summary>
        public MultidimensionalArray LocalVertices {
            get {
                if (NumberOfLocalVertices == 0) {
                    return null;
                } else {
                    if (localVerticesAsArray == null
                        || localVerticesAsArray.GetLength(0) != NumberOfLocalVertices) {
                        localVerticesAsArray = new double[NumberOfLocalVertices * SpatialDimension];
                        for (int i = 0; i < NumberOfLocalVertices; i++) {
                            double[] vertex = listOfVertices[i].Coordinates;
                            for (int j = 0; j < SpatialDimension; j++) {
                                localVerticesAsArray[i * SpatialDimension + j] = vertex[j];
                            }
                        }
                    }
                    return MultidimensionalArray.CreateWrapper(localVerticesAsArray, NumberOfLocalVertices, SpatialDimension);
                }
            }
        }

        /// <summary>
        /// Registers a new vertex in the set if it is not already present in
        /// this set or one of its ancestors.
        /// </summary>
        /// <param name="vertex">
        /// An array containing the spatial coordinates of the vertex.
        /// </param>
        /// <returns>
        /// The unique index of the vertex <paramref name="vertex"/> (i.e., the
        /// index is unique within this set and its ancestors). This might be a
        /// new index (if the vertex was previously unknown) or an existing
        /// index (if the vertex is already present in this set or one of its
        /// ancestors).
        /// </returns>
        public int RegisterVertex(double[] vertex) {
            Vertex boxedVertex = new Vertex(vertex);
            int index = GetGlobalVertexIndex(ref boxedVertex);
            if (index < 0) {
                listOfVertices.Add(boxedVertex);

                int localIndex = listOfVertices.Count - 1;
                normToVertexIndexMap.Add(boxedVertex.Norm, localIndex);

                return IndexOffset + localIndex;
            } else {
                return index;
            }
        }

        /// <summary>
        /// Checks whether the given index <paramref name="globalVertexIndex"/>
        /// is associated with a vertex in this set (i.e., ancestors are
        /// <b>not</b> considered)
        /// </summary>
        /// <param name="globalVertexIndex">
        /// The index in question.
        /// </param>
        /// <returns>
        /// True, if the vertex identified by
        /// <paramref name="globalVertexIndex"/> is contained in this set.
        /// Otherwise, false is returned.
        /// </returns>
        public bool VertexIsContainedInLocalSet(int globalVertexIndex) {
            int localIndex = GetLocalFromGlobalVertexIndex(globalVertexIndex);
            return (localIndex >= 0) && (localIndex < listOfVertices.Count);
        }

        /// <summary>
        /// Transforms the given index <paramref name="globalVertexIndex"/>
        /// into a local vertex index (i.e., an index in the range
        /// [0, <see cref="NumberOfLocalVertices"/>]
        /// </summary>
        /// <param name="globalVertexIndex"></param>
        /// <returns>
        /// The local index in the range [0, <see cref="NumberOfLocalVertices"/>]
        /// associated with <paramref name="globalVertexIndex"/>.
        /// </returns>
        public int GetLocalFromGlobalVertexIndex(int globalVertexIndex) {
            return globalVertexIndex - IndexOffset;
        }

        /// <summary>
        /// Retrieves a vertex by its global id.
        /// </summary>
        /// <param name="globalVertexIndex">
        /// The global id of the vertex to be retrieved.
        /// </param>
        /// <returns>
        /// The spatial coordinates of the vertex with global id
        /// <paramref name="globalVertexIndex"/>.
        /// </returns>
        public double[] GetVertex(int globalVertexIndex) {
            int localIndex = GetLocalFromGlobalVertexIndex(globalVertexIndex);
            if (localIndex < 0) {
                if (Parrent == null) {
                    throw new Exception("Unknown vertex index " + globalVertexIndex);
                } else {
                    return Parrent.GetVertex(globalVertexIndex);
                }
            } else {
                return listOfVertices[localIndex].Coordinates;
            }
        }

        /// <summary>
        /// Tries to determine the global vertex index associated with the
        /// given vertex.
        /// </summary>
        /// <param name="vertex">
        /// The vertex in question.
        /// </param>
        /// <returns>
        /// The global vertex id of <paramref name="vertex"/> if it is
        /// contained in this set or one of its ancestors. Otherwise, -1 will
        /// be returned.
        /// </returns>
        private int GetGlobalVertexIndex(ref Vertex vertex) {
            Debug.Assert(
                vertex.Coordinates.Length == SpatialDimension,
                "Wrong spatial dimension of vertex");

            int index = GetLocalVertexIndex(ref vertex);
            if (index < 0) {
                if (Parrent != null) {
                    index = Parrent.GetGlobalVertexIndex(ref vertex);
                }
                return index;
            } else {
                return index + IndexOffset;
            }
        }

        /// <summary>
        /// Tries to determine the local vertex index associated with the given
        /// vertex.
        /// </summary>
        /// <param name="vertex">
        /// The vertex in question.
        /// </param>
        /// <returns>
        /// The local vertex id of <paramref name="vertex"/> if it is
        /// contained in this set. Otherwise, -1 will be returned.
        /// </returns>
        private int GetLocalVertexIndex(ref Vertex vertex) {
            int index = -1;

            foreach (int candidateIndex in normToVertexIndexMap.ValuesInRange(vertex.Norm - EPS, vertex.Norm + EPS)) {
                Vertex candidateVertex = listOfVertices[candidateIndex];

                bool isEqual = true;
                for (int i = 0; i < vertex.Coordinates.Length; i++) {
                    isEqual &= Math.Abs(vertex.Coordinates[i] - candidateVertex.Coordinates[i]) <= EPS;
                }

                if (isEqual) {
                    index = candidateIndex;
                }
            }
            return index;
        }

        /// <summary>
        /// The offset when transforming between local und global coordinates.
        /// Obviously, this is equal to the number of vertices contained in all
        /// ancestors of this set.
        /// </summary>
        private int IndexOffset {
            get {
                int offset = 0;
                if (Parrent != null) {
                    offset += Parrent.listOfVertices.Count;
                    offset += Parrent.IndexOffset;
                }
                return offset;
            }
        }

        /// <summary>
        /// A vertex (i.e., a point). Mainly exists because it caches the norm
        /// of the vertex.
        /// </summary>
        private struct Vertex {

            /// <summary>
            /// <see cref="Coordinates"/>.
            /// </summary>
            private double[] coordinates;

            /// <summary>
            /// <see cref="Norm"/>.
            /// </summary>
            private double norm;

            /// <summary>
            /// Constructs a new vertex with spatial coordinates
            /// <paramref name="coordinates"/> and calculates the
            /// <see cref="Norm"/> of the vertex.
            /// </summary>
            /// <param name="coordinates"></param>
            public Vertex(double[] coordinates) {
                this.coordinates = coordinates;
                norm = 0.0;
                for (int i = 0; i < coordinates.Length; i++) {
                    norm += coordinates[i] * coordinates[i];
                }
            }

            /// <summary>
            /// The spatial coordinates of the vertex.
            /// </summary>
            public double[] Coordinates {
                get {
                    return coordinates;
                }
            }

            /// <summary>
            /// Some norm of the coordinates of this vertex.
            /// </summary>
            /// <remarks>
            /// Currently this is the squared L2-norm because we can calculate
            /// it efficiently and do not care about the real norm.
            /// </remarks>
            public double Norm {
                get {
                    return norm;
                }
            }
        }

        /// <summary>
        /// Compares two vertices based on their norm. In order to do so, it
        /// consumes their indices into <see cref="listOfVertices"/> and uses
        /// the cached norm to impose the ordering. Is used for the ordering of
        /// <see cref="normToVertexIndexMap"/>.
        /// </summary>
        private class VertexNormComparer : IComparer<int> {

            /// <summary>
            /// A list containing the vertices to be compared.
            /// </summary>
            private List<Vertex> listOfVertices;

            /// <summary>
            /// Constructs a comparer for vertices in
            /// <paramref name="listOfVertices"/>.
            /// </summary>
            /// <param name="listOfVertices">
            /// Vertices to be compared.
            /// </param>
            public VertexNormComparer(List<Vertex> listOfVertices) {
                this.listOfVertices = listOfVertices;
            }

            #region IComparer<int> Members

            /// <summary>
            /// Compares the vertices identified by their local indices (i.e.,
            /// the indicies into <see cref="listOfVertices"/>).
            /// </summary>
            /// <param name="x">First operand</param>
            /// <param name="y">Second operand</param>
            /// <returns>
            /// <see cref="IComparer{Double}.Compare"/>
            /// </returns>
            public int Compare(int x, int y) {
                return listOfVertices[x].Norm.CompareTo(listOfVertices[y].Norm);
            }

            #endregion
        }

        /// <summary>
        /// Variant of <see cref="SortedList"/> (which, in fact, is a
        /// list-backed sorted dictionary rather than a list). That allows
        /// multiple values per key and provides an efficient method for
        /// iterating over the values in a range of keys (see 
        /// <see cref="ValuesInRange"/>).
        /// </summary>
        /// <typeparam name="TKey">
        /// The type of the keys of the dictionary.
        /// </typeparam>
        /// <typeparam name="TValue">
        /// The type of the value of a single entry.
        /// </typeparam>
        public class OrderedMultiDictionary<TKey, TValue>
            : IDictionary<TKey, ICollection<TValue>> where TKey : IComparable {

            /// <summary>
            /// The main part of the class which actually holds the data.
            /// </summary>
            private SortedList<TKey, ICollection<TValue>> storage =
                new SortedList<TKey, ICollection<TValue>>();

            /// <summary>
            /// Varient of <see cref="Add(TKey, ICollection{TValue})"/> which
            /// allows adding a single value (instead of a list of values).
            /// </summary>
            /// <param name="key">
            /// The descired dictionary key.
            /// </param>
            /// <param name="value">
            /// The value to be stored.
            /// </param>
            public void Add(TKey key, TValue value) {
                if (storage.ContainsKey(key)) {
                    storage[key].Add(value);
                } else {
                    storage[key] = new List<TValue>();
                    storage[key].Add(value);
                }
            }

            /// <summary>
            /// Provides an efficient iterator for all values whose keys are in
            /// the range [<paramref name="fromKey"/>,<paramref name="toKey"/>]
            /// (borders excluded).
            /// </summary>
            /// <param name="fromKey">
            /// The start of the range.
            /// </param>
            /// <param name="toKey">
            /// The end of the range.
            /// </param>
            /// <returns>
            /// In iterator containing all values from (key, value) for which 
            /// <paramref name="fromKey"/> &lt; key &lt; <paramref name="toKey"/>
            /// holds.
            /// </returns>
            public IEnumerable<TValue> ValuesInRange(TKey fromKey, TKey toKey) {
                IList<TKey> keys = storage.Keys;

                // Binaray search for start key
                int firstIndex = 0;
                int lastIndex = keys.Count - 1;

                int currentIndex = -1;
                while (firstIndex <= lastIndex) {
                    int midIndex = (firstIndex + lastIndex) / 2;
                    int comparisonResult = keys[midIndex].CompareTo(fromKey);

                    if (comparisonResult < 0) {
                        firstIndex = midIndex + 1;
                    } else if (comparisonResult > 0) {
                        lastIndex = midIndex - 1;

                        // We are looking for the first key where key > fromKey. We
                        // already such an element fulfilling this condition
                        // -> store it even though we might find something better
                        currentIndex = midIndex;
                    } else {
                        currentIndex = midIndex;
                        break;
                    }
                }

                if (currentIndex < 0 || currentIndex >= keys.Count) {
                    yield break;
                }

                // yield until key > toKey
                while (currentIndex < keys.Count && keys[currentIndex].CompareTo(toKey) < 0) {
                    foreach (TValue value in storage[keys[currentIndex]]) {
                        yield return value;
                    }
                    currentIndex++;
                }
            }

            #region IDictionary<TKey,ICollection<TValue>> Members

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.Add"/>
            /// </summary>
            /// <param name="key">
            /// <see cref="IDictionary{TKey,TValueCollection}.Add"/>
            /// </param>
            /// <param name="values">
            /// <see cref="IDictionary{TKey,TValueCollection}.Add"/>
            /// </param>
            /// <remarks>
            /// Leaves are currently stored as standard
            /// <see cref="List{TValue}"/>s. This may change in the future.
            /// </remarks>
            public void Add(TKey key, ICollection<TValue> values) {
                if (storage.ContainsKey(key)) {
                    foreach (TValue value in values) {
                        storage[key].Add(value);
                    }
                } else {
                    storage[key] = new List<TValue>(values);
                }
            }

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.ContainsKey"/>
            /// </summary>
            /// <param name="key">
            /// <see cref="IDictionary{TKey,TValueCollection}.ContainsKey"/>
            /// </param>
            /// <returns>
            /// <see cref="IDictionary{TKey,TValueCollection}.ContainsKey"/>
            /// </returns>
            public bool ContainsKey(TKey key) {
                return storage.ContainsKey(key);
            }

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.Keys"/>
            /// </summary>
            public ICollection<TKey> Keys {
                get {
                    return storage.Keys;
                }
            }

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.Remove"/>
            /// <param name="key">
            /// <see cref="IDictionary{TKey,TValueCollection}.Remove"/>
            /// </param>
            /// </summary>
            public bool Remove(TKey key) {
                return storage.Remove(key);
            }

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.TryGetValue"/>
            /// <param name="key">
            /// <see cref="IDictionary{TKey,TValueCollection}.TryGetValue"/>
            /// </param>
            /// <param name="value">
            /// <see cref="IDictionary{TKey,TValueCollection}.TryGetValue"/>
            /// </param>
            /// </summary>
            public bool TryGetValue(TKey key, out ICollection<TValue> value) {
                return TryGetValue(key, out value);
            }

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.Values"/>
            /// </summary>
            public ICollection<ICollection<TValue>> Values {
                get {
                    return storage.Values;
                }
            }

            /// <summary>
            /// <see cref="IDictionary{TKey,TValueCollection}.this"/>
            /// <param name="key">
            /// <see cref="IDictionary{TKey,TValueCollection}.this"/>
            /// </param>
            /// </summary>
            public ICollection<TValue> this[TKey key] {
                get {
                    return storage[key];
                }
                set {
                    storage[key] = value;
                }
            }

            #endregion

            #region ICollection<KeyValuePair<TKey,ICollection<TValue>>> Members

            /// <summary>
            /// <see cref="ICollection{TKeyTValuePair}.Add"/>
            /// </summary>
            /// <param name="item">
            /// <see cref="ICollection{TKeyTValuePair}.Add"/>
            /// </param>
            public void Add(KeyValuePair<TKey, ICollection<TValue>> item) {
                storage.Add(item.Key, item.Value);
            }

            /// <summary>
            /// <see cref="ICollection{TKeyTValuePair}.Clear"/>
            /// </summary>
            public void Clear() {
                storage.Clear();
            }

            /// <summary>
            /// <see cref="ICollection{TKeyTValuePair}.Contains"/>
            /// </summary>
            /// <param name="item">
            /// <see cref="ICollection{TKeyTValuePair}.Contains"/>
            /// </param>
            /// <returns>
            /// True, if all elements in <paramref name="item"/>.Value are
            /// associated with the key <paramref name="item"/>.Key in this
            /// dictionary. Otherwise, false is returned.
            /// </returns>
            public bool Contains(KeyValuePair<TKey, ICollection<TValue>> item) {
                if (this.ContainsKey(item.Key)) {
                    var storedValues = this[item.Key];
                    return item.Value.Intersect(storedValues).Count() == item.Value.Count();
                } else {
                    return false;
                }
            }

            /// <summary>
            /// <see cref="ICollection{TKeyTValuePair}.CopyTo"/>
            /// </summary>
            /// <param name="array">
            /// <see cref="ICollection{TKeyTValuePair}.CopyTo"/>
            /// </param>
            /// <param name="arrayIndex">
            /// <see cref="ICollection{TKeyTValuePair}.CopyTo"/>
            /// </param>
            public void CopyTo(KeyValuePair<TKey, ICollection<TValue>>[] array, int arrayIndex) {
                storage.ToArray().CopyTo(array, arrayIndex);
            }

            /// <summary>
            /// <see cref="ICollection{TKeyTValuePair}.Count"/>
            /// </summary>
            public int Count {
                get {
                    return storage.Count;
                }
            }

            /// <summary>
            /// <see cref="ICollection{TKeyTValuePair}.IsReadOnly"/>
            /// </summary>
            public bool IsReadOnly {
                get {
                    return false;
                }
            }

            /// <summary>
            /// Tries to remove all values in <paramref name="item"/>.Value
            /// from the entry associated with <paramref name="item"/>.Key.
            /// </summary>
            /// <param name="item">
            /// <see cref="ICollection{TKeyTValuePair}.Remove"/>
            /// </param>
            /// <returns>
            /// True, if at least one item has been removed. Otherwise, false
            /// will be returned.
            /// </returns>
            public bool Remove(KeyValuePair<TKey, ICollection<TValue>> item) {
                if (storage.ContainsKey(item.Key)) {
                    bool removedSomething = false;
                    foreach (TValue value in item.Value) {
                        removedSomething |= storage[item.Key].Remove(value);
                    }
                    return removedSomething;
                } else {
                    return false;
                }
            }

            #endregion

            #region IEnumerable<KeyValuePair<TKey,ICollection<TValue>>> Members

            /// <summary>
            /// <see cref="IEnumerable{TKeyTValuePair}.GetEnumerator"/>
            /// </summary>
            /// <returns>
            /// <see cref="IEnumerable{TKeyTValuePair}.GetEnumerator"/>
            /// </returns>
            public IEnumerator<KeyValuePair<TKey, ICollection<TValue>>> GetEnumerator() {
                return storage.GetEnumerator();
            }

            #endregion

            #region IEnumerable Members

            /// <summary>
            /// <see cref="IEnumerable.GetEnumerator"/>
            /// </summary>
            /// <returns>
            /// <see cref="IEnumerable.GetEnumerator"/>
            /// </returns>
            IEnumerator IEnumerable.GetEnumerator() {
                return storage.GetEnumerator();
            }

            #endregion
        }
    }
}
