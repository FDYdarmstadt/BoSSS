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
using System.IO;
using ilPSP.Utils;

namespace ilPSP {

    /// <summary>
    /// Extension methods for <see cref="IEnumerable{T}"/>
    /// </summary>
    public static class IEnumerableExtensions {

        /// <summary>
        /// Checks whether the given sequence is either null or empty.
        /// </summary>
        /// <typeparam name="TSource">The type of the source elements</typeparam>
        /// <param name="source">The sequence in question</param>
        /// <returns>
        /// True, if <paramref name="source"/> is either null or empty; False,
        /// otherwise.
        /// </returns>
        public static bool IsNullOrEmpty<TSource>(this IEnumerable<TSource> source) {
            if (source == null) {
                return true;
            }

            return !source.Any();
        }

        /// <summary>
        /// Checks whether all elements of the given sequence are uniform (or
        /// homogeneous) with respect to the given <paramref name="selector"/>.
        /// That is, it checks whether <paramref name="selector"/>(element)
        /// returns equal values for all elements of <paramref name="source"/>.
        /// </summary>
        /// <typeparam name="TSource">
        /// The type of the elements of <paramref name="source"/>.
        /// </typeparam>
        /// <typeparam name="TCompare">
        /// The of the values that are compared.
        /// </typeparam>
        /// <param name="source">
        /// A generic list of elements.
        /// </param>
        /// <param name="selector">
        /// A function that maps each element of <paramref name="source"/> onto
        /// the values that should be checked for homogeneity.
        /// </param>
        /// <returns>
        /// True, if <paramref name="source"/> is either null, empty or
        /// contains only values for which <paramref name="selector"/> returns
        /// equal values.
        /// </returns>
        public static bool IsUniform<TSource, TCompare>(this IEnumerable<TSource> source, Func<TSource, TCompare> selector) {
            if (source.IsNullOrEmpty()) {
                return true;
            }

            using (IEnumerator<TSource> enumerator = source.GetEnumerator()) {
                enumerator.MoveNext();
                TCompare referenceValue = selector(enumerator.Current);

                while (enumerator.MoveNext()) {
                    TCompare currentValue = selector(enumerator.Current);
                    if (!currentValue.Equals(referenceValue)) {
                        return false;
                    }
                }

                return true;
            }
        }

        /// <summary>
        /// Determines the index of the element with the maximum value with
        /// respect to the comparer of <typeparamref name="TCompare"/>.
        /// </summary>
        /// <typeparam name="TSource">
        /// The type of the source elements of the source sequence.
        /// </typeparam>
        /// <typeparam name="TCompare">
        /// The type of the elements to be compared. Its comparer defines what
        /// "maximum" means in this context.
        /// </typeparam>
        /// <param name="source">
        /// A generic sequence of elements.
        /// </param>
        /// <param name="selector">
        /// A function that maps each element of <paramref name="source"/> onto
        /// the values that should be compared to find the maximum.
        /// </param>
        /// <returns>
        /// The index of the element with the maximum value with respect to the
        /// comparer of <typeparamref name="TCompare"/>.
        /// </returns>
        public static int IndexOfMax<TSource, TCompare>(this IEnumerable<TSource> source, Func<TSource, TCompare> selector)
            where TCompare : IComparable<TCompare> {
            if (source.IsNullOrEmpty()) {
                return -1;
            }

            int maxIndex = 0;
            using (IEnumerator<TSource> enumerator = source.GetEnumerator()) {
                enumerator.MoveNext();

                TSource t = enumerator.Current;
                if (!enumerator.MoveNext()) {
                    return maxIndex;
                }

                TCompare maxValue = selector(t);
                TCompare currentValue;
                int i = 1;
                do {
                    if ((currentValue = selector(enumerator.Current)).CompareTo(maxValue) > 0) {
                        maxValue = currentValue;
                        maxIndex = i;
                    }
                    i++;
                } while (enumerator.MoveNext());
            }

            return maxIndex;
        }

        /// <summary>
        /// Determines the index of the element with the minimum value with
        /// respect to the comparer of <typeparamref name="TCompare"/>.
        /// </summary>
        /// <typeparam name="TSource">
        /// The type of the source elements of the source sequence.
        /// </typeparam>
        /// <typeparam name="TCompare">
        /// The type of the elements to be compared. Its comparer defines what
        /// "minimum" means in this context.
        /// </typeparam>
        /// <param name="source">
        /// A generic sequence of elements.
        /// </param>
        /// <param name="selector">
        /// A function that maps each element of <paramref name="source"/> onto
        /// the values that should be compared to find the minimum.
        /// </param>
        /// <returns>
        /// The index of the element with the minimum value with respect to the
        /// comparer of <typeparamref name="TCompare"/>.
        /// </returns>
        public static int IndexOfMin<TSource, TCompare>(this IEnumerable<TSource> source, Func<TSource, TCompare> selector)
            where TCompare : IComparable<TCompare> {
            if (source.IsNullOrEmpty()) {
                return -1;
            }

            int minIndex = 0;
            using (IEnumerator<TSource> enumerator = source.GetEnumerator()) {
                enumerator.MoveNext();

                TSource t = enumerator.Current;
                if (!enumerator.MoveNext()) {
                    return minIndex;
                }

                TCompare minValue = selector(t);
                TCompare currentValue;
                int i = 1;
                do {
                    if ((currentValue = selector(enumerator.Current)).CompareTo(minValue) < 0) {
                        minValue = currentValue;
                        minIndex = i;
                    }
                    i++;
                } while (enumerator.MoveNext());
            }

            return minIndex;
        }

        /// <summary>
        /// Returns the element at the index returned by
        /// <see cref="IndexOfMax"/>.
        /// </summary>
        /// <typeparam name="TSource"></typeparam>
        /// <typeparam name="TCompare"></typeparam>
        /// <param name="source"></param>
        /// <param name="selector"></param>
        /// <returns></returns>
        public static TSource ElementAtMax<TSource, TCompare>(this IEnumerable<TSource> source, Func<TSource, TCompare> selector)
            where TCompare : IComparable<TCompare> {
            if (source.IsNullOrEmpty()) {
                return default(TSource);
            }

            int maxIndex = source.IndexOfMax(selector);
            return source.ElementAt(maxIndex);
        }

        /// <summary>
        /// Returns the element at the index returned by
        /// <see cref="IndexOfMin"/>.
        /// </summary>
        /// <typeparam name="TSource"></typeparam>
        /// <typeparam name="TCompare"></typeparam>
        /// <param name="source"></param>
        /// <param name="selector"></param>
        /// <returns></returns>
        public static TSource ElementAtMin<TSource, TCompare>(this IEnumerable<TSource> source, Func<TSource, TCompare> selector)
            where TCompare : IComparable<TCompare> {
            if (source.IsNullOrEmpty()) {
                return default(TSource);
            }

            int minIndex = source.IndexOfMin(selector);
            return source.ElementAt(minIndex);
        }

        /// <summary>
        /// The index of an element in an enumeration, based on
        /// <see cref="object.ReferenceEquals"/>-comparison.
        /// </summary>
        /// <returns>
        /// The index of the first occurrence of <paramref name="elm"/>
        /// in <paramref name="List"/>, or -1, if not found.
        /// </returns>
        public static int IndexOf<T>(this IEnumerable<T> List, T elm) {
            return IndexOf(List, elm, (a, b) => object.ReferenceEquals(a, b));
        }

        /// <summary>
        /// The index of an element in an enumeration, based on the custom
        /// comparison <paramref name="comparer"/>.
        /// </summary>
        /// <returns>
        /// The index of the first occurrence of <paramref name="element"/> in
        /// <paramref name="seq"/>, or -1, if not found.
        /// </returns>
        public static int IndexOf<T>(this IEnumerable<T> seq, T element, Func<T, T, bool> comparer) {
            int i = 0;
            foreach (var e in seq) {
                if (comparer(e, element))
                    return i;
                i++;
            }
            return -1;
        }

        /// <summary>
        /// Determines the index of the entry in <paramref name="seq"/>,
        /// where <paramref name="condition"/> evaluates to true
        /// </summary>
        /// <param name="seq">
        /// The sequence to be searched
        /// </param>
        /// <param name="condition">
        /// A condition that is applied to each element of
        /// <paramref name="seq"/>.
        /// </param>
        /// <returns>
        /// The index of the entry in <paramref name="seq"/>,
        /// where <paramref name="condition"/> evaluates to true
        /// </returns>
        public static int IndexWhere<T>(this IEnumerable<T> seq, Func<T, bool> condition) {
            int cnt = 0;
            int ret = -1;
            int ifound = 0;
            foreach (var t in seq) {
                if (condition(t)) {
                    ret = cnt;
                    ifound++;
                }
                cnt++;
            }
            if (ifound > 1)
                throw new ArgumentException("Condition is fulfilled " + ifound + " times.");

            return ret;
        }


        /// <summary>
        /// Determines the index of the first entry in <paramref name="seq"/>,
        /// where <paramref name="condition"/> evaluates to true
        /// </summary>
        /// <param name="seq">
        /// The sequence to be searched
        /// </param>
        /// <param name="condition">
        /// A condition that is applied to each element of
        /// <paramref name="seq"/>.
        /// </param>
        /// <returns>
        /// The index of the first entry in <paramref name="seq"/>,
        /// where <paramref name="condition"/> evaluates to true
        /// </returns>
        public static int FirstIndexWhere<T>(this IEnumerable<T> seq, Func<T, bool> condition) {
            int cnt = 0;
            foreach (var t in seq) {
                if (condition(t))
                    return cnt;
                cnt++;
            }
            return -1;
        }

        /// <summary>
        /// Determines the index of the last entry in <paramref name="seq"/>,
        /// where <paramref name="condition"/> evaluates to true
        /// </summary>
        /// <param name="seq">
        /// The sequence to be searched
        /// </param>
        /// <param name="condition">
        /// A condition that is applied to each element of
        /// <paramref name="seq"/>.
        /// </param>
        /// <returns>
        /// The index of the last entry in <paramref name="seq"/>,
        /// where <paramref name="condition"/> evaluates to true
        /// </returns>
        public static int LastIndexWhere<T>(this IEnumerable<T> seq, Func<T, bool> condition) {
            int r = -1;
            int cnt = 0;
            foreach (var t in seq) {
                if (condition(t))
                    r = cnt;
                cnt++;
            }
            return r;
        }

        /// <summary>
        /// Variant of 'Contains' with custom equality comparer
        /// </summary>
        /// <typeparam name="T">
        /// Type of the sequence elements
        /// </typeparam>
        /// <typeparam name="R">
        /// Type of the element we're looking for
        /// </typeparam>
        /// <param name="seq">
        /// The sequence of elements to be searched
        /// </param>
        /// <param name="element">
        /// The element we're looking for
        /// </param>
        /// <param name="comparer">
        /// A comparer that defines a measure for equality between the elements
        /// of <paramref name="seq"/> and <paramref name="element"/>
        /// </param>
        /// <returns>
        /// True, if <paramref name="element"/> is contained in
        /// <paramref name="seq"/>; false otherwise.
        /// </returns>
        public static bool Contains<T, R>(this IEnumerable<T> seq, R element, Func<T, R, bool> comparer) {
            return seq.Any(s => comparer(s, element));
        }

        /// <summary>
        /// Variant of 'Contains' that uses
        /// <see cref="object.ReferenceEquals"/> to compare objects.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the sequence elements
        /// </typeparam>
        /// <param name="seq">
        /// The sequence to be searched
        /// </param>
        /// <param name="element">
        /// The element we're looking for
        /// </param>
        /// <returns>
        /// True, if <paramref name="element"/> is contained in
        /// <paramref name="seq"/>; false otherwise.
        /// </returns>
        public static bool ContainsExactly<T>(this IEnumerable<T> seq, T element) {
            return seq.Any(s => object.ReferenceEquals(s, element));
        }

        /// <summary>
        /// Applies the given <paramref name="action"/> to each element of the
        /// given sequence of elements.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the element in the sequence.
        /// </typeparam>
        /// <param name="source">
        /// The sequence of items the <paramref name="action"/> should be
        /// applied to.
        /// </param>
        /// <param name="action">
        /// The action to be performed on each element.
        /// </param>
        public static void ForEach<T>(this IEnumerable<T> source, Action<T> action) {
            if (source == null) {
                return;
            }

            using (var enumerator = source.GetEnumerator()) {
                while (enumerator.MoveNext()) {
                    action(enumerator.Current);
                }
            }
        }

        /// <summary>
        /// Applies the given <paramref name="action"/> to each element of the
        /// given sequence of elements.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the element in the sequence.
        /// </typeparam>
        /// <param name="source">
        /// The sequence of items the <paramref name="action"/> should be
        /// applied to.
        /// </param>
        /// <param name="action">
        /// The action to be performed on each element.
        /// </param>
        public static void ForEach<T>(this IEnumerable<T> source, Action<int, T> action) {
            if (source == null) {
                return;
            }
            int i = 0;
            foreach (T x in source) {
                action(i, x);
                i++;
            }
        }

        /// <summary>
        /// Concatenates a single element to a given sequence (i.e., there's no
        /// need to convert <paramref name="element"/> first).
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="sequence"></param>
        /// <param name="element"></param>
        /// <returns>
        /// The sequence [<paramref name="sequence"/>, <paramref name="element"/> ]
        /// </returns>
        public static IEnumerable<T> Concat<T>(this IEnumerable<T> sequence, T element) {
            if (sequence == null) {
                yield break;
            }

            var enumerator = sequence.GetEnumerator();
            while (enumerator.MoveNext()) {
                yield return enumerator.Current;
            }

            yield return element;
        }

        /// <summary>
        /// converts an enumeration of strings <paramref name="s"/> into a single string, 
        /// where all entries are separated by separated by <paramref name="separator"/>
        /// </summary>
        static public string CatStrings(this IEnumerable<string> s, string separator) {
            StringWriter ret = new StringWriter();
            var t = s.ToArray();
            for (int i = 0; i < t.Length; i++) {
                ret.Write(t[i]);
                if (i < t.Length - 1)
                    ret.Write(separator);
            }
            return ret.ToString();
        }

        /// <summary>
        /// adds a set of objects to a collection.
        /// </summary>
        public static void AddRange<T>(this ICollection<T> target, IEnumerable<T> ObjectsToAdd) {
            ObjectsToAdd.ForEach(o => target.Add(o));
        }

        /// <summary>
        /// adds a set of objects to a collection.
        /// </summary>
        public static void AddRange<T>(this ICollection<T> target, params T[] ObjectsToAdd) {
            ObjectsToAdd.ForEach(o => target.Add(o));
        }

        /// <summary>
        /// converts a BitArray into a boolean array.
        /// </summary>
        public static bool[] ToBoolArray(this BitArray bitArray) {
            bool[] ret = new bool[bitArray.Length];
            for (int k = 0; k < ret.Length; k++)
                ret[k] = bitArray[k];
            return ret;
        }

        /// <summary>
        /// Sets the entries of some BitArray according to a boolean array.
        /// </summary>
        public static void FromBoolArray(this BitArray bitArray, bool[] A) {
            if(A.Length != bitArray.Length)
                throw new ArgumentException("Arrays must have the same length.");
            for(int i = 0; i < A.Length; i++) {
                bitArray[i] = A[i];
            }
        }


        /// <summary>
        /// for each item in <paramref name="keyCollection"/>, the mapping <paramref name="f"/> is executed and the returned key/value pairs are 
        /// collected in a dictionary.
        /// </summary>
        public static Dictionary<T, V> ToDictionary<S, T, V>(this ICollection<S> keyCollection, Func<S, KeyValuePair<T, V>> f) {
            var R = new Dictionary<T, V>();
            foreach (var key in keyCollection) {
                var kv = f(key);
                R.Add(kv.Key, kv.Value);
            }
            return R;
        }

   

        /// <summary>
        /// True, if all elements in <paramref name="A"/> are also in <paramref name="B"/>.
        /// </summary>
        public static bool IsSubsetOf<T>(this IEnumerable<T> A, IEnumerable<T> B) { 
            foreach (var a in A) {
                if (!B.Contains(a))
                    return false;
            }
            return true;
        }

        /// <summary>
        /// true, if all elements in <paramref name="A"/> are also in <paramref name="B"/> and vice-versa
        /// </summary>
        public static bool SetEquals<T>(this IEnumerable<T> A, IEnumerable<T> B) {
            return (A.IsSubsetOf<T>(B) && B.IsSubsetOf<T>(A));
        }
        
        /// <summary>
        /// Wraps each entity into an <see cref="SmartEnumerable{T}.Entry"/>
        /// and returns the result as an enumerable. See
        /// <see cref="SmartEnumerable{T}"/>.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="entities"></param>
        /// <returns></returns>
        public static SmartEnumerable<T> AsSmartEnumerable<T>(this IEnumerable<T> entities) {
            return new SmartEnumerable<T>(entities);
        }

        /// <summary>
        /// Macht Array aus Zweier-Tupel.
        /// </summary>
        public static T[] TupleToArray<T>(this Tuple<T, T> t) {
            return new T[] { t.Item1, t.Item2 };
        }

        /// <summary>
        /// Creates a set from an enumeration.
        /// </summary>
        public static ISet<T> ToSet<T>(this IEnumerable<T> e) {
            return new HashSet<T>(e);
        }

        /// <summary>
        /// Tests if a string is empty or contains only whitespaces.
        /// </summary>
        public static bool IsEmptyOrWhite(this string s) {
            if (s == null)
                return true;

            int L = s.Length;
            for (int l = 0; l < L; l++) {
                if (!char.IsWhiteSpace(s[l]))
                    return false;
            }

            return true;
        }

        /// <summary>
        /// Returns the original <paramref name="sequence"/> except for
        /// elements that are equal to (<see cref="object.Equals(object)"/>) to
        /// to the <paramref name="excludedItem"/>.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="sequence"></param>
        /// <param name="excludedItem"></param>
        /// <returns></returns>
        public static IEnumerable<T> Except<T>(this IEnumerable<T> sequence, T excludedItem) {
            return sequence.Where(item => !item.Equals(excludedItem));
        }
    }
}
