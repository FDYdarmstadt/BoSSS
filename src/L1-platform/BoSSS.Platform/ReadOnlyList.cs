/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Platform {
#pragma warning disable 1591
    public interface IReadOnlyList<out T> : IEnumerable<T> {

        int Count {
            get;
        }

        T this[int index] {
            get;
        }
    }

    public class ReadOnlyList<T> : IReadOnlyList<T> {

        private IList<T> list;

        public ReadOnlyList() {
            list = new List<T>();
        }

        public ReadOnlyList(int capacity) {
            list = new List<T>(capacity);
        }

        public ReadOnlyList(IEnumerable<T> items) {
            this.list = items.ToList();
        }

        #region IReadonlyList<T> Members

        public int Count {
            get {
                return list.Count;
            }
        }

        public T this[int index] {
            get {
                return list[index];
            }
        }

        #endregion

        #region IEnumerable<T> Members

        public IEnumerator<T> GetEnumerator() {
            return list.GetEnumerator();
        }

        #endregion

        #region IEnumerable Members

        IEnumerator IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }

        #endregion
    }
#pragma warning restore 1591
}
