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
using System.Globalization;
using System.IO;
using System.Linq;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution {

    /// <summary>
    /// By this class, arbitrary objects (e.g. residuals of solvers, query results) can be saved (<see cref="LogValue"/>)
    /// and are related to certain composite keys.
    /// Finally, the object for keys and values can be converted into a table <see cref="FormatTable"/>.
    /// </summary>
    /// <remarks>
    /// This class is designed to work as a state-machine:
    /// An internal 'current' key is defined by subsequent calls to <see cref="UpdateKey"/>.
    /// If a certain value is logged by calling <see cref="LogValue"/>, it is linked to the current key. 
    /// </remarks>
    public class QueryResultTable {

        /// <summary>
        /// contructor.
        /// </summary>
        public QueryResultTable() {
            this.CurrentKey = new CompositeKey();
            this.CurrentKeyHistory = new List<CompositeKey>();
            this.CurrentKeyHistory.Add(this.CurrentKey);
        }

        /// <summary>
        /// This class represents a key in the log, resp. a row in the table.
        /// Each key may consist of several sub-keys which can be arbitrary objects, e.g. iteration numbers,
        /// timestep numbers, parameter cases.
        /// </summary>
        sealed public class CompositeKey : ICloneable {

            private readonly Dictionary<string, object> SubkeyValues = new Dictionary<string, object>();

            /// <summary>
            /// You'll never guess it.
            /// </summary>
            public object Clone() {
                var r = new CompositeKey();
                r.SubkeyValues.AddRange(this.SubkeyValues);
                return r;
            }

            /// <summary>
            /// the value of a specific sub-key
            /// </summary>
            /// <param name="kn">the sub-key name</param>
            public object this[string kn] {
                get {
                    object kv;
                    SubkeyValues.TryGetValue(kn, out kv);
                    return kv;
                }
            }

            /// <summary>
            /// Modifies one sub-key in this composite key.
            /// </summary>
            /// <param name="SubkeyName">name of the sub-key that should be altered</param>
            /// <param name="SubkeyValue">new value of the sub-key</param>
            /// <returns>
            /// If the new value set for the sub-key is equal to the current value of the sub-key, this object itselv is returned;
            /// otherwise, a clone with an altered sub-key is returned;
            /// </returns>
            public CompositeKey Modify(string SubkeyName, object SubkeyValue) {
                CompositeKey ret = this;

                if (SubkeyValue == null) {
                    if (this[SubkeyName] != null) {
                        ret = this.CloneAs();
                        Debug.Assert(ret.SubkeyValues.Keys.Contains(SubkeyName));
                        ret.SubkeyValues.Remove(SubkeyName);
                    }
                } else if (this.SubkeyValues.Keys.Contains(SubkeyName)) {
                    if (!SubkeyName.Equals(this[SubkeyName])) {
                        ret = this.CloneAs();
                        ret.SubkeyValues[SubkeyName] = (SubkeyValue is ICloneable) ? ((ICloneable)SubkeyValue).Clone() : SubkeyValue;
                    }
                } else {
                    ret = this.CloneAs();
                    ret.SubkeyValues[SubkeyName] = (SubkeyValue is ICloneable) ? ((ICloneable)SubkeyValue).Clone() : SubkeyValue;
                }

                return ret;
            }

            /// <summary>
            /// the names of all keys present in 
            /// </summary>
            public IEnumerable<string> AllKeyNames {
                get {
                    return this.SubkeyValues.Keys;
                }
            }

            /// <summary>
            /// comparison of two composite keys
            /// </summary>
            public override bool Equals(object obj) {
                if (object.ReferenceEquals(this, obj)) {
                    return true;
                }

                CompositeKey otha = obj as CompositeKey;
                if (otha == null) {
                    return false;
                }

                if (this.SubkeyValues.Count != otha.SubkeyValues.Count
                    || this.SubkeyValues.GetHashCode() != otha.SubkeyValues.GetHashCode()) {
                    return false;
                }

                foreach (var kv in this.SubkeyValues) {
                    object o;
                    if (!otha.SubkeyValues.TryGetValue(kv.Key, out o)) {
                        return false;
                    } else {
                        if (!o.Equals(kv.Value)) {
                            return false;
                        }
                    }
                }

                return true;
            }

            /// <summary>
            /// h.c.
            /// </summary>
            public override int GetHashCode() {
                return this.SubkeyValues.GetHashCode();
            }
        }

        /// <summary>
        /// The content of the table, sorted first by value names, second by keys; <br/>
        /// 1st level key: value name 
        /// 2nd level key (value of inner dictionary): key status 
        /// value: the logged value
        /// </summary>
        Dictionary<string, Dictionary<CompositeKey, object>> Values = new Dictionary<string, Dictionary<CompositeKey, object>>();

        /// <summary>
        /// Returns the logged values for a given name
        /// </summary>
        /// <param name="valueName">The name of the value</param>
        /// <returns>The logged values</returns>
        public IEnumerable<KeyValuePair<CompositeKey, object>> this[string valueName] {
            get {
                return Values.ContainsKey(valueName) ? Values[valueName] : null;
            }
        }

        CompositeKey CurrentKey;

        /// <summary>
        /// whenever the current key is modified, it is pushed here; (gemurxe!)
        /// </summary>
        internal List<CompositeKey> CurrentKeyHistory;

        /// <summary>
        /// Modifies the current key, see <see cref="CompositeKey.Modify"/>.
        /// </summary>
        public void UpdateKey(string KeyName, object KeyValue) {
            CurrentKey = CurrentKey.Modify(KeyName, KeyValue);
            this.CurrentKeyHistory.Add(CurrentKey);
        }

        /// <summary>
        /// Logs some value (an arbitrary object, e.g. a residual norm) and links it to the currently set key.
        /// </summary>
        /// <param name="ValueName">Name for the value (column name)</param>
        /// <param name="ValueValue">object which represents the value.</param>
        public void LogValue(string ValueName, object ValueValue) {
            Dictionary<CompositeKey, object> dings;
            if (!Values.TryGetValue(ValueName, out dings)) {
                dings = new Dictionary<CompositeKey, object>();
                Values.Add(ValueName, dings);
            }

            if (dings.ContainsKey(this.CurrentKey)) {
                dings[this.CurrentKey] = ValueValue;
            } else {
                dings.Add(CurrentKey, ValueValue);
            }
        }

        /// <summary>
        /// Returns a collection of all sub-key names used in the lifetime of this object.
        /// </summary>
        public IEnumerable<string> AllSubkeyNames {
            get {
                HashSet<string> all = new HashSet<string>();

                foreach (var v in Values.Values) {
                    foreach (var w in v.Keys) {
                        foreach (string s in w.AllKeyNames) {
                            all.Add(s);
                        }
                    }
                }


                return all.ToArray();
            }
        }

        /// <summary>
        /// Returns a collection of all value names used in the lifetime of this object.
        /// </summary>
        public IEnumerable<string> AllValueNames {
            get {
                return Values.Keys;
            }
        }

        /// <summary>
        /// Returns a collection of all keys defined in the lifetime of this object.
        /// </summary>
        public IEnumerable<CompositeKey> AllKeys {
            get {
                HashSet<CompositeKey> all = new HashSet<CompositeKey>();

                foreach (var v in Values.Values) {
                    foreach (CompositeKey k in v.Keys) {
                        all.Add(k);
                    }
                }

                return all;
            }
        }

        /// <summary>
        /// sorts all logged values in a table
        /// </summary>
        /// <param name="KeyColumnNames">
        /// The column names for the key table, i.e. names of different keys.
        /// the index correlates with the 2nd index, i.e. the column index of <paramref name="KeyTable"/>
        /// </param>
        /// <param name="KeyTable">
        /// The values of all keys that were put into the logger:
        /// 1st/row index: all composite keys
        /// 2nd/column index: the components of the composite key (sub-key);
        /// </param>
        /// <param name="ValueColumnNames">
        /// column names for the value table, correlates with the 2nd index of <paramref name="ValueTable"/>;
        /// </param>
        /// <param name="ValueTable">
        /// The values of all looged values:
        /// 1st/row index: all composite keys
        /// 2nd/column index: the components of the composite key (sub-key);
        /// </param>
        /// <param name="RowFilter">
        /// select a certain subset of keys, resp. table rows.
        /// </param>
        public void FormatTable(out string[] KeyColumnNames, out string[] ValueColumnNames, out object[,] KeyTable, out object[,] ValueTable, IEnumerable<CompositeKey> RowFilter = null) {
            KeyColumnNames = this.AllSubkeyNames.ToArray();
            ValueColumnNames = this.AllValueNames.ToArray();

            CompositeKey[] keys;
            if (RowFilter == null) {
                keys = this.AllKeys.ToArray();
            } else {
                keys = this.AllKeys.Intersect(RowFilter).ToArray();
            }

            KeyTable = new object[keys.Length, KeyColumnNames.Length];
            ValueTable = new object[keys.Length, ValueColumnNames.Length];

            Dictionary<string, int> KeyColumnNamesIdx = new Dictionary<string, int>();
            int cnt = 0;
            foreach (var s in KeyColumnNames) {
                KeyColumnNamesIdx.Add(s, cnt);
                cnt++;
            }

            Dictionary<string, int> ValueColumnNamesIdx = new Dictionary<string, int>();
            cnt = 0;
            foreach (var s in ValueColumnNames) {
                ValueColumnNamesIdx.Add(s, cnt);
                cnt++;
            }

            Dictionary<CompositeKey, int> KeyRowIdx = new Dictionary<CompositeKey, int>();
            cnt = 0;
            foreach (var k in keys) {
                KeyRowIdx.Add(k, cnt);
                foreach (string keyName in k.AllKeyNames) {
                    int KeyCol = KeyColumnNamesIdx[keyName];
                    KeyTable[cnt, KeyCol] = k[keyName];
                }

                cnt++;
            }

            foreach (var kv1 in this.Values) {
                string ValName = kv1.Key;
                int ValCol = ValueColumnNamesIdx[ValName];

                foreach (var kv2 in kv1.Value) {
                    CompositeKey ky = kv2.Key;
                    object val = kv2.Value;

                    int RowIdx;
                    if (KeyRowIdx.TryGetValue(ky, out RowIdx)) {
                        ValueTable[RowIdx, ValCol] = val;
                    }
                }
            }
        }


        /// <summary>
        /// writes a tabular representation of the logged values to a string
        /// </summary>
        public override string ToString() {
            var stw = new StringWriter();
            this.WriteToStream(stw);
            return stw.ToString();
        }

        /// <summary>
        /// converts the table obtained by <see cref="FormatTable"/> into csv-text.
        /// </summary>
        /// <param name="s">speeratior in the text</param>
        /// <param name="sepString">string used to speerate columns in the text output</param>
        /// <param name="RowFilter">
        /// select a certain subset of keys, resp. table rows.
        /// </param>
        public void WriteToStream(TextWriter s, string sepString = "\t", IEnumerable<CompositeKey> RowFilter = null) {
            string[] KeyColumnNames;
            string[] ValueColumnNames;
            object[,] KeyTable;
            object[,] ValueTable;
            FormatTable(out KeyColumnNames, out ValueColumnNames, out KeyTable, out ValueTable, RowFilter);


            object[,] fullTable =
                ArrayTools.CatVert<object>(
                    ArrayTools.CatHoriz<object>(KeyColumnNames.ToRowVec(), ValueColumnNames.ToRowVec()),
                    ArrayTools.CatHoriz<object>(KeyTable, ValueTable));

            for (int i = 0; i < fullTable.GetLength(0); i++) {
                for (int j = 0; j < fullTable.GetLength(1); j++) {
                    var o = fullTable[i, j];

                    if (o == null) {
                        s.Write("null");
                    } else if (o is double) {
                        s.Write(((double)o).ToString(NumberFormatInfo.InvariantInfo));
                    } else if (o is float) {
                        s.Write(((float)o).ToString(NumberFormatInfo.InvariantInfo));
                    } else {
                        s.Write(o.ToString());
                    }

                    if (j + 1 < fullTable.GetLength(1))
                        s.Write(sepString);
                }
                s.WriteLine();
            }

        }

        /// <summary>
        /// see <see cref="WriteToStream"/>
        /// </summary>
        public void ToTextFile(string name) {
            using (var s = new StreamWriter(name)) {
                WriteToStream(s);
                s.Flush();
                s.Close();
            }
        }
    }
}
