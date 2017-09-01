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
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using Newtonsoft.Json;
using ilPSP;


namespace BoSSS.Foundation.IO {
    public partial class SessionInfo {

        [OnDeserialized]
        private void AfterDeserialisation(StreamingContext context) {
            if (this.m_KeysAndQueries == null)
                this.m_KeysAndQueries = new KeysDict();
            this.m_KeysAndQueries.m_Owner = this;
        }

        /// <summary>
        /// Custom (very-low-performance) dictionary for the <see cref="SessionInfo.KeysAndQueries"/>-property,
        /// triggers the serialization of the session info if something is changed.
        /// </summary>
        [Serializable]
        [DataContract]
        class KeysDict : IDictionary<string, object> {

            [DataMember]
            public SessionInfo m_Owner;

            [DataMember]
            string[] m_Keys = new string[0];

            [DataMember]
            object[] m_Values = new object[0];

            private void TestSerializibility(object o) {
                try {
                    using (var ms = new MemoryStream()) {
                        ms.Position = 0;


                        JsonSerializer formi = new JsonSerializer() {
                            NullValueHandling = NullValueHandling.Ignore,
                            TypeNameHandling = TypeNameHandling.Auto,
                            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                            ReferenceLoopHandling = ReferenceLoopHandling.Ignore

                        };

                        JsonWriter jWrt = new JsonTextWriter(new StreamWriter(ms));
                        formi.Serialize(jWrt, o);
                        jWrt.Flush();

                        ms.Position = 0;
                        JsonReader jRed = new JsonTextReader(new StreamReader(ms));
                        object backAgain = formi.Deserialize(jRed, o.GetType());

                        if (!(o.Equals(backAgain))) {
                            throw new ArgumentException("Value '" + o.ToString() + "' is not serializeable (not equal before and after serialization)");
                        }
                    }
                } catch (Exception e) {
                    throw new ArgumentException("Value '" + o.ToString() + "' is not serializeable (during serialization or de-serialization, the following exception occured: " + e.GetType().Name + ", " + e.Message +")");
                }
            }
            
            public object this[string key] {
                get {
                    
                    if (key == null)
                        throw new ArgumentNullException();
                    int idx = m_Keys.IndexOf(key,(a, b) => a.Equals(b));
                    if (idx < 0)
                        throw new KeyNotFoundException();

                    return m_Values[idx];
                }
                set {
                    if (key == null)
                        throw new ArgumentNullException();
                    TestSerializibility(value);

                    int idx = m_Keys.IndexOf(key, (a,b) => a.Equals(b));
                    if (idx < 0) {
                        Array.Resize(ref this.m_Keys, this.m_Keys.Length + 1);
                        Array.Resize(ref this.m_Values, this.m_Values.Length + 1);

                        this.m_Keys[this.m_Keys.Length - 1] = key;
                        this.m_Values[this.m_Values.Length - 1] = value;
                    } else {
                        m_Values[idx] = value;
                    }
                    if(m_Owner != null)
                        m_Owner.Save();
                }
            }

            public int Count {
                get {
                    return m_Keys.Length;
                }
            }

            public bool IsReadOnly {
                get {
                    return false;
                }
            }

            public ICollection<string> Keys {
                get {
                    return m_Keys.ToList().AsReadOnly();
                }
            }

            public ICollection<object> Values {
                get {
                    return m_Values.ToList().AsReadOnly();
                }
            }

            public void Add(KeyValuePair<string, object> item) {
                this.Add(item.Key, item.Value);
            }

            public void Add(string key, object value) {
                if (key == null)
                    throw new ArgumentNullException();
                int idx = m_Keys.IndexOf(key, (a, b) => a.Equals(b));
                if (idx >= 0)
                    throw new ArgumentException("Key '" + key + "' is already contained in this dictionary.");

                TestSerializibility(value);

                Array.Resize(ref this.m_Keys, this.m_Keys.Length + 1);
                Array.Resize(ref this.m_Values, this.m_Values.Length + 1);
                
                this.m_Keys[this.m_Keys.Length - 1] = key;
                this.m_Values[this.m_Values.Length - 1] = value;

                if (m_Owner != null)
                    m_Owner.Save();
            }

            public void Clear() {
                m_Values = new object[0];
                m_Keys = new string[0];
                if (m_Owner != null)
                    m_Owner.Save();
            }

            bool ICollection<KeyValuePair<string, object>>.Contains(KeyValuePair<string, object> item) {
                int idxK = m_Keys.IndexOf(item.Key, (a, b) => a.Equals(b));
                if(idxK < 0)
                    return false;
                return m_Values[idxK].Equals(item.Value);
            }

            public bool ContainsKey(string key) {
                return (m_Keys.IndexOf(key, (a, b) => a.Equals(b)) >= 0);
            }

            public void CopyTo(KeyValuePair<string, object>[] array, int arrayIndex) {
                for (int i = 0; i < m_Keys.Length; i++) {
                    array[i + arrayIndex] = new KeyValuePair<string, object>(m_Keys[i], m_Values[i]);
                }
            }

            IEnumerator<KeyValuePair<string, object>> IEnumerable<KeyValuePair<string, object>>.GetEnumerator() {
                return this.ToList().GetEnumerator();
            }

            bool ICollection<KeyValuePair<string, object>>.Remove(KeyValuePair<string, object> item) {
                int idxK = m_Keys.IndexOf(item.Key, (a, b) => a.Equals(b));
                if (idxK < 0)
                    return false;
                return m_Values[idxK].Equals(item.Value);
            }

            public bool Remove(string key) {
                if (key == null)
                    throw new ArgumentNullException();
                int idx = m_Keys.IndexOf(key, (a, b) => a.Equals(b));
                if (idx < 0) {
                    return false;
                } else {
                    for (int i = idx; i < m_Keys.Length - 1; i++) {
                        this.m_Keys[i] = this.m_Keys[i + 1];
                        this.m_Values[i] = this.m_Values[i + 1];
                    }
                    Array.Resize(ref this.m_Keys, this.m_Keys.Length - 1);
                    Array.Resize(ref this.m_Values, this.m_Values.Length - 1);
                    if (m_Owner != null)
                        m_Owner.Save();
                    return true;
                }
            }

            public bool TryGetValue(string key, out object value) {
                if (key == null)
                    throw new ArgumentNullException();
                int idx = m_Keys.IndexOf(key, (a, b) => a.Equals(b));
                if (idx < 0) {
                    value = null;
                    return false;
                } else {
                    value = m_Values[idx];
                    return true;
                }
            }

            IEnumerator IEnumerable.GetEnumerator() {
                return this.ToList().GetEnumerator();
            }
        }
    }
}
