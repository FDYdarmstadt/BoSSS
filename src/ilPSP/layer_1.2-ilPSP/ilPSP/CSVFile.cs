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
using System.IO;
using System.Linq;
using System.Text;

namespace ilPSP.Utils {
    
    /// <summary>
    /// Utils for IO of csv (comma-separated-vector) files
    /// </summary>
    public static class CSVFile {

        /// <summary>
        /// Writes a 'table' in a csv format to a file named <paramref name="filename"/>.
        /// </summary>
        public static void SaveToCSVFile<T>(this IDictionary<string, T> table, string filename, FileMode fm = FileMode.Create, char ColSep = '\t', bool writeHeader = true)
            where T : System.Collections.IEnumerable //    
        {
            using(var txt = new StreamWriter(new FileStream(filename, fm), new UTF8Encoding())) {
                WriteCSVToStream(table, txt, ColSep, writeHeader);
                txt.Flush();
            }
        }

        /*

        static void Test() {
            IDictionary<string, string[]> A = null;
            Dictionary<string, IEnumerable<string>> B = null;

            var R1 = M<string[], string>(A);
            //M<IEnumerable<string>, string>(R1);
            //M<string[], string>(A);

            var R2 = M2(A);
            var R3 = M2(R2);

        }

        static IDictionary<string, System.Collections.IEnumerable> M2<V>(IDictionary<string, V> tableA)
            where V : System.Collections.IEnumerable
        {
            return null;
        }



        static IDictionary<string, IEnumerable<T>> M<V,T>(IDictionary<string, V> tableA)
            where V : System.Collections.IEnumerable
        {
            //IDictionary<string, V> R = new Dictionary<string, V>();

            

            return null;
        }
         */


        /// <summary>
        /// Merge two tables verically.
        /// </summary>
        public static IDictionary<string, IList<V>> VerticalCat<T, V>(IDictionary<string, T> tableA, IDictionary<string, T> tableB, bool ColumnMerge = true)
            where T : IEnumerable<V> //
        {
            Dictionary<string, IList<V>> Ret = new Dictionary<string, IList<V>>();

            int LenA = -1; // number of rows in table A
            foreach(var kv in tableA) { // add all columns of table A to output table
                if(LenA < 0) {
                    LenA = kv.Value.Count();
                } else {
                    if(LenA != kv.Value.Count()) {
                        throw new ArgumentException("Illegal table A: column '" + kv.Value + "' has different length than the first column.");
                    }
                }
                Ret.Add(kv.Key, kv.Value.ToList());
            }
            if(LenA < 0)
                LenA = 0;


            int LenB = -1;
            foreach(var kv in tableB) {
                if(LenB < 0) {
                    LenB = kv.Value.Count();
                } else {
                    if(LenB != kv.Value.Count()) {
                        throw new ArgumentException("Illegal table B: column '" + kv.Value + "' has different length than the first column.");
                    }
                }
                
                if(ColumnMerge == false) {
                    if(!Ret.ContainsKey(kv.Key)) {
                        throw new ArgumentException("Tables A and B have different columns: unable to find column '" + kv.Key + "' from table B in table A. Maybe wanna set ColumnMerge to true?");
                    }
                } else {
                    if(!Ret.ContainsKey(kv.Key)) {
                        Ret.Add(kv.Key, new List<V>());
                        for(int j = 0; j < LenB; j++)
                            Ret[kv.Key].Add(default(V));
                    }
                }

                foreach(V vau in kv.Value)
                    Ret[kv.Key].Add(vau);
            }


            return Ret;

            //Dictionary<string, List<V>> Ret2 = new Dictionary<string, List<V>>();
            //foreach(var kv in Ret) {
            //    Ret2.Add(kv.Key, kv.Value);                
            //}
            ////Dictionary<string, T> Ret2 = new Dictionary<string, T>();
            ////foreach(var kv in Ret) {
            ////    Ret2.Add(kv.Key, kv.Value);
            ////}

            //return Ret2;
        }



        /// <summary>
        /// Writes a 'table' in a csv format to a stream <paramref name="txt"/>.
        /// </summary>
        public static void WriteCSVToStream<T>(this IDictionary<string, T> table, TextWriter txt, char ColSep = '\t', bool writeHeader = true)
            where T : System.Collections.IEnumerable //
        {
            string[] cols = table.Keys.ToArray();
            if(cols.Length <= 0)
                return;

      
            int L = table[cols[0]].Count();

            foreach(var kv in table) {
                if(kv.Value.Count() != L)
                    throw new NotSupportedException(string.Format("Unable to format vectors of different length as csv: vector {0} has {1} entries, vector {2} has {3} entries.",cols[0], L, kv.Key, kv.Value.Count()));
            }
             
            int I = cols.Length, IM = cols.Length - 1;
            if(writeHeader) {
                for(int i = 0; i < I; i++) {
                    txt.Write(cols[i]);
                    if(i < IM) 
                        txt.Write(ColSep);
                }
                txt.WriteLine();
            }

            System.Collections.IEnumerator[] CC = cols.Select(cn => table[cn].GetEnumerator()).ToArray();

            for(int l = 0; l < L; l++) {
                for(int i = 0; i < I; i++) {
                    bool mn = CC[i].MoveNext();
                    Debug.Assert(mn);
                    object obj = CC[i].Current;

                    if(obj == null) {
                        txt.Write("NULL");
                    } else if(obj is double) {
                        double? dd = (obj) as double?;
                        double ddd = (dd != null ? (double)dd : double.NaN);
                        txt.Write(ddd.ToString(System.Globalization.NumberFormatInfo.InvariantInfo));
                    } else if(obj is float) {
                        float? ff = (obj) as float?;
                        float fff = (ff != null ? (float)ff : float.NaN);
                        txt.Write(fff.ToString(System.Globalization.NumberFormatInfo.InvariantInfo));
                    } else {
                        txt.Write(obj.ToString());
                    }
                    if(i < IM) 
                        txt.Write(ColSep);
                }
                txt.WriteLine();
            }

        }

        /// <summary>
        /// Reads a 'table' in a csv format from a file named <paramref name="filename"/>.
        /// </summary>
        public static void ReadFromCSVFile(this IDictionary<string, IEnumerable<double>> table, string filename, char ColSep = '\t') {
            using(var txt = new StreamReader(filename, new UTF8Encoding())) {
                ReadCSVFromStream(table, txt, ColSep);
            }
        }

        /// <summary>
        /// Reads a 'table' in a csv format from <paramref name="txt"/>.
        /// </summary>
        public static void ReadCSVFromStream(this IDictionary<string, IEnumerable<double>> table, TextReader txt, char ColSep = '\t') {
            if(table.Count > 0)
                throw new ArgumentException("Expecting an empty table");
            
            
            char[] _colSep = new char[] { ColSep };
            string line;
            line = txt.ReadLine();
            if(line == null)
                return;
            string[] ColNames = line.Split(_colSep, StringSplitOptions.RemoveEmptyEntries);
            int I = ColNames.Length;

            List<double>[] ColCont = new List<double>[I];
            for(int i  =0; i < I; i++) {
                ColCont[i] = new List<double>();
            }

            int iLine = 1;
            for(line = txt.ReadLine(); line != null; line = txt.ReadLine()) {
                iLine++;
                string[] vals = line.Split(_colSep, StringSplitOptions.RemoveEmptyEntries);
                if(vals.Length != I)
                    throw new IOException(string.Format("Errors in CSV in line {0}; expecting {1} value(s), but found {2}.", iLine, I, vals.Length));
                
                for(int i = 0; i < I; i++) {
                    ColCont[i].Add(double.Parse(vals[i], System.Globalization.NumberFormatInfo.InvariantInfo));
                }
            }

            for(int i = 0; i < I; i++) {
                table.Add(ColNames[i], ColCont[i].ToArray());
            }
        }

        /// <summary>
        /// Reads a 'table' in a csv format from a file named <paramref name="filename"/>.
        /// </summary>
        public static void ReadFromCSVFile(this IDictionary<string, IEnumerable<string>> table, string filename, char ColSep = '\t') {
            using(var txt = new StreamReader(filename, new UTF8Encoding())) {
                ReadCSVFromStream(table, txt, ColSep);
            }
        }

        
        /// <summary>
        /// Reads a 'table' in a csv format from a file named <paramref name="filename"/>.
        /// </summary>
        public static void ReadFromCSVFile(this IDictionary<string, System.Collections.IEnumerable> table, string filename, char ColSep = '\t', params Func<string, object>[] parsers) {
            using(var txt = new StreamReader(filename, new UTF8Encoding())) {
                ReadCSVFromStream(table, txt, ColSep, parsers);
            }
        }

        /// <summary>
        /// Reads a 'table' in a csv format from <paramref name="txt"/>.
        /// </summary>
        public static void ReadCSVFromStream(this IDictionary<string, IEnumerable<string>> table, TextReader txt, char ColSep = '\t') {
            if(table.Count > 0)
                throw new ArgumentException("Expecting an empty table");


            char[] _colSep = new char[] { ColSep };
            string line;
            line = txt.ReadLine();
            if(line == null)
                return;
            string[] ColNames = line.Split(_colSep, StringSplitOptions.RemoveEmptyEntries);
            int I = ColNames.Length;

            List<string>[] ColCont = new List<string>[I];
            for(int i  =0; i < I; i++) {
                ColCont[i] = new List<string>();
            }

            int iLine = 1;
            for(line = txt.ReadLine(); line != null; line = txt.ReadLine()) {
                iLine++;
                string[] vals = line.Split(_colSep, StringSplitOptions.RemoveEmptyEntries);
                if(vals.Length != I)
                    throw new IOException(string.Format("Errors in CSV in line {0}; expecting {1} value(s), but found {2}.", iLine, I, vals.Length));

                for(int i = 0; i < I; i++) {
                    ColCont[i].Add(vals[i]);
                }
            }

            for(int i = 0; i < I; i++) {
                table.Add(ColNames[i], ColCont[i].ToArray());
            }
        }

         /// <summary> 
        /// Finds the most derived common base class of all the provided types, or System.Object if there is no common base class.
        /// </summary>
        static Type CommonBaseClass(System.Collections.IEnumerable enu) {
            HashSet<Type> types = new HashSet<Type>();
            foreach(object o in enu) {
                if(o != null)
                    types.Add(o.GetType());

            }

            if(types.Count > 0)
                return CommonBaseClass(types.ToArray());
            else
                return typeof(object);
        }

        /// <summary>
        /// Gets the class hierarchy of the provided type, in order of derivation, e.g. : (System.Object,CustomBaseType,CustomConcreteType,...), or the singleton of System.Object type if the provided type is an interface or null .
        /// </summary>
        static IEnumerable<Type> ClassHierarchy(this Type type) {
            if (type == null || type.IsInterface) type = typeof(object);
            var stack = new Stack<Type>();
            do {
                stack.Push(type);
                type = type.BaseType;
            } while (type != null);
            return stack;

        }


        /// <summary> 
        /// Finds the most derived common base class of all the provided types, or System.Object if there is no common base class.
        /// </summary>
        static Type CommonBaseClass(params Type[] types) {
            if (ReferenceEquals(types, null)) return typeof(object);
            types = types.Where(x => !ReferenceEquals(x, null)).Distinct().ToArray();
            switch (types.Length) {
                case 0: return typeof(object);
                case 1: return types[0].IsInterface ? typeof(object) : types[0];
                default:
                IEnumerable<IEnumerable<Type>> hierarchies = types.Select(ClassHierarchy).OrderBy(x => x.Count());
                Queue<Type> smallest = new Queue<Type>(hierarchies.First().Reverse());
                hierarchies = hierarchies.Skip(1);
                do {
                    int maxPossible = smallest.Count;
                    hierarchies = hierarchies.Select(each => each.Take(maxPossible));
                    Type candidate = smallest.Dequeue();
                    if (hierarchies.All(each => each.Last() == candidate))
                        return candidate;
                } while (smallest.Count > 1);
                return typeof(object);
            }
        }

        
        /// <summary>
        /// Reads a 'table' in a csv format from <paramref name="txt"/>.
        /// </summary>
        public static void ReadCSVFromStream(this IDictionary<string, System.Collections.IEnumerable> table, TextReader txt, char ColSep = '\t', params Func<string, object>[] parsers) {
            if(table.Count > 0)
                throw new ArgumentException("Expecting an empty table");


            char[] _colSep = new char[] { ColSep };
            string line;
            line = txt.ReadLine();
            if(line == null)
                return;
            string[] ColNames = line.Split(_colSep, StringSplitOptions.RemoveEmptyEntries);
            int I = ColNames.Length;

            List<object>[] ColCont = new List<object>[I];
            for(int i = 0; i < I; i++) {
                ColCont[i] = new List<object>();
            }

            int iLine = 1;
            for(line = txt.ReadLine(); line != null; line = txt.ReadLine()) {
                iLine++;
                string[] vals = line.Split(_colSep, StringSplitOptions.RemoveEmptyEntries);
                if(vals.Length != I)
                    throw new IOException(string.Format("Errors in CSV in line {0}; expecting {1} value(s), but found {2}.", iLine, I, vals.Length));

                for(int i = 0; i < I; i++) {
                    ColCont[i].Add(parsers[i](vals[i]));
                }
            }

            Array[] CastedColumns = new Array[I];
            for(int i = 0; i < I; i++) {
                var Col_i = ColCont[i];

                Type T = CommonBaseClass(Col_i);
                Array CastCol = Array.CreateInstance(T, Col_i.Count);
                for(int l = 0; l < Col_i.Count; l++)
                    CastCol.SetValue(Col_i[l], l);

                CastedColumns[i] = CastCol;
            }


            for(int i = 0; i < I; i++) {
                table.Add(ColNames[i], CastedColumns[i]);
            }
        }
    }
}
