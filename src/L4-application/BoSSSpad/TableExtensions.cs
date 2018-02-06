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

using BoSSS.Foundation.IO;
using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {
    /// <summary>
    /// Collection of extension methods to use <see cref="System.Data.DataTable"/>-objects in the BoSSSpad.
    /// </summary>
    public static class TableExtensions {

        /// <summary>
        /// The keys and queries (see <see cref="ISessionInfo.KeysAndQueries"/>) of all sessions in an enumeration
        /// <paramref name="sessions"/> in one table.
        /// </summary>
        public static DataTable GetSessionTable(this IEnumerable<ISessionInfo> sessions, Tuple<string, Func<ISessionInfo, object>>[] AdditionalColums) {

            Dictionary<string, object[]> Ret = new Dictionary<string, object[]>();

            for (int iSess = 0; iSess < sessions.Count(); iSess++) {
                var SS = sessions.ElementAt(iSess);
                var kq = SS.KeysAndQueries.ToList();

                // add additional columns
                kq.Add(new KeyValuePair<string, object>("Session", SS));
                kq.Add(new KeyValuePair<string, object>("RegularTerminated", !SS.Tags.Contains(BoSSS.Solution.Application.NOT_TERMINATED_TAG)));

                if(AdditionalColums != null) {
                    foreach(var t in AdditionalColums) {
                        object val;
                        try {
                            val = t.Item2(SS);
                        } catch(Exception) {
                            val = null;
                        }
                        kq.Add(new KeyValuePair<string, object>(t.Item1, val));
                    }
                }


                // convert to table
                foreach (var kv in kq) {
                    string ColumnName = kv.Key;
                    object ValueInCol = kv.Value;

                    if (Ret.ContainsKey(ColumnName)) {
                        object[] Column = Ret[ColumnName];
                        Debug.Assert(Column.Length == iSess);
                        Array.Resize(ref Column, Column.Length + 1);
                        Column[iSess] = ValueInCol;
                        Ret[ColumnName] = Column;
                    } else {
                        object[] newColumn = new object[iSess + 1];
                        newColumn[iSess] = ValueInCol;
                        Ret.Add(ColumnName, newColumn);
                    }
                }

                foreach (var kv in Ret.ToArray()) {
                    string ColumnName = kv.Key;
                    object[] Column = kv.Value;

                    Debug.Assert((Column.Length == iSess) || (Column.Length == iSess + 1));

                    if (Column.Length == iSess) {
                        Array.Resize(ref Column, Column.Length + 1);
                        Ret[ColumnName] = Column;
                    }
                }

            }


            return Ret.ToDataTable();
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


        static Type GetColType(object[] col) {
            return CommonBaseClass(col.Where(o => o != null).Select(o => o.GetType()).ToArray());
        }

        /// <summary>
        /// Converts a table which is stored as a dictionary into a <see cref="DataTable"/>-object.
        /// </summary>
        /// <param name="t">
        /// Input table:
        ///  - keys: column names
        ///  - values: each column of the table; all column-arrays are expected to have the same length.
        /// </param>
        public static DataTable ToDataTable(this IDictionary<string, object[]> t) {
            // check arguments
            // ===============
            if (t.Count <= 0)
                return new DataTable();
            int L = t.First().Value.Length;
            foreach (var kv in t) {
                if (kv.Value.Length != L)
                    throw new ArgumentException(string.Format("Input dictionary does not represent a table: column \"{0}\" has a different lenght than column \"{1}\".", kv.Key, t.First().Key));
            }


            DataTable R = new DataTable();

            foreach (KeyValuePair<string, object[]> col in t) {

                // Create first column and add to the DataTable.
                DataColumn column = new DataColumn();
                column.DataType = GetColType(col.Value);
                column.ColumnName = col.Key;
                column.Caption = col.Key;
                column.ReadOnly = false;
                column.Unique = false;
                column.DefaultValue = null;

                R.Columns.Add(column);
            }

            for (int i = 0; i < L; i++) {
                DataRow row = R.NewRow();

                foreach (KeyValuePair<string, object[]> col in t) {
                    row[col.Key] = col.Value[i] != null ? col.Value[i] : DBNull.Value;
                }

                R.Rows.Add(row);
            }

            return R;
        }


        /// <summary>
        /// The inverse of <see cref="ToDataTable(IDictionary{string, object[]})"/>
        /// </summary>
        /// <param name="tab"></param>
        /// <returns>
        /// Table as dictionary:
        ///  - keys: column names
        ///  - values: each column of the table; all column-arrays are expected to have the same length.
        /// </returns>
        public static IDictionary<string,object[]> ToDictionary(this DataTable tab) {
            Dictionary<string, object[]> R = new Dictionary<string, object[]>();

            string[] ColNames = tab.GetColumnNames();
            foreach(string ColNmn in ColNames) {
                object[] col = tab.GetColumn<object>(ColNmn);

                R.Add(ColNmn, col);
            }

            return R;
        }



        /// <summary>
        /// Serializes an entire table into a string.
        /// </summary>
        public static string Serialize(this DataTable tab) {
            JsonSerializer formatter = new JsonSerializer() {
                TypeNameHandling = TypeNameHandling.All,
                Formatting = Formatting.Indented,
                ReferenceLoopHandling = ReferenceLoopHandling.Error,
                NullValueHandling = NullValueHandling.Ignore,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor
            };
                
            string[] ColNameS = tab.GetColumnNames();
            Array[] Cols = ColNameS.Select(name => tab.GetColumn(name)).ToArray();

            var intermedTable = new Tuple<string[], Array[]>(ColNameS, Cols);

            using(var tw = new StringWriter()) {
                //tw.WriteLine(this.GetType().AssemblyQualifiedName);
                using(JsonWriter writer = new JsonTextWriter(tw)) {  // Alternative: binary writer: BsonWriter
                    formatter.Serialize(writer, intermedTable);
                }

                string Ret = tw.ToString();
                return Ret;
            }
        }
        
        /// <summary>
        /// De-serializes an entire table into a string.
        /// </summary>
        public static DataTable Deserialize(string data) {

            // part 1: load data
            // =================

            Tuple<string[], Array[]> intermedTable;
            {
                JsonSerializer formatter = new JsonSerializer() {
                    TypeNameHandling = TypeNameHandling.All,
                    Formatting = Formatting.Indented,
                    ReferenceLoopHandling = ReferenceLoopHandling.Error,
                    NullValueHandling = NullValueHandling.Ignore,
                    ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor
                };



                using (var tr = new StringReader(data)) {
                    using (JsonReader reader = new JsonTextReader(tr)) {
                        intermedTable = formatter.Deserialize<Tuple<string[], Array[]>>(reader);
                    }
                }
            }

            // part 2: convert & return
            // ========================
            {

                DataTable R = new DataTable();

                if (intermedTable.Item1.Length != intermedTable.Item2.Length)
                    throw new IOException();
                int NoOfCol = intermedTable.Item1.Length;

                int L = 0;
                for(int iCol = 0; iCol < NoOfCol; iCol++) {
                    string name = intermedTable.Item1[iCol];
                    Array colData = intermedTable.Item2[iCol];
                    if(iCol == 0) {
                        L = colData.Length;
                    } else {
                        if(L != colData.Length) {
                            throw new IOException();
                        }
                    }


                    // Create column and add to the DataTable.
                    DataColumn column = new DataColumn();
                    column.DataType = colData.GetType().GetElementType();
                    column.ColumnName = name;
                    column.Caption = name;
                    column.ReadOnly = false;
                    column.Unique = false;
                    column.DefaultValue = null;

                    R.Columns.Add(column);
                }

                for (int i = 0; i < L; i++) { // loop over rows
                    DataRow row = R.NewRow();
                    for (int iCol = 0; iCol < NoOfCol; iCol++) { // loop over columns
                        object dataItem = intermedTable.Item2[iCol].GetValue(i);
                        string name = intermedTable.Item1[iCol];
                        
                        row[name] = dataItem != null ? dataItem : DBNull.Value;

                    }
                    R.Rows.Add(row);
                }
                return R;
            }

        }


        /// <summary>
        /// Extracts one column from a table, and tries to convert the content into type <typeparamref name="T"/>.
        /// </summary>
        public static Array GetColumn(this DataTable tab, string ColName) {
            DataColumn col = tab.Columns[ColName];
            Type colType = col.DataType;

            object[] ColContent = tab.GetColumn<object>(ColName);
            int L = ColContent.Length;
            Array R = Array.CreateInstance(colType, L);

            for(int i = 0; i < L; i++) {
                R.SetValue(ColContent[i], i);
            }

            return R;
        }

        /// <summary>
        /// Extracts one column from a table, and tries to convert the content into type <typeparamref name="T"/>.
        /// </summary>
        public static T[] GetColumn<T>(this DataTable tab, string ColumnName) {
            T[] R = new T[tab.Rows.Count];

            //tab.Columns[ColumnName].
            for (int i = 0; i < R.Length; i++) {
                object val = tab.Rows[i][ColumnName];
                if (val.Equals(DBNull.Value)) {
                    R[i] = default(T);
                } else {
                    //? default(T) : ((T)val);
                    try {
                        if (typeof(T) == typeof(int)) {
                            R.SetValue(Convert.ToInt32(val), i);
                        } else if (typeof(T) == typeof(short)) {
                            R.SetValue(Convert.ToInt16(val), i);
                        } else if (typeof(T) == typeof(long)) {
                            R.SetValue(Convert.ToInt64(val), i);
                        } else if (typeof(T) == typeof(ushort)) {
                            R.SetValue(Convert.ToUInt16(val), i);
                        } else if (typeof(T) == typeof(uint)) {
                            R.SetValue(Convert.ToUInt32(val), i);
                        } else if (typeof(T) == typeof(ulong)) {
                            R.SetValue(Convert.ToUInt64(val), i);
                        } else if (typeof(T) == typeof(byte)) {
                            R.SetValue(Convert.ToByte(val), i);
                        } else if (typeof(T) == typeof(sbyte)) {
                            R.SetValue(Convert.ToSByte(val), i);
                        } else if (typeof(T) == typeof(double)) {
                            R.SetValue(Convert.ToDouble(val), i);
                        } else if (typeof(T) == typeof(char)) {
                            R.SetValue(Convert.ToChar(val), i);
                        } else if (typeof(T) == typeof(decimal)) {
                            R.SetValue(Convert.ToDecimal(val), i);
                        } else if (typeof(T) == typeof(bool)) {
                            R.SetValue(Convert.ToBoolean(val), i);
                        } else if (typeof(T) == typeof(float)) {
                            R.SetValue(Convert.ToSingle(val), i);
                        } else if (typeof(T) == typeof(string)) {
                            R.SetValue(Convert.ToString(val), i);
                        } else if (typeof(T) == typeof(DateTime)) {
                            R.SetValue(Convert.ToDateTime(val), i);
                        } else {
                            R[i] = (T)val;
                        }
                    } catch (InvalidCastException) {
                        R[i] = default(T);
                    }
                }
            }

            return R;
        }

        /// <summary>
        /// A string array, listing for each column in table <paramref name="tab"/> its type and name.
        /// </summary>
        public static string[] GetColumnInfo(this DataTable tab) {
            List<string> R = new List<string>();
            for (int iCol = 0; iCol < tab.Columns.Count; iCol++) {
                DataColumn col = tab.Columns[iCol];
                R.Add(string.Format("{0}:\"{1}\"", col.DataType.Name, col.ColumnName));
            }
            return R.ToArray();
        }

        /// <summary>
        /// The names of all columns in table <paramref name="tab"/>.
        /// </summary>
        public static string[] GetColumnNames(this DataTable tab) {
            List<string> R = new List<string>();
            for (int iCol = 0; iCol < tab.Columns.Count; iCol++) {
                DataColumn col = tab.Columns[iCol];
                R.Add(col.ColumnName);
            }
            return R.ToArray();
        }

        /// <summary>
        /// The data type of all columns in table <paramref name="tab"/>.
        /// </summary>
        public static Type[] GetColumnTypes(this DataTable tab) {
            List<Type> R = new List<Type>();
            for (int iCol = 0; iCol < tab.Columns.Count; iCol++) {
                DataColumn col = tab.Columns[iCol];
                R.Add(col.DataType);
            }
            return R.ToArray();
        }

        /// <summary>
        /// Extracts rows by a selector.
        /// </summary>
        /// <param name="Tab">The original table.</param>
        /// <param name="Selector">
        /// Maps (RowIndex, Row) to either true or false; those rows, for which this function evaluates to 
        /// true will be taken into the result table.
        /// </param>
        /// <param name="IgnoreOnException">
        /// If <paramref name="Selector"/> throws an exception in a certain row, the row will not be selected/not selected based on the value
        /// of <paramref name="IgnoreOnException"/>.
        /// </param>
        public static DataTable ExtractRows(this DataTable Tab, Func<int, IDictionary<string, object>, bool> Selector, bool IgnoreOnException = true) {
            // check arguments
            // ===============

            DataTable R = new DataTable();

            List<string> ColNames = new List<string>();
            foreach (DataColumn orgCol in Tab.Columns) {

                // Create first column and add to the DataTable.
                DataColumn column = new DataColumn();
                column.DataType = orgCol.DataType;
                column.ColumnName = orgCol.ColumnName;
                column.Caption = orgCol.Caption;
                column.ReadOnly = orgCol.ReadOnly;
                column.Unique = orgCol.Unique;
                column.DefaultValue = orgCol.DefaultValue;
                ColNames.Add(orgCol.ColumnName);

                R.Columns.Add(column);
            }

            int L = Tab.Rows.Count;
            int J = Tab.Columns.Count;
            for (int i = 0; i < L; i++) {
                DataRow orgRow = Tab.Rows[i];
                Dictionary<string, object> orgRowAsDict = new Dictionary<string, object>();
                foreach (string ColName in ColNames) {
                    object obj_ColName = orgRow[ColName];
                    if (obj_ColName == DBNull.Value) {
                        orgRowAsDict.Add(ColName, null);
                    } else {
                        orgRowAsDict.Add(ColName, obj_ColName);
                    }
                }

                bool take;
                try {
                    take = Selector(i, orgRowAsDict);
                } catch (Exception e) {
                    Console.WriteLine("Exception in the selection test of row {0}: {1}, Message: {2}.", i, e.GetType().Name, e.Message);
                    take = !IgnoreOnException;
                }

                if (take) {
                    DataRow rowNew = R.NewRow();
                    for (int j = 0; j < J; j++) {
                        rowNew[j] = orgRow[j];
                    }

                    R.Rows.Add(rowNew);
                }
            }

            Func<int, int> G = delegate (int i) {
                int j = i * 2;
                return j;
            };

            return R;
        }

        /// <summary>
        /// Extracts columns by name.
        /// </summary>
        public static DataTable ExtractColumns(this DataTable Tab, params string[] ColumnNames) {
            DataTable R = new DataTable();

            // check args
            // ==========
            foreach (string colName in ColumnNames) {
                var orgCol = Tab.Columns[colName];
                if (orgCol == null)
                    throw new ArgumentException("Column \"" + colName + "\" does not exist.");

                // Create first column and add to the DataTable.
                DataColumn column = new DataColumn();
                column.DataType = orgCol.DataType;
                column.ColumnName = orgCol.ColumnName;
                column.Caption = orgCol.Caption;
                column.ReadOnly = orgCol.ReadOnly;
                column.Unique = orgCol.Unique;
                column.DefaultValue = orgCol.DefaultValue;

                R.Columns.Add(column);
            }

            // extract columns
            // ===============
            int L = Tab.Rows.Count;
            for (int i = 0; i < L; i++) {
                DataRow rowSrc = Tab.Rows[i];
                DataRow rowDst = R.NewRow();

                foreach (string colName in ColumnNames) {
                    rowDst[colName] = rowSrc[colName];
                }

                R.Rows.Add(rowDst);
            }

            return R;
        }

        /// <summary>
        /// Prints table <paramref name="Tab"/> to the console.
        /// </summary>
        public static void Print(this DataTable Tab) {
            WriteCSVToStream(Tab, Console.Out, ' ', true, true, true);
        }

        /// <summary>
        /// Prints table <paramref name="Tab"/> to a CSV-file.
        /// </summary>
        public static void ToCSVFile(this DataTable Tab,
            string filename, FileMode fm = FileMode.Create,
            char ColSep = '\t', bool EnforceEqualColumns = false, bool writeHeader = true, bool writeRowIdx = false) {

            using (var txt = new StreamWriter(new FileStream(filename, fm), new UTF8Encoding())) {
                WriteCSVToStream(Tab, txt, ColSep, EnforceEqualColumns, writeHeader, writeRowIdx);
                txt.Flush();
            }
        }

        static string WriteWhite(int I) {
            using (var t = new StringWriter()) {
                for (int i = 0; i < I; i++) {
                    t.Write(" ");
                }
                return t.ToString();
            }
        }

        /// <summary>
        /// Prints table <paramref name="table"/> to a text-writer <paramref name="txt"/>.
        /// </summary>
        public static void WriteCSVToStream(this DataTable table, TextWriter txt,
            char ColSep, bool EnforceEqualColumns, bool writeHeader, bool writeRowIdx) {

            if (EnforceEqualColumns && ColSep != ' ')
                throw new ArgumentException("Setting 'EnforceEqualColumns' true requires the column separator to be the space character.");

            int NoOfCols = table.Columns.Count;
            int NoOfRows = table.Rows.Count;
            if (NoOfCols <= 0)
                return;
            string[] cols = new string[NoOfCols];
            for (int iCol = 0; iCol < cols.Length; iCol++)
                cols[iCol] = table.Columns[iCol].ColumnName;

            // convert ot text
            // ===============

            int ColOffset = writeRowIdx ? 1 : 0;
            int RowOffset = writeHeader ? 1 : 0;
            string[,] TextArray = new string[RowOffset + NoOfRows, ColOffset + NoOfCols];
            if (writeHeader && writeRowIdx)
                TextArray[0, 0] = "";

            if (writeHeader) {
                for (int iCol = 0; iCol < NoOfCols; iCol++) {
                    TextArray[0, iCol + ColOffset] = table.Columns[iCol].ColumnName;
                }
            }

            if (writeRowIdx) {
                for (int iRow = 0; iRow < NoOfRows; iRow++) {
                    TextArray[iRow + RowOffset, 0] = string.Format("{0}:", iRow);
                }
            }

            for (int iRow = 0; iRow < NoOfRows; iRow++) {
                DataRow row = table.Rows[iRow];
                for (int iCol = 0; iCol < NoOfCols; iCol++) {
                    object obj = row[iCol];

                    if (obj == null) {
                        TextArray[iRow + RowOffset, iCol + ColOffset] = "NULL";
                    } else if (obj == DBNull.Value) {
                        TextArray[iRow + RowOffset, iCol + ColOffset] = "NULL";
                    } else if (obj is double) {
                        double? dd = (obj) as double?;
                        double ddd = (dd != null ? (double)dd : double.NaN);
                        TextArray[iRow + RowOffset, iCol + ColOffset] = (ddd.ToString(System.Globalization.NumberFormatInfo.InvariantInfo));
                    } else if (obj is float) {
                        float? ff = (obj) as float?;
                        float fff = (ff != null ? (float)ff : float.NaN);
                        TextArray[iRow + RowOffset, iCol + ColOffset] = (fff.ToString(System.Globalization.NumberFormatInfo.InvariantInfo));
                    } else {
                        TextArray[iRow + RowOffset, iCol + ColOffset] = (obj.ToString());
                    }
                }
            }

            // column separation
            // =================

            if (EnforceEqualColumns) {
                int[] ColWidths = new int[TextArray.GetLength(1)];
                for (int iRow = 0; iRow < TextArray.GetLength(0); iRow++) {
                    for (int iCol = 0; iCol < TextArray.GetLength(1); iCol++) {
                        ColWidths[iCol] = Math.Max(ColWidths[iCol], TextArray[iRow, iCol].Length);
                    }
                }

                for (int iCol = 0; iCol < ColWidths.Length; iCol++) {
                    ColWidths[iCol] += 1;
                }

                for (int iRow = 0; iRow < TextArray.GetLength(0); iRow++) {
                    for (int iCol = 0; iCol < TextArray.GetLength(1); iCol++) {
                        int WhiteWidth = ColWidths[iCol] - TextArray[iRow, iCol].Length;
                        if (WhiteWidth <= 0)
                            throw new ApplicationException();
                        TextArray[iRow, iCol] = TextArray[iRow, iCol] + WriteWhite(WhiteWidth);
                    }
                }


            } else {
                for (int iRow = 0; iRow < TextArray.GetLength(0); iRow++) {
                    for (int iCol = 0; iCol < TextArray.GetLength(1) - 1; iCol++) {
                        TextArray[iRow, iCol] = TextArray[iRow, iCol] + ColSep;
                    }
                }
            }

            // output to text writer
            // =====================

            for (int iRow = 0; iRow < TextArray.GetLength(0); iRow++) {
                for (int iCol = 0; iCol < TextArray.GetLength(1); iCol++) {
                    txt.Write(TextArray[iRow, iCol]);
                    if (iCol == TextArray.GetLength(1) - 1 && iRow < TextArray.GetLength(0) - 1)
                        txt.WriteLine();
                }
            }
        }

        /// <summary>
        /// Creates an xy-plot form a table
        /// </summary>
        /// <param name="Tab"></param>
        /// <param name="ColName_ForXValues">
        /// Column name, where the values for the x-axis are taken.
        /// </param>
        /// <param name="ColName_ForYValues"></param>
        /// <param name="RowSelector">
        /// Selects, which table row will end up in which graph, resp. data group (<see cref="Plot2Ddata.dataGroups"/>).
        /// If the returned name is null, or if an exception is thrown, the respective data row will not be included in any graph.
        /// </param>
        /// <returns></returns>
        public static Plot2Ddata ToPlot(this DataTable Tab,
            string ColName_ForXValues, string ColName_ForYValues,
            PlotRowSelector RowSelector) {
            return Tab.ToPlot(
                delegate (int iRow, IDictionary<string, object> Row, out string GraphName, out BoSSS.Solution.Gnuplot.PlotFormat GraphFormat, out double xValue, out double yValue) {
                    RowSelector(iRow, Row, out GraphName, out GraphFormat);
                    xValue = Convert.ToDouble(Row[ColName_ForXValues]);
                    yValue = Convert.ToDouble(Row[ColName_ForYValues]);
                });
        }

        /// <summary>
        /// Creates an xy-plot form a table
        /// </summary>
        /// <param name="Tab"></param>
        /// <param name="RowSelector">
        /// Selects, which table row will end up in which graph, resp. data group (<see cref="Plot2Ddata.dataGroups"/>).
        /// If the returned name is null, or if an exception is thrown, the respective data row will not be included in any graph.
        /// </param>
        /// <returns></returns>
        public static Plot2Ddata ToPlot(this DataTable Tab, 
            PlotRowSelectorEx RowSelector) {
            

            Plot2Ddata ret = new Plot2Ddata();

            // loop over table rows
            // ====================
            string[] ColNames = Tab.GetColumnNames();
            int L = Tab.Rows.Count;
            int J = Tab.Columns.Count;
            for (int i = 0; i < L; i++) {
                DataRow orgRow = Tab.Rows[i];
                Dictionary<string, object> orgRowAsDict = new Dictionary<string, object>();
                foreach (string ColName in ColNames) {
                    object obj_ColName = orgRow[ColName];
                    if (obj_ColName == DBNull.Value) {
                        orgRowAsDict.Add(ColName, null);
                    } else {
                        orgRowAsDict.Add(ColName, obj_ColName);
                    }
                }

                string groupName;
                Solution.Gnuplot.PlotFormat graphFormat;
                double xValue;
                double yValue;
                try {
                    RowSelector(i, orgRowAsDict, out groupName, out graphFormat, out xValue, out yValue);
                } catch (Exception e) {
                    Console.WriteLine("Exception in the selection test of row {0}: {1}, Message: {2}.", i, e.GetType().Name, e.Message);
                    groupName = null;
                    graphFormat = null;
                    xValue = 0.0;
                    yValue = 0.0;
                }

                if (groupName != null) {
                    //double xValue = Convert.ToDouble(orgRowAsDict[ColName_ForXValues]);
                    //double yValue = Convert.ToDouble(orgRowAsDict[ColName_ForYValues]);

                    Plot2Ddata.XYvalues xyGroup = Array.Find(ret.dataGroups, xyG => xyG.Name.Equals(groupName));
                    if(xyGroup == null) {
                        xyGroup = new Plot2Ddata.XYvalues(groupName);
                        ArrayTools.AddToArray(xyGroup, ref ret.dataGroups);
                    }

                    ArrayTools.AddToArray(xValue, ref xyGroup.Abscissas);
                    ArrayTools.AddToArray(yValue, ref xyGroup.Values);

                    xyGroup.Format = graphFormat;
                }
            }

            // sort data
            // =========
            foreach(var xyGroup in ret.dataGroups) {
                Array.Sort(xyGroup.Abscissas, xyGroup.Values);
            }


            // return
            // ======
            return ret;
        }



        /// <summary>
        /// Creates an xy-plot form a table
        /// </summary>
        /// <param name="Tab"></param>
        /// <param name="ColName_ForXValues">
        /// Column name, where the values for the x-axis are taken.
        /// </param>
        /// <param name="ColName_ForYValues"></param>
        /// <param name="ColName_GroupSelection">
        /// Selects, which table row will end up in which graph, resp. data group (<see cref="Plot2Ddata.dataGroups"/>).
        /// </param>
        /// <returns></returns>
        public static Plot2Ddata ToPlot(this DataTable Tab, 
            string ColName_ForXValues, string ColName_ForYValues, 
            params string[] ColName_GroupSelection) {
            

            Plot2Ddata ret = new Plot2Ddata();

            // loop over table rows
            // ====================
            string[] ColNames = Tab.GetColumnNames();
            int L = Tab.Rows.Count;
            int J = Tab.Columns.Count;
            for (int i = 0; i < L; i++) {
                DataRow orgRow = Tab.Rows[i];
                Dictionary<string, object> orgRowAsDict = new Dictionary<string, object>();
                foreach (string ColName in ColNames) {
                    object obj_ColName = orgRow[ColName];
                    if (obj_ColName == DBNull.Value) {
                        orgRowAsDict.Add(ColName, null);
                    } else {
                        orgRowAsDict.Add(ColName, obj_ColName);
                    }
                }

                string groupName = "";
                try {
                    //groupName = GroupSelector(i, orgRowAsDict);

                    for(int iS = 0; iS < ColName_GroupSelection.Length; iS++) {
                        groupName += ColName_GroupSelection[iS] + orgRow[ColName_GroupSelection[iS]].ToString();
                        if(iS < ColName_GroupSelection.Length - 1)
                            groupName += "--";
                    }

                } catch (Exception e) {
                    Console.WriteLine("Exception in the selection test of row {0}: {1}, Message: {2}.", i, e.GetType().Name, e.Message);
                    groupName = null;
                }

                if (groupName != null) {
                    double xValue = Convert.ToDouble(orgRowAsDict[ColName_ForXValues]);
                    double yValue = Convert.ToDouble(orgRowAsDict[ColName_ForYValues]);

                    Plot2Ddata.XYvalues xyGroup = Array.Find(ret.dataGroups, xyG => xyG.Name.Equals(groupName));
                    if(xyGroup == null) {
                        xyGroup = new Plot2Ddata.XYvalues(groupName);
                        ArrayTools.AddToArray(xyGroup, ref ret.dataGroups);
                    }

                    ArrayTools.AddToArray(xValue, ref xyGroup.Abscissas);
                    ArrayTools.AddToArray(yValue, ref xyGroup.Values);
                }
            }

            // sort data
            // =========
            foreach(var xyGroup in ret.dataGroups) {
                Array.Sort(xyGroup.Abscissas, xyGroup.Values);
            }


            // return
            // ======
            return ret;
        }
    }

    /// <summary>
    /// Used by <see cref="TableExtensions.ToPlot(DataTable, string, string, PlotRowSelector)"/>, to select data rows which are put into a plot.
    /// </summary>
    /// <param name="iRow">
    /// Row index.
    /// </param>
    /// <param name="Row">
    /// One row of the data table, 
    /// </param>
    /// <param name="GraphName">
    /// On exit, the name of the graph; rows with the same name will be in one graph.
    /// *if null, the respective row is ignored.*
    /// </param>
    /// <param name="GraphFormat">
    /// Symbol, Color, dash type, etc.
    /// </param>
    public delegate void PlotRowSelector(int iRow, IDictionary<string, object> Row, out string GraphName, out BoSSS.Solution.Gnuplot.PlotFormat GraphFormat);
    
    /// <summary>
    /// Used by <see cref="TableExtensions.ToPlot(DataTable, PlotRowSelectorEx)"/>, to select data rows which are put into a plot.
    /// </summary>
    /// <param name="iRow">
    /// Row index.
    /// </param>
    /// <param name="Row">
    /// One row of the data table, 
    /// </param>
    /// <param name="GraphName">
    /// On exit, the name of the graph; rows with the same name will be in one graph.
    /// *if null, the respective row is ignored.*
    /// </param>
    /// <param name="GraphFormat">
    /// Symbol, Color, dash type, etc.
    /// </param>
    /// <param name="xValue">
    /// On exit, the x-value for the graph.
    /// </param>
    /// <param name="yValue">
    /// On exit, the y-value for the graph.
    /// </param>
    public delegate void PlotRowSelectorEx(int iRow, IDictionary<string, object> Row, out string GraphName, out BoSSS.Solution.Gnuplot.PlotFormat GraphFormat, out double xValue, out double yValue);
}
