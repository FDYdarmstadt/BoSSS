using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;
using System.Text;

namespace BoSSS.Solution {
    public static class SolverExceptionLogger {

        public static void SaveException<T>(Exception e, Application<T> app) where T : AppControl, new() {

            if (app == null || app.IsDisposed) return; // object is disposed of, do nothing
                                                       // i'm actually not sure if this can happen
                                                       // some testing with a using(){} around a statement like here did not give like a null reference exception...
                                                       // apparently the object.Dispose() is called, but it is not collected and the reference in the event handler is still valid

            if (app.Control.savetodb && !app.CurrentSessionInfo.SuccessfulTermination) {
                //Console.WriteLine("**saving exception to database ...");

                SaveException(e, app.CurrentSessionInfo);
                // app.SaveToDatabase(-1, -1.0); this would be nice, but atm can cause deadlocks
            }

            //Console.WriteLine("**logging exception ...");
            MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int rank);
            Console.Error.WriteLine(("").PadRight(50, '='));
            Console.Error.WriteLine();
            Console.Error.WriteLine($"MPI rank {rank}: {e.GetType().Name } :");
            Console.Error.WriteLine(e.Message);
            Console.Error.WriteLine(e.StackTrace);
            Console.Error.WriteLine();
            Console.Error.WriteLine(("").PadRight(50, '='));
            Console.Error.WriteLine();

            ilPSP.Environment.StdOut.Flush();
            ilPSP.Environment.StdErr.Flush();
        }

        private static void SaveException(Exception e, SessionInfo s) {

            var dbDrv = s.Database.Controller.DBDriver;
            SingleExceptionDataSet.LogException(dbDrv, s, e);
        }

        /// <summary>
        /// Load all exceptions, logged for the selected sessions
        /// </summary>
        /// <param name="sessions"></param>
        /// <returns></returns>
        public static ExceptionDataSet LoadExceptions(IEnumerable<ISessionInfo> sessions) {

            // Create dataset with unique correlation between session and exception
            var set = new ExceptionDataSet();

            string[] exceptions = new string[0];
            foreach (var s in sessions) {
                string dir = DatabaseDriver.GetSessionDirectory(s);
                exceptions = Directory.GetFiles(dir, "Exception.*.txt").Select(f => Path.GetFullPath(f)).ToArray();

                foreach (var exception in exceptions) {
                    var eset = SingleExceptionDataSet.LoadException(exception);
                    set.AddLoggedException(eset);
                }
            }            
            return set;
        }

        internal class SingleExceptionDataSet : DataSet {

            SingleExceptionDataSet() {
                this.DataSetName = "Exception";
                CreateExceptionTable();

                void CreateExceptionTable() {
                    this.Tables.Add(new DataTable("Exception"));
                    DataColumn column;

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "SessionID";
                    column.ReadOnly = true;
                    column.Unique = true;
                    this.Tables["Exception"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(int);
                    column.ColumnName = "MPI-Rank";
                    this.Tables["Exception"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "ExceptionType";
                    this.Tables["Exception"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Method";
                    this.Tables["Exception"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Stacktrace";
                    this.Tables["Exception"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Message";
                    this.Tables["Exception"].Columns.Add(column);
                }
            }

            internal static void LogException(IDatabaseDriver drv, SessionInfo s, Exception e) {

                SingleExceptionDataSet set = new SingleExceptionDataSet();

                MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out int rank);

                // this is only called once, so every dataset has only one table with one row.
                // we use the dataset anyway, for simpler unification of all these outputs later on
                DataRow row;
                row = set.Tables["Exception"].NewRow();
                row["SessionID"] = s.ID.ToString();
                row["MPI-Rank"] = rank;
                row["ExceptionType"] = e.GetType().FullName;
                row["Method"] = $"{e.TargetSite.DeclaringType.FullName}.{ e.TargetSite.Name}";
                row["Stacktrace"] = e.StackTrace;
                row["Message"] = e.Message;
                set.Tables["Exception"].Rows.Add(row);

                set.WriteXml(drv.GetNewLogStream(s, "Exception"), XmlWriteMode.WriteSchema);
            }

            internal static SingleExceptionDataSet LoadException(string loggedException) {
                SingleExceptionDataSet set = new SingleExceptionDataSet();
                set.ReadXml(loggedException, XmlReadMode.ReadSchema);
                return set;
            }
        }

        public class ExceptionDataSet : DataSet {

            internal ExceptionDataSet() {
                this.DataSetName = "ExceptionsDB";
                CreateExceptionTable();
                CreateSessionTable();
                CreateRelation();

                void CreateExceptionTable() {
                    this.Tables.Add(new DataTable("Exceptions"));

                    DataColumn column;

                    column = new DataColumn();
                    column.DataType = typeof(int);
                    column.ColumnName = "ExceptionID";
                    column.ReadOnly = true;
                    column.Unique = true;
                    column.AutoIncrement = true;
                    this.Tables["Exceptions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "ExceptionType";
                    this.Tables["Exceptions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Method";
                    this.Tables["Exceptions"].Columns.Add(column);
                }

                void CreateSessionTable() {
                    this.Tables.Add(new DataTable("Sessions"));

                    DataColumn column;

                    column = new DataColumn();
                    column.DataType = typeof(int);
                    column.ColumnName = "ID";
                    column.ReadOnly = true;
                    column.Unique = true;
                    column.AutoIncrement = true;
                    this.Tables["Sessions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "SessionID";
                    this.Tables["Sessions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(int);
                    column.ColumnName = "MPI-Rank";
                    this.Tables["Sessions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "ExceptionType";
                    this.Tables["Sessions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(int);
                    column.ColumnName = "ExceptionID";
                    this.Tables["Sessions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Stacktrace";
                    this.Tables["Sessions"].Columns.Add(column);

                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Message";
                    this.Tables["Sessions"].Columns.Add(column);
                }

                void CreateRelation() {
                    DataColumn parentColumn = this.Tables["Exceptions"].Columns["ExceptionID"];
                    DataColumn childColumn = this.Tables["Sessions"].Columns["ExceptionID"];
                    DataRelation relation = new DataRelation("Session2Exception", parentColumn, childColumn);
                    this.Relations.Add(relation);
                }
            }

            public void SaveToHtml(string path) {
                foreach (DataTable tab in this.Tables) {
                    SaveToHTML(tab, tab.TableName + ".html");
                }

                void SaveToHTML(DataTable table, string filename) {
                    string fullPath;
                    if (path == null) {
                        fullPath = Path.Combine(Path.GetFullPath(Directory.GetCurrentDirectory()), filename);
                    } else {
                        fullPath = Path.Combine(Path.GetFullPath(path), filename);
                    }

                    File.WriteAllText(fullPath, ConvertDataTableToHTML(table));
                }

                string ConvertDataTableToHTML(DataTable dt) {
                    if (dt.Rows.Count == 0) return ""; // enter code here

                    StringBuilder builder = new StringBuilder();
                    builder.Append("<html>");
                    builder.Append("<head>");
                    builder.Append("<title>");
                    builder.Append(dt.TableName);
                    builder.Append("</title>");
                    builder.Append("</head>");
                    builder.Append("<body>");
                    builder.Append("<table border='1px' cellpadding='5' cellspacing='0' ");
                    builder.Append("style='border: solid 1px Silver; font-size: x-small;'>");
                    builder.Append("<tr align='left' valign='top'>");
                    foreach (DataColumn c in dt.Columns) {
                        builder.Append("<td align='left' valign='top'><b>");
                        builder.Append(System.Web.HttpUtility.HtmlEncode(c.ColumnName));
                        builder.Append("</b></td>");
                    }
                    builder.Append("</tr>");
                    foreach (DataRow r in dt.Rows) {
                        builder.Append("<tr align='left' valign='top'>");
                        foreach (DataColumn c in dt.Columns) {
                            builder.Append("<td align='left' valign='top'>");
                            builder.Append(System.Web.HttpUtility.HtmlEncode(r[c.ColumnName]));
                            builder.Append("</td>");
                        }
                        builder.Append("</tr>");
                    }
                    builder.Append("</table>");
                    builder.Append("</body>");
                    builder.Append("</html>");

                    return builder.ToString();
                }
            }

            internal void AddLoggedException(SingleExceptionDataSet exception) {

                var tab = exception.Tables["Exception"];

                foreach (DataRow row in tab.Rows) {
                    DataRow rrow;
                    rrow = this.Tables["Exceptions"].Select($"ExceptionType = '{row["ExceptionType"]}' AND Method = '{row["Method"]}'").SingleOrDefault();
                    if (rrow == null) {
                        rrow = this.Tables["Exceptions"].NewRow();
                        rrow["ExceptionType"] = row["ExceptionType"];
                        rrow["Method"] = row["Method"];
                        this.Tables["Exceptions"].Rows.Add(rrow);
                    }
                    int exceptionId = (int)rrow["ExceptionID"];

                    rrow = this.Tables["Sessions"].NewRow();
                    rrow["SessionID"] = row["SessionID"];
                    rrow["MPI-Rank"] = row["MPI-Rank"];
                    rrow["ExceptionType"] = row["ExceptionType"];
                    rrow["ExceptionID"] = exceptionId;
                    rrow["Stacktrace"] = row["Stacktrace"];
                    rrow["Message"] = row["Message"];
                    this.Tables["Sessions"].Rows.Add(rrow);
                }
            }
        }
    }
}
