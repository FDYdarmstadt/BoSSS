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

using BoSSS.Application.BoSSSpad;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Tracing;
using ilPSP.Utils;
using Mono.CSharp;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Xml;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Extensions methods for <see cref="ISessionInfo"/>.
    /// </summary>
    public static class ISessionInfoExtensions {

        /// <summary>
        /// Copies a session to another database.
        /// </summary>
        /// <param name="session">The session to be copied.</param>
        /// <param name="targetDB">The target database.</param>
        public static void Copy(this ISessionInfo session, IDatabaseInfo targetDB) {
            session.Database.Controller.CopySession(session, targetDB);
        }

        /// <summary>
        /// Moves a session to another database.
        /// </summary>
        /// <param name="session">The session to be moved.</param>
        /// <param name="targetDB">The target database.</param>
        public static void Move(this ISessionInfo session, IDatabaseInfo targetDB) {
            session.Database.Controller.MoveSession(session, targetDB);
        }

        /// <summary>
        /// Deletes a session from its database.
        /// </summary>
        /// <param name="session">The session to be deleted.</param>
        /// <param name="force">If true, the user will not be asked for confirmation.</param>
        public static void Delete(this ISessionInfo session, bool force = false) {
            if (force == true) {
                session.Database.Controller.DeleteSession(session);
                Console.WriteLine("Session " + session.ID.ToString() + " deleted.");
            } else {
                Console.WriteLine("Session: " + session.ToString());
                Console.Write("Do you really want to delete this session? [y/n]: ");
                string line = Console.ReadLine();

                if (line.ToLower().Equals("y")) {
                    session.Database.Controller.DeleteSession(session);
                    Console.WriteLine("Session " + session.ID.ToString() + " deleted.");
                } else {
                    Console.WriteLine("Session delete canceled.");
                }
            }
        }

        /// <summary>
        /// Deletes all sessions inside a collection. The user is not asked to
        /// confirm this action once, and not for every individual element in
        /// the collection.
        /// </summary>
        /// <param name="sessions">The entities to be deleted.</param>
        public static void DeleteAll(this IEnumerable<ISessionInfo> sessions) {
            if (sessions.IsNullOrEmpty()) {
                Console.WriteLine("Given collection is empty; nothing to delete");
                return;
            }

            Console.WriteLine("Sessions to delete:");
            foreach (ISessionInfo session in sessions) {
                Console.WriteLine(session.ToString());
            }
            Console.Write("Do you really want to delete these sessions? [y/n]: ");
            string line = Console.ReadLine();

            if (line.ToLower().Equals("y")) {
                foreach (ISessionInfo session in sessions) {
                    session.Database.Controller.DeleteSession(session);
                }
                Console.WriteLine("Sessions successfully deleted.");
            } else {
                Console.WriteLine("Session delete canceled.");
            }
        }

        /// <summary>
        /// Opens the directory where the export for the selected
        /// <paramref name="session"/> are stored in the explorer.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Obviously, this only works in Windows environments.
        /// </remarks>
        public static void OpenExportDirectory(this ISessionInfo session) {
            Process.Start(Utils.GetExportDirectory(session));
        }

        /// <summary>
        /// Prints the directory where the exports for the selected
        /// <paramref name="session"/> are stored to the console.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static void PrintExportDirectory(this ISessionInfo session) {
            Console.WriteLine(Utils.GetExportDirectory(session));
        }

        /// <summary>
        /// Opens the database directory where the files for the selected
        /// <paramref name="session"/> are stored in the explorer.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Obviously, this only works in Windows environments.
        /// </remarks>
        public static void OpenSessionDirectory(this ISessionInfo session) {
            Process.Start(DatabaseDriver.GetSessionDirectory(session));
        }

        /// <summary>
        /// Prints the directory where the files for the selected
        /// <paramref name="session"/> are stored to the console.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static void PrintSessionDirectory(this ISessionInfo session) {
            Console.WriteLine(DatabaseDriver.GetSessionDirectory(session));
        }

        /// <summary>
        /// Opens the specified <paramref name="TextFile"/>
        /// of the selected <paramref name="session"/>.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <param name="TextFile">
        /// The text file to be opened.
        /// </param>
        public static void OpenTextFile(this ISessionInfo session, string TextFile) {
            if (TextFile.Contains(".txt"))
                Process.Start(Path.Combine(DatabaseDriver.GetSessionDirectory(session), TextFile));
            else
                Process.Start(Path.Combine(DatabaseDriver.GetSessionDirectory(session), TextFile + ".txt"));
        }

        /// <summary>
        /// Convenience interface to create a
        /// <see cref="SessionExportInstruction"/>
        /// </summary>
        /// <param name="session"></param>
        /// <returns></returns>
        public static SessionExportInstruction Export(this ISessionInfo session) {
            return new SessionExportInstruction(session);
        }

        /// <summary>
        /// Reads a residual log file.
        /// </summary>
        /// <param name="session">
        /// The session for which the residuals are to be read.
        /// </param>
        /// <param name="norm">
        /// Identifier of the norm to be read from. Use e.g. "L2_abs" if you want
        /// to plot the residuals from residual-L2_abs.txt.
        /// </param>
        /// <param name="stride">
        /// Stride length, i.e. every n-th time-step will be read..
        /// </param>
        /// <param name="variables">
        /// The names of the variables to plot the read out.
        /// </param>
        public static ResidualLog Residuals(
            this ISessionInfo session, string norm, int stride = 1, params string[] variables) {

            return new ResidualLog(session, norm, variables, stride);
        }

        /// <summary>
        /// Reads queryResults.txt file.
        /// </summary>
        /// <param name="session">
        /// The session which queryResults should be read.
        /// </param>
        /// <param name="path">Path to the Query File</param>
        /// <returns>
        /// Dictionary of QueryResults.
        /// Key: QueryId
        /// Value: QueryValue
        /// </returns>
        public static IDictionary<string, double> QueryResults(this ISessionInfo session, string path = "queryResults.txt") {
            Dictionary<string, double> QueryDictionary = new Dictionary<string, double>();

            string QueryPath = Path.Combine(DatabaseDriver.GetSessionDirectory(session), path);

            using (StreamReader reader = new StreamReader(QueryPath)) {
                string[] QueryHeader = reader.ReadLine().Trim().Split(new char[] { '\t' });
                string[] QueryResults = reader.ReadLine().Trim().Split(new char[] { '\t' });

                if (QueryHeader.Length != QueryResults.Length) {
                    throw new ArgumentException("Broken queryResults.txt file.");
                }

                for (int i = 0; i < QueryHeader.Length; i++) {
                    QueryDictionary.Add(QueryHeader[i], Convert.ToDouble(QueryResults[i]));
                }
            }

            return QueryDictionary;
        }


        /// <summary>
        /// Reads the queryResult of a specific <paramref name="QueryId"/>.
        /// </summary>
        /// <param name="session">
        /// The session for which the query should be read.
        /// </param>
        /// <param name="QueryId">
        /// The Id of the query which should be read.
        /// </param>
        /// <param name="path">
        /// Name of query results file
        /// </param>
        /// <returns>
        /// The queryResult.
        /// </returns>
        public static double QueryResults(this ISessionInfo session, string QueryId, string path = "queryResults.txt") {
            double res;

            IDictionary<string, double> QueryDictionary = session.QueryResults();

            if (QueryDictionary.Count > 0) {
                try {
                    res = QueryDictionary[QueryId];
                } catch (KeyNotFoundException knf) {
                    Console.WriteLine(knf.Message);
                    res = double.NaN;
                }
            } else {
                res = double.NaN;
            }

            return res;
        }


        /// <summary>
        /// Lists all files in the session directory;
        /// </summary>
        /// <returns></returns>
        public static IEnumerable<string> FilesInSessionDir(this ISessionInfo session, string searcPattern = null) {
            if (searcPattern != null)
                return Directory.GetFiles(DatabaseDriver.GetSessionDirectory(session), searcPattern);
            else
                return Directory.GetFiles(DatabaseDriver.GetSessionDirectory(session));
        }

        /// <summary>
        /// True if the session is part of some parameter study.
        /// </summary>
        public static bool IsInParameterStudy(this ISessionInfo session) {
            return session.FilesInSessionDir("ParameterStudy*").Any();
        }

        /// <summary>
        /// Reads tabulated text files, and tries to interpret the content as floating-point values.
        /// </summary>
        /// <param name="session">
        /// The session for which the text file should be read.
        /// </param>
        /// <param name="TextFile">
        /// The text file, which should be read.
        /// </param>
        /// <param name="SepChars">
        /// A list of chars that are used as column separators in the given
        /// text file. If none are given, tabs and semicolons will be assumed
        /// </param>
        /// <returns>
        /// Returns the table as a dictionary.
        /// Key: The header of the table, i.e. column specifier.
        /// Value: The table values of the specific column.
        /// </returns>
        public static IDictionary<string, IList<double>> ReadTabulatedTextFileAsDoubles(this ISessionInfo session, string TextFile, params char[] SepChars) {
            if (SepChars == null || SepChars.Length <= 0)
                SepChars = new char[] { '\t', ';' };

            var raw = ReadTabulatedTextFileAsStrings(session, TextFile, SepChars);
            Dictionary<string, IList<double>> ret = new Dictionary<string, IList<double>>();

            foreach (var kv in raw) {
                IList<string> ColumnStrings = kv.Value;

                double[] ColumnDoubles = ColumnStrings.Select(delegate (string s) {
                    double ResDouble;
                    TimestepNumber ResTimestep;
                    if (Double.TryParse(s, out ResDouble))
                        return ResDouble;
                    else if (TimestepNumber.TryParse(s, out ResTimestep))
                        return (double)(ResTimestep[ResTimestep.Length - 1]);
                    else
                        throw new ArgumentException("Found not supported data format in table, see table entry '" + s + "'.");

                }).ToArray();
            }

            return ret;
        }

        /// <summary>
        /// Loads the profiling information for a session
        /// </summary>
        /// <param name="session"></param>
        /// <returns>
        /// An array of profiling trees, one for each MPI rank; th index into the returned array corresponds with the MPI rank.
        /// </returns>
        public static MethodCallRecord[] GetProfiling(this ISessionInfo session) {
            // find
            string sessDir = DatabaseDriver.GetSessionDirectory(session);
            string[] TextFils = Directory.GetFiles(sessDir, "profiling_bin.*.txt");
            if (TextFils.Count() <= 0)
                throw new IOException("Unable to find profiling information.");

            // sort according to process rank
            int[] Ranks = new int[TextFils.Length];
            for (int i = 0; i < Ranks.Length; i++) {
                var parts = TextFils[i].Split(new string[] { "profiling_bin.", ".txt" }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length <= 0)
                    throw new IOException("Unable to determine file rank from path '" + parts[i] + "'.");
                Ranks[i] = int.Parse(parts.Last());
                if (Ranks[i] < 0)
                    throw new IOException("Unable to determine file rank from path '" + parts[i] + "'.");
            }

            int MPISize = session.ComputeNodeNames.Count;
            if(MPISize != Ranks.Max() + 1) {
                Console.WriteLine("WARNING: mismatch between number of MPI ranks (" + MPISize + ") for session and max rank of profiling information (" + (Ranks.Max() + 1) + ").");
            }

            // load 
            var R = new MethodCallRecord[Ranks.Max() + 1];
            for(int i = 0; i < Ranks.Length; i++) {
                int rnk = Ranks[i];

                var f = TextFils[i];
                var JSON = File.ReadAllText(f);
                var mcr = MethodCallRecord.Deserialize(JSON);

                if (R[rnk] != null)
                    throw new IOException("It seems profiling info was written more than once for MPI rank " + rnk + ".");

                R[rnk] = mcr;
            }

            // return
            return R;
        }

        /// <summary>
        /// Reads tabulated text files.
        /// </summary>
        /// <param name="session">
        /// The session for which the text file should be read.
        /// </param>
        /// <param name="TextFile">
        /// The text file, which should be read.
        /// </param>
        /// <param name="SepChars">
        /// A list of chars that are used as column separators in the given
        /// text file. If none are given, tabs and semicolons will be assumed
        /// </param>
        /// <returns>
        /// Returns the table as a dictionary.
        /// Key: The header of the table, i.e. column specifier.
        /// Value: The table values of the specific column.
        /// </returns>
        public static IDictionary<string, IList<string>> ReadTabulatedTextFileAsStrings(this ISessionInfo session, string TextFile, params char[] SepChars) {
            if (SepChars == null || SepChars.Length <= 0)
                SepChars = new char[] { '\t', ';', ',', ' ' };
            IDictionary<string, IList<string>> Dict = new Dictionary<string, IList<string>>();

            //string TextFilePath;
            //if(TextFile.EndsWith(".txt"))
            //    TextFilePath = Path.Combine(DatabaseDriver.GetSessionDirectory(session), TextFile);
            //else
            //    TextFilePath = Path.Combine(DatabaseDriver.GetSessionDirectory(session), TextFile + ".txt");

            string sessDir = DatabaseDriver.GetSessionDirectory(session);
            string[] TextFils = Directory.GetFiles(sessDir, TextFile);
            if (TextFils.Count() <= 0)
                throw new ArgumentException("Cannot find file '" + TextFile + "' in session directory '" + sessDir + "'");
            if (TextFils.Count() > 1) {
                throw new ArgumentException("Multiple files matching '" + TextFile + "' in session directory '" + sessDir + "'");
            }
            string TextFilePath = TextFils[0];


            using (StreamReader reader = new StreamReader(TextFilePath)) {
                // Read header of table
                string[] Header = reader.ReadLine().Split(SepChars, StringSplitOptions.RemoveEmptyEntries);

                // Init list for results
                IList<string>[] Results = new List<string>[Header.Length];
                for (int col = 0; col < Header.Length; col++)
                    Results[col] = new List<string>();

                // Read table line by line
                string CurrentLine;
                do {
                    CurrentLine = reader.ReadLine();
                    if (CurrentLine != null) {
                        string[] Items = CurrentLine.Split(SepChars, StringSplitOptions.RemoveEmptyEntries);

                        if (Items.Length != Header.Length)
                            throw new ArgumentException("Given text file is not a table.");

                        for (int iCol = 0; iCol < Results.Length; iCol++) {
                            Results[iCol].Add(Items[iCol]);
                        }

                    }
                } while (CurrentLine != null);

                // Store results in dictionary
                for (int col = 0; col < Header.Length; col++) {
                    Dict.Add(Header[col], Results[col]);
                }
            }

            return Dict;
        }

        /// <summary>
        /// Reads text-file tables from multiple sessions and merges them.
        /// </summary>
        public static IDictionary<string, IList<string>> ReadTabulatedTextFileAsStrings(this IEnumerable<ISessionInfo> sessions, string TextFile, params char[] SepChars) {
            IDictionary<string, IList<string>> R = null;

            int i = 0;
            foreach (var s in sessions) {
                IDictionary<string, IList<string>> W = ReadTabulatedTextFileAsStrings(s, TextFile, SepChars);
                if (i == 0) {
                    R = W;
                } else {
                    R = CSVFile.VerticalCat<IList<string>, string>(R, W);
                }
                i++;
            }

            return R;
        }


        /// <summary>
        /// Adds one or more tags to a session.
        /// </summary>
        public static void AddTags(this ISessionInfo session, params string[] newTags) {
            string[] oldTags = session.Tags.ToArray();
            string[] sessTags = new string[oldTags.Length + newTags.Length];
            oldTags.CopyTo(sessTags, 0);
            newTags.CopyTo(sessTags, oldTags.Length);
            session.Tags = sessTags;
        }

        /// <summary>
        /// Removes a tag from a session.
        /// </summary>
        /// <param name="session">The session in question.</param>
        /// <param name="tag">The tag to delete.</param>
        public static void RemoveTag(this ISessionInfo session, string tag) {
            IList<string> sessTags = session.Tags.ToList();
            sessTags.Remove(tag);
            session.Tags = sessTags;
        }

        /// <summary>
        /// Returns the single session that matches the given
        /// <paramref name="name"/>.
        /// </summary>
        /// <param name="sessions">
        /// The list of sessions to be queried.
        /// </param>
        /// <param name="name">
        /// The name of the session in question.
        /// </param>
        /// <returns>
        /// A single session that matches <paramref name="name"/>, if it
        /// exists. Otherwise, an error will be thrown.
        /// </returns>
        public static ISessionInfo Find(this IEnumerable<ISessionInfo> sessions, string name) {
            var result = sessions.WithName(name);
            if (result.Count() == 0) {
                return null;
            } else if (result.Count() > 1) {
                throw new ArgumentException(String.Format(
                    "Found multiple sessions with name '{0}'", name));
            }

            return result.Single();
        }

        /// <summary>
        /// Finds all sessions that match the given <paramref name="name"/>.
        /// </summary>
        /// <param name="sessions">
        /// The list of sessions to be queried.
        /// </param>
        /// <param name="name">
        /// The name of the sessions in question.
        /// </param>
        /// <returns>
        /// A sequence of sessions matching with a name that matches
        /// <paramref name="name"/>.
        /// </returns>
        public static IEnumerable<ISessionInfo> WithName(this IEnumerable<ISessionInfo> sessions, string name) {
            return sessions.Where(s => s.Name.ToLowerInvariant().Equals(name.ToLowerInvariant()));
        }

        /// <summary>
        /// Finds all sessions that match the given <paramref name="guid"/>.
        /// </summary>
        /// <param name="sessions">
        /// The list of sessions to be queried.
        /// </param>
        /// <param name="guid">
        /// The guid of the sessions in question.
        /// </param>
        /// <returns>
        /// A sequence of sessions matching with an id that matches
        /// <paramref name="guid"/>.
        /// </returns>
        public static IEnumerable<ISessionInfo> WithGuid(this IEnumerable<ISessionInfo> sessions, PartialGuid guid) {
            return sessions.Where(s => s.ID == guid);
        }

        /// <summary>
        /// Returns the single session that matches the given
        /// <paramref name="guid"/>.
        /// </summary>
        /// <param name="sessions">
        /// The list of sessions to be queried.
        /// </param>
        /// <param name="guid">
        /// The guid of the session in question.
        /// </param>
        /// <returns>
        /// A single session that matches <paramref name="guid"/>, if it
        /// exists. Otherwise, an error will be thrown.
        /// </returns>
        public static ISessionInfo FindByGuid(this IEnumerable<ISessionInfo> sessions, PartialGuid guid) {
            return sessions.WithGuid(guid).Single();
        }

        /// <summary>
        /// Finds all sessions that match _all_ the given tags, see <see cref="ISessionInfo.Tags"/>.
        /// The search ignores upper-/lowercase.
        /// </summary>
        /// <param name="sessions">The list of sessions to be queried.</param>
        /// <param name="tags">The tags in question.</param>
        /// <returns>A collection of matching sessions.</returns>
        public static IEnumerable<ISessionInfo> WithTag(this IEnumerable<ISessionInfo> sessions, params string[] tags) {
            return sessions.Where(s =>
                    tags.All(tag =>
                        s.Tags
                        .Select(sessTag => sessTag.ToLowerInvariant())
                        .Contains(tag.ToLowerInvariant())));
        }

        /// <summary>
        /// Finds all sessions associated with the given
        /// <paramref name="projectName"/> via
        /// <see cref="ISessionInfo.ProjectName"/>
        /// </summary>
        /// <param name="sessions">The list of sessions to be queried.</param>
        /// <param name="projectName">The name of the relevant project</param>
        /// <returns>A collection of matching sessions.</returns>
        public static IEnumerable<ISessionInfo> WithProject(this IEnumerable<ISessionInfo> sessions, string projectName) {
            return sessions.Where(s => s.ProjectName == projectName);
        }

        /// <summary>
        /// Converts a list of sessions into a <see cref="Plot2Ddata"/> based on
        /// the information stored in the last time-step of the given
        /// <paramref name="sessions"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="ySelector">
        /// Selector for the relevant data at the data points.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and <paramref name="ySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> xSelector, Func<ITimestepInfo, double> ySelector) {
            return sessions.Select(s => s.Timesteps.Last()).ToDataSet(xSelector, ySelector);
        }

        /// <summary>
        /// Converts a list of sessions into a <see cref="Plot2Ddata"/> based on
        /// the information stored in the last time-step of the given
        /// <paramref name="sessions"/>, while grouping the results using
        /// <paramref name="groupKeySelector"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="ySelector">
        /// Selector for the relevant data at the data points.
        /// </param>
        /// <param name="groupKeySelector">
        /// A function defining a group id to each key-value pair
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and <paramref name="ySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> xSelector, Func<ITimestepInfo, double> ySelector, Func<ITimestepInfo, string> groupKeySelector) {
            return sessions.Select(s => s.Timesteps.Last()).ToDataSet(xSelector, ySelector, groupKeySelector);
        }

        /// <summary>
        /// Converts a list of sessions into a <see cref="Plot2Ddata"/> based on
        /// the information stored in the last time-step of the given
        /// <paramref name="sessions"/> and the results of a query named
        /// <paramref name="queryName"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="queryName">
        /// Name of the query whose results will be used as error measure.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and the results of the query named
        /// <paramref name="queryName"/>.
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> xSelector, string queryName) {
            return sessions.Select(s => s.Timesteps.Last()).ToDataSet(xSelector, queryName);
        }

        /// <summary>
        /// Converts a list of sessions into a <see cref="Plot2Ddata"/> based on
        /// the information stored in the last time-step of the given
        /// <paramref name="sessions"/>, while grouping the results using the
        /// DG degree of the field identified by
        /// <paramref name="groupFieldName"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="ySelector">
        /// Selector for the relevant data at the data points.
        /// </param>
        /// <param name="groupFieldName">
        /// A name of a DG field present in all <paramref name="sessions"/>
        /// whose DG degree will be used as a grouping function.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and <paramref name="ySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> xSelector, Func<ITimestepInfo, double> ySelector, string groupFieldName) {
            return sessions.Select(s => s.Timesteps.Last()).ToDataSet(xSelector, ySelector, groupFieldName);
        }

        /// <summary>
        /// Converts a list of sessions into a <see cref="Plot2Ddata"/> based on
        /// the information stored in the last time-step of the given
        /// <paramref name="sessions"/> and the results of a query named
        /// <paramref name="queryName"/>, while grouping the results using
        /// <paramref name="groupKeySelector"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="queryName">
        /// Name of the query whose results will be used as error measure.
        /// </param>
        /// <param name="groupKeySelector">
        /// A function defining a group id for each key-value pair
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and the results of the query named
        /// <paramref name="queryName"/>, grouped by means of
        /// <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> xSelector, string queryName, Func<ITimestepInfo, string> groupKeySelector) {
            return sessions.Select(s => s.Timesteps.Last()).ToDataSet(xSelector, queryName, groupKeySelector);
        }

        /// <summary>
        /// Converts a list of sessions into a <see cref="Plot2Ddata"/> based on
        /// the information stored in the last time-step of the given
        /// <paramref name="sessions"/> and the results of a query named
        /// <paramref name="queryName"/>, while grouping the results using the
        /// DG degree of the field identified by
        /// <paramref name="groupingFieldName"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions
        /// </param>
        /// <param name="xSelector">
        /// Selector for the data points.
        /// </param>
        /// <param name="queryName">
        /// Name of the query whose results will be used as error measure.
        /// </param>
        /// <param name="groupingFieldName">
        /// A name of a DG field present in all <paramref name="sessions"/>
        /// whose DG degree will be used as a grouping function.
        /// </param>
        /// <returns>
        /// A new <see cref="Plot2Ddata"/> filled with the data extracted via
        /// <paramref name="xSelector"/> and the results of the query named
        /// <paramref name="queryName"/>, grouped by means of the DG degree of
        /// the field named <paramref name="groupingFieldName"/>
        /// </returns>
        public static Plot2Ddata ToDataSet(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> xSelector, string queryName, string groupingFieldName) {
            return sessions.Select(s => s.Timesteps.Last()).ToDataSet(xSelector, queryName, groupingFieldName);
        }

        /// <summary>
        /// Extracts information about the grid convergence of the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// the given <paramref name="errorFunctional"/>.
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a grid convergence study. The
        /// DG degree of the relevant fields should be the same for all
        /// time-steps.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of
        /// <paramref name="errorFunctional"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> errorFunctional) {
            return sessions.Select(s => s.Timesteps.Last()).ToGridConvergenceData(errorFunctional);
        }

        /// <summary>
        /// Extracts information about the grid convergence from the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// the given <paramref name="errorFunctional"/>, while grouping the
        /// results using <paramref name="groupKeySelector"/>
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a grid convergence study. The
        /// DG degree of the relevant fields should be the same for all
        /// time-steps.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <param name="groupKeySelector">
        /// A grouping function.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of
        /// <paramref name="errorFunctional"/>, grouped by means of
        /// <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> errorFunctional, Func<ITimestepInfo, string> groupKeySelector) {
            return sessions.Select(s => s.Timesteps.Last()).ToGridConvergenceData(errorFunctional, groupKeySelector);
        }

        /// <summary>
        /// Extracts information about the grid convergence of the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// <paramref name="errorFunctional"/> and groups the results according
        /// to the DG degree of a field named
        /// <paramref name="groupFieldName"/>.
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a grid convergence study.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <param name="groupFieldName">
        /// A name of a relevant DG field which is present in all
        /// <paramref name="sessions"/> and whose DG degree will be used as a
        /// group key.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/>, the
        /// values are determined via the logarithm of
        /// <paramref name="errorFunctional"/> and the results are grouped with
        /// respect to the DG degree of a field named
        /// <paramref name="groupFieldName"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> errorFunctional, string groupFieldName) {
            return sessions.Select(s => s.Timesteps.Last()).ToGridConvergenceData(errorFunctional, groupFieldName);
        }

        /// <summary>
        /// Extracts information about the grid convergence from the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// the results of the query named <paramref name="queryName"/>.
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a grid convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used an error measure.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of the results of a query
        /// named <paramref name="queryName"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ISessionInfo> sessions, string queryName) {
            return sessions.Select(s => s.Timesteps.Last()).ToGridConvergenceData(queryName);
        }

        /// <summary>
        /// Extracts information about the grid convergence from the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// the results of the query named <paramref name="queryName"/>, while
        /// grouping the results using the DG degree of the field identified
        /// by <paramref name="groupingFieldName"/>.
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a grid convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used as an error measure.
        /// </param>
        /// <param name="groupingFieldName">
        /// A name of a DG field present in all <paramref name="sessions"/>
        /// whose DG degree will be used as a grouping function.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/>, the
        /// values are determined via the logarithm of the results of a query
        /// named <paramref name="queryName"/> and the results are grouped with
        /// respect to the DG degree of a field named
        /// <paramref name="groupingFieldName"/>.
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ISessionInfo> sessions, string queryName, string groupingFieldName) {
            return sessions.Select(s => s.Timesteps.Last()).ToGridConvergenceData(queryName, groupingFieldName);
        }

        /// <summary>
        ///  Extracts information about the grid convergence from the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// the results of the query named <paramref name="queryName"/>, while
        /// grouping the results using <paramref name="groupKeySelector"/>
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a grid convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used as an error measure.
        /// </param>
        /// <param name="groupKeySelector">
        /// A grouping function.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="IGridInfoExtensions.GetMeshSize"/> and the
        /// values are determined via the logarithm of the results of a query
        /// named <paramref name="queryName"/>, grouped by means of
        /// <paramref name="groupKeySelector"/>
        /// </returns>
        public static Plot2Ddata ToGridConvergenceData(this IEnumerable<ISessionInfo> sessions, string queryName, Func<ITimestepInfo, string> groupKeySelector) {
            return sessions.Select(s => s.Timesteps.Last()).ToGridConvergenceData(queryName, groupKeySelector);
        }

        /// <summary>
        /// Extracts information about the time convergence of the last
        /// time-steps of the given <paramref name="sessions"/> with respect to
        /// the given <paramref name="errorFunctional"/>
        /// </summary>
        /// <param name="sessions">
        /// A set of sessions representing a time convergence study.
        /// </param>
        /// <param name="errorFunctional">
        /// A function defining the error for a given
        /// <see cref="ITimestepInfo"/>.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="ITimestepInfoExtensions.GetTimeStepSize"/>
        /// and the values are determined via the logarithm of
        /// <paramref name="errorFunctional"/>.
        /// </returns>
        public static Plot2Ddata ToTimeConvergenceData(this IEnumerable<ISessionInfo> sessions, Func<ITimestepInfo, double> errorFunctional) {
            return sessions.Select(s => s.Timesteps.Last()).ToTimeConvergenceData(errorFunctional);
        }

        /// <summary>
        /// Extracts information about the time convergence of the given 
        /// <paramref name="sessions"/> with respect to the results of the
        /// query named <paramref name="queryName"/>
        /// </summary>
        /// <param name="sessions">
        /// A set of time-steps representing a time convergence study.
        /// </param>
        /// <param name="queryName">
        /// The name of a query whose results will be used as an error measure.
        /// </param>
        /// <returns>
        /// A <see cref="Plot2Ddata"/> where the abscissas are given by the
        /// logarithm of <see cref="ITimestepInfoExtensions.GetTimeStepSize"/>
        /// and the values are determined via the logarithm of the results of a
        /// query named <paramref name="queryName"/>
        /// </returns>
        public static Plot2Ddata ToTimeConvergenceData(this IEnumerable<ISessionInfo> sessions, string queryName) {
            return sessions.Select(s => s.Timesteps.Last()).ToTimeConvergenceData(queryName);
        }

        /// <summary>
        /// Estimates the errors of the fields with name
        /// <paramref name="fieldName"/> in the given last time-steps of the
        /// given <paramref name="sessions"/> by computing the errors with
        /// respect to the solution on the finest corresponding grid by making
        /// use of
        /// <see cref="BoSSS.Solution.Statistic.DGFieldComparison.ComputeErrors"/>.
        /// The result is then grouped according to the polynomial degree of
        /// field <paramref name="fieldName"/>.
        /// </summary>
        /// <param name="sessions">
        /// The sessions containing the fields whose errors should be
        /// estimated.
        /// </param>
        /// <param name="fieldName">
        /// The name of the DG field whose error should be estimated.
        /// </param>
        /// <returns>
        /// A data set containing information about the grid resolution and the
        /// corresponding errors with respect to the finest corresponding grid,
        /// grouped by the polynomial degree. Obviously, the time-step
        /// associated with the finest grid for each polynomial degree has an
        /// estimated error of zero (by definition) and is thus excluded from
        /// the result.
        /// </returns>
        public static Plot2Ddata ToEstimatedGridConvergenceData(this IEnumerable<ISessionInfo> sessions, string fieldName) {
            return sessions.Select(s => s.Timesteps.Last()).ToEstimatedGridConvergenceData(fieldName);
        }

        /// <summary>
        /// Tries to loads the control file of the given
        /// <paramref name="session"/> in the new REPL format. Note: Use
        /// <see cref="InteractiveBase.LoadAssembly"/> to load the
        /// corresponding solver assembly to be able to use solver-specific
        /// sub-classes of <see cref="AppControl"/>.
        /// </summary>
        /// <param name="session">
        /// The session whose configuration file should be loaded
        /// </param>
        /// <returns>
        /// The properly initialized configuration object
        /// </returns>
        /// <param name="t">
        /// Optional type of the control object to be instantiated. The type
        /// <see cref="AppControl"/> should work in all cases, but will not
        /// give access to solver-specific configuration options.
        /// </param>
        public static AppControl GetControl(this ISessionInfo session, Type t = null) {
        
            string sessionDir = DatabaseDriver.GetSessionDirectory(session);
            string path_script = Path.Combine(sessionDir, "Control-script.txt");
            string path_obj = Path.Combine(sessionDir, "Control-obj.txt");

            if (File.Exists(path_obj)) {
                string ctrlfileContent = File.ReadAllText(path_obj);

                return AppControl.Deserialize(ctrlfileContent);

            } else if (File.Exists(path_script)) {

                string ctrlfileContent = File.ReadAllText(path_script);

                AppControl.FromCode(ctrlfileContent, t, out AppControl ctrl, out AppControl[] ctrl_paramstudy);

                if (ctrl != null) {
                    return ctrl;
                } else if (ctrl_paramstudy != null) {
                    // We did a parameter study -> extract the correct control object from the list
                    int.TryParse(session.KeysAndQueries["id:pstudy_case"].ToString(), out int parameterStudyCase);

                    if (parameterStudyCase < 0 || parameterStudyCase >= ctrl_paramstudy.Count()) {
                        throw new Exception(
                            "Parameter study case index out of range. This should not have happened.");
                    }

                    return ctrl_paramstudy.ElementAt(parameterStudyCase);
                } else {
                    //throw new Exception(string.Format(
                    //    "Invalid control instruction: unable to cast the last"
                    //        + " result of the control file/cs-script of type {0} to type {1}",
                    //    controlObj.GetType().FullName,
                    //    typeof(T).FullName));
                    throw new NotSupportedException("unknown state.");
                }
            } else {
                throw new IOException("Unable to find control object (" + path_obj + ") or control script (" + path_script + ").");
            }
        }

        /// <summary>
        /// Compares the control files of <paramref name="session"/> and
        /// <paramref name="other"/> and displays the result using kdiff
        /// </summary>
        /// <param name="session"></param>
        /// <param name="other"></param>
        /// <returns></returns>
        public static string Diff(this ISessionInfo session, ISessionInfo other) {
            string sessionDir = DatabaseDriver.GetSessionDirectory(session);
            string controlFile = Path.Combine(sessionDir, "Control.txt");

            if (!File.Exists(controlFile)) {
                throw new FileNotFoundException(String.Format(
                    "Could not find a control file for session {0}", session.ID));
            }

            string otherControlFile = Path.Combine(
                DatabaseDriver.GetSessionDirectory(other), "Control.txt");
            if (!File.Exists(otherControlFile)) {
                throw new FileNotFoundException(string.Format(
                    "Could not find Control.txt for session {0}; either"
                        + " the control file is missing or the session is"
                        + " using the old control file version",
                    other.ID));
            }

            ProcessStartInfo psi = new ProcessStartInfo();
            string programFiles = System.Environment.GetFolderPath(System.Environment.SpecialFolder.ProgramFilesX86);
            psi.FileName = Path.Combine(Path.Combine(programFiles, "KDiff3"), "kdiff3.exe");
            psi.Arguments = String.Format("\"{0}\" \"{1}\"", controlFile, otherControlFile);
            psi.UseShellExecute = false;

            using (Process process = Process.Start(psi)) {
                return "Opening files in external application (kdiff3)";
            }
        }

        /// <summary>
        /// Retrieves an estimate for the run time of the given session
        /// </summary>
        /// <param name="session">
        /// The session of interest
        /// </param>
        /// <param name="firstIndex">
        /// optional Index for the first time-step, if the first x time-steps
        /// should be left out
        /// </param>
        /// <param name="lastIndex">
        /// optional Index for the last time-step, if the last x time-steps
        /// should be left out
        /// </param>
        /// <returns>
        /// The time elapsed between the first (<paramref name="firstIndex"/>)
        /// and the last (<paramref name="lastIndex"/>) time-step of
        /// <paramref name="session"/>
        /// </returns>
        public static TimeSpan GetApproximateRunTime(this ISessionInfo session, int firstIndex = -1, int lastIndex = -1) {
            var orderedTimesteps = session.Timesteps.WithoutSubSteps().OrderBy(t => t.TimeStepNumber);
            ITimestepInfo first = (firstIndex == -1) ? orderedTimesteps.First() : orderedTimesteps.Pick(firstIndex);
            ITimestepInfo last = (lastIndex == -1) ? orderedTimesteps.Last() : orderedTimesteps.Pick(lastIndex);
            return last.CreationTime - first.CreationTime;
        }

        /// <summary>
        /// Retrieves the average computing time (in seconds) for a single
        /// time-step in the given session.
        /// </summary>
        /// <param name="session">
        /// The session of interest.
        /// </param>
        /// <param name="firstIndex">
        /// optional Index for the first time-step, if the first x time-steps should be left out
        /// </param>
        /// <param name="lastIndex">
        /// optional Index for the last time-step, if the last x time-steps should be left out
        /// </param>
        /// <returns>
        /// The time elapsed between the first(firstIndex) and the last (lastIndex) time-step of
        /// <paramref name="session"/>, divided by the number of time-steps
        /// </returns>
        public static double GetAverageComputingTimePerTimestep(this ISessionInfo session, int firstIndex = -1, int lastIndex = -1) {
            var orderedTimesteps = session.Timesteps.WithoutSubSteps().OrderBy(t => t.TimeStepNumber);
            ITimestepInfo first = (firstIndex == -1) ? orderedTimesteps.First() : orderedTimesteps.Pick(firstIndex);
            ITimestepInfo last = (lastIndex == -1) ? orderedTimesteps.Last() : orderedTimesteps.Pick(lastIndex);
            return (last.CreationTime - first.CreationTime).TotalSeconds
                / (last.TimeStepNumber.MajorNumber - first.TimeStepNumber.MajorNumber);
        }

        /// <summary>
        /// Retrieves the average CPU (in seconds) for a single time-step in
        /// the given session.
        /// </summary>
        /// <param name="session">
        /// The session of interest.
        /// </param>
        /// <param name="firstIndex">
        /// Index for the first time-step, if the first x time-steps
        /// should be left out
        /// </param>
        /// <param name="lastIndex">
        /// Index for the last time-step, if the last x time-steps
        /// should be left out
        /// </param>
        /// <returns>
        /// The time elapsed between the first and the last time-step of
        /// <paramref name="session"/>, divided by the number of time-steps and
        /// multiplied by the number of processes
        /// </returns>
        public static double GetAverageCPUTimePerTimestep(this ISessionInfo session, int firstIndex = -1, int lastIndex = -1) {
            return GetAverageComputingTimePerTimestep(session, firstIndex, lastIndex) *
                session.ComputeNodeNames.Count;
        }

        /// <summary>
        /// Summarizes performance data about the given
        /// <paramref name="sessions"/> in a data set. In particular, this
        /// function considers the mesh size vs. average CPU time per time-step
        /// (<see cref="GetAverageCPUTimePerTimestep"/>).
        /// </summary>
        /// <param name="sessions">
        /// The sessions of interest.
        /// </param>
        /// <param name="groupKeySelector">
        /// Optional: A grouping function
        /// </param>
        /// <param name="firstIndex">
        /// Optional: Index for the first time-step, if the first x time-steps
        /// should be left out
        /// </param>
        /// <param name="lastIndex">
        /// Optional: Index for the last time-step, if the last x time-steps
        /// should be left out
        /// </param>
        /// <returns>
        /// For each session: <see cref="IGridInfoExtensions.GetMeshSize"/> vs.
        /// <see cref="GetAverageCPUTimePerTimestep"/> with a logarithmic
        /// scaling for both axes.
        /// </returns>
        public static Plot2Ddata ToPerformanceData(this IEnumerable<ISessionInfo> sessions, Func<ISessionInfo, string> groupKeySelector = null, int firstIndex = -1, int lastIndex = -1) {
            if (groupKeySelector == null) {
                groupKeySelector = (s => "allGroups");
            }

            return new Plot2Ddata(sessions.GroupBy(s => groupKeySelector(s)).
                Select(g => new KeyValuePair<string, double[][]>(
                    g.Key,
                    new double[][] {
                        g.Select(s => s.GetGrids().Single().GetMeshSize()).ToArray(),
                        g.Select(s => s.GetAverageCPUTimePerTimestep(firstIndex, lastIndex)).ToArray()
                    })).ToArray()).WithLogX().WithLogY();
        }

        /// <summary>
        /// Summarizes the information about the parallel speed up of the
        /// considered sessions.
        /// </summary>
        /// <param name="sessions">
        /// The sessions whose speedup shall be analyzed.
        /// </param>
        /// <param name="groupKeySelector">
        /// Optional: Grouping of the sessions, if multiple cases are considered.
        /// </param>
        /// <param name="firstIndex">
        /// Optional: Index for the first time-step, if the first x time-steps
        /// should be left out
        /// </param>
        /// <param name="lastIndex">
        /// Optional: Index for the last time-step, if the last x time-steps
        /// should be left out
        /// </param>
        /// <returns>
        /// A log-log data set where the x-values are given by the number of
        /// processes and the y-values are given by the speed-up with respect
        /// to the run with the fewest processes within
        /// <paramref name="sessions"/>. This data set automatically contains
        /// a data row displaying the optimal speed-up.
        /// </returns>
        public static Plot2Ddata ToSpeedUpData(this IEnumerable<ISessionInfo> sessions, Func<ISessionInfo, string> groupKeySelector = null, int firstIndex = -1, int lastIndex = -1) {
            var sortedSessions = sessions.OrderBy(s => s.ComputeNodeNames.Count());
            double minNodes = sortedSessions.First().ComputeNodeNames.Count();
            double minTime = sortedSessions.First().GetAverageComputingTimePerTimestep();

            if (groupKeySelector == null) {
                groupKeySelector = (s => "allGroups");
            }

            var data = sortedSessions.GroupBy(s => groupKeySelector(s)).
                Select(g => new KeyValuePair<string, double[][]>(
                    g.Key,
                    new double[][] {
                        g.Select(s => (double)s.ComputeNodeNames.Count()).ToArray(),
                        g.Select(s => minTime / s.GetAverageComputingTimePerTimestep(firstIndex, lastIndex)).ToArray()
                    }));
            var idealData = sortedSessions.GroupBy(s => groupKeySelector(s)).
                Select(g => new KeyValuePair<string, double[][]>(
                    "ideal",
                    new double[][] {
                        g.Select(s => (double)s.ComputeNodeNames.Count()).ToArray(),
                        g.Select(s => s.ComputeNodeNames.Count() / minNodes).ToArray()
                    }));

            return new Plot2Ddata(data.Concat(idealData).ToArray()).WithLogX().WithLogY();
        }

        /// <summary>
        /// Summarizes the information about the parallel efficiency of the
        /// considered sessions.
        /// </summary>
        /// <param name="sessions">
        /// The sessions whose efficiency shall be analyzed.
        /// </param>
        /// <param name="groupKeySelector">
        /// Optional: Grouping of the sessions, if multiple cases are considered.
        /// </param>
        /// <param name="firstIndex">
        /// Optional: Index for the first time-step, if the first x time-steps
        /// should be left out
        /// </param>
        /// <param name="lastIndex">
        /// Optional: Index for the last time-step, if the last x time-steps
        /// should be left out
        /// </param>
        /// <returns>
        /// A semi-logarithmic data set where the x-values are given by the
        /// number of processes and the y-values are given by the efficiency.
        /// That is, the actual speed-up compared to the nominal speed-up with
        /// respect to the run using the fewest processes within
        /// <paramref name="sessions"/>.
        /// </returns>
        public static Plot2Ddata ToEfficiencyData(this IEnumerable<ISessionInfo> sessions, Func<ISessionInfo, string> groupKeySelector = null, int firstIndex = -1, int lastIndex = -1) {
            var sortedSessions = sessions.OrderBy(s => s.ComputeNodeNames.Count());
            double minNodes = sortedSessions.First().ComputeNodeNames.Count();
            double minTime = sortedSessions.First().GetAverageComputingTimePerTimestep();

            Func<ISessionInfo, double> efficiency = s =>
                minTime / s.GetAverageComputingTimePerTimestep(firstIndex, lastIndex) *
                minNodes / s.ComputeNodeNames.Count();

            if (groupKeySelector == null) {
                groupKeySelector = (s => "allGroups");
            }

            return new Plot2Ddata(sortedSessions.GroupBy(s => groupKeySelector(s)).
                Select(g => new KeyValuePair<string, double[][]>(
                    g.Key,
                    new double[][] {
                        g.Select(s => (double)s.ComputeNodeNames.Count()).ToArray(),
                        g.Select(s => efficiency(s)).ToArray()
                    })).ToArray()).WithLogX();
        }


        /// <summary>
        /// Retrieves all sessions from <paramref name="sessions"/> who have
        /// been restarted from a session whose id matches the given (partial)
        /// <paramref name="id"/>.
        /// </summary>
        /// <param name="sessions">
        /// A list of sessions.
        /// </param>
        /// <param name="id">
        /// A partial guid of a session from which the retrieved sessions have
        /// been restarted
        /// </param>
        /// <returns>
        /// All sessions where <see cref="ISessionInfo.RestartedFrom"/> equals
        /// the given (partial) <paramref name="id"/>.
        /// </returns>
        public static IEnumerable<ISessionInfo> RestartedFrom(this IEnumerable<ISessionInfo> sessions, PartialGuid id) {
            return sessions.Where(s => id.Equals(s.RestartedFrom));
        }

        /// <summary>
        /// Extracts information about the run index within a parameter study
        /// from the description of a session.
        /// </summary>
        /// <param name="session">
        /// The session in question.
        /// </param>
        /// <returns>
        /// An integer greater than 0 zero if the information about the run
        /// could be extracted; -1 otherwise.
        /// </returns>
        public static int GetRun(this ISessionInfo session) {
            string desc = session.Description;
            int start = desc.LastIndexOf("run ") + 4;
            int end = desc.LastIndexOf(" of ");

            if (start < 0 || end < 0) {
                return -1;
            }

            return int.Parse(desc.Substring(start, end - start));
        }

        /// <summary>
        /// Filters all sessions which carry the <see cref="BoSSS.Solution.Application{T}.NOT_TERMINATED_TAG"/>-tag.
        /// </summary>
        public static IEnumerable<ISessionInfo> RunningOrCrashed(this IEnumerable<ISessionInfo> sessions) {
            return sessions.Where(S => S.Tags.Contains(BoSSS.Solution.Application.NOT_TERMINATED_TAG)).ToArray();
        }

        /// <summary>
        /// Filters all sessions that have been written on the last <paramref name="NoPrevDays"/> days;
        /// </summary>
        /// <param name="sessions"></param>
        /// <param name="NoPrevDays">
        /// 0: only today's sessions, <br/>
        /// 1: today and yesterday, etc.
        /// </param>
        /// <returns></returns>
        public static IEnumerable<ISessionInfo> LastDays(this IEnumerable<ISessionInfo> sessions, int NoPrevDays = 0) {
            var Today = DateTime.Now;
            var maxAge = new TimeSpan((NoPrevDays + 1) * 24, 0, 0);
            return sessions.Where(S => (Today - S.CreationTime) <= maxAge).ToArray();
        }

        /// <summary>
        /// Collects all tags in a set of sessions.
        /// </summary>
        public static IEnumerable<string> CollectTags(this IEnumerable<ISessionInfo> sessions) {
            HashSet<string> R = new HashSet<string>();
            foreach (var S in sessions) {
                foreach (var tag in S.Tags) {
                    if (!R.Contains(tag)) {
                        R.Add(tag);
                    }
                }
            }
            return R.ToArray();
        }

        /// <summary>
        /// Retrieves the polynomial Ansatz order of the DG field with name
        /// <paramref name="fieldName"/>.
        /// </summary>
        /// <param name="session">
        /// The considered session
        /// </param>
        /// <param name="fieldName">
        /// The name of the DG field.
        /// </param>
        /// <returns>
        /// The polynomial Ansatz order of field <paramref name="fieldName"/>.
        /// </returns>
        public static int GetOrder(this ISessionInfo session, string fieldName) {
            DGField field = session.Timesteps.First().Fields.Find(fieldName);
            if (field == null) {
                throw new Exception(String.Format(
                    "Could not find field {0} in the given session", fieldName));
            }

            return field.Basis.Degree;
        }

        /// <summary>
        /// Computes the total number of DOF for DG field
        /// <paramref name="fieldName"/> for the given
        /// <paramref name="session"/>. Note that the result is only meaningful
        /// if the number of DOF is constant over time.
        /// </summary>
        /// <param name="session">
        /// The considered session
        /// </param>
        /// <param name="fieldName">
        /// The name of the DG field.
        /// </param>
        /// <returns>
        /// The total number of degrees of freedom for DG field
        /// <paramref name="fieldName"/>.
        /// </returns>
        public static int GetDOF(this ISessionInfo session, string fieldName) {
            int order = session.GetOrder(fieldName);

            int dofPerCell = 1;
            int D = session.Timesteps.First().Grid.SpatialDimension;
            int faculty = 1;
            for (int d = 1; d <= D; d++) {
                dofPerCell *= (order + d);
                faculty *= d;
            }
            dofPerCell = dofPerCell / faculty;

            int cellCount = session.Timesteps.First().Grid.NumberOfCells;
            return dofPerCell * cellCount;
        }

        /// <summary>
        /// Retrieves the ordered list of ancestors of the given session. That
        /// is, recursively gathers the sessions from which the given
        /// <paramref name="session"/> was restarted (cf.
        /// <see cref="ISessionInfo.RestartedFrom"/>), and so on. 
        /// </summary>
        /// <param name="session">
        /// The session whose ancestors shall be retrieved.
        /// </param>
        /// <returns>
        /// The ordered list sessions from which the given
        /// <paramref name="session"/> originates.
        /// </returns>
        public static IEnumerable<ISessionInfo> Ancestors(this ISessionInfo session) {
            if (session == null) {
                yield break;
            }

            Guid currentId = session.RestartedFrom;
            while (currentId != Guid.Empty) {
                var matches = session.Database.Sessions.Where(s => s.ID == currentId);
                if (matches.Count() == 0) {
                    Console.WriteLine(
                        "Aborting search because given session depends on session"
                            + " {0} which is not present in the current database.",
                        currentId);
                    yield break;
                } else if (matches.Count() == 1) {
                    ISessionInfo match = matches.Single();
                    yield return match;
                    currentId = match.RestartedFrom;
                } else {
                    throw new Exception(String.Format(
                        "Found multiple sessions with identical ID {0}."
                            + " This should not have happened.",
                        currentId));
                }
            }
        }

        /// <summary>
        /// Returns <paramref name="session"/> with all its
        /// <see cref="Ancestors(ISessionInfo)"/> as an ordered list
        /// </summary>
        /// <param name="session">
        /// The session whose ancestors shall be retrieved.
        /// </param>
        /// <returns>
        /// [session, [session.RestartedFrom, [session.RestartedFrom.RestartedFrom, ...] ] ]
        /// </returns>
        public static IEnumerable<ISessionInfo> AncestorsAndSelf(this ISessionInfo session) {
            if (session == null) {
                yield break;
            }

            yield return session;
            foreach (var ancestor in session.Ancestors()) {
                yield return ancestor;
            }
        }

        /// <summary>
        /// The number of cores/processes that where used for this simulation
        /// </summary>
        public static int NumberOfCores(this ISessionInfo session) {
            return session.ComputeNodeNames.Count();
        }

        /// <summary>
        /// Function to evaluate results of calculations of FixedCylinder, OscillatingCylinder, ParticleInShear and Particle in Gravity
        /// </summary>
        /// <param name="sessions"></param>
        /// <param name="tags">
        /// 1st: Testcase
        /// 2nd: polynomial degree
        /// 3rd: time discretization
        /// </param>
        public static void IBMPaperEvaluateByTag(this IEnumerable<ISessionInfo> sessions, params string[] tags) {

            //var list = sessions.FindByTag(tags[0]);
            var list = sessions.Where(s => tags.All(tag => s.Tags.Select(sessTag => sessTag.ToLowerInvariant()).Contains(tag.ToLowerInvariant())));

            string sessionIDs;

            sessionIDs = "";

            for (int i = 0; i < list.Count(); i++) {
                sessionIDs = string.Concat(sessionIDs, "'" + list.Pick(i).ID.ToString() + "';");
            }

            string sessionDescriptions = "";
            for (int i = 0; i < list.Count(); i++) {
                sessionDescriptions = string.Concat(sessionDescriptions, "'" + list.Pick(i).Name.ToString() + "';");
            }

            /*            for (int i = 0; i < list.Count(); i++)
                        {
                            var tagList = list.Pick(i + 1).Tags;
                            string nameByTag = "";
                            for (int j = 0; j < list.Count(); j++)
                            {
                                nameByTag = string.Concat(nameByTag, tagList.Pick(j + 1).ToString() + "_");
                            }
                            sessionDescriptions = string.Concat(sessionDescriptions, "'" + nameByTag + "';");
                        }*/

            // Which evaluation case
            int k;

            switch (tags[0]) {
                case ("IBMCylinderFlow"):
                    k = 1;
                    break;

                case ("OscillatingCylinder"):
                    k = 2;
                    break;

                case ("ParticleInShearFlow"):
                    k = 3;
                    break;

                case ("ParticleUnderGravity"):
                    k = 4;
                    break;

                default:
                    throw new NotImplementedException("This testcase can not be evaluated, please check if you've got the correct tag");

            }

            string path = list.First().Database.Path;

            Console.WriteLine("Calling MATLAB...");

            using (BatchmodeConnector bmc = new BatchmodeConnector(path)) {
                bmc.Cmd("Guids = {0} {1} {2}", "[", sessionIDs, "]");
                bmc.Cmd("Descriptions = {0} {1} {2}", "{", sessionDescriptions, "}");
                bmc.Cmd("IBMPaperPostProcessing(Guids,Descriptions,'" + path + "'," + k + ")");
                bmc.Execute(true);
            }

            Console.WriteLine("...Evaluation done");
        }

        /*
        /// <summary>
        /// Calls EvaluatePerformance and plots the DataSets.
        /// </summary>
        /// <param name="sessions"> List of sessions of the same problem but different MPIs </param>
        /// <param name="methods"> Array of methods to be evaluated. If methods == null, the 10 most expensive methods will be taken. </param>
        /// <param name="exclusive"> Boolean that defines if exclusive or inclusive times will be calculated. Methods will still be chosen by exclusive times. </param>
        /// <param name="solver"> String that indicates the solver. Up to now only implemented for IBM_Solver and CNS. </param>
        public static void EvaluatePerformanceAndPlot(this IEnumerable<ISessionInfo> sessions, string[] methods = null, bool exclusive = true, string solver = "IBM_Solver", bool weakScaling = false)
        {
            Plot2Ddata[] data = sessions.EvaluatePerformance(methods,exclusive);
            int numberDataSets = data.Length;
            int numberSessions = sessions.Count();

            // Plotting of all methods' execution times over processors and their ideal curves using Gnuplot
            for (int i = 0; i < numberDataSets/2; i++) {
                Gnuplot gp = new Gnuplot();
                gp.SetMultiplot(1, 2);
                gp.SetSubPlot(0, 0);
                gp.SetXLabel("Processors");
                if (exclusive)
                {
                    gp.SetYLabel("Exlusive times [s]");
                } else
                {
                    gp.SetYLabel("Inclusive times [s]");
                }
                gp.Cmd("set terminal wxt noraise");
                gp.Cmd("set grid xtics ytics");

                int lineColor = 0;
                foreach (var group in data[i].dataGroups)
                {
                    gp.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                        new PlotFormat(lineColor: ((LineColors)(++lineColor)), pointType: ((PointTypes)4), pointSize: 1.5, Style: Styles.LinesPoints));
                }
                gp.WriteDeferredPlotCommands();
                gp.SetSubPlot(0, 1);
                gp.SetXLabel("Processors");
                gp.SetYLabel("Speedup");
                gp.Cmd("set terminal wxt noraise");
                gp.Cmd("set grid xtics ytics");

                lineColor = 0;
                foreach (var group in data[i+numberDataSets/2].dataGroups)
                {
                    gp.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                        new PlotFormat(lineColor: ((LineColors)(++lineColor)), pointType: ((PointTypes)4), pointSize: 1.5, Style: Styles.LinesPoints));
                }
                gp.WriteDeferredPlotCommands();
                gp.Execute();
            }
        }
        */

        /// <summary>
        /// Calculates performance times from profiling_bins for each session for specified methods. Writes out a table of the most expensive and (of those) worst scaling functions. 
        /// Returns data of convergence and speedup for each method over number of MPIs
        /// </summary>
        /// <param name="sessions"> List of sessions of the same problem but different MPIs </param>
        /// <param name="methods"> Array of methods to be evaluated. If methods == null, the 10 most expensive methods will be taken. </param>
        /// <param name="exclusive"> Boolean that defines if exclusive or inclusive times will be calculated. Methods will still be chosen by exclusive times. </param>
        /// <returns>
        /// Returns an array of DataSets, where the first half contains the convergence data for every method and the second half the speedup data.
        /// </returns>
        public static Plot2Ddata[] EvaluatePerformance(this IEnumerable<ISessionInfo> sessions, string[] methods = null, bool exclusive = true/*, string solver = "IBM_Solver",*/, bool weakScaling = false) {
        // <param name="solver"> String that indicates the solver. Up to now only implemented for IBM_Solver and CNS. </param>
            string path = sessions.Pick(0).Database.Path;
            string mainMethod = "*.RunSolverOneStep"; // use wildcard!
            //switch (solver) {
            //    case "IBM_Solver":
            //        mainMethod = "BoSSS.Application.IBM_Solver.IBM_SolverMain.RunSolverOneStep";
            //        break;
            //    case "CNS":
            //        mainMethod = "CNS.Program`1.RunSolverOneStep";
            //        break;
            //    default:
            //        throw new ApplicationException("Main method not defined for this solver yet");
            //}



            // Change maxNumberMethods to how many methods you want considered if no  methods specified
            int maxNumberMethods = 8;
            double[] fraction = new double[maxNumberMethods];
            int idx = sessions.IndexOfMax(s => s.ComputeNodeNames.Count());

            var mcr = sessions.Pick(idx).GetProfiling();

            // Find methods if none given
            if (methods == null) {
                
                var findMainMethod = mcr[0].FindChild(mainMethod);
                IOrderedEnumerable<CollectionReport> mostExpensive;
                
                if (findMainMethod != null) {
                    mostExpensive = findMainMethod.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveTimeFractionOfRoot);
                } else {
                    mostExpensive = mcr[0].CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveTimeFractionOfRoot);
                }

                methods = new string[maxNumberMethods];
                for (int i = 0; i < maxNumberMethods; i++) {
                    methods[i] = mostExpensive.Pick(i).Name;
                }
            }
            int numberMethods = methods.Length;
            string[] methodCalls = new string[numberMethods];

            // Initialise variables
            int numberSessions = sessions.Count();
            Plot2Ddata[] data = new Plot2Ddata[2 * numberMethods];
            double[][] times = new double[numberSessions][];
            int[] processors = new int[numberSessions];

            // Iterate over sessions
            for (int i = 0; i < numberSessions; i++) {
                // Get number of processors and save for later
                int fileCount = (from file in Directory.EnumerateFiles(@path + "\\sessions\\" + sessions.Pick(i).ID, "profiling_bin.*", SearchOption.AllDirectories)
                                 select file).Count();
                int numberProcessors = fileCount;
                processors[i] = numberProcessors;

                double[] maxTime = new double[numberMethods];

                // Iterate over MPIs
                for (int j = 0; j < numberProcessors; j++) {
                    MethodCallRecord value;
                    // Iterate over methods
                    for (int k = 0; k < numberMethods; k++) {
                        // Get execution time of current method for current processor
                        double[] tempTime = new double[numberMethods];
                        double[] tempFractions = new double[numberMethods];
                        int occurence = methods.Take(k+1).Where(x => x.Equals(methods[k])).Count();

                        value = mcr[j].FindChild(mainMethod);
                        if (value == null) {
                            value = mcr[j];
                        }
                        if (exclusive) {
                            tempTime[k] = value.FindChildren(methods[k]).OrderByDescending(s => s.TimeExclusive.TotalSeconds).Pick(occurence-1).TimeExclusive.TotalSeconds;
                            if (i == idx) {
                                IEnumerable<MethodCallRecord> calls = value.FindChildren(methods[k]).OrderByDescending(s => s.ExclusiveTimeFractionOfRoot);
                                double maxValue = calls.Pick(occurence-1).ExclusiveTimeFractionOfRoot;
                                int maxIndex = calls.Select(s => s.ExclusiveTimeFractionOfRoot).ToList().IndexOf(maxValue);
                                tempFractions[k] = maxValue;
                                MethodCallRecord correctCall = calls.Pick(maxIndex);
                                IEnumerable<MethodCallRecord> neighbourCalls = calls.Except(correctCall);
                                if (maxValue > fraction[k]) {
                                    methodCalls[k] = getUniqueParentName(correctCall, neighbourCalls) + " (" + occurence + "/" + calls.Count() + ")";
                                }
                                
                            }
                        } else {
                            tempTime[k] = value.FindChildren(methods[k]).OrderByDescending(s => s.TimeSpentInMethod.TotalSeconds).Pick(occurence-1).TimeSpentInMethod.TotalSeconds;
                            if (i == idx) {
                                IEnumerable<MethodCallRecord> calls = value.FindChildren(methods[k]).OrderByDescending(s => s.TimeFractionOfRoot);
                                double maxValue = calls.Pick(occurence-1).TimeFractionOfRoot;
                                int maxIndex = calls.Select(s => s.TimeFractionOfRoot).ToList().IndexOf(maxValue);
                                tempFractions[k] = maxValue;
                                MethodCallRecord correctCall = calls.Pick(maxIndex);
                                IEnumerable<MethodCallRecord> neighbourCalls =calls.Except(correctCall);
                                if (maxValue > fraction[k]) {
                                    methodCalls[k] = getUniqueParentName(correctCall, neighbourCalls) + " (" + occurence + "/" + calls.Count() + ")";
                                }
                            }
                        }
                        // Only save execution time if it is the highest value of all processor times
                        if (tempTime[k] > maxTime[k]) {
                            maxTime[k] = tempTime[k];
                        }
                        if (tempFractions[k] > fraction[k]) {
                            fraction[k] = tempFractions[k];
                        }
                    }
                }
                times[i] = maxTime;
            }
            Array.Sort(processors, times);

            KeyValuePair<string, double>[] methodRegressionPair = new KeyValuePair<string, double>[numberMethods];
            KeyValuePair<string, double>[] methodFractionPair = new KeyValuePair<string, double>[numberMethods];
            KeyValuePair<string, double>[] callsFractionPair = new KeyValuePair<string, double>[numberMethods];
            // Create DataSets and ideal curves
            for (int i = 0; i < numberMethods; i++) {
                // Calculation of ideal curves
                double[] ideal = new double[numberSessions];
                double[] idealSpeedUp = new double[numberSessions];
                double startIdeal = times.Pick(0)[i];
                for (int j = 0; j < numberSessions; j++) {
                    if (weakScaling) {
                        ideal[j] = startIdeal;
                        idealSpeedUp[j] = 0;
                    } else {
                        ideal[j] = Math.Pow(0.5, j) * startIdeal;
                        idealSpeedUp[j] = processors[j];
                    }
                }

                double[] speedUpTimes = new double[numberSessions];

                var timeArray = times.Select(t => t.Pick(i));
                if (weakScaling) {
                    speedUpTimes = timeArray.Select(x =>( x - startIdeal ) / startIdeal).ToArray();
                } else {
                    speedUpTimes = timeArray.Select(x => startIdeal * processors[0] / x).ToArray();
                }

                // Create DataRows for convergence and speedup with actual and ideal curve
                KeyValuePair<string, double[][]>[] dataRowsConvergence = new KeyValuePair<string, double[][]>[2];
                KeyValuePair<string, double[][]>[] dataRowsSpeedup = new KeyValuePair<string, double[][]>[2];
                double[] doubleProcessors = processors.Select(Convert.ToDouble).ToArray();

                dataRowsConvergence[0] = new KeyValuePair<string, double[][]>(methods[i] + " (" + methodCalls[i].Split('(').Last(), new double[][] { doubleProcessors, times.Select(s => s[i]).ToArray() });
                dataRowsConvergence[1] = new KeyValuePair<string, double[][]>("ideal", new double[][] { doubleProcessors, ideal });
                dataRowsSpeedup[0] = new KeyValuePair<string, double[][]>(methods[i] + " (" + methodCalls[i].Split('(').Last(), new double[][] { doubleProcessors, speedUpTimes });
                dataRowsSpeedup[1] = new KeyValuePair<string, double[][]>("ideal", new double[][] { doubleProcessors, idealSpeedUp });

                // Create DataSets from DataRows
                data[i] = new Plot2Ddata(dataRowsConvergence);
                data[i + numberMethods] = new Plot2Ddata(dataRowsSpeedup);
                methodRegressionPair[i] = new KeyValuePair<string, double>(methods[i], Math.Min(data.Skip(numberMethods).Pick(i).Regression().Pick(0).Value, data.Skip(numberMethods).Pick(i).Regression().Pick(1).Value));
                methodFractionPair[i] = new KeyValuePair<string, double>(methods[i], fraction[i]);
                callsFractionPair[i] = new KeyValuePair<string, double>(methodCalls[i], fraction[i]);
            }

            // Use slope of actual speedup curve to sort methods and DataSets by "worst scaling"
            if (weakScaling) {
                methodRegressionPair = methodRegressionPair.OrderByDescending(t => t.Value).ToArray();
            } else {
                methodRegressionPair = methodRegressionPair.OrderBy(t => t.Value).ToArray();
            }
            methodFractionPair = methodFractionPair.OrderByDescending(t => t.Value).ToArray();
            callsFractionPair = callsFractionPair.OrderByDescending(t => t.Value).ToArray();
            double[] regressions = methodRegressionPair.Select(s => s.Value).ToArray();
            double[] regressions2 = regressions;
            string[] methods2 = methodFractionPair.Select(s => s.Key).ToArray();
            double[] fractions2 = methodFractionPair.Select(s => s.Value).ToArray();
            string[] methodCalls2 = callsFractionPair.Select(s => s.Key).ToArray();
            string[] sortedMethods = methodRegressionPair.Select(s => s.Key).ToArray();

            // Write out the most expensive functions and the worst scaling functions
            Console.WriteLine("\n Most expensive functions");
            Console.WriteLine("============================");
            for (int i = 0; i < numberMethods; i++) {
                Console.WriteLine("Rank " + i + ": " + methods2[i] + " (" + methodCalls2[i].Split('(').Last());
                Console.WriteLine("\t Time fraction of root: " + fractions2[i].ToString("p3") + "\t in " + methodCalls2[i].Split('(').First());
            }
            if (weakScaling) {
                Console.WriteLine("\n ======= WEAK SCALING ========");
            } else {
                Console.WriteLine("\n ======= STRONG SCALING ========");
            }
            Console.WriteLine("\n Sorted by worst scaling");
            Console.WriteLine("============================");
            for (int i = 0; i < numberMethods; i++) {
                Console.WriteLine("Rank " + i + ": " + sortedMethods[i]);
                Console.WriteLine("\t speedup slope: " + regressions[i].ToString("N3"));
            }

            return data;
        }



        private static string getUniqueParentName(MethodCallRecord correctCall, IEnumerable<MethodCallRecord> neighbourCalls) {
            while (neighbourCalls.Count() > 0) {
                foreach (MethodCallRecord m in neighbourCalls) {
                    if (m.ParrentCall.Name != correctCall.ParrentCall.Name) {
                        neighbourCalls = neighbourCalls.Except(m);
                    }
                }
                correctCall = correctCall.ParrentCall;
                if (neighbourCalls != null) {
                    neighbourCalls = neighbourCalls.Select(c => c.ParrentCall);
                }
            }
            return correctCall.Name;
        }

        /// <summary>
        /// Plots circularity and rise velocity over time if a  "BenchmarkQuantities_RisingBubble.txt" exists.
        /// </summary>
        /// <param name="sess"></param> List of sessions to be evaluated
        public static void EvalRisingBubble(this IEnumerable<ISessionInfo> sess) {
            int numberSessions = sess.Count();
            double[][] times = new double[numberSessions][];
            double[][] circularities = new double[numberSessions][];
            double[][] riseVelocities = new double[numberSessions][];

            // Read all data
            for (int j = 0; j < numberSessions; j++) {
                string path = @sess.Pick(j).Database.Path + "\\sessions\\" + sess.Pick(j).ID + "\\BenchmarkQuantities_RisingBubble.txt";
                string[] lines = File.ReadAllLines(path);
                double[] time = new double[lines.Length - 1];
                double[] circularity = new double[lines.Length - 1];
                double[] riseVelocity = new double[lines.Length - 1];

                for (int i = 0; i < lines.Length - 1; i++) {
                    time[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                    circularity[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[4]);
                    riseVelocity[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[5]);
                }
                times[j] = time;
                circularities[j] = circularity;
                riseVelocities[j] = riseVelocity;

            }
            // Build DataSets
            KeyValuePair<string, double[][]>[] dataRowsCircularity = new KeyValuePair<string, double[][]>[numberSessions];
            KeyValuePair<string, double[][]>[] dataRowsRiseVelocity = new KeyValuePair<string, double[][]>[numberSessions];
            for (int i = 0; i < numberSessions; i++) {
                dataRowsCircularity[i] = new KeyValuePair<string, double[][]>(sess.Pick(i).Name, new double[][] { times[i], circularities[i] });
                dataRowsRiseVelocity[i] = new KeyValuePair<string, double[][]>(sess.Pick(i).Name, new double[][] { times[i], riseVelocities[i] });
            }
            Plot2Ddata Time_Circularity = new Plot2Ddata(dataRowsCircularity);
            Plot2Ddata Time_riseVelocity = new Plot2Ddata(dataRowsRiseVelocity);

            // Plot circularity
            int lineColor = 0;
            PlotFormat format = new PlotFormat(lineColor: ((LineColors)(++lineColor)));
            Gnuplot gp = new Gnuplot(baseLineFormat: format);
            gp.SetXLabel("Time");
            gp.SetYLabel("Circularity");
            gp.Cmd("set grid xtics ytics");
            foreach (var group in Time_Circularity.dataGroups) {
                gp.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                    new PlotFormat(lineColor: ((LineColors)(++lineColor))));
            }
            gp.WriteDeferredPlotCommands();
            gp.Execute();

            // Plot rise velocity
            lineColor = 0;
            Gnuplot gp2 = new Gnuplot(baseLineFormat: format);
            gp2.SetXLabel("Time");
            gp2.SetYLabel("Rise velocity");
            gp2.Cmd("set grid xtics ytics");
            foreach (var group in Time_riseVelocity.dataGroups) {
                gp2.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                    new PlotFormat(lineColor: ((LineColors)(++lineColor))));
            }
            gp2.WriteDeferredPlotCommands();
            gp2.Execute();
        }
    }
}
