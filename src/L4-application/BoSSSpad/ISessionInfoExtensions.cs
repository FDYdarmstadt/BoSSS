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
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Xml;
using MathNet.Numerics.Interpolation;
using static BoSSS.Solution.Gnuplot.Plot2Ddata;
using BoSSS.Solution.Statistic;
using BoSSS.Application.XdgPoisson3;
using BoSSS.Solution;

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
        public static string PrintSessionDirectory(this ISessionInfo session) {
            string dir = DatabaseDriver.GetSessionDirectory(session);
            Console.WriteLine(dir);
            return dir;
        }

        /// <summary>
        /// Returns the directory where the files for the selected
        /// <paramref name="session"/> are stored to the console.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        public static string GetSessionDirectory(this ISessionInfo session) {
            string dir = DatabaseDriver.GetSessionDirectory(session);
            return dir;
        }

        /// <summary>
        /// Opens the database directory where the files for the selected
        /// <paramref name="session"/> are stored in the explorer.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <param name="Clientidx">
        /// index into <see cref="BatchProcessorConfig.AllQueues"/>, <see cref="BoSSSshell.ExecutionQueues"/>.
        /// </param>
        public static void OpenDeployDirectory(this ISessionInfo session, int Clientidx = -1) {
            string deployname = session.DeployPath.Split('/', '\\').Last();
            var bpc = BatchProcessorConfig.LoadOrDefault();
            int idx = -1;

            if(Clientidx < 0) {
                for(int i = 0; i < bpc.AllQueues.Length; i++) 
                    Console.WriteLine(i + " : " + bpc.AllQueues[i].ToString());
                Console.WriteLine("Choose Client-Config:");
                string usrinput = Console.ReadLine();
                idx = Convert.ToInt32(usrinput);
            } else {
                idx = Clientidx;
            }

            if(idx > bpc.AllQueues.Length)
                throw new IndexOutOfRangeException();

            var bpclnt = bpc.AllQueues[idx];
            string deploydir = bpclnt.DeploymentBaseDirectory;
            string deploypath = Path.Combine(deploydir, deployname);

            if(!deploydir.IsEmptyOrWhite() && Directory.Exists(deploypath))
                Process.Start(deploypath);
            else
                Console.WriteLine("Attempt failed. Search directory via DeployPath-Attribute");
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
            this ISessionInfo session, string norm = null, int stride = 1, params string[] variables) {

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

            // dbg_launch();

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

                ret.Add(kv.Key, ColumnDoubles);
            }

            return ret;
        }

        /// <summary>
        /// Loads the profiling information for a session
        /// </summary>
        /// <param name="session"></param>
        /// <param name="Ranks">
        /// A selection of MPI ranks; there is a separate profiling log for each MPI rank; if null or empty, all profiling for all MPI ranks are returned.
        /// </param>
        /// <returns>
        /// An array of profiling trees, one for each MPI rank; the index into the returned array corresponds with <paramref name="Ranks"/>.
        /// </returns>
        public static OnlineProfiling[] GetProfiling(this ISessionInfo session, params int[] Ranks) {
            // find
            string sessDir = DatabaseDriver.GetSessionDirectory(session);
            string[] TextFiles = Directory.GetFiles(sessDir, "profiling_bin.*.txt");
            if (TextFiles.Count() <= 0)
                throw new IOException("Unable to find profiling information.");

            // sort according to process rank
            int MPISize = session.ComputeNodeNames.Count;
            if (Ranks == null || Ranks.Length <= 0)
                Ranks = MPISize.ForLoop(rnk => rnk);
            string[] TextFilesSorted = new string[Ranks.Length];
            for (int i = 0; i < Ranks.Length; i++) {
                if (Ranks[i] >= MPISize || Ranks[i] < 0)
                    throw new ArgumentException($"Illegal MPI rank specified (Rank: {Ranks[i]}). MPI size is {MPISize}, ranks are excepted in the range of 0 to {MPISize - 1} (including); Session: {session} @ {session.GetSessionDirectory()}");
                TextFilesSorted[i] = TextFiles.SingleOrDefault(FileName => FileName.Contains("." + Ranks[i] + "."));
                if (TextFilesSorted[i] == null)
                    throw new IOException($"Unable to find profiling file for MPI rank {Ranks[i]}; Session: {session} @ {session.GetSessionDirectory()}");

                //var parts = TextFils[i].Split(new string[] { "profiling_bin.", ".txt" }, StringSplitOptions.RemoveEmptyEntries);
                //if (parts.Length <= 0)
                //    throw new IOException("Unable to determine file rank from path '" + parts[i] + "'.");
                //Ranks[i] = int.Parse(parts.Last());
                                
                //if (Ranks[i] < 0)
                //    throw new IOException("Unable to determine file rank from path '" + parts[i] + "'.");
            }

            
            // load 
            var R = new OnlineProfiling[Ranks.Max() + 1];
            for(int i = 0; i < Ranks.Length; i++) {
                int rnk = Ranks[i];

                var f = TextFilesSorted[i];
                var JSON = File.ReadAllText(f);
                var mcr = OnlineProfiling.Deserialize(JSON);

                if (R[rnk] != null)
                    throw new IOException("It seems profiling info was written more than once for MPI rank " + rnk + ".");

                R[rnk] = mcr;
            }

            // return
            return R;
        }

        /// <summary>
        /// Gets a profiling tree for a given MPI-rank, of a session
        /// </summary>
        /// <param name="session"></param>
        /// <param name="rank"></param>
        /// <returns></returns>
        public static OnlineProfiling GetProfilingOfRank(this ISessionInfo session, int rank) {

            // check if within bounds
            int MPISize = session.ComputeNodeNames.Count;
            if (MPISize < rank) {
                Console.WriteLine("WARNING: the searched rank (" + rank + ") is greater then MPI-size (" + MPISize + ").");
            }

            // check if exists
            string sessDir = DatabaseDriver.GetSessionDirectory(session);
            string pathf = String.Concat(sessDir, @"\profiling_bin.", rank, ".txt");
            string namef = String.Concat("profiling_bin.", rank, ".txt");
            if (!File.Exists(pathf))
                throw new IOException("Unable to locate '" + pathf + "'.");

            // check if unique
            var many = Directory.GetFiles(sessDir, namef);
            if (many.Length > 1)
                throw new ArgumentException("profiling is not unique; got " + many.ToConcatString("", ", ", ";"));

            // load 
            var f = pathf;
            var JSON = File.ReadAllText(f);
            var mcr = OnlineProfiling.Deserialize(JSON);

            return mcr;
        }

        private static void PrintImbalance(Dictionary<string, (double RelInbalance, double Imbalance, int CallCount)> dictImbalances, int printcnt) {
            var mostimbalance = dictImbalances.OrderByDescending(im => im.Value.RelInbalance);
            int i = 1;
            var wrt = Console.Out;
            foreach (var kv in mostimbalance) {
                wrt.Write("#" + i + ": ");
                wrt.WriteLine(string.Format(
                "'{0}': {1} calls, {2:F3}% / {3:0.##E-00} sec. runtime exclusive",
                    kv.Key,
                    kv.Value.Item3,
                    kv.Value.Item1,
                    kv.Value.Item2));
                if (i == printcnt) return;
                i++;
            }
            Console.Out.Flush();
        }

        /// <summary>
        /// Prints the methods with the highest Imbalance of runtime over MPI-ranks,
        /// which reflects load imbalance. Runtime of FuncTraces are taken.
        /// Uses <see cref="ISessionInfoExtensions.GetProfiling"/> internally, which means this can be expensive.
        /// NOTE: this considered no idle time within methods!
        /// </summary>
        /// <param name="SI"></param>
        /// <param name="printcnt"></param>
        public static void PrintTotalImbalance(this ISessionInfo SI, int printcnt = 0) {
            var dictImbalances = MethodCallRecordExtension.GetFuncImbalance(SI.GetProfiling().Select(p => p.RootCall).ToArray());
            PrintImbalance(dictImbalances, printcnt);
        }

        /// <summary>
        /// Prints the methods with the highest Imbalance within MPI blocking methods,
        /// which reflects communication delay.
        /// Uses <see cref="ISessionInfoExtensions.GetProfiling"/> internally, which means this can be expensive.
        /// NOTE: Nonblocking MPI-Methods are not included! Refer to this as a hint rather beeing accurate
        /// If <see cref="PrintTotalImbalance"/> and this show same tendency, it is likley,
        /// that <see cref="PrintTotalImbalance"/> is dominated by communication delays.
        /// </summary>
        /// <param name="SI"></param>
        /// <param name="printcnt"></param>
        public static void PrintMPIImbalance(this ISessionInfo SI, int printcnt = 0) {
            var dictImbalances = MethodCallRecordExtension.GetMPIImbalance(SI.GetProfiling().Select(p => p.RootCall).ToArray());
            PrintImbalance(dictImbalances, printcnt);
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


            using (StreamReader reader = new StreamReader(new FileStream(TextFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))) {
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
        /// true, if the BoSSS app which wrote the session is successfully terminated (w.o. exceptions, etc.)
        /// </summary>
        public static bool SuccessfulTermination(this ISessionInfo session) {
            return !session.Tags.Contains(SessionInfo.NOT_TERMINATED_TAG);
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
        /// <see cref="BoSSS.Solution.Statistic.DGFieldComparison.ComputeErrors(IList{IEnumerable{DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]}, NormType)"/>.
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
        /// <param name="xAxis_Is_hOrDof">
        /// - true: the x-axis (<see cref="Plot2Ddata.XYvalues.Abscissas"/>) is the grid resolution \f$ h \f$
        /// - false: the x-axis (<see cref="Plot2Ddata.XYvalues.Abscissas"/>) is the number of degrees-of-freedom
        /// </param>     
        /// <param name="normType">
        /// H1, L2, etc.
        /// </param>
        /// <returns>
        /// A data set containing information about the grid resolution and the
        /// corresponding errors with respect to the finest corresponding grid,
        /// grouped by the polynomial degree. Obviously, the time-step
        /// associated with the finest grid for each polynomial degree has an
        /// estimated error of zero (by definition) and is thus excluded from
        /// the result.
        /// </returns>
        public static Plot2Ddata ToEstimatedGridConvergenceData(this IEnumerable<ISessionInfo> sessions, string fieldName, bool xAxis_Is_hOrDof = true, NormType normType = NormType.L2_approximate) {
            ISessionInfo[] _session = sessions.ToArray();
            ITimestepInfo[] _timesteps = sessions.Select(s => s.Timesteps.Last()).ToArray();
            //Debugger.Launch();
            return _timesteps.ToEstimatedGridConvergenceData(fieldName, xAxis_Is_hOrDof, normType);
        }

        /// <summary>
        /// Tries to loads the control file of the given
        /// <paramref name="session"/> in the new REPL format. 
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
        /// Filters all sessions which carry the <see cref="SessionInfo.NOT_TERMINATED_TAG"/>-tag.
        /// </summary>
        public static IEnumerable<ISessionInfo> RunningOrCrashed(this IEnumerable<ISessionInfo> sessions) {
            return sessions.Where(S => S.Tags.Contains(SessionInfo.NOT_TERMINATED_TAG)).ToArray();
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
        /// Computes the total number of DOF for DG or XDG field for the first timestep.
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
        /// The total number of degrees of freedom for DG/XDG field
        /// <paramref name="fieldName"/>.
        /// </returns>
        public static int GetDOF(this ISessionInfo session, string fieldName) {

            int cellCount = session.Timesteps.First().Grid.NumberOfCells;
            int order = session.GetOrder(fieldName);
            var targetfield = session.Timesteps.First().Fields.Find(fieldName);

            int dofPerCell = 1;
            int D = session.Timesteps.First().Grid.SpatialDimension;
            int faculty = 1;
            for (int d = 1; d <= D; d++) {
                dofPerCell *= (order + d);
                faculty *= d;
            }
            dofPerCell = dofPerCell / faculty;

            if (targetfield.GetType() == typeof(XDGField)) {
                var phi = session.Timesteps.First().Fields.ElementAt(0);
                var LevSet = new LevelSet(phi.Basis, "LevelSet");
                LevSet.Acc(1.0, phi);
                var LsTrk = new LevelSetTracker((BoSSS.Foundation.Grid.Classic.GridData)phi.GridDat, XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B" }, LevSet);
                int numCC = LsTrk.Regions.GetCutCellMask().Count();
                cellCount += numCC;
            }
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

        
        /// <summary>
        /// Calls EvaluatePerformance and plots the DataSets.
        /// </summary>
        /// <param name="sessions"> List of sessions of the same problem but different MPIs </param>
        /// <param name="methods"> Array of methods to be evaluated. If methods == null, the 10 most expensive methods will be taken. </param>
        /// <param name="exclusive"> Boolean that defines if exclusive or inclusive times will be calculated. Methods will still be chosen by exclusive times. </param>
        /// <param name="solver"> String that indicates the solver. Up to now only implemented for IBM_Solver and CNS. </param>
        /// <param name="weakScaling"></param>
        public static void EvaluatePerformanceAndPlot(this IEnumerable<ISessionInfo> sessions, string[] methods = null, bool exclusive = true, string solver = "IBM_Solver", bool weakScaling = false)
        {
            Plot2Ddata[] data = sessions.EvaluatePerformance(methods,exclusive,weakScaling);
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
        

        /// <summary>
        /// Calculates performance times from profiling_bins for each session for specified methods. Writes out a table of the most expensive and (of those) worst scaling functions. 
        /// Returns data of convergence and speedup for each method over number of MPIs
        /// </summary>
        /// <param name="sessions"> List of sessions of the same problem but different MPIs </param>
        /// <param name="methods"> Array of methods to be evaluated. If methods == null, the 10 most expensive methods will be taken. </param>
        /// <param name="exclusive"> Boolean that defines if exclusive or inclusive times will be calculated. Methods will still be chosen by exclusive times. </param>
        /// <param name="weakScaling"></param>
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

            var profiling = sessions.Pick(idx).GetProfiling();

            // Find methods if none given
            if (methods == null) {
                
                var findMainMethod = profiling[0].RootCall.FindChild(mainMethod);
                IOrderedEnumerable<CollectionReport> mostExpensive;
                
                if (findMainMethod != null) {
                    mostExpensive = findMainMethod.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveTimeFractionOfRoot);
                } else {
                    mostExpensive = profiling[0].RootCall.CompleteCollectiveReport().OrderByDescending(cr => cr.ExclusiveTimeFractionOfRoot);
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

                // was missing
                profiling = sessions.Pick(i).GetProfiling();

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

                        value = profiling[j].RootCall.FindChild(mainMethod);
                        if (value == null) {
                            value = profiling[j].RootCall;
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
        /// imports the specified log file data 
        /// </summary>
        /// <param name="sess"> List of sessions to be evaluated </param>
        /// <param name="logName"> which log values to be evaluated </param>
        /// <param name="evalName"></param>
        /// <param name="keyName"></param>
        /// <returns></returns>
        public static List<Plot2Ddata> ReadLogDataForXNSE(this List<ISessionInfo> sess, string logName, string evalName = null, string keyName = null) {


            string[] values;
            switch (logName) {
                case Application.XNSE_Solver.PhysicalBasedTestcases.WaveLikeLogging.LogfileName: {
                        values = new string[] { "#timestep", "time", "magnitude", "real", "imaginary" };
                        break;
                    }
                case Application.XNSE_Solver.PhysicalBasedTestcases.Dropletlike.LogfileName: {
                        values = new string[] { "#timestep", "time", "semi axis x", "semi axis y", "area", "perimeter" };
                        break;
                    }
                case Application.XNSE_Solver.PhysicalBasedTestcases.RisingBubble2DBenchmarkQuantities.LogfileName: {
                        values = new string[] { "#timestep", "time", "area", "center of mass - x", "center of mass - y", "circularity", "rise velocity" };
                        break;
                    }
                case Application.XNSE_Solver.PhysicalBasedTestcases.MovingContactLineLogging.LogfileName: {
                        values = new string[] { "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle" };
                        break;
                    }
                case Application.XNSFE_Solver.PhysicalBasedTestcases.EvaporationLogging.LogfileName: {
                    values = new string[] { "#timestep", "time", "interfacePosition", "meanInterfaceVelocity", "meanMassFlux" };
                    break;
                }
                case Application.XNSFE_Solver.PhysicalBasedTestcases.StefanProblemBenchmarkQuantities.LogfileName: {
                    values = new string[] { "#timestep", "time", "interface-x-pos-min", "interface-x-pos-max", "mass-vapor", "mass-liquid", "massflux-interface", "massflux-outlet" };
                    break;
                }
                case Application.XNSFE_Solver.PhysicalBasedTestcases.MassfluxLogging.LogfileName: {
                    values = new string[] { "#timestep", "time", "mass-liq", "mass-vap", "mass-total", "masschange-evap", "masschange-vapor", "masschange-liquid", "masschange-total", "interface length" };
                    break;
                }
                default:
                    throw new ArgumentException("No specified LogFormat");
            }


            List<Plot2Ddata> plotData = new List<Plot2Ddata>();

            int numberSessions = sess.Count();
            int numberValues = values.Count();      
            for (int vIdx = 2; vIdx < numberValues; vIdx++) {       

                double[][] times = new double[numberSessions][];
                double[][] valueDatas = new double[numberSessions][];

                // Read all data
                for (int j = 0; j < numberSessions; j++) {
                    string path = Path.Combine(sess.Pick(j).Database.Path, "sessions", sess.Pick(j).ID.ToString(), logName + ".txt");
                    string[] lines = File.ReadAllLines(path);

                    if (sess.Pick(j).RestartedFrom == Guid.Empty) { 
                   
                        double[] time = new double[lines.Length - 1];
                        double[] valueData = new double[lines.Length - 1];

                        for (int i = 0; i < lines.Length - 1; i++) {
                            time[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                            valueData[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[vIdx]);
                        }
                        times[j] = time;
                        valueDatas[j] = valueData;

                    } else {

                        string pathR = @sess.Pick(j).Database.Path + "\\sessions\\" + sess.Pick(j).RestartedFrom + logName;
                        string[] linesR = File.ReadAllLines(pathR);

                        int len = (lines.Length - 1) + (linesR.Length - 1);
                        double[] time = new double[len];
                        double[] valueData = new double[len];
                        int iL = 0;
                        for (int i = 0; i < linesR.Length - 1; i++) {
                            time[iL] = Convert.ToDouble(linesR[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                            valueData[iL] = Convert.ToDouble(linesR[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[vIdx]);
                            iL++;
                        }
                        for (int i = 0; i < lines.Length - 1; i++) {
                            time[iL] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                            valueData[iL] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[vIdx]);
                            iL++;
                        }

                        // remove doubled time steps 
                        List<double> rTime = new List<double>();
                        List<double> rValDat = new List<double>();
                        rTime.Add(time[len-1]);
                        rValDat.Add(valueData[len-1]);
                        for(int i = len-2; i >= 0; i--) {
                            if (time[i] < rTime.Last()) {
                                rTime.Add(time[i]);
                                rValDat.Add(valueData[i]);
                            }
                        }
                        rTime.Reverse();
                        rValDat.Reverse();

                        times[j] = rTime.ToArray();
                        valueDatas[j] = rValDat.ToArray();

                    }
                }

                // Build DataSet
                KeyValuePair<string, double[][]>[] dataRowsValue = new KeyValuePair<string, double[][]>[numberSessions];
                for (int i = 0; i < numberSessions; i++) {
                    string sessName;
                    if (evalName == null || keyName == null)
                        sessName = (sess.Pick(i).Name).Replace("_", "-");
                    else
                        sessName = evalName + (Convert.ToDouble(sess.Pick(i).KeysAndQueries[keyName])).ToString();

                    dataRowsValue[i] = new KeyValuePair<string, double[][]>(sessName, new double[][] { times[i], valueDatas[i] });
                }
                Console.WriteLine("Element at {0}: time vs {1}", vIdx - 2, values[vIdx]);
                plotData.Add(new Plot2Ddata(dataRowsValue));
            }

            return plotData;

        }

        /// <summary>
        /// special purpose method, most likely legacy stuff
        /// </summary>
        public static List<Plot2Ddata>[] ReadLogDataForMovingContactLine(this IEnumerable<ISessionInfo> sess) {

            string[] values = new string[] { "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle" };

            // check number of contact lines
            string path = @sess.Pick(0).Database.Path + "\\sessions\\" + sess.Pick(0).ID + "\\ContactAngle.txt";
            string[] lines = File.ReadAllLines(path);
            int numCL = 0;
            for (int i = 1; i <= 4; i++) {       // max number of contact lines should be 4
                int ts = (int)Convert.ToDouble(lines[i].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[0]);
                if (ts == 0)
                    numCL++;
            }

            Console.WriteLine("number of contact lines: {0}", numCL);

            List<Plot2Ddata>[] plotDataCL = new List<Plot2Ddata>[numCL];

            int numberSessions = sess.Count();
            int numberValues = values.Count();
            for (int c = 0; c < numCL; c++) {

                List<Plot2Ddata> plotData = new List<Plot2Ddata>();

                for (int vIdx = 2; vIdx < numberValues; vIdx++) {       // considered are only the values over time

                    double[][] times = new double[numberSessions][];
                    double[][] valueDatas = new double[numberSessions][];

                    // Read all data
                    for (int j = 0; j < numberSessions; j++) {
                        path = @sess.Pick(j).Database.Path + "\\sessions\\" + sess.Pick(j).ID + "\\ContactAngle.txt";
                        lines = File.ReadAllLines(path);
                        double[] time = new double[(lines.Length - 1)/numCL];
                        double[] valueData = new double[(lines.Length - 1)/numCL];

                        int ind = 0;
                        for (int i = c; i < lines.Length - 1; i+=numCL) {
                            time[ind] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                            valueData[ind] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[vIdx]);
                            ind++;
                        }
                        times[j] = time;
                        valueDatas[j] = valueData;
                    }

                    // Build DataSet
                    KeyValuePair<string, double[][]>[] dataRowsValue = new KeyValuePair<string, double[][]>[numberSessions];
                    for (int i = 0; i < numberSessions; i++) {
                        dataRowsValue[i] = new KeyValuePair<string, double[][]>(sess.Pick(i).Name, new double[][] { times[i], valueDatas[i] });
                    }
                    if(c == 0)
                        Console.WriteLine("Element at {0}: time vs {1}", vIdx - 2, values[vIdx]);

                    plotData.Add(new Plot2Ddata(dataRowsValue));
                }

                plotDataCL[c] = plotData;
            }

            return plotDataCL;

        }

        /// <summary>
        /// 
        /// </summary>
        public static List<Plot2Ddata> LogDataToConvergenceData(List<Plot2Ddata> LogData, double[] abscissa, double[] _refAbs = null, double[] _refVal = null) {

            List<Plot2Ddata> convData = new List<Plot2Ddata>();

            foreach (var p2d in LogData) {

                int numSess = p2d.dataGroups.Length;
                int i0;

                double[] refAbs;
                double[] refVal;
                if (_refAbs != null && _refVal != null) {
                    if (abscissa.Length != numSess)
                        throw new ArgumentException("wrong length of abscissa");
                    refAbs = _refAbs;
                    refVal = _refVal;
                    i0 = 0;
                } else {
                    if (abscissa.Length != numSess - 1)
                        throw new ArgumentException("wrong length of abscissa");
                    // set reference data [0] (log data in ascending order)
                    refAbs = p2d.dataGroups[0].Abscissas;
                    refVal = p2d.dataGroups[0].Values;
                    i0 = 1;
                }
                int numRefVal = refVal.Length;

                double[] l1Norm = new double[numSess - i0];
                double[] l2Norm = new double[numSess - i0];
                double[] linfNorm = new double[numSess - i0];
                KeyValuePair<string, double[][]>[] dataRowsValue = new KeyValuePair<string, double[][]>[3];
                //foreach (var datgrp in p2d.dataGroups.Skip(1)) {
                for (int i = i0; i < numSess; i++) {

                    double[] abs = p2d.dataGroups[i].Abscissas;
                    double[] val = p2d.dataGroups[i].Values;
                    int numVal = val.Length;

                    double[] diff = new double[numRefVal];
                    if (numVal != numRefVal) {
                        if (numRefVal < numVal)
                            throw new ArgumentException("reference data should have at least the same length as comparison data");
                        // interpolate solution
                        //LinearSplineInterpolation LinSpline = new LinearSplineInterpolation();
                        //LinSpline.Initialize(abs, val);
                        LinearSpline LinSpline = LinearSpline.InterpolateSorted(abs, val);

                        val = new double[numRefVal];
                        for (int p = 0; p < numRefVal; p++) {
                            val[p] = LinSpline.Interpolate(refAbs[p]);
                        }
                    }
                    diff = refVal.Zip(val, (r, v) => Math.Abs(v - r)).ToArray();

                    l1Norm[i - i0] = diff.Sum() / refVal.Sum();
                    l2Norm[i - i0] = (diff.L2NormPow2() / refVal.L2NormPow2()).Sqrt();
                    linfNorm[i - i0] = diff.Max() / refVal.Max();

                    //double[] diff = new double[numVal];
                    //double[] refValCom = new double[numVal];
                    //if (numVal != numRefVal) {
                    //    if (numRefVal < numVal)
                    //        throw new ArgumentException("reference data should have at least the same length as comparison data");
                    //    // interpolate solution
                    //    LinearSplineInterpolation LinSpline = new LinearSplineInterpolation();
                    //    LinSpline.Initialize(refAbs, refVal);

                    //    for (int p = 0; p < numVal; p++) {
                    //        refValCom[p] = LinSpline.Interpolate(abs[p]);
                    //    }
                    //}
                    //diff = refValCom.Zip(val, (r, v) => Math.Abs(v - r)).ToArray();

                    //l1Norm[i - i0] = diff.Sum() / refValCom.Sum();
                    //l2Norm[i - i0] = (diff.L2NormPow2() / refValCom.L2NormPow2()).Sqrt();
                    //linfNorm[i - i0] = diff.Max() / refValCom.Max();
                }

                dataRowsValue[0] = new KeyValuePair<string, double[][]>("l_1 error norm", new double[][] { abscissa, l1Norm });
                dataRowsValue[1] = new KeyValuePair<string, double[][]>("l_2 error norm", new double[][] { abscissa, l2Norm });
                dataRowsValue[2] = new KeyValuePair<string, double[][]>("l_{inf} error norm", new double[][] { abscissa, linfNorm });

                convData.Add(new Plot2Ddata(dataRowsValue).WithLogX().WithLogY());

            }

            return convData;
        }

        /// <summary>
        /// ???
        /// </summary>
        public static Tuple<List<double>, List<double>> ComputeOscillationProperties(XYvalues waveData, bool twoPiPeriodic) {

            List<double> maxValues = new List<double>();
            List<double> periods = new List<double>();

            double[] times = waveData.Abscissas;
            double[] values = waveData.Values;

            double max = values[0];
            maxValues.Add(max);
            periods.Add(times[0]);

            double min = max;
            bool searchMin = true;
            for (int i = 0; i < times.Length; i++) {
                double val = values[i];
                if (searchMin) {
                    if (val < min) {
                        min = val;
                    } else {
                        max = min;
                        searchMin = false;
                    }
                }
                if (!searchMin) {
                    if (val > max) {
                        max = val;
                    } else {
                        maxValues.Add(max);
                        periods.Add(times[i]);
                        min = max;
                        searchMin = true;
                    }
                }
            }

            //return new Tuple<List<double>, List<double>>(periods, maxValues);

            List<double> frequencies = new List<double>();
            List<double> dampingRates = new List<double>();
            for (int i = 0; i < maxValues.Count()-1; i++) {
                double period = periods[i + 1] - periods[i];
                if (period <= 1.0e-8)
                    continue;
                double freq = twoPiPeriodic ? 1.0 / period : 1.0 / (2.0 * period);
                double damp = Math.Log(maxValues[i + 1] / maxValues[i]) / period;
                frequencies.Add(freq);
                dampingRates.Add(damp);
            }

            return new Tuple<List<double>, List<double>>(frequencies, dampingRates);

        }

        /// <summary>
        /// ???
        /// </summary>
        public static void CheckForEnergyLogging(this IEnumerable<ISessionInfo> pSessions) {

            int numberSessions = pSessions.Count();
            for (int j = 0; j < numberSessions; j++) {
                Console.WriteLine("Session: {0}", pSessions.Pick(j).ID);
                string path = @pSessions.Pick(j).Database.Path + "\\sessions\\" + pSessions.Pick(j).ID + "\\Energy.txt";

                string header;
                try {
                    header = File.ReadAllLines(path)[0];
                } catch {
                    Console.WriteLine("no energy file available");
                    continue;
                }

                int energycount = header.Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries).Length;

                for (int i = 0; i < energycount; i++) {
                    string energyName = header.Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[i];
                    Console.WriteLine("Element at {0}: {1}", i, energyName);
                }
            }

        }


        /// <summary>
        /// Check for condition number loggings 
        /// </summary>
        /// <param name="pSessions"></param>
        /// <returns>An array of dictionaries, where each dictionary represents a session, with keys as column names (string) and values as a list of doubles.</returns>
        public static Dictionary<Guid, Dictionary<string, List<double>>> CheckForCondLogging(this IEnumerable<ISessionInfo> pSessions)
        {
            string[] allColumnNames = new string[] { ""};
            int numberSessions = pSessions.Count();

            //the so-called database the first key: session Id, second key: column name, value: list of entries
            Dictionary<Guid, Dictionary<string, List<double>>> logsForAllSessions = new Dictionary<Guid, Dictionary<string, List<double>>>(numberSessions);
            
            for (int j = 0; j < numberSessions; j++){
                ISessionInfo currentSession = pSessions.Pick(j);
                //Initiate the "database", it should suffice the need
                Dictionary<string, List<double>> logs = new Dictionary<string, List<double>>();
                logsForAllSessions[currentSession.ID] = logs;

                Console.WriteLine("Session: {0}", currentSession.ID);
                string path = @currentSession.Database.Path + "\\sessions\\" + currentSession.ID + "\\CondNumbers.txt";
                string[] lines;
                string header;

                try{ //reading
                    lines = File.ReadAllLines(path);
                    header = lines[0];
                } catch{
                    Console.WriteLine("no logging file available");
                    continue;
                }

                //Get column names
                var columnsNames = header.Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries);
                if (!allColumnNames.SequenceEqual(columnsNames)) {
                    allColumnNames = columnsNames;
                    allColumnNames.ForEach(c => Console.Write(c + ", "));
                    Console.WriteLine("");
                }
                //Initiate a list for each column and add to the database
                for (int i = 0; i < columnsNames.Length; i++){
                    string columnName = columnsNames[i];
                    var list = new List<double>();
                    logs.Add(columnName, list);  
                }

                // loop for each line in the log
                for (int k = 1; k < lines.Length; k++){
                    for (int i = 0; i < columnsNames.Length; i++){
                        string currentColumn = columnsNames[i];
                        string value = lines[k].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[i];
                        double valueDouble;

                        //Check if it is NaN or Infinity
                        if (Double.TryParse(value, out valueDouble)){
                            logs[currentColumn].Add(valueDouble);
                        } else {                    
                            Console.WriteLine($"The value '{value}' could not be converted to a double for line {k} at column {i} in session {currentSession.Name}.");
                        }
                    }
                }

            }
            return logsForAllSessions;
        }

        /// <summary>
        /// Plots selected energy over time if an  "Energy.txt" exists.
        /// </summary>
        /// <param name="pSessions">List of sessions to be evaluated</param> 
        /// <param name="energytype"> Energytypes to be plotted, can be partial</param>
        /// <param name="singlePlot"></param>
        public static Plot2Ddata[] EvalEnergy(this IEnumerable<ISessionInfo> pSessions, string[] energytype, bool singlePlot) {
            int numberSessions = pSessions.Count();
            //List<Gnuplot> Plots = new List<Gnuplot>();

            Plot2Ddata[] Time_Energies = new Plot2Ddata[energytype.Length];

            // Cycle over all given energytypes
            for (int g = 0; g < energytype.Length; g++) {
                // Create evaluation variables
                int validSessions = 0;
                int energypos = -1;
                List<double[]> times = new List<double[]>();
                List<double[]> energies = new List<double[]>();

                // Read the data from the sessions
                for (int j = 0; j < numberSessions; j++) {
                    string path = @pSessions.Pick(j).Database.Path + "\\sessions\\" + pSessions.Pick(j).ID + "\\Energy.txt";
                    string[] lines = File.ReadAllLines(path);

                    int energycount = lines[0].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries).Length;

                    List<string> energynames = new List<string>();
                    for (int i = 0; i < energycount; i++) {
                        energynames.Add(lines[0].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[i]);
                    }
                    energypos = energynames.FindIndex(s => s.Contains(energytype[g]));

                    double[] time = new double[lines.Length - 1];
                    double[] energy = new double[lines.Length - 1];

                    if (energypos != -1) {
                        energytype[g] = energynames[energypos];
                        //Console.WriteLine("Found " + energytype[g] + " in " + pSessions.Pick(j).Name);
                        for (int k = 0; k < lines.Length - 1; k++) {
                            time[k] = Convert.ToDouble(lines[k + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                            energy[k] = Convert.ToDouble(lines[k + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[energypos]);
                        }
                        validSessions++;
                    } else {
                        Console.WriteLine(pSessions.Pick(j).Name + " does not contain an Energytype with key '" + energytype[g] + "', maybe try a different spelling/capitalization");
                    }
                    times.Add(time);
                    energies.Add(energy);
                }

                // Build DataSets
                KeyValuePair<string, double[][]>[] dataRowsEnergy = new KeyValuePair<string, double[][]>[validSessions];
                for (int i = 0; i < validSessions; i++) {
                    dataRowsEnergy[i] = new KeyValuePair<string, double[][]>(pSessions.Pick(i).Name, new double[][] { times[i], energies[i] });
                }
                Time_Energies[g] = new Plot2Ddata(dataRowsEnergy);

            }

            int numplt = (singlePlot) ? 1 : energytype.Length;
            for (int g = 0; g < numplt; g++) {
                // Plot energy
                int e0 = (singlePlot) ? 0 : g;
                int eL = (singlePlot) ? energytype.Length : g + 1;
                for (int e = e0; e < eL; e++) {
                    int lineColor = 0;
                    PlotFormat format = new PlotFormat(lineColor: ((LineColors)(++lineColor)));
                    Gnuplot gp = new Gnuplot(baseLineFormat: format);
                    gp.SetXLabel("Time");
                    gp.SetYLabel(energytype[g]);
                    gp.Cmd("set grid xtics ytics");
                    foreach (var group in Time_Energies[e].dataGroups) {
                        gp.PlotXY(group.Abscissas, group.Values, group.Name.Split().Last(),
                            new PlotFormat(lineColor: ((LineColors)(++lineColor))));
                    }
                    gp.WriteDeferredPlotCommands();
                    gp.Execute();
                }
            }

            return Time_Energies;

        }

        /// <summary>
        /// Standard output of some session
        /// </summary>
        public static string GetStdout(this ISessionInfo sess, int rank = 0) {
            var StdoutFile = sess.FilesInSessionDir("stdout." + rank + ".txt").FirstOrDefault();
            if(StdoutFile == null || !File.Exists(StdoutFile)) {
                Console.Error.WriteLine("Missing stdout file for session: " + sess);
                return "";
            }

            using (FileStream stream = File.Open(StdoutFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) {
                using (StreamReader reader = new StreamReader(stream)) {
                    string stdout = reader.ReadToEnd();
                    return stdout;
                }
            }
        }

        /// <summary>
        /// Standard output of some session
        /// </summary>
        public static string GetStderr(this ISessionInfo sess, int rank = 0) {
            var StdoutFile = sess.FilesInSessionDir("stderr." + rank + ".txt").FirstOrDefault();
            if (StdoutFile == null || !File.Exists(StdoutFile)) {
                Console.Error.WriteLine("Missing stderr file for session: " + sess);
                return "";
            }

            using (FileStream stream = File.Open(StdoutFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) {
                using (StreamReader reader = new StreamReader(stream)) {
                    string stdout = reader.ReadToEnd();
                    return stdout;
                }
            }
        }

        /// <summary>
        /// total memory (aka. sum) over all MPI ranks over time
        /// </summary>
        static public Plot2Ddata GetMPItotalMemory(this ISessionInfo sess, int LineOffset = 0, int MaxLines = 10000) {
            var ana = new SessionMemtrace(new DirectoryInfo(sess.GetSessionDirectory()), LineOffset, MaxLines);
            var ret = ana.GetMPItotalMemory();

            ret.Title = "Total memory of session " + sess;


            return ret;        
        }

        /// <summary>
        /// minimum, average and maximum memory allocations over all MPI ranks over time
        /// </summary>
        static public Plot2Ddata GetMinAvgMaxMemory(this ISessionInfo sess, int LineOffset = 0, int MaxLines = 10000) {
            var ana = new SessionMemtrace(new DirectoryInfo(sess.GetSessionDirectory()), LineOffset, MaxLines);
            var ret = ana.GetMinAvgMaxMemPlot();
            ret.Title = "Memory of session " + sess.ID;
            return ret;
        }

        /// <summary>
        /// Returns the memory instrumentation for a session (if available),
        /// combined from files `memory.mpi_rank.txt` in the session directory
        /// </summary>
        public static SessionMemtrace GetMemtrace(this ISessionInfo sess, int LineOffset = 0, int MaxLines = 10000) {
            var ret = new SessionMemtrace(new DirectoryInfo(sess.GetSessionDirectory()), LineOffset, MaxLines);
            return ret;
        }


        /// <summary>
        /// Reports the largest memory-allocating routines in descending order
        /// </summary>
        static public (int TimelineIndex, double Megs, string Name)[] ReportLargestAllocators(this ISessionInfo sess, int LineOffset = 0, int MaxLines = 10000) {
            var ana = new SessionMemtrace(new DirectoryInfo(sess.GetSessionDirectory()), LineOffset, MaxLines);
            return ana.ReportLargestAllocators();
        }


        /// <summary>
        /// total memory (aka. sum) over all MPI ranks over time
        /// </summary>
        static public Plot2Ddata GetMPItotalMemory(this IEnumerable<ISessionInfo> sessS, int LineOffset = 0, int MaxLines = 10000) {
            var ana = new SessionsComparisonMemtrace(
                sessS.Select(sess => new DirectoryInfo(sess.GetSessionDirectory())).ToArray(),
                LineOffset, MaxLines);

            int L = ana.NoOfTimeEntries;

            var ret = new Plot2Ddata();
            int ColorCounter = 0;
            foreach (var tt in ana.GetTotalMemoryMegs()) {

                var fmt = new PlotFormat();
                fmt.SetLineColorFromIndex(ColorCounter); ColorCounter++;
                fmt.WithStyle(Styles.Lines);

                ret.AddDataGroup(new XYvalues(
                    $"Tot Mem [MegB] at {tt.MPISz} cores",
                    L.ForLoop(i => (double)i),
                    tt.TotalMem),
                    fmt);
            }

            ret.Title = "Total memory of sessions at " + ana.MPIsizeOfRuns.ToConcatString("", ", ", "") + " MPI cores";


            return ret;
        }


        /// <summary>
        /// Reports the largest differences in memory allocation between the multiple runs
        /// </summary>
        static public (int TimelineIndex, double Imbalance, double[] AllocMegs, string Name)[] ReportLargestAllocatorImbalance(this IEnumerable<ISessionInfo> sessS, int LineOffset = 0, int MaxLines = 10000) {
            var ana = new SessionsComparisonMemtrace(
               sessS.Select(sess => new DirectoryInfo(sess.GetSessionDirectory())).ToArray(),
               LineOffset, MaxLines);
            return ana.ReportLargestAllocatorImbalance();
        }

        /// <summary>
        /// 
        /// </summary>
        public static void PlotData(Plot2Ddata pltDat, string xLabel, string yLabel, bool convData = false) {

            int lineColor = 0;
            PlotFormat format = new PlotFormat(lineColor: ((LineColors)(++lineColor)));
            Gnuplot gp = new Gnuplot(baseLineFormat: format);
            gp.SetXLabel(xLabel);
            gp.SetYLabel(yLabel);
            gp.Cmd("set grid xtics ytics");
            foreach (var group in pltDat.dataGroups) {
                if (convData) {
                    gp.PlotXY(group.Abscissas, group.Values, group.Name,
                        new PlotFormat(lineColor: ((LineColors)(++lineColor))), logX: true, logY: true);
                } else {
                    gp.PlotXY(group.Abscissas, group.Values, group.Name,
                        new PlotFormat(lineColor: ((LineColors)(++lineColor))));
                }
            }
            gp.WriteDeferredPlotCommands();
            gp.Execute();

        }


        /// <summary>
        /// Plots interface position/velocity and evaporative mass flux over time if a  "Evaporation.txt" exists.
        /// </summary>
        /// <param name="pSessions"> List of sessions to be evaluated </param>
        public static void EvalEvaporationData(this IEnumerable<ISessionInfo> pSessions) {

            int numberSessions = pSessions.Count();
            double[][] Stimesteps = new double[numberSessions][];
            double[][] SpositionI = new double[numberSessions][];
            double[][] SvelocityI = new double[numberSessions][];
            double[][] SevapMassI = new double[numberSessions][];

            // Read all data
            for (int j = 0; j < numberSessions; j++) {
                string path = pSessions.Pick(j).Database.Path + "\\sessions\\" + pSessions.Pick(j).ID + "\\Evaporation.txt";
                string[] lines = File.ReadAllLines(path);
                double[] timesteps = new double[lines.Length - 1];
                double[] positionI = new double[lines.Length - 1];
                double[] velocityI = new double[lines.Length - 1];
                double[] evapMassI = new double[lines.Length - 1];

                for (int i = 0; i < lines.Length - 1; i++) {
                    timesteps[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[1]);
                    positionI[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[2]);
                    velocityI[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[3]);
                    evapMassI[i] = Convert.ToDouble(lines[i + 1].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[4]);
                }
                Stimesteps[j] = timesteps;
                SpositionI[j] = positionI;
                SvelocityI[j] = velocityI;
                SevapMassI[j] = evapMassI;
            }

            // Build DataSets
            KeyValuePair<string, double[][]>[] dataRowsPositionI = new KeyValuePair<string, double[][]>[numberSessions];
            KeyValuePair<string, double[][]>[] dataRowsVelocityI = new KeyValuePair<string, double[][]>[numberSessions];
            KeyValuePair<string, double[][]>[] dataRowsEvapMassI = new KeyValuePair<string, double[][]>[numberSessions];
            for (int i = 0; i < numberSessions; i++) {
                dataRowsPositionI[i] = new KeyValuePair<string, double[][]>(pSessions.Pick(i).Name, new double[][] { Stimesteps[i], SpositionI[i] });
                dataRowsVelocityI[i] = new KeyValuePair<string, double[][]>(pSessions.Pick(i).Name, new double[][] { Stimesteps[i], SvelocityI[i] });
                dataRowsEvapMassI[i] = new KeyValuePair<string, double[][]>(pSessions.Pick(i).Name, new double[][] { Stimesteps[i], SevapMassI[i] });
            }
            Plot2Ddata Time_PositionI = new Plot2Ddata(dataRowsPositionI);
            Plot2Ddata Time_VelocityI = new Plot2Ddata(dataRowsVelocityI);
            Plot2Ddata Time_EvapMassI = new Plot2Ddata(dataRowsEvapMassI);

            // Plot interface position
            int lineColor = 0;
            PlotFormat format = new PlotFormat(lineColor: ((LineColors)(++lineColor)));
            Gnuplot gp = new Gnuplot(baseLineFormat: format);
            gp.SetXLabel("time");
            gp.SetYLabel("interface postion");
            gp.Cmd("set grid xtics ytics");
            foreach (var group in Time_PositionI.dataGroups) {
                gp.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                    new PlotFormat(lineColor: ((LineColors)(++lineColor))));
            }
            gp.WriteDeferredPlotCommands();
            gp.Execute();

            // Plot interface velocity
            lineColor = 0;
            Gnuplot gp2 = new Gnuplot(baseLineFormat: format);
            gp2.SetXLabel("time");
            gp2.SetYLabel("interface velocity");
            gp2.Cmd("set grid xtics ytics");
            foreach (var group in Time_VelocityI.dataGroups) {
                gp2.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                    new PlotFormat(lineColor: ((LineColors)(++lineColor))));
            }
            gp2.WriteDeferredPlotCommands();
            gp2.Execute();

            // Plot evaporative mass flux
            lineColor = 0;
            Gnuplot gp3 = new Gnuplot(baseLineFormat: format);
            gp3.SetXLabel("time");
            gp3.SetYLabel("evaporative mass flux");
            gp3.Cmd("set grid xtics ytics");
            foreach (var group in Time_EvapMassI.dataGroups) {
                gp3.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
                    new PlotFormat(lineColor: ((LineColors)(++lineColor))));
            }
            gp3.WriteDeferredPlotCommands();
            gp3.Execute();

        }



        // <summary>
        // Plots the temperature profile if a  "Evaporation.txt" exists.
        // </summary>
        //public static void PlotTemperatureProfileAt(this ISessionInfo pSession, int[] timestepIndex) {

        //    int numberTimesteps = timestepIndex.Count();
        //    double[][] tsTemperatureP = new double[numberTimesteps][];

        //    double L = (double)pSession.KeysAndQueries["AdditionalParameters[0]"];
        //    int len = 0;

        //    // Read all data
        //    string path = pSession.Database.Path + "\\sessions\\" + pSession.ID + "\\Evaporation.txt";
        //    string[] lines = File.ReadAllLines(path);
        //    for (int j = 0; j < numberTimesteps; j++) {

        //        string[] tsData = lines[timestepIndex[j]].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries);
        //        len = tsData.Count();

        //        double[] TemperatureP = new double[len - 5];
        //        for (int i = 5; i < len; i++) {
        //            TemperatureP[i - 5] = Convert.ToDouble(tsData[i]);
        //        }
        //        tsTemperatureP[j] = TemperatureP;
        //    }

        //    double[] profile = GenericBlas.Linspace(0, L, len - 5);

        //    // Build DataSets
        //    KeyValuePair<string, double[][]>[] dataRowsTemperatureP = new KeyValuePair<string, double[][]>[numberTimesteps];
        //    for (int i = 0; i < numberTimesteps; i++) {
        //        dataRowsTemperatureP[i] = new KeyValuePair<string, double[][]>(pSession.Name, new double[][] { profile, tsTemperatureP[i] });
        //    }
        //    Plot2Ddata Profile_Temperature = new Plot2Ddata(dataRowsTemperatureP);

        //    // Plot interface position
        //    int lineColor = 0;
        //    PlotFormat format = new PlotFormat(lineColor: ((LineColors)(++lineColor)));
        //    Gnuplot gp = new Gnuplot(baseLineFormat: format);
        //    gp.SetXLabel("profile");
        //    gp.SetYLabel("temperature");
        //    gp.Cmd("set grid xtics ytics");
        //    foreach (var group in Profile_Temperature.dataGroups) {
        //        gp.PlotXY(group.Abscissas, group.Values, group.Name.Split('.').Last(),
        //            new PlotFormat(lineColor: ((LineColors)(++lineColor))));
        //    }
        //    gp.WriteDeferredPlotCommands();
        //    gp.Execute();

        //}



    }
}
