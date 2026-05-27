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
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace BoSSS.PlotGenerator {

    /// <summary>
    /// Container for the main entry point
    /// </summary>
    public partial class PlotApplication {

        private static int rank;

        /// <summary>
        /// Main entry point. Uses <see cref="PlotApplication"/> to plot all fields
        /// referenced in the config file supplied as a command line argument
        /// </summary>
        /// <param name="args">
        /// Command line arguments. Currently, only the first entry is used
        /// which has to contain the path to the configuration file containing
        /// the plot information
        /// </param>
        static void Main(string[] args) {
            bool mustFinalizeMPI;

            var t = typeof(BoSSS.Application.Rheology.RheologyTimestepInfo);


            mustFinalizeMPI = BoSSS.Solution.Application.InitMPI();
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);

            WriteLine(" --------------------------------------------");
            WriteLine("| Process BoSSS.PlotGenerator.exe started   |");
            WriteLine("| Time: " + DateTime.Now.ToString() + "                  |");
            WriteLine(" --------------------------------------------");
            WriteLine("");

            //Debugger.Break();

            if (args.Length == 0 || args.Any(IsHelpRequest)) {
                WriteUsage();
                FinishMPI(mustFinalizeMPI);
                return;
            }

            WriteLine("args[0]: " + args[0]);

            FieldStateConfiguration config;
            string outputPath;

            if (IsXmlInput(args[0])) {
                outputPath = Path.GetDirectoryName(args[0]);

                // deserialize plotConfig.xml 
                WriteLine("Reading plotConfig.xml ...");
                config = FieldStateConfiguration.LoadConfig(args[0]);
            } else {
                WriteLine("Reading command line arguments ...");
                config = CreateConfigFromCommandLine(args, out outputPath);
            }

            if (!string.IsNullOrWhiteSpace(outputPath) && !Directory.Exists(outputPath)) {
                Directory.CreateDirectory(outputPath);
            }

            // In the next 5 lines of code you find all necessary data to create
            // plot files
            string[] basePaths = config.BasePaths;
            Guid sessionGuid = config.SessionGuid;
            FieldStateConfiguration.ExportFormats exportFormat = config.ExportFormat;
            List<string> fieldNameList = config.FieldNames ?? new List<string>();
            List<TimestepNumber> timeStepList = config.TimeSteps ?? new List<TimestepNumber>();

            // a few debug lines
            WriteLine("BasePaths[0]      : " + basePaths[0]);
            WriteLine("SessionGuid       : " + sessionGuid.ToString());
            WriteLine("ExportFormat      : " + exportFormat.ToString());
            WriteLine("Count(FieldNames) : " + (fieldNameList.Count == 0 ? "all" : fieldNameList.Count.ToString()));
            WriteLine("Count(TimeSteps)  : " + (timeStepList.Count == 0 ? "all" : timeStepList.Count.ToString()));
            WriteLine("output to         : " + outputPath);

            WriteLine("plotting...");

            PlotApplication app = new PlotApplication(config, outputPath);
            app.Run();

            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            FinishMPI(mustFinalizeMPI);
        }

        private static void WriteLine(string message) {
            if (rank == 0) {
                Console.WriteLine(message);
            }
        }

        private static bool IsXmlInput(string argument) {
            return argument.EndsWith(".xml", StringComparison.OrdinalIgnoreCase);
        }

        private static bool IsHelpRequest(string argument) {
            return argument == "-h" || argument == "--help" || argument == "/?";
        }

        private static void FinishMPI(bool mustFinalizeMPI) {
            if (mustFinalizeMPI) {
                BoSSS.Solution.Application.FinalizeMPI();
            }
        }

        private static FieldStateConfiguration CreateConfigFromCommandLine(string[] args, out string outputPath) {
            string sessionDirectory = Path.GetFullPath(args[0]);
            Guid sessionGuid = InferSessionGuid(sessionDirectory);
            string databasePath = InferDatabasePath(sessionDirectory);

            CommandLineOptions options = CommandLineOptions.Parse(args.Skip(1));
            if (options.HelpRequested) {
                WriteUsage();
            }

            IDatabaseInfo database = DatabaseInfo.Open(databasePath);
            ISessionInfo session = database.Controller.GetSessionInfo(sessionGuid);
            outputPath = options.OutputPath ?? Utils.GetExportDirectory(session);

            FieldStateConfiguration config = new FieldStateConfiguration();
            config.BasePaths = new string[] { databasePath };
            config.SessionGuid = sessionGuid;
            config.FieldNames = ParseFieldNames(options.Fields);
            config.TimeSteps = ParseTimesteps(options.Timesteps, session);
            config.SuperSampling = options.SuperSampling;
            config.PlotShadowFields = options.PlotShadowFields;
            config.ReconstructionType = options.ReconstructionType;
            config.GhostLevel = options.GhostLevel;
            config.NumberOfProcesses = 1;
            config.ExportFormat = options.ExportFormat;

            return config;
        }

        private static Guid InferSessionGuid(string sessionDirectory) {
            if (!Directory.Exists(sessionDirectory)) {
                throw new ArgumentException("Session directory does not exist: " + sessionDirectory);
            }

            string sessionGuidText = new DirectoryInfo(sessionDirectory).Name;
            if (!Guid.TryParse(sessionGuidText, out Guid sessionGuid)) {
                throw new ArgumentException("Expected the session directory name to be a GUID: " + sessionDirectory);
            }

            return sessionGuid;
        }

        private static string InferDatabasePath(string sessionDirectory) {
            DirectoryInfo sessionDir = new DirectoryInfo(sessionDirectory);
            DirectoryInfo sessionsDir = sessionDir.Parent;
            if (sessionsDir == null || !sessionsDir.Name.Equals("sessions", StringComparison.OrdinalIgnoreCase)) {
                throw new ArgumentException("Expected a session directory below a 'sessions' directory: " + sessionDirectory);
            }

            DirectoryInfo databaseDir = sessionsDir.Parent;
            if (databaseDir == null) {
                throw new ArgumentException("Could not infer database directory from: " + sessionDirectory);
            }

            return databaseDir.FullName;
        }

        private static List<string> ParseFieldNames(string fields) {
            if (string.IsNullOrWhiteSpace(fields) || fields.Equals("all", StringComparison.OrdinalIgnoreCase)) {
                return new List<string>();
            }

            return SplitList(fields).Distinct().ToList();
        }

        private static List<TimestepNumber> ParseTimesteps(string spec, ISessionInfo session) {
            if (string.IsNullOrWhiteSpace(spec) || spec.Equals("all", StringComparison.OrdinalIgnoreCase)) {
                return new List<TimestepNumber>();
            }

            TimestepNumber[] available = session.Timesteps.Select(tsi => tsi.TimeStepNumber).OrderBy(ts => ts).ToArray();
            if (available.Length == 0) {
                throw new ArgumentException("Session contains no timesteps: " + session.ID);
            }

            SortedSet<TimestepNumber> selected = new SortedSet<TimestepNumber>();
            foreach (string token in SplitList(spec)) {
                AddTimesteps(token, available, selected);
            }

            if (selected.Count == 0) {
                throw new ArgumentException("Timestep selection did not match any timesteps: " + spec);
            }

            return selected.ToList();
        }

        private static void AddTimesteps(string token, TimestepNumber[] available, SortedSet<TimestepNumber> selected) {
            if (token.Equals("first", StringComparison.OrdinalIgnoreCase)) {
                selected.Add(available.First());
                return;
            }

            if (token.Equals("last", StringComparison.OrdinalIgnoreCase)) {
                selected.Add(available.Last());
                return;
            }

            if (token.Contains("..")) {
                string[] range = token.Split(new string[] { ".." }, StringSplitOptions.None);
                AddRange(range, 1, token, available, selected);
                return;
            }

            string[] colonRange = token.Split(':');
            if (colonRange.Length == 2 || colonRange.Length == 3) {
                int stride = 1;
                if (colonRange.Length == 3 && !int.TryParse(colonRange[2], NumberStyles.Integer, CultureInfo.InvariantCulture, out stride)) {
                    throw new ArgumentException("Could not parse timestep range stride in '" + token + "'");
                }

                AddRange(colonRange.Take(2).ToArray(), stride, token, available, selected);
                return;
            }

            selected.Add(ParseTimestepNumber(token));
        }

        private static void AddRange(string[] range, int stride, string token, TimestepNumber[] available, SortedSet<TimestepNumber> selected) {
            if (range.Length != 2 || string.IsNullOrWhiteSpace(range[0]) || string.IsNullOrWhiteSpace(range[1])) {
                throw new ArgumentException("Could not parse timestep range '" + token + "'");
            }

            if (stride <= 0) {
                throw new ArgumentException("Timestep range stride must be positive in '" + token + "'");
            }

            TimestepNumber first = ParseTimestepNumber(range[0]);
            TimestepNumber last = ParseTimestepNumber(range[1]);
            int ordinal = 0;

            foreach (TimestepNumber timestep in available.Where(ts => ts.CompareTo(first) >= 0 && ts.CompareTo(last) <= 0)) {
                if (ordinal % stride == 0) {
                    selected.Add(timestep);
                }
                ordinal++;
            }
        }

        private static TimestepNumber ParseTimestepNumber(string text) {
            if (TimestepNumber.TryParse(text, out TimestepNumber timestep)) {
                return timestep;
            }

            throw new ArgumentException("Could not parse timestep number '" + text + "'");
        }

        private static IEnumerable<string> SplitList(string text) {
            return text.Split(new char[] { ',', ';' }, StringSplitOptions.RemoveEmptyEntries)
                .Select(s => s.Trim())
                .Where(s => s.Length > 0);
        }

        private static void WriteUsage() {
            WriteLine("Usage:");
            WriteLine("  BoSSS.PlotGenerator <plotConfig.xml>");
            WriteLine("  BoSSS.PlotGenerator <database>/sessions/<session-guid> [options]");
            WriteLine("");
            WriteLine("Options:");
            WriteLine("  --supersampling <n>       Plot supersampling level, default 0");
            WriteLine("  --timesteps <spec>        all, first, last, comma list, ranges 0..10 or 0:10[:stride]");
            WriteLine("  --fields <names>          comma-separated field names, default all");
            WriteLine("  --format <format>         TecPlot, CSV, Curve, BinaryStream, Cgns, Hdf5");
            WriteLine("  --output <path>           output directory, default BoSSSpad session export directory");
            WriteLine("  --reconstruction <mode>   None or Average");
            WriteLine("  --ghost <0|1>             ghost level, default 0");
            WriteLine("  --shadow-fields           plot XDG fields specieswise");
        }

        private class CommandLineOptions {
            public string Timesteps;
            public string Fields;
            public string OutputPath;
            public int SuperSampling;
            public bool PlotShadowFields;
            public bool HelpRequested;
            public FieldStateConfiguration.ReconstructionTypes ReconstructionType = FieldStateConfiguration.ReconstructionTypes.None;
            public FieldStateConfiguration.GhostLevels GhostLevel = FieldStateConfiguration.GhostLevels.Zero;
            public FieldStateConfiguration.ExportFormats ExportFormat = FieldStateConfiguration.ExportFormats.TecPlot;

            public static CommandLineOptions Parse(IEnumerable<string> args) {
                CommandLineOptions options = new CommandLineOptions();
                string pendingKey = null;

                foreach (string rawArg in args) {
                    string arg = rawArg.Trim();
                    if (arg.Length == 0) {
                        continue;
                    }

                    if (pendingKey != null) {
                        options.ApplyOption(pendingKey, arg);
                        pendingKey = null;
                        continue;
                    }

                    if (IsHelpRequest(arg)) {
                        options.HelpRequested = true;
                        continue;
                    }

                    if (arg.StartsWith("--", StringComparison.Ordinal)) {
                        string option = arg.Substring(2);
                        int separator = option.IndexOf('=');
                        if (separator >= 0) {
                            options.ApplyOption(option.Substring(0, separator), option.Substring(separator + 1));
                        } else if (options.IsBooleanOption(option)) {
                            options.ApplyOption(option, "true");
                        } else {
                            pendingKey = option;
                        }
                    } else if (arg.StartsWith("-", StringComparison.Ordinal)) {
                        string option = ExpandShortOption(arg.Substring(1));
                        if (options.IsBooleanOption(option)) {
                            options.ApplyOption(option, "true");
                        } else {
                            pendingKey = option;
                        }
                    } else if (options.Timesteps == null) {
                        options.Timesteps = arg;
                    } else {
                        options.Timesteps += "," + arg;
                    }
                }

                if (pendingKey != null) {
                    throw new ArgumentException("Missing value for command line option '--" + pendingKey + "'");
                }

                return options;
            }

            private static string ExpandShortOption(string option) {
                switch (option) {
                    case "s":
                        return "supersampling";
                    case "t":
                        return "timesteps";
                    case "f":
                        return "fields";
                    case "o":
                        return "output";
                    default:
                        return option;
                }
            }

            private bool IsBooleanOption(string key) {
                return NormalizeKey(key) == "shadowfields";
            }

            private void ApplyOption(string key, string value) {
                switch (NormalizeKey(key)) {
                    case "supersampling":
                        SuperSampling = int.Parse(value, CultureInfo.InvariantCulture);
                        if (SuperSampling < 0) {
                            throw new ArgumentException("Supersampling must be non-negative.");
                        }
                        break;
                    case "timesteps":
                    case "timestep":
                        Timesteps = value;
                        break;
                    case "fields":
                    case "field":
                        Fields = value;
                        break;
                    case "format":
                        ExportFormat = ParseEnum<FieldStateConfiguration.ExportFormats>(value, key);
                        break;
                    case "output":
                    case "out":
                        OutputPath = Path.GetFullPath(value);
                        break;
                    case "reconstruction":
                        ReconstructionType = ParseEnum<FieldStateConfiguration.ReconstructionTypes>(value, key);
                        break;
                    case "ghost":
                    case "ghostlevel":
                        GhostLevel = ParseGhostLevel(value);
                        break;
                    case "shadowfields":
                        PlotShadowFields = bool.Parse(value);
                        break;
                    default:
                        throw new ArgumentException("Unknown command line option '--" + key + "'");
                }
            }

            private static T ParseEnum<T>(string value, string key) where T : struct {
                if (Enum.TryParse(value, true, out T result)) {
                    return result;
                }

                throw new ArgumentException("Could not parse value '" + value + "' for option '--" + key + "'");
            }

            private static FieldStateConfiguration.GhostLevels ParseGhostLevel(string value) {
                if (value == "0") {
                    return FieldStateConfiguration.GhostLevels.Zero;
                }

                if (value == "1") {
                    return FieldStateConfiguration.GhostLevels.One;
                }

                return ParseEnum<FieldStateConfiguration.GhostLevels>(value, "ghost");
            }

            private static string NormalizeKey(string key) {
                return key.Replace("-", "").Replace("_", "").ToLowerInvariant();
            }
        }
    }
}
