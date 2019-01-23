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
using System.Diagnostics;
using System.IO;
using BoSSS.Foundation.IO;
using Mono.CSharp;
using System.Collections.Generic;
using System.Reflection;
using BoSSS.Solution.Gnuplot;
using System.Linq;
using ilPSP;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Extends/replaces the standard commands provided by
    /// <see cref="InteractiveBase"/> by some BoSSSPad-specific stuff.
    /// </summary>
    public class InteractiveShell : InteractiveBase {

        /// <summary>
        /// Provides a help text
        /// </summary>
        new public static string help {
            get {
                try {
                    string dbeCommandOverviewDocPath = Path.Combine(
                        Utils.GetBoSSSInstallDir(),
                        "doc",
                        "BoSSSPad_Command_Overview.pdf");
                    System.Diagnostics.Process.Start(dbeCommandOverviewDocPath);
                    return "Displaying external help document...";
                } catch (Exception e) {
                    return "Displaying external help failed ( " + e.GetType().Name + ": " + e.Message + ")";
                }
            }
        }

        /// <summary>
        /// <see cref="help"/>
        /// </summary>
        public static string Help {
            get {
                return help;
            }
        }

        /// <summary>
        /// Holds the last exception that has been thrown during the execution.
        /// </summary>
        public static Exception LastError {
            get;
            set;
        }

        /// <summary>
        /// Holds the result of the last operation
        /// </summary>
        public static object LastResult {
            get {
                return ans;
            }
        }

        /// <summary>
        /// Holds the result of the last operation
        /// </summary>
        public static object ans {
            get;
            internal set;
        }

        /// <summary>
        /// Opens the folder containing config files like the DBE.xml
        /// </summary>
        public static void OpenConfigDirectory() {
            string dbeXmlPath = Path.Combine(Utils.GetBoSSSUserSettingsPath(), "etc");
            Process.Start(dbeXmlPath);
        }

        /// <summary>
        /// Opens the DBE.xml
        /// </summary>
        public static void OpenConfigFile() {
            string dbeXmlPath = Path.Combine(
                Utils.GetBoSSSUserSettingsPath(), "etc", "DBE.xml");
            Process.Start(dbeXmlPath);
        }

        /// <summary>
        /// Saves the current interactive session as a worksheet that can be
        /// loaded by the worksheet edition of the BoSSSPad
        /// </summary>
        /// <param name="path"></param>
        public static void SaveSessionAsWorksheet(string path) {
            ReadEvalPrintLoop.SaveSessionAsWorksheet(path);
        }

        /// <summary>
        /// Clears the console window.
        /// </summary>
        public static void Clear() {
            Console.Clear();
        }

        private static WorkflowMgm m_WorkflowMgm;

        /// <summary>
        /// Link to the workflow-management facility
        /// </summary>
        public static WorkflowMgm WorkflowMgm {
            get {
                if (m_WorkflowMgm == null)
                    m_WorkflowMgm = new WorkflowMgm();
                return m_WorkflowMgm;
            }
        }

        /// <summary>
        /// Reset states which survive a restart of the interpreter.
        /// </summary>
        static internal void Reset() {
            databases = new IDatabaseInfo[0];
            m_WorkflowMgm = null;
        }



        /// <summary>
        /// All the databases; the workflow-management (see <see cref="WorkflowMgm"/>) must have access to those.
        /// </summary>
        public static IList<IDatabaseInfo> databases;

        /// <summary>
        /// path to the default BoSSS database directory for the current user
        /// </summary>
        static public string GetDefaultDatabaseDir() {
            string basepath = null;
            basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");

            if (basepath.IsEmptyOrWhite())
                basepath = System.Environment.GetEnvironmentVariable("HOME");

            string path = Path.Combine(basepath, "default_bosss_db");

            return path;
        }


        /// <summary>
        /// Similar to <see cref="OpenOrCreateDatabase"/>, but uses a default database
        /// </summary>
        static public IDatabaseInfo OpenOrCreateDefaultDatabase() {
            string path = GetDefaultDatabaseDir();
            return OpenOrCreateDatabase(path);
        }
        
        /// <summary>
        /// Opens a database at a specific path, resp. creates one if the 
        /// </summary>
        static public IDatabaseInfo OpenOrCreateDatabase(string dbDir) {
            if (Directory.Exists(dbDir)) {
                if (!DatabaseUtils.IsValidBoSSSDatabase(dbDir)) {
                    throw new ArgumentException("Directory '" + dbDir + "' exists, but is not a valid BoSSS database.");
                }
                Console.WriteLine("Opening existing database '" + dbDir + "'.");
            } else {
                DatabaseUtils.CreateDatabase(dbDir);
                Console.WriteLine("Creating database '" + dbDir + "'.");
            }

            var dbi = new DatabaseInfo(dbDir);

            List<IDatabaseInfo> mod_databases = new List<IDatabaseInfo>();
            if (databases != null) {
                mod_databases.AddRange(databases);
            }
            mod_databases.Add(dbi);
            databases = mod_databases.ToArray();

            return dbi;
        }

        static internal Document CurrentDoc = null;

        /// <summary>
        /// Extracts the source code of some function, which can be used as an initial value or boundary condition.
        /// </summary>
        /// <param name="f">
        /// Must be the reference to a static method of a static class.
        /// </param>
        static public BoSSS.Solution.Control.Formula GetFormulaObject(Func<double[], double> f) {
            return GetFormulaObject(f, false);
        }

        /// <summary>
        /// Extracts the source code of some function, which can be used as an initial value or boundary condition.
        /// </summary>
        /// <param name="f">
        /// Must be the reference to a static method of a static class.
        /// </param>
        static public BoSSS.Solution.Control.Formula GetFormulaObject(Func<double[], double, double> f) {
            return GetFormulaObject(f, true);
        }

        private static Solution.Control.Formula GetFormulaObject(System.Delegate f, bool timedep) {
            if (CurrentDoc == null) {
                throw new NotSupportedException("Only supported when a bws-document is present (GUI or batch mode).");
            }
            if (f == null)
                throw new ArgumentNullException();
            Assembly SearchedAssembly = f.Method.DeclaringType.Assembly;

            if (SearchedAssembly == null)
                throw new ApplicationException("Unable to find some assembly for delegate.");

            string AssemblyCode = null;
            foreach (var entry in CurrentDoc.CommandAndResult) {
                if (SearchedAssembly.Equals(entry.AssemblyProduced)) {
                    AssemblyCode = entry.Command;
                }
            }
            if (AssemblyCode == null) {
                throw new ApplicationException("Unable to find code of " + SearchedAssembly.FullName);
            } else {
                //Console.WriteLine("Found code:");
                //Console.WriteLine(AssemblyCode);
            }

            var fo = new Solution.Control.Formula(f.Method.DeclaringType.Name + "." + f.Method.Name, timedep, AssemblyCode);

            return fo;
        }

        static internal string _CurrentDocFile = null;

        /// <summary>
        /// Path of the current file.
        /// </summary>
        static public string CurrentDocFile {
            get {
                if (_CurrentDocFile == null) {
                    Console.WriteLine("No current document/not saved yet.");
                    return null;
                }

                if (!File.Exists(_CurrentDocFile)) {
                    Console.WriteLine("Document path '{0}' seems non-existent.", _CurrentDocFile);
                }

                return _CurrentDocFile;
            }
        }

        /// <summary>
        /// Directory where the current file is stored.
        /// </summary>
        static public string CurrentDocDir {
            get {
                string f = CurrentDocFile;
                if (f == null)
                    return null;
                return Path.GetDirectoryName(f);
            }
        }

        /// <summary>
        /// Simple plotting interface
        /// </summary>
        /// <returns>Output of <see cref="GnuplotExtensions.PlotNow(Gnuplot)"/></returns>
        static public object Plot(IEnumerable<double> X1, IEnumerable<double> Y1, string Name1 = null, string Format1 = null,
            IEnumerable<double> X2 = null, IEnumerable<double> Y2 = null, string Name2 = null, string Format2 = null,
            IEnumerable<double> X3 = null, IEnumerable<double> Y3 = null, string Name3 = null, string Format3 = null,
            IEnumerable<double> X4 = null, IEnumerable<double> Y4 = null, string Name4 = null, string Format4 = null,
            IEnumerable<double> X5 = null, IEnumerable<double> Y5 = null, string Name5 = null, string Format5 = null,
            IEnumerable<double> X6 = null, IEnumerable<double> Y6 = null, string Name6 = null, string Format6 = null,
            IEnumerable<double> X7 = null, IEnumerable<double> Y7 = null, string Name7 = null, string Format7 = null,
            bool logX = false, bool logY = false) {

            using (var gp = new Gnuplot()) {


                IEnumerable<double>[] Xs = new[] { X1, X2, X3, X4, X5, X6, X7 };
                IEnumerable<double>[] Ys = new[] { Y1, Y2, Y3, Y4, Y5, Y6, Y7 };
                string[] Ns = new string[] { Name1, Name2, Name3, Name4, Name5, Name6, Name7 };
                string[] Fs = new string[] { Format1, Format2, Format3, Format4, Format5, Format6, Format7 };

                for (int i = 0; i < 7; i++) {
                    if (Ys[i] != null) {
                        var f1 = new PlotFormat();
                        if (Fs[i] != null)
                            f1.FromString(Fs[i]);
                        gp.PlotXY(Xs[i], Ys[i], title: Ns[i], format: f1, logX: logX, logY: logY);
                    }
                }

                return gp.PlotNow();
            }
        }


        /// <summary>
        /// Driver interface for the <see cref="BoSSS.Solution.Tecplot.Tecplot"/> functionality.
        /// </summary>
        static public void Tecplot(string filename, params BoSSS.Foundation.DGField[] flds) {
            if (flds == null || flds.Length <= 0) {
                Console.WriteLine("No DG fields specified - not writing any output files.");
                return;
            }

            int susamp = Math.Max(flds.Max(f => f.Basis.Degree) + 1, 3);


            Tecplot(filename, 0.0, susamp, flds);
        }


        /// <summary>
        /// Driver interface for the <see cref="BoSSS.Solution.Tecplot.Tecplot"/> functionality.
        /// </summary>
        static public void Tecplot(string filename, double time, int supersampling, params BoSSS.Foundation.DGField[] flds) {
            if (flds == null || flds.Length <= 0) {
                Console.WriteLine("No DG fields specified - not writing any output files.");
                return;
            }

            if (supersampling > 3) {
                Console.WriteLine("Plotting with a supersampling greater than 3 is deactivated because it would very likely exceed this machines memory.");
                Console.WriteLine("Higher supersampling values are supported by external plot application.");
                supersampling = 3;
            }

            string directory = Path.GetDirectoryName(filename);
            string FullPath;
            if (directory == null || directory.Length <= 0) {
                directory = CurrentDocDir ?? "";
                FullPath = Path.Combine(directory, filename);
            } else {
                FullPath = filename;
            }

            Console.WriteLine("Writing output file {0}...", FullPath);


            BoSSS.Solution.Tecplot.Tecplot.PlotFields(flds, FullPath, time, supersampling);

        }

    }
}

