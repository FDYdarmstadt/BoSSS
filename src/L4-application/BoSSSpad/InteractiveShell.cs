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
        /// Opens a database at a specific path, resp. creates one if the 
        /// </summary>
        /// <param name="dbDir"></param>
        /// <returns></returns>
        static public IDatabaseInfo OpenOrCreateDatabase(string dbDir) {
            if(Directory.Exists(dbDir)) {
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
            if(databases != null) {
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
                Console.WriteLine("Found code:");
                Console.WriteLine(AssemblyCode);
            }

            var fo = new Solution.Control.Formula(f.Method.DeclaringType.Name + "." + f.Method.Name, timedep, AssemblyCode);

            return fo;
        }
    }
}
