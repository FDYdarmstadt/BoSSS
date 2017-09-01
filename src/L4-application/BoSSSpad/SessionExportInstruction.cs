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
using ilPSP;
using BoSSS.PlotGenerator;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Configuration options for the plotting of a <see cref="ISessionInfo"/>.
    /// </summary>
    public class SessionExportInstruction : ExportInstruction {

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is
        /// <see cref="FieldStateConfiguration.ReconstructionTypes.None"/>.
        /// </summary>
        private FieldStateConfiguration.ReconstructionTypes reconstructionType =
            FieldStateConfiguration.ReconstructionTypes.None;

        /// <summary>
        /// Retriever for the names of the fields to be plotted
        /// </summary>
        private Lazy<List<string>> fieldNames;

        /// <summary>
        /// Retriever for the indices of the time steps to be plotted.
        /// </summary>
        private Lazy<List<TimestepNumber>> timeSteps;

        /// <summary>
        /// Creates a plot instruction for the plotting of
        /// <paramref name="session"/>.
        /// </summary>
        /// <param name="session">
        /// The session to be plotted.
        /// </param>
        public SessionExportInstruction(ISessionInfo session) {
            this.Session = session;
            this.fieldNames = new Lazy<List<string>>(
                () => session.Timesteps
                    .SelectMany(ts => ts.FieldInitializers
                        .Select(fi => fi.Identification))
                        .Distinct().ToList());
            this.timeSteps = new Lazy<List<TimestepNumber>>(
                () => session.Timesteps.Select(t => t.TimeStepNumber).ToList());
        }

        /// <summary>
        /// The session to be plotted.
        /// </summary>
        public ISessionInfo Session {
            get;
            set;
        }

        /// <summary>
        /// The names of the fields to be plotted. Default is all fields.
        /// </summary>
        public List<string> FieldNames {
            get {
                return fieldNames.Value;
            }
            set {
                fieldNames = new Lazy<List<string>>(() => value);
            }
        }

        /// <summary>
        /// The indices of the time-steps to be plotted. Default is all
        /// time-steps.
        /// </summary>
        public List<TimestepNumber> TimeSteps {
            get {
                return timeSteps.Value;
            }
            set {
                timeSteps = new Lazy<List<TimestepNumber>>(() => value);
            }
        }

        /// <summary>
        /// The type of reconstruction to be used. Default is
        /// <see cref="FieldStateConfiguration.ReconstructionTypes.None"/>.
        /// </summary>
        public FieldStateConfiguration.ReconstructionTypes ReconstructionType {
            get {
                return reconstructionType;
            }
            set {
                reconstructionType = value;
            }
        }

        private FieldStateConfiguration CreateConfiguration() {
            FieldStateConfiguration fsConfig = new FieldStateConfiguration();
            fsConfig.BasePaths = new string[] { Session.Database.Path };
            fsConfig.SessionGuid = Session.ID;
            fsConfig.FieldNames = FieldNames;
            fsConfig.TimeSteps = TimeSteps;
            fsConfig.SuperSampling = SuperSampling;
            fsConfig.ReconstructionType = reconstructionType;
            fsConfig.GhostLevel = GhostLevels;
            fsConfig.NumberOfProcesses = NumberOfProcesses;
            fsConfig.ExportFormat = ExportFormat;
            return fsConfig;
        }

        private string PlotDirPath {
            get {
                if (AlternativeDirectoryName == null) {
                    return Path.Combine(
                        Utils.GetExportOutputPath(),
                        "sessions",
                        Session.ID.ToString());
                } else {
                    return AlternativeDirectoryName;
                }
            }
        }

        /// <summary>
        /// Starts the export in an external process window.
        /// </summary>
        public override string YouMust() {
            Console.Write("Starting export process... ");
            PrepareProcess(CreateConfiguration(), PlotDirPath).Start();
            Console.WriteLine("Data will be written to the following directory:");
            return PlotDirPath;
        }

        /// <summary>
        /// Starts the export while suppressing the output and opens the first
        /// exported file in an external viewer.
        /// </summary>
        public void AndOpenYouMust() {
            Process process = PrepareProcess(CreateConfiguration(), PlotDirPath);
            process.StartInfo.UseShellExecute = false;
            process.StartInfo.RedirectStandardInput = true;
            process.StartInfo.RedirectStandardOutput = true;

            Console.Write("Starting export (this may take a while)...");
            process.Start();

            // Try to find out when process stops writing data (without
            // resorting to threads...)
            bool probablyFinished = false;
            while (true) {
                Thread.Sleep(300);
                if (process.StandardOutput.Peek() == -1) {
                    if (probablyFinished) {
                        break;
                    } else {
                        probablyFinished = true;
                    }
                } else {
                    process.StandardOutput.ReadLine();
                    probablyFinished = false;
                }
            }

            // Automated key-press
            process.StandardInput.WriteLine(true);
            if (process.WaitForExit(2000)) {
                Console.WriteLine("Done.");
            } else {
                process.Kill();
                Console.WriteLine("The process did not end as expected and has been aborted");
            }

            if (process.ExitCode != 0) {
                return;
            }

            DirectoryInfo dir = new DirectoryInfo(PlotDirPath);
            FileInfo[] matchingFiles = dir.GetFiles("state_" + TimeSteps.First() + "*");

            if (matchingFiles.IsNullOrEmpty()) {
                throw new Exception(
                    "Could not find plot file. This should not have happened");
            }

            FileInfo file = matchingFiles.OrderBy(f => f.CreationTime).Last();

            Console.WriteLine("Opening file in external viewer.");
            process = new Process();
            process.StartInfo.FileName = Path.Combine(file.DirectoryName, file.Name);
            process.Start();
        }

        /// <summary>
        /// Synonym for <see cref="AndOpenYouMust"/>.
        /// </summary>
        public void DoAndOpen() {
            AndOpenYouMust();
        }
    }

    /// <summary>
    /// Specialized components of the fluent interface defined in
    /// <see cref="ExportInstructionExtensions"/> for <see cref="SessionExportInstruction"/>
    /// </summary>
    public static class SessionExportInstructionExtensions {

        /// <summary>
        /// Modifies <see cref="SessionExportInstruction.FieldNames"/>
        /// </summary>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="fieldNames">
        /// The names of the fields to be plotted.
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static SessionExportInstruction WithFields(
            this SessionExportInstruction instruction, params string[] fieldNames) {
            instruction.FieldNames = fieldNames.ToList();
            return instruction;
        }

        /// <summary>
        /// Modifies <see cref="SessionExportInstruction.TimeSteps"/>
        /// </summary>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="timeSteps">
        /// The time-steps to be plotted
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static SessionExportInstruction WithTimesteps(
            this SessionExportInstruction instruction, params TimestepNumber[] timeSteps) {
            instruction.TimeSteps = timeSteps.ToList();
            return instruction;
        }

        /// <summary>
        /// Modifies <see cref="SessionExportInstruction.TimeSteps"/>
        /// </summary>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="FirstTimestep">
        /// The first time-step (<see cref="TimestepNumber.MajorNumber"/>) to be plotted
        /// </param>
        /// <param name="LastTimestep">
        /// The last time-step (<see cref="TimestepNumber.MajorNumber"/>) to be plotted
        /// </param>
        /// <returns></returns>
        public static SessionExportInstruction WithTimesteps(
            this SessionExportInstruction instruction, int FirstTimestep, int LastTimestep) {
                List<TimestepNumber> timeSteps = new List<TimestepNumber>();
                for (int i = FirstTimestep; i <= LastTimestep; i++)
                    timeSteps.Add(i);
                instruction.TimeSteps = timeSteps;
                return instruction;
        }

        /// <summary>
        /// Modifies <see cref="SessionExportInstruction.ReconstructionType"/>
        /// </summary>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static SessionExportInstruction WithReconstruction(this SessionExportInstruction instruction) {
            instruction.ReconstructionType = FieldStateConfiguration.ReconstructionTypes.Average;
            return instruction;
        }
    }
}
