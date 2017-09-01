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

using BoSSS.Platform;
using BoSSS.PlotGenerator;
using System;
using System.Diagnostics;
using System.IO;
using System.Reflection;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Base class for all plotting instructions defining some common
    /// configuration options an their default values.
    /// </summary>
    public abstract class ExportInstruction {

        /// <summary>
        /// Name of the MPI executable.
        /// </summary>
        public const string MPI_EXECUTABLE = "mpiexec";

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is zero.
        /// </summary>
        public int SuperSampling {
            get;
            set;
        }

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is none.
        /// </summary>
        public FieldStateConfiguration.GhostLevels GhostLevels {
            get;
            set;
        }

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is one.
        /// </summary>
        private uint numberOfProcesses = 1;

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is one.
        /// </summary>
        public uint NumberOfProcesses {
            get {
                return numberOfProcesses;
            }
            set {
                numberOfProcesses = value;
            }
        }

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is Tecplot.
        /// </summary>
        private FieldStateConfiguration.ExportFormats exportFormat =
            FieldStateConfiguration.ExportFormats.TecPlot;

        /// <summary>
        /// <see cref="FieldStateConfiguration"/>. Default is Tecplot.
        /// </summary>
        public FieldStateConfiguration.ExportFormats ExportFormat {
            get {
                return exportFormat;
            }
            set {
                exportFormat = value;
            }
        }

        /// <summary>
        /// By default, the exported files are saved in a default directory
        /// specified by classes deriving from this class. If this is field
        /// not null, its value is taken as an alternative directory name.
        /// </summary>
        public string AlternativeDirectoryName {
            get;
            set;
        }

        /// <summary>
        /// Does the actual plotting using the plot generator
        /// </summary>
        /// <param name="config">
        /// The field state configuration to be plotted
        /// </param>
        /// <param name="outputPath">
        /// Path where the output should be placed
        /// </param>
        protected Process PrepareProcess(FieldStateConfiguration config, string outputPath) {
            if (Directory.Exists(outputPath) == false) {
                Directory.CreateDirectory(outputPath);
            }
            string plotConfigPath = Path.Combine(outputPath, "plotConfig.xml");
            FieldStateConfiguration.Serialize(plotConfigPath, config);

            // we are expecting 'BoSSS.PlotGen.exe' to be in the same dir as DBE.exe
            Assembly a = System.Reflection.Assembly.GetEntryAssembly();
            string plotGenPath = Path.Combine(
                System.IO.Path.GetDirectoryName(a.Location),
                "BoSSS.PlotGenerator.exe");

            if (!File.Exists(plotGenPath)) {
                throw new Exception(plotGenPath + " could not be found.");
            }

            Process plotProcess = new Process();
            plotProcess.StartInfo.WindowStyle = ProcessWindowStyle.Minimized;
            if (config.NumberOfProcesses == 1) {
                plotProcess.StartInfo.FileName = plotGenPath;
            } else {
                plotProcess.StartInfo.FileName = MPI_EXECUTABLE;
                plotProcess.StartInfo.Arguments = string.Format(
                    " -n {0} {1} ", config.NumberOfProcesses, plotGenPath);
            }
            plotProcess.StartInfo.Arguments += "\"" + plotConfigPath + "\"";

            return plotProcess;
        }

        /// <summary>
        /// Starts the export using the configuration represented by this
        /// object.
        /// </summary>
        /// <returns>
        /// Export path
        /// </returns>
        /// <remarks>
        /// If you are curious about the naming of this method, see
        /// <see cref="ExportInstructionExtensions"/>.
        /// </remarks>
        public abstract string YouMust();

        /// <summary>
        /// See <see cref="YouMust()"/>
        /// </summary>
        /// <returns>
        /// Export path
        /// </returns>
        public string Do() {
            return YouMust();
        }
    }

    /// <summary>
    /// Creates a fluent interface to create a plot instruction. That is, it
    /// allows you to chain the setting of configuration options.
    /// <example>
    /// Using
    /// <code>
    /// Plot.Sessions($session).WithSuperSampling(2).WitMPI(2).YouMust()
    /// </code>
    /// will export the $session with two times super-sampling and using two
    /// processors.
    /// </example>
    /// </summary>
    public static class ExportInstructionExtensions {

        /// <summary>
        /// Activates super-sampling in <paramref name="instruction"/>.
        /// </summary>
        /// <typeparam name="T">
        /// The sub-type of <see cref="ExportInstruction"/> to be used.
        /// </typeparam>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="superSampling">
        /// <see cref="FieldStateConfiguration"/>
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static T WithSupersampling<T>(this T instruction, int superSampling)
            where T : ExportInstruction {
            instruction.SuperSampling = superSampling;
            return instruction;
        }

        /// <summary>
        /// Activates ghost cells in <paramref name="instruction"/>
        /// </summary>
        /// <typeparam name="T">
        /// The sub-type of <see cref="ExportInstruction"/> to be used.
        /// </typeparam>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static T WithGhostCells<T>(this T instruction)
            where T : ExportInstruction {
            instruction.GhostLevels = FieldStateConfiguration.GhostLevels.One;
            return instruction;
        }

        /// <summary>
        /// Activates ghost cells in <paramref name="instruction"/>
        /// </summary>
        /// <typeparam name="T">
        /// The sub-type of <see cref="ExportInstruction"/> to be used.
        /// </typeparam>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="numberOfProcesses">
        /// <see cref="FieldStateConfiguration"/>
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static T WithMPI<T>(this T instruction, uint numberOfProcesses)
            where T : ExportInstruction {
            instruction.NumberOfProcesses = numberOfProcesses;
            return instruction;
        }

        /// <summary>
        /// Changes the format in <paramref name="instruction"/>
        /// </summary>
        /// <typeparam name="T">
        /// The sub-type of <see cref="ExportInstruction"/> to be used.
        /// </typeparam>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="format">
        /// <see cref="FieldStateConfiguration"/>
        /// </param>
        /// <returns>
        /// <paramref name="instruction"/>
        /// </returns>
        public static T WithFormat<T>(this T instruction, string format)
            where T : ExportInstruction {
            instruction.ExportFormat = Enum<FieldStateConfiguration.ExportFormats>.Parse(format, true);
            return instruction;
        }

        /// <summary>
        /// Changes <see cref="ExportInstruction.AlternativeDirectoryName"/> in
        /// the given <paramref name="instruction"/>.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="instruction">
        /// The instruction to be modified
        /// </param>
        /// <param name="targetDirectory">
        /// The new target directory of the export
        /// </param>
        /// <returns>
        /// The modified <paramref name="instruction"/>
        /// </returns>
        public static T To<T>(this T instruction, string targetDirectory)
            where T : ExportInstruction {
            instruction.AlternativeDirectoryName = targetDirectory;
            return instruction;
        }
    }
}
