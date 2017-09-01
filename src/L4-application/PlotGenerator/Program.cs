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
using System.IO;

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
            
            ilPSP.Environment.Bootstrap(args, BoSSS.Solution.Application.GetBoSSSInstallDir(), out mustFinalizeMPI);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);

            WriteLine(" --------------------------------------------");
            WriteLine("| Process BoSSS.PlotGenerator.exe started   |");
            WriteLine("| Time: " + DateTime.Now.ToString() + "                  |");
            WriteLine(" --------------------------------------------");
            WriteLine("");

            //Debugger.Break();

            // path to the plotConfig.xml
            WriteLine("args[0]: " + args[0]);

            string outputPath = Path.GetDirectoryName(args[0]);

            // deserialize plotConfig.xml 
            WriteLine("Reading plotConfig.xml ...");
            FieldStateConfiguration config = FieldStateConfiguration.LoadConfig(args[0]);

            // In the next 5 lines of code you find all necessary data to create
            // plot files
            string[] basePaths = config.BasePaths;
            Guid sessionGuid = config.SessionGuid;
            FieldStateConfiguration.ExportFormats exportFormat = config.ExportFormat;
            List<string> fieldNameList = config.FieldNames;
            List<TimestepNumber> timeStepList = config.TimeSteps;

            // a few debug lines
            WriteLine("BasePaths[0]      : " + basePaths[0]);
            WriteLine("SessionGuid       : " + sessionGuid.ToString());
            WriteLine("ExportFormat      : " + exportFormat.ToString());
            WriteLine("Count(FieldNames) : " + fieldNameList.Count.ToString());
            WriteLine("Count(TimeSteps)  : " + timeStepList.Count.ToString());
            WriteLine("output to         : " + outputPath);

            WriteLine("plotting...");

            PlotApplication app = new PlotApplication(config, outputPath);
            app.Run();

            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            if (mustFinalizeMPI) {
                csMPI.Raw.mpiFinalize();
            }
        }

        private static void WriteLine(string message) {
            if (rank == 0) {
                Console.WriteLine(message);
            }
        }
    }
}
