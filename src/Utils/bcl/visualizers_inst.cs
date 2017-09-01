using System;
using System.Collections.Generic;
using System.IO;

namespace bcl.visualizers_inst {

    class Program : IProgram {

        private static readonly Dictionary<string, Dictionary<string, string[]>> toolVersionFilesMap =
            new Dictionary<string, Dictionary<string, string[]>> {
            { "MatrixVisualizer",
                new Dictionary<string, string[]>() {
                    { "2015", new string[] {
                        @"MatrixVisualizer.dll", @"MatrixVisualizerVS15.dll" } },
                    { "2017", new string[] {
                        @"MatrixVisualizer.dll", @"MatrixVisualizerVS17.dll" } }
                }
            }
        };

        #region IProgram Members

        public void Execute() {
            string myDocumentsPath = Environment.GetFolderPath(
                Environment.SpecialFolder.MyDocuments);
            
            MyEnvironment env = new MyEnvironment();
            string binariesPath = Path.Combine(env.BOSSS_BIN.FullName, "Release");

            foreach (var entry in toolVersionFilesMap) {
                string tool = entry.Key;
                foreach (var versionFilesPair in entry.Value) {
                    DirectoryInfo dest = new DirectoryInfo(Path.Combine(
                        myDocumentsPath,
                        "Visual Studio " + versionFilesPair.Key,
                        "Visualizers"));
                    if (!dest.Exists) {
                        // VS version not installed
                        continue;
                    }

                    Console.WriteLine(
                        "Installing visualizer '{0}' to '{1}'",
                        tool,
                        dest.FullName);
                    foreach (string fileName in versionFilesPair.Value) {
                        FileInfo file = new FileInfo(Path.Combine(
                            binariesPath, fileName));
                        if (!file.Exists) {
                            throw new Exception(String.Format(
                                "Could not find required file '{0}' in directory '{1}'.",
                                fileName,
                                binariesPath));
                        }

                        file.CopyTo(Path.Combine(dest.FullName, fileName), true);
                    }
                }
            }
        }

        public void DecodeArgs(string[] args) {
            if (args.Length != 0) {
                throw new UserInputException("The tool visualizers-inst command has no arguments.");
            }
        }

        public void PrintUsage() {
            Console.WriteLine("Usage: bcl visualizers-inst");
            Console.WriteLine();
        }

        #endregion
    }
}
