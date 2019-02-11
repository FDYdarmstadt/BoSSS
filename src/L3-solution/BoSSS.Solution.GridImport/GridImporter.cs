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
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.GridImport.Gambit;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.GridImport {

    /// <summary>
    /// Class which contains driver routines for grid import.
    /// </summary>
    public static class GridImporter {

        /// <summary>
        /// Enum to link file type and importer class
        /// </summary>
        public enum ImporterTypes {

            /// <summary>
            /// <see cref="BoSSS.Solution.GridImport.Cgns"/>
            /// </summary>
            CGNS,

            /// <summary>
            /// <see cref="BoSSS.Solution.GridImport.Gambit.GambitNeutral"/>
            /// </summary>
            Gambit,

            /// <summary>
            /// <see cref="BoSSS.Solution.GridImport.Gmsh"/>
            /// </summary>
            Gmsh
        }

        private static readonly IDictionary<ImporterTypes, string[]> fileEndingMapping = new Dictionary<ImporterTypes, string[]>() {
            { ImporterTypes.CGNS, new[] { ".cgns", ".adf", ".hdf" } },
            { ImporterTypes.Gambit, new[] { ".neu" } },
            { ImporterTypes.Gmsh, new[] { ".msh" } }
        };

        /// <summary>
        /// Determines the type of grid file from the file ending.
        /// </summary>
        public static ImporterTypes GetImporterType(string fileName) {
            FileInfo fi = new FileInfo(fileName);
            foreach (var entry in fileEndingMapping) {
                if (entry.Value.Contains(fi.Extension.ToLowerInvariant())) {
                    return entry.Key;
                }
            }

            throw new ArgumentException("Unknown mesh file type", "fileName");
        }

        /// <summary>
        /// Loads a BoSSS grid from an grid file; the file type (see <see cref="ImporterTypes"/>) are determined by the file ending.
        /// </summary>
        public static GridCommons Import(string fileName) {
            using (var tr = new FuncTrace()) {
                int myrank;
                int size;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out myrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                ImporterTypes importerType = default(ImporterTypes);
                if (myrank == 0) {
                    importerType = GetImporterType(fileName);
                }
                importerType = importerType.MPIBroadcast(0);


                IGridImporter importer;
                {
                    tr.Info(string.Format("Loading {0} file '{1}'...", importerType.ToString(), fileName));

                    using (new BlockTrace("Import", tr)) {

                        switch (importerType) {
                            case ImporterTypes.Gambit:
                                if (size > 1)
                                    throw new NotSupportedException("Not supported in parallel mode");
                                GambitNeutral gn = new GambitNeutral(fileName);
                                if (gn.BoSSSConversionNeccessary()) {
                                    gn = gn.ToLinearElements();
                                }

                                importer = gn;
                                break;

                            case ImporterTypes.CGNS:
                                if (size > 1)
                                    throw new NotSupportedException("Not supported in parallel mode");
                                importer = new Cgns(fileName);
                                break;

                            case ImporterTypes.Gmsh:
                                importer = new Gmsh(fileName);
                                break;

                            default:
                                throw new NotImplementedException();
                        }
                    }

                    tr.Info("Converting to BoSSS grid ...");
                }

                GridCommons grid;
                using (new BlockTrace("Conversion", tr)) {
                    grid = importer.GenerateBoSSSGrid();
                }

                return grid;
            }
        }
    }
}