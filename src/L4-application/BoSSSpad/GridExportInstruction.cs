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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using BoSSS.Solution.Tecplot;
using ilPSP;
using System;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Configuration options for the plotting of a <see cref="IGridInfo"/>.
    /// </summary>
    public class GridExportInstruction : ExportInstruction {

        /// <summary>
        /// The grid to be plotted
        /// </summary>
        public IGridInfo Grid;

        /// <summary>
        /// Information about the grid.
        /// </summary>
        private GridData GridDat;

        /// <summary>
        /// Database driver
        /// </summary>
        private IDatabaseDriver dbDriver;

        /// <summary>
        /// Creates a plotting instruction with default configuration options.
        /// </summary>
        /// <param name="grid">
        /// The grid to be plotted.
        /// </param>
        public GridExportInstruction(IGridInfo grid) {
            this.Grid = grid;

            // because hacks work. (TODO)
            IFileSystemDriver fsdrv = new StandardFsDriver(grid.Database.Path);
            dbDriver = new DatabaseDriver(fsdrv);
            var gridComm = dbDriver.LoadGrid(grid.ID, grid.Database);
            GridDat = (GridData)(gridComm.iGridData);
        }

        /// <summary>
        /// Starts the export.
        /// </summary>
        public override string YouMust() {
            Console.Write("Starting export process... ");
            //// BoSSS.PlotGenerator currently does not support plots without a
            //// session...
            //FieldStateConfiguration fsConfig = new FieldStateConfiguration();
            //fsConfig.BasePaths = new string[] { Grid.Database.Path };
            //fsConfig.SessionGuid = Guid.Empty;
            //fsConfig.FieldNames = new List<string>();
            //fsConfig.TimeSteps = new List<TimestepNumber>();
            //fsConfig.SuperSampling = SuperSampling;
            //fsConfig.ReconstructionType = FieldStateConfiguration.ReconstructionTypes.None;
            //fsConfig.GhostLevel = GhostLevels;
            //fsConfig.NumberOfProcesses = NumberOfProcesses;
            //fsConfig.ExportFormat = ExportFormat;

            //string sRootpath = Utils.GetExportOutputPath();
            //string plotDirPath = Path.Combine(
            //    sRootpath, "grids", Grid.ID.ToString());

            //Perform(fsConfig, plotDirPath);

            //
            // Temporary hack
            //

            string plotDirPath = AlternativeDirectoryName;
            if (AlternativeDirectoryName == null) {
                plotDirPath = Path.Combine(
                    Utils.GetExportOutputPath(),
                    "grids",
                    Grid.ID.ToString());
            }

            // create the subdirectory if necessary
            if (!Directory.Exists(plotDirPath)) {
                Directory.CreateDirectory(plotDirPath);
            }


            // Calculate sum of edge tags for each cell (for bc debugging)
            int[] tagSum = new int[GridDat.Cells.NoOfLocalUpdatedCells];
            for (int e = 0; e < GridDat.Edges.Count; e++) {
                byte tag = GridDat.Edges.EdgeTags[e];
                if (tag != 0) {
                    Debug.Assert(GridDat.Edges.CellIndices[e, 1] < 0);

                    int cell = GridDat.Edges.CellIndices[e, 0];
                    if (cell < GridDat.Cells.NoOfLocalUpdatedCells) {
                        tagSum[cell] += tag;
                    }
                }
            }

            DGField tagSumField = new SinglePhaseField(new Basis(GridDat, 0), "edgeTagSum");
            for (int i = 0; i < tagSum.Length; i++) {
                tagSumField.SetMeanValue(i, tagSum[i]);
            }

            DGField[] partitioningFields = new DGField[GridDat.Grid.PredefinedGridPartitioning.Count];
            int fieldIndex = 0;
            foreach (var namePartitioningPair in GridDat.Grid.PredefinedGridPartitioning) {
                Partitioning currentPartitioning = GridDat.Grid.CellPartitioning;
                var cellToRankMap = dbDriver.LoadVector<int>(namePartitioningPair.Value.Guid, ref currentPartitioning);

                DGField partitioningField = new SinglePhaseField(
                    new Basis(GridDat, 0),
                    "partitioning_" + namePartitioningPair.Key);

                for (int i = 0; i < cellToRankMap.Count(); i++) {
                    partitioningField.SetMeanValue(i, cellToRankMap[i]);
                }

                partitioningFields[fieldIndex] = partitioningField;
                fieldIndex++;
            }

            Tecplot.PlotFields(
                partitioningFields.Concat(tagSumField).ToArray(),
                Path.Combine(plotDirPath, "grid"),
                0.0,
                SuperSampling);

            Console.WriteLine(" Data will be written to the following directory:");
            return plotDirPath;
        }
    }
}
