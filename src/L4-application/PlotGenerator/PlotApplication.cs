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
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.ASCIIExport;
using BoSSS.Solution.Tecplot;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace BoSSS.PlotGenerator {

    /// <summary>
    /// Actual implementation of the plot generation using the BoSSS API for
    /// convenience reasons
    /// </summary>
    public partial class PlotApplication {

        /// <summary>
        /// No Real functionality; this just forces the compiler to link the XDG dll.
        /// </summary>
        public XDGField ForceXDGLinking {
            get {
                return new XDGField(null);
            }
        }

        /// <summary>
        /// The format of the plot file
        /// </summary>
        private FieldStateConfiguration.ExportFormats m_format;

        /// <summary>
        /// The destination of the plot file
        /// </summary>
        private string m_path;

        /// <summary>
        /// De-serialized configuration options
        /// </summary>
        private FieldStateConfiguration m_config;

        /// <summary>
        /// Session information
        /// </summary>
        private ISessionInfo m_info;

        ///// <summary>
        ///// The file system driver
        ///// </summary>
        //private StandardFsDriver m_FsDriver;

        /// <summary>
        /// Handle for the plot driver
        /// </summary>
        private PlotDriver m_PlotDriver;

        /// <summary>
        /// Indicates whether jumps should be shown or whether values should be
        /// averaged at congruent nodes
        /// </summary>
        private bool m_showJumps;

        /// <summary>
        /// Indicates whether ghost overlap should be shown
        /// </summary>
        private bool m_GhostZone;

        /// <summary>
        /// Saves the <paramref name="config"/> and the <paramref name="path"/>
        /// and determines the reconstruction type (<see cref="m_showJumps"/>) and
        /// ghost level (<see cref="m_GhostZone"/>)
        /// from the <paramref name="config"/>.
        /// </summary>
        /// <param name="config">Configuration options</param>
        /// <param name="path">The output path</param>
        public PlotApplication(FieldStateConfiguration config, string path) {
            m_config = config;
            m_path = path;
            m_format = config.ExportFormat;

            switch (config.ReconstructionType) {
                case FieldStateConfiguration.ReconstructionTypes.None:
                    m_showJumps = true;
                    break;

                case FieldStateConfiguration.ReconstructionTypes.Average:
                    m_showJumps = false;
                    break;

                default:
                    throw new ApplicationException("Unknown reconstruction type");
            }

            switch (config.GhostLevel) {
                case FieldStateConfiguration.GhostLevels.One:
                    m_GhostZone = true;
                    break;

                case FieldStateConfiguration.GhostLevels.Zero:
                    m_GhostZone = false;
                    break;

                default:
                    throw new ApplicationException("Unknown ghost level");
            }

            m_Database = GetDatabase(); //init db-info
            this.GridDat = null;
        }

        private IDatabaseInfo m_Database;

        /// <summary>
        /// The driver of the database containing the data to be plotted.
        /// </summary>
        public IDatabaseDriver DBDriver {
            get {
                return m_Database.Controller.DBDriver;
            }
        }

        /// <summary>
        /// the grid
        /// </summary>
        public IGridData GridDat {
            get;
            private set;
        }

        /// <summary>
        /// Override the original run method defined by the super-class because
        /// we don't want to run an actual calculation. Plots the data for
        /// every selected time-step.
        /// </summary>
        protected void Run() {
            ISessionInfo session = m_Database.Controller.GetSessionInfo(m_config.SessionGuid);
            var timesteps = session.Timesteps.Where(
                tsi => this.m_config.TimeSteps.Contains(tsi.TimeStepNumber)).ToArray();
            int TotCnt = timesteps.Count();

            int process = this.DBDriver.MyRank + 1;
            int processCount = this.DBDriver.Size;
            this.m_info = this.m_Database.Controller.GetSessionInfo(this.m_config.SessionGuid);


            GC.Collect();
            for (int i = 0; i < timesteps.Length; i++) {
                var ts = timesteps[i];
                double physTime = ts.PhysicalTime;
                TimestepNumber timestepNo = ts.TimeStepNumber;

                if (this.GridDat == null || !this.GridDat.GridID.Equals(ts.Grid.ID)) {
                    WriteMessage(process, processCount, "Loading grid ...");
                    GridCommons grid = DBDriver.LoadGrid(ts.Grid.ID, m_Database);
                    this.GridDat = this.m_Database.Controller.GetInitializationContext(ts).GridData;
                    WriteMessage(process, processCount, "   Number of cells: " + this.GridDat.iLogicalCells.NoOfLocalUpdatedCells);
                    this.CreatePlotter();
                }


                
                WriteMessage(process, processCount, "Loading timestep ... (" + (i+1) + " of " + TotCnt + ")");
                var fields = DBDriver.LoadFields(ts, this.GridDat, this.m_config.FieldNames).ToList();
                WriteMessage(process, processCount, "Loaded timestep " + timestepNo + ". Plotting...");

                //{
                //    Console.WriteLine("computing vorticity...");
                //    DGField velX = fields.Single(f => f.Identification == "VelocityX");
                //    DGField velY = fields.Single(f => f.Identification == "VelocityY");

                //    DGField vortZ = velX.CloneAs();
                //    vortZ.Identification = "Vorticity";
                //    vortZ.Clear();
                //    vortZ.DerivativeByFlux(1.0, velY, 0);
                //    vortZ.DerivativeByFlux(-1.0, velX, 1);
                //    Console.WriteLine("done.");

                //    fields.Add(vortZ);
                //}
                PlotCurrentState(fields, physTime, timestepNo);

                double perc = Math.Round(100.0 * (double)(i+1) / (double)TotCnt, 1);
                WriteMessage(process, processCount, "Finished timestep (" + perc + "% of timesteps done)");

                // Free memory if possible
                timesteps[i] = null;
                GC.Collect();
            }
        }

        private void WriteMessage(int process, int processCount, string message) {
            Console.WriteLine("Process {0}/{1}: {2}", process, processCount, message);
        }

        /// <summary>
        /// Returns a "lightweight" object that stores information about the current database.
        /// </summary>
        /// <returns></returns>
        protected IDatabaseInfo GetDatabase() {
            IDatabaseInfo db = new DatabaseInfo(m_config.BasePaths[0]);
            return db;
        }

        /// <summary>
        /// Initializes the plot driver
        /// </summary>
        protected void CreatePlotter() {
            switch (m_format) {
                case FieldStateConfiguration.ExportFormats.TecPlot:
                    m_PlotDriver = new Tecplot(this.GridDat, m_showJumps, m_GhostZone, (uint)m_config.SuperSampling);
                    break;

                case FieldStateConfiguration.ExportFormats.Cgns:
                    throw new NotImplementedException();
                //m_PlotDriver = new CgnsExport(this.GridDat, false, m_showJumps, m_GhostZone, (uint)m_config.SuperSampling);
                //break;

                case FieldStateConfiguration.ExportFormats.Hdf5:
                    throw new NotImplementedException();
                //m_PlotDriver = new CgnsExport(this.GridDat, true, m_showJumps, m_GhostZone, (uint)m_config.SuperSampling);
                //break;

                case FieldStateConfiguration.ExportFormats.CSV:
                    m_PlotDriver = new CSVExportDriver(this.GridDat, m_showJumps, (uint)m_config.SuperSampling);
                    break;

                case FieldStateConfiguration.ExportFormats.Curve:
                    m_PlotDriver = new CurveExportDriver(this.GridDat, m_showJumps, (uint)m_config.SuperSampling);
                    break;

                case FieldStateConfiguration.ExportFormats.BinaryStream:
                    m_PlotDriver = new BinaryStreamExportDriver(this.GridDat, m_showJumps, m_GhostZone, (uint)m_config.SuperSampling);
                    break;

                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Actual plotting routine which uses
        /// <see cref="PlotDriver.PlotFields(string, double, IEnumerable{DGField})"/>
        /// to plot a single time-step
        /// </summary>
        /// <param name="fields"></param>
        /// <param name="physTime"></param>
        /// <param name="timestepNo"></param>
        protected void PlotCurrentState(IEnumerable<DGField> fields, double physTime, TimestepNumber timestepNo) {
            m_PlotDriver.PlotFields(
                Path.Combine(m_path, "state_") + timestepNo,
                physTime,
                fields);
        }
    }
}
