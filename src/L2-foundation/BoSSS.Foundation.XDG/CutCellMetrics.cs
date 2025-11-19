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

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.Intersecting;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using IntersectingQuadrature;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {


    public interface ICutCellMetrics {
        /// <summary>
        /// Quadrature order used to obtain the metrics
        /// </summary>
        int CutCellQuadratureOrder { get; }

        /// <summary>
        /// available species
        /// </summary>
        IEnumerable<SpeciesId> SpeciesList { get; }


        /// <summary>
        /// Volume of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        Dictionary<SpeciesId, MultidimensionalArray> CutCellVolumes { get; }


        /// <summary>
        /// Surface (Level-Set plus cut edges) of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        Dictionary<SpeciesId, MultidimensionalArray> CellSurface { get; }
    }


    /// <summary>
    /// HMF-based computation of cut-cell volumes and cut-edge areas.
    /// </summary>
    public class CutCellMetrics : ICutCellMetrics {

        /// <summary>
        /// owner object.
        /// </summary>
        public XDGSpaceMetrics XDGSpaceMetrics {
            get;
            private set;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        internal CutCellMetrics(XDGSpaceMetrics owner) {
            XDGSpaceMetrics = owner;
            ComputeNonAgglomeratedMetrics();
        }

        /// <summary>
        /// The quadrature order used for computing cell volumes and edge areas.
        /// </summary>
        public int CutCellQuadratureOrder {
            get {
                return XDGSpaceMetrics.CutCellQuadOrder;
            }
        }

        /// <summary>
        /// The kind of HMF which is be used for computing cell volumes.
        /// </summary>
        public CutCellQuadratureMethod HMFvariant {
            get {
                return XDGSpaceMetrics.CutCellQuadratureType;
            }
        }

        /// <summary>
        /// All species for which agglomeration is available.
        /// </summary>
        public IEnumerable<SpeciesId> SpeciesList {
            get {
                return XDGSpaceMetrics.SpeciesList;
            }
        }

        /// <summary>
        /// Volume of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutCellVolumes {
            get;
            private set;
        }

        /// <summary>
        /// Area of the interface, resp. level-set-surface of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> InterfaceArea {
            get;
            private set;
        }

        /// <summary>
        /// Area of non-agglomerated cut edges.
        /// - key: species
        /// - index: edge index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutEdgeAreas {
            get;
            private set;
        }

        /// <summary>
        /// Surface (Level-Set plus cut edges) of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CellSurface {
            get;
            private set;
        }

        /// <summary>
        /// Writes all the edge rules as vtp files
        /// </summary>
        public void WriteEdgeRulesToVtp() {
            var schH = new XQuadSchemeHelper(XDGSpaceMetrics);
            var gd = XDGSpaceMetrics.GridDat;
            SpeciesId[] species = this.SpeciesList.ToArray();

            for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                SpeciesId spc = species[iSpc];

                var edgeScheme = schH.GetEdgeQuadScheme(spc);
                var chunRulePairList = edgeScheme.Compile(gd, this.CutCellQuadratureOrder);

                chunRulePairList.ToVtpFilesEdge(gd, "edgeQuadFor" + spc.ToString(XDGSpaceMetrics.Tracker) + HMFvariant);
            }
        }

        /// <summary>
        /// Writes all the volume rules as vtp files
        /// </summary>
        public void WriteVolumeRulesToVtp() {
            var schH = new XQuadSchemeHelper(XDGSpaceMetrics);
            var gd = XDGSpaceMetrics.GridDat;
            SpeciesId[] species = this.SpeciesList.ToArray();

            for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                SpeciesId spc = species[iSpc];

                var volScheme = schH.GetVolumeQuadScheme(spc);
                var chunRulePairList = volScheme.Compile(gd, this.CutCellQuadratureOrder);

                chunRulePairList.ToVtpFilesCell(gd, "volQuadFor" + spc.ToString(XDGSpaceMetrics.Tracker) + HMFvariant);
            }
        }

        /// <summary>
        /// Writes all the surface rules (level set=0) as vtp files
        /// </summary>
        public void WriteSurfaceRulesToVtp() {
            var schH = new XQuadSchemeHelper(XDGSpaceMetrics);
            var gd = XDGSpaceMetrics.GridDat;
            SpeciesId[] species = this.SpeciesList.ToArray();

            if(species.Length > 0) {
                var AllSpc = XDGSpaceMetrics.TotalSpeciesList;
                var requiredSpecies = XDGSpaceMetrics.SpeciesList;

                // loop over all possible pairs of species
                for(int iSpcA = 0; iSpcA < AllSpc.Count - 1; iSpcA++) {
                    var SpeciesA = AllSpc[iSpcA];
                    var SpeciesADom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesA);
                    int iLocalSpcA = requiredSpecies.IndexOf(SpeciesA);

                    // Standard XDG case, where level sets are always
                    // an interface between species
                    for(int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                        var SpeciesB = AllSpc[iSpcB];
                        int iLocalSpcB = requiredSpecies.IndexOf(SpeciesB);

                        if(iLocalSpcA > -1 || iLocalSpcB > -1) {
                            var SpeciesBDom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesB);
                            var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);
                            int NoOfLs = XDGSpaceMetrics.NoOfLevelSets;
                            for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                                if(schH.SpeciesAreSeparatedByLevSet(iLevSet, SpeciesA, SpeciesB)) {
                                    var LsDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
                                    var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                    CellQuadratureScheme SurfIntegration = schH.GetLevelSetQuadScheme(iLevSet, SpeciesA, IntegrationDom);
                                    var chunRulePairList = SurfIntegration.Compile(gd, this.CutCellQuadratureOrder);
                                    chunRulePairList.ToVtpFilesCell(gd, "surfQuadFor" + SpeciesA.ToString(XDGSpaceMetrics.Tracker) + "-" + SpeciesB.ToString(XDGSpaceMetrics.Tracker) + HMFvariant);
                                }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Computes Cell-volumes and edge areas before agglomeration.
        /// </summary>
        void ComputeNonAgglomeratedMetrics() {
            using(var tr = new FuncTrace()) {
                MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD);

                var gd = XDGSpaceMetrics.GridDat;
                int JE = gd.iLogicalCells.Count;
                int J = gd.iLogicalCells.NoOfLocalUpdatedCells;
                int D = gd.SpatialDimension;
                int EE = gd.Edges.Count;
                SpeciesId[] species = this.SpeciesList.ToArray();
                int NoSpc = species.Count();
                int[,] E2C = gd.iLogicalEdges.CellIndices;

                var schH = new XQuadSchemeHelper(XDGSpaceMetrics);
                //var schH = this.XDGSpaceMetrics.XQuadSchemeHelper;

                // collect all per-cell-metrics in the same MultidimArry, for MPI-exchange (only 1 exchange instead of three, saving some overhead)
                // 1st index: cell
                // 2nd index: species
                // 3rd index: 0 for interface surface per cell, 1 for cut-cell-volume, 2 for cut-cell surface
                double[] vec_cellMetrics = new double[JE * NoSpc * 3];
                MultidimensionalArray cellMetrics = MultidimensionalArray.CreateWrapper(vec_cellMetrics, JE, NoSpc, 3);

                this.CutEdgeAreas = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CutCellVolumes = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.InterfaceArea = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CellSurface = new Dictionary<SpeciesId, MultidimensionalArray>();

                //var schS = new List<CellQuadratureScheme>();
                //var rulz = new List<ICompositeQuadRule<QuadRule>>();

                tr.Info("Checkpoint1");

                // edges and volumes
                // =================
                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var cellVol = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1);
                    SpeciesId spc = species[iSpc];

                    // compute cut edge area
                    // ---------------------
                    MultidimensionalArray edgArea = MultidimensionalArray.Create(EE);
                    this.CutEdgeAreas.Add(spc, edgArea);

                    var edgeScheme = schH.GetEdgeQuadScheme(spc);
                    var edgeRule = edgeScheme.Compile(gd, this.CutCellQuadratureOrder);

                    BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                        [1], gd,
                        edgeRule,
                        _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                        {
                            EvalResult.SetAll(1.0);
                        },
                        _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                        {
                            for(int i = 0; i < Length; i++) {
                                int iEdge = i + i0;
                                Debug.Assert(edgArea[iEdge] == 0);
                                edgArea[iEdge] = ResultsOfIntegration[i, 0];
                                Debug.Assert(!(double.IsNaN(edgArea[iEdge]) || double.IsInfinity(edgArea[iEdge])));
                            }
                        }).Execute();

                    // sum up edges for surface
                    // ------------------------

                    tr.Info("Checkpoint1.1");

                    var cellSurf = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 2);


                    for(int e = 0; e < EE; e++) {

                        double a = edgArea[e];
                        int jCell0 = E2C[e, 0];
                        int jCell2 = E2C[e, 1];
                        cellSurf[jCell0] += a;
                        if(jCell2 >= 0)
                            cellSurf[jCell2] += a;

                    }

                    tr.Info("Checkpoint1.2");

                    // compute cut cell volumes
                    // ------------------------

                    tr.Info("Species: " + XDGSpaceMetrics.Tracker.GetSpeciesName(spc));
                    var volScheme = schH.GetVolumeQuadScheme(spc);
                    tr.Info("Checkpoint1.3");
                    var volRule = volScheme.Compile(gd, this.CutCellQuadratureOrder);
                    tr.Info("Checkpoint1.4");

                    BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                        new int[] { 1 }, gd,
                        volRule,
                        _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                        {
                            EvalResult.SetAll(1.0);
                        },
                        _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                        {
                            for(int i = 0; i < Length; i++) {
                                int jCell = i + i0;
                                Debug.Assert(cellVol[jCell] == 0);
                                cellVol[jCell] = ResultsOfIntegration[i, 0];
                                //Console.WriteLine("cellVol " + cellVol[jCell] + " of cell " + jCell);
                                Debug.Assert(!(double.IsNaN(cellVol[jCell]) || double.IsInfinity(cellVol[jCell])));
                            }
                        }).Execute();
                }

                tr.Info("Checkpoint2");

                // interface surface
                // =================

                // loop over surfaces
                if(species.Length > 0) {
                    var AllSpc = XDGSpaceMetrics.TotalSpeciesList;
                    var requiredSpecies = XDGSpaceMetrics.SpeciesList;

                    // loop over all possible pairs of species
                    for(int iSpcA = 0; iSpcA < AllSpc.Count - 1; iSpcA++) {
                        var SpeciesA = AllSpc[iSpcA];
                        var SpeciesADom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesA);
                        int iLocalSpcA = requiredSpecies.IndexOf(SpeciesA);

                        // Standard XDG case, where level sets are always
                        // an interface between species
                        for(int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                            var SpeciesB = AllSpc[iSpcB];
                            int iLocalSpcB = requiredSpecies.IndexOf(SpeciesB);

                            if(iLocalSpcA > -1 || iLocalSpcB > -1) {
                                var SpeciesBDom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesB);
                                var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);
                                int NoOfLs = XDGSpaceMetrics.NoOfLevelSets;
                                for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                                    if(schH.SpeciesAreSeparatedByLevSet(iLevSet, SpeciesA, SpeciesB)) {
                                        var LsDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
                                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                        CellQuadratureScheme SurfIntegration = schH.GetLevelSetQuadScheme(iLevSet, SpeciesA, IntegrationDom);
                                        var rule = SurfIntegration.Compile(gd, this.CutCellQuadratureOrder);

                                        BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                                            [1], gd,
                                            rule,
                                            _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                                                EvalResult.SetAll(1.0);
                                            },
                                            _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                for(int i = 0; i < Length; i++) {
                                                    int jCell = i + i0;
                                                    if(iLocalSpcA > -1) {
                                                        cellMetrics[jCell, iLocalSpcA, 0] += ResultsOfIntegration[i, 0];
                                                        Debug.Assert(!(double.IsNaN(cellMetrics[jCell, iLocalSpcA, 0]) || double.IsInfinity(cellMetrics[jCell, iLocalSpcA, 0])));
                                                    }
                                                    if(iLocalSpcB > -1) {
                                                        cellMetrics[jCell, iLocalSpcB, 0] += ResultsOfIntegration[i, 0];
                                                        Debug.Assert(!(double.IsNaN(cellMetrics[jCell, iLocalSpcB, 0]) || double.IsInfinity(cellMetrics[jCell, iLocalSpcB, 0])));
                                                    }
                                                }
                                            }).Execute();
                                    }
                                }
                            }
                        }
                    }
                }

                tr.Info("Checkpoint3");

                // MPI exchange & store
                // ====================

                if(species.Length > 0) {
#if DEBUG
                    int NoOfSpc = species.Length;
                    var cellMetricsB4 = cellMetrics.ExtractSubArrayShallow(new[] { 0, 0, 0 }, new[] { J - 1, NoOfSpc - 1, 1 }).CloneAs();
#endif
                    vec_cellMetrics.MPIExchange(gd);
#if DEBUG
                    var cellMetricsComp = cellMetrics.ExtractSubArrayShallow(new[] { 0, 0, 0 }, new[] { J - 1, NoOfSpc - 1, 1 }).CloneAs();
                    cellMetricsComp.Acc(-1.0, cellMetricsB4);
                    Debug.Assert(cellMetricsComp.L2Norm() == 0.0);
#endif
                }
                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var spc = species[iSpc];
                    this.InterfaceArea.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 0).CloneAs());
                    this.CutCellVolumes.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1).CloneAs());
                    this.CellSurface.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 2).CloneAs());

                    this.CellSurface[spc].Acc(1.0, this.InterfaceArea[spc]);
                }
            }
        }

        /// <summary>
        /// Total length (i.e., D-2 dimensional measure) of cut line (i.e., intersection of D-1 dimensional level-set surface with D-1 dimensional cell boundary) for each non-agglomerated cut cell.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutLineLength {
            get {
                if(m_CutLineLength == null) {
                    ComputeNonAgglomeratedCutLineMetrics();
                }
                return m_CutLineLength;
            }
        }

        Dictionary<SpeciesId, MultidimensionalArray> m_CutLineLength;

        /// <summary>
        /// Total length (i.e., D-2 dimensional measure) of 
        /// level-set intersection (i.e., intersection of the D-1 dimensional level-set 0 with the d-1 dimensional level-set 1) 
        /// for each non-agglomerated cut cell.
        /// - key: species
        /// - index: cell index
        /// See also: <see cref="XQuadSchemeHelper.GetContactLineQuadScheme"/>
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> IntersectionLength {
            get {
                if(m_IntersectionLength == null) {
                    ComputeNonAgglomeratedCutLineMetrics();
                }
                return m_IntersectionLength;
            }
        }

        Dictionary<SpeciesId, MultidimensionalArray> m_IntersectionLength;


        /// <summary>
        /// Total length (i.e., D-2 dimensional measure) of cut line (i.e., intersection of D-1 dimensional level-set surface with D-1 dimensional cell boundary) for each edge.
        /// for each edge.
        /// - key: species
        /// - index: edge index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutLineLengthEdge {
            get {
                if(m_CutLineLengthEdge == null) {
                    ComputeNonAgglomeratedCutLineMetrics();
                }
                return m_CutLineLengthEdge;
            }
        }

        Dictionary<SpeciesId, MultidimensionalArray> m_CutLineLengthEdge;



        void ComputeNonAgglomeratedCutLineMetrics() {
            using(var tr = new FuncTrace()) {
                MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD);

                var gd = XDGSpaceMetrics.GridDat;
                int JE = gd.iLogicalCells.Count;
                int J = gd.iLogicalCells.NoOfLocalUpdatedCells;
                int D = gd.SpatialDimension;
                int EE = gd.Edges.Count;
                SpeciesId[] species = this.SpeciesList.ToArray();
                int NoSpc = species.Count();
                int NoOfLevelSets = this.XDGSpaceMetrics.NoOfLevelSets;
                int[,] E2C = gd.iLogicalEdges.CellIndices;

                //var schH = new XQuadSchemeHelper(XDGSpaceMetrics);
                var schH = this.XDGSpaceMetrics.XQuadSchemeHelper;

                // collect all per-cell-metrics in the same MultidimArry, for MPI-exchange (only 1 exchange instead of three, saving some overhead)
                // 1st index: cell
                // 2nd index: species
                double[] vec_cellMetrics = new double[JE * NoSpc * 2];
                MultidimensionalArray cellMetrics = MultidimensionalArray.CreateWrapper(vec_cellMetrics, JE, NoSpc, 2);

                this.m_CutLineLength = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.m_IntersectionLength = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.m_CutLineLengthEdge = new Dictionary<SpeciesId, MultidimensionalArray>();

                tr.Info("Checkpoint1");

                // edges
                // =====


                // compute cut line measure
                // ------------------------
                if(species.Length > 0) {
                    var AllSpc = XDGSpaceMetrics.TotalSpeciesList;
                    var requiredSpecies = XDGSpaceMetrics.SpeciesList;

                    // loop over all possible pairs of species, pt. 1 (loop over first species)...
                    for(int iSpcA = 0; iSpcA < AllSpc.Count - 1; iSpcA++) {
                        var SpeciesA = AllSpc[iSpcA];
                        var SpeciesADom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesA);
                        int iLocalSpcA = requiredSpecies.IndexOf(SpeciesA);

                        // loop over all possible pairs of species, pt. 2 (loop over second species)....
                        for(int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                            var SpeciesB = AllSpc[iSpcB];
                            int iLocalSpcB = requiredSpecies.IndexOf(SpeciesB);

                            if(iLocalSpcA > -1 || iLocalSpcB > -1) {
                                var SpeciesBDom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesB);
                                var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);
                                int NoOfLs = XDGSpaceMetrics.NoOfLevelSets;

                                // loop over all level sets...
                                for(int iLevelSet = 0; iLevelSet < NoOfLs; iLevelSet++) {
                                    if(schH.SpeciesAreSeparatedByLevSet(iLevelSet, SpeciesA, SpeciesB)) {
                                        var LsDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevelSet);
                                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                        var cutLineScheme = schH.Get_SurfaceElement_EdgeQuadScheme(SpeciesA, SpeciesB, iLevelSet);
                                        var cutLineRule = cutLineScheme.Compile(gd, this.CutCellQuadratureOrder);

                                        MultidimensionalArray edgMeas = MultidimensionalArray.Create(EE);
                                        BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                                            [1], gd,
                                            cutLineRule,
                                            _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                                            {
                                                EvalResult.SetAll(1.0);
                                            },
                                            _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                                            {
                                                for(int i = 0; i < Length; i++) {
                                                    int iEdge = i + i0;
                                                    Debug.Assert(edgMeas[iEdge] == 0);
                                                    edgMeas[iEdge] = ResultsOfIntegration[i, 0];
                                                    Debug.Assert(!(double.IsNaN(edgMeas[iEdge]) || double.IsInfinity(edgMeas[iEdge])));
                                                }
                                            }).Execute();

                                        foreach(var spc in new SpeciesId[] { SpeciesA, SpeciesB }) {
                                            string spc_name = XDGSpaceMetrics.Tracker.GetSpeciesName(spc);
                                            string spc_nameA = XDGSpaceMetrics.Tracker.GetSpeciesName(SpeciesA);
                                            string spc_nameB = XDGSpaceMetrics.Tracker.GetSpeciesName(SpeciesB);

                                            if(!m_CutLineLengthEdge.ContainsKey(spc)) {
                                                m_CutLineLengthEdge.Add(spc, MultidimensionalArray.Create(EE));
                                            }

                                            m_CutLineLengthEdge[spc].Acc(1.0, edgMeas);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    SpeciesId spc = species[iSpc];

                    // sum up edges for surface
                    // ------------------------

                    tr.Info("Checkpoint1.1");
                    {
                        var edgMeas = m_CutLineLengthEdge[spc];
                        var cutLineInCell = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 0);

                        for(int e = 0; e < EE; e++) {
                            double a = edgMeas[e];
                            int jCell0 = E2C[e, 0];
                            int jCell2 = E2C[e, 1];
                            cutLineInCell[jCell0] += a;
                            if(jCell2 >= 0)
                                cutLineInCell[jCell2] += a;

                        }
                    }

                    // compute intersection line (between level-set 0 and 1)
                    // -----------------------------------------------------

                    var intersectLine = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1);

                    if(this.XDGSpaceMetrics.NoOfLevelSets > 1) {
                        var intersectScheme = schH.GetContactLineQuadScheme(spc, 0, 1);
                        var intersectRule = intersectScheme.Compile(gd, this.CutCellQuadratureOrder);

                        BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                            [1], gd,
                            intersectRule,
                            _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                            {
                                EvalResult.SetAll(1.0);
                            },
                            _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                            {
                                for(int i = 0; i < Length; i++) {
                                    int jCell = i + i0;
                                    Debug.Assert(intersectLine[jCell] == 0);
                                    intersectLine[jCell] = ResultsOfIntegration[i, 0];
                                    //Console.WriteLine("cellVol " + cellVol[jCell] + " of cell " + jCell);
                                    Debug.Assert(!(double.IsNaN(intersectLine[jCell]) || double.IsInfinity(intersectLine[jCell])));
                                }
                            }).Execute();
                    } else {
                        intersectLine.Clear();
                    }

                    tr.Info("Checkpoint2");
                }


                // MPI exchange & store
                // ====================

                if(species.Length > 0) {
#if DEBUG
                    int NoOfSpc = species.Length;
                    var cellMetricsB4 = cellMetrics.ExtractSubArrayShallow(new[] { 0, 0, 0 }, new[] { J - 1, NoOfSpc - 1, 1 }).CloneAs();
#endif
                    vec_cellMetrics.MPIExchange(gd);
#if DEBUG
                    var cellMetricsComp = cellMetrics.ExtractSubArrayShallow(new[] { 0, 0, 0 }, new[] { J - 1, NoOfSpc - 1, 1 }).CloneAs();
                    cellMetricsComp.Acc(-1.0, cellMetricsB4);
                    Debug.Assert(cellMetricsComp.L2Norm() == 0.0);
#endif
                }
                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var spc = species[iSpc];
                    this.m_CutLineLength.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 0).CloneAs());
                    this.m_IntersectionLength.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1).CloneAs());
                }
            }
        }
    }
}

