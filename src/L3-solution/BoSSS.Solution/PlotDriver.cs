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
using System.Runtime.InteropServices;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution {

    /// <summary>
    /// Abstract base class for all drivers implementing a specific output
    /// format for plotting purposes.
    /// </summary>
    abstract public class PlotDriver {

        /// <summary>
        /// Constructs a new PlotDriver
        /// </summary>
        /// <param name="context">The omnipresent context</param>
        /// <param name="showJumps">
        /// If true, ghost cell values will be calculated and additional ghost
        /// zone information will be exported (i.e. ghost cells will be
        /// included in the plots)
        /// </param>
        /// <param name="showGhostCells">
        /// If true, the real DG data (including discontinuities) should be
        /// exported. Otherwise, the mean value of all values in one node
        /// should be calculated
        /// </param>
        /// <param name="superSampling">
        /// how often one computational cell should be subdivided;
        /// <see cref="ZoneDriver.superSampling"/>
        /// </param>
        /// <param name="__sgrd">
        /// </param>
        protected PlotDriver(GridData context, bool showJumps, bool showGhostCells, uint superSampling, SubGrid __sgrd) {

            var RefElms = context.Grid.RefElements;
            ZoneDrivers = new ZoneDriver[RefElms.Length];
            for (int iKref = 0; iKref < RefElms.Length; iKref++) {
                SubGrid ZoneSgrd;
                if (__sgrd == null) {
                    ZoneSgrd = context.GetRefElementSubGrid(iKref);
                } else {
                    var B = context.GetRefElementSubGrid(iKref);
                    ZoneSgrd = new SubGrid(B.VolumeMask.Intersect(__sgrd.VolumeMask));
                }

                if (ZoneSgrd.GlobalNoOfCells <= 0)
                    continue;

                //Debug.Assert(false, "break: " + context.MyRank); 
                ZoneDrivers[iKref] = CreateZoneDriver(context, iKref, showJumps, showGhostCells, superSampling, ZoneSgrd);
            }

            this.gridData = context;
        }

        /// <summary>
        /// %
        /// </summary>
        protected GridData gridData;

        /// <summary>
        /// %
        /// </summary>
        abstract protected ZoneDriver CreateZoneDriver(GridData context, int iKref, bool showJumps, bool showGhostCells, uint superSampling, SubGrid __sgrd);

        /// <summary>
        /// index: correlates with reference element index
        /// </summary>
        ZoneDriver[] ZoneDrivers;

        /// <summary>
        /// opens an output file
        /// </summary>
        abstract protected void OpenFile(string filename, IEnumerable<string> fieldNames);

        /// <summary>
        /// Closes the currently open output file 
        /// </summary>
        abstract protected void CloseFile();

        /// <summary>
        /// plots one zone, i.e. all elements that have the same cell type, i.e. the same reference element
        /// </summary>
        abstract public class ZoneDriver {

            /// <summary>
            /// Border vertices;
            /// </summary>
            private List<int> borderVertices = new List<int>();

            /// <summary>
            /// Subdivisions per cell;
            /// </summary>
            protected int subdivisionsPerCell;

            /// <summary>
            /// Number of all external cells;
            /// </summary>
            private int globalNoOfExternalCells;

            /// <summary>
            /// Vertices already processed?;
            /// </summary>
            private bool[] verticesAlreadyProcessed;

            /// <summary>
            /// Sampling modified indices into the vertices-array (<see cref="vertices"/>);
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>
            ///   1st index: local cell index.
            ///   </item>
            ///   <item>
            ///   2nd index: vertex index for the <see cref="RefElement"/>;
            ///   </item>
            /// </list>
            /// </remarks>
            protected int[,] connectivity;

            /// <summary>
            /// Inter process generated smoothed field values
            /// </summary>
            protected double[] smoothedResult;

            /// <summary>
            /// Generated jumped field values
            /// </summary>
            protected double[] notSmoothedResult;

            /// <summary>
            /// Inter process shared cell indices
            /// </summary>
            private long[] sharedCellIndices;

            /// <summary>
            /// Dimension of the grid
            /// </summary>
            protected int dimension;

            /// <summary>
            /// Number of super sampled local cells
            /// </summary>
            private int NewNoOfLocalCells;

            /// <summary>
            /// Number of super sampled cells
            /// </summary>
            private int NewNoOfCells;

            /// <summary>
            /// Inter process shared vertices multiplicity
            /// </summary>
            private byte[,] sharedVerticesMultiplicity;

            /// <summary>
            /// Maximum NoOfCells (see <see cref="NoOfCells"/>) - NoOfLocalCells
            /// </summary>
            private int maxLocalNoOfExternalCells;


            /// <summary>
            /// total number of cells (locally updated plus external)
            /// </summary>
            protected int NoOfCells;

            /// <summary>
            /// total number of used cells
            /// </summary>
            protected int totalCells;

            /// <summary>
            /// total number of used vertices
            /// </summary>
            protected int totalVertices;

            /// <summary>
            /// Number of used cells without ghost cells
            /// </summary>
            protected int innerCells;

            /// <summary>
            /// Number of used vertices without ghost cells
            /// </summary>
            protected int innerVertices;

            /// <summary>
            /// Number of used ghost cells
            /// </summary>
            protected int ghostCells;

            /// <summary>
            /// Number of used ghost vertices
            /// </summary>
            protected int ghostVertices;

            /// <summary>
            /// The omnipresent context
            /// </summary>
            protected GridData context;

            /// <summary>
            /// If true, the real DG data (including discontinuities) should be
            /// exported. Otherwise, the mean value of all values in one node
            /// should be calculated
            /// </summary>
            protected bool showJumps;

            /// <summary>
            /// If true, mesh overlap is calculated
            /// </summary>
            protected bool ghostZone;

            /// <summary>
            /// Controls the number of subdivisions of a single grid cell
            /// (<see cref="RefElement.SubdivisionTreeNode"/>) in the output. If 0, the original
            /// grid should be exported
            /// </summary>
            /// <example>
            /// Let A be a grid <see cref="RefElement"/> with N vertices. Then, for
            /// superSampling=0, any sub-class will export N node-centered values
            /// for A. For superSampling=1, A will be subdivided into N sub-cells
            /// which leads to total of N^N node-centered values for A. Thus, for
            /// superSampling=M, N^M values per cell will be evaluated.
            /// Note:
            /// If <see cref="showJumps"/>=true, exactly N^M values will be
            /// exported. On the other hand, if <see cref="showJumps"/>=false, only
            /// only values at the nodes of A plus values at distinct nodes of the
            /// children should exported. E.g. in case of a square (-> N=4) with
            /// showJumps=false and superSampling=1 (-> M=1), 9 values will be
            /// exported (4 original vertices + 4 new vertices on the edges + 1 new
            /// vertex in the cell center)
            /// </example>
            protected uint superSampling;

            /// <summary>
            /// Global vertices coordinates. See
            /// <see cref="GridData.TransformLocal2Global(MultidimensionalArray, MultidimensionalArray, int)"/>
            /// for the definition of the indices (parameter
            /// "GlobalVerticesOut").
            /// </summary>
            protected MultidimensionalArray verticeCoordinates;

            /// <summary>
            /// indices into the multiplicity array;
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>
            ///   1st index: local cell index (including local ghost cell index)
            ///   </item>
            /// </list>
            /// </remarks>
            private byte[] multiplicity;

            /// <summary>
            /// indices into the number of external cells;
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>
            ///   1st index: process index
            ///   </item>
            /// </list>
            /// </remarks>
            private int[] noOfExternalCells;

            /// <summary>
            /// indices into the <see cref="vertices"/>-array;
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>
            ///   1st index: local cell index.
            ///   </item>
            ///   <item>
            ///   2nd index: vertex index for the <see cref="RefElement"/>;
            ///   </item>
            /// </list>
            /// </remarks>
            private int[,] cellVertices;

            /// <summary>
            /// all vertices/nodes of the gird;
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>1st index: vertex index</item>
            ///   <item>2nd index: spatial dimension</item>
            /// </list>
            /// </remarks>
            protected double[,] vertices;

            /// <summary>
            /// Leaves of the subdivision tree which hold the results of the
            /// <see cref="superSampling"/> divisions of the grid cells.
            /// </summary>
            protected RefElement.SubdivisionTreeNode[] subdivisionTreeLeaves;

            /// <summary>
            /// The number of vertices one cell consists of (depends on the grid
            /// type AND on <see cref="superSampling"/>)
            /// </summary>
            protected int verticesPerCell;

            /// <summary>
            /// The local vertices' coordinates
            /// <see cref="RefElement.SubdivisionTreeNode.GlobalVertice"/>
            /// </summary>
            private NodeSet localVerticeCoordinates;

            /// <summary>
            /// If running parallel
            /// </summary>
            private bool parallelMode;

            /// <summary>
            /// the part of the computational grid that is plotted
            /// </summary>
            protected SubGrid sgrd;

            /// <summary>
            /// the reference element (element type) for the zone
            /// </summary>
            protected RefElement Zone_Element;

            /// <summary>
            /// Constructs a new PlotDriver
            /// </summary>
            /// <param name="context">The omnipresent context</param>
            /// <param name="showJumps">
            /// If true, ghost cell values will be calculated and additional ghost
            /// zone information will be exported (i.e. ghost cells will be
            /// included in the plots)
            /// </param>
            /// <param name="showGhostCells">
            /// If true, the real DG data (including discontinuities) should be
            /// exported. Otherwise, the mean value of all values in one node
            /// should be calculated
            /// </param>
            /// <param name="superSampling">
            /// how often one computational cell should be subdivided;
            /// <see cref="superSampling"/>
            /// </param>
            /// <param name="__sgrd">
            /// </param>
            /// <param name="iKref">
            /// reference element index, <see cref="GridCommons.GetRefElement"/>;
            /// </param>
            protected ZoneDriver(GridData context, int iKref, bool showJumps, bool showGhostCells, uint superSampling, SubGrid __sgrd) {
                using (new FuncTrace()) {

                    // record args
                    // ===========
                    this.context = context;
                    this.showJumps = showJumps;
                    this.ghostZone = showGhostCells;
                    this.superSampling = superSampling;
                    if (__sgrd == null)
                        // safe some if's in subsequence
                        throw new ArgumentNullException();
                    this.sgrd = __sgrd;
                    this.Zone_Element = context.Grid.GetRefElement(iKref);

                    // default values
                    // ==============

                    RefElement.SubdivisionTreeNode subdiv = this.Zone_Element.GetSubdivisionTree((int)superSampling);
                    subdivisionTreeLeaves = subdiv.GetLeaves();
                    verticesPerCell = subdiv.GlobalVertice.GetLength(0);

                    int NoOfLocalCells = sgrd.LocalNoOfCells;
                    NoOfCells = NoOfLocalCells + sgrd.NoOfGhostCells;

                    subdivisionsPerCell = subdivisionTreeLeaves.Length;
                    NewNoOfLocalCells = NoOfLocalCells * subdivisionsPerCell;
                    NewNoOfCells = NoOfCells * subdivisionsPerCell;
                    dimension = context.Grid.SpatialDimension;

                    // vertices per cell, local coordinates (from the leaves of the subdivision)
                    localVerticeCoordinates = new NodeSet(this.Zone_Element, subdiv.GlobalVertice);

                    // create vertices of the grid
                    // ===========================


                    parallelMode = NoOfCells != NoOfLocalCells;
                    if (parallelMode && !showJumps) {
                        totalCells = NewNoOfCells;
                        verticeCoordinates = MultidimensionalArray.Create(NoOfCells, verticesPerCell, dimension);
                        //context.GridDat.TransformLocal2Global(localVerticeCoordinates, verticeCoordinates, 0, NoOfCells, 0);
                    } else {
                        totalCells = NewNoOfLocalCells;
                        NoOfCells = NoOfLocalCells;
                        NewNoOfCells = NewNoOfLocalCells;
                        verticeCoordinates = MultidimensionalArray.Create(context.Cells.Count, verticesPerCell, dimension);
                        //context.GridDat.TransformLocal2Global(localVerticeCoordinates, verticeCoordinates, 0, NoOfLocalCells, 0);
                    }

                    double[] hmin = new double[verticeCoordinates.GetLength(0)];
                    int cnt = 0;
                    foreach (var cnk in this.sgrd.VolumeMask.GetEnumerableWithExternal()) {
                        context.TransformLocal2Global(localVerticeCoordinates, cnk.i0, cnk.Len,
                            verticeCoordinates, cnt);

                        for (int j = cnk.i0; j < cnk.JE; j++) {
                            hmin[cnt] = context.Cells.h_min[j];
                            cnt++;
                        }
                    }


                    innerCells = NewNoOfLocalCells;
                    ghostCells = NewNoOfCells - NewNoOfLocalCells;


                    //if (parallelMode || superSampling > 0) {
                    //    // eliminate duplicate vertices, build connectivity
                    //    // ++++++++++++++++++++++++++++++++++++++++++++++++

                    InitializeVertice2(
                        verticeCoordinates,
                        out cellVertices,
                        out vertices,
                        hmin);

                    //(GridData grdDat, out int[][] cellVertices, out MultidimensionalArray vertice, double[] h)
                    //} else {
                    //    // grid vertices and cell connectivity can be taken directly from BoSSS
                    //    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    //    cellVertices = context.GridDat.CellVertices;
                    //    vertices = context.GridDat.Vertice;
                    //}


                    if (superSampling > 0) {
                        connectivity = new int[totalCells, this.Zone_Element.NoOfVertices];
                        for (int i = 0; i < NoOfLocalCells; i++) {
                            for (int j = 0; j < subdivisionsPerCell; j++) {
                                int cell = i * subdivisionsPerCell + j;
                                int[] indices = subdivisionTreeLeaves[j].GlobalVerticeInd;

                                for (int k = 0; k < indices.Length; k++) {
                                    connectivity[cell, k] = cellVertices[i, indices[k]];
                                }
                            }
                        }
                    } else {
                        connectivity = cellVertices;
                    }

                    innerVertices = totalVertices = vertices.GetLength(0);

                    if (showJumps && parallelMode) {

                        bool[] outerBorderVertices = new bool[totalVertices];

                        for (int cell = NoOfLocalCells; cell < NoOfCells; cell++) {
                            for (int vertexInCell = 0; vertexInCell < verticesPerCell; vertexInCell++) {
                                int vertex = cellVertices[cell, vertexInCell];

                                if (!outerBorderVertices[vertex]) {
                                    borderVertices.Add(vertex);
                                    outerBorderVertices[vertex] = true;
                                }
                            }
                        }
                        ghostVertices = borderVertices.ToArray().Length;
                        innerVertices = totalVertices - ghostVertices;
                        borderVertices.Clear();
                    }

                    if (!showJumps) {
                        // ++++++++++++++++++++++++++
                        // average at cell boundaries
                        // ++++++++++++++++++++++++++

                        multiplicity = new byte[totalVertices];

                        bool[] outerBorderVertices = new bool[totalVertices];

                        for (int cell = 0; cell < NoOfCells; cell++) {
                            for (int vertexInCell = 0; vertexInCell < verticesPerCell; vertexInCell++) {
                                int vertex = cellVertices[cell, vertexInCell];

                                if (cell >= NoOfLocalCells && !outerBorderVertices[vertex]) {
                                    borderVertices.Add(vertex);
                                    outerBorderVertices[vertex] = true;
                                }
                                multiplicity[vertex]++;
                            }
                        }

                        ghostVertices = borderVertices.ToArray().Length;
                        innerVertices = totalVertices - ghostVertices;

                        if (!parallelMode || !showGhostCells && parallelMode) {
                            if (!showGhostCells) {
                                totalCells = innerCells;
                            }
                            borderVertices.Clear();
                            return;
                        }

                        verticesAlreadyProcessed = new bool[totalVertices];

                        noOfExternalCells = new int[1];
                        noOfExternalCells[0] = NoOfCells - NoOfLocalCells;

                        int[] result = new int[context.MpiSize];

                        noOfExternalCells = GatherGlobalNumberOfGhostCells(noOfExternalCells, result);

                        maxLocalNoOfExternalCells = -1;

                        for (int i = 0; i < noOfExternalCells.Length; i++) {
                            if (maxLocalNoOfExternalCells < noOfExternalCells[i]) {
                                maxLocalNoOfExternalCells = noOfExternalCells[i];
                            }
                        }

                        long[] globalIndicesExternals = context.Parallel.GlobalIndicesExternalCells;
                        long[] indices = new long[maxLocalNoOfExternalCells * context.MpiSize];

                        for (int i = 0; i < indices.Length; i++) {
                            indices[i] = -1;
                        }

                        long[] mySharedCellIndices = new long[maxLocalNoOfExternalCells];

                        for (int i = NoOfLocalCells; i < NoOfCells; i++) {
                            int cellIndex = i - NoOfLocalCells;

                            mySharedCellIndices[cellIndex] = globalIndicesExternals[cellIndex];
                        }

                        sharedCellIndices = GatherGlobalSharedCellIndices(mySharedCellIndices, indices);

                        globalNoOfExternalCells = sharedCellIndices.Length;

                        byte[,] mySharedVerticesMultiplicity = new byte[maxLocalNoOfExternalCells, verticesPerCell];

                        for (int i = NoOfLocalCells; i < NoOfCells; i++) {
                            for (int j = 0; j < verticesPerCell; j++) {
                                mySharedVerticesMultiplicity[i - NoOfLocalCells, j] = multiplicity[cellVertices[i, j]];
                            }
                        }

                        byte[,] verticesMultiplicity = new byte[globalNoOfExternalCells, verticesPerCell];

                        sharedVerticesMultiplicity = GatherGlobalSharedVerticesMultiplicity(mySharedVerticesMultiplicity, verticesMultiplicity);

                        for (int proc = 0; proc < noOfExternalCells.Length; proc++) {
                            int start = maxLocalNoOfExternalCells * proc;
                            int end = maxLocalNoOfExternalCells * proc + noOfExternalCells[proc];

                            for (int i = 0; i < totalVertices; i++) {
                                verticesAlreadyProcessed[i] = false;
                            }

                            for (int i = start; i < end; i++) {
                                long localCell = sharedCellIndices[i] - context.Grid.CellPartitioning.i0;

                                if (localCell < 0 || localCell > NoOfLocalCells) {
                                    continue;
                                }
                                for (int j = 0; j < verticesPerCell; j++) {
                                    int localVertex = cellVertices[localCell, j];

                                    if (verticesAlreadyProcessed[localVertex]) {
                                        continue;
                                    }

                                    verticesAlreadyProcessed[localVertex] = true;

                                    if (borderVertices.Contains(localVertex)) {
                                        if (multiplicity[localVertex] < sharedVerticesMultiplicity[i, j]) {
                                            multiplicity[localVertex] = sharedVerticesMultiplicity[i, j];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// computes <see cref="ZoneDriver.vertices"/> and
            /// <see cref="ZoneDriver.cellVertices"/>
            /// </summary>
            /// <param name="celVtx">
            /// input:
            /// vertices of each cell, unmerged;
            /// 1st index: cell index;
            /// 2nd index: vertex index within cell;
            /// 3rd index: spatial direction;
            /// </param>
            /// <param name="vertice">
            /// output: the merged vertices, i.e. each vertex from the <paramref name="celVtx"/>
            /// occurs only once.
            /// 1st index: vertex index;
            /// 2nd index: spatial direction
            /// </param>
            /// <param name="cellVertices">
            /// output: indices into the <paramref name="vertice"/>-array.
            /// 1st index: cell index;
            /// 2nd index: vertex index within cell;
            /// </param>
            /// <param name="h">
            /// cell measure
            /// </param>
            public static void InitializeVertice2(MultidimensionalArray celVtx, out int[,] cellVertices, out double[,] vertice, double[] h) {
                using (var tr = new FuncTrace()) {

                    // fk: bemerkung: sollte einegermaßen skalieren;
                    // Test am 10sept10: Gitter mit 400x400 dauert cd 8 sec., Gitter mit 800x400 dauert 16 sekunden

                    int J_update = celVtx.GetLength(0);
                    int D = celVtx.GetLength(2);
                    int J = J_update;


                    int NV = celVtx.GetLength(1);


                    // collect vertices
                    // ================

                    MultidimensionalArray Verts = MultidimensionalArray.Create(J * NV, D);

                    // loop over locally updated cells...
                    int cnt = 0;
                    for (int j = 0; j < J_update; j++) {

                        // loop over cell vertices...
                        for (int iv = 0; iv < NV; iv++) {
                            for (int d = 0; d < D; d++)
                                Verts[cnt, d] = celVtx[j, iv, d];
                            cnt++;
                        }
                    }

                    // build tree
                    // ==========

                    BoundingBox bb = new BoundingBox(Verts);
                    bb.ExtendByFactor(0.005);
                    int[] Perm = new int[Verts.GetLength(0)];
                    PointLocalization locTree = new PointLocalization(Verts, bb, Perm);
                    Verts = locTree.Points;

                    // eliminate duplicate points
                    // ==========================

                    double[] pt = new double[D];
                    int N = J * NV;
                    List<int> foundPoints = new List<int>();
                    int[] AliasPts = new int[J * NV];
                    ArrayTools.SetAll(AliasPts, int.MinValue);
                    int NewVertice = 0;
                    List<int> VerticeTmp = new List<int>();


                    using (new BlockTrace("duplicate Point elimination", tr)) {
                        for (int j = 0; j < J; j++) {
                            for (int _n = 0; _n < NV; _n++) {
                                int n = j * NV + _n;

                                if (AliasPts[n] >= 0)
                                    continue; // point is already assigned.

                                for (int d = 0; d < D; d++)
                                    pt[d] = Verts[n, d];

                                //if (n == 95) 
                                //    Console.Write("");

                                double eps = h[n / NV] * 1.0e-6;

                                locTree.FindNearPoints(foundPoints, eps, pt);
                                if (foundPoints.Count < 1)
                                    throw new ApplicationException("error in algorithm");
                                if (!foundPoints.Contains(n)) {
                                    throw new ApplicationException("error in algorithm");
                                }

                                VerticeTmp.Add(n);

                                for (int k = 0; k < foundPoints.Count; k++) {
                                    AliasPts[foundPoints[k]] = NewVertice;
                                }

                                NewVertice++;
                            }
                        }
                    }

                    // test
                    // ====
                    for (int i = 0; i < AliasPts.Length; i++)
                        if (AliasPts[i] < 0)
                            throw new ApplicationException("error in alg");

                    // store results
                    // =============
                    int[] PermInv = new int[Perm.Length];
                    for (int i = 0; i < PermInv.Length; i++)
                        PermInv[Perm[i]] = i;

                    cellVertices = new int[J, NV];
                    for (int j = 0; j < J; j++) {
                        for (int _n = 0; _n < NV; _n++) {
                            int n = j * NV + _n;
                            cellVertices[j, _n] = AliasPts[PermInv[n]];
                        }
                    }

                    vertice = new double[VerticeTmp.Count, D];
                    for (int i = 0; i < NewVertice; i++) {
                        int n = VerticeTmp[i];
                        for (int d = 0; d < D; d++) {
                            vertice[i, d] = Verts[n, d];
                        }
                    }
                }
            }

            /// <summary>
            /// Evaluate a field at every single point of the (possibly sub-
            /// divided) grid.
            /// </summary>
            /// <param name="Evaluator">The function to be evaluated</param>
            /// <returns>
            /// A list of results. See
            /// <see cref="ScalarFunctionEx"/>
            /// for information about the indexing.
            /// </returns>
            private MultidimensionalArray SampleField(ScalarFunctionEx Evaluator) {
                MultidimensionalArray result = MultidimensionalArray.Create(NoOfCells, verticesPerCell);

                IEnumerable<Chunk> enumerable;
                if (NoOfCells == sgrd.VolumeMask.NoOfItemsLocally) {
                    enumerable = sgrd.VolumeMask;
                } else {
                    enumerable = sgrd.VolumeMask.GetEnumerableWithExternal();
                }


                int iKref = Array.IndexOf(context.Grid.RefElements, this.Zone_Element);

                foreach (var cnk in enumerable) {
                    Evaluator(cnk.i0, cnk.Len, localVerticeCoordinates,
                        result.ExtractSubArrayShallow(new int[] { cnk.i0, 0 }, new int[] { cnk.JE - 1, verticesPerCell - 1 }));
                }
                return result;
            }

            /// <summary>
            /// Evaluation of a function <paramref name="field"/> on the plotting grid
            /// and optional averaging.
            /// </summary>
            /// <remarks>
            /// The output is stored in members <see cref="smoothedResult"/> (if <see cref="showJumps"/>==false)
            /// or <see cref="notSmoothedResult"/> (if <see cref="showJumps"/>==true).
            /// </remarks>
            protected void SampleField(ScalarFunctionEx field, bool showJumps) {
                MultidimensionalArray originalResult = SampleField(field);

                if (showJumps) {
                    notSmoothedResult = originalResult.Storage;
                    return;
                }

                smoothedResult = new double[totalVertices];
                for (int i = 0; i < NoOfCells; i++) {
                    for (int j = 0; j < verticesPerCell; j++) {
                        smoothedResult[cellVertices[i, j]] += originalResult[i, j];
                    }
                }

                // If not in parallel, take a short cut
                if (!parallelMode || (!showJumps && !ghostZone && parallelMode)) {
                    for (int i = 0; i < totalVertices; i++) {
                        smoothedResult[i] /= multiplicity[i];
                    }
                    return;
                }

                int NoOfLocalCells = sgrd.LocalNoOfCells;
                double[,] mySharedVerticesSmoothedResult = new double[maxLocalNoOfExternalCells, verticesPerCell];
                for (int i = NoOfLocalCells; i < NoOfCells; i++) {
                    int cellIndex = i - NoOfLocalCells;

                    for (int j = 0; j < verticesPerCell; j++) {
                        mySharedVerticesSmoothedResult[cellIndex, j] = smoothedResult[cellVertices[i, j]];
                    }
                }

                double[,] result = new double[globalNoOfExternalCells, verticesPerCell];

                double[,] sharedVerticesSmoothedResult = GatherGlobalSharedVerticesSmoothedResult(mySharedVerticesSmoothedResult, result);

                for (int proc = 0; proc < noOfExternalCells.Length; proc++) {
                    int start = maxLocalNoOfExternalCells * proc;
                    int end = maxLocalNoOfExternalCells * proc + noOfExternalCells[proc];

                    for (int i = 0; i < totalVertices; i++) {
                        verticesAlreadyProcessed[i] = false;
                    }

                    for (int i = start; i < end; i++) {
                        long localCell = sharedCellIndices[i] - context.CellPartitioning.i0;

                        if (localCell < 0 || localCell > NoOfLocalCells) {
                            continue;
                        }
                        for (int j = 0; j < verticesPerCell; j++) {
                            int localVertex = cellVertices[localCell, j];

                            if (verticesAlreadyProcessed[localVertex]) {
                                continue;
                            }

                            verticesAlreadyProcessed[localVertex] = true;

                            if (borderVertices.Contains(localVertex)) {
                                double newValue = sharedVerticesSmoothedResult[i, j];
                                double oldValue = smoothedResult[localVertex];

                                if (oldValue * oldValue < newValue * newValue) {
                                    smoothedResult[localVertex] = newValue;
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < totalVertices; i++) {
                    smoothedResult[i] /= multiplicity[i];
                }
            }

            /// <summary>
            /// Gathers the number of _local_ cells that share a given vertex
            /// (on the partition boundary) to all processes. If the multiplicity
            /// is zero, the vertex is not actually located at the boundary but
            /// just happens to belong to a cell that has at least one vertex
            /// located on this boundary
            /// </summary>
            /// <param name="myMultiplicity">
            /// The number of _local_ cells that share a given vertex
            /// </param>
            /// /// <param name="multiplicity">
            /// The number of all global cells that share a given global vertex
            /// </param>
            /// <returns>
            /// [multiplicityProc1, ..., multiplicityProcN] where
            /// |multiplicityProcI| = |sharedIndicesProcI| x |NoOfVerticesPerCell|
            /// (see <see cref="GatherGlobalSharedCellIndices"/> for the definition
            /// of sharedIndicesProcI)
            /// </returns>
            /// <remarks>
            /// The cell can be identified uniquely using the result of
            /// <see cref="GatherGlobalSharedCellIndices"/>
            /// </remarks>
            private byte[,] GatherGlobalSharedVerticesMultiplicity(byte[,] myMultiplicity, byte[,] multiplicity) {
                GCHandle sendPin = GCHandle.Alloc(myMultiplicity, GCHandleType.Pinned);
                GCHandle receivePin = GCHandle.Alloc(multiplicity, GCHandleType.Pinned);

                IntPtr sendBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(myMultiplicity, 0);
                IntPtr receiveBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(multiplicity, 0);
                csMPI.Raw.Allgather(
                    sendBuffer,
                    myMultiplicity.Length,
                    csMPI.Raw._DATATYPE.BYTE,
                    receiveBuffer,
                    myMultiplicity.Length,
                    csMPI.Raw._DATATYPE.BYTE,
                    csMPI.Raw._COMM.WORLD);

                sendPin.Free();
                receivePin.Free();

                return multiplicity;
            }

            /// <summary>
            /// Gathers the sum of the field values of all _local_ cells that share
            /// a given vertex (on the partition boundary) to all processes.
            /// </summary>
            /// <param name="myResult">
            /// The sum of the field values of all _local_ cells that share a given
            /// vertex
            /// </param>
            /// <param name="result">
            /// The sum of the field values over all processes
            /// </param>
            /// <returns>
            /// [sumProc1, ..., sumProcN] where
            /// |sumProcI| = |sharedIndicesProcI| x |NoOfVerticesPerCell|
            /// (see <see cref="GatherGlobalSharedCellIndices"/> for the definition
            /// of sharedIndicesProcI)
            /// </returns>
            /// <remarks>
            /// The cell can be identified uniquely using the result of
            /// <see cref="GatherGlobalSharedCellIndices"/>
            /// </remarks>
            private double[,] GatherGlobalSharedVerticesSmoothedResult(double[,] myResult, double[,] result) {
                GCHandle sendPin = GCHandle.Alloc(myResult, GCHandleType.Pinned);
                GCHandle receivePin = GCHandle.Alloc(result, GCHandleType.Pinned);

                IntPtr sendBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(myResult, 0);
                IntPtr receiveBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(result, 0);
                csMPI.Raw.Allgather(
                    sendBuffer,
                    myResult.Length,
                    csMPI.Raw._DATATYPE.DOUBLE,
                    receiveBuffer,
                    myResult.Length,
                    csMPI.Raw._DATATYPE.DOUBLE,
                    csMPI.Raw._COMM.WORLD);

                sendPin.Free();
                receivePin.Free();

                return result;
            }

            /// <summary>
            /// Gathers global cell indices of shared (external) cells to all
            /// processes. Mainly needed to define an ordering for the values the
            /// shared multiplicities (see
            /// <see cref="GatherGlobalSharedVerticesMultiplicity"/>) and the
            /// shared results (see
            /// <see cref="GatherGlobalSharedVerticesMultiplicity"/>)
            /// </summary>
            /// <param name="myIndices">
            /// The global indices of external cells of this process
            /// </param>
            /// <param name="indices">
            /// The global indices of external cells of all processes
            /// </param>
            /// <returns>
            /// [sharedIndicesProc1, ..., sharedIndicesProcN]
            /// </returns>
            private long[] GatherGlobalSharedCellIndices(long[] myIndices, long[] indices) {

                GCHandle sendPin = GCHandle.Alloc(myIndices, GCHandleType.Pinned);
                GCHandle receivePin = GCHandle.Alloc(indices, GCHandleType.Pinned);

                IntPtr sendBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(myIndices, 0);
                IntPtr receiveBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(indices, 0);
                csMPI.Raw.Allgather(
                    sendBuffer,
                    myIndices.Length,
                    csMPI.Raw._DATATYPE.LONG_LONG_INT,
                    receiveBuffer,
                    myIndices.Length,
                    csMPI.Raw._DATATYPE.LONG_LONG_INT,
                    csMPI.Raw._COMM.WORLD);

                sendPin.Free();
                receivePin.Free();

                return indices;
            }

            /// <summary>
            /// Determine number of ghost cells on each process
            /// </summary>
            /// <param name="myResult">
            /// The number of ghost cells in the current process
            /// </param>
            /// <param name="result">
            /// The number of ghost cells in all processes
            /// </param>
            /// <returns>
            ///  Determines the number of ghost cells on each process and returns an array with process number indexing
            /// </returns>
            private int[] GatherGlobalNumberOfGhostCells(int[] myResult, int[] result) {
                GCHandle sendPin = GCHandle.Alloc(myResult, GCHandleType.Pinned);
                GCHandle receivePin = GCHandle.Alloc(result, GCHandleType.Pinned);

                IntPtr sendBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(myResult, 0);
                IntPtr receiveBuffer = Marshal.UnsafeAddrOfPinnedArrayElement(result, 0);
                csMPI.Raw.Allgather(
                    sendBuffer,
                    1,
                    csMPI.Raw._DATATYPE.INT,
                    receiveBuffer,
                    1,
                    csMPI.Raw._DATATYPE.INT,
                    csMPI.Raw._COMM.WORLD);

                sendPin.Free();
                receivePin.Free();

                return result;
            }

            /// <summary>
            /// Gets the global cell indices of shared (external) cells to all
            /// processes. Mainly needed to define an ordering for the values the
            /// shared multiplicities (see
            /// <see cref="GatherGlobalSharedVerticesMultiplicity"/>) and the
            /// shared results (see
            /// <see cref="GatherGlobalSharedVerticesMultiplicity"/>)
            /// </summary>
            /// <param name="indices">
            /// The global indices of external cells of all processes
            /// </param>
            /// <param name="globalIndexValue">
            /// The given global index value
            /// </param>
            /// <returns>
            /// The index to the given global index value
            /// </returns>
            private int GetGlobalSharedCellIndices(long[] indices, long globalIndexValue) {
                for (int i = 0; i < indices.Length; i++) {
                    if (indices[i] == globalIndexValue) {
                        return i;
                    }
                }
                return 0;
            }

            /// <summary>
            /// Implement this methods by exporting/plotting all functions contained
            /// in <paramref name="fieldsToPlot"/> while regarding the settings
            /// <see cref="superSampling"/> and <see cref="showJumps"/>.
            /// </summary>
            /// <param name="ZoneName">
            /// The name of the zone
            /// </param>
            /// <param name="time">
            /// The solution time represented by <paramref name="fieldsToPlot"/>
            /// </param>
            /// <param name="fieldsToPlot">
            /// A list of mappings between field names and the fields which
            /// should be plotted
            /// </param>
            abstract public void PlotZone(string ZoneName, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> fieldsToPlot);
        }


        /// <summary>
        /// Generates file name
        /// </summary>
        /// <param name="fileNameBase">The base of the file name</param>
        /// <returns>
        /// A file name which does not already exist in the current working
        /// directory
        /// </returns>
        protected virtual string GenerateFileName(string fileNameBase) {
            int size = gridData.MpiSize;

            string MPI_prefix = "";
            if (size > 1) {
                MPI_prefix = "MPI" + (gridData.MpiRank + 1) + "of" + size + ".";
            }

            string filename = MPI_prefix + fileNameBase + "." + FileEnding;
            int cnt = 2;
            while (File.Exists(filename)) {
                filename = fileNameBase + "." + cnt + MPI_prefix + "." + FileEnding;
                cnt++;
            }

            return filename;
        }

        /// <summary>
        /// Implement this methods by exporting/plotting all fields contained
        /// in <paramref name="fieldsToPlot"/> while regarding the settings
        /// <see cref="ZoneDriver.superSampling"/> and
        /// <see cref="ZoneDriver.showJumps"/>.
        /// </summary>
        /// <param name="fileNameBase">
        /// The base name of the file (if any) that should be created. If
        /// called in parallel, rank information will be appended. If a file
        /// with the given name already exists, the resulting file name will
        /// contain an additional counter.
        /// </param>
        /// <param name="time">
        /// The solution time represented by <paramref name="fieldsToPlot"/>
        /// </param>
        /// <param name="fieldsToPlot">
        /// A list of fields which should be plotted
        /// </param>
        virtual public void PlotFields(string fileNameBase, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> fieldsToPlot) {
            this.OpenFile(this.GenerateFileName(fileNameBase), fieldsToPlot.Select(x => x.Item1));

            for (int i = 0; i < this.ZoneDrivers.Length; i++) {
                if (ZoneDrivers[i] != null) {
                    ZoneDrivers[i].PlotZone(
                        "Zone_" + gridData.Grid.GetRefElement(i).GetType().Name,
                        time,
                        fieldsToPlot);
                }
            }

            this.CloseFile();
        }

        /// <summary>
        /// see <see cref="PlotFields(string, string, double, IEnumerable{DGField})"/>.
        /// </summary>
        virtual public void PlotFields(string fileNameBase, double time, IEnumerable<DGField> fieldsToPlot) {
            this.PlotFields(fileNameBase, time, fieldsToPlot.Select(x => new Tuple<string, ScalarFunctionEx>(x.Identification, x.Evaluate)));
        }

        /// <summary>
        /// see <see cref="PlotFields(string, string, double, IEnumerable{DGField})"/>.
        /// </summary>
        virtual public void PlotFields(string fileNameBase, double time, params DGField[] fieldsToPlot) {
            this.PlotFields(fileNameBase, time, (IEnumerable<DGField>)fieldsToPlot);
        }

        /// <summary>
        /// Implement this method by returning the file ending that should be
        /// appended at the end of a file name.
        /// </summary>
        abstract protected string FileEnding {
            get;
        }
    }
}

