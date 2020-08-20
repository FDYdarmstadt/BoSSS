using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System.Linq;
using MPI.Wrappers.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP;
using System.Diagnostics;
using BoSSS.Solution;

namespace BoSSS.Solution.Tecplot
{
    /// <summary>
    /// Tecplot plotting. Aggregated Cells are depicted as lines on a dedicated mesh.
    /// </summary>
    public class AggregatedTecplot : Tecplot
    {
        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public AggregatedTecplot(AggregationGridData gridData, uint superSampling) : base(gridData, superSampling)
        {
        }

        /// <summary>
        /// creates a zone driver for Tecplot. Aggregated Cells are depicted as lines on a dedicated mesh.
        /// </summary>
        protected override PlotDriver.ZoneDriver CreateZoneDriver(IGridData context, int iKref, bool showJumps, bool showGhostCells, uint superSampling, CellMask __mask)
        {
            return new AggregationZone((AggregationGridData)context, iKref, showJumps, showGhostCells, superSampling, __mask);
        }

        class AggregationZone : TecplotZone
        {
            readonly AggregationGridData gridData;

            UnsafeTECIO m_TECPLOT = new UnsafeTECIO();

            public AggregationZone(AggregationGridData gridData, int iKref, bool showJumps, bool ghostZone, uint superSampling, CellMask mask = null) 
                : base(gridData, iKref, showJumps, ghostZone, superSampling, mask)
            {
                this.gridData = gridData;
            }

            /// <summary>
            /// Plots all fields in <paramref name="fieldsToPlot"/> into a Tecplot file, see
            /// <see cref="PlotDriver.ZoneDriver.PlotZone"/>.
            /// Aggregated Cells are depicted as a line Zone with a dedicated mesh.
            /// </summary>
            public override void PlotZone(string ZoneName, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> fieldsToPlot)
            {
                base.PlotZone(ZoneName, time, fieldsToPlot);
                PlotAsFELineSegZone(time, fieldsToPlot.Count(), gridData, m_TECPLOT);
            }

            /// <summary>
            /// The beloved spaghetti style and flavor
            /// </summary>
            /// <param name="time"></param>
            /// <param name="fieldCount"></param>
            /// <param name="gridData"></param>
            unsafe static void PlotAsFELineSegZone(double time, int fieldCount, AggregationGridData gridData, UnsafeTECIO m_TECPLOT)
            {
                int[,] edge2LogicalCell = gridData.iGeomEdges.LogicalCellIndices;
                int[][] cellVerts = gridData.iGeomCells.CellVertices;
                int[,] edge2GeomCell = gridData.iGeomEdges.CellIndices;
                byte[,] edge2Face = gridData.iGeomEdges.FaceIndices;
                byte[] edgeTags = gridData.iGeomEdges.EdgeTags;

                //Find Line Connections: Each Periodic BoundaryEdge Twice
                List<(int Vertex0, int Vertex1)> connection = new List<(int, int)>(gridData.iLogicalEdges.Count);
                for (int geomEdgeIndex = 0; geomEdgeIndex < gridData.iGeomEdges.Count; ++geomEdgeIndex)
                {
                    int logicalCell1 = edge2LogicalCell[geomEdgeIndex, 0];
                    int logicalCell2 = edge2LogicalCell[geomEdgeIndex, 1];
                    if (logicalCell1 != logicalCell2)
                    {
                        if (edge2GeomCell[geomEdgeIndex, 0] > -1)
                        {
                            AddEdgeToConnection(geomEdgeIndex, 0);
                            if (edgeTags[geomEdgeIndex] >= GridCommons.FIRST_PERIODIC_BC_TAG) //Periodic Edges
                            {
                                AddEdgeToConnection(geomEdgeIndex, 1);
                            }
                        }
                        else
                        {
                            AddEdgeToConnection(geomEdgeIndex, 1);
                        }
                    }
                }

                void AddEdgeToConnection(int edgeIndex, int sideOfEdge)
                {
                    int geomCell = edge2GeomCell[edgeIndex, sideOfEdge];
                    int faceIndice = edge2Face[edgeIndex, sideOfEdge];

                    int numberOfFaces = gridData.iGeomCells.GetRefElement(0).NoOfFaces;
                    int[] geomVertices = cellVerts[geomCell];
                    int vertex0 = geomVertices[faceIndice] + 1;
                    int vertex1 = geomVertices[(faceIndice + 1) % numberOfFaces] + 1;
                    connection.Add((vertex0, vertex1));
                }

                //Setup Zone
                MultidimensionalArray vertices = ((GridData)gridData.ParentGrid).Vertices.Coordinates;
                int numberOfVertices = vertices.GetLength(0);
                int numberOfElements = connection.Count;
                IntPtr ptrZoneTitle = Marshal.StringToHGlobalAnsi("Grid");
                int zoneTypeIndex = 1; //FELineSeg
                int KMax = 0, ICellMax = 0, JCellMax = 0, KCellMax = 0, NFConn = 0, FNMode = 0;
                int IsBlock = 1;
                int StrandID = 0, ParentZone = 0, ShrConn = 0;
                m_TECPLOT.TECZNE110(ptrZoneTitle,
                    ref zoneTypeIndex,
                    ref numberOfVertices,
                    ref numberOfElements,
                    ref KMax,
                    ref ICellMax,
                    ref JCellMax,
                    ref KCellMax,
                    ref time,
                    ref StrandID,
                    ref ParentZone,
                    ref IsBlock,
                    ref NFConn,
                    ref FNMode,
                    null,           // No passive variables
                    null,           // All variables node-centered
                    null,           // No shared variables
                    ref ShrConn);
                Marshal.FreeHGlobal(ptrZoneTitle);

                //Set Vertices
                int VIsDouble = 1;
                double[] globalCoordinates = new double[numberOfVertices];
                for (int d = 0; d < gridData.SpatialDimension; d++)
                {
                    for (int vertex = 0; vertex < numberOfVertices; vertex++)
                    {
                        globalCoordinates[vertex] = vertices[vertex, d];
                    }

                    int n = globalCoordinates.Length;
                    m_TECPLOT.TECDAT110(ref n, globalCoordinates, ref VIsDouble);
                }

                //Set Values
                double[] values = new double[numberOfVertices];
                for (int i = 0; i < fieldCount; ++i)
                {
                    m_TECPLOT.TECDAT110(ref numberOfVertices, values, ref VIsDouble);
                }

                //Set Connection
                int[,] connectionArray = new int[numberOfElements, 2];
                for (int i = 0; i < numberOfElements; ++i)
                {
                    connectionArray[i, 0] = connection[i].Vertex0;
                    connectionArray[i, 1] = connection[i].Vertex1;
                }
                fixed (int* ptr = connectionArray)
                {
                    m_TECPLOT.TECNOD110(ptr);
                }
            }
        }
    }
}
