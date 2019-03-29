using System;
using System.Collections.Generic;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Connectors.Matlab;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi
{
    static class MatlabMeshGenerator
    {
        class MatlabRidgeComparer : EqualityComparer<MatlabRidge>
        {
            public override bool Equals(MatlabRidge x, MatlabRidge y)
            {
                return ((x.Start == y.Start && x.End == y.End) || (x.Start == y.End && x.End == y.Start));
            }

            public override int GetHashCode(MatlabRidge obj)
            {
                //http://eternallyconfuzzled.com/tuts/algorithms/jsw_tut_hashing.aspx 
                //If x == y hash must be hash(x) = hash(y)
                int start = obj.End > obj.Start ? obj.Start : obj.End; 
                int end = obj.End < obj.Start ? obj.Start : obj.End;

                int hash = 17;
                hash = hash * 23 + start.GetHashCode();
                hash = hash * 23 + end.GetHashCode();

                return hash;
            }
        }

        class MatlabRidge
        {
            public int Start;
            public int End;

            public MatlabRidge(int i_start, int i_end)
            {
                Start = i_start;
                End = i_end;
            }
        }

        class OrderedMatlabMesh : SimpleIdMesh
        {
            public OrderedMatlabMesh(int[][] VocellVertexIndex, List<Vertex> vertexCoordinates)
            {
                cells = CreateCells(VocellVertexIndex, vertexCoordinates);
                vertices = vertexCoordinates;
            }

            static List<Cell> CreateCells(int[][] vertIndex, List<Vertex> vertCoords)
            {
                Dictionary<MatlabRidge, Edge> ridgeMap = new Dictionary<MatlabRidge, Edge>(vertCoords.Count,
                                                                                             new MatlabRidgeComparer());
                List<Cell> cells = new List<Cell>(vertIndex.Length);
                for (int i = 0; i < vertIndex.Length; ++i)
                {
                    int[] cell = vertIndex[i];
                    Cell voronoiCell = new Cell { ID = i };
                    Edge[] Ridges = new Edge[cell.Length];
                    Vertex[] Vertices = new Vertex[cell.Length];


                    //Collect Vertices
                    for (int j = 0; j < cell.Length; ++j)
                    {
                        Vertices[j] = vertCoords[cell[j]];
                    }
                    //Correct Orientation of Cells
                    bool infCell = false;
                    for (int j = 0; j < cell.Length; ++j)
                    {
                        if (cell[j] == 0)
                        {
                            infCell = true;
                            break;
                        }
                    }
                    if (!infCell)
                    {
                        int sign = CheckOrientation(Vertices);
                        if (sign > 0)
                        {
                            // nop
                        }
                        else if (sign < 0)
                        {
                            Array.Reverse(Vertices);
                            Array.Reverse(cell);
                        }
                        else
                        {
                            throw new ArithmeticException("got indefinite polygon from matlab");
                        }
                    }

                    //Collect Ridges
                    for (int j = 0; j < cell.Length; ++j)
                    {
                        int j_plus_1 = (j + 1) % cell.Length;
                        MatlabRidge matlabRidge = new MatlabRidge(cell[j], cell[j_plus_1]);
                        Edge voronoiRidge = new Edge
                        {
                            Cell = voronoiCell,
                            Start = vertCoords[cell[j]],
                            End = vertCoords[cell[j_plus_1]]
                        };
                        Edge twinRidge;
                        if (ridgeMap.TryGetValue(matlabRidge, out twinRidge))
                        {
                            voronoiRidge.Twin = twinRidge;
                            twinRidge.Twin = voronoiRidge;
                        }
                        else
                        {
                            ridgeMap.Add(matlabRidge, voronoiRidge);
                        }
                        Ridges[j] = voronoiRidge;

                    }

                    //Assemble Cell
                    voronoiCell.Vertices = Vertices;
                    voronoiCell.Edges = Ridges;
                    cells.Add(voronoiCell);
                }
                return cells;
            }

            static int CheckOrientation(Vertex[] vertices)
            {
                int L = vertices.Length;

                double[] signs = new double[L - 2];

                bool AllPos = true;
                bool AllNeg = true;

                for (int iTri = 0; iTri < L - 2; iTri++)
                { // loop over triangles of voronoi cell
                    int iV0 = 0;
                    int iV1 = iTri + 1;
                    int iV2 = iTri + 2;

                    Vector V0 = vertices[iV0].Position;
                    Vector V1 = vertices[iV1].Position;
                    Vector V2 = vertices[iV2].Position;

                    Vector D1 = V1 - V0;
                    Vector D2 = V2 - V0;

                    signs[iTri] = D1.CrossProduct2D(D2);

                    AllPos = AllPos && (signs[iTri] > 0);
                    AllNeg = AllNeg && (signs[iTri] < 0);
                }

                if (AllNeg == AllPos)
                    return 0;
                //throw new ArithmeticException("Indefinite polygon");

                if (AllPos)
                    return 1;

                if (AllNeg)
                {
                    //Polygon = Polygon.Reverse().ToArray();
                    //iVtx = iVtx.Reverse().ToArray();
                    return -1;
                }

                throw new ArithmeticException("Indefinite polygon");
            }
        }

        public static IntersectionMesh CreateMesh(MultidimensionalArray nodes)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            IntersectionMesh mesh;
            using (var Matlab = new BatchmodeConnector())
            {
                Matlab.PutMatrix(nodes, "Nodes");

                // compute Voronoi diagramm
                Matlab.Cmd("[V, C] = voronoin(Nodes);");

                // output (export from matlab)
                int[][] VocellVertexIndex = new int[nodes.NoOfRows][];
                Matlab.GetStaggeredIntArray(VocellVertexIndex, "C");
                Matlab.GetMatrix(null, "V");

                // run matlab
                Matlab.Execute(false);

                // Vertices: List<Vertex>
                //---------------------------
                MultidimensionalArray matlabVertexCoordinates = (MultidimensionalArray)(Matlab.OutputObjects["V"]);
                List<Vertex> vertexCoordinates = new List<Vertex>(matlabVertexCoordinates.NoOfRows);
                for (int i = 0; i < matlabVertexCoordinates.NoOfRows; ++i)
                {
                    Vertex tmp = new Vertex();
                    tmp.Position = new Vector(matlabVertexCoordinates.GetRow(i));
                    tmp.ID = i;
                    vertexCoordinates.Add(tmp);
                }

                //Cells: int[][]
                //--------------------------------------------
                // correct indices (1-based index to 0-based index)
                
                for (int i = 0; i < VocellVertexIndex.Length; ++i)
                {
                    for (int k = 0; k < VocellVertexIndex[i].Length; k++)
                    {
                        VocellVertexIndex[i][k]--;
                    }
                }

                stopwatch.Stop();
                Console.WriteLine(stopwatch.ElapsedMilliseconds);
                //Convert to mesh data types
                //--------------------------------------------
                mesh = new IntersectionMesh(new OrderedMatlabMesh( VocellVertexIndex, vertexCoordinates));
            }

            return mesh;
        }

        static void ChooseBatchMode()
        {
            if (System.Environment.MachineName.ToLowerInvariant().EndsWith("rennmaschin")
          //|| System.Environment.MachineName.ToLowerInvariant().Contains("jenkins")
          )
            {
                // This is Florians Laptop;
                // he is to poor to afford MATLAB, so he uses OCTAVE
                BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;
                //BatchmodeConnector.MatlabExecuteable = "C:\\cygwin64\\bin\\bash.exe";
            }
        }
 
    }
}
