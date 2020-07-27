using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.MICHMesher
{
    class VariableCell<T> : MeshCell<T>
    {
        public VariableCell()
        {
            edgeList = new List<Edge<T>>();
        }

        public bool boundary = false;

        List<Edge<T>> edgeList;

        public void Insert(Edge<T> a)
        {
            edgeList.Add(a);
        }

        public void Init()
        {
            if (!boundary)
            {
                Edges = GetEdges();
                Vertices = GetVertices();
            }
            else
            {
                //Empty Ridges
                Edges = edgeList.ToArray();
                for (int i = 0; i < Edges.Length; ++i)
                {
                    Edges[i].IsBoundary = true;
                }
            }

        }

        //This should be possible in nlogn time!!111
        Edge<T>[] GetEdges()
        {
            for (int i = 0; i < edgeList.Count() - 2; ++i)
            {
                Edge<T> c_i = edgeList[i];
                for (int j = i + 1; j < edgeList.Count(); ++j)
                {
                    Edge<T> c_j = edgeList[j];
                    if (c_i.End.ID == c_j.Start.ID)
                    {
                        //Ridge is in order
                        if (j != i + 1)
                        {
                            Edge<T> tmp = edgeList[i + 1];
                            edgeList[i + 1] = c_j;
                            edgeList[j] = tmp;
                        }
                        break;
                    }
                }
            }
            return edgeList.ToArray();
        }

        Vertex[] GetVertices()
        {
            Vertex[] verts = new Vertex[Edges.Length];
            for (int i = 0; i < Edges.Length; ++i)
            {
                verts[i] = Edges[i].Start;
            }
            return verts;
        }

    }

}
