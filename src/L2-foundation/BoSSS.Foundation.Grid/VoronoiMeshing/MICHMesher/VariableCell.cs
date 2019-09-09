using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class VariableCell<T> : MeshCell<T>
    {
        public VariableCell()
        {
            ridgeList = new List<Edge<T>>();
        }

        public bool boundary = false;

        List<Edge<T>> ridgeList;

        public void Insert(Edge<T> a)
        {
            ridgeList.Add(a);
        }

        public void Init()
        {
            if (!boundary)
            {
                Edges = GetRidges();
                Vertices = GetVertices();
            }
            else
            {
                //Empty Ridges
                Edges = ridgeList.ToArray();
                for (int i = 0; i < Edges.Length; ++i)
                {
                    Edges[i].IsBoundary = true;
                }
            }

        }

        //This should be possible in nlogn time!!111
        Edge<T>[] GetRidges()
        {
            //Remove zero Ridges
            for (int i = 0; i < ridgeList.Count(); ++i)
            {
                Edge<T> c_i = ridgeList[i];
                if (c_i.End.ID == c_i.Start.ID)
                {
                    ridgeList.RemoveAt(i);
                }
            }
            for (int i = 0; i < ridgeList.Count() - 2; ++i)
            {
                Edge<T> c_i = ridgeList[i];
                for (int j = i + 1; j < ridgeList.Count(); ++j)
                {
                    Edge<T> c_j = ridgeList[j];
                    if (c_i.End.ID == c_j.Start.ID)
                    {
                        //Ridge is in order
                        if (j != i + 1)
                        {
                            Edge<T> tmp = ridgeList[i + 1];
                            ridgeList[i + 1] = c_j;
                            ridgeList[j] = tmp;
                        }
                        break;
                    }
                }
            }
            return ridgeList.ToArray();
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
