using System;
using System.Collections.Generic;
using ilPSP;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.MICHMesher
{
    static class MIConvexHullMeshGenerator
    {
        public static IDMesh<T> CreateMesh<T>(IList<T> nodes)
            where T : ILocatable, new()
        {
            ResetDataIDCounters<T>();
            MICHVertex<T>[] startNodes = Wrap(nodes);
            IDMesh<T> mesh = MICMesher<T>.Create(startNodes);
            return mesh;
        }

        static MICHVertex<T>[] Wrap<T>(IList<T> nodes)
            where T : ILocatable
        {
            MICHVertex<T>[] startNodes = new MICHVertex<T>[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                startNodes[i] = new MICHVertex<T>(nodes[i].Position, nodes[i]);
            }
            return startNodes;
        }

        static void ResetDataIDCounters<T>()
        {
            MICHVertex<T>.Clear();
            MICHDelaunayCell<T>.Clear();
        }
    }
}
