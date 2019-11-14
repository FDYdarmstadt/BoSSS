using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.MICHMesher
{
    class MICHVertex<T> : MIConvexHull.IVertex
    {
        public T Node;

        static int ID_Counter = 0;
        public static void Clear()
        {
            ID_Counter = 0;
        }
        public readonly int ID;

        public MICHVertex(Vector position, T node)
        {
            Position = position;
            Node = node;
            ID = ID_Counter;
            ++ID_Counter;
        }

        public double[] Position { get; set; }
    }
}
