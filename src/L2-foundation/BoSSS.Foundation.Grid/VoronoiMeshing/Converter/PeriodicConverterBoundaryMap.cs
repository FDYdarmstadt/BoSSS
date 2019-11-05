using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    class PeriodicBoundaryConverterMap
    {
        public IDictionary<byte, AffineTrafo> PeriodicTrafos;

        public IDictionary<int, bool> PeriodicInverseMap;

        public bool IsPeriodicInverse(int boundaryEdgeNumber)
        {
            if (PeriodicInverseMap.TryGetValue(boundaryEdgeNumber, out bool isInverse))
            {
                return isInverse;
            }
            else
            {
                throw new Exception("BoundaryEdgeNumber is not registerd as periodic boundary.");
            }
        }

        public AffineTrafo GetTransformation(byte edgeTag)
        {
            if (PeriodicTrafos.TryGetValue(edgeTag, out AffineTrafo trafo))
            {
                return trafo;
            }
            else
            {
                throw new Exception("BoundaryEdgeNumber is not registerd as periodic boundary.");
            }
        }
    }
}
