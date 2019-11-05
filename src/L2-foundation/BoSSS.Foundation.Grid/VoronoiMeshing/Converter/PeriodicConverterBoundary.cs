using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    class PeriodicBoundaryConverter
    {
        readonly EdgePairer edgePairer;

        readonly PeriodicBoundaryConverterMap boundaryMap;

        public PeriodicBoundaryConverter(
            VoronoiBoundary boundary,
            PeriodicMap map)
        {
            edgePairer = new EdgePairer();
            this.boundaryMap = CreatePeriodicBoundaryMap
                <SortedList<byte, AffineTrafo>, LinkedListDictionary<int, bool>>(
                    boundary.EdgeTags, 
                    map.PeriodicBoundaryMap, 
                    map.PeriodicTransformationMap);
        }

        static PeriodicBoundaryConverterMap CreatePeriodicBoundaryMap<TtrafoDictionary, TinverseDictionary>(
            byte[] edgeTags,
            IDictionary<int, int> periodicBoundaryMap,
            IDictionary<int, Transformation> periodicTransformations)
            where TtrafoDictionary : IDictionary<byte, AffineTrafo>, new()
            where TinverseDictionary : IDictionary<int, bool>, new()
        {
            var periodicInverseMap = ExtractPeriodicInverseMap<TinverseDictionary>(
                    periodicBoundaryMap);
            var periodicTrafos = ExtractPeriodicTrafos<TtrafoDictionary>(
                    edgeTags,
                    periodicInverseMap,
                    periodicTransformations);
            var boundaryMap = new PeriodicBoundaryConverterMap()
            {
                PeriodicInverseMap = periodicInverseMap,
                PeriodicTrafos = periodicTrafos
            };
            return boundaryMap;
        }

        static T ExtractPeriodicInverseMap<T>(IDictionary<int, int> periodicBoundaryMap)
            where T : IDictionary<int, bool>, new()
        {
            var periodicInverseMap = new T();
            foreach (var entry in periodicBoundaryMap)
            {
                int edgeNumber = entry.Key;
                int pairedEdgeNumber = entry.Value;
                if (periodicInverseMap.ContainsKey(pairedEdgeNumber))
                {
                    periodicInverseMap.Add(edgeNumber, false);
                }
                else
                {
                    periodicInverseMap.Add(edgeNumber, true);
                }
            }
            return periodicInverseMap;
        }

        static T ExtractPeriodicTrafos<T>(
            byte[] edgeTags,
            IDictionary<int, bool> periodicInverseMap,
            IDictionary<int, Transformation> periodicTransformations)
            where T : IDictionary<byte, AffineTrafo>, new()
        {
            var periodicTrafos = new T();
            foreach (var indiceInversePair in periodicInverseMap)
            {
                int boundaryEdgeNumber = indiceInversePair.Key;
                bool isPeriodicInverse = indiceInversePair.Value;

                if (!isPeriodicInverse)
                {
                    byte edgeTag = edgeTags[boundaryEdgeNumber];
                    AffineTrafo trafo = Convert(periodicTransformations[boundaryEdgeNumber]);
                    periodicTrafos.Add(edgeTag, trafo);
                }
            }
            return periodicTrafos;
        }

        public void RegisterPeriodicBoundariesTo(GridCommons grid)
        {
            Debug.Assert(boundaryMap.PeriodicTrafos is SortedList<byte, AffineTrafo>);
            var periodicTrafos = (SortedList<byte, AffineTrafo>)boundaryMap.PeriodicTrafos;
            //Sort by edgeTag
            foreach (var pair in periodicTrafos)
            {
                AffineTrafo trafo = pair.Value;
                byte edgeTag = pair.Key;
                byte edgeTagInGrid = grid.AddPeriodicEdgeTrafo(trafo);
                Debug.Assert(edgeTag == edgeTagInGrid);
            }
        }

        public SortedList<byte, AffineTrafo> GetPeriodicTrafos()
        {
            Debug.Assert(boundaryMap.PeriodicTrafos is SortedList<byte, AffineTrafo>);
            var periodicTrafos = (SortedList<byte, AffineTrafo>)boundaryMap.PeriodicTrafos;
            return periodicTrafos;
        }

        public void SetPeriodicData(
            int tagIndice,
            CellFaceTag[] tags,
            int boundaryNumber,
            long globalCellID,
            int edgeStartID,
            int edgeEndID)
        {
            tags[tagIndice].ConformalNeighborship = true;
            tags[tagIndice].PeriodicInverse = boundaryMap.IsPeriodicInverse(boundaryNumber);
            if (tags[tagIndice].PeriodicInverse)
            {
                edgePairer.SetNeighborCell(edgeStartID, tagIndice, tags, globalCellID);
            }
            else
            {
                edgePairer.SetNeighborCell(edgeEndID, tagIndice, tags, globalCellID);
            }
        }

        static AffineTrafo Convert(Transformation transformation)
        {
            var trafo = new AffineTrafo(transformation.Matrix, transformation.AffineTransformation);
            return trafo;
        }
    }
}
