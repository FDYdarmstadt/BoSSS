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

        readonly byte[] edgeTags; 

        public PeriodicBoundaryConverter(
            VoronoiBoundary boundary,
            PeriodicMap map)
        {
            edgePairer = new EdgePairer();
            edgeTags = ExtractEdgeTags(boundary, map);
            this.boundaryMap = 
                CreatePeriodicBoundaryMap<SortedList<byte, AffineTrafo>, LinkedListDictionary<int, bool>>(
                    edgeTags, 
                    map);
        }

        public void Clear()
        {
            edgePairer.Clear();
        }

        static byte[] ExtractEdgeTags(VoronoiBoundary boundary, PeriodicMap map)
        {
            byte[] voronoiEdgeTags = boundary.EdgeTags;
            byte[] periodicCornerEdgeTags = GetPeriodicCornerEdgeTags(
                map.PeriodicCornerCorrelation,
                map.PeriodicBoundaryCorrelation);
            byte[] edgeTags = Concat(voronoiEdgeTags, periodicCornerEdgeTags);
            return edgeTags;
        }

        static T[] Concat<T>(T[] a, T[] b)
        {
            T[] both = new T[a.Length + b.Length];
            a.CopyTo(both, 0);
            b.CopyTo(both, a.Length);
            return both;
        }

        static byte[] GetPeriodicCornerEdgeTags(
            IDictionary<Corner, int> periodicCornerMap,
            IDictionary<int, int> periodicBoundaryCorrelation)
        {
            int periodicCornerCount = periodicCornerMap.Count;
            int allBoundariesCount = periodicBoundaryCorrelation.Count;
            int nonCornerEdgeTagCount = allBoundariesCount - periodicCornerCount;

            byte[] newEdgeTags = new byte[periodicCornerCount];
            int i = 0;
            foreach (var corner in periodicCornerMap)
            {
                int cornerIndice = corner.Value - nonCornerEdgeTagCount;
                if (cornerIndice < (periodicCornerCount / 2))
                {
                    byte edgeTag = (byte)(GridCommons.FIRST_PERIODIC_BC_TAG + nonCornerEdgeTagCount / 2 + i);
                    newEdgeTags[cornerIndice] = edgeTag;
                    periodicBoundaryCorrelation.TryGetValue(corner.Value, out int pairedEdge);
                    newEdgeTags[pairedEdge - nonCornerEdgeTagCount] = edgeTag;
                    ++i;
                }
                else
                {
                    break;
                }
            }
            return newEdgeTags;
        }

        static PeriodicBoundaryConverterMap CreatePeriodicBoundaryMap<TtrafoDictionary, TinverseDictionary>(
            byte[] edgeTags,
            PeriodicMap map)
            where TtrafoDictionary : IDictionary<byte, AffineTrafo>, new()
            where TinverseDictionary : IDictionary<int, bool>, new()
        {
            var periodicInverseMap = ExtractPeriodicInverseMap<TinverseDictionary>(map.PeriodicBoundaryCorrelation);

            
            var periodicTrafos = ExtractPeriodicTrafos<TtrafoDictionary>(
                edgeTags,
                periodicInverseMap,
                map.PeriodicBoundaryTransformations);

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
            where T : IDictionary<byte, AffineTrafo>, new() //
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

        public void RegisterPeriodicBoundariesTo(GridCommons grid) {
            Debug.Assert(boundaryMap.PeriodicTrafos is SortedList<byte, AffineTrafo>);
            var periodicTrafos = boundaryMap.PeriodicTrafos;
            //Sort by edgeTag
            foreach (var pair in periodicTrafos) {
                AffineTrafo trafo = pair.Value;
                byte edgeTag = pair.Key;
                byte edgeTagInGrid = grid.AddPeriodicEdgeTrafo(trafo);
                Debug.Assert(edgeTag == edgeTagInGrid);

                if (!grid.EdgeTagNames.ContainsKey(edgeTagInGrid)) {
                    grid.EdgeTagNames.Add(edgeTagInGrid, $"periodic corner {edgeTagInGrid}");
                }
            }
        }

        public void SetPeriodicData(
            int tagIndice,
            CellFaceTag[] tags,
            int boundaryNumber,
            long globalCellID,
            int edgeID,
            int edgeTwinID)
        {
            tags[tagIndice].ConformalNeighborship = true;
            tags[tagIndice].PeriodicInverse = boundaryMap.IsPeriodicInverse(boundaryNumber);
            if (tags[tagIndice].PeriodicInverse)
            {
                edgePairer.SetNeighborCell(edgeID, tagIndice, tags, globalCellID);
            }
            else
            {
                edgePairer.SetNeighborCell(edgeTwinID, tagIndice, tags, globalCellID);
            }
        }

        static AffineTrafo Convert(Transformation transformation)
        {
            var trafo = new AffineTrafo(transformation.Matrix, transformation.AffineTransformation);
            return trafo;
        }

        public byte GetEdgeTagOf(int boundaryEdgeNumber)
        {
            return edgeTags[boundaryEdgeNumber];
        }
    }
}
