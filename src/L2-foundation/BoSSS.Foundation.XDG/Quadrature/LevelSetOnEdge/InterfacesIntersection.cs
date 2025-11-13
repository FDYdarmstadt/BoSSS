using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature.LevelSetOnEdge {



    /// <summary>
    /// 
    /// </summary>
    internal class InterfacesIntersection : IQuadRuleFactory<QuadRule> {


        /// <summary>
        /// Returns all cells which are cut by two level-sets and one of them coincides with some face of the respective cell
        /// </summary>
        static public CellMask ComputeCellMask(LevelSetTracker.LevelSetRegions Regions, int iLevelSet, int jLevelSet, int iKref, CellMask baseMask) {
            var CoincidFaces = Regions.LevSetCoincidingFaces;
            var gdat = Regions.GridDat;
            if(baseMask.MaskType != MaskType.Geometrical) {
                throw new ArgumentException("expecting a geometrical mask", nameof(baseMask));
            } 

            if(CoincidFaces == null) {
                return CellMask.GetEmptyMask(gdat, MaskType.Geometrical);
            }

            int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
            BitArray bitMask = new BitArray(J);

            foreach(int j in baseMask.ItemEnum) {
                if(gdat.iGeomCells.GetRefElementIndex(j) != iKref)
                    continue;

                if(CoincidFaces[j] == null)
                    continue;

                foreach(var t in CoincidFaces[j]) {
                    if(t.iLevSet == iLevelSet || t.iLevSet == jLevelSet)
                        bitMask[j] = true;
                }

            }


            return new CellMask(gdat, bitMask, MaskType.Geometrical);
        }



        public InterfacesIntersection(RefElement Kref, LevelSetTracker.LevelSetRegions __Regions, int __LevelSet1_Index, int __LevelSet2_Index,
            IQuadRuleFactory<CellBoundaryQuadRule> __bndyRule_LevelSet1, IQuadRuleFactory<CellBoundaryQuadRule> __bndyRule_LevelSet2,
            MultidimensionalArray cellVol4Species
            ) { 
            RefElement = Kref;
            this.Regions = __Regions;
            if(!gdat.iGeomCells.RefElements.Contains(Kref))
                throw new ArgumentException("expecting a cell reference element", nameof(Kref));
            LevelSet1_Index = __LevelSet1_Index;
            LevelSet2_Index = __LevelSet2_Index;
            BndyRule_LevelSet1 = __bndyRule_LevelSet1;
            BndyRule_LevelSet2 = __bndyRule_LevelSet2;
            CellVol4Species = cellVol4Species;
        }

        IGridData gdat => Regions.GridDat;

        public RefElement RefElement {
            get;
            private set;
        }

        LevelSetTracker.LevelSetRegions Regions;

        readonly int LevelSet1_Index;
        readonly int LevelSet2_Index;

        readonly IQuadRuleFactory<CellBoundaryQuadRule> BndyRule_LevelSet1;
        readonly IQuadRuleFactory<CellBoundaryQuadRule> BndyRule_LevelSet2;

        readonly MultidimensionalArray CellVol4Species;


        public int[] GetCachedRuleOrders() {
            return [];
        }

        /// <summary>
        /// Extracts the rule from the <paramref name="BndyRule"/>;
        /// If the respective <paramref name="iFace"/> has an empty quadrature, returns null;
        /// </summary>
        QuadRule GetIntersectionRuleImpl(int jCell, IQuadRuleFactory<CellBoundaryQuadRule> BndyRule, int iOrder, int iFace) {
            if(!UseCell(jCell, iFace))
                return null;

            var crp = BndyRule.GetQuadRuleSet(new CellMask(gdat, Chunk.GetSingleElementChunk(jCell), MaskType.Geometrical), iOrder);
            var boundaryRule = crp.First().Rule;

            int K = boundaryRule.NumbersOfNodesPerFace[iFace];
            if(K == 0)
                return null;

            int n0 = 0;
            for(int iF = 0; iF < iFace; iF++)
                n0 += boundaryRule.NumbersOfNodesPerFace[iF];
            
            int D = this.RefElement.SpatialDimension;
            var rule = QuadRule.CreateBlank(this.RefElement, K, D, false);
            rule.Weights.Set(boundaryRule.Weights.ExtractSubArrayShallow([n0], [n0 + K - 1]));
            if(rule.Weights.Sum() <= 0.0)
                return null;
            rule.OrderOfPrecision = boundaryRule.OrderOfPrecision;
            rule.Nodes.Set(boundaryRule.Nodes.ExtractSubArrayShallow([n0, 0], [n0 + K - 1, D - 1]));
            rule.Nodes.LockForever();

            return rule;
        }


        



        IEnumerable<int> GetLogicalNeighbors(int jGeom, int iFace) {
            int jLog = this.gdat.GetLogicalCellIndex(jGeom);
            byte[,] edge2face = this.gdat.iGeomEdges.FaceIndices;
            int[,] edge2cells = this.gdat.iGeomEdges.CellIndices;

            var GeomEdges = this.gdat.GetGeometricalEdgesForCells(jGeom); // get edges for cell `jGeom`

            int InOrOut(int iE) {
                if(edge2cells[iE, 0] == jGeom)
                    return 0;
                if(edge2cells[iE, 1] == jGeom)
                    return 1;
                throw new ArgumentException();
            }


            var GeomNeighs = new HashSet<int>();
            foreach(int iGeomEdge in GeomEdges) { // loop over edges...
                int inOut = InOrOut(iGeomEdge);
                int outIn = inOut == 0 ? 1 : 0;
                if(edge2face[iGeomEdge, inOut] == iFace) { // face matches
                    GeomNeighs.Add(edge2cells[iGeomEdge,outIn]);
                }
            }
            return GeomNeighs;
        }



        bool UseCell(int jCell, int iFace) {

            double NeighVol = 0;            
            foreach(int jNeigh in gdat.GetLogicalCellIndices(GetLogicalNeighbors(jCell, iFace))) {
                NeighVol += CellVol4Species[jNeigh];
            }
            double ThisVol = CellVol4Species[jCell];

            return ThisVol >= NeighVol;
        }




        /// <summary>
        /// driver routine
        /// </summary>
        QuadRule GetIntersectionRuleAsym(int jCell, int iLc_Coincide, int iLs2, int order) {
            if(Regions.LevSetCoincidingFaces == null)
                return null;
            if(Regions.LevSetCoincidingFaces[jCell] == null)
                return null;

            IQuadRuleFactory<CellBoundaryQuadRule> BndyRule;
            if(iLs2 == LevelSet1_Index)
                BndyRule = BndyRule_LevelSet1;
            else if(iLs2 == LevelSet2_Index)
                BndyRule = BndyRule_LevelSet2;
            else
                throw new ArgumentException();

            foreach(var t in Regions.LevSetCoincidingFaces[jCell]) {
                if(t.iLevSet == iLc_Coincide) {
                   return GetIntersectionRuleImpl(jCell, BndyRule, order, t.iFace);
                }
            }
            return null;
        }








        /// <summary>
        /// driver routine
        /// </summary>
        QuadRule GetIntersectionRule(int jCell, int order) {
            var qr = GetIntersectionRuleAsym(jCell, LevelSet1_Index, LevelSet2_Index, order);
            if(qr != null)
                return qr;
            var qr2 = GetIntersectionRuleAsym(jCell, LevelSet2_Index, LevelSet1_Index, order);
            if(qr2 != null)
                return qr2;

            var Empty = QuadRule.CreateBlank(RefElement, 1, RefElement.SpatialDimension);
            Empty.Nodes.LockForever();
            Empty.OrderOfPrecision = order;
            return Empty;
        }


        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if((!(mask is CellMask)) || (mask.MaskType != MaskType.Geometrical)) {
                throw new ArgumentException("expecting a geometrical cell mask", nameof(mask));
            }

            var ret = new List<IChunkRulePair<QuadRule>>();
            foreach(int jCell in mask.ItemEnum) {
                var qr = GetIntersectionRule(jCell, order);
                ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr));
            }

            return ret;
        }
    }
}
