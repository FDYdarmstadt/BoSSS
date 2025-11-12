using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Configuration;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature.LevelSetOnEdge {
    internal class InterfaceIntegration4DoubleCut : InterfaceIntegration {
        public InterfaceIntegration4DoubleCut(RefElement _RefElement, LevelSetTracker.LevelSetData levelSetData, IQuadRuleFactory<QuadRule> __edgeSchemeForTrimmingLevelSet, EdgeMask __dblCutEdges)
            : base(_RefElement, levelSetData) //
        {
            grddat = levelSetData.GridDat;
            EdgeSchemeForTrimmingLevelSet = __edgeSchemeForTrimmingLevelSet;
            if(__dblCutEdges.MaskType != MaskType.Geometrical)
                throw new ArgumentException("expecting a geometrical mask", nameof(__dblCutEdges));
            DoubleCutEdgeMask = __dblCutEdges.GetBitMask();
        }

        readonly IGridData grddat;

        readonly IQuadRuleFactory<QuadRule> EdgeSchemeForTrimmingLevelSet;

        readonly BitArray DoubleCutEdgeMask;

        protected override QuadRule GetFaceRule(int jCell, int iFace, int intOrder) {
            int iEdge = grddat.GetEdgesForFace(jCell, iFace, out int InOrOut, out _);
            if(DoubleCutEdgeMask[iEdge]) {
                var Kref = grddat.iGeomCells.GetRefElement(jCell);

                // get edge rule
                var ruleSet = EdgeSchemeForTrimmingLevelSet.GetQuadRuleSet(new EdgeMask(grddat, Chunk.GetSingleElementChunk(iEdge), MaskType.Geometrical), intOrder);
                var edgeRule = ruleSet.First().Rule;

                // transform rule to cell
                int iTrafo = grddat.iGeomEdges.Edge2CellTrafoIndex[iEdge, InOrOut];
                var trafo_e2c = grddat.iGeomEdges.Edge2CellTrafos[iTrafo];
                var gramian = (trafo_e2c.Matrix.TransposeTo() * trafo_e2c.Matrix).Determinant().Sqrt();
                if(gramian <= 0)
                    throw new ArithmeticException();

                // transform from cell to face
                var VolNodes = trafo_e2c.Transform(edgeRule.Nodes);
                var trafo_c2f = Kref.GetInverseFaceTrafo(iFace);
                var faceRule = QuadRule.CreateBlank(Kref.FaceRefElement, edgeRule.NoOfNodes, Kref.FaceRefElement.SpatialDimension, true);
                trafo_c2f.Transform(VolNodes, faceRule.Nodes);

                // scale the weights (scale will be 1.0, but for snake of correctness?)
                faceRule.Weights.Set(edgeRule.Weights);
                faceRule.Weights.Scale(gramian * (1.0 / Kref.FaceTrafoGramianSqrt[iFace]));
                faceRule.OrderOfPrecision = edgeRule.OrderOfPrecision;
                faceRule.Nodes.LockForever();

                // return
                return faceRule;
            } else {
                return base.GetFaceRule(jCell, iFace, intOrder);
            }
        }


    }
}
