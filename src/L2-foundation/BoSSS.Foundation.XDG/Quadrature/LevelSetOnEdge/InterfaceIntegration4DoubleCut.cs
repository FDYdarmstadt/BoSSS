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
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature.LevelSetOnEdge {
    internal class InterfaceIntegration4DoubleCut : InterfaceIntegration {
        public InterfaceIntegration4DoubleCut(RefElement _RefElement, LevelSetTracker.LevelSetData levelSetData, 
            IQuadRuleFactory<QuadRule> __edgeSchemeForTrimmingLevelSet, EdgeMask __dblCutEdges,
            LevelSetTracker.LevelSetData[] __levelSetDataS, LevelSetSignCode[] __allSignCodes)
            : base(_RefElement, levelSetData) //
        {
            grddat = levelSetData.GridDat;
            EdgeSchemeForTrimmingLevelSet = __edgeSchemeForTrimmingLevelSet;
            if(__dblCutEdges.MaskType != MaskType.Geometrical)
                throw new ArgumentException("expecting a geometrical mask", nameof(__dblCutEdges));
            DoubleCutEdgeMask = __dblCutEdges.GetBitMask();
            LevelSetDataS = __levelSetDataS;
            AllSignCodes = __allSignCodes;

            IndexOfCoincidingLevelSet = levelSetData.LevelSetIndex;
            if(__levelSetDataS.IndexOf(levelSetData, object.ReferenceEquals) != IndexOfCoincidingLevelSet) {
                throw new ArgumentException("expecting the `LevelSetIndex` to match", nameof(__levelSetDataS));
            }
        }

        readonly int IndexOfCoincidingLevelSet;

        readonly IGridData grddat;

        readonly IQuadRuleFactory<QuadRule> EdgeSchemeForTrimmingLevelSet;

        readonly BitArray DoubleCutEdgeMask;

        readonly LevelSetTracker.LevelSetData[] LevelSetDataS;
        readonly LevelSetSignCode[] AllSignCodes;

        bool IsSpeciesPresentAtNodes(NodeSet Ns, int jCell) {
            int NoOfLevSets = LevelSetDataS.Length;
            var LevelSetValues = MultidimensionalArray.Create(NoOfLevSets, Ns.NoOfNodes);

            var Xglobal = LevelSetDataS[0].GridDat.GlobalNodes.GetValue_Cell(Ns, jCell);
            if(Xglobal[0, 0].Abs() <= 0.001 && Xglobal[0, 1] < 0)
                Console.WriteLine("dsjlkjdskljskjcjckjxhjdhhjdsacsahdchbjb");


            for(int iLevSet = 0; iLevSet < NoOfLevSets; iLevSet++) {
                if(iLevSet != IndexOfCoincidingLevelSet) {
                    LevelSetValues.ExtractSubArrayShallow(iLevSet, -1).Set(
                        LevelSetDataS[iLevSet].GetLevSetValues(Ns, jCell, 1).ExtractSubArrayShallow(0, -1)
                        );
                } else {
                    // level-set is exactly 0, nothing to do, since buffer is freshly allocated.
                }
            }

            //var codes = new LevelSetSignCode[NoOfLevSets];
            bool completelyEmty = true, completelyFull = true;
            for(int k = 0; k < Ns.NoOfNodes; k++) {
                var code1 = LevelSetSignCode.ComputeLevelSetBytecode(LevelSetValues.GetColumn(k)); code1.SetSign(IndexOfCoincidingLevelSet, -1);
                var code2 = LevelSetSignCode.ComputeLevelSetBytecode(LevelSetValues.GetColumn(k)); code2.SetSign(IndexOfCoincidingLevelSet, +1);
                bool nonvoid = Array.IndexOf(AllSignCodes, code1) >= 0 || Array.IndexOf(AllSignCodes, code2) >= 0;

                completelyFull &= nonvoid;
                completelyEmty &= (!nonvoid);
            }

            if(completelyFull == completelyEmty)
                throw new ArithmeticException("cannot decide");

            return completelyFull;
        }


        protected override QuadRule GetFaceRule(int jCell, int iFace, int intOrder) {
            int iEdge = grddat.GetEdgesForFace(jCell, iFace, out int InOrOut, out _);
            if(DoubleCutEdgeMask[iEdge]) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // the face/edge is by the trimming level-set;
                // => the edge rule of the trimming level-set is actually the right quadrature for the coinciding level-set
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                var Kref = base.RefElement;
                Debug.Assert(Kref == grddat.iGeomCells.GetRefElement(jCell));

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
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // the face/edge which coincides with the level-set is either 
                // * completely full 
                //     -- or --
                // * completely empty
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                var testNode = base.RefElement.GetFaceCenter(iFace);

                if(IsSpeciesPresentAtNodes(testNode, jCell)) {
                    return base.GetFaceRule(jCell, iFace, intOrder);
                } else {
                    var empty = QuadRule.CreateBlank(base.RefElement.FaceRefElement, 1, base.RefElement.FaceRefElement.SpatialDimension);
                    empty.Nodes.LockForever();
                    empty.OrderOfPrecision = intOrder;
                    return empty;
                }
            }
        }


    }
}
