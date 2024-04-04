using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;

namespace BoSSS.Foundation.XDG.Quadrature {
    public static class AlgoimFactories {
        #region Edge rules


        public static IQuadRuleFactory<CellBoundaryQuadRule> AlgoimGaussRule_EdgeVolume3D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder) {

            ISayeGaussEdgeRule rule = new SayeGaussRule_EdgeCube(
                _lsData,
                RootFinder,
                SayeGaussRule_Cube.QuadratureMode.PositiveVolume);
            return new SayeGaussEdgeRuleFactory(rule);
        }

        #endregion
    }

    //QuadRule RuleToRuleThemAll = QuadRule.CreateEmpty(refElement, count, spatialDim);
    //RuleToRuleThemAll.Nodes = new NodeSet(refElement, nodes, true);
    //RuleToRuleThemAll.Weights = weights;
    //        return RuleToRuleThemAll;



    class AlgoimEdgeRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {
        public RefElement RefElement => throw new NotImplementedException();


        public AlgoimEdgeRuleFactory() { }

        public int[] GetCachedRuleOrders() {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            throw new NotImplementedException();
        }
    }

}


