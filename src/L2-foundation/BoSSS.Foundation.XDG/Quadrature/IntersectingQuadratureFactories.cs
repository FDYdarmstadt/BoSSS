using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System;
using System.Collections.Generic;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    internal static class IntersectingQuadratureFactories {

        static Symbol ToSymbol(JumpTypes sign) {
            switch (sign) {
                case JumpTypes.OneMinusHeaviside:
                return Symbol.Minus;
                case JumpTypes.Heaviside:
                return Symbol.Plus;
                default:
                throw new NotSupportedException();
            }
        }

        public static IQuadRuleFactory<QuadRule> Volume(LevelSetData levSet0, JumpTypes jmp0, LevelSetData levSet1, JumpTypes jmp1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory( 
                new GlobalCellMapping(levSet0.GridDat),
                levSet0,
                ToSymbol(jmp0),
                levSet1,
                ToSymbol(jmp1)
            );
            return factory;
        }

        public static IQuadRuleFactory<QuadRule> Surface(LevelSetData levSet0, LevelSetData levSet1, JumpTypes jmp1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory(
                new GlobalCellMapping(levSet0.GridDat),
                levSet0,
                Symbol.Zero,
                levSet1,
                ToSymbol(jmp1)
            );
            return factory;
        }
        
        public static IQuadRuleFactory<QuadRule> Intersection(LevelSetData levSet0, LevelSetData levSet1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory(
                new GlobalCellMapping(levSet0.GridDat),
                levSet0,
                Symbol.Zero,
                levSet1,
                Symbol.Zero
            );
            return factory;
        }
        
        public static IQuadRuleFactory<QuadRule> Edge(LevelSetData levSet0, JumpTypes jmp0, LevelSetData levSet1, JumpTypes jmp1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory(
                new GlobalEdgeMapping(levSet0.GridDat),
                levSet0,
                ToSymbol(jmp0),
                levSet1,
                ToSymbol(jmp1)
            );
            return factory;
        }

        public static IQuadRuleFactory<QuadRule> EdgePoint(LevelSetData levSet0, LevelSetData levSet1, JumpTypes jmp1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory(
                new GlobalEdgeMapping(levSet0.GridDat),
                levSet1,
                ToSymbol(jmp1),
                levSet0,
                Symbol.Zero
            );
            return factory;
        }
    }

    

}
