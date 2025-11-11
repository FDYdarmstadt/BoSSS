using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    internal class IntersectingRuleFactory : IQuadRuleFactory<QuadRule> {

        readonly IFunctionMap map;

        readonly LevelSetData alpha;
        readonly Symbol signAlpha;
        readonly LevelSetData beta;
        readonly Symbol signBeta;

        readonly IQuadrater finder;

        public RefElement RefElement => map.Domain;

        public IntersectingRuleFactory(IFunctionMap translater, LevelSetData alpha, Symbol signAlpha, LevelSetData beta, Symbol signBeta) {
            this.map = translater;
            this.alpha = alpha;
            this.signAlpha = signAlpha;
            this.beta = beta;
            this.signBeta = signBeta;
            finder = IntersectingQuadrature.Methods.Create(5);
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            List<ChunkRulePair<QuadRule>> rules = new List<ChunkRulePair<QuadRule>>();

            foreach(int j in mask.ItemEnum) {
                try {
                    QuadRule rule = GetQuadRule(j, order);
                    if(rule.NoOfNodes == 0) {
                        rule = QuadRule.CreateBlank(RefElement, 1, RefElement.SpatialDimension);
                        rule.OrderOfPrecision = order;
                        rule.Nodes.LockForever();
                    }
                    rule.OrderOfPrecision = order;
                    rules.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(j), rule));
                } catch(Exception e) {

                    var grid = mask.GridData.Grid as GridCommons;
                    if(grid != null) {
                        Console.Error.WriteLine("Cell: ");
                        Console.Error.WriteLine(grid.Cells[j].ToString());

                        var lsAlpah = (alpha.LevelSet as DGField);
                        Console.Error.WriteLine("DG coordinates of level set alpha: " + lsAlpah.Coordinates.GetRow(j).ToConcatString("[", ", ", "]"));

                        var lsBeta = (beta.LevelSet as DGField);
                        Console.Error.WriteLine("DG coordinates of level set beta: " + lsBeta.Coordinates.GetRow(j).ToConcatString("[", ", ", "]"));
                    }

                    throw e;
                }
            }
            return rules;
        }

        QuadRule GetQuadRule(int j, int order) {
            HyperRectangle domain = map.Codomain(j);
            IScalarFunction a = map.MapFromDomainToCodomain(alpha, j);
            IScalarFunction b = map.MapFromDomainToCodomain(beta, j);
            (int nodeCount, int subdivisions) = Convert(order);

            QuadratureRule ruleQ;
            try {
                ruleQ = finder.FindRule(a, signAlpha, b, signBeta, domain, nodeCount, subdivisions);
            } catch(Exception e) {
                Console.WriteLine(e);
                throw new ArithmeticException($"IntersectingQuadrature failed for input nodeCount={nodeCount}, subdivisions={subdivisions}, map={map} ");
            }

            QuadRule q = map.MapFromCodomainToDomain(ruleQ, j);
            return q;
        }

        static (int nodeCount, int subdivisions) Convert(int order) {
            return (Math.Min(32, Math.Max(1, order)), order / 32);
        }
    }
}
