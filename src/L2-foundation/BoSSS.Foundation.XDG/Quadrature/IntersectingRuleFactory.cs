using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using static BoSSS.Foundation.XDG.LevelSetTracker;


namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {

    /*
    public static class Report {
        public static void R() {
            Console.WriteLine("-----------------------------------------------");
            Console.WriteLine($"all calls: {IntersectingRuleFactory.calls.Elapsed.TotalSeconds:g7}");
            Console.WriteLine($"  opt val: {GlobalCellFunctionOptimized.stwV.Elapsed.TotalSeconds:g7} \t(calls: {GlobalCellFunctionOptimized.EvalCounterV})");
            Console.WriteLine($" opt grad: {GlobalCellFunctionOptimized.stwG.Elapsed.TotalSeconds:g7} \t(calls: {GlobalCellFunctionOptimized.EvalCounterG})");
            Console.WriteLine($" opt hess: {GlobalCellFunctionOptimized.stwH.Elapsed.TotalSeconds:g7} \t(calls: {GlobalCellFunctionOptimized.EvalCounterH})");
            Console.WriteLine($" opt ctor: {GlobalCellFunctionOptimized.stwC.Elapsed.TotalSeconds:g7}");

            Console.WriteLine($"  ref val: {GlobalCellFunction.stwV.Elapsed.TotalSeconds:g7} \t(calls: {GlobalCellFunction.EvalCounterV})");
            Console.WriteLine($" ref grad: {GlobalCellFunction.stwG.Elapsed.TotalSeconds:g7} \t(calls: {GlobalCellFunction.EvalCounterG})");
            Console.WriteLine($" ref hess: {GlobalCellFunction.stwH.Elapsed.TotalSeconds:g7} \t(calls: {GlobalCellFunction.EvalCounterH})");

            Console.WriteLine("-----------------------------------------------");
        }
    }
    */

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
                QuadRule rule;
                bool restore = ilPSP.Environment.StdOut.surpressStream0;
                try {
                    ilPSP.Environment.StdOut.surpressStream0 = true;
                    rule = GetQuadRule(j, order);
                    if(rule.NoOfNodes == 0) {
                        rule = QuadRule.CreateBlank(RefElement, 1, RefElement.SpatialDimension);
                        rule.OrderOfPrecision = order;
                        rule.Nodes.LockForever();
                    }
                } catch(Exception e) {
                    Console.Error.WriteLine("Intersecting Quadrature, cell " + j + ": " + e.Message + ", defaulting to empty rule");

                    /*


                    var grid = mask.GridData.Grid as GridCommons;
                    if(grid != null) {

                        XQuadFactoryHelper.Plot();


                        Console.Error.WriteLine("Cell: ");
                        Console.Error.WriteLine(grid.Cells[j].ToString());

                        var lsAlpah = (alpha.LevelSet as DGField);
                        Console.Error.WriteLine("DG coordinates of level set alpha: " + lsAlpah.Coordinates.GetRow(j).ToConcatString("[", ", ", "]"));

                        var lsBeta = (beta.LevelSet as DGField);
                        Console.Error.WriteLine("DG coordinates of level set beta: " + lsBeta.Coordinates.GetRow(j).ToConcatString("[", ", ", "]"));
                    }

                    */
                    //throw e;

                    rule = QuadRule.CreateBlank(RefElement, 1, RefElement.SpatialDimension);
                    rule.OrderOfPrecision = order;
                    rule.Nodes.LockForever();
                } finally {
                    ilPSP.Environment.StdOut.surpressStream0 = restore;
                }
                rule.OrderOfPrecision = order;
                rules.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(j), rule));

            }
            return rules;
        }

        //public static Stopwatch calls = new Stopwatch();

        QuadRule GetQuadRule(int j, int order) {
            HyperRectangle domain = map.Codomain(j);
            //calls.Start();
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
            //calls.Stop();
            return q;
        }

        static (int nodeCount, int subdivisions) Convert(int order) {
            return (Math.Min(32, Math.Max(1, order)), order / 32);
        }


        
    }
}
