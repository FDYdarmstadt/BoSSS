using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System;
using System.Collections.Generic;
using System.Linq;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    public sealed class IntersectingSubcell {
        public double[] Center;
        public double[] Diameters;
        public Func<double[], double> LevelSet0;
        public Func<double[], double> LevelSet1;
    }

    
    class DelegateIntersectingRuleFactory : IQuadRuleFactory<QuadRule> {
        readonly RefElement Kref;
        readonly Func<int, IEnumerable<IntersectingSubcell>> subcellProvider;
        readonly Symbol sign0;
        readonly Symbol sign1;
        readonly IQuadrater finder;

        public DelegateIntersectingRuleFactory(RefElement Kref, Func<int, IEnumerable<IntersectingSubcell>> subcellProvider, Symbol sign0, Symbol sign1) {
            this.Kref = Kref;
            this.subcellProvider = subcellProvider;
            this.sign0 = sign0;
            this.sign1 = sign1;
            finder = IntersectingQuadrature.Methods.Create(5);
        }

        public RefElement RefElement {
            get {
                return Kref;
            }
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            var rules = new List<ChunkRulePair<QuadRule>>();
            foreach(int j in mask.ItemEnum) {
                rules.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(j), GetQuadRule(j, order)));
            }
            return rules;
        }

        QuadRule GetQuadRule(int jCell, int order) {
            var nodes = new List<double[]>();
            var weights = new List<double>();
            (int nodeCount, int subdivisions) = Convert(order);

            foreach(var subcell in subcellProvider(jCell)) {
                if(subcell == null)
                    continue;
                if(IsDegenerate(subcell))
                    continue;

                HyperRectangle domain = new HyperRectangle(Kref.SpatialDimension);
                domain.Dimension = Kref.SpatialDimension;
                for(int d = 0; d < Kref.SpatialDimension; d++) {
                    domain.Center[d] = subcell.Center[d];
                    domain.Diameters[d] = subcell.Diameters[d];
                }

                QuadratureRule subRule = finder.FindRule(
                    new DelegateScalarFunction(Kref.SpatialDimension, subcell.LevelSet0),
                    sign0,
                    new DelegateScalarFunction(Kref.SpatialDimension, subcell.LevelSet1),
                    sign1,
                    domain,
                    nodeCount,
                    subdivisions);

                for(int i = 0; i < subRule.Count; i++) {
                    double[] node = new double[Kref.SpatialDimension];
                    for(int d = 0; d < Kref.SpatialDimension; d++)
                        node[d] = subRule[i].Point[d];
                    nodes.Add(node);
                    weights.Add(subRule[i].Weight);
                }
            }

            if(nodes.Count == 0) {
                var empty = QuadRule.CreateBlank(Kref, 1, Kref.SpatialDimension);
                empty.OrderOfPrecision = order;
                empty.Nodes.LockForever();
                return empty;
            }

            QuadRule rule = QuadRule.CreateBlank(Kref, nodes.Count, Kref.SpatialDimension);
            rule.OrderOfPrecision = order;
            for(int i = 0; i < nodes.Count; i++) {
                for(int d = 0; d < Kref.SpatialDimension; d++)
                    rule.Nodes[i, d] = nodes[i][d];
                rule.Weights[i] = weights[i];
            }
            rule.Nodes.LockForever();
            return rule;
        }

        static bool IsDegenerate(IntersectingSubcell subcell) {
            if(subcell.Center == null || subcell.Diameters == null || subcell.LevelSet0 == null || subcell.LevelSet1 == null)
                throw new ArgumentException("Invalid intersecting subcell.");
            return subcell.Diameters.Any(d => Math.Abs(d) <= 1.0e-14);
        }

        static (int nodeCount, int subdivisions) Convert(int order) {
            return (Math.Min(32, Math.Max(1, order)), order / 32);
        }
    }

    class DelegateScalarFunction : IScalarFunction {
        readonly Func<double[], double> function;

        public DelegateScalarFunction(int dimension, Func<double[], double> function) {
            M = dimension;
            this.function = function;
        }

        public int M {
            get;
        }

        public double Evaluate(Tensor1 x) {
            return function(ToArray(x));
        }

        public (double evaluation, Tensor1 gradient) EvaluateAndGradient(Tensor1 x) {
            double[] point = ToArray(x);
            double evaluation = function(point);
            Tensor1 gradient = Tensor1.Zeros(M);
            double h = 1.0e-6;
            for(int d = 0; d < M; d++) {
                double old = point[d];
                point[d] = old + h;
                double plus = function(point);
                point[d] = old - h;
                double minus = function(point);
                point[d] = old;
                gradient[d] = (plus - minus) / (2.0 * h);
            }
            return (evaluation, gradient);
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            var (evaluation, gradient) = EvaluateAndGradient(x);
            Tensor2 hessian = Tensor2.Zeros(M);
            double[] point = ToArray(x);
            double h = 1.0e-5;
            double center = function(point);

            for(int i = 0; i < M; i++) {
                double oldI = point[i];
                point[i] = oldI + h;
                double plus = function(point);
                point[i] = oldI - h;
                double minus = function(point);
                point[i] = oldI;
                hessian[i, i] = (plus - 2.0 * center + minus) / (h * h);

                for(int j = i + 1; j < M; j++) {
                    double oldJ = point[j];
                    point[i] = oldI + h;
                    point[j] = oldJ + h;
                    double pp = function(point);
                    point[j] = oldJ - h;
                    double pm = function(point);
                    point[i] = oldI - h;
                    point[j] = oldJ + h;
                    double mp = function(point);
                    point[j] = oldJ - h;
                    double mm = function(point);
                    point[i] = oldI;
                    point[j] = oldJ;
                    double mixed = (pp - pm - mp + mm) / (4.0 * h * h);
                    hessian[i, j] = mixed;
                    hessian[j, i] = mixed;
                }
            }

            return (evaluation, gradient, hessian);
        }

        double[] ToArray(Tensor1 x) {
            double[] point = new double[M];
            for(int d = 0; d < M; d++)
                point[d] = x[d];
            return point;
        }
    }



}
