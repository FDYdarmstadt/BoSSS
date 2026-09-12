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
    /*
    public sealed class SubCell {
        public double[] Center;
        public double[] Diameters;
        public Func<double[], double> LevelSet0;
        public Func<double[], double[]> LevelSet0Gradient;
        public Func<double[], double[,]> LevelSet0Hessian;
        public Func<double[], double> LevelSet1;
        public Func<double[], double[]> LevelSet1Gradient;
        public Func<double[], double[,]> LevelSet1Hessian;
    }


    class DelegateIntersectingRuleFactory : IQuadRuleFactory<QuadRule> {
        readonly RefElement Kref;
        readonly Func<int, IEnumerable<SubCell>> subcellProvider;
        readonly Symbol sign0;
        readonly Symbol sign1;
        readonly IQuadrater finder;

        public DelegateIntersectingRuleFactory(RefElement Kref, Func<int, IEnumerable<SubCell>> subcellProvider, Symbol sign0, Symbol sign1) {
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

            var subCells = subcellProvider(jCell);
            foreach(var subcell in subCells) {
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
                    new DelegateScalarFunction(Kref.SpatialDimension, subcell.LevelSet0, subcell.LevelSet0Gradient, subcell.LevelSet0Hessian),
                    sign0,
                    new DelegateScalarFunction(Kref.SpatialDimension, subcell.LevelSet1, subcell.LevelSet1Gradient, subcell.LevelSet1Hessian),
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

        static bool IsDegenerate(SubCell subcell) {
            if(subcell.Center == null || subcell.Diameters == null
                || subcell.LevelSet0 == null || subcell.LevelSet0Gradient == null || subcell.LevelSet0Hessian == null
                || subcell.LevelSet1 == null || subcell.LevelSet1Gradient == null || subcell.LevelSet1Hessian == null)
                throw new ArgumentException("Invalid intersecting subcell.");
            return subcell.Diameters.Any(d => Math.Abs(d) <= 1.0e-14);
        }

        static (int nodeCount, int subdivisions) Convert(int order) {
            return (Math.Min(32, Math.Max(1, order)), order / 32);
        }
    }


    /// <summary>
    /// just a fucking wrapper around another wrapper
    /// </summary>
    class DelegateScalarFunction : IScalarFunction {
        readonly Func<double[], double> function;
        readonly Func<double[], double[]> gradientFunction;
        readonly Func<double[], double[,]> hessianFunction;

        public DelegateScalarFunction(int dimension, Func<double[], double> function, Func<double[], double[]> gradientFunction, Func<double[], double[,]> hessianFunction) {
            M = dimension;
            this.function = function ?? throw new ArgumentNullException(nameof(function));
            this.gradientFunction = gradientFunction ?? throw new ArgumentNullException(nameof(gradientFunction));
            this.hessianFunction = hessianFunction ?? throw new ArgumentNullException(nameof(hessianFunction));
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
            double[] grad = gradientFunction(point);
            if(grad == null || grad.Length != M)
                throw new InvalidOperationException($"Gradient delegate must return a vector of length {M}.");
            Tensor1 gradient = Tensor1.Zeros(M);
            for(int d = 0; d < M; d++)
                gradient[d] = grad[d];
            return (evaluation, gradient);
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            double[] point = ToArray(x);
            double evaluation = function(point);

            double[] grad = gradientFunction(point);
            if(grad == null || grad.Length != M)
                throw new InvalidOperationException($"Gradient delegate must return a vector of length {M}.");
            Tensor1 gradient = Tensor1.Zeros(M);
            for(int d = 0; d < M; d++)
                gradient[d] = grad[d];

            double[,] hess = hessianFunction(point);
            if(hess == null || hess.GetLength(0) != M || hess.GetLength(1) != M)
                throw new InvalidOperationException($"Hessian delegate must return a {M}x{M} matrix.");
            Tensor2 hessian = Tensor2.Zeros(M);
            for(int i = 0; i < M; i++) {
                for(int j = 0; j < M; j++)
                    hessian[i, j] = hess[i, j];
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


    */
}
