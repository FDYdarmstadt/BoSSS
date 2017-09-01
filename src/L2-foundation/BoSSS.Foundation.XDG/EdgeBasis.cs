using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Foundation.Grid;
using System.Diagnostics;

namespace BoSSS.Foundation.XDG.Quadrature.Unknown {
    
    public class EdgeBasis : Basis {

        private int localEdge;

        private EdgeNodeCache edgeNodeCache;

        private MonomialCache monomialCache;

        private ValueCache valueCache;

        private GradientCache gradientCache;

        public EdgeBasis(Context context, int degree, int localEdge)
            : base(context, context.Grid.GridSimplex.EdgeSimplex, degree) {
            ConstructorCommon(localEdge);
        }

        public EdgeBasis(Context context, IEnumerable<Polynomial> polynomials, int localEdge)
            : base(context, context.Grid.GridSimplex.EdgeSimplex, polynomials) {
            ConstructorCommon(localEdge);
        }

        private void ConstructorCommon(int localEdge) {
            this.localEdge = localEdge;
            edgeNodeCache = new EdgeNodeCache(Context.NSC, this, localEdge);
            monomialCache = new MonomialCache(Context.NSC, this);
            valueCache = new ValueCache(Context.NSC, this);
            gradientCache = new GradientCache(Context.NSC, this);
        }

        public override MultidimensionalArray Evaluate(int NodeSetIndex) {
            return valueCache.GetValues(NodeSetIndex);
        }

        public override MultidimensionalArray EvaluateGradient(int NodeSetIndex) {
            return gradientCache.GetValues(NodeSetIndex);
        }

        public override MultidimensionalArray Evaluate2ndDeriv(int NodeSetIndex) {
            throw new NotImplementedException();
        }

        public override bool IsSubBasis(Basis other) {
            EdgeBasis otherEdgeBasis = other as EdgeBasis;
            if (otherEdgeBasis == null) {
                return false;
            }

            if (otherEdgeBasis.localEdge != this.localEdge) {
                return false;
            }

            return base.IsSubBasis(other);
        }

        public override bool Equals(object obj) {
            EdgeBasis otherEdgeBasis = obj as EdgeBasis;
            if (otherEdgeBasis == null) {
                return false;
            }

            if (otherEdgeBasis.localEdge != this.localEdge) {
                return false;
            }

            return base.Equals(obj);
        }

        public override int GetHashCode() {
            return base.GetHashCode() ^ localEdge;
        }

        class EdgeNodeCache : NodeSetController.ValueCacheBase<MultidimensionalArray> {

            private EdgeBasis edgeBasis;

            private int localEdge;

            public EdgeNodeCache(NodeSetController owner, EdgeBasis edgeBasis, int localEdge)
                : base(owner) {
                this.edgeBasis = edgeBasis;
                this.localEdge = localEdge;
            }

            protected override MultidimensionalArray Evaluate(MultidimensionalArray Nodes) {
                Debug.Assert(
                    edgeBasis.Context.Grid.GridSimplex.EdgeSimplex == edgeBasis.Simplex,
                    "An edge basis must live on the edge of a simplex");

                MultidimensionalArray edgeNodes = MultidimensionalArray.Create(
                    Nodes.GetLength(0), edgeBasis.Simplex.SpatialDimension);
                edgeBasis.Context.Grid.GridSimplex.VolumeToEdgeCoordinates(
                    localEdge, Nodes, edgeNodes);
                return edgeNodes;
            }
        }

        class MonomialCache : NodeSetController.ValueCacheBase<MultidimensionalArray> {

            private EdgeBasis edgeBasis;

            private MultidimensionalArray edgeNodes;

            public MonomialCache(NodeSetController owner, EdgeBasis edgeBasis)
                : base(owner) {
                this.edgeBasis = edgeBasis;
            }

            public override MultidimensionalArray GetValues(int NodeSetIndex) {
                edgeNodes = edgeBasis.edgeNodeCache.GetValues(NodeSetIndex);
                MultidimensionalArray result = base.GetValues(NodeSetIndex);
                edgeNodes = null;

                return result;
            }

            protected override MultidimensionalArray Evaluate(MultidimensionalArray Nodes) {
                return Polynomial.GetMonomials(
                    edgeNodes, edgeBasis.Simplex.SpatialDimension, edgeBasis.Degree);
            }
        }

        class ValueCache : NodeSetController.ValueCacheBase<MultidimensionalArray> {

            private EdgeBasis edgeBasis;

            private MultidimensionalArray edgeNodes;

            private MultidimensionalArray monomials;

            public ValueCache(NodeSetController owner, EdgeBasis edgeBasis)
                : base(owner) {
                this.edgeBasis = edgeBasis;
            }

            public override MultidimensionalArray GetValues(int NodeSetIndex) {
                edgeNodes = edgeBasis.edgeNodeCache.GetValues(NodeSetIndex);
                monomials = edgeBasis.monomialCache.GetValues(NodeSetIndex);
                MultidimensionalArray result = base.GetValues(NodeSetIndex);
                edgeNodes = null;
                monomials = null;

                return result;
            }

            protected override MultidimensionalArray Evaluate(MultidimensionalArray Nodes) {
                MultidimensionalArray ret = MultidimensionalArray.Create(
                    edgeNodes.GetLength(0), edgeBasis.Polynomials.Length);

                for (int m = 0; m < edgeBasis.Polynomials.Length; m++) {
                    edgeBasis.Polynomials[m].Evaluate(
                        ret.ExtractSubArrayShallow(-1, m),
                        edgeNodes,
                        monomials);
                }

                return ret;
            }
        }

        class GradientCache : NodeSetController.ValueCacheBase<MultidimensionalArray> {

            private EdgeBasis edgeBasis;

            private MultidimensionalArray edgeNodes;

            private MultidimensionalArray monomials;

            internal GradientCache(NodeSetController owner, EdgeBasis edgeBasis)
                : base(owner) {
                this.edgeBasis = edgeBasis;
            }

            public override MultidimensionalArray GetValues(int NodeSetIndex) {
                edgeNodes = edgeBasis.edgeNodeCache.GetValues(NodeSetIndex);
                monomials = edgeBasis.monomialCache.GetValues(NodeSetIndex);
                MultidimensionalArray result = base.GetValues(NodeSetIndex);
                edgeNodes = null;
                monomials = null;

                return result;
            }

            protected override MultidimensionalArray Evaluate(MultidimensionalArray Nodes) {
                int noOfNodes = edgeNodes.GetLength(0);
                int D = edgeNodes.GetLength(1);
                MultidimensionalArray ret = MultidimensionalArray.Create(
                    noOfNodes, edgeBasis.Polynomials.Length, D);
                for (int i = 0; i < edgeBasis.Polynomials.Length; i++) {
                    for (int j = 0; j < D; j++) {
                        if (edgeBasis.DerivativePolynomials[i, j].Coeff.Any(d => d != 0.0)) {
                            edgeBasis.DerivativePolynomials[i, j].Evaluate(
                                ret.ExtractSubArrayShallow(-1, i, j),
                                edgeNodes,
                                monomials);
                        }
                    }
                }

                return ret;
            }
        }
    }
}
