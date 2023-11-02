using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System;
using System.Collections.Generic;
using System.Text;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using IntersectingQuadrature;
using IntersectingQuadrature.TensorAnalysis;
using BoSSS.Foundation.Grid.Classic;
using System.Data;
using System.Reflection;
using ilPSP;
using NUnit.Framework.Internal.Execution;
using System.Runtime.CompilerServices;
using System.ComponentModel;

namespace BoSSS.Foundation.XDG.Quadrature {
    internal static class IntersectingFactories {

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
                new VolumeMapping(levSet0.GridDat.SpatialDimension),
                levSet0,
                ToSymbol(jmp0),
                levSet1,
                ToSymbol(jmp1)
            );
            return factory;
        }

        public static IQuadRuleFactory<QuadRule> Surface(LevelSetData levSet0, LevelSetData levSet1, JumpTypes jmp1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory(
                new VolumeMapping(levSet0.GridDat.SpatialDimension),
                levSet0,
                Symbol.Zero,
                levSet1,
                ToSymbol(jmp1)
            );
            return factory;
        }
        
        /*
        public IQuadRuleFactory<QuadRule> Intersection(int levSetIndex0, int levSetIndex1) {
            IntersectingRuleFactory factory = new IntersectingRuleFactory(
                new VolumeMapping(spatialDimension),
                levelSets[levSetIndex0],
                Symbol.Zero,
                levelSets[levSetIndex1],
                Symbol.Zero
            );
            return factory;
        }

        public IQuadRuleFactory<QuadRule> Edge(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1) {
            throw new NotImplementedException();
        }

        public IQuadRuleFactory<QuadRule> EdgePoint(int levSetIndex0, int levSetIndex1, JumpTypes jmp1) {
            throw new NotImplementedException();
        }
        */
    }

    internal class IntersectingRuleFactory : IQuadRuleFactory<QuadRule> {

        readonly IFunctionMap map;

        readonly LevelSetData alpha;
        readonly Symbol signAlpha;
        readonly LevelSetData beta;
        readonly Symbol signBeta;

        readonly Quadrater finder;

        public RefElement RefElement => map.Domain;

        public IntersectingRuleFactory( IFunctionMap translater,
                                        LevelSetData alpha,
                                        Symbol signAlpha,
                                        LevelSetData beta,
                                        Symbol signBeta) {
            this.map = translater;
            this.alpha = alpha;
            this.signAlpha = signAlpha;
            this.beta = beta;
            this.signBeta = signBeta;
            finder = new Quadrater();
        }

        public int[] GetCachedRuleOrders() {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            List<ChunkRulePair<QuadRule>> rules = new List<ChunkRulePair<QuadRule>>();
            foreach (Chunk chunk in mask) {
                for (int i = chunk.i0; i < chunk.JE; ++i) {
                    QuadRule rule = GetQuadRule(i, order);
                    rules.Add(new ChunkRulePair<QuadRule>(chunk, rule));
                }
            }
            return rules;
        }

        QuadRule GetQuadRule(int jCell, int order) {
            
            IScalarFunction a = map.MapToCodomain( alpha, jCell);
            IScalarFunction b = map.MapToCodomain( beta, jCell);
            (int nodeCount, int subdivisions) = Convert(order);
            HyperRectangle domain = map.Codomain( jCell);
            QuadratureRule ruleQ = finder.FindRule(a, signAlpha, b, signBeta, domain, nodeCount, subdivisions);
            QuadRule q = map.MapToDomain(ruleQ, jCell);
            return q;
        }

        static (int nodeCount, int subdivisions) Convert(int order) {
            return (1,0);
            //return (Math.Min(5, order), order/5);
        }
    }

    interface IFunctionMap {

        RefElement Domain { get; }

        HyperRectangle Codomain( int jCell);

        QuadRule MapToDomain(QuadratureRule rule, int jCell);

        IScalarFunction MapToCodomain(LevelSetData levelSetData, int jCell);
    }

    class VolumeMapping : IFunctionMap {

        int dim;

        readonly HyperRectangle codomain;

        public RefElement Domain { get; private set; }

        public VolumeMapping(int dim) {
            switch (dim) {
                case 1:
                Domain = Line.Instance;
                break;
                case 2:
                Domain = Square.Instance;
                break;
                case 3:
                Domain = Cube.Instance;
                break;
                default:
                throw new NotImplementedException();
            }
            codomain = new UnitCube(dim);
            this.dim = dim;
        }

        public HyperRectangle Codomain(int jCell) {
            return codomain;
        }

        public QuadRule MapToDomain(QuadratureRule rule, int jCell) {
            QuadRule q = QuadRule.CreateEmpty(Domain, rule.Count, dim);
            for(int i = 0; i < rule.Count; ++i) {
                for(int d = 0; d < dim; ++d) {
                    q.Nodes[i, d] = rule[i].Point[d];
                }
                q.Weights[i] = rule[i].Weight;
            }
            q.Nodes.LockForever();
            return q;
        }

        public IScalarFunction MapToCodomain(LevelSetData levelSetData, int jCell) {
            return new ReferenceElementFunction(Domain, levelSetData, jCell);
        }
    }

    class ReferenceElementFunction : IScalarFunction {

        readonly LevelSetData levelSet;

        readonly RefElement refElement;

        readonly int jCell;

        public ReferenceElementFunction(RefElement refElement, LevelSetData levelSetData, int jCell) {
            this.jCell = jCell;
            this.M = levelSetData.GridDat.SpatialDimension;
            this.levelSet = levelSetData;
            this.refElement = refElement;
        }

        public int M { get; private set; }

        public double Evaluate(Tensor1 x) {
            NodeSet X = ToNodeSet(x);
            X.LockForever();
            MultidimensionalArray evaluation = levelSet.GetLevSetValues(X, jCell,1);
            return evaluation[0, 0];
        }

        public (double evaluation, Tensor1 gradient) EvaluateAndGradient(Tensor1 x) {
            NodeSet X = ToNodeSet(x);
            X.LockForever();
            MultidimensionalArray evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            MultidimensionalArray gradient = levelSet.GetLevelSetReferenceGradients(X, jCell, 1);
            gradient.Scale(1.0 / levelSet.GridDat.Cells.JacobiDet[jCell]);
            return (evaluation[0, 0], ToTensor1(gradient));
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            NodeSet X = ToNodeSet(x);
            X.LockForever();
            MultidimensionalArray evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            MultidimensionalArray gradient = levelSet.GetLevelSetReferenceGradients(X, jCell, 1);
            gradient.Scale(1.0 / levelSet.GridDat.Cells.JacobiDet[jCell]);
            MultidimensionalArray hessian = levelSet.GetLevelSetReferenceHessian(X, jCell, 1);
            hessian.Scale(1.0 / levelSet.GridDat.Cells.JacobiDet[jCell]);
            return (evaluation[0, 0], ToTensor1(gradient), ToTensor2(hessian));
        }

        NodeSet ToNodeSet(Tensor1 x) {
            MultidimensionalArray X = MultidimensionalArray.Create( 1, x.M);
            for (int i = 0; i < x.M; ++i) {
                X[0, i] = x[i];
            }
            NodeSet n = new NodeSet(refElement, X, true);
            return n;
        }

        static Tensor1 ToTensor1(MultidimensionalArray x) {
            Tensor1 X = Tensor1.Zeros(x.Length);
            for(int i= 0; i < X.M; ++i) {
                X[i] = x[0,0,i];
            }
            return X;
        }

        static Tensor2 ToTensor2(MultidimensionalArray x) {
            Tensor2 X = Tensor2.Zeros(x.GetLength(x.Dimension - 1));
            for (int i = 0; i < X.M; ++i) {
                for (int j = 0; j < X.M; ++j) {
                    X[i,j] = x[0, 0, i, j];
                }
            }
            return X;
            
        }
    }

}
