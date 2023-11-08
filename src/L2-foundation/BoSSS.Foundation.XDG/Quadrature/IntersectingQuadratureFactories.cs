using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using System;
using System.Collections.Generic;
using System.Text;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using BoSSS.Foundation.Grid.Classic;
using System.Data;
using System.Reflection;
using ilPSP;
using NUnit.Framework.Internal.Execution;
using System.Runtime.CompilerServices;
using System.ComponentModel;
using System.Collections;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG.Quadrature {
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
                levSet0,
                Symbol.Zero,
                levSet1,
                ToSymbol(jmp1)
            );
            return factory;
        }
    }

    internal class IntersectingRuleFactory : IQuadRuleFactory<QuadRule> {

        readonly IFunctionMap map;

        readonly LevelSetData alpha;
        readonly Symbol signAlpha;
        readonly LevelSetData beta;
        readonly Symbol signBeta;

        readonly Quadrater finder;

        public RefElement RefElement => map.Domain;

        public IntersectingRuleFactory( IFunctionMap translater, LevelSetData alpha, Symbol signAlpha, LevelSetData beta, Symbol signBeta) {
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
                for (int j = chunk.i0; j < chunk.JE; ++j) {
                    QuadRule rule = GetQuadRule(j, order);
                    if (rule.NoOfNodes == 0) {
                        rule = QuadRule.CreateEmpty(RefElement, 1, RefElement.SpatialDimension);
                        rule.Nodes.LockForever();
                    }
                    rules.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(j), rule));
                }
            }
            return rules;
        }

        QuadRule GetQuadRule(int j, int order) {
            HyperRectangle domain = map.Codomain(j);
            IScalarFunction a = map.MapFromDomainToCodomain( alpha, j);
            IScalarFunction b = map.MapFromDomainToCodomain( beta, j);
            (int nodeCount, int subdivisions) = Convert(order);
            QuadratureRule ruleQ = finder.FindRule(a, signAlpha, b, signBeta, domain, nodeCount, subdivisions);
            QuadRule q = map.MapFromCodomainToDomain(ruleQ, j);
            return q;
        }

        static (int nodeCount, int subdivisions) Convert(int order) {
            return (Math.Min(32, order), order/32);
        }
    }

    interface IFunctionMap {

        RefElement Domain { get; }

        HyperRectangle Codomain( int j);

        QuadRule MapFromCodomainToDomain(QuadratureRule rule, int j);

        IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int j);
    }

    class GlobalEdgeMapping : IFunctionMap {

        class Selection {

            int except;

            public Selection(int except) {
                this.except = except;
            }

            public int FromSubIndexToIndex(int i) {
                if (i < except) {
                    return i;
                } else {
                    return ++i;
                }
            }
        }

        GridData grid;

        public RefElement Domain { get; private set; }
        
        RefElement CellDomain;

        public GlobalEdgeMapping(GridData grid) {
            this.grid = grid;
            int dim = grid.SpatialDimension;

            switch (dim - 1) {
                case 1:
                Domain = Line.Instance;
                CellDomain = Square.Instance;
                break;
                case 2:
                Domain = Square.Instance;
                CellDomain = Cube.Instance;
                break;
                default:
                throw new NotImplementedException();
            }
            center = MultidimensionalArray.Create(1,dim);
        }

        MultidimensionalArray center;

        int EmptyDim(MultidimensionalArray center) {
            double max = double.MinValue;
            int emptyDim = -1;
            for (int i = 0; i < center.GetLength(1); ++i) {
                if (Math.Abs(center[0, i]) > max) {
                    emptyDim = i;
                    max = Math.Abs(center[0, i]);
                }
            }
            return emptyDim;
        }

        public HyperRectangle Codomain(int jEdge) {
            HyperRectangle codomain = new HyperRectangle(Domain.SpatialDimension);
            codomain.Dimension = Domain.SpatialDimension;

            int jCell = grid.Edges.CellIndices[jEdge, 0];
            byte iFace = grid.Edges.FaceIndices[jEdge, 0];
            NodeSet faceCenter = CellDomain.GetFaceCenter(iFace);
            int emptyDim = EmptyDim(faceCenter);

            grid.TransformLocal2Global(faceCenter, center, jCell);
            Selection edgeToCell = new Selection(emptyDim);

            for ( int i = 0; i < Domain.SpatialDimension; ++i) {
                int j = edgeToCell.FromSubIndexToIndex(i);
                codomain.Center[i] = center[0, j];
                codomain.Diameters[i] = grid.Cells.Transformation[jCell, j, j] * 2;
            }
            return codomain;
        }
        
        public QuadRule MapFromCodomainToDomain(QuadratureRule rule, int jEdge) {
            
            MultidimensionalArray globalNodes = MultidimensionalArray.Create(rule.Count, Domain.SpatialDimension);
            for (int i = 0; i < rule.Count; ++i) {
                for (int d = 0; d < Domain.SpatialDimension; ++d) {
                    globalNodes[i, d] = rule[i].Point[d];
                }
            }

            AffineTrafo t = FromCodomainToDomain(jEdge);
            QuadRule q = QuadRule.CreateEmpty(Domain, rule.Count, Domain.SpatialDimension);
            t.Transform(globalNodes, q.Nodes);

            double jacobianDeterminant = t.Matrix.Determinant();
            for (int i = 0; i < rule.Count; ++i) {
                q.Weights[i] = rule[i].Weight * jacobianDeterminant;
            }
            q.Nodes.LockForever();
            return q;
        }
        
        AffineTrafo FromCodomainToDomain(int jEdge) {
            AffineTrafo t = new AffineTrafo(Domain.SpatialDimension);
            HyperRectangle codomain = Codomain(jEdge);
            for (int i = 0; i < Domain.SpatialDimension; ++i) {
                t.Affine[i] = -codomain.Center[i] * 2.0 / codomain.Diameters[i];
                t.Matrix[i, i] = 2.0 / codomain.Diameters[i];
            }
            return t;
        }

        public IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int jEdge) {
            int jCell = grid.Edges.CellIndices[jEdge, 0];
            IScalarFunction globalLevelSet = new GlobalCellFunction(levelSet, jCell);
            IVectorFunction T = FromDomainToEmbeddedCodomain(jEdge);
            IScalarFunction alpha = new ScalarComposition(globalLevelSet, T);
            return alpha;
        }

        IVectorFunction FromDomainToEmbeddedCodomain(int jEdge) {
            int jCell = grid.Edges.CellIndices[jEdge, 0];
            byte iFace = grid.Edges.FaceIndices[jEdge, 0];
            NodeSet faceCenter = CellDomain.GetFaceCenter(iFace);
            int emptyDim = EmptyDim(faceCenter);

            grid.TransformLocal2Global(faceCenter, center, jCell);
            Selection edgeToCell = new Selection(emptyDim);

            Tensor1 affine = Tensor1.Zeros(CellDomain.SpatialDimension);
            Tensor2 m = Tensor2.Zeros(CellDomain.SpatialDimension, Domain.SpatialDimension);

            for (int i = 0; i < Domain.SpatialDimension; ++i) {
                int j = edgeToCell.FromSubIndexToIndex(i);
                affine[j] = center[0, j];
                m[j, i] = grid.Cells.Transformation[jCell, j, j] * 2;
            }
            affine[emptyDim] = center[0, emptyDim];

            LinearVectorPolynomial t = new LinearVectorPolynomial(affine, m);
            return t;
        }
    }

    class GlobalCellMapping : IFunctionMap {

        int dim;
        GridData grid;

        public RefElement Domain { get; private set; }

        public GlobalCellMapping(GridData grid) {
            this.grid = grid;
            dim = grid.SpatialDimension;
            
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
        }

        public HyperRectangle Codomain(int jCell) {
            HyperRectangle codomain = new HyperRectangle(grid.SpatialDimension);
            codomain.Dimension = grid.SpatialDimension;
            for (int i= 0; i < codomain.Dimension; ++i) {
                codomain.Center[i] = grid.Cells.CellCenter[jCell, i];
                codomain.Diameters[i] = grid.Cells.Transformation[jCell, i, i] * 2;
            }
            return codomain;
        }

        public QuadRule MapFromCodomainToDomain(QuadratureRule rule, int jCell) {
            MultidimensionalArray globalNodes = MultidimensionalArray.Create(rule.Count, dim);
            for (int i = 0; i < rule.Count; ++i) {
                for (int d = 0; d < dim; ++d) {
                    globalNodes[i, d] = rule[i].Point[d];
                }
            }
            
            QuadRule q = QuadRule.CreateEmpty(Domain, rule.Count, dim);
            grid.TransformGlobal2Local(globalNodes, q.Nodes, jCell, null);
            double jacobianDeterminant = grid.Cells.JacobiDet[jCell];
            for (int i = 0; i < rule.Count; ++i) {
                q.Weights[i] = rule[i].Weight / jacobianDeterminant;
            }
            q.Nodes.LockForever();
            return q;
        }

        public IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int jCell) {
            return new GlobalCellFunction(levelSet, jCell);
        }
    }

    class GlobalCellFunction : IScalarFunction {

        MultidimensionalArray evaluation;
        MultidimensionalArray gradient;
        MultidimensionalArray hessian;

        readonly ILevelSet levelSet;
        readonly GridData grid;
        readonly int jCell;

        public GlobalCellFunction(LevelSetData levelSet, int jCell) {
            this.jCell = jCell;
            this.M = levelSet.GridDat.SpatialDimension;
            this.grid = levelSet.GridDat;
            this.levelSet = levelSet.LevelSet;
            evaluation = MultidimensionalArray.Create(1, 1);
            gradient = MultidimensionalArray.Create(1, 1, M);
            hessian = MultidimensionalArray.Create(1, 1, M, M);
        }

        //Spatial Dimension
        public int M { get; private set; }

        public double Evaluate(Tensor1 x) {
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            X.LockForever();
            levelSet.Evaluate( jCell, 1, X, evaluation);
            return evaluation[0, 0];
        }

        public (double evaluation, Tensor1 gradient) EvaluateAndGradient(Tensor1 x) {
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            X.LockForever();
            levelSet.Evaluate(jCell, 1, X, evaluation);
            levelSet.EvaluateGradient(jCell, 1, X, gradient);
            return (evaluation[0, 0], ToTensor1(gradient));
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            X.LockForever();
            levelSet.Evaluate(jCell, 1, X, evaluation);
            levelSet.EvaluateGradient(jCell, 1, X, gradient);
            levelSet.EvaluateHessian(jCell, 1, X, hessian);
            return (evaluation[0, 0], ToTensor1(gradient), ToTensor2(hessian));
        }

        NodeSet FromGlobalToReferenceNodeSet(Tensor1 x) {
            MultidimensionalArray XGlobal = MultidimensionalArray.Create(1, x.M);
            for (int i = 0; i < x.M; ++i) {
                XGlobal[0, i] = x[i];
            }
            MultidimensionalArray XReference = MultidimensionalArray.Create(1, x.M);
            grid.TransformGlobal2Local(XGlobal, XReference, jCell, null);
            NodeSet n = new NodeSet(grid.Cells.GetRefElement(jCell), XReference, true);
            return n;
        }

        static Tensor1 ToTensor1(MultidimensionalArray x) {
            Tensor1 X = Tensor1.Zeros(x.Length);
            for (int i = 0; i < X.M; ++i) {
                X[i] = x[0, 0, i];
            }
            return X;
        }

        static Tensor2 ToTensor2(MultidimensionalArray x) {
            Tensor2 X = Tensor2.Zeros(x.GetLength(x.Dimension - 1));
            for (int i = 0; i < X.M; ++i) {
                for (int j = 0; j < X.M; ++j) {
                    X[i, j] = x[0, 0, i, j];
                }
            }
            return X;

        }
    }
}
