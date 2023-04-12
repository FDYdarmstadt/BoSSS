using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using log4net.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using IntersectingQuadrature;
using IntersectingQuadrature.TensorAnalysis;
using BoSSS.Foundation.XDG.Quadrature.Subdivision;
using NUnit.Framework.Interfaces;
using System.Drawing;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Platform;
using System.ComponentModel;

namespace BoSSS.Foundation.XDG.Quadrature
{
    public class MultiLevelSetSayeFactoryCreator
    {
        LevelSetCombination[] phis;

        GridData grid;

        LevelSetData[] levelSets;

        public MultiLevelSetSayeFactoryCreator(LevelSetData[] levelSets)
        {
            this.levelSets = levelSets;
            LevelSet levelSet0 = (LevelSet)levelSets[0].LevelSet;
            LevelSet levelSet1 = (LevelSet)levelSets[1].LevelSet;
            grid = levelSets[0].GridDat;
            phis = new LevelSetCombination[]
            {
                new LevelSetCombination( new CombinedID
                    {
                        LevSet0 = 0,
                        Jmp0 = JumpTypes.Heaviside,
                        LevSet1 = 1,
                        Jmp1 = JumpTypes.Heaviside
                    }, 
                    levelSet0, 
                    levelSet1),
                new LevelSetCombination( new CombinedID
                    {
                        LevSet0 = 0,
                        Jmp0 = JumpTypes.OneMinusHeaviside,
                        LevSet1 = 1,
                        Jmp1 = JumpTypes.Heaviside
                    },
                    levelSet0,
                    levelSet1),
                new LevelSetCombination( new CombinedID
                    {
                        LevSet0 = 0,
                        Jmp0 = JumpTypes.OneMinusHeaviside,
                        LevSet1 = 1,
                        Jmp1 = JumpTypes.OneMinusHeaviside
                    },
                    levelSet0,
                    levelSet1),
                new LevelSetCombination( new CombinedID
                    {
                        LevSet0 = 0,
                        Jmp0 = JumpTypes.Heaviside,
                        LevSet1 = 1,
                        Jmp1 = JumpTypes.OneMinusHeaviside
                    },
                    levelSet0,
                    levelSet1),
            };

        }
        
        int GetQuadOrder(int? quadorder)
        {
            int qO = 0;
            if (quadorder == null)
            {
                qO = Math.Max(((LevelSet)levelSets[0].LevelSet).Basis.Degree, ((LevelSet)levelSets[1].LevelSet).Basis.Degree);
            }
            else
            {
                qO = (int)quadorder;
            }
            return qO;
        }

        LevelSetCombination FindPhi(CombinedID iD)
        {
            foreach(LevelSetCombination phi in phis)
            {
                if (phi.Equals(iD))
                {
                    return phi;
                }
            }
            throw new Exception("Factory not found");
        }

        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1)
        {
            CombinedID id = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = jmp0,
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };
            LevelSetCombination lscomb = FindPhi(id);
            return new DoubleSayeQuadratureFactory(new DoubleSayeEdgeScheme(levelSets, lscomb), levelSets);

        }

        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0,
            int levSetIndex1,
            JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
        {
            CombinedID id0 = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = JumpTypes.Heaviside,
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };

            LevelSetCombination lscomb = FindPhi(id0);
            return new DoubleSayeQuadratureFactory(new DoubleSayeSurfaceScheme(levelSets, lscomb), levelSets);
        }

        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1)
        {
            CombinedID id = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = jmp0,
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };
            LevelSetCombination lscomb = FindPhi(id);

            return new DoubleSayeQuadratureFactory(new DoubleSayeVolumeScheme(levelSets, lscomb),levelSets);
        }

        public IQuadRuleFactory<QuadRule> GetEdgePointRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory) {

            CombinedID id = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = jmp1, // does not matter here, we need the id for the special edge detector
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };

            //var edgeScheme = new BruteForceEdgePointScheme(Phi, Phi2, SqrtGram);
            //return new BruteForceQuadratureFactory(backupFactory, levelSets, id, edgeScheme, 200) ;
            throw new NotImplementedException();
        }

        public IQuadRuleFactory<QuadRule> GetIntersectionFactory(int levSetIndex0, int levSetIndex1, IQuadRuleFactory<QuadRule> backupFactory) {

            //void Phi(int cell, NodeSet nodes, MultidimensionalArray result) {
            //    LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
            //    levelSet0.Evaluate(cell, 1, nodes, result);

            //    LevelSet levelSet1 = (LevelSet)levelSets[levSetIndex1].LevelSet;
            //    MultidimensionalArray result1 = MultidimensionalArray.Create(1, result.Length);
            //    levelSet1.Evaluate(cell, 1, nodes, result1);

            //    for(int i = 0; i < result.Length; ++i) {
            //        result[0, i] = Math.Max(Math.Abs(result[0, i]), Math.Abs(result1[0, i]));
            //    }
            //}

            //double Det(int cell) {
            //    MultidimensionalArray inverseJacobian = grid.JacobianDeterminat.GetValue_Cell(Square.Instance.Center, cell, 1);
            //    double g = 1 / inverseJacobian[0, 0];
            //    return g;
            //}

            //CombinedID id = new CombinedID {
            //    LevSet0 = levSetIndex0,
            //    Jmp0 = JumpTypes.Heaviside, // does not matter here, we need the id for the special edge detector
            //    LevSet1 = levSetIndex1,
            //    Jmp1 = JumpTypes.Heaviside
            //};

            //var surfaceScheme = new BruteForceZeroScheme(Phi, Det);
            throw new NotImplementedException();
            //return new BruteForceQuadratureFactory(backupFactory, levelSets, id, surfaceScheme, 400);
        }

    }

    /// <summary>
    /// Just a hack to enable special handling of cells where levelset=1 lies on an edge
    /// </summary>
    public static class SayeSettingsOverride {

        /// <summary>
        /// switch to enable special handling for cells, where the levelset=1 lies on an edge
        /// </summary>
        public static bool doubleCutCellOverride = false;
    }

    class DoubleSayeQuadratureFactory : IQuadRuleFactory<QuadRule>
    {
        DoubleSayeBaseScheme scheme;

        //MultiLevelSetOnEdgeDetector detector;

        LevelSetData[] data;

        public DoubleSayeQuadratureFactory(DoubleSayeBaseScheme scheme, LevelSetTracker.LevelSetData[] data)
        {
            this.scheme = scheme;
            this.data = data;
        }

        //public DoubleSayeQuadratureFactory(LevelSetTracker.LevelSetData[] data, CombinedID id, DoubleSayeBaseScheme scheme) : this(scheme,data)
        //{
        //    detector = new MultiLevelSetOnEdgeDetector(data, id);
        //}

        public RefElement RefElement => scheme.GetRefElement();

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public virtual IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            List<ChunkRulePair<QuadRule>> rule = new List<ChunkRulePair<QuadRule>>();
            //adds a quadrature rule for every element in the mask
            foreach (Chunk chunk in mask)
            {
                for (int i = chunk.i0; i < chunk.JE; ++i)
                {
                    Chunk singleChunk = Chunk.GetSingleElementChunk(i);
                    ChunkRulePair<QuadRule> pair = new ChunkRulePair<QuadRule>(singleChunk, scheme.GetQuadRule(i,order));
                    rule.Add(pair);
                }
            }
            return rule;
        }

        //unused
        public HyperRectangle CellToGlobalHR(int j, GridData gdat)
        {
            var eC1 = gdat.Edges.CellIndices[j,0];
            var eC2 = gdat.Edges.CellIndices[j,1];
            // get spatial dimension
            int D = gdat.SpatialDimension;
            // get the cell center
            var cC = gdat.iGeomCells.GetCenter(j);
            //create tensor for diameters
            var dC = Tensor1.Zeros(D);
            //get the indicies for the cell verticies
            var vC = gdat.iGeomCells.CellVertices[j];
            // will store min max of x,y,z 
            var mmV = new double[D, 2];
            // set default values to +,- infinity
            for (int d = 0; d < D; d++)
            {
                mmV[d, 0] = Double.MaxValue;
                mmV[d, 1] = Double.MinValue;
            }
            //loop over all vertex indices in order to obtain min/max 
            for (int i = 0; i < vC.Length; i++)
            {
                // get the coordinates of the vertex
                var coord = gdat.Vertices.Coordinates.ExtractSubArrayShallow(vC[i], -1);

                // assign if bigger/smaller than current minimum
                for (int d = 0; d < D; d++)
                {
                    if (mmV[d, 0] > coord[d])
                    {
                        mmV[d, 0] = coord[d];
                    }
                    if (mmV[d, 1] < coord[d])
                    {
                        mmV[d, 1] = coord[d];
                    }
                    // update the diameter accordingly
                    dC[d] = mmV[d, 1] - mmV[d, 0];
                }

            }

            var HR = new HyperRectangle(D);
            if (D == 2)
            {
                HR.Center = Tensor1.Vector(cC.x, cC.y);
            }
            else
            {
                HR.Center = Tensor1.Vector(cC.x, cC.y, cC.z);
            }
            HR.Diameters = dC;


            return HR;
        }


    }
    internal interface ISchemeWO
    {
        RefElement ReferenceElement { get; }

        QuadRule GetQuadRule(int j,int order);

        void Initialize(int resolution);
    }
    internal abstract class DoubleSayeBaseScheme: ISchemeWO
    {
        public LevelSetData[] data;
        
        public LevelSetCombination lscomb;
        public int D ;

        RefElement ISchemeWO.ReferenceElement => GetRefElement();

        public abstract RefElement GetRefElement();

        public DoubleSayeBaseScheme( LevelSetData[] data, LevelSetCombination lscomb)
        {
            this.data = data;
            
            this.lscomb = lscomb;
            this.D = ((LevelSet)data[0].LevelSet).Basis.GridDat.SpatialDimension;
        }
        
        public int GetNoOfQuadNodes(int quadorder)
        {
            return quadorder + 1;
        }
        public QuadRule GetQuadRule(int j,int order)
        {
            //Get Phi Evaluators
            (var phi0, var phi1) = GetPhiEval(j);

            //Get NoOfQuadNodes
            int noOfNodes = GetNoOfQuadNodes(order);

            // get a quadrule finder
            Quadrater finder = new Quadrater();

            //get the cell on which we integrate, here depending on the scheme the right integration domain is chosen
            HyperRectangle cell = GetCell();

            //get the symbols
            (var s1, var s2) = GetSymbols(lscomb);

            // Get The Quadrature rule From Intersectingquadrature
            QuadratureRule ruleQ = finder.FindRule(phi0, s1, phi1, s2, cell, noOfNodes);

            // Creates a QuadRule Object <- The One BoSSS uses
            QuadRule rule = QuadRule.CreateEmpty(GetRefElement(), ruleQ.Count, D, true);

            //loop over all quadrature Nodes
            for (int i = 0; i < ruleQ.Count; ++i)
            {
                QuadratureNode qNode = ruleQ[i];
                for (int d = 0; d < D; d++)
                {
                    rule.Nodes[i, d] = qNode.Point[d];
                }
                rule.Weights[i] = qNode.Weight;
            }
            rule.Nodes.LockForever();
            return rule;

        }
        /// <summary>
        /// This method reutrn a Hyperrectangle that is exactly the domain we want to integrate
        /// </summary>
        /// <returns></returns>
        public abstract HyperRectangle GetCell();

        /// <summary>
        /// This method reutrn a Hyperrectangle that is exactly the domain we want to integrate
        /// </summary>
        /// <returns></returns>
        public abstract (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j);

        //does nothing so far
        public void Initialize(int resolution)
        {
            
        }
        /// <summary>
        /// utility function to get the Symbols used by Intersecting Quadrature
        /// </summary>
        /// <param name="lscomb"></param>
        /// <returns></returns>
        public (Symbol s1, Symbol s2) GetSymbols(LevelSetCombination lscomb)
        {
            double sign0 = lscomb.sign0;
            double sign1 = lscomb.sign1;

            if (sign0 == -1.0 && sign1 == 1.0)
            {
                return (Symbol.Minus, Symbol.Plus);
            }
            else if (sign0 == -1.0 && sign1 == -1.0)
            {
                return (Symbol.Minus, Symbol.Minus);
            }
            else if (sign0 == 1.0 && sign1 == 1.0)
            {
                return (Symbol.Plus, Symbol.Plus);
            }
            else
            {
                return (Symbol.Plus, Symbol.Minus);
            }
        }    }
    internal class DoubleSayeVolumeScheme : DoubleSayeBaseScheme
    {
        public DoubleSayeVolumeScheme(LevelSetData[] data, LevelSetCombination lscomb) : base( data, lscomb)
        {
            
        }
        public override RefElement GetRefElement()
        {
            if (D == 1)
            {
                return Line.Instance;
            }else if (D == 2)
            {
                return Square.Instance;
            }
            else if (D == 3)
            {
                return Cube.Instance;
            }
            else
            {
                throw new ArgumentException();
            }
        }
        public override HyperRectangle GetCell()
        {
            return new UnitCube(D);
        }
        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            return (new lSEval((LevelSet)data[0].LevelSet,j), new lSEval((LevelSet)data[1].LevelSet,j));
        }
    }
    internal class DoubleSayeEdgeScheme : DoubleSayeBaseScheme
    {
        public DoubleSayeEdgeScheme( LevelSetData[] data, LevelSetCombination lscomb) : base(data, lscomb)
        {
        }
        public override RefElement GetRefElement()
        {
            if (D == 1)
            {
                return Grid.RefElements.Point.Instance;
            }
            else if (D == 2)
            {
                return Line.Instance;
            }
            else if (D == 3)
            {
                return Square.Instance;
            }
            else
            {
                throw new ArgumentException();
            }
        }
        public override HyperRectangle GetCell()
        {
            //define the Cell
            HyperRectangle cell = new UnitCube(D-1);
            return cell;
        }
        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            return (new lSEvalEdge((LevelSet)data[0].LevelSet, j), new lSEvalEdge((LevelSet)data[1].LevelSet, j));
        }

    }
    internal class DoubleSayeSurfaceScheme : DoubleSayeBaseScheme
    {
        public DoubleSayeSurfaceScheme(LevelSetData[] data, LevelSetCombination lscomb) : base(data,  lscomb)
        {
        }
        public override RefElement GetRefElement()
        {
            if (D == 1)
            {
                return Line.Instance;
            }
            else if (D == 2)
            {
                return Square.Instance;
            }
            else if (D == 3)
            {
                return Cube.Instance;
            }
            else
            {
                throw new ArgumentException();
            }
        }
        public override HyperRectangle GetCell()
        {
            //define the Cell
            HyperRectangle cell = new UnitCube(D);
            cell.Dimension = D - 1;
            return cell;
        }
        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            return (new lSEval((LevelSet)data[0].LevelSet, j), new lSEval((LevelSet)data[1].LevelSet, j));
        }
    }

    internal class lSEvalEdge : IScalarFunction
    {
        //levelSet which is evaluated by this class
        public LevelSet m_levelSet;
        //cell the Level Set is evaluated on
        public int j;
        int IScalarFunction.M => m_levelSet.Basis.GridDat.SpatialDimension-1;

        public lSEvalEdge(LevelSet levelSet, int j)
        {
            m_levelSet = levelSet;
            this.j = j;
        }
        double IScalarFunction.Evaluate(Tensor1 x)
        {
            (NodeSet NS, int i0) = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //var inp = MultidimensionalArray.Create(1, m_levelSet.Basis.GridDat.SpatialDimension - 1);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            var inp = MultidimensionalArray.Create(1, 1);
            m_levelSet.EvaluateEdge(i0, 1, NS,inp, ev);
            return inp[0];
        }

        (double evaluation, Tensor1 gradient) IScalarFunction.EvaluateAndGradient(Tensor1 x)
        {
            (NodeSet NS, int i0) = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var inp = MultidimensionalArray.Create(1, 1);
            var ev = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
            MultidimensionalArray gradIn = MultidimensionalArray.Create(1, 1, x.M);
            m_levelSet.EvaluateEdge(i0, 1, NS, inp, ev,GradientIN: grad ,GradientOT:gradIn);
            return (inp[0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)));
        }

        (double evaluation, Tensor1 gradient, Tensor2 hessian) IScalarFunction.EvaluateAndGradientAndHessian(Tensor1 x)
        {
            (NodeSet NS, int i0) = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var inp = MultidimensionalArray.Create(1, 1);
            var ev = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
            MultidimensionalArray gradIn = MultidimensionalArray.Create(1, 1, x.M);
            m_levelSet.EvaluateEdge(i0, 1, NS, inp, ev, GradientIN: grad, GradientOT: gradIn);

            //no Hessian Evaluation on Edges
            MultidimensionalArray hess = MultidimensionalArray.Create(1, 1, x.M, x.M);
            //m_levelSet.EvaluateHessian(i0, 1, NS, hess);

            return (inp[0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)), lSEvalUtil.ToTensor2(hess.ExtractSubArrayShallow(0, 0, -1, -1)));

        }
    }
    internal class lSEval : IScalarFunction
    {
        //levelSet which is evaluated by this class
        public LevelSet m_levelSet;
        //cell the Level Set is evaluated on
        public int j;
        int IScalarFunction.M => m_levelSet.Basis.GridDat.SpatialDimension;

        public lSEval(LevelSet levelSet, int j)
        {
            m_levelSet = levelSet;
            this.j = j;
        }

        double IScalarFunction.Evaluate(Tensor1 x)
        {
            (NodeSet NS, int i0) = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            m_levelSet.Evaluate(i0, 1, NS, ev);
            return ev[0];
        }

        (double evaluation, Tensor1 gradient) IScalarFunction.EvaluateAndGradient(Tensor1 x)
        {
            (NodeSet NS, int i0) = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            m_levelSet.Evaluate(i0, 1, NS, ev);

            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
            m_levelSet.EvaluateGradient(i0, 1, NS, grad);

            return (ev[0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)));

        }

        (double evaluation, Tensor1 gradient, Tensor2 hessian) IScalarFunction.EvaluateAndGradientAndHessian(Tensor1 x)
        {
            (NodeSet NS, int i0) = lSEvalUtil.NSFromTensor(x,m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            m_levelSet.Evaluate(i0, 1, NS, ev);


            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
            MultidimensionalArray hess = MultidimensionalArray.Create(1, 1, x.M, x.M);
            m_levelSet.EvaluateGradient(i0, 1, NS, grad);
            m_levelSet.EvaluateHessian(i0, 1, NS, hess);

            return (ev[0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)), lSEvalUtil.ToTensor2(hess.ExtractSubArrayShallow(0, 0, -1, -1)));

        }
        //(double evaluation, Tensor1 gradient) EvaluateAndGradientGlobal(Tensor1 x)
        //{
        //    (NodeSet NS, int i0) = NSFromGlobalTensor(x);
        //    //Evaluates the LevelSet using the NodeSet and the Cell Number
        //    var ev = MultidimensionalArray.Create(1, 1);
        //    m_levelSet.Evaluate(i0, 1, NS, ev);

        //    MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
        //    m_levelSet.EvaluateGradient(i0, 1, NS, grad);

        //    return (ev[0], ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)));

        //}
        //double EvaluateGlobal(Tensor1 x)
        //{
        //    (NodeSet NS, int i0) = NSFromGlobalTensor(x);
        //    //Evaluates the LevelSet using the NodeSet and the Cell Number
        //    var ev = MultidimensionalArray.Create(1, 1);
        //    m_levelSet.Evaluate(i0, 1, NS, ev);
        //    return ev[0];

        //}
    }
    internal static class lSEvalUtil {
        
        public static Vector TensorToVector(Tensor1 x)
        {
            Vector point;
            if (x.M == 2)
            {
                point = new Vector(x[0], x[1], x[2]);
            }
            else
            {
                point = new Vector(x[0], x[1], x[2]);
            }
            return point;
        }
        //public Tensor1 ToTensor(IEnumerable x)
        //{
        //    Tensor1 point;
        //    if (x.Count() == 1)
        //    {
        //        point = Tensor1.Vector(x[0]);
        //    }
        //    else if (x.Count() == 2)
        //    {
        //        point = Tensor1.Vector(x[0], x[1]);
        //    }
        //    else
        //    {
        //        point = Tensor1.Vector(x[0], x[1], x[2]);
        //    }
        //    return point;
        //}
        public static Tensor1 ToTensor1(MultidimensionalArray x)
        {
            Tensor1 vec;
            if (x.Length == 1)
            {
                vec = Tensor1.Vector(x[0]);
            }
            else if (x.Length == 2)
            {
                vec = Tensor1.Vector(x[0], x[1]);
            }
            else
            {
                vec = Tensor1.Vector(x[0], x[1], x[2]);
            }
            return vec;
        }
        public static Tensor2 ToTensor2(MultidimensionalArray x)
        {
            Tensor2 mat = Tensor2.Zeros(x.Lengths[0], x.Lengths[0]);
            for (int i = 0; i < x.Lengths[0]; i++)
            {
                for (int j = 0; j < x.Lengths[1]; j++)
                {
                    mat[i, j] = x[i, j];
                }
            }

            return mat;
        }
        public static MultidimensionalArray ToMultArray(Tensor1 x)
        {
            MultidimensionalArray vec = MultidimensionalArray.Create(x.M);
            for (int i = 0; i < x.M; i++)
            {
                vec[i] = x[i];
            }
            return vec;
        }
        public static MultidimensionalArray ToMultArray(Tensor2 x)
        {
            MultidimensionalArray vec = MultidimensionalArray.Create(x.M, x.N);
            for (int i = 0; i < x.M; i++)
            {
                for (int j = 0; j < x.N; j++)
                {
                    vec[i, j] = x[i, j];
                }
            }
            return vec;
        }

        public static (NodeSet NS, int i0) NSFromTensor(Tensor1 x, LevelSet m_levelSet)
        {
            //transform into multarray
            MultidimensionalArray GlobIn = MultidimensionalArray.Create(1, x.M);
            for (int i = 0; i < x.M; i++)
            {
                GlobIn[0, i] = x[i];
            }

            //output array for local coordinates
            //var LocOut = MultidimensionalArray.Create(1, 1, GlobIn.Lengths[1]);
            //gives the cell the point is located
            m_levelSet.GridDat.LocatePoint(GlobIn.ExtractSubArrayShallow(0, -1).To1DArray(), out long id, out long i0, out bool IsInside, out bool Onm_levelSetProces);
            //gives the local Coordinates
            //m_levelSet.GridDat.TransformGlobal2Local(GlobIn, LocOut, (int)i0, 1, 0);

            //gives the NodeSet using the local Coords
            var NS = new NodeSet(m_levelSet.GridDat.iGeomCells.RefElements[0], GlobIn, true);

            return (NS, (int)i0);
        }
        public static (NodeSet NS, int i0) NSFromGlobalTensor(Tensor1 x, LevelSet m_levelSet)
        {
            //transform into multarray
            MultidimensionalArray GlobIn = MultidimensionalArray.Create(1, x.M);
            for (int i = 0; i < x.M; i++)
            {
                GlobIn[0, i] = x[i];
            }

            //output array for local coordinates
            var LocOut = MultidimensionalArray.Create(1, 1, GlobIn.Lengths[1]);
            //gives the cell the point is located
            m_levelSet.GridDat.LocatePoint(GlobIn.ExtractSubArrayShallow(0, -1).To1DArray(), out long id, out long i0, out bool IsInside, out bool Onm_levelSetProces);
            //gives the local Coordinates
            m_levelSet.GridDat.TransformGlobal2Local(GlobIn, LocOut, (int)i0, 1, 0);

            //gives the NodeSet using the local Coords
            var NS = new NodeSet(m_levelSet.GridDat.iGeomCells.RefElements[0], LocOut.ExtractSubArrayShallow(0, -1, -1), true);

            return (NS, (int)i0);
        }
    }
}
