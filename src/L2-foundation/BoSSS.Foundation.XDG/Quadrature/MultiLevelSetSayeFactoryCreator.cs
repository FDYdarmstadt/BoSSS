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

        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
        {
            //CombinedID id = new CombinedID
            //{
            //    LevSet0 = levSetIndex0,
            //    Jmp0 = jmp0,
            //    LevSet1 = levSetIndex1,
            //    Jmp1 = jmp1
            //};
            //LevelSetCombination phi = FindPhi(id);
            //var edgeScheme = new BruteForceEdgeScheme(phi.EvaluateEdge);           

            //return new BruteForceQuadratureFactory(backupFactory, this.levelSets, id, edgeScheme, 400);
            throw new NotImplementedException();    
        }

        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0,
            int levSetIndex1,
            JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
        {
            //CombinedID id0 = new CombinedID
            //{
            //    LevSet0 = levSetIndex0,
            //    Jmp0 = JumpTypes.Heaviside,
            //    LevSet1 = levSetIndex1,
            //    Jmp1 = jmp1
            //};
            
            //LevelSetCombination phi0 = FindPhi(id0);
            
            //Vector Scales(int cell)
            //{
            //    MultidimensionalArray inverseJacobian = grid.InverseJacobian.GetValue_Cell(Square.Instance.Center, cell, 1).ExtractSubArrayShallow(0, 0, -1, -1);
            //    return new Vector(inverseJacobian[0, 0], inverseJacobian[1, 1]);
            //}

            //void Gradient(int cell, NodeSet nodes, MultidimensionalArray result)
            //{
            //    LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
            //    levelSet0.EvaluateGradient(cell, 1, nodes, result);
            //}

            //void Phi0(int cell, NodeSet nodes, MultidimensionalArray result)
            //{
            //    LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
            //    levelSet0.Evaluate(cell, 1, nodes, result);
            //}

            //var surfaceScheme = new BruteForceSurfaceScheme(phi0.Evaluate, Phi0, Scales, Gradient);
            //return new BruteForceQuadratureFactory(backupFactory, levelSets, id0, surfaceScheme, 400);
            throw new NotImplementedException();
        }

        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
        {
            //CombinedID id = new CombinedID
            //{
            //    LevSet0 = levSetIndex0,
            //    Jmp0 = jmp0,
            //    LevSet1 = levSetIndex1,
            //    Jmp1 = jmp1
            //};
            //LevelSetCombination phi = FindPhi(id);
            //var volumeScheme = new BruteForceVolumeScheme(phi.Evaluate);

            //return new BruteForceQuadratureFactory(backupFactory, this.levelSets, id, volumeScheme, 200);
            throw new NotImplementedException();
        }

        public IQuadRuleFactory<QuadRule> GetEdgePointRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory) {

            //void Phi(int j0, NodeSet x, MultidimensionalArray resultIn, MultidimensionalArray resultOut) {
            //    LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
            //    levelSet0.EvaluateEdge(j0, 1, x, resultIn, resultOut);
            //}

            //void Phi2(int j0, NodeSet x, MultidimensionalArray resultIn, MultidimensionalArray resultOut) {
            //    LevelSet levelSet1 = (LevelSet)levelSets[levSetIndex1].LevelSet;
            //    levelSet1.EvaluateEdge(j0, 1, x, resultIn, resultOut);

            //    for (int i = 0; i < resultIn.Length; ++i) {
            //        resultIn[0, i] = LevelSetCombination.ToDouble(jmp1) * resultIn[0, i];
            //        resultOut[0, i] = LevelSetCombination.ToDouble(jmp1) * resultOut[0, i];
            //    }
            //}

            //double SqrtGram(int edge) {
            //    var g = grid.Edges.SqrtGramian[edge];
            //    g = 1 / g;
            //    return g;
            //}

            //CombinedID id = new CombinedID {
            //    LevSet0 = levSetIndex0,
            //    Jmp0 = jmp1, // does not matter here, we need the id for the special edge detector
            //    LevSet1 = levSetIndex1,
            //    Jmp1 = jmp1
            //};

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
        IScheme scheme;

        MultiLevelSetOnEdgeDetector detector;

        int quadorder;

        public DoubleSayeQuadratureFactory(IScheme scheme, int quadorder)
        {
            this.scheme = scheme;
            this.quadorder = quadorder;
        }

        public DoubleSayeQuadratureFactory(LevelSetTracker.LevelSetData[] data, CombinedID id, IScheme scheme, int quadorder) : this(scheme, quadorder)
        {
            detector = new MultiLevelSetOnEdgeDetector(data, id);
        }

        public RefElement RefElement => scheme.ReferenceElement;

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public virtual IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {

            //Copied from BruteForceScheme
            scheme.Initialize(quadorder);

            List<ChunkRulePair<QuadRule>> rule = new List<ChunkRulePair<QuadRule>>();
            foreach (Chunk chunk in mask)
            {
                for (int i = chunk.i0; i < chunk.JE; ++i)
                {
                    Chunk singleChunk = Chunk.GetSingleElementChunk(i);
                    ChunkRulePair<QuadRule> pair = new ChunkRulePair<QuadRule>(singleChunk, this.GetQuadRule(i));
                    rule.Add(pair);
                }
            }
            return rule;
        }
        public int GetNoOfQuadNodes(int quadorder)
        {
            return quadorder + 1;
        }
        public HyperRectangle CellToHR(int j, GridData gdat)
        {
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

        public (Symbol s1, Symbol s2) GetSymbols()
        {
            int sign0 = -1;
            int sign1 = 1;
            if (sign0 == -1 && sign1 == 1)
            {
                return (Symbol.Minus, Symbol.Plus);
            }
            else if (sign0 == -1 && sign1 == -1)
            {
                return (Symbol.Minus, Symbol.Minus);
            }
            else if (sign0 == 1 && sign1 == 1)
            {
                return (Symbol.Plus, Symbol.Plus);
            }
            else
            {
                return (Symbol.Plus, Symbol.Minus);
            }
        }

        public QuadRule GetQuadRule(int j, LevelSet phi0, LevelSet phi1, LevelSetCombination lscomb)
        {

            //Get NoOfQuadNodes
            int noOfNodes = GetNoOfQuadNodes(quadorder);

            // get a quadrule finder
            Quadrater finder = new Quadrater();

            //define the Cell
            HyperRectangle cell = CellToHR(j, (GridData)phi0.GridDat);

            //get the symbols
            (var s1, var s2) = GetSymbols();

            // Get The Quadrature rule From Intersectingquadrature
            QuadratureRule ruleQ = finder.FindRule(phi0, s1, phi1, s2, cell, noOfNodes);

            // Creates a QuadRule Object <- The One BoSSS uses
            QuadRule rule = QuadRule.CreateEmpty(Square.Instance, ruleQ.Count, phi0.Basis.GridDat.SpatialDimension, true);

            //loop over all quadrature Nodes
            for (int i = 0; i < ruleQ.Count; ++i)
            {
                QuadratureNode qNode = ruleQ[i];
                for (int d = 0; d < phi0.Basis.GridDat.SpatialDimension; d++)
                {
                    rule.Nodes[i, d] = qNode.Point[d];
                }
                rule.Weights[i] = qNode.Weight;
            }
            rule.Nodes.LockForever();
            return rule;
        }


    }
}
