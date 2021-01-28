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

namespace BoSSS.Foundation.XDG.Quadrature
{
    class LevelSetCombination
    {
        CombinedID ID;

        LevelSet levelSet0;

        LevelSet levelSet1;

        double sign0;

        double sign1;

        public LevelSetCombination(CombinedID iD, LevelSet levelSet0, LevelSet levelSet1)
        {
            this.ID = iD;
            sign0 = ToDouble(iD.Jmp0);
            sign1 = ToDouble(iD.Jmp1);
            this.levelSet0 = levelSet0;
            this.levelSet1 = levelSet1;
        }

        public void Evaluate(int j0, NodeSet x, MultidimensionalArray result)
        {
            levelSet0.Evaluate(j0, 1, x, result);
            MultidimensionalArray result1 = MultidimensionalArray.Create(1, result.Length);
            levelSet1.Evaluate(j0, 1, x, result1);

            for(int i = 0; i < result.Length; ++i)
            {
                result[0,i] = Math.Min(sign0 * result[0,i], sign1 * result1[0,i]);
            }
        }

        public void EvaluateEdge(int j0, NodeSet x, MultidimensionalArray resultIn, MultidimensionalArray resultOut)
        {
            levelSet0.EvaluateEdge(j0, 1, x, resultIn, resultOut);
            MultidimensionalArray result1 = MultidimensionalArray.Create(1, resultIn.Length);
            MultidimensionalArray result2 = MultidimensionalArray.Create(1, resultOut.Length);
            levelSet1.EvaluateEdge(j0, 1, x, result1, result2);

            for (int i = 0; i < resultIn.Length; ++i)
            {
                resultIn[0, i] = Math.Min(sign0 * resultIn[0, i], sign1 * result1[0, i]);
                resultOut[0, i] = Math.Min(sign0 * resultOut[0, i], sign1 * result2[0, i]);
            }
        }

        double ToDouble(JumpTypes type)
        {
            if (type == JumpTypes.Heaviside)
            {
                return 1.0;
            }
            else if (type == JumpTypes.OneMinusHeaviside)
            {
                return -1.0;
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        public bool Equals(CombinedID otherID)
        {
            return ID.Equals(otherID);
        }
    }

    class MultiLevelSetBruteForceQuadratureFactory
    {
        LevelSetCombination[] phis;

        GridData grid;

        LevelSetData[] levelSets;

        public MultiLevelSetBruteForceQuadratureFactory(LevelSetData[] levelSets)
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

        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1)
        {
            CombinedID id = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = jmp0,
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };
            LevelSetCombination phi = FindPhi(id);
            var edgeScheme = new BruteForceEdgeScheme(phi.EvaluateEdge);

            return new BruteForceQuadratureFactory(edgeScheme);
        }

        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0,
            int levSetIndex1,
            JumpTypes jmp1)
        {
            CombinedID id0 = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = JumpTypes.Heaviside,
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };
            
            LevelSetCombination phi0 = FindPhi(id0);
            
            Vector Scales(int cell)
            {
                MultidimensionalArray jacobian = grid.Jacobian.GetValue_Cell(Square.Instance.Center, cell, 1).ExtractSubArrayShallow(0, 0, -1, -1);
                return new Vector(1 / jacobian[0, 0], 1 / jacobian[1, 1]);
            }

            void Gradient(int cell, NodeSet nodes, MultidimensionalArray result)
            {
                LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
                levelSet0.EvaluateGradient(cell, 1, nodes, result);
            }

            void Phi0(int cell, NodeSet nodes, MultidimensionalArray result)
            {
                LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
                levelSet0.Evaluate(cell, 1, nodes, result);
            }

            var surfaceScheme = new BruteForceSurfaceScheme(phi0.Evaluate, Phi0, Scales, Gradient);
            return new BruteForceQuadratureFactory(surfaceScheme);
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
            LevelSetCombination phi = FindPhi(id);
            var volumeScheme = new BruteForceVolumeScheme(phi.Evaluate);
            return new BruteForceQuadratureFactory(volumeScheme);
        }
    }

    class BruteForceQuadratureFactory : IQuadRuleFactory<QuadRule>
    {
        IScheme scheme;

        public BruteForceQuadratureFactory(IScheme scheme)
        {
            this.scheme = scheme;
        }

        public RefElement RefElement => scheme.ReferenceElement;

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public virtual IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            int resolution = 99; // (int) Math.Pow(2, order);
            scheme.Initialize(resolution);
            List<ChunkRulePair<QuadRule>> rule = new List<ChunkRulePair<QuadRule>>();
            foreach(Chunk chunk in mask)
            {
                for(int i = chunk.i0; i < chunk.JE; ++i)
                {
                    Chunk singleChunk = Chunk.GetSingleElementChunk(i);
                    ChunkRulePair<QuadRule> pair = new ChunkRulePair<QuadRule>(singleChunk, scheme.GetQuadRule(i));
                    rule.Add(pair);
                }
            }
            return rule;
        }
    }
}
