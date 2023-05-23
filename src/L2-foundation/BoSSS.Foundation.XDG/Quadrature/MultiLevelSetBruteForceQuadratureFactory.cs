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
        public CombinedID ID;

        LevelSet levelSet0;

        LevelSet levelSet1;

        public double sign0 { get; set; }

        public double sign1 { get; set; }

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

        public static double ToDouble(JumpTypes type)
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

        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
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

            return new BruteForceQuadratureFactory(backupFactory, this.levelSets, id, edgeScheme, 400);
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
            
            LevelSetCombination phi0 = FindPhi(id0);
            
            Vector Scales(int cell)
            {
                MultidimensionalArray inverseJacobian = grid.InverseJacobian.GetValue_Cell(Square.Instance.Center, cell, 1).ExtractSubArrayShallow(0, 0, -1, -1);
                return new Vector(inverseJacobian[0, 0], inverseJacobian[1, 1]);
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
            return new BruteForceQuadratureFactory(backupFactory, levelSets, id0, surfaceScheme, 400);
        }

        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
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

            return new BruteForceQuadratureFactory(backupFactory, this.levelSets, id, volumeScheme, 200);
        }

        public IQuadRuleFactory<QuadRule> GetEdgePointRuleFactory(int levSetIndex0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory) {

            void Phi(int j0, NodeSet x, MultidimensionalArray resultIn, MultidimensionalArray resultOut) {
                LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
                levelSet0.EvaluateEdge(j0, 1, x, resultIn, resultOut);
            }

            void Phi2(int j0, NodeSet x, MultidimensionalArray resultIn, MultidimensionalArray resultOut) {
                LevelSet levelSet1 = (LevelSet)levelSets[levSetIndex1].LevelSet;
                levelSet1.EvaluateEdge(j0, 1, x, resultIn, resultOut);

                for (int i = 0; i < resultIn.Length; ++i) {
                    resultIn[0, i] = LevelSetCombination.ToDouble(jmp1) * resultIn[0, i];
                    resultOut[0, i] = LevelSetCombination.ToDouble(jmp1) * resultOut[0, i];
                }
            }

            double SqrtGram(int edge) {
                var g = grid.Edges.SqrtGramian[edge];
                g = 1 / g;
                return g;
            }

            CombinedID id = new CombinedID {
                LevSet0 = levSetIndex0,
                Jmp0 = jmp1, // does not matter here, we need the id for the special edge detector
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };

            var edgeScheme = new BruteForceEdgePointScheme(Phi, Phi2, SqrtGram);
            return new BruteForceQuadratureFactory(backupFactory, levelSets, id, edgeScheme, 200) ;
        }

        public IQuadRuleFactory<QuadRule> GetIntersectionFactory(int levSetIndex0, int levSetIndex1, IQuadRuleFactory<QuadRule> backupFactory) {

            void Phi(int cell, NodeSet nodes, MultidimensionalArray result) {
                LevelSet levelSet0 = (LevelSet)levelSets[levSetIndex0].LevelSet;
                levelSet0.Evaluate(cell, 1, nodes, result);

                LevelSet levelSet1 = (LevelSet)levelSets[levSetIndex1].LevelSet;
                MultidimensionalArray result1 = MultidimensionalArray.Create(1, result.Length);
                levelSet1.Evaluate(cell, 1, nodes, result1);
                
                for(int i = 0; i < result.Length; ++i) {
                    result[0, i] = Math.Max(Math.Abs(result[0, i]), Math.Abs(result1[0, i]));
                }
            }

            double Det(int cell) {
                MultidimensionalArray inverseJacobian = grid.JacobianDeterminat.GetValue_Cell(Square.Instance.Center, cell, 1);
                double g = 1 / inverseJacobian[0, 0];
                return g;
            }

            CombinedID id = new CombinedID {
                LevSet0 = levSetIndex0,
                Jmp0 = JumpTypes.Heaviside, // does not matter here, we need the id for the special edge detector
                LevSet1 = levSetIndex1,
                Jmp1 = JumpTypes.Heaviside
            };

            var surfaceScheme = new BruteForceZeroScheme(Phi, Det);
            return new BruteForceQuadratureFactory(backupFactory, levelSets, id, surfaceScheme, 400);
        }

    }

    /// <summary>
    /// Just a hack to enable special handling of cells where levelset=1 lies on an edge
    /// </summary>
    public static class BruteForceSettingsOverride {

        /// <summary>
        /// switch to enable special handling for cells, where the levelset=1 lies on an edge
        /// </summary>
        public static bool doubleCutCellOverride = false;
    }

    class BruteForceQuadratureFactory : IQuadRuleFactory<QuadRule>
    {
        IScheme scheme;

        IQuadRuleFactory<QuadRule> fallBackFactory;
        MultiLevelSetOnEdgeDetector detector;

        int resolution; 

        public BruteForceQuadratureFactory(IScheme scheme, int resolution = 200)
        {
            this.scheme = scheme;
            this.resolution = resolution;
        }

        public BruteForceQuadratureFactory(IQuadRuleFactory<QuadRule> fallBackFactory, LevelSetTracker.LevelSetData[] data, CombinedID id, IScheme scheme, int resolution = 200) : this(scheme, resolution) {
            this.fallBackFactory = fallBackFactory;
            detector = new MultiLevelSetOnEdgeDetector(data, id);
        }

        public RefElement RefElement => scheme.ReferenceElement;

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public virtual IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            scheme.Initialize(resolution);

            List<ChunkRulePair<QuadRule>> rule = new List<ChunkRulePair<QuadRule>>();
            foreach(Chunk chunk in mask)
            {
                for(int i = chunk.i0; i < chunk.JE; ++i)
                {
                    Chunk singleChunk = Chunk.GetSingleElementChunk(i);
                    if (BruteForceSettingsOverride.doubleCutCellOverride) {
                        if (detector != null) {
                            ExecutionMask singleMask;
                            bool special;
                            bool active;
                            if (mask is CellMask cm) {
                                special = detector.IsSpecialCell(i);
                                active = detector.IsActiveCell(i);
                                singleMask = new CellMask(cm.GridData, singleChunk, MaskType.Geometrical);
                            } else {
                                EdgeMask em = (EdgeMask)mask;
                                special = detector.IsSpecialEdge(i);
                                active = detector.IsActiveEdge(i);
                                singleMask = new EdgeMask(em.GridData, singleChunk, MaskType.Geometrical);
                            }

                            // detected a "special" cell, treat differently
                            if (special) {
                                QuadRule specialRule;
                                specialRule = GetSpecialQuadRule(i, singleMask, order, active);

                                // add to List and jump to next chunk
                                ChunkRulePair<QuadRule> pairSpecial = new ChunkRulePair<QuadRule>(singleChunk, specialRule);
                                rule.Add(pairSpecial);
                                continue;
                            }
                        }
                    }
                    ChunkRulePair<QuadRule> pair = new ChunkRulePair<QuadRule>(singleChunk, scheme.GetQuadRule(i));
                    rule.Add(pair);
                }
            }
            return rule;
        }

        private QuadRule GetSpecialQuadRule(int i, ExecutionMask singleMask, int order, bool activeChunk) {
            QuadRule specialRule;

            switch (scheme) {
                case BruteForceZeroScheme b:
                    // to avoid deadlocks? construct here
                    specialRule = ConstructSpecialQuadRule(i, singleMask, order);

                    // skip cells not in active phase, could otherwise lead to problems, when there are kinks in the levelset (i.e. different contact angle from right and left)
                    if (activeChunk) {
                        // do nothing, already built
                    } else {
                        specialRule = QuadRule.CreateEmpty(b.ReferenceElement, 1, b.ReferenceElement.SpatialDimension);
                        specialRule.Nodes.LockForever();
                    }
                    break;
                case BruteForceEdgePointScheme b:
                    specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;
                    break;
                case BruteForceEdgeScheme b: 
                    {
                        specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;

                        // remove edge scheme from "true" coinciding edge
                        int jCell = ((GridData)singleMask.GridData).Edges.CellIndices[i, 0];
                        int specialFace = detector.GetSpecialFace(jCell);
                        int iFace = ((GridData)singleMask.GridData).Edges.FaceIndices[i, 0];
                        bool coinciding = specialFace == iFace;
                        if (coinciding) {
                            specialRule = QuadRule.CreateEmpty(specialRule.RefElement, 1, specialRule.RefElement.SpatialDimension);
                            specialRule.Nodes.LockForever();
                        }
                        break;
                    }
                case BruteForceSurfaceScheme b:
                    if (fallBackFactory.RefElement == scheme.ReferenceElement) {
                        // skip cells not in active phase
                        if (activeChunk) {
                            specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;
                        } else {
                            specialRule = QuadRule.CreateEmpty(b.ReferenceElement, 1, b.ReferenceElement.SpatialDimension);
                            specialRule.Nodes.LockForever();
                        }
                    } else {
                        int iFace = detector.GetSpecialFace(i);
                        int iEdge = ((GridData)singleMask.GridData).Cells.GetEdgesForFace(i, iFace, out int InOrOut, out int[] FurtherEdges);

                        int J = ((GridData)singleMask.GridData).Cells.NoOfLocalUpdatedCells;
                        int OtherCell = ((GridData)singleMask.GridData).Edges.CellIndices[iEdge, InOrOut == 0 ? 1 : 0];
                        long ThisCellGlob = i + ((GridData)singleMask.GridData).CellPartitioning.i0;
                        long OtherCellGlob = OtherCell < J ? OtherCell + ((GridData)singleMask.GridData).CellPartitioning.i0 : ((GridData)singleMask.GridData).Parallel.GlobalIndicesExternalCells[OtherCell - J];
                        bool ThisConform = InOrOut == 0 ? ((GridData)singleMask.GridData).Edges.IsEdgeConformalWithCell1(iEdge) : ((GridData)singleMask.GridData).Edges.IsEdgeConformalWithCell2(iEdge);
                        bool OtherConform = InOrOut == 0 ? ((GridData)singleMask.GridData).Edges.IsEdgeConformalWithCell2(iEdge) : ((GridData)singleMask.GridData).Edges.IsEdgeConformalWithCell1(iEdge);

                        // to avoid deadlocks? construct here
                        specialRule = ConstructSpecialQuadRule(i, singleMask, order);

                        // take conformal cell or that one with globally lower index
                        if (ThisConform & OtherConform) {
                            if(ThisCellGlob < OtherCellGlob) {
                                // do nothing, already built
                            } else {
                                specialRule = QuadRule.CreateEmpty(b.ReferenceElement, 1, b.ReferenceElement.SpatialDimension);
                                specialRule.Nodes.LockForever();
                            }
                        } else if (ThisConform & !OtherConform) {
                            // do nothing, already built
                        } else if (!ThisConform & OtherConform) {
                            specialRule = QuadRule.CreateEmpty(b.ReferenceElement, 1, b.ReferenceElement.SpatialDimension);
                            specialRule.Nodes.LockForever();
                        } else {
                            throw new NotSupportedException(String.Format("Error in cell {0}, {1}: Only one cell should have the hanging node.", i, OtherCell));
                        }                        
                    }
                    break;
                case BruteForceVolumeScheme b:
                    // skip cells not in active phase
                    if (activeChunk) {
                        specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;
                    } else {
                        specialRule = QuadRule.CreateEmpty(b.ReferenceElement, 1, b.ReferenceElement.SpatialDimension);
                        specialRule.Nodes.LockForever();
                    }
                    break;
                default:
                    throw new NotImplementedException();
            }            
            return specialRule;
        }

        private QuadRule ConstructSpecialQuadRule(int i, ExecutionMask singleMask, int order) {
            QuadRule specialRule;
            if (singleMask.GetType() != typeof(CellMask))
                throw new NotSupportedException();

            // determine face
            int iFace = detector.GetSpecialFace(i);
            int iEdge = ((GridData)singleMask.GridData).Cells.GetEdgesForFace(i, iFace, out int InOrOut, out int[] FurtherEdges);
            // select all edges, relevant in cells with hanging nodes
            singleMask = new EdgeMask(singleMask.GridData, Chunk.GetSingleElementChunk(iEdge), MaskType.Geometrical);
            if (FurtherEdges != null) {
                specialRule = QuadRule.CreateEmpty(scheme.ReferenceElement, 1, scheme.ReferenceElement.SpatialDimension);
                specialRule.Nodes.LockForever();

                foreach (var edg in FurtherEdges) {
                    singleMask = singleMask.Union(new EdgeMask(singleMask.GridData, Chunk.GetSingleElementChunk(edg), MaskType.Geometrical));
                }
            }
            // construct edgerule
            var specialRule_t = fallBackFactory.GetQuadRuleSet(singleMask, order);
            List<NodeSet> specialNodes = new List<NodeSet>();
            List<MultidimensionalArray> specialWeights = new List<MultidimensionalArray>();
            int NoOfNodes = 0;
            foreach (var crp in specialRule_t) {
                var rule_t = crp.Rule;
                foreach (int edg in crp.Chunk.Elements) {
                    int trf;
                    if (singleMask.GridData.iGeomEdges.CellIndices[edg, 0] == i) {
                        trf = singleMask.GridData.iGeomEdges.Edge2CellTrafoIndex[edg, 0];
                    } else {
                        trf = singleMask.GridData.iGeomEdges.Edge2CellTrafoIndex[edg, 1];
                    }
                    // transform edgerule to volume rule
                    var nodes_t = rule_t.Nodes.GetVolumeNodeSet(singleMask.GridData, trf, false);
                    var weights_t = rule_t.Weights;

                    // scale accordingly!, for a volume rule generated through an edge rule, this is length of linerefelem / length of edge perpendicular to rule edge
                    double scale = singleMask.GridData.iGeomCells.GetRefElement(i).Volume / singleMask.GridData.iGeomCells.GetCellVolume(i) * singleMask.GridData.iGeomEdges.SqrtGramian[edg];
                    weights_t.Scale(scale);

                    specialNodes.Add(nodes_t);
                    specialWeights.Add(weights_t);
                    NoOfNodes += nodes_t.NoOfNodes;
                }
            }
            MultidimensionalArray specialNodes_t = MultidimensionalArray.Create(NoOfNodes, singleMask.GridData.SpatialDimension);
            MultidimensionalArray specialWeights_t = MultidimensionalArray.Create(NoOfNodes);
            int offset = 0;
            for (int j = 0; j < specialNodes.Count; j++) {
                int count = specialNodes[j].Lengths[0];
                specialNodes_t.ExtractSubArrayShallow(new int[] { offset, 0 }, new int[] { offset + count - 1, singleMask.GridData.SpatialDimension - 1 }).Acc(1.0, specialNodes[j]);
                specialWeights_t.ExtractSubArrayShallow(new int[] { offset }, new int[] { offset + count - 1 }).Acc(1.0, specialWeights[j]);
                offset += count;
            }
            specialRule = new QuadRule();
            specialRule.Nodes = new NodeSet(scheme.ReferenceElement, specialNodes_t, true);
            specialRule.Weights = specialWeights_t;

            return specialRule;
        }
    }
}
