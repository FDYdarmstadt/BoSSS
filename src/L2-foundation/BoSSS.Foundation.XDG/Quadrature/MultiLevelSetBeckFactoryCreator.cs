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
using System.Security.Cryptography.X509Certificates;
using ilPSP.LinSolvers.monkey;
using System.Data;

namespace BoSSS.Foundation.XDG.Quadrature
{
    public class MultiLevelSetBeckFactoryCreator
    {
        LevelSetCombination[] phis;

        GridData grid;

        LevelSetData[] levelSets;

        public MultiLevelSetBeckFactoryCreator(LevelSetData[] levelSets)
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

        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, IQuadRuleFactory<QuadRule> backupFactory)
        {
            CombinedID id = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = jmp0,
                LevSet1 = levSetIndex1,
                Jmp1 = jmp1
            };
            LevelSetCombination lscomb = new LevelSetCombination(id,
                (LevelSet)levelSets[levSetIndex0].LevelSet,
                (LevelSet)levelSets[levSetIndex1].LevelSet);
            return new BeckQuadratureFactory(backupFactory, new BeckEdgeScheme(levelSets, lscomb,false), levelSets,id);
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
            LevelSetCombination lscomb = new LevelSetCombination(id0, 
                (LevelSet) levelSets[levSetIndex0].LevelSet,
                (LevelSet) levelSets[levSetIndex1].LevelSet);
            lscomb.sign0 = 0;
            return new BeckQuadratureFactory(backupFactory, new BeckSurfaceScheme(levelSets, lscomb,true), levelSets,id0);
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
            LevelSetCombination lscomb = FindPhi(id);

            return new BeckQuadratureFactory(backupFactory, new BeckVolumeScheme(levelSets, lscomb,true),levelSets,id);
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

        public IQuadRuleFactory<QuadRule> GetIntersectionFactory( int levSetIndex0, int levSetIndex1, IQuadRuleFactory<QuadRule> backupFactory) {

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

            CombinedID id0 = new CombinedID
            {
                LevSet0 = levSetIndex0,
                Jmp0 = JumpTypes.Heaviside, // does not matter here, we need the id for the special edge detector
                LevSet1 = levSetIndex1,
                Jmp1 = JumpTypes.Heaviside
            };

            LevelSetCombination lscomb = new LevelSetCombination(id0,
                (LevelSet)levelSets[levSetIndex0].LevelSet,
                (LevelSet)levelSets[levSetIndex1].LevelSet);
            lscomb.sign0 = 0;
            lscomb.sign1 = 0;
            return new BeckQuadratureFactory(backupFactory, new BeckSurfaceScheme(levelSets, lscomb, true), levelSets, id0);
            //var surfaceScheme = new BruteForceZeroScheme(Phi, Det);
            //throw new NotImplementedException();
            //return new BruteForceQuadratureFactory(backupFactory, levelSets, id, surfaceScheme, 400);
        }

    }

    /// <summary>
    /// Just a hack to enable special handling of cells where levelset=1 lies on an edge
    /// </summary>
    public static class BeckSettingsOverride {

        /// <summary>
        /// switch to enable special handling for cells, where the levelset=1 lies on an edge
        /// </summary>
        public static bool doubleCutCellOverride = false;
    }

    class BeckQuadratureFactory : IQuadRuleFactory<QuadRule>
    {
        BeckBaseScheme scheme;

        MultiLevelSetOnEdgeDetector detector;

        LevelSetData[] data;

        IQuadRuleFactory<QuadRule> fallBackFactory;

        public BeckQuadratureFactory(IQuadRuleFactory<QuadRule> fallBackFactory, BeckBaseScheme scheme, LevelSetTracker.LevelSetData[] data, CombinedID id)
        {
            this.fallBackFactory = fallBackFactory;
            this.scheme = scheme;
            this.data = data;
            detector = new MultiLevelSetOnEdgeDetector(data, id);
        }

        //public BeckQuadratureFactory(LevelSetTracker.LevelSetData[] data, CombinedID id, BeckBaseScheme scheme) : this(scheme,data)
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
                    if (detector != null)
                    {
                        ExecutionMask singleMask;
                        bool special;
                        bool active;
                        if (mask is CellMask cm)
                        {
                            special = detector.IsSpecialCell(i);
                            active = detector.IsActiveCell(i);
                            singleMask = new CellMask(cm.GridData, singleChunk, MaskType.Geometrical);
                        }
                        else
                        {
                            EdgeMask em = (EdgeMask)mask;
                            special = detector.IsSpecialEdge(i);
                            active = detector.IsActiveEdge(i);
                            singleMask = new EdgeMask(em.GridData, singleChunk, MaskType.Geometrical);
                        }

                        // detected a "special" cell, treat differently
                        if (special)
                        {
                            QuadRule specialRule;
                            specialRule = GetSpecialQuadRule(i, singleMask, order, active);

                            // add to List and jump to next chunk
                            ChunkRulePair<QuadRule> pairSpecial = new ChunkRulePair<QuadRule>(singleChunk, specialRule);
                            rule.Add(pairSpecial);
                            continue;
                        }
                    }
                    ChunkRulePair<QuadRule> pair = new ChunkRulePair<QuadRule>(singleChunk, scheme.GetQuadRule(i, order));
                    rule.Add(pair);
                }
            }
            return rule;
        }

        private QuadRule GetSpecialQuadRule(int i, ExecutionMask singleMask, int order, bool activeChunk)
        {
            QuadRule specialRule;

            switch (scheme)
            {
                case BeckZeroScheme b:
                    // to avoid deadlocks? construct here
                    specialRule = ConstructSpecialQuadRule(i, singleMask, order);

                    // skip cells not in active phase, could otherwise lead to problems, when there are kinks in the levelset (i.e. different contact angle from right and left)
                    if (activeChunk)
                    {
                        // do nothing, already built
                    }
                    else
                    {
                        specialRule = QuadRule.CreateEmpty(b.GetRefElement(), 1, b.GetRefElement().SpatialDimension);
                        specialRule.Nodes.LockForever();
                    }
                    break;
                case BeckEdgePointScheme b:
                    specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;
                    break;
                case BeckEdgeScheme b:
                    {
                        specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;

                        // remove edge scheme from "true" coinciding edge
                        int jCell = ((GridData)singleMask.GridData).Edges.CellIndices[i, 0];
                        int specialFace = detector.GetSpecialFace(jCell);
                        int iFace = ((GridData)singleMask.GridData).Edges.FaceIndices[i, 0];
                        bool coinciding = specialFace == iFace;
                        if (coinciding)
                        {
                            specialRule = QuadRule.CreateEmpty(specialRule.RefElement, 1, specialRule.RefElement.SpatialDimension);
                            specialRule.Nodes.LockForever();
                        }
                        break;
                    }
                case BeckSurfaceScheme b:
                    if (fallBackFactory.RefElement == scheme.GetRefElement())
                    {
                        // skip cells not in active phase
                        if (activeChunk)
                        {
                            specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;
                        }
                        else
                        {
                            specialRule = QuadRule.CreateEmpty(b.GetRefElement(), 1, b.GetRefElement().SpatialDimension);
                            specialRule.Nodes.LockForever();
                        }
                    }
                    else
                    {
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
                        if (ThisConform & OtherConform)
                        {
                            if (ThisCellGlob < OtherCellGlob)
                            {
                                // do nothing, already built
                            }
                            else
                            {
                                specialRule = QuadRule.CreateEmpty(b.GetRefElement(), 1, b.GetRefElement().SpatialDimension);
                                specialRule.Nodes.LockForever();
                            }
                        }
                        else if (ThisConform & !OtherConform)
                        {
                            // do nothing, already built
                        }
                        else if (!ThisConform & OtherConform)
                        {
                            specialRule = QuadRule.CreateEmpty(b.GetRefElement(), 1, b.GetRefElement().SpatialDimension);
                            specialRule.Nodes.LockForever();
                        }
                        else
                        {
                            throw new NotSupportedException(String.Format("Error in cell {0}, {1}: Only one cell should have the hanging node.", i, OtherCell));
                        }
                    }
                    break;
                case BeckVolumeScheme b:
                    // skip cells not in active phase
                    if (activeChunk)
                    {
                        specialRule = fallBackFactory.GetQuadRuleSet(singleMask, order).Single().Rule;
                    }
                    else
                    {
                        specialRule = QuadRule.CreateEmpty(b.GetRefElement(), 1, b.GetRefElement().SpatialDimension);
                        specialRule.Nodes.LockForever();
                    }
                    break;
                default:
                    throw new NotImplementedException();
            }
            return specialRule;
        }

        private QuadRule ConstructSpecialQuadRule(int i, ExecutionMask singleMask, int order)
        {
            QuadRule specialRule;
            if (singleMask.GetType() != typeof(CellMask))
                throw new NotSupportedException();

            // determine face
            int iFace = detector.GetSpecialFace(i);
            int iEdge = ((GridData)singleMask.GridData).Cells.GetEdgesForFace(i, iFace, out int InOrOut, out int[] FurtherEdges);
            // select all edges, relevant in cells with hanging nodes
            singleMask = new EdgeMask(singleMask.GridData, Chunk.GetSingleElementChunk(iEdge), MaskType.Geometrical);
            if (FurtherEdges != null)
            {
                specialRule = QuadRule.CreateEmpty(scheme.GetRefElement(), 1, scheme.GetRefElement().SpatialDimension);
                specialRule.Nodes.LockForever();

                foreach (var edg in FurtherEdges)
                {
                    singleMask = singleMask.Union(new EdgeMask(singleMask.GridData, Chunk.GetSingleElementChunk(edg), MaskType.Geometrical));
                }
            }
            // construct edgerule
            var specialRule_t = fallBackFactory.GetQuadRuleSet(singleMask, order);
            List<NodeSet> specialNodes = new List<NodeSet>();
            List<MultidimensionalArray> specialWeights = new List<MultidimensionalArray>();
            int NoOfNodes = 0;
            foreach (var crp in specialRule_t)
            {
                var rule_t = crp.Rule;
                foreach (int edg in crp.Chunk.Elements)
                {
                    int trf;
                    if (singleMask.GridData.iGeomEdges.CellIndices[edg, 0] == i)
                    {
                        trf = singleMask.GridData.iGeomEdges.Edge2CellTrafoIndex[edg, 0];
                    }
                    else
                    {
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
            for (int j = 0; j < specialNodes.Count; j++)
            {
                int count = specialNodes[j].Lengths[0];
                specialNodes_t.ExtractSubArrayShallow(new int[] { offset, 0 }, new int[] { offset + count - 1, singleMask.GridData.SpatialDimension - 1 }).Acc(1.0, specialNodes[j]);
                specialWeights_t.ExtractSubArrayShallow(new int[] { offset }, new int[] { offset + count - 1 }).Acc(1.0, specialWeights[j]);
                offset += count;
            }
            specialRule = new QuadRule();
            specialRule.Nodes = new NodeSet(scheme.GetRefElement(), specialNodes_t, true);
            specialRule.Weights = specialWeights_t;

            return specialRule;

        }
    }
    internal interface ISchemeWO
    {
        RefElement ReferenceElement { get; }

        QuadRule GetQuadRule(int j,int order);

        void Initialize(int resolution);
    }
    internal abstract class BeckBaseScheme : ISchemeWO
    {
        public LevelSetData[] data;

        public LevelSetCombination lscomb;
        public int D;
        public bool isGlobalMode = false;


        RefElement ISchemeWO.ReferenceElement => GetRefElement();

        public abstract RefElement GetRefElement();

        public BeckBaseScheme(LevelSetData[] data, LevelSetCombination lscomb, bool isGlobalMode=false)
        {
            this.data = data;
            this.isGlobalMode= isGlobalMode;

            this.lscomb = lscomb;
            this.D = ((LevelSet)data[0].LevelSet).Basis.GridDat.SpatialDimension;
        }

        public (int noOfNodes,int subdiv) GetNoOfQuadNodes(int quadorder)
        {
            if (quadorder == 0)
            {
                return (1, 0);
            }
            //we need at least so much nodes
            int neededNodes = (int)(quadorder+1) / 2;
            //max noOfNodes supported is 4
            int neededSubdiv = (int) (neededNodes) / 4;
            if(neededSubdiv > 0)
            {
                return (4, neededSubdiv+0);
            }
            else
            {
                return (neededNodes, 0);
            }
        }
        public QuadRule GetQuadRule(int j, int order)
        {
            //get GridData 
            var gdat = (GridData)((LevelSet)this.data[0].LevelSet).Basis.GridDat;

            //Get Phi Evaluators
            (var phi0, var phi1) = GetPhiEval(j);

            //Get NoOfQuadNodes
            (int noOfNodes, int subdiv) = GetNoOfQuadNodes(order);

            //Get a quadrule finder
            Quadrater finder = new Quadrater();

            //Get the cell on which we integrate, here depending on the scheme the right integration domain is chosen
            HyperRectangle cell = GetCell(j);

            //get the symbols
            (var s1, var s2) = GetSymbols(lscomb);
            // Get The Quadrature rule From Intersectingquadrature
            QuadratureRule ruleQ = (default);
            
            try
            {
                //find the quad rule
                ruleQ = finder.FindRule(phi0, s1, phi1, s2, cell, noOfNodes, subdiv);
            }
            catch(Exception ex)
            {
                throw ex;
            }
            

            // Creates a QuadRule Object <- The One BoSSS uses
            QuadRule rule = QuadRule.CreateEmpty(GetRefElement(), ruleQ.Count, cell.Dimension, true);
            rule.OrderOfPrecision = order;

            //transfer from Beck into BoSSS structure
            //loop over all quadrature Nodes
            for (int i = 0; i < ruleQ.Count; ++i)
            {
                QuadratureNode qNode = ruleQ[i];
                PointsToNodes(rule.Nodes.ExtractSubArrayShallow(i,-1), qNode, cell.Dimension, j);
                rule.Weights[i] = qNode.Weight;
            }
            
            //do some scaling of the Nodes
            if (isGlobalMode) //used for Surface and Volume Rules
            {
                //Tensor1 root0 = Tensor1.Vector(0.5);
                //Tensor1 root1 = Tensor1.Vector(-0.5);

                //var ret0 = phi0.Evaluate(root0);
                //var ret1 = phi1.Evaluate(root1);

                //We transform the Nodes, which are in global coordainates into Local onese
                var NodesOut = MultidimensionalArray.Create(1,rule.Nodes.Lengths[0], rule.Nodes.Lengths[1]);
                gdat.TransformGlobal2Local(rule.Nodes, NodesOut, j, 1, 0);
                rule.Nodes.Set(NodesOut.ExtractSubArrayShallow(0,-1,-1));
                rule.Nodes.LockForever();

                if (rule.Nodes.Lengths[0] != 0)
                {
                    //We need to scale the weights, as they will be multiplied by the determinant of the jacobian

                    //Get the Scaling
                    ScaleWeights(rule,j);

                    
                }
            }
            else
            {
                //Get the Scaling
                double scaling = GetScaling(j);
                rule.Weights.Scale(1 / scaling);
                rule.Nodes.LockForever();
            }
            
            //rule.Nodes.SaveToTextFile($"quadNodes_{this.ToString()}_{lscomb.ID.LevSet0}{s1.ToString()}_{lscomb.ID.LevSet1}{s2.ToString()}.txt");
            return rule;

        }

        public abstract void ScaleWeights(QuadRule rule, int j);
        public abstract void PointsToNodes(MultidimensionalArray multidimensionalArray, QuadratureNode qNode, int d, int j);

        /// <summary>
        /// Creates a Hyperrectangle corresponding to a Grid Cell
        /// </summary>
        /// <param name="j">Cell No.</param>
        /// <param name="gdat">Grid </param>
        /// <returns></returns>
        public HyperRectangle CellToGlobalHR(int j, GridData gdat)
        {
            var eC1 = gdat.Edges.CellIndices[j, 0];
            var eC2 = gdat.Edges.CellIndices[j, 1];
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
            var HR = new UnitCube(D);
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
        /// <summary>
        /// Creates a HyperRectangle From an Edge
        /// </summary>
        /// <param name="j">Edge No.</param>
        /// <param name="gdat">Grid</param>
        /// <returns></returns>
        public HyperRectangle EdgeToGlobalHR(int j, GridData gdat)
        {
            var jNeighCell = gdat.Edges.CellIndices[j, 0];
            var e2C = gdat.iGeomEdges.Edge2CellTrafos[gdat.Edges.Edge2CellTrafoIndex[j, 0]];

            var center = new Vector(gdat.SpatialDimension - 1);

            // create the Hyperrectangle
            var ret = new UnitCube(gdat.SpatialDimension-1);

            //get the Coordinate of the Center
            var centerInCell = e2C.Transform(center);

            //Helper Function to transform a Vector into a MA of lengts[1,x.Dim]
            MultidimensionalArray VectorToMA(Vector x)
            {
                MultidimensionalArray MA = MultidimensionalArray.Create(1, x.Dim);
                for (int i = 0; i < x.Dim; i++)
                {
                    MA[0, i] = x[i];
                }
                return MA;
            }

            MultidimensionalArray LocalVerticesIn = VectorToMA(centerInCell);
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.Create(1, LocalVerticesIn.Lengths[1]);


            gdat.TransformLocal2Global(LocalVerticesIn, GlobalVerticesOut, jNeighCell);

            //asign it
            //Pick the relevant entries of the gradient
            int count = 0;
            for (int d = 0; d < gdat.SpatialDimension; d++)
            {
                if (centerInCell[d] == 0)
                {
                    ret.Center[count] = GlobalVerticesOut.ExtractSubArrayShallow(0, -1)[d];
                    count++;
                }
            }
            //ret.Center = lSEvalUtil.ToTensor1(GlobalVerticesOut.ExtractSubArrayShallow(0, -1));

            var corner = new Vector(gdat.SpatialDimension - 1);
            //Get Diameters from information of Corners
            for (int iD = 0; iD < ret.Dimension; iD++)
            {
                //get first corner
                corner[iD] = 1;
                var cornerInCell = e2C.Transform(corner);
                LocalVerticesIn = VectorToMA(cornerInCell);
                MultidimensionalArray corner1 = MultidimensionalArray.Create(1, LocalVerticesIn.Lengths[1]);
                gdat.TransformLocal2Global(LocalVerticesIn, corner1, jNeighCell);

                //get second corner
                corner[iD] = -1;
                cornerInCell = e2C.Transform(corner);
                LocalVerticesIn = VectorToMA(cornerInCell);
                MultidimensionalArray corner2 = MultidimensionalArray.Create(1, LocalVerticesIn.Lengths[1]);
                gdat.TransformLocal2Global(LocalVerticesIn, corner2, jNeighCell);

                //substract from each other
                corner1.Acc(-1, corner2);

                //set diameter
                ret.Diameters[iD] = corner1.AbsSum(); //only one entry should be non zero.

                //reset
                corner[iD] = 0;
            }
            return ret;

        }
        /// <summary>
        /// This method reutrn a Hyperrectangle that is exactly the domain we want to integrate
        /// </summary>
        /// <returns></returns>
        public abstract HyperRectangle GetCell(int j);

        /// <summary>
        /// gives the scaling
        /// </summary>
        /// <param name="j"></param>
        /// <returns></returns>
        public abstract double GetScaling(int j);

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

            Symbol s1 = (sign0 == -1.0) ? Symbol.Minus : Symbol.Zero;
            s1 = (sign0 == 1.0) ? Symbol.Plus : s1;

            Symbol s2 = (sign1 == -1.0) ? Symbol.Minus : Symbol.Zero;
            s2 = (sign1 == 1.0) ? Symbol.Plus : s2;

            return (s1, s2);
        }
    }

    internal class BeckVolumeScheme : BeckBaseScheme
    {
        public BeckVolumeScheme(LevelSetData[] data, LevelSetCombination lscomb, bool isGlobalMode = false) : base( data, lscomb, isGlobalMode)
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
        public override HyperRectangle GetCell(int j)
        {
            if (base.isGlobalMode)
            {
                return CellToGlobalHR(j, (GridData) ((LevelSet)data[lscomb.ID.LevSet0].LevelSet).GridDat);
            }
            else
            {
                return new UnitCube(D);
            }
            
        }
        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            if (base.isGlobalMode)
            {
                return (new lSEvalGlobal((LevelSet)data[lscomb.ID.LevSet0].LevelSet, j), new lSEvalGlobal((LevelSet)data[lscomb.ID.LevSet1].LevelSet, j));
            }
            else
            {
                return (new lSEvalLocal((LevelSet)data[lscomb.ID.LevSet0].LevelSet, j), new lSEvalLocal((LevelSet)data[lscomb.ID.LevSet1].LevelSet, j));
            }
            
        }

        public override double GetScaling(int j)
        {
            return 1;
        }

        public override void ScaleWeights(QuadRule rule, int j)
        {
            var gdat = ((LevelSet)data[lscomb.ID.LevSet0].LevelSet).GridDat;
            var jacDet = gdat.JacobianDeterminat.GetValue_Cell(rule.Nodes, j, 1);
            for (int iWeight = 0; iWeight < rule.Weights.Lengths[0]; iWeight++)
            {
                //rule.Weights[iWeight] = rule.Weights[iWeight] / jacDet[0, iWeight];
                rule.Weights[iWeight] = rule.Weights[iWeight]/ jacDet[0, iWeight];
            }
        }

        public override void PointsToNodes(MultidimensionalArray rNode, QuadratureNode qNode, int D, int j)
        {
            for (int d = 0; d < D; d++)
            {
                rNode[d] = qNode.Point[d];
            }
        }
    }
    internal class BeckEdgeScheme : BeckBaseScheme
    {
        public BeckEdgeScheme( LevelSetData[] data, LevelSetCombination lscomb, bool isGlobalMode = false) : base(data, lscomb, isGlobalMode)
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
        public override void PointsToNodes(MultidimensionalArray rNode, QuadratureNode qNode, int D,int j)
        {
            for (int d = 0; d < D; d++)
            {
                rNode[d] = qNode.Point[d];
            }
            //var gdat = (GridData)((LevelSet)data[0].LevelSet).Basis.GridDat;
            //var jNeighCell = gdat.Edges.CellIndices[j, 0];
            //var e2C = gdat.iGeomEdges.Edge2CellTrafos[gdat.Edges.Edge2CellTrafoIndex[j, 0]];
            
            //var trafoNode= e2C.Transform(lSEvalUtil.TensorToVector(qNode.Point));
            //var globTrafoNode = lSEvalUtil.VectorToMA(trafoNode);
            //if (isGlobalMode)
            //{
            //    gdat.TransformLocal2Global(lSEvalUtil.VectorToMA(trafoNode), globTrafoNode, jNeighCell);
            //}
            
            //for (int d = 0; d < gdat.SpatialDimension; d++)
            //{
            //    rNode[d] = globTrafoNode[0,d];
            //}
        }
        public override HyperRectangle GetCell( int j)
        {
            if (isGlobalMode)
            {
                return EdgeToGlobalHR(j, (GridData)((LevelSet)data[lscomb.ID.LevSet0].LevelSet).GridDat);
            }
            else
            {
                return new UnitCube(D-1);
            }
        }
        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            return (new lSEvalEdge((LevelSet)data[lscomb.ID.LevSet0].LevelSet, j), new lSEvalEdge((LevelSet)data[lscomb.ID.LevSet1].LevelSet, j));
        }

        public override double GetScaling(int j)
        {
            return 1;
        }
        public override void ScaleWeights(QuadRule rule, int j)
        {
            var gdat = ((LevelSet)data[lscomb.ID.LevSet0].LevelSet).GridDat;
            var jacDet = gdat.JacobianDeterminat.GetValue_Cell(rule.Nodes, j, 1);
            for (int iWeight = 0; iWeight < rule.Weights.Lengths[0]; iWeight++)
            {
                rule.Weights[iWeight] = rule.Weights[iWeight] / jacDet[0, iWeight];
            }
        }
    }
    internal class BeckSurfaceScheme : BeckBaseScheme
    {
        public BeckSurfaceScheme(LevelSetData[] data, LevelSetCombination lscomb, bool isGlobalMode = false) : base(data, lscomb, isGlobalMode)
        {
        }
        public override void PointsToNodes(MultidimensionalArray rNode, QuadratureNode qNode, int D, int j)
        {
            for (int d = 0; d < D; d++)
            {
                rNode[d] = qNode.Point[d];
            }
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
        public override HyperRectangle GetCell(int j)
        {
            if (isGlobalMode)
            {
                return CellToGlobalHR(j, (GridData)((LevelSet)data[lscomb.ID.LevSet0].LevelSet).GridDat);
            }
            else
            {
                return new UnitCube(D);
            }
        }
        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            if(isGlobalMode)
            {
                return (new lSEvalGlobal((LevelSet)data[lscomb.ID.LevSet0].LevelSet, j), new lSEvalGlobal((LevelSet)data[lscomb.ID.LevSet1].LevelSet, j));
            }
            else
            {
                return (new lSEvalLocal((LevelSet)data[lscomb.ID.LevSet0].LevelSet, j), new lSEvalLocal((LevelSet)data[lscomb.ID.LevSet1].LevelSet, j));
            }
            
        }

        public override double GetScaling(int j)
        {
            if (isGlobalMode)
            {
                return 1;
            }
            else
            {
                var ls = (LevelSet)data[lscomb.ID.LevSet0].LevelSet;
                var grd = ls.GridDat;
                var ret = grd.Jacobian.GetValue_Cell(GetRefElement().Center, j, 1);
                return ret[0, 0, 0];
            }
            
        }
        public override void ScaleWeights(QuadRule rule, int j)
        {
            var gdat = ((LevelSet)data[lscomb.ID.LevSet0].LevelSet).GridDat;
            var jacDet = gdat.JacobianDeterminat.GetValue_Cell(rule.Nodes, j, 1);
            for (int iWeight = 0; iWeight < rule.Weights.Lengths[0]; iWeight++)
            {
                rule.Weights[iWeight] = rule.Weights[iWeight] / jacDet[0, iWeight];
            }
        }
    }

    internal class BeckZeroScheme : BeckBaseScheme
    {
        public BeckZeroScheme(LevelSetData[] data, LevelSetCombination lscomb, bool isGlobalMode = false) : base(data, lscomb, isGlobalMode)
        {
        }

        public override HyperRectangle GetCell(int j)
        {
            throw new NotImplementedException();
        }

        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            throw new NotImplementedException();
        }

        public override RefElement GetRefElement()
        {
            throw new NotImplementedException();
        }

        public override double GetScaling(int j)
        {
            throw new NotImplementedException();
        }

        public override void PointsToNodes(MultidimensionalArray multidimensionalArray, QuadratureNode qNode, int d, int j)
        {
            throw new NotImplementedException();
        }

        public override void ScaleWeights(QuadRule rule, int j)
        {
            throw new NotImplementedException();
        }
    }
    internal class BeckEdgePointScheme : BeckBaseScheme
    {
        public BeckEdgePointScheme(LevelSetData[] data, LevelSetCombination lscomb, bool isGlobalMode = false) : base(data, lscomb, isGlobalMode)
        {
        }

        public override HyperRectangle GetCell(int j)
        {
            throw new NotImplementedException();
        }

        public override (IScalarFunction phi0, IScalarFunction phi1) GetPhiEval(int j)
        {
            throw new NotImplementedException();
        }

        public override RefElement GetRefElement()
        {
            throw new NotImplementedException();
        }

        public override double GetScaling(int j)
        {
            throw new NotImplementedException();
        }

        public override void PointsToNodes(MultidimensionalArray multidimensionalArray, QuadratureNode qNode, int d, int j)
        {
            throw new NotImplementedException();
        }

        public override void ScaleWeights(QuadRule rule, int j)
        {
            throw new NotImplementedException();
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
            NodeSet NS = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //var inp = MultidimensionalArray.Create(1, m_levelSet.Basis.GridDat.SpatialDimension - 1);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            var inp = MultidimensionalArray.Create(1, 1);
            if (x.M == m_levelSet.Basis.GridDat.SpatialDimension)
            {
                m_levelSet.Evaluate(j, 1, NS, ev);
            }
            else
            {
                m_levelSet.EvaluateEdge(j, 1, NS, ev, inp);
            }
            return ev[0, 0];
        }

        (double evaluation, Tensor1 gradient) IScalarFunction.EvaluateAndGradient(Tensor1 x)
        {
            int D = m_levelSet.Basis.GridDat.SpatialDimension;
            NodeSet NS = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var inp = MultidimensionalArray.Create(1, 1);
            var ev = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray gradIn = MultidimensionalArray.Create(1, 1, D);

            MultidimensionalArray gradOut = MultidimensionalArray.Create(D - 1);
            MultidimensionalArray hessOut = MultidimensionalArray.Create(D - 1, D - 1);
            if (x.M == D)
            {
                m_levelSet.Evaluate(j, 1, NS, ev);
                m_levelSet.EvaluateGradient(j, 1, NS, grad);
                return (ev[0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)));
            }
            else
            {
                m_levelSet.EvaluateEdge(j, 1, NS, ev, inp, GradientIN: grad, GradientOT: gradIn);
                (gradOut, hessOut) = DecideForEntries(grad);
                return (ev[0, 0], lSEvalUtil.ToTensor1(gradOut));
            }
            
            
        }

        (double evaluation, Tensor1 gradient, Tensor2 hessian) IScalarFunction.EvaluateAndGradientAndHessian(Tensor1 x)
        {
            int D = m_levelSet.Basis.GridDat.SpatialDimension;
            NodeSet NS = lSEvalUtil.NSFromTensor(x, m_levelSet);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var inp = MultidimensionalArray.Create(1, 1);
            var ev = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray gradIn = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray hess = MultidimensionalArray.Create(1, 1,D,D);
            MultidimensionalArray gradOut = MultidimensionalArray.Create(D-1);
            MultidimensionalArray hessOut = MultidimensionalArray.Create(D - 1, D - 1);
            if (x.M == m_levelSet.Basis.GridDat.SpatialDimension)
            {
                m_levelSet.Evaluate(j, 1, NS, ev);
                m_levelSet.EvaluateGradient(j, 1, NS, grad);
                m_levelSet.EvaluateHessian(j, 1, NS, hess);
                return (ev[0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0,0,-1)), lSEvalUtil.ToTensor2(hess.ExtractSubArrayShallow(0, 0, -1,-1)));
            }
            else
            {
                m_levelSet.EvaluateEdge(j, 1, NS, ev, inp, GradientIN: grad, GradientOT: gradIn);
                (gradOut, hessOut) = DecideForEntries(grad);
                return (ev[0, 0], lSEvalUtil.ToTensor1(gradOut), lSEvalUtil.ToTensor2(hessOut));
            }


            
            

            //no Hessian Evaluation on Edges
            //m_levelSet.EvaluateHessian(i0, 1, NS, hess);

            

        }

        private (MultidimensionalArray gradOut, MultidimensionalArray hessOut) DecideForEntries(MultidimensionalArray grad)
        {
            //output of gradient must be dim-1
            var gradOut = MultidimensionalArray.Create(grad.Lengths[2] - 1);
            //Dimensions must agree, will stay zero as no Hessian is provided on edges
            var hessOut = MultidimensionalArray.Create(grad.Lengths[2] - 1, grad.Lengths[2] - 1);
            
            // we use the Center of the edge in Local Coordinates to decide which Entries are relevant
            var gdat = (GridData)this.m_levelSet.Basis.GridDat;
            var jNeighCell = gdat.Edges.CellIndices[j, 0];
            var e2C = gdat.iGeomEdges.Edge2CellTrafos[gdat.Edges.Edge2CellTrafoIndex[j, 0]];


            var center = new Vector(gdat.SpatialDimension - 1);
            //get the Coordinate of the Center
            var centerInCell = e2C.Transform(center);


            //Pick the relevant entries of the gradient
            int count = 0;
            for (int d = 0; d < gdat.SpatialDimension; d++)
            {
                if (centerInCell[d] == 0)
                {
                    gradOut[count] = grad[0,0,d];
                    count++;
                }
            }
            return (gradOut, hessOut);
        }
    }
    internal class lSEvalLocal : lSEvalBase
    {
        public lSEvalLocal(LevelSet levelSet, int j) : base(levelSet, j)
        {
        }

        public override NodeSet GetNodeSet(Tensor1 x)
        {
            return lSEvalUtil.NSFromTensor(x, m_levelSet);
        }
    }

    internal abstract class lSEvalBase : IScalarFunction
    {
        //levelSet which is evaluated by this class
        public LevelSet m_levelSet;
        //cell the Level Set is evaluated on
        public int j;
        int IScalarFunction.M => m_levelSet.Basis.GridDat.SpatialDimension;

        public lSEvalBase(LevelSet levelSet, int j)
        {
            m_levelSet = levelSet;
            this.j = j;
        }

        double IScalarFunction.Evaluate(Tensor1 x)
        {
            NodeSet NS = GetNodeSet(x);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            m_levelSet.Evaluate(j, 1, NS, ev);
            return ev[0, 0];
        }

        public abstract NodeSet GetNodeSet(Tensor1 x);

        (double evaluation, Tensor1 gradient) IScalarFunction.EvaluateAndGradient(Tensor1 x)
        {
            NodeSet NS = GetNodeSet(x);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            m_levelSet.Evaluate(j, 1, NS, ev);

            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
            m_levelSet.EvaluateGradient(j, 1, NS, grad);

            return (ev[0, 0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)));

        }

        (double evaluation, Tensor1 gradient, Tensor2 hessian) IScalarFunction.EvaluateAndGradientAndHessian(Tensor1 x)
        {

            NodeSet NS = GetNodeSet(x);
            //Evaluates the LevelSet using the NodeSet and the Cell Number
            var ev = MultidimensionalArray.Create(1, 1);
            m_levelSet.Evaluate(j, 1, NS, ev);


            MultidimensionalArray grad = MultidimensionalArray.Create(1, 1, x.M);
            MultidimensionalArray hess = MultidimensionalArray.Create(1, 1, x.M, x.M);
            m_levelSet.EvaluateGradient(j, 1, NS, grad);
            m_levelSet.EvaluateHessian(j, 1, NS, hess);

            return (ev[0, 0], lSEvalUtil.ToTensor1(grad.ExtractSubArrayShallow(0, 0, -1)), lSEvalUtil.ToTensor2(hess.ExtractSubArrayShallow(0, 0, -1, -1)));

        }
    }

    internal class lSEvalGlobal : lSEvalBase
    {
        public lSEvalGlobal(LevelSet levelSet, int j):base(levelSet, j) 
        {
        }
        public override NodeSet GetNodeSet(Tensor1 x)
        {
            
            return lSEvalUtil.NSFromGlobalTensor(x, m_levelSet, j);
        }

    }
    internal static class lSEvalUtil {
        
        public static Vector TensorToVector(Tensor1 x)
        {
            Vector point;
            if (x.M == 1)
            {
                point = new Vector(x[0]);
            }
            else if(x.M==2) {
                point = new Vector(x[0], x[1]);
            }
            else
            {
                point = new Vector(x[0], x[1], x[2]);
            }
            return point;
        }
        //Helper Function to transform a Vector into a MA of lengts[1,x.Dim]
        public static MultidimensionalArray VectorToMA(Vector x)
        {
            MultidimensionalArray MA = MultidimensionalArray.Create(1, x.Dim);
            for (int i = 0; i < x.Dim; i++)
            {
                MA[0, i] = x[i];
            }
            return MA;
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

        public static NodeSet NSFromTensor(Tensor1 x, LevelSet m_levelSet)
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
            //m_levelSet.GridDat.LocatePoint(GlobIn.ExtractSubArrayShallow(0, -1).To1DArray(), out long id, out long i0, out bool IsInside, out bool Onm_levelSetProces);
            //gives the local Coordinates
            //m_levelSet.GridDat.TransformGlobal2Local(GlobIn, LocOut, (int)i0, 1, 0);


            //gives the NodeSet using the local Coords
            if(x.M == m_levelSet.Basis.GridDat.SpatialDimension)
            {
                return new NodeSet(m_levelSet.GridDat.iGeomCells.RefElements[0], GlobIn, true);

            }
            else if(x.M == m_levelSet.Basis.GridDat.SpatialDimension-1)
            {
                return new NodeSet(m_levelSet.GridDat.iGeomEdges.EdgeRefElements[0], GlobIn, true);

            }
            else
            {
                return new NodeSet(Grid.RefElements.Point.Instance, GlobIn, true);
            }
        }
        public static NodeSet NSFromGlobalTensor(Tensor1 x, LevelSet m_levelSet, int j)
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
            //m_levelSet.GridDat.LocatePoint(GlobIn.ExtractSubArrayShallow(0, -1).To1DArray(), out long id, out long i0, out bool IsInside, out bool Onm_levelSetProces);
            //gives the local Coordinates
            m_levelSet.GridDat.TransformGlobal2Local(GlobIn, LocOut, j, 1, 0);

            //gives the NodeSet using the local Coords
            var NS = new NodeSet(m_levelSet.GridDat.iGeomCells.RefElements[0], LocOut.ExtractSubArrayShallow(0, -1, -1), true);

            return NS;
        }
    }
}
