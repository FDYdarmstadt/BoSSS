// Ignore Spelling: Algoim

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using MPI.Wrappers.Utils;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.Text;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineAndPointQuadratureFactory;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;

namespace BoSSS.Foundation.XDG.Quadrature {


    public class AlgoimFactories {


        public IQuadRuleFactory<QuadRule> GetSurfaceFactory() {
            return new Factory() {
                m_Owner = this  };
        }

        public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            return new Factory() {
                m_Owner = this
            };
        }

        //This would return an factory object with the configuration 
        //input: level set, tolerance etc.
        public AlgoimFactories(LevelSetTracker.LevelSetData ls, RefElement e) {
            refElement = e;
            lsData = ls;
        }

        LevelSetTracker.LevelSetData lsData;
        RefElement refElement;

        #region Edge rules


        class Factory : IQuadRuleFactory<QuadRule> {
            internal AlgoimFactories m_Owner;

            public RefElement RefElement => m_Owner.refElement;

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                            throw new NotImplementedException();

            }
            
//            {
//                if (!(mask is CellMask))
//                    throw new ArgumentException("Expecting a cell mask.");
//                if (mask.MaskType != MaskType.Geometrical)
//                    throw new ArgumentException("Expecting a geometrical mask.");


//                //int InternalSurfaceOrder = m_Owner.OrderToInternalOrder(RequestedOrder);
//                //int InternalVolumeOrder = InternalSurfaceOrder - 1;

//                //int FiledOrder;
//                //Debug.Assert(object.ReferenceEquals(this.Rules, m_Owner.m_VolumeRules) || object.ReferenceEquals(this.Rules, m_Owner.m_SurfaceRules));
//                //if (object.ReferenceEquals(this.Rules, m_Owner.m_VolumeRules))
//                //    // I'm a volume rule factory
//                //    FiledOrder = InternalVolumeOrder;
//                //else
//                //    // I'm a surface rule factory
//                //    FiledOrder = InternalSurfaceOrder;

//#if DEBUG
//                if (mask.Except(m_Owner.MaxGrid).NoOfItemsLocally > 0)
//                    throw new NotSupportedException("'mask' must be a subset of the cut cells, for my reference element.");
//#endif
//                if (!Rules.ContainsKey(FiledOrder))
//                    m_Owner.GetQuadRuleSet_Internal(InternalSurfaceOrder);

//                if (mask.NoOfItemsLocally == m_Owner.MaxGrid.NoOfItemsLocally) {
//                    // aggressive
//                    return Rules[FiledOrder];
//                } else {
//                    var Rule = Rules[FiledOrder];

//                    int L = mask.NoOfItemsLocally, H = Rule.Length;
//                    var Ret = new ChunkRulePair<QuadRule>[L];
//                    int h = 0;
//                    //for (int jsub = 0; jsub < L; jsub++) {
//                    //    int jCell = jsub2jcell[jsub];
//                    int jsub = 0;
//                    foreach (int jCell in mask.ItemEnum) {

//                        Debug.Assert(Rule[h].Chunk.Len == 1);

//                        while (jCell > Rule[h].Chunk.i0) {
//                            h++;
//                        }

//                        Debug.Assert(jCell == Rule[h].Chunk.i0);
//                        Ret[jsub] = Rule[h];
//#if DEBUG
//                        Ret[jsub].Rule.Weights.CheckForNanOrInf();
//                        Ret[jsub].Rule.Nodes.CheckForNanOrInf();
//#endif
//                        jsub++;
//                    }
//                    Debug.Assert(jsub == L);

//                    return Ret;
//                }

            
        }

    }




        #endregion
    



    //QuadRule RuleToRuleThemAll = QuadRule.CreateEmpty(refElement, count, spatialDim);
    //RuleToRuleThemAll.Nodes = new NodeSet(refElement, nodes, true);
    //RuleToRuleThemAll.Weights = weights;
    //        return RuleToRuleThemAll;



    class AlgoimEdgeRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {
        public RefElement RefElement => throw new NotImplementedException();


        public AlgoimEdgeRuleFactory() { }

        public int[] GetCachedRuleOrders() {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            throw new NotImplementedException();
        }



        //public static IQuadRuleFactory<CellBoundaryQuadRule> AlgoimGaussRule_EdgeVolume3D(
        //    LevelSetTracker.LevelSetData _lsData,
        //    IRootFindingAlgorithm RootFinder) {

        //    ISayeGaussEdgeRule rule = new SayeGaussRule_EdgeCube(
        //        _lsData,
        //        RootFinder,
        //        SayeGaussRule_Cube.QuadratureMode.PositiveVolume);
        //    return new SayeGaussEdgeRuleFactory(rule);
        //}

    }

}


