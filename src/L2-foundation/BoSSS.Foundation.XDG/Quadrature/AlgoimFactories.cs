using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using MPI.Wrappers.Utils;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineAndPointQuadratureFactory;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;

namespace BoSSS.Foundation.XDG.Quadrature {

    //public sealed class Algoim : DynLibLoader {

        //    public Algoim() :
        //base(new string[] { "algoim.dll", "algoim_seq.so" },
        //      new string[2][][],
        //      new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.BoSSS_Prefix },
        //      Helper(), //new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix, PlatformID.Unix, PlatformID.Unix },
        //      new int[] { -1, -1 }) { }

        //}
    //}

    public static class AlgoimFactories {

        //This would return an factory object with the configuration 
        //input: level set, tolerance etc.
        public static void GetFactories() {



        }

        #region Edge rules


        class Factory {


            //public IQuadRuleFactory<QuadRule> GetSurfaceFactory() {
            //    return new IQuadRuleFactory<QuadRule>();
            //}

            //public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            //    return new IQuadRuleFactory<QuadRule>();
            //}

        }




        #endregion
    }



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



        public static IQuadRuleFactory<CellBoundaryQuadRule> AlgoimGaussRule_EdgeVolume3D(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder) {

            ISayeGaussEdgeRule rule = new SayeGaussRule_EdgeCube(
                _lsData,
                RootFinder,
                SayeGaussRule_Cube.QuadratureMode.PositiveVolume);
            return new SayeGaussEdgeRuleFactory(rule);
        }

    }

}


