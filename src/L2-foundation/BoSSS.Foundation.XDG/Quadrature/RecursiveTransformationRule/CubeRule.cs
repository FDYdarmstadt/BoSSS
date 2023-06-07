using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG.Quadrature.RecursiveMappingRule {
    class CubeRule<T> {

        //SayeIntegrand<T, SayeArgument<T>> spaceDecomposer;

        //IList<Transformation> trafos;

        public CubeRule(LevelSetTracker.LevelSetData lsData, HMF.LineSegment.IRootFindingAlgorithm rooter) {
            throw new NotImplementedException();
            //spaceDecomposer = new SayeGaussRule_Cube(lsData, rooter, SayeGaussRule_Cube.QuadratureMode.NegativeVolume);
        }

        TreeNode<T> DecomposeSpace() {
            throw new NotImplementedException();
            //spaceDecomposer
        }

        IList<Transformation> CreateTransformations() {
            throw new NotImplementedException();
        }
        
    }

    class Transformation {
        public QuadRule Transform(QuadRule rule) {
            throw new NotImplementedException();
        }
    }
}
