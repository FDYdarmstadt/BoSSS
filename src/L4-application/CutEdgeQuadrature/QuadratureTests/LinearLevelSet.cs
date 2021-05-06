using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CutEdgeQuadrature.QuadratureTests {
    class LinearLevelSet : IEdgeQuadratureTest3D {

        public LinearLevelSet() {
            
        }

        double Phi(double x, double y, double z) {
            return z + 0.1 + 0.15 * y;
        }

        public _3D LevelSet => Phi;

        public double EdgeArea => 1.8;

        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant => XQuadFactoryHelper.MomentFittingVariants.Saye;
    }
}
