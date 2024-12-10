using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.XQuadFactoryHelperBase;

namespace CutEdgeQuadrature {

    interface IEdgeQuadratureTest2D : IEdgeQuadratureTest {

        _2D LevelSet { get; }
    }

    interface IEdgeQuadratureTest3D : IEdgeQuadratureTest{

        _3D LevelSet { get; }
    }

    interface IEdgeQuadratureTest {

        double EdgeArea { get; }

        MomentFittingVariants MomentFittingVariant { get; set; }
    }
}
