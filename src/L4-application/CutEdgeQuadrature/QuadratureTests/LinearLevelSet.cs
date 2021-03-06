﻿using BoSSS.Foundation.XDG;
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
            return z + 0.1 + 0.1 * y;
        }

        public _3D LevelSet => Phi;

        public double EdgeArea => 1.8;

        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant => XQuadFactoryHelper.MomentFittingVariants.Saye;
    }

    class SkewLinearLevelSet : IEdgeQuadratureTest3D {

        public SkewLinearLevelSet() {

        }

        double Phi(double x, double y, double z) {
            return z + x + y;
        }

        public _3D LevelSet => Phi;

        public double EdgeArea => 3.5;

        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant => XQuadFactoryHelper.MomentFittingVariants.Saye;
    }

    class QuadraticLevelSet : IEdgeQuadratureTest3D {

        public QuadraticLevelSet() {

        }

        double Phi(double x, double y, double z) {
            return z + 0.8 * (1.0/3.0 - y * y);
        }

        public _3D LevelSet => Phi;

        public double EdgeArea => 2.0;

        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant => XQuadFactoryHelper.MomentFittingVariants.Saye;
    }

    class CubicLevelSet : IEdgeQuadratureTest3D {

        public CubicLevelSet() {

        }

        double Phi(double x, double y, double z) {
            return z + 0.8 * (- y * y * y);
        }

        public _3D LevelSet => Phi;

        public double EdgeArea => 2.0;

        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant => XQuadFactoryHelper.MomentFittingVariants.Saye;
    }
}
