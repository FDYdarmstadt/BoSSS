using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {
    class ExtensionSIP : SipViscosity_GradU {

        public ExtensionSIP(double _penalty_safety, int iComp, int D, IncompressibleBoundaryCondMap bcmap, ViscosityOption _ViscosityMode, double constantViscosityValue) 
            : base(_penalty_safety, iComp, D, bcmap,_ViscosityMode, constantViscosityValue) {

        }

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }

    }

    class ExtensionPressureGradient : PressureGradientLin_d {

        int _d;

        public ExtensionPressureGradient(int _d, IncompressibleBoundaryCondMap bcmap) : base(_d, bcmap) {
            this._d = _d;
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return 0.0;
        }
    }

    class ExtensionDivergenceFlux : Divergence_DerivativeSource_Flux {
        public ExtensionDivergenceFlux(int component, IncompressibleBoundaryCondMap bcmap) : base(component, bcmap) {
            
        }

        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            FluxInCell = 0.0;
        }
    }
}
