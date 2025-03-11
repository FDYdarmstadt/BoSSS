using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {
    class ExtensionSIP : SipViscosity_GradU {

        public ExtensionSIP(double _penalty_safety, int iComp, int D, IncompressibleBoundaryCondMap bcmap, ViscosityOption _ViscosityMode, double constantViscosityValue, bool useBCMap = false) 
            : base(_penalty_safety, iComp, D, bcmap,_ViscosityMode, constantViscosityValue) {
            this.m_useBCMap = useBCMap;
        }

        bool m_useBCMap;

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            if (m_useBCMap) {
                return base.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
            } else {
                return 0;
            }
        }
    }

    class ExtensionPressureGradient : PressureGradientLin_d {

        int _d;
        bool m_useBCMap;

        public ExtensionPressureGradient(int _d, IncompressibleBoundaryCondMap bcmap, bool useBCMap = false) : base(_d, bcmap) {
            this._d = _d;
            this.m_useBCMap = useBCMap;
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            if (m_useBCMap) {
                return base.BorderEdgeFlux(ref inp, Uin);
            } else {
                return 0;
            }
        }
    }

    class ExtensionDivergenceFlux : Divergence_DerivativeSource_Flux {

        bool m_useBCMap;

        public ExtensionDivergenceFlux(int component, IncompressibleBoundaryCondMap bcmap, bool useBCMap = false) : base(component, bcmap) {
            this.m_useBCMap = useBCMap;
        }

        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            if (m_useBCMap) {
                base.BorderEdgeFlux_(ref inp, Uin, out FluxInCell);
            } else {
                FluxInCell = 0;
            }
        }
    }
}
