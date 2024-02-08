using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace XESF.Fluxes {
    public class CentralFlux_XDG_Interface : ILevelSetForm, IEquationComponentChecking, ISupportsJacobianComponent {
        #region ILevelSetForm members
        public int LevelSetIndex {
            get;
            private set;
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return CompressibleEnvironment.PrimalArgumentOrdering;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
        #endregion
        public bool IgnoreVectorizedImplementation => false;

        public CentralFlux_XDG_Interface(IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int component = int.MinValue, int levelSetIndex = 0, string posSpecies = "R", string negSpecies = "L") {
            //this.levelSetTracker = levelSetTracker;

            switch(flux) {
                case FluxVariables.Density:
                this.bulkFlux = new CentralDensityFlux(boundaryConditionMap, material);
                break;
                case FluxVariables.Momentum:
                this.bulkFlux = new CentralMomentumFlux(boundaryConditionMap, material, component);
                break;
                case FluxVariables.Energy:
                this.bulkFlux = new CentralEnergyFlux(boundaryConditionMap, material);
                break;
                default:
                throw new NotImplementedException();
            }

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");
            return new IEquationComponent[] {
                new LevelSetFormDifferentiator(this,SpatialDimension)
            };
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return this.bulkFlux.InnerEdgeForm(ref inp, _uIN, _uOUT, _Grad_uIN, _Grad_uOUT, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
        }

        protected readonly CentralFlux bulkFlux;

        
    }
}
