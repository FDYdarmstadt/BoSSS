using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;

namespace XESF.Fluxes {
    public abstract class RHFlux : IEdgeForm, ISupportsJacobianComponent {

        /// <summary>
        /// <see cref="HLLCDensityFlux.HLLCDensityFlux"/>
        /// </summary>
        protected readonly IBoundaryConditionMap boundaryMap;

        protected readonly Material material;

        protected string speciesName;

        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        public RHFlux(IBoundaryConditionMap boundaryMap, Material material) {
            this.boundaryMap = boundaryMap;
            this.material = material;
        }
        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
        public TermActivationFlags InnerEdgeTerms {
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

        public double[] GetExternalState(ref CommonParamsBnd inp, double[] Uin) {
            double[] Uout = new double[Uin.Length];

            var xdgBoudaryMap = this.boundaryMap as XDGCompressibleBoundaryCondMap;
            var boundaryMap = this.boundaryMap as CompressibleBoundaryCondMap;
            if(xdgBoudaryMap == null && boundaryMap == null)
                throw new NotSupportedException("This type of boundary condition map is not supported.");

            //get EdgeTag
            byte EdgeTag = inp.EdgeTag;

            // Get boundary condition on this edge
            BoundaryCondition boundaryCondition;
            if(xdgBoudaryMap != null)
                boundaryCondition = xdgBoudaryMap.GetBoundaryConditionForSpecies(EdgeTag, this.speciesName);
            else
                boundaryCondition = boundaryMap.GetBoundaryCondition(EdgeTag);

            //create state vector from Uin
            Vector momentum = new Vector(inp.D);
            for(int i = 0; i < inp.D; i++) {
                momentum[i] = Uin[i + 1];
            }
            var stateIn = new StateVector(material, Uin[0], momentum, Uin[inp.D + 1]);

            //Call state evaluation
            var stateOut = boundaryCondition.GetBoundaryState(inp.time, inp.X, inp.Normal, stateIn);

            //Creat Uout from stateOut
            Uout[0] = stateOut.Density;
            for(int i = 0; i < inp.D; i++) {
                Uout[1 + i] = stateOut.Momentum[i];
            }
            Uout[inp.D + 1] = stateOut.Energy;

            return Uout;
        }
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {


            double[] Uout = GetExternalState(ref inp, Uin);
            double[,] _Grad_Uout = _Grad_Uin.CloneAs();
            double Vout = 0;
            double[] _Grad_Vout = new double[_Grad_Vin.Length];

            var inpIn = new CommonParams();
            inpIn.Normal = inp.Normal;
            inpIn.X = inp.X;
            inpIn.time = inp.time;
            inpIn.iEdge = inp.iEdge;
            inpIn.GridDat = inp.GridDat;
            inpIn.EdgeTag = 0;
            



            return InnerEdgeForm(ref inpIn, Uin, Uout, _Grad_Uin, _Grad_Uout, Vin,Vout, _Grad_Vin, _Grad_Vout); //<- is this even possible somehow?
            
        }
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");

            return new IEquationComponent[] {
                new EdgeFormDifferentiator(this, SpatialDimension)
            };
        }

        public abstract Vector GetFluxVector(int D, double[] U);

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            //double[] U = new double[_uIN.Length];
            //for(int i = 0; i < U.Length; i++) {
            //    U[i] = 0.5*(_uIN[i] + _uOUT[i]);
            //}
            //double[] GradU = new double[_uIN.Length];
            //double[] V = new double[_uIN.Length];
            //double[] GradV = new double[_uIN.Length];

            int D = inp.D;
            return (GetFluxVector(D, _uIN) * inp.Normal - GetFluxVector(D, _uOUT) * inp.Normal) * (_vIN - _vOUT);
            //return GetFluxVector(D, U) * inp.Normal * (_vIN - _vOUT);
        }
    }
}
