using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using System.Security.Cryptography.X509Certificates;

namespace XESTSF.SPFluxes {
    public static class SPFluxConverter {

        //maybe we only need to adjust the dimension? efficency
        public static CommonParams GetRedInp(CommonParams inp) {
            var inpDM1 = new CommonParams();
            inpDM1.time = inp.X[inp.D - 1];
            inpDM1.X = new double[inp.D - 1];
            inpDM1.Normal = new double[inp.D - 1];
            for(int i = 0; i < inp.D - 1; i++) {
                inpDM1.X[i] = inp.X[i];
                inpDM1.Normal[i] = inp.Normal[i];
            }
            inpDM1.EdgeTag = inp.EdgeTag;
            inpDM1.iEdge = inp.iEdge;
            inpDM1.jCellIn = inp.jCellIn;
            inpDM1.jCellOut = inp.jCellOut;
            inpDM1.Parameters_OUT = inp.Parameters_OUT;
            inpDM1.GridDat = inp.GridDat;
            inpDM1.Parameters_IN = inp.Parameters_IN;

            return inpDM1;
        }
        public static CommonParamsVol GetRedInp(CommonParamsVol inp) {
            var inpDM1 = new CommonParamsVol();
            inpDM1.time = inp.Xglobal[inp.D - 1];
            inpDM1.Xglobal = new double[inp.D - 1];
            for(int i = 0; i < inp.D - 1; i++) {
                inpDM1.Xglobal[i] = inp.Xglobal[i];
            }
            inpDM1.Parameters = inp.Parameters;
            inpDM1.GridDat = inp.GridDat;
            inpDM1.jCell=inp.jCell;

            return inpDM1;
        }
        public static CommonParamsBnd GetRedInp(CommonParamsBnd inp) {
            var inpDM1 = new CommonParamsBnd();
            inpDM1.time = inp.X[inp.D - 1];
            inpDM1.X = new double[inp.D - 1];
            inpDM1.Normal = new double[inp.D - 1];
            for(int i = 0; i < inp.D - 1; i++) {
                inpDM1.X[i] = inp.X[i];
                inpDM1.Normal[i] = inp.Normal[i];
            }
            inpDM1.EdgeTag = inp.EdgeTag;
            inpDM1.iEdge = inp.iEdge;
            inpDM1.GridDat = inp.GridDat;
            inpDM1.Parameters_IN = inp.Parameters_IN;

            return inpDM1;
        }

    }
    public class SPBulkFlux : IEdgeForm, IVolumeForm, ISupportsJacobianComponent,ISpeciesFilter {

        public IEdgeForm m_EdgeSpaceFlux;
        public IVolumeForm m_VolumeSpaceFlux;
        public int i_comp;
        public IBoundaryConditionMap boundaryMap;
        public SPBulkFlux(IBoundaryConditionMap boundaryMap,IEdgeForm EdgeSpaceFlux,IVolumeForm VolumeSpaceFlux, int i_comp,SpeciesId speciesId, string spcName,Material material) {
            m_EdgeSpaceFlux = EdgeSpaceFlux;
            m_VolumeSpaceFlux = VolumeSpaceFlux;
            this.i_comp = i_comp;
            this.speciesId = speciesId;
            this.ValidSpecies = spcName;
            this.material = material;
            this.boundaryMap = boundaryMap;
        }
        public TermActivationFlags BoundaryEdgeTerms => m_EdgeSpaceFlux.BoundaryEdgeTerms;

        public TermActivationFlags InnerEdgeTerms => m_EdgeSpaceFlux.InnerEdgeTerms;

        public IList<string> ArgumentOrdering => m_EdgeSpaceFlux.ArgumentOrdering;

        public IList<string> ParameterOrdering => m_EdgeSpaceFlux.ParameterOrdering;

        public TermActivationFlags VolTerms => m_VolumeSpaceFlux.VolTerms;

        public string ValidSpecies {
            get;
            private set;
        }

        private SpeciesId speciesId;
        protected readonly Material material;
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
                boundaryCondition = xdgBoudaryMap.GetBoundaryConditionForSpecies(EdgeTag, this.ValidSpecies);
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
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            //TODO: Think about wether this makes truly sense.....
            var inpDM1 = SPFluxConverter.GetRedInp(inp);

            double[] Uout = GetExternalState(ref inpDM1, _uA);
            //Create a D-1 inp object
            var SpaceFlux = m_EdgeSpaceFlux.BoundaryEdgeForm(ref inpDM1, _uA, _Grad_uA, _vA, _Grad_vA);
            var nt = inp.Normal[inp.D - 1];

            if(nt >= 0) {
                return SpaceFlux + nt * _uA[i_comp] * _vA;
            } else {
                return SpaceFlux + nt * Uout[i_comp] *_vA; 
            }

        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");

            return new IEquationComponent[] {
                new EdgeFormDifferentiator(this, SpatialDimension),
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }
        //public double[] ExtractNTArray(double[] inp) {
        //    var output= new double[inp.Length - 1];
        //    for(int i=0;i<inp.Length-1;i++) { 

        //    }
        //}

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            //Create a D-1 inp object
            var inpDM1 = SPFluxConverter.GetRedInp(inp);
            var SpaceFlux = m_EdgeSpaceFlux.InnerEdgeForm(ref inpDM1, _uIN, _uOUT,_Grad_uIN, _Grad_uOUT, _vIN,_vOUT, _Grad_vIN, _Grad_vOUT);
            var nt = inp.Normal[inp.D-1];
            double timeflux = 0;

            if(nt >= 0) {
                timeflux = nt * _uIN[i_comp];
            } else {
                timeflux = nt * _uOUT[i_comp];
            }
            timeflux = nt * 0.5 * (_uIN[i_comp] + _uOUT[i_comp]);
            return SpaceFlux + timeflux * (_vIN - _vOUT);

        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            var cpvDM1 = SPFluxConverter.GetRedInp(cpv);
            var GradVDM1 = new double[GradV.Length - 1];
            for (int i = 0; i < GradVDM1.Length; i++) {
                GradVDM1[i] = GradV[i];
            }
            var SpaceFlux = m_VolumeSpaceFlux.VolumeForm(ref cpvDM1, U, GradU, V, GradVDM1);


            return SpaceFlux - U[i_comp] * GradV[cpv.D-1];
        }
    }
    public class SPInterfaceFlux : ILevelSetForm, ISupportsJacobianComponent {

        public ILevelSetForm m_InterfaceSpaceFlux;
        public int i_comp;
        public SPInterfaceFlux(ILevelSetForm InterfaceSpaceFlux, int i_comp) {
            m_InterfaceSpaceFlux = InterfaceSpaceFlux;
            this.i_comp = i_comp;
        }
        public int LevelSetIndex => m_InterfaceSpaceFlux.LevelSetIndex;

        public string PositiveSpecies => m_InterfaceSpaceFlux.PositiveSpecies;

        public string NegativeSpecies => m_InterfaceSpaceFlux.NegativeSpecies;

        public TermActivationFlags LevelSetTerms => m_InterfaceSpaceFlux.LevelSetTerms;

        public IList<string> ArgumentOrdering => m_InterfaceSpaceFlux.ArgumentOrdering;

        public IList<string> ParameterOrdering => m_InterfaceSpaceFlux.ParameterOrdering;

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");
            return new IEquationComponent[] {
                new LevelSetFormDifferentiator(this,SpatialDimension)
            };
        }

        
        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            var inpDM1 = SPFluxConverter.GetRedInp(inp);
            var SpaceFlux = m_InterfaceSpaceFlux.InnerEdgeForm(ref inpDM1, _uIN,_uOUT, _Grad_uIN, _Grad_uOUT, _vIN,_vOUT,_Grad_vIN,_Grad_vOUT);
            var nt = inp.Normal[inp.D - 1];
            double timeflux = 0;
            if(nt >= 0) {
                timeflux = nt * _uIN[i_comp];
            } else {
                timeflux = nt * _uOUT[i_comp];
            }
            timeflux = nt * 0.5 * (_uIN[i_comp] + _uOUT[i_comp]);
            return SpaceFlux + timeflux * (_vIN - _vOUT);
        }
    }
}
