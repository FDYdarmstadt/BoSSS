using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
 
    /// <summary>
    /// Transport term 
    /// ```math
    ///    \textrm{div} ( \rho u \vec{v} )
    /// ```
    /// using a local Lax-Friedrichs form; Version to use in combination with <see cref="BoSSS.Solution.Control.NonLinearSolverCode.Newton"/>, see also <see cref="ISpatialOperator.LinearizationHint"/>
    /// Here $` \vec{v} `$ is a velocity field, 
    /// $` \rho `$ is a constant density and 
    /// $` u `$ is a property which should be transported;
    /// </summary>
    /// <remarks>
    /// All boundaries are a impermeable walls for now
    /// </remarks>
    class NonLinearConvectionForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;

        protected double rho;

        string[] variableNames;

        string[] parameternames;

        protected int D;

        double penaltyScale = 1;

        IncompressibleMultiphaseBoundaryCondMap boundaryMap;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="speciesName"></param>
        /// <param name="variableName">
        /// variable name of transport property
        /// </param>
        /// <param name="velocity">
        /// variable names for velocity component.
        /// </param>
        /// <param name="d"></param>
        /// <param name="rho"></param>
        public NonLinearConvectionForm(string speciesName, string variableName, string[] velocity, int d, double rho) {
            this.speciesName = speciesName;
            this.variableNames = velocity.Cat(variableName);
            this.D = velocity.Length;
            //this.d = d;
            this.rho = rho;
            this.parameternames = new string[] { };
        }

        public NonLinearConvectionForm(string speciesName, string variableName, string[] velocity, int d, double rho, IncompressibleMultiphaseBoundaryCondMap boundaryMap) 
            : this(speciesName, variableName, velocity, d, rho) {
            this.boundaryMap = boundaryMap;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.UxGradV | TermActivationFlags.GradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => parameternames;

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public string ValidSpecies => speciesName;

        public virtual double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {
            //return 0.0; // solid wall

            Vector VelocityIn = new Vector(_uIN, 0, D);
            if(boundaryMap != null) {
                IncompressibleBcType edgType = boundaryMap.EdgeTag2Type[inp.EdgeTag];
                Vector vDirichlet = new Vector(D);
                for(int i = 0; i< D; ++i) {
                    vDirichlet[i] = boundaryMap.bndFunction[variableNames[i]][inp.EdgeTag](inp.X, inp.time);
                }
                double uDirichlet = boundaryMap.bndFunction[variableNames[D]][inp.EdgeTag](inp.X, inp.time);
                switch(edgType) {
                    case IncompressibleBcType.FreeSlip:
                    case IncompressibleBcType.Wall:
                        return rho * uDirichlet * (vDirichlet * inp.Normal) * (_vIN);
                    case IncompressibleBcType.Velocity_Inlet:
                    case IncompressibleBcType.Outflow:
                    case IncompressibleBcType.Pressure_Outlet:
                        return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN); // outflow
                    default:
                    throw new NotImplementedException();
                }
            } else {
                if(VelocityIn * inp.Normal >= 0) {
                    return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN); // outflow
                } else {
                    return 0.0 * (_vIN); // inflow
                }
            }
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double r = 0.0;

            Vector VelocityIn = new Vector(_uIN, 0, D);
            Vector VelocityOt = new Vector(_uOT, 0, D);
            Vector VelocityAvg = 0.5 * (VelocityIn + VelocityOt);
            //*
            double penalty = rho * 0.5 * Math.Abs(VelocityAvg * inp.Normal) * (_uIN[D] - _uOT[D]) * (_vIN - _vOUT);
            penalty *= penaltyScale;
            return rho * 0.5 * (_uIN[D] + _uOT[D]) * (VelocityAvg * inp.Normal) * (_vIN - _vOUT) + penalty;
            //*/
            // Upwinding:
            /*
            if(VelocityAvg*inp.Normal >= 0) {
                return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN - _vOUT);
            } else {
                return rho * _uOT[D] * (VelocityOt * inp.Normal) * (_vIN - _vOUT);
            }
            //*/
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;

            for(int dim = 0; dim < D; dim++)
                acc += U[dim] * U[D] * GradV[dim];
            acc *= rho;
            return -acc;
        }
    }


    class SolidNonLinearConvectionForm : NonLinearConvectionForm {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="speciesName"></param>
        /// <param name="variableName">
        /// variable name of transport property
        /// </param>
        /// <param name="velocity">
        /// variable names for velocity component.
        /// </param>
        /// <param name="d"></param>
        /// <param name="rho"></param>
        public SolidNonLinearConvectionForm(string speciesName, string variableName, string[] velocity, int d, double rho)
            : base(speciesName, variableName, velocity, d, rho) {
        }

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {
            //return 0.0; // solid wall
            Vector VelocityIn = new Vector(_uIN, 0, D);
            if(VelocityIn * inp.Normal >= 0) {
                return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN); // outflow
            } else {
                return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN);
            }
        }
    }
}
