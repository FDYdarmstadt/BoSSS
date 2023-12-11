using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
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
    /// using a local Lax-Friedrichs form; Version to use in combination with <see cref="BoSSS.Solution.Control.NonLinearSolverCode.Newton"/>, see also <see cref="IDifferentialOperator.LinearizationHint"/>
    /// Here $` \vec{v} `$ is a velocity field, 
    /// $` \rho `$ is a constant density and 
    /// $` u `$ is a property which should be transported;
    /// </summary>
    /// <remarks>
    /// All boundaries are a impermeable walls for now
    /// </remarks>
    class NonLinearConvectionForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;

        double rho;

        string[] variableNames;

        string[] parameternames;

        int D;


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

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {
            //return 0.0; // solid wall

            Vector VelocityIn = new Vector(_uIN, 0, D);

// Upwinding:
            if(VelocityIn*inp.Normal >= 0) {
                return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN); // outflow
            } else {
                return 0.0 * (_vIN); // inflow
            }

            /*
            double r = 0.0;
            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            //r += _uIN[d] * (_uIN[D+0] * inp.Normal[0] + _uIN[D + 1] * inp.Normal[1]);
            //if(D == 3) {
            //    r += _uIN[d] * _uIN[D + 2] * inp.Normal[2];
            //}
            r += _uIN[D] + (VelocityIn * inp.Normal);

            // Calculate dissipative part
            // ==========================


            double LambdaIn = LambdaConvection.GetLambda(VelocityIn, inp.Normal);
            double uJump = _uIN[D];
            r += LambdaIn * uJump;

            r *= rho;
            return r * (_vIN);
            */
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

            // Upwinding:
            if(VelocityAvg*inp.Normal >= 0) {
                return rho * _uIN[D] * (VelocityIn * inp.Normal) * (_vIN - _vOUT);
            } else {
                return rho * _uOT[D] * (VelocityOt * inp.Normal) * (_vIN - _vOUT);
            }

            /*

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            //r += _uIN[d] * (_uIN[D+0] * inp.Normal[0] + _uIN[D + 1] * inp.Normal[1]);
            //r += _uOUT[d] * (_uOUT[D+0] * inp.Normal[0] + _uOUT[D+1] * inp.Normal[1]);
            //if(D == 3) {
            //    r += _uIN[d] * _uIN[D+2] * inp.Normal[2] + _uOUT[d] * _uOUT[D+2] * inp.Normal[2];
            //}
            r += _uIN[D] * (VelocityIn * inp.Normal)*0.5;
            r += _uOT[D] * (VelocityOt * inp.Normal)*0.5;

            // Calculate dissipative part
            // ==========================


            double LambdaIn = LambdaConvection.GetLambda(VelocityIn, inp.Normal);
            double LambdaOut = LambdaConvection.GetLambda(VelocityOt, inp.Normal);
            double Lambda = Math.Max(LambdaIn, LambdaOut);
            //double uJump = _uIN[D + d] - _uOUT[D + d];
            double uJump = _uIN[D] - _uOT[D];
            r += Lambda * uJump*0.5;

            r *= rho;
            return r * (_vIN - _vOUT);
            */
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;

            //Vector Velocity = new Vector(U, 0, this.D);

            for(int dim = 0; dim < D; dim++)
                acc += U[dim] * U[D] * GradV[dim];
            acc *= rho;
            return -acc;
        }
    }
}
