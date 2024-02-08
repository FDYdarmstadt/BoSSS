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
    public class LinearConvectionForm : IVolumeForm, IEdgeForm, ISpeciesFilter {

        string speciesName;

        double rho;

        string[] variableNames;

        string[] parameternames;

        int D;

        public LinearConvectionForm(string speciesName, string variableName, int D, double rho) {
            this.speciesName = speciesName;
            this.variableNames = new string[] { variableName };
            this.D = D; 
            this.rho = rho;
            this.parameternames = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D).Cat(
                BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
        }

        public TermActivationFlags VolTerms { 
            get { return TermActivationFlags.UxGradV | TermActivationFlags.GradV;}
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
            
            double r = 0.0;
            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += _uIN[0] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            if(D == 3) {
                r += _uIN[0] * inp.Parameters_IN[2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            Vector VelocityMeanIn = new Vector(D);
            for(int d = 0; d < D; ++d) {
                VelocityMeanIn[d] = inp.Parameters_IN[D + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double uJump = _uIN[0];
            r += LambdaIn * uJump;

            r *= rho;
            return r * (_vIN);
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double r = 0.0;

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += _uIN[0] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            r += _uOUT[0] * (inp.Parameters_OUT[0] * inp.Normal[0] + inp.Parameters_OUT[1] * inp.Normal[1]);
            if(D == 3) {
                r += _uIN[0] * inp.Parameters_IN[2] * inp.Normal[2] + _uOUT[0] * inp.Parameters_OUT[2] * inp.Normal[2];
            }
            
            // Calculate dissipative part
            // ==========================

            Vector VelocityMeanIn = new Vector(D);
            Vector VelocityMeanOut = new Vector(D);
            for(int d = 0; d < D; ++d) {
                VelocityMeanIn[d] = inp.Parameters_IN[D + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[D + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal);
            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = _uIN[0] - _uOUT[0];
            r += Lambda * uJump;

            r *= rho;
            r *= 0.5;
            return r * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector flux = new Vector(D);
            for(int d = 0; d < D; ++d) {
                flux[d] = U[0] * cpv.Parameters[d];
            }
            double acc = 0;
            for(int d = 0; d < D; d++)
                acc += flux[d] * GradV[d];
            acc *= rho;
            return -acc;            
        }
    }
}
