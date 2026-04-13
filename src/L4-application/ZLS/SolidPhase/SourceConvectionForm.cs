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
    class SourceConvectionForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;

        protected double rho;

        string[] variableNames;

        string[] parameternames;

        protected int D;


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
        public SourceConvectionForm(string speciesName, string variableName, string[] velocity, double rho) {
            this.speciesName = speciesName;
            this.variableNames = velocity.Cat(variableName);
            this.D = velocity.Length;
            //this.d = d;
            this.rho = rho;
            this.parameternames = new string[] { };
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV | TermActivationFlags.UxV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => parameternames;

        public string ValidSpecies => speciesName;

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;

            for(int i = 0; i < D; i++)
                acc += U[i] * GradU[D, i]; 
            acc *= rho;
            return acc * V;
        }

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.None;

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return 0;
            double acc = 0;

            for(int i = 0; i < D; i++)
                acc += (_uIN[D] * _uIN[i] - _uOUT[D] * _uOUT[i]) * inp.Normal[i];
            acc *= rho;
            return -acc * 0.5 * (_vIN + _vOUT);
        }

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {
            return 0;
            double acc = 0;

            for(int i = 0; i < D; i++)
                acc += (_uIN[D] * _uIN[i]) * inp.Normal[i];
            acc *= rho;
            return -acc * (_vIN);
        }
    }
}
