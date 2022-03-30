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

    //All boundaries are a Wall for now
    public class ParameterTransportForm : IVolumeForm, ISupportsJacobianComponent, IEdgeForm, ISpeciesFilter {

        string speciesName;

        double rho;

        string[] parameternames;

        int D;

        public ParameterTransportForm(string speciesName, string fieldName, string[] velocityNames, int D, double rho) {
            this.speciesName = speciesName;
            this.D = D;
            this.rho = rho;
            this.parameternames = new string[] { fieldName }.Cat(velocityNames);
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradV; }
        }

        public IList<string> ArgumentOrdering => Array.Empty<string>();

        public IList<string> ParameterOrdering => parameternames;

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.V; }
        }

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {


            return 0.0; // inflow
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double r = 0.0;

            Vector VelocityIn = new Vector(inp.Parameters_IN, 1, D);
            Vector VelocityOt = new Vector(inp.Parameters_OUT, 1, D);
            Vector VelocityAvg = 0.5 * (VelocityIn + VelocityOt);
            double a = (inp.Parameters_IN[0] - inp.Parameters_OUT[0]) * Math.Abs(VelocityAvg * inp.Normal);
            return rho * (0.5 * (inp.Parameters_IN[0] + inp.Parameters_OUT[0]) * (VelocityAvg * inp.Normal) + 0) * (_vIN - _vOUT);
            }

            return rho * ( 0.5 *  (inp.Parameters_IN[0]+ inp.Parameters_OUT[0]) * (VelocityAvg * inp.Normal) 
                +(inp.Parameters_IN[0] - inp.Parameters_OUT[0]) * 1 * (VelocityAvg * inp.Normal)) * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            Vector Up = new Vector(cpv.Parameters, 1, D);
            for(int i = 0; i < D; i++)
                acc += Up[i] * GradV[i];
            acc *= cpv.Parameters[0];
            acc *= rho;
            return -acc;
        }
    }

    class ParameterTransportBoundaryForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] parameterNames;
        int m_D;

        public ParameterTransportBoundaryForm(string parameterName, string[] velocityNames, 
            double rho, int _D, string FluidSpc, string SolidSpecies, int iLevSet) {
            m_D = _D;
            m_iLevSet = iLevSet;
            m_SolidSpecies = SolidSpecies;
            m_FluidSpc = FluidSpc;
            m_rho = rho;
            this.parameterNames = new string[] { parameterName }.Cat(velocityNames);
        }

        public IList<string> ArgumentOrdering {
            get { return Array.Empty<string>(); }
        }

        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public string PositiveSpecies {
            get { return m_SolidSpecies; }
        }

        /// <summary>
        /// Species ID of the fluid;
        /// </summary>
        public string NegativeSpecies {
            get { return m_FluidSpc; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get { return parameterNames; }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIN, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIN, double[] Grad_vOUT) {

            Vector VelocityOt = new Vector(inp.Parameters_OUT, 1, m_D);

            return m_rho * (inp.Parameters_OUT[0]) * (VelocityOt * inp.Normal) * (-vOut);
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
