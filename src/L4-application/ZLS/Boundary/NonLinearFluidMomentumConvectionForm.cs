using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {

    class NonLinearFluidMomentumConvectionForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variableNames;
        string[] parameterNames;
        int m_D;

        public NonLinearFluidMomentumConvectionForm(string variableName, string[] velocityName, double rho, int _D, int iLevSet, string FluidSpc, string SolidSpecies) {
            m_D = _D;
            m_iLevSet = iLevSet;
            m_SolidSpecies = SolidSpecies;
            m_FluidSpc = FluidSpc;
            m_rho = rho;
            variableNames = new string[] { variableName }.Cat(velocityName);
            parameterNames = new string[] { };
        }

        public IList<string> ArgumentOrdering {
            get { return variableNames; }
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
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return parameterNames; }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIN, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIN, double[] Grad_vOUT) {
            double r = 0.0;
            Vector VelocityOt = new Vector(uOut, 1, m_D);
            Vector VelocityIn = new Vector(uIn, 1, m_D);
            return m_rho * uIn[0] * (VelocityIn * inp.Normal) * (vIn);
            /*
            if(VelocityOt * inp.Normal >= 0) {
                return m_rho * uIn[0] * (VelocityOt * inp.Normal) * (vIn);
            } else {
                return m_rho * uOut[0] * (VelocityOt * inp.Normal) * (vIn);
            }
            */
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}