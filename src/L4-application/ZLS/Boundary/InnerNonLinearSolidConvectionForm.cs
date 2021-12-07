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

    class InnerNonLinearSolidConvectionForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variableNames;
        string[] parameterNames;
        int D;
        int d;

        public InnerNonLinearSolidConvectionForm(string[] variableNames, string[] velocityNames, double rho, int d, int iLevSet, string FluidSpc, string SolidSpecies) {
            D = velocityNames.Length;
            m_iLevSet = iLevSet;
            m_SolidSpecies = SolidSpecies;
            m_FluidSpc = FluidSpc;
            m_rho = rho;
            this.variableNames = velocityNames.Cat( variableNames);
            parameterNames = new string[] { };
            this.d = d;
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
            Vector VelocityIn = new Vector(uIn, 0, D);
            Vector VelocityOt = new Vector(uOut, 0, D);
            Vector VelocityAvg = 0.5 * (VelocityIn + VelocityOt);

            r += uOut[D + d] * (VelocityAvg * inp.Normal);

            return m_rho * r * (-vOut);

            // Upwinding:
            /*
            if (VelocityAvg * inp.Normal >= 0)
            {
                return m_rho * uOut[D+d] * (VelocityIn * inp.Normal) * ( - vOut);
            }
            else
            {
                return m_rho * uOut[D+d] * (VelocityOt * inp.Normal) * (- vOut);
            }
            */
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class NonLinearSolidConvectionForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variableNames;
        string[] parameterNames;
        int D;
        int d;

        public NonLinearSolidConvectionForm(string[] variableNames, string[] velocityNames, double rho, int d, int iLevSet, string FluidSpc, string SolidSpecies) {
            D = velocityNames.Length;
            m_iLevSet = iLevSet;
            m_SolidSpecies = SolidSpecies;
            m_FluidSpc = FluidSpc;
            m_rho = rho;
            this.variableNames = velocityNames.Cat(variableNames);
            parameterNames = new string[] { };
            this.d = d;
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
            Vector VelocityIn = new Vector(uIn, 0, D);
            Vector VelocityOt = new Vector(uOut, 0, D);
            Vector VelocityAvg = 0.5 * (VelocityIn + VelocityOt);

            double penalty = m_rho * 0.5 * Math.Abs(VelocityAvg * inp.Normal) * (uIn[D+d] - uOut[D+d]) * (vIn-vOut);
            return m_rho * 0.5 * (uOut[D + d] + uIn[D+d]) * (VelocityAvg * inp.Normal) * (vIn - vOut) + penalty;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class NonLinearSolidMomentumConvectionForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variableNames;
        string[] parameterNames;
        int m_D;

        public NonLinearSolidMomentumConvectionForm(string variableName, string[] velocityName, double rho, int _D, int iLevSet, string FluidSpc, string SolidSpecies) {
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

            Vector VelocityIn = new Vector(uIn, 1, m_D);
            Vector VelocityOt = new Vector(uOut, 1, m_D);
            Vector VelocityAvg = 0.5 * (VelocityIn + VelocityOt);

            double penalty = 0.5 * m_rho * Math.Abs(VelocityAvg * inp.Normal) * (uIn[0] - uOut[0]) * (- vOut);
            return m_rho * uIn[0] * (VelocityAvg * inp.Normal) * ( - vOut) + penalty;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}