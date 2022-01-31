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

    class FluidMomentumConvectionForm : ILevelSetForm {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variableNames;
        string[] parameterNames;
        int m_D;

        public FluidMomentumConvectionForm(string variableName, double rho, int _D, int iLevSet, string FluidSpc, string SolidSpecies) {
            m_D = _D;
            m_iLevSet = iLevSet;
            m_SolidSpecies = SolidSpecies;
            m_FluidSpc = FluidSpc;
            m_rho = rho; 
            variableNames = new string[] { variableName };
            parameterNames = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(m_D).Cat(
                BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(m_D)); ;
        }

        public IList<string> ArgumentOrdering {
            get { return variableNames;}
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

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += uIn[0] * (inp.Parameters_IN[0] * inp.Normal[0] + inp.Parameters_IN[1] * inp.Normal[1]);
            r += uOut[0] * (inp.Parameters_OUT[0] * inp.Normal[0] + inp.Parameters_OUT[1] * inp.Normal[1]);
            if(m_D == 3) {
                r += uIn[0] * (inp.Parameters_IN[2] + inp.Parameters_OUT[2]) * inp.Normal[2];
            }

            Vector VelocityMeanIn = new Vector(m_D);
            Vector VelocityMeanOut = new Vector(m_D);
            for(int d = 0; d < m_D; ++d) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_D + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_D + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal);
            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = uIn[0] - uOut[0];
            r += Lambda * uJump;

            r *= 0.5 * m_rho;
            return r * (vIn);
        }
    }
}