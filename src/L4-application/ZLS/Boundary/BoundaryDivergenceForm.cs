using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class BoundaryDivergenceForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variables;
        int d;
        
        public BoundaryDivergenceForm(string variableName, int d, int iLevSet, string FluidSpc, string SolidSpecies) {
            this.d = d;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            variables = new string[] { variableName };
        }


        public IList<string> ArgumentOrdering {
            get {
                return variables;
            }
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
            get {
                return null;
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] Grad_pA, double[,] Grad_pB, double _vIN, double _vOUT, double[] Grad_vA, double[] Grad_vB) {
            double flux = -0.5 * (_uIN[0] - _uOUT[0]) * inp.Normal[d];
            flux *= (_vIN + _vOUT);
            return flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
