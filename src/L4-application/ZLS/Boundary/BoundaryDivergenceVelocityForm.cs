using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class BoundaryDivergenceVelocityForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variables;
        int d;
        
        public BoundaryDivergenceVelocityForm(int d, int iLevSet, string FluidSpc, string SolidSpecies) {
            this.d = d;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            variables = new string[] { BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d)};
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
            double flux = (_uIN[0] - _uOUT[0]) *0.5* (_vIN + _vOUT ) * inp.Normal[d];
            return flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class BoundaryDivergenceDisplacementForm : ILevelSetForm, ISupportsJacobianComponent
    {
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variables;
        int d;

        public BoundaryDivergenceDisplacementForm(int d, int iLevSet, string FluidSpc, string SolidSpecies)
        {
            this.d = d;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            variables = new string[] { BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d), VariableNames.DisplacementComponent(d)};
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

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] Grad_pA, double[,] Grad_pB, double _vIN, double _vOUT, double[] Grad_vA, double[] Grad_vB)
        {
            double flux = -_uOUT[1] *  _vOUT * inp.Normal[d];
            flux += _uIN[0] * _vIN * inp.Normal[d];
            return flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new IEquationComponent[] { this };
        }
    }
}
