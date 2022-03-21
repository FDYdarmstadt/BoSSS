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
        double scale;

        public BoundaryDivergenceDisplacementForm(int d, int iLevSet, string FluidSpc, string SolidSpecies, double scale)
        {
            this.d = d;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            variables = new string[] { BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d), VariableNames.DisplacementComponent(d)};
            this.scale = scale;
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
            //double flux = -_uOUT[1] *  _vOUT * inp.Normal[d];
            double flux = (_uIN[0]- _uOUT[0]) * _vIN * inp.Normal[d];
            //flux += scale *(_uIN[0]- _uOUT[0]) * -_vOUT * inp.Normal[d];
            return flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new IEquationComponent[] { this };
        }
    }

    class ConvectionDivergenceBoundaryForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variables;
        int D;
        double scale;

        public ConvectionDivergenceBoundaryForm(int D, int iLevSet, string FluidSpc, string SolidSpecies) {
            this.D = D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            variables = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D).Cat(VariableNames.DisplacementVector(D));
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
                return TermActivationFlags.GradUxV|TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] Grad_pA, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] Grad_vA, double[] Grad_vB) {
            double acc = 0;
            for(int i = 0; i < D; i++) {
                double t = 0;
                for(int j = 0; j < D; ++j) {
                    t = (_uOUT[j]) * ( _Grad_uOUT[D + i, j]);
                }
                acc += t * inp.Normal[i]; ;
            }
            return -acc * (- _vOUT);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}
