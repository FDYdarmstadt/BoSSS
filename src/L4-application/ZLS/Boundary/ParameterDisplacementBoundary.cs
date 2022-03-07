using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class ParameterDisplacementBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public ParameterDisplacementBoundary(LevelSetTracker LsTrkr, string fluidSpecies, string solidSpecies, int d, int D, double artificialViscosity) {
            codomainName = ZwoLevelSetSolver.EquationNames.DisplacementEvolutionComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddComponent(new InterfaceConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), 1.0, d, D, 1, fluidSpecies, solidSpecies));
            AddParameter(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D));

            if(artificialViscosity != 0.0) {
                AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, d, D, 1, artificialViscosity, 0,0));
            }
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }

    class InterfaceConvectionForm : ILevelSetForm, ISupportsJacobianComponent {
        int m_iLevSet;
        double m_rho;
        string m_FluidSpc;
        string m_SolidSpecies;
        string[] variableNames;
        string[] parameterNames;
        int m_D;
        int m_d;

        public InterfaceConvectionForm(string[] variableNames, double rho, int d, int _D, int iLevSet, string FluidSpc, string SolidSpecies) {
            m_D = _D;
            m_d = d;
            m_iLevSet = iLevSet;
            m_SolidSpecies = SolidSpecies;
            m_FluidSpc = FluidSpc;
            m_rho = rho;
            this.variableNames = variableNames;
            parameterNames = ZwoLevelSetSolver.VariableNames.Displacement0Vector(_D);
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

            Vector VelocityOt = new Vector(uOut, 0, m_D);

            return m_rho * (inp.Parameters_OUT[m_d]) * (VelocityOt * inp.Normal) *(-vOut);
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
