using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IBM_Solver {
 
   
    public class ConvectionAtIB : ILevelSetForm, ISupportsJacobianComponent {
        public ConvectionAtIB(LevelSetTracker LsTrk, int _d, int _D, double fluidDensity, bool UseMovingMesh) {
            m_LsTrk = LsTrk;
            m_D = _D;
            m_d = _d;
            fDensity = fluidDensity;
            m_UseMovingMesh = UseMovingMesh;
        }

        int m_D;
        int m_d;
        double fDensity;
        bool m_UseMovingMesh;
        LevelSetTracker m_LsTrk;

        // Use Fluxes as in Bulk Convection
        LinearizedConvection NegFlux;


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
           
        public double InnerEdgeForm(ref CommonParams cp, 
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(m_D == cp.D);

            if (m_UseMovingMesh) {
                return 0.0;
            } else {
                Vector Uin = new Vector(U_Neg);
                return (Uin * cp.Normal) * Uin[m_d] * fDensity * v_Neg;
            }

        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new NotImplemntedClass(this) };
        }

        class NotImplemntedClass : ILevelSetForm {
            public NotImplemntedClass(ConvectionAtIB __owner) {
                m_owner = __owner;
            }
            ConvectionAtIB m_owner;

            public int LevelSetIndex => m_owner.LevelSetIndex;

            public SpeciesId PositiveSpecies => m_owner.PositiveSpecies;

            public SpeciesId NegativeSpecies => m_owner.NegativeSpecies;

            public TermActivationFlags LevelSetTerms => m_owner.LevelSetTerms;

            public IList<string> ArgumentOrdering => m_owner.ArgumentOrdering;

            public IList<string> ParameterOrdering => m_owner.ParameterOrdering;

            public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
                throw new NotImplementedException();
            }
        }


    }

}
