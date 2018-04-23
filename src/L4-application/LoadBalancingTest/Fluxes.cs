using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Application.LoadBalancingTest {

    /// <summary>
    /// fluss fuer du/dx; (Ableitung nach 1. Raumrichtung), bulk-Phase;
    /// </summary>
    class DxFlux : LinearFlux, IEquationComponentSpeciesNotification {

        LevelSetTracker m_LsTrk;
        double alpha_A;
        double alpha_B;

        public DxFlux(LevelSetTracker _LsTrk, double _alpha_A, double _alpha_B) {
            m_LsTrk = _LsTrk;
            alpha_A = _alpha_A;
            alpha_B = _alpha_B;
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
        }

        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            return -alpha * Uin[0] * inp.Normale[0];
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return -alpha * 0.5 * (Uin[0] + Uout[0]) * inp.Normale[0];
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            output[0] = -alpha * U[0];
        }

        double alpha;

        public void NowIntegratingBulk(string speciesName, SpeciesId SpcId) {
            switch(speciesName) {
                case "A": alpha = alpha_A; break;
                case "B": alpha = alpha_B; break;
                default: throw new NotImplementedException();
            }
        }
    }

    /// <summary>
    /// Flux for du/dx at the interface;
    /// </summary>
    class LevSetFlx : ILevelSetForm {

        LevelSetTracker m_LsTrk;
        double alpha_A;
        double alpha_B;

        public LevSetFlx(LevelSetTracker _LsTrk, double _alpha_A, double _alpha_B) {
            m_LsTrk = _LsTrk;
            alpha_A = _alpha_A;
            alpha_B = _alpha_B;
        }
        
        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
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

        public double LevelSetForm(ref CommonParamsLs inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double FlxNeg = -U_Neg[0] * inp.n[0] * alpha_A;
            double FlxPos = -U_Pos[0] * inp.n[0] * alpha_B;
            return (FlxNeg * vA - FlxPos * vB);
        }

        public int LevelSetIndex {
            get {
                return 0;
            }
        }

        public SpeciesId PositiveSpecies {
            get { 
                return m_LsTrk.GetSpeciesId("B"); 
            }
        }

        public SpeciesId NegativeSpecies {
            get { 
                return m_LsTrk.GetSpeciesId("A"); 
            }
        }
    }
}
