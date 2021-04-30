using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IntersectingLevelSetTest {
    class LevelSetJumpFlux  : LinearFlux {

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" };
            }
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return Uin[0] * inp.Normal[0];
        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return 0.5 * (Uin[0] + Uout[0]) * inp.Normal[0];
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            output[0] = U[0];
        }
    }

    class LevSetJump_AB : LevSetFlx {

        public LevSetJump_AB() : base() { }

        public override double InnerEdgeForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double FlxNeg = U_Neg[0] * inp.Normal[0];
            double FlxPos = U_Pos[0] * inp.Normal[0];

            return (FlxNeg + FlxPos) / 2 * (vA - vB);

        }

        public override int LevelSetIndex {
            get { return 0; }
        }

        public override string PositiveSpecies {
            get {
                return "A";
            }
        }

        public override string NegativeSpecies {
            get {
                return "B";
            }
        }
    }

    class LevSetJump_CB : LevSetFlx {
        public LevSetJump_CB() : base() { }

        public override double InnerEdgeForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double FlxNeg = U_Neg[0] * inp.Normal[0];
            double FlxPos = U_Pos[0] * inp.Normal[0];

            return (FlxNeg + FlxPos) / 2 * (vA - vB);
        }

        public override int LevelSetIndex {
            get { return 1; }
        }

        public override string PositiveSpecies {
            get { return "C"; }
        }

        public override string NegativeSpecies {
            get { return "B"; }
        }
    }

    class LevSetJump_CA : LevSetFlx {
        public LevSetJump_CA() : base() { }

        public override double InnerEdgeForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double FlxNeg = U_Neg[0] * inp.Normal[0];
            double FlxPos = U_Pos[0] * inp.Normal[0];

            return (FlxNeg + FlxPos) / 2 * (vA - vB);
        }

        public override int LevelSetIndex {
            get { return 1; }
        }

        public override string PositiveSpecies {
            get { return "C"; }
        }

        public override string NegativeSpecies {
            get { return "A"; }
        }
    }


}
