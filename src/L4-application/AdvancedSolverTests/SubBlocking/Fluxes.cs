/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;

namespace AdvancedSolverTests {

    /// <summary>
    /// This flux is formally wrong. Its ment to get rid of secondary diagonals, but filling main diagonal of LES
    /// </summary>
    class SourceTest : LinearSource {

        public SourceTest(string varname, double factor) {
            m_varname = varname;
            m_factor = factor;
        }

        string m_varname;
        double m_factor;

        public override IList<string> ArgumentOrdering {
            get { return new string[] { m_varname }; }
        }

        protected override double Source(double[] x, double[] parameters, double[] U) {
            //return U[0] -(x[0] * x[0]+ x[1] * x[1]- 2 *x[0]*x[1]+3*x[0]+ 4*x[1] - 5);
            return U[0] - 5;
        }


    }

    /// <summary>
    /// fluss fuer du/dx; (Ableitung nach 1. Raumrichtung), bulk-Phase;
    /// </summary>
    class DxFlux : LinearFlux {

        public DxFlux(string varname, double factor) {
            m_varname = varname;
            m_factor = factor;
        }

        string m_varname;
        double m_factor;

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { m_varname };
            }
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return Uin[0] * inp.Normal[0] * m_factor;
        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return 0.5 * (Uin[0] + Uout[0]) * inp.Normal[0] * m_factor;
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            output[0] = U[0] * m_factor;
        }
    }

    /// <summary>
    /// Fluss fuer du/dx; (Ableitung nach 1. Raumrichtung), common parts for both level-sets;
    /// </summary>
    class LevSetFlx : ILevelSetForm {

        //protected LevelSetTracker m_LsTrk;

        public LevSetFlx(string varname, double factor) {
            //m_LsTrk = _LsTrk;
            m_varname = varname;
            m_factor = factor;
        }

        double m_factor;
        string m_varname;

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { m_varname };
            }
        }

        //public TermActivationFlags LevelSetTerms {
        //    get {
        //        return TermActivationFlags.UxV;
        //    }
        //}

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double Flx = 0.5 * (U_Pos[0] + U_Neg[0]) * inp.Normal[0];
            return Flx * vA - Flx * vB * m_factor;
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get {
                return "A";
            }
        }

        public string NegativeSpecies {
            get {
                return "B";
            }
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.None; }
        }
    }


    public class XLaplaceBCs {

        /// <summary>
        /// Function which determines which part of the domain boundary is of Dirichlet type
        /// </summary>
        public Func<CommonParamsBnd, bool> IsDirichlet;

        /// <summary>
        /// Dirichlet boundary value
        /// </summary>
        public Func<CommonParamsBnd, double> g_Diri;

        /// <summary>
        /// Neumann boundary value
        /// </summary>
        public Func<CommonParamsBnd, double> g_Neum;
    }

    /// <summary>
    /// Laplace operator in the bulk.
    /// </summary>
    public class XLaplace_Bulk : BoSSS.Solution.NSECommon.SIPLaplace, IEquationComponentSpeciesNotification, IEquationComponentCoefficient {

        public XLaplace_Bulk(double _muA, double _muB, double __penatly_baseFactor, string n)
            : base(__penatly_baseFactor, n) {
            muA = _muA;
            muB = _muB;
            base.m_alpha = 1.0;
            var BC = new XLaplaceBCs();
            BC.g_Diri = ((CommonParamsBnd inp) => 0.0);
            BC.IsDirichlet = (inp => true);
            this.boundaries = BC;
            this.m_Mode = XLaplace_Interface.Mode.SIP;
        }

        double muA;
        double muB;

        XLaplace_Interface.Mode m_Mode;


        public void SetParameter(string speciesName) {
            switch (speciesName) {
                case "A": species_Mu = muA; otherSpecies_Mu = muB; break;
                case "B": species_Mu = muB; otherSpecies_Mu = muA; break;
                default: throw new ArgumentException("Unknown species.");
            }

        }

        double species_Mu;
        double otherSpecies_Mu;

        public override double Nu(double[] x, double[] parameters, int jCell) {
            return this.species_Mu;
        }

        XLaplaceBCs boundaries;

        protected override double g_Diri(ref CommonParamsBnd inp) {
            return boundaries.g_Diri(inp);
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return boundaries.IsDirichlet(inp);
        }
        protected override double g_Neum(ref CommonParamsBnd inp) {
            return boundaries.g_Neum(inp);
        }

        double GetPenalty(ref CommonParams inp) {
            double mu = base.GetPenalty(inp.jCellIn, inp.jCellOut);

            double penalty_muFactor;
            switch (this.m_Mode) {
                case XLaplace_Interface.Mode.SIP:
                //case XLaplace_Interface.Mode.VolumeScaled:
                case XLaplace_Interface.Mode.SWIP:
                    penalty_muFactor = species_Mu;
                    break;

                default:
                    throw new NotImplementedException();
            }

            return mu * penalty_muFactor;
        }



        override public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;


            double muA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);
            double muB = this.Nu(inp.X, inp.Parameters_OUT, inp.jCellOut);

            double scaleIN, scaleOT;
            ComputeScaling(ref inp, out scaleIN, out scaleOT);


            for (int d = 0; d < inp.D; d++) {
                Acc += (scaleIN * muA * _Grad_uA[0, d] + scaleOT * muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term
                Acc += (scaleIN * muA * _Grad_vA[d] + scaleOT * muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];  // symmetry term
            }
            Acc *= this.m_alpha;


            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * this.GetPenalty(ref inp); // penalty term

            return Acc;
        }

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            return base.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
        }


        void ComputeScaling(ref CommonParams inp, out double scaleIN, out double scaleOT) {
            switch (this.m_Mode) {
                case XLaplace_Interface.Mode.SIP:
                case XLaplace_Interface.Mode.SWIP: {
                    scaleIN = 0.5;
                    scaleOT = 0.5;
                    return;
                }

                default:
                    throw new NotImplementedException();
            }
        }

        override public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);
            //this.m_LenScales = cs.CellLengthScales;
        }
    }

    /// <summary>
    /// Laplace operator at the interface
    /// </summary>
    public class XLaplace_Interface : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        public enum Mode {
            SWIP,
            SIP
        }

        //protected LevelSetTracker m_LsTrk;

        public XLaplace_Interface(double _muA, double _muB, double __penatly_baseFactor, string varname) {
            //this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penatly_baseFactor = __penatly_baseFactor;
            this.m_mode = XLaplace_Interface.Mode.SIP;
            this.m_varname = varname;

        }



        protected double muA;
        protected double muB;
        protected double penatly_baseFactor;
        protected Mode m_mode;
        string m_varname;


        public virtual double InnerEdgeForm(ref CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;
            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            int D = N.Length;
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double omega_A, omega_B;
            ComputeScaling(ref inp, out omega_A, out omega_B);

            double Ret = 0.0;
            Ret += (muA * omega_A * Grad_uA_xN + muB * omega_B * Grad_uB_xN) * (vA - vB);
            Ret += (muA * omega_A * Grad_vA_xN + muB * omega_B * Grad_vB_xN) * (uA[0] - uB[0]);

            //
            Ret -= GetPenalty(ref inp) * (uA[0] - uB[0]) * (vA - vB);

            return Ret;
        }

        MultidimensionalArray NegLengthScaleS, PosLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            double _p = DomainDGdeg.Max();
            double _D = csA.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_deg = Math.Max(penalty_deg_tri, penalty_deg_sqr);

        }

        /// <summary>
        /// penalty degree multiplier
        /// </summary>
        double m_penalty_deg;

        protected double GetPenalty(ref CommonParams inp) {

            double PosCellLengthScale = PosLengthScaleS[inp.jCellOut];
            double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

            double penaltySizeFactor_A = 1.0 / NegCellLengthScale;
            double penaltySizeFactor_B = 1.0 / PosCellLengthScale;
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);


            double penalty_muFactor;
            switch (this.m_mode) {
                case Mode.SWIP:
                    penalty_muFactor = (2.0 * muA * muB) / (muA + muB);
                    break;


                case Mode.SIP:
                    //case Mode.VolumeScaled:
                    penalty_muFactor = Math.Max(Math.Abs(muA), Math.Abs(muB)) * Math.Sign(muA);
                    break;

                default:
                    throw new NotImplementedException();
            }

            double mu = this.penatly_baseFactor * penaltySizeFactor * penalty_muFactor * m_penalty_deg;
            if (mu.IsNaNorInf())
                throw new ArithmeticException("NAN/INF in penalty param");
            return mu;
        }

        //List<int> cellElo = new List<int>();


        protected void ComputeScaling(ref CommonParams inp, out double scaleIN, out double scaleOT) {
            Debug.Assert(Math.Sign(muA) == Math.Sign(muB));

            switch (this.m_mode) {
                case Mode.SWIP: {
                    scaleIN = muB / (muA + muB);
                    scaleOT = muA / (muA + muB);
                    return;
                }

                case Mode.SIP: {
                    // Konventionell:
                    scaleIN = 0.5;
                    scaleOT = 0.5;
                    return;
                }

                default:
                    throw new NotImplementedException();
            }
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { m_varname }; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.None; }
        }

    }

}
