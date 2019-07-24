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
using System.Text;
using System.Threading.Tasks;

using ilPSP;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;

namespace BoSSS.Solution.XheatCommon {


    public class ConductivityAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ConductivityAtLevelSet(LevelSetTracker lstrk, double _kA, double _kB, double _penalty, double _Tsat) {
            this.kA = _kA;
            this.kB = _kB;
            this.penalty = _penalty;
            this.Tsat = _Tsat;

            m_LsTrk = lstrk;

        }

        double kA;
        double kB;

        double penalty;

        double Tsat;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            Debug.Assert(!(double.IsInfinity(hCutCellMin) || double.IsNaN(hCutCellMin)));

            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            double Ret = 0.0;


            //Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret -= (kA * Grad_uA_xN) * (vA - 0);                           // consistency term
            Ret -= (kB * Grad_uB_xN) * (0 - vB);                           // consistency term

            //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - uB[0]);                     // symmetry term
            Ret -= (kA * Grad_vA_xN) * (uA[0] - Tsat);                     // symmetry term
            Ret -= (kB * Grad_vB_xN) * (Tsat - uB[0]);                     // symmetry term

            //Ret += (penalty / hCutCellMin) * (uA[0] - uB[0]) * (vA - vB) * (Math.Abs(kA) > Math.Abs(kB) ? kA : kB); // penalty term
            Ret += (2.0*penalty / hCutCellMin) * (uA[0] - Tsat) * (vA - 0) * kA; // penalty term
            Ret += (2.0*penalty / hCutCellMin) * (Tsat - uB[0]) * (0 - vB) * kB; // penalty term



            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }



        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }


    }

    /// <summary>
    /// 
    /// </summary>
    public class HeatFluxDivergencetAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm {

        int m_D;

        LevelSetTracker m_LsTrk;

        public HeatFluxDivergencetAtLevelSet(LevelSetTracker lstrk, bool _DiriCond) {

            this.m_D = lstrk.GridDat.SpatialDimension;
            m_LsTrk = lstrk;

            //this.kAsqrt = Math.Sqrt(_kA);
            //this.kBsqrt = Math.Sqrt(_kB);

            this.DirichletCond = _DiriCond;
        }


        bool DirichletCond;

        /// <summary>
        /// 
        /// </summary>
        public double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double uAxN = GenericBlas.InnerProd(U_Neg, cp.n);
            double uBxN = GenericBlas.InnerProd(U_Pos, cp.n);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict = uBxN;
            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict = uAxN;

            double FlxNeg = (DirichletCond) ? uAxN : Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = (DirichletCond) ? uBxN : Flux(uBxN_fict, uBxN);  // flux on B-side


            return FlxNeg * vA - FlxPos * vB;

        }

        private double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in + UxN_out);
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.HeatFluxVector(m_D); }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
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

    }


    public class AuxiliaryStabilizationFormAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm {

        int m_D;

        LevelSetTracker m_LsTrk;

        public AuxiliaryStabilizationFormAtLevelSet(LevelSetTracker lstrk, bool _DiriCond) {

            this.m_D = lstrk.GridDat.SpatialDimension;
            m_LsTrk = lstrk;

            this.DirichletCond = _DiriCond;
        }

        bool DirichletCond;


        /// <summary>
        /// 
        /// </summary>
        public double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double uAxN = GenericBlas.InnerProd(U_Neg, cp.n);
            double uBxN = GenericBlas.InnerProd(U_Pos, cp.n);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict = uBxN;
            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict = uAxN;

            double FlxNeg = (DirichletCond) ? 0.0 : Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = (DirichletCond) ? 0.0 : Flux(uBxN_fict, uBxN);  // flux on B-side


            return -(FlxNeg * vA - FlxPos * vB);

        }

        static double Flux(double UxN_in, double UxN_out) {
            return (UxN_in - UxN_out);
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.HeatFluxVector(m_D); }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
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

    }


    public class TemperatureGradientAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm {

        int m_d;

        LevelSetTracker m_LsTrk;

        public TemperatureGradientAtLevelSet(int _d, LevelSetTracker lstrk, double _kA, double _kB, bool _DiriCond, double _Tsat) {

            this.m_d = _d;
            m_LsTrk = lstrk;

            this.kA = _kA;
            this.kB = _kB;

            this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;
        }

        double kA;
        double kB;

        bool DirichletCond;
        double Tsat;

        /// <summary>
        /// 
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double Acc = (DirichletCond) ? Tsat : 0.5 * (uB[0] + uA[0]);

            return -Acc * (kB * vB - kA * vA) * inp.n[m_d];
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }


    public class TemperatureStabilizationFormAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm {

        int m_d;

        LevelSetTracker m_LsTrk;

        public TemperatureStabilizationFormAtLevelSet(int _d, LevelSetTracker lstrk, bool _DiriCond, double _Tsat) {

            this.m_d = _d;
            m_LsTrk = lstrk;

            this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;
        }

        bool DirichletCond;
        double Tsat;

        /// <summary>
        /// 
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            //return (uA[0] - uB[0]) * inp.n[m_d] * (vA - vB);

            double Acc = 0.0;

            if (DirichletCond) {
                Acc += 2.0 * (uA[0] - Tsat) * inp.n[m_d] * (vA - 0.0);
                Acc += 2.0 * (Tsat - uB[0]) * inp.n[m_d] * (0.0 - vB);
            } else {
                Acc += (uA[0] - uB[0]) * inp.n[m_d] * (vA - vB);
            }

            return Acc;
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }

}
