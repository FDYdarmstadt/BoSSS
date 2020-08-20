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
using System.Collections;

namespace BoSSS.Solution.XheatCommon {


    public class ConductivityAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ConductivityAtLevelSet(LevelSetTracker lstrk, double _kA, double _kB, double _penalty, double _Tsat) {
            this.kA = _kA;
            this.kB = _kB;
            this.m_penalty_base = _penalty;
            this.m_D = lstrk.GridDat.SpatialDimension;

            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;

            m_LsTrk = lstrk;

        }

        double kA;
        double kB;

        double penalty;
        int m_D;

        //bool DirichletCond;
        double Tsat;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(inp.jCellIn == inp.jCellOut);

            double[] N = inp.Normal;
            //double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];

            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == m_D);
            Debug.Assert(Grad_uB.GetLength(1) == m_D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            //double PosCellLengthScale = PosLengthScaleS[inp.jCellOut];
            //double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

            //double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            //Debug.Assert(!(double.IsInfinity(hCutCellMin) || double.IsNaN(hCutCellMin)));

            //if (hCutCellMin <= 1.0e-10 * hCellMin)
            //    // very small cell -- clippling
            //    hCutCellMin = hCellMin;

            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);
            double wPenalty = (Math.Abs(kA) > Math.Abs(kB)) ? kA : kB;

            double Ret = 0.0;

            if (!evapMicroRegion[inp.jCellIn]) {
                Ret -= (kA * Grad_uA_xN) * (vA - 0);                           // consistency term
                Ret -= (kB * Grad_uB_xN) * (0 - vB);                           // consistency term
                Ret -= (kA * Grad_vA_xN) * (uA[0] - Tsat);                     // symmetry term
                Ret -= (kB * Grad_vB_xN) * (Tsat - uB[0]);                     // symmetry term

                //Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - vB);
                //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - Tsat);
                //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (Tsat - uB[0]);

                Ret += (uA[0] - Tsat) * (vA - 0) * 2.0 * pnlty * kA; // penalty term
                Ret += (Tsat - uB[0]) * (0 - vB) * 2.0 * pnlty * kB; // penalty term

                //Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - vB);                           // consistency term
                //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - uB[0]);                     // symmetry term
                //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * Tsat;

                //Ret += (uA[0] - uB[0]) * (vA - vB) * pnlty * wPenalty; // penalty term
                //Ret -= Tsat * (vA - vB) * pnlty * wPenalty;

            } else {
                //Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - vB);                           // consistency term
                //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - uB[0]);                     // symmetry term

                //Ret += (penalty / hCutCellMin) * (uA[0] - uB[0]) * (vA - vB) * (Math.Abs(kA) > Math.Abs(kB) ? kA : kB); // penalty term
            }


            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        //MultidimensionalArray PosLengthScaleS;
        //MultidimensionalArray NegLengthScaleS;

        //public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

        //    NegLengthScaleS = csA.CellLengthScales;
        //    PosLengthScaleS = csB.CellLengthScales;

        //    if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
        //        evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        //}

        BitArray evapMicroRegion;


        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double Penalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];
            double penaltySizeFactor_B = 1.0 / PosLengthScaleS[jCellOut];

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        /// <param name="csA"></param>
        /// <param name="csB"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];
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
                return TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }


    public class HeatFluxAtLevelSet : EvaporationAtLevelSet {

        public HeatFluxAtLevelSet(int _D, LevelSetTracker _LsTrk, ThermalParameters thermParams, double _sigma) 
            : base(_D, _LsTrk, thermParams, _sigma) {

        }


        public override double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, 
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(cp.jCellIn == cp.jCellOut);

            if (!evapMicroRegion[cp.jCellIn]) {

                return 0.0;

            } else {

                double q = ComputeHeatFlux(cp.Parameters_IN, cp.Parameters_OUT, cp.Normal, cp.jCellIn);

                double FlxNeg = -0.5 * q;
                double FlxPos = +0.5 * q;
                
                double Ret = FlxNeg * vA - FlxPos * vB;

                return Ret;

            }


        }

    }


    /// <summary>
    /// 
    /// </summary>
    public class HeatFluxDivergencetAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        int m_D;

        LevelSetTracker m_LsTrk;

        public HeatFluxDivergencetAtLevelSet(LevelSetTracker lstrk) {

            this.m_D = lstrk.GridDat.SpatialDimension;
            m_LsTrk = lstrk;

            //this.kAsqrt = Math.Sqrt(_kA);
            //this.kBsqrt = Math.Sqrt(_kB);

           // this.DirichletCond = _DiriCond;
        }


        //bool DirichletCond;

        /// <summary>
        /// 
        /// </summary>
        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);
            double uBxN = GenericBlas.InnerProd(U_Pos, cp.Normal);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict = uBxN;
            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict = uAxN;

            double FlxNeg = (!evapMicroRegion[cp.jCellIn]) ? uAxN : Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = (!evapMicroRegion[cp.jCellOut]) ? uBxN : Flux(uBxN_fict, uBxN);  // flux on B-side


            return FlxNeg * vA - FlxPos * vB;

        }

        private double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in + UxN_out);
        }


        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        BitArray evapMicroRegion;

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
        public double InnerEdgeForm(ref CommonParams cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);
            double uBxN = GenericBlas.InnerProd(U_Pos, cp.Normal);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict = uBxN;
            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict = uAxN;

            double FlxNeg = (DirichletCond) ? 0.0 : -Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = (DirichletCond) ? 0.0 : Flux(uBxN_fict, uBxN);  // flux on B-side


            return (FlxNeg * vA - FlxPos * vB);

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
                if (DirichletCond)
                    return TermActivationFlags.V;
                else
                    return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }


    public class TemperatureGradientAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        int m_d;

        LevelSetTracker m_LsTrk;

        public TemperatureGradientAtLevelSet(int _d, LevelSetTracker lstrk, double _kA, double _kB, double _Tsat) {

            this.m_d = _d;
            m_LsTrk = lstrk;

            this.kA = _kA;
            this.kB = _kB;

            //this.DirichletCond = _DiriCond;
            this.Tsat = _Tsat;
        }

        double kA;
        double kB;

        //bool DirichletCond;
        double Tsat;

        /// <summary>
        /// 
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(inp.jCellIn == inp.jCellOut);

            double FlxNeg = 0.0;
            double FlxPos = 0.0;
            if (!evapMicroRegion[inp.jCellIn]) {
                double Avg = Tsat;
                FlxNeg += kA * Avg; // + 0.5 * uB[0] * (kA - kB);
                FlxPos += kB * Avg; // + 0.5 * uA[0] * (kA - kB);
            } else {
                double Avg = 0.5 * (uB[0] + uA[0]);
                FlxNeg += kA * Avg; // + 0.5 * uB[0] * (kA - kB);
                FlxPos += kB * Avg; // + 0.5 * uA[0] * (kA - kB);
            }

            return (FlxNeg * vA - FlxPos * vB) * inp.Normal[m_d];

            //double Acc = (DirichletCond) ? Tsat : 0.5 * (uB[0] + uA[0]);
            //return Acc * (kA * vA - kB * vB) * inp.n[m_d];
        }


        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        BitArray evapMicroRegion;


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

        public TemperatureStabilizationFormAtLevelSet(int _d, LevelSetTracker lstrk, double _kA, double _kB, bool _DiriCond, double _Tsat) {

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
        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            //return (uA[0] - uB[0]) * inp.n[m_d] * (vA - vB);

            double Acc = 0.0;

            if (DirichletCond) {
                Acc += 2.0 * (uA[0] - Tsat) * inp.Normal[m_d] * (vA - 0.0);
                Acc += 2.0 * (Tsat - uB[0]) * inp.Normal[m_d] * (0.0 - vB);
            } else {
                Acc += (kA * uA[0] - kB * uB[0]) * inp.Normal[m_d] * (vA - vB);
                //Acc += (uA[0] - uB[0]) * inp.Normal[m_d] * (vA - vB);
            }

            return -Acc;
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
                if(DirichletCond)
                    return TermActivationFlags.UxV | TermActivationFlags.V;
                else
                    return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

    }

}
