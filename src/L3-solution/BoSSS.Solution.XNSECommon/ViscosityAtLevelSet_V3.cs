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
using System.Linq;
using System.Text;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;
using System.Collections;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.XNSECommon.Operator.Viscosity {

    /// <summary>
    /// Coupling of the fully symmetric viscous stress, i.e.
    /// ```math
    ///    - \mu \left( \nabla \vec{u} + \nabla \vec{u}^T \right) .
    /// ```
    /// Must be used together with <see cref="ViscosityInSpeciesBulk_GradUTerm"/> **and** <see cref="ViscosityInSpeciesBulk_GradUtranspTerm"/> in the bulk.
    /// </summary>
    public class ViscosityAtLevelSet_FullySymmetric : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        //LevelSetTracker m_LsTrk;

        public ViscosityAtLevelSet_FullySymmetric(int SpaceDim, double _muA, double _muB, double _penalty, int _component, 
            bool _freeSurface = false, double _sI = 0.0) {
            //this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.m_penalty_base = _penalty;
            this.component = _component;
            this.m_D = SpaceDim;
            this.freeSurface = _freeSurface;
            this.sI = _sI;
            this.beta = sI.IsInfinity() ? 0.0 : (muA + muB)/(2 * sI);
        }

        double muA;
        double muB;
        double sI;
        double beta;
        int component;
        int m_D;

        bool freeSurface;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;
            //double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];

            int D = N.Length;
            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double[] Grad_uA_xN = new double[D], Grad_uB_xN = new double[D];
            double Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    Grad_uA_xN[dd] += Grad_uA[dd, d] * N[d];
                    Grad_uB_xN[dd] += Grad_uB[dd, d] * N[d];
                }
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }
            double Ret = 0.0;

            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);


            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);

            double wA;
            double wB;
            double wPenalty;
            wA = 0.5;
            wB = 0.5;
            wPenalty = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;


            if (!freeSurface && sI == 0.0) {
                Ret -= (wA * muA * Grad_uA_xN[component] + wB * muB * Grad_uB_xN[component]) * (vA - vB);                           // consistency term
                Ret -= (wA * muA * Grad_vA_xN + wB * muB * Grad_vB_xN) * (uA[component] - uB[component]);     // symmetry term
                //Ret += (penalty / hCutCellMin) * (uA[component] - uB[component]) * (vA - vB) * wPenalty; // penalty term
                Ret += (uA[component] - uB[component]) * (vA - vB) * pnlty * wPenalty; // penalty term
                // Transpose Term
                for (int i = 0; i < D; i++) {
                    Ret -= (wA * muA * Grad_uA[i, component] + wB * muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                    Ret -= (wA * muA * Grad_vA[i] + wB * muB * Grad_vB[i]) * (uA[i] - uB[i]) * N[component];  // symmetry term
                }
            } else if(!freeSurface && sI > 0.0){
                // handle normal and tangential components seperately

                // GradU in normal direction
                for(int dN = 0; dN < D; dN++) {
                    for(int dD = 0; dD < D; dD++) {
                        // consistency
                        Ret -=  (wA * muA * inp.Normal[dN] * Grad_uA[dN, dD] * inp.Normal[dD] + wB * muB * inp.Normal[dN] * Grad_uB[dN, dD] * inp.Normal[dD]) * ((vA - vB) * inp.Normal[component]);
                        // symmetry
                        Ret -= (wA * muA * inp.Normal[component] * Grad_vA[dD] * inp.Normal[dD] + wB * muB * inp.Normal[component] * Grad_vB[dD] * inp.Normal[dD]) * (uA[dN] - uB[dN]) * inp.Normal[dN];
                    }
                    // penalty
                    Ret += ((uA[dN] - uB[dN]) * inp.Normal[dN]) * ((vA - vB) * inp.Normal[component]) * pnlty * wPenalty;
                }

                // GradU Transpose
                for(int dN = 0; dN < D; dN++) {
                    for(int dD = 0; dD < D; dD++) {
                        // consistency
                        Ret -= (wA * muA * inp.Normal[dN] * Grad_uA[dD, dN] * inp.Normal[dD] + wB * muB * inp.Normal[dN] * Grad_uB[dD, dN] * inp.Normal[dD]) * ((vA - vB) * inp.Normal[component]);
                        Ret -= (wA * muA * inp.Normal[dN] * Grad_vA[dN] * inp.Normal[component] + wB * muB * inp.Normal[dN] * Grad_vB[dN] * inp.Normal[component]) * (uA[dD] - uB[dD]) * inp.Normal[dD] ;
                    }
                    // penalty
                    // Ret += ((uA[dN] - uB[dN]) * inp.Normal[dN]) * ((vA - vB) * inp.Normal[component]) * pnlty * wPenalty;
                }

                double[,] P = new double[D, D];
                for(int d1 = 0; d1 < D; d1++) {
                    for(int d2 = 0; d2 < D; d2++) {
                        double nn = inp.Normal[d1] * inp.Normal[d2];
                        if(d1 == d2) {
                            P[d1, d2] = 1 - nn;
                        } else {
                            P[d1, d2] = -nn;
                        }
                    }
                }

                // tangential dissipation force term aka slip
                for(int d1 = 0; d1 < D; d1++) {
                    for(int d2 = 0; d2 < D; d2++) {
                        Ret += (beta * P[d1, d2] * (uA[d2] - uB[d2])) * (P[d1, component] * (vA - vB));
                    }
                }
            } else {
                //free slip
                for (int d = 0; d < D; d++) {
                    Ret -= N[d] * (wA * muA * Grad_uA_xN[d] + wB * muB * Grad_uB_xN[d]) * (vA - vB) * N[component];                           // consistency term
                    Ret -= N[component] * (wA * muA * Grad_vA_xN + wB * muB * Grad_vB_xN) * (uA[d] - uB[d]) * N[d];     // symmetry term
                    //Ret += (penalty / hCutCellMin) * (uA[d] - uB[d]) * N[d] * (vA - vB) * N[component] * wPenalty; // penalty term
                    Ret += (uA[d] - uB[d]) * N[d] * (vA - vB) * N[component] * pnlty * wPenalty; // penalty term
                }
                // Transpose Term
                for (int dN = 0; dN < D; dN++) {
                    for (int dD = 0; dD < D; dD++) {
                        Ret -= N[dN] * (wA * muA * Grad_uA[dD, dN] + wB * muB * Grad_uB[dD, dN]) * N[dD] * (vA - vB) * N[component];  // consistency term
                        Ret -= N[dN] * (wA * muA * Grad_vA[dN] + wB * muB * Grad_vB[dN]) * N[component] * (uA[dD] - uB[dD]) * N[dD];  // symmetry term
                    }
                }
            }
            return Ret;
        }



        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
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

            double scaledPenalty = penaltySizeFactor * m_penalty * m_penalty_base;
            if(scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

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
        }

        //private static bool rem = true;

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            //var DerivEdg = new LevelSetFormDifferentiator(this, SpatialDimension);
            //return new IEquationComponent[] { DerivEdg };
            return new IEquationComponent[] { this };
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }


        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.GradUxGradV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null;  }
        }

     
    }    

    //public class ViscosityAtLevelSet_Slip : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

    //    LevelSetTracker m_LsTrk;

    //    public ViscosityAtLevelSet_Slip(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component) {
    //        this.m_LsTrk = lstrk;
    //        this.muA = _muA;
    //        this.muB = _muB;
    //        this.penalty = _penalty;
    //        this.component = _component;
    //        this.m_D = lstrk.GridDat.SpatialDimension;
    //    }

    //    double muA;
    //    double muB;
    //    double penalty;
    //    int component;
    //    int m_D;


    //    /// <summary>
    //    /// 
    //    /// </summary>
    //    public double InnerEdgeForm(ref CommonParams inp,
    //        double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
    //        double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
    //        double[] N = inp.n;
    //        double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

    //        int D = N.Length;
    //        Debug.Assert(this.ArgumentOrdering.Count == D);
    //        Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
    //        Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
    //        Debug.Assert(Grad_uA.GetLength(1) == D);
    //        Debug.Assert(Grad_uB.GetLength(1) == D);


    //        double Ret = 0.0;

    //        double PosCellLengthScale = PosLengthScaleS[inp.jCell];
    //        double NegCellLengthScale = NegLengthScaleS[inp.jCell];

    //        double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
    //        if(hCutCellMin <= 1.0e-10 * hCellMin)
    //            // very small cell -- clippling
    //            hCutCellMin = hCellMin;

    //        Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
    //        Debug.Assert(uB.Length == this.ArgumentOrdering.Count);

    //        double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;

    //        if(NegSlipLengths[inp.jCell] == 0.0 && NegSlipLengths[inp.jCell] == 0.0) {

    //            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
    //            for(int d = 0; d < D; d++) {
    //                Grad_uA_xN += Grad_uA[component, d] * N[d];
    //                Grad_uB_xN += Grad_uB[component, d] * N[d];
    //                Grad_vA_xN += Grad_vA[d] * N[d];
    //                Grad_vB_xN += Grad_vB[d] * N[d];
    //            }

    //            Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
    //            Ret -= 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * (uA[component] - uB[component]);     // symmetry term
    //            Ret += (penalty / hCutCellMin) * (uA[component] - uB[component]) * (vA - vB) * muMax; // penalty term
    //                                                                                                  // Transpose Term
    //            for(int i = 0; i < D; i++) {
    //                Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
    //                Ret -= 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * (uA[i] - uB[i]) * N[component];  // symmetry term
    //            }

    //        } else {

    //            for(int dN = 0; dN < m_D; dN++) {
    //                for(int dD = 0; dD < m_D; dD++) {
    //                    Ret -= 0.5 * N[dN] * (muA * Grad_uA[dD, dN] + muB * Grad_uB[dN, dD]) * N[dD] * (vA - vB) * N[component];     // consistency term
    //                    Ret -= 0.5 * N[dN] * (muA * Grad_vA[dN] + muB * Grad_vB[dD]) * N[component] * (uA[dD] - uB[dN]) * N[dD];    // symmetry term
    //                                                                                                                                // transposed term
    //                    Ret -= 0.5 * N[dN] * (muA * Grad_uA[dN, dD] + muB * Grad_uB[dN, dD]) * N[dD] * (vA - vB) * N[component];    // consistency term
    //                    Ret -= 0.5 * N[component] * (muA * Grad_vA[dD] + muB * Grad_vB[dD]) * N[dD] * (uA[dN] - uB[dN]) * N[dN];    // symmetry term
    //                }
    //                Ret += (penalty / hCutCellMin) * (uA[dN] - uB[dN]) * (vA - vB) * N[dN] * muMax;     // penalty term
    //            }
    //        }

    //        return Ret;
    //    }


    //    MultidimensionalArray PosLengthScaleS;
    //    MultidimensionalArray NegLengthScaleS;

    //    MultidimensionalArray PosSlipLengths;
    //    MultidimensionalArray NegSlipLengths;

    //    public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
    //        NegLengthScaleS = csA.CellLengthScales;
    //        PosLengthScaleS = csB.CellLengthScales;

    //        NegSlipLengths = (MultidimensionalArray)csA.UserDefinedValues["SlipLengths"];
    //        PosSlipLengths = (MultidimensionalArray)csB.UserDefinedValues["SlipLengths"];
    //    }


    //    public int LevelSetIndex {
    //        get { return 0; }
    //    }

    //    public IList<string> ArgumentOrdering {
    //        get { return VariableNames.VelocityVector(this.m_D); }
    //    }

    //    public SpeciesId PositiveSpecies {
    //        get { return m_LsTrk.GetSpeciesId("B"); }
    //    }

    //    public SpeciesId NegativeSpecies {
    //        get { return m_LsTrk.GetSpeciesId("A"); }
    //    }

    //    public TermActivationFlags LevelSetTerms {
    //        get {
    //            return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
    //        }
    //    }

    //    public IList<string> ParameterOrdering {
    //        get { return null; }
    //    }


    //}

}
