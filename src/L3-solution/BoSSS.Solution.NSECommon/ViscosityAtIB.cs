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

namespace BoSSS.Solution.NSECommon.Operator.Viscosity {

    public class ViscosityAtIB : BoSSS.Foundation.XDG.ILevelSetComponent {

        LevelSetTracker m_LsTrk;

        public ViscosityAtIB(int _d, int _D, LevelSetTracker t, double penalty, double _muA, Func<double, double>[] _uLevSet, Func<double, double>[] _wLevSet, double particleRadius) {

            m_penalty = penalty;
            this.m_LsTrk = t;
            this.muA = _muA;
            this.component = _d;
            this.uLevSet = _uLevSet;
            this.wLevSet = _wLevSet;
            this.m_D = _D;
            this.pRadius = particleRadius;
        }

        int component;
        int m_D;
        double pRadius;

        Func<double, double>[] uLevSet, wLevSet;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        double muA;

        double m_penalty;

        /// <summary>
        /// penalty param
        /// </summary>
        double penalty(double hmin) {
            double µ = m_penalty / hmin;
            Debug.Assert(!(double.IsNaN(µ) || double.IsInfinity(µ)));
            return µ;
        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            //Debug.Assert(!double.IsNaN(inp.PosCellLengthScale));
            //double hCutCellMin = Math.Min(inp.PosCellLengthScale, inp.NegCellLengthScale);
            double hCutCellMin = inp.NegCellLengthScale; // for IBM, there is no positve species!
            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;
            double _penalty = penalty(hCutCellMin);


            int D = N.Length;
            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[component, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
            }


            // Evaluate the complete velocity as a sum of translation and angular velocity
            double uAFict;
            double Ret = 0.0;


            // 3D for IBM_Solver
            if (inp.x.Length == 3) {
                
                Ret -= Grad_uA_xN * (vA);                           // consistency term
                Ret -= Grad_vA_xN * (uA[component] - 0);     // symmetry term
                Ret += _penalty * (uA[component] - 0) * (vA); // penalty term

                Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
                return Ret * muA;
            }


            if (component == 0) {
                uAFict = (uLevSet[component])(inp.time) + pRadius * wLevSet[0](inp.time) * -inp.n[1];
            } else {
                uAFict = (uLevSet[component])(inp.time) + pRadius * wLevSet[0](inp.time) * inp.n[0];
            }


            Ret -= Grad_uA_xN * (vA);                           // consistency term
            Ret -= Grad_vA_xN * (uA[component] - uAFict);     // symmetry term
            Ret += _penalty * (uA[component] - uAFict) * (vA); // penalty term

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret * muA;
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
    }

    /*
    
    public class ViscosityAtIB2 
    {

        public ViscosityAtIB2(int _d, int _D, LevelSetTracker t, double penalty,
            double _Re, Func<double, double>[] _uLevSet)
            : base(t)
        {
            m_penalty = penalty;
            this.Re = _Re;
            this.m_d = _d;
            this.uLevSet = _uLevSet;
            this.m_D = _D;
        }

        int m_d;
        int m_D;

        Func<double, double>[] uLevSet;

        /// <summary>
        /// Reynolds in species A
        /// </summary>
        double Re;

        double m_penalty;

        /// <summary>
        /// penalty param
        /// </summary>
        double penalty(double hmin)
        {
            double µ = m_penalty / hmin;
            Debug.Assert(!(double.IsNaN(µ) || double.IsInfinity(µ)));
            return µ;
        }

        public override void DerivativVar_NodeFunc(ref Linear2ndDerivativeCouplingFlux.CommonParams inp)
        {
        }


        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp,
            double[] U_A, double[] U_B, double[,] GradU_A, double[,] GradU_B)
        {

            int D = GradU_A.GetLength(1);
            Debug.Assert(D == this.m_D);
            int _d = m_d;//T.d;
                         //T.PassParams(cp, U_A, U_B, GradU_A, GradU_B);

            // Extract values from parameters
            // ==============================

            unsafe
            {
                fixed (double* pGradU_A = GradU_A, pGradU_B = GradU_B, pN = cp.n)
                {

                    // u
                    Debug.Assert(GradU_A.GetLength(0) == this.ArgumentOrdering.Count);
                    Debug.Assert(GradU_B.GetLength(0) == this.ArgumentOrdering.Count);

                    double uA = U_A[_d];
                    double uB = U_B[_d];


                    //// interface velocity
                    //double s = cp.ParamsPos[D];
                    //Debug.Assert(s == cp.ParamsNeg[D], "Interface velocity must be continuous across the level-set.");
                    double s = double.NaN;

                    // Computation of fictitious values: transform values from A to B and vice-versa
                    // =============================================================================

                    // fictitious u
                    double uAFict;
                    uAFict = (uLevSet[_d])(cp.time);
                    //TransformU(_d, D, s, cp.n, cp.x, U_A, U_B, out uAFict, out uBFict);

                    // ficticious Gradient times Normal
                    double grUA;
                    TransformGradU(_d, D, pN, pGradU_A, pGradU_B, out grUA);

                    // compute the fluxes
                    // ==================

                    // penalty
                    double hCutCellMin = Math.Min(cp.PosCellLengthScale, cp.NegCellLengthScale);
                    if (hCutCellMin <= 1.0e-10 * cp.hCellMin)
                        // very small cell -- clippling
                        hCutCellMin = cp.hCellMin;
                    double _penalty = penalty(hCutCellMin);

                    // kinematic viscosity for both phases, resp. dynamic viscosity for the u_p_2 -- case

                    // Finally, compute the fluxes
                    FlxNeg = 0;
                    FlxNeg += _penalty * (uA - uAFict) * 1/Re; //VORZEICHEN UMGEDREHT
                    FlxNeg -= grUA * 1/Re; //VORZEICHEN UMGEDREHT
                    FlxPos = 0;
                }
            }
        }

        static unsafe void TransformGradU(
            int _d, int D, double* n, double* GradU_A, double* GradU_B,
            out double grUA)
        {
            Debug.Assert(_d < D);

            // Gradient U * normal
            grUA = 0;
            double* pGradU_A_d = GradU_A + _d * D;
            double* pGradU_B_d = GradU_B + _d * D;
            for (int d1 = 0; d1 < D; d1++)
            {
               grUA += pGradU_A_d[d1] * n[d1];
            }
        }


        unsafe static double TransformGradU_2D(int d, double muA, double muB, double* N, double NormalSign, double* GradUA)
        {
            //double uxA = GradUA[0, 0];
            //double uyA = GradUA[0, 1];
            //double vxA = GradUA[1, 0];
            //double vyA = GradUA[1, 1];
            //Debug.Assert(GradUA.GetLength(1) == 2);
            double uxA = *GradUA;
            double uyA = *(GradUA + 1);
            double vxA = *(GradUA + 2);
            double vyA = *(GradUA + 3);
            Debug.Assert(Math.Abs(NormalSign) == 1.0);


            double NX = N[0] * NormalSign;
            double NY = N[1] * NormalSign;
            double NXp2 = NX * NX;
            double NYp2 = NY * NY;

            switch (d)
            {
                case 0:
                    {
                        return (NormalSign / muB) * (
                              uxA * (2 * NYp2 * NX * (-muB + muA))
                            + uyA * (-NY * (-2 * muB * NXp2 + 2 * muA * NXp2 - muA))
                            + vxA * (-NY * (muB + 2 * muA * NXp2 - muA - 2 * muB * NXp2))
                            + vyA * (NX * (muB - 2 * muB * NXp2 + 2 * muA * NXp2 - 2 * muA))
                            );
                    }

                case 1:
                    {
                        return (NormalSign / muB) * (
                               uxA * (-NY * (2 * muA * NXp2 - 2 * muB * NXp2 + muB))
                             + uyA * (NX * (muB + 2 * muA * NXp2 - muA - 2 * muB * NXp2))
                             + vxA * (NX * (-2 * muB * NXp2 - muA + 2 * muB + 2 * muA * NXp2))
                             + vyA * (2 * NY * NXp2 * (-muB + muA))
                            );

                    }

                default: throw new ArgumentException();
            }

        }



        /// <summary>
        /// velocity transformation according to tangential no-slip.
        /// </summary>
        double TransformU_Noslip(int d, int D, double[] n, double[] Ufremd)
        {
            Debug.Assert(n.Length == D);
            Debug.Assert(d >= 0 && d < D);

            double acc = 0;
            for (int d1 = 0; d1 < D; d1++)
                acc += n[d1] * (Ufremd[d1]);


            return Ufremd[d] - acc * n[d];
        }


        void TransformU(int d, int D, double s, double[] n, double[] X, double[] uA, double[] uB, out double uAfict_d, out double uBfict_d)
        {
            Debug.Assert(d < D);

            //if (!T.MaterialInterface) {
            //     double uAxN = BLAS.ddot(D, n, 1, uA, 1);
            //     double uBxN = BLAS.ddot(D, n, 1, uB, 1);

            //     double rhoJmp = m_rhoB - m_rhoA;


            //    // normal component: using the continuity in normal direction
            //    double uAxN_fict, uBxN_fict;

            //    uAxN_fict = (1 / m_rhoA) * (m_rhoB * uBxN + (-s) * rhoJmp); // transform from species B to A: we call this the "A-fictitious" value
            //    uBxN_fict = (1 / m_rhoB) * (m_rhoA * uAxN - (-s) * rhoJmp); // transform from species A to B: we call this the "B-fictitious" value


            //    // Variant: no constraint on normal component
            //    //double uAxN_fict = uAxN;
            //    //double uBxN_fict = uBxN;

            //    // return
            //    uBfict_d = uBxN_fict*n[d] + TransformU_Noslip(d, D, n, uA);
            //    uAfict_d = uAxN_fict*n[d] + TransformU_Noslip(d, D, n, uB);
            //} else {
            uBfict_d = uA[d];
            uAfict_d = uB[d];
            //}
        }



        public override void PrimalVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos)
        {
            //EquationAndVarMode varMode = T.varmode;
            int m_D = this.m_D;
            int m_d = this.m_d;       

            // extract velocity component
            double uNeg = U_Neg[m_d];
            double uPos = U_Pos[m_d];


            // interface velocity
            //double s = cp.ParamsPos[m_D];
            //Debug.Assert(s == cp.ParamsNeg[m_D], "Interface velocity must be continuous across the level-set.");
            double s = double.NaN;

            // Fluxes
            FlxNeg = (uLevSet[m_d])(cp.time);     //VORZEICHEM UMGEDREHT
            FlxPos = 0;    
        }

        public override void FluxPotential(out double G, double[] U)
        {
            //G = 0;
            G = U[this.m_d]; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        }

        public override void Nu(out double NuNeg, out double NuPos, ref CommonParamsLs cp)
        {
            NuPos = 1;
            NuNeg = 1;
        }

        public override int LevelSetIndex
        {
            get { return 0; }
        }

        public override IList<string> ArgumentOrdering
        {
            get
            {
                return VariableNames.VelocityVector(m_D);
            }
        }

        public override IList<string> ParameterOrdering
        {
            get
            {
                return null;
            }
        }

        public override SpeciesId PositiveSpecies
        {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public override SpeciesId NegativeSpecies
        {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }
    }
    */
}