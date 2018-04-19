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
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Platform;

namespace BoSSS.Solution.XNSECommon.Operator.Viscosity {
    
    /*
    /// <summary>
    /// Explicit transformation; works only in 2D;
    /// gives at least a symmetric form for equal viscosities.
    /// </summary>
    public class ViscosityAtLevelSet_Explicit : Linear2ndDerivativeCouplingFlux {

        public ViscosityAtLevelSet_Explicit(int _d, int _D, LevelSetTracker t, double penalty,
            double _muA, double _muB)
            : base(t) {
            m_penalty = penalty;
            this.m_muA = _muA;
            this.m_muB = _muB;
            this.m_D = _D;
            this.m_d = _d;
        }

        int m_d;
        int m_D;

        /// <summary>
        /// Viscosity in species A
        /// </summary>
        double m_muA;

        /// <summary>
        /// Viscosity in species B
        /// </summary>
        double m_muB;

        double m_penalty;

        /// <summary>
        /// penalty param
        /// </summary>
        double penalty(double hmin) {
            double µ = m_penalty/hmin;
            Debug.Assert(!(double.IsNaN(µ) || double.IsInfinity(µ)));
            return µ;
        }

        public override void DerivativVar_NodeFunc(ref Linear2ndDerivativeCouplingFlux.CommonParams inp) {
        }


        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp,
            double[] U_A, double[] U_B, double[,] GradU_A, double[,] GradU_B) {
       
            int D = GradU_A.GetLength(1);
            Debug.Assert(D == this.m_D);
            int _d = m_d;//T.d;
            //T.PassParams(cp, U_A, U_B, GradU_A, GradU_B);
            
            // Extract values from parameters
            // ==============================

            unsafe {
                fixed (double* pGradU_A = GradU_A, pGradU_B = GradU_B, pN = cp.n) {

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
                    double uBFict;
                    TransformU(_d, D, s, cp.n, cp.x, U_A, U_B, out uAFict, out uBFict);

                    // ficticious Gradient times Normal
                    double grUA, grUB;
                    double grUBFict, grUAFict;
                    TransformGradU(_d, D, m_muA, m_muB, pN, pGradU_A, pGradU_B, out grUA, out grUB, out grUBFict, out grUAFict);

#if DEBUG
                    {
                        double grUA_check = 0.0;
                        double grUB_check = 0.0;

                        Debug.Assert(GradU_A.GetLength(0) == D);
                        Debug.Assert(GradU_A.GetLength(1) == D);
                        Debug.Assert(GradU_B.GetLength(0) == D);
                        Debug.Assert(GradU_B.GetLength(1) == D);

                        for (int d1 = 0; d1 < D; d1++) {
                            grUA_check += GradU_A[_d, d1]*cp.n[d1];
                            grUB_check += GradU_B[_d, d1]*cp.n[d1];
                        }

                        Debug.Assert(Math.Abs(grUA_check - grUA) < 1.0e-10);
                        Debug.Assert(Math.Abs(grUB_check - grUB) < 1.0e-10);
                    }
#endif

                    // compute the fluxes
                    // ==================

                    // penalty
                    double hCutCellMin = Math.Min(cp.PosCellLengthScale, cp.NegCellLengthScale);
                    if(hCutCellMin <= 1.0e-10 * cp.hCellMin)
                        // very small cell -- clippling
                        hCutCellMin = cp.hCellMin;
                    double _penalty = penalty(hCutCellMin);

                    // kinematic viscosity for both phases, resp. dynamic viscosity for the u_p_2 -- case
                    double nuA, nuB;
                    Nu(out nuA, out nuB, ref cp);

                    // Finally, compute the fluxes
                    FlxNeg = 0;
                    FlxNeg -= _penalty*(uA - uAFict)*nuA; //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    FlxPos = 0;
                    FlxPos -= _penalty*(uBFict - uB)*nuB; //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    FlxNeg += 0.5 * (grUA * nuA + grUAFict * nuA);
                    FlxPos += 0.5 * (grUBFict * nuB + grUB * nuB);
                }
            }
        }

        static unsafe void TransformGradU(
            int _d, int D, double muA, double muB, double* n, double* GradU_A, double* GradU_B,
            out double grUA, out double grUB,
            out double grUBFict, out double grUAFict) {
            Debug.Assert(_d < D); 

            // Gradient U * normal
            grUA = 0;
            grUB = 0;
            double* pGradU_A_d = GradU_A + _d*D;
            double* pGradU_B_d = GradU_B + _d*D;
            for (int d1 = 0; d1 < D; d1++) {
                grUA += pGradU_A_d[d1]*n[d1];
                grUB += pGradU_B_d[d1]*n[d1];
            }

            {
                // new style (august 13) - fully explicit transformation
                // NO dependence on pressure !
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                
                switch (D) {

                    case 2:
                    grUBFict = TransformGradU_2D(_d, muA, muB, n, +1.0, GradU_A);
                    grUAFict = TransformGradU_2D(_d, muB, muA, n, -1.0, GradU_B);
                    break;
                    default:
                    throw new NotImplementedException();
                }
                
            }
        }


        unsafe static double TransformGradU_2D(int d, double muA, double muB, double* N, double NormalSign, double* GradUA) {
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


            double NX = N[0]*NormalSign;
            double NY = N[1]*NormalSign;
            double NXp2 = NX*NX;
            double NYp2 = NY*NY;

            switch (d) {
                case 0: {
                    return (NormalSign/muB)*(
                          uxA*(2*NYp2*NX*(-muB+muA))
                        + uyA*(-NY*(-2*muB*NXp2+2*muA*NXp2-muA))
                        + vxA*(-NY*(muB+2*muA*NXp2-muA-2*muB*NXp2))
                        + vyA*(NX*(muB-2*muB*NXp2+2*muA*NXp2-2*muA))
                        );
                }

                case 1: {
                    return (NormalSign/muB)*(
                           uxA*(-NY*(2*muA*NXp2-2*muB*NXp2+muB))
	                     + uyA*(NX*(muB+2*muA*NXp2-muA-2*muB*NXp2))
	                     + vxA*(NX*(-2*muB*NXp2-muA+2*muB+2*muA*NXp2))
	                     + vyA*(2*NY*NXp2*(-muB+muA))
                        );

                }
                
                default: throw new ArgumentException();
            }

        }



        /// <summary>
        /// velocity transformation according to tangential no-slip.
        /// </summary>
        double TransformU_Noslip(int d, int D, double[] n, double[] Ufremd) {
            Debug.Assert(n.Length == D);
            Debug.Assert(d >= 0 && d < D);

            double acc = 0;
            for( int d1 = 0; d1 < D; d1++)
                acc += n[d1]*(Ufremd[d1]);


            return Ufremd[d] - acc*n[d];
        }


        void TransformU(int d, int D, double s, double[] n, double[] X, double[] uA, double[] uB, out double uAfict_d, out double uBfict_d) {
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
            double[] U_Neg, double[] U_Pos) {

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

            // fictitious u
            double uNegFict;
            double uPosFict;
            TransformU(m_d, m_D, s, cp.n, cp.x, U_Neg, U_Pos, out uNegFict, out uPosFict);

            // Fluxes
            FlxNeg = 0.5*(uNeg + uNegFict);    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            FlxPos = 0.5*(uPosFict + uPos);    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        }

        public override void FluxPotential(out double G, double[] U) {
            G = U[this.m_d]; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        }

        public override void Nu(out double NuNeg, out double NuPos, ref CommonParamsLs cp) {
            NuPos = -this.m_muB;
            NuNeg = -this.m_muA;
        }

        public override int LevelSetIndex {
            get { return 0; }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D); 
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public override SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public override SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }
    }
    */
}
