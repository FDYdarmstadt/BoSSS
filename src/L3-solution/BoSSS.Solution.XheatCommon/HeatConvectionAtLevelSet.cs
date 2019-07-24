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
using ilPSP.Utils;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using System.Collections;

namespace BoSSS.Solution.XheatCommon {


    public class HeatConvectionAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public HeatConvectionAtLevelSet(int _D, LevelSetTracker LsTrk, double _capA, double _capB, double _LFFA, double _LFFB, 
            ThermalMultiphaseBoundaryCondMap _bcmap, bool _movingmesh, bool _DiriCond, double _Tsat) {

            m_D = _D;

            m_LsTrk = LsTrk;

            //MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, _LFFA, double.NaN, LsTrk);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"));
            PosFlux = new HeatConvectionInBulk(_D, _bcmap, _capA, _capB, double.NaN, _LFFB, LsTrk);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"));


            DirichletCond = _DiriCond;
            Tsat = _Tsat;

            capA = _capA;
            capB = _capB;
            LFFA = _LFFA;
            LFFB = _LFFB;

        }

        //bool MaterialInterface;
        int m_D;

        double capA;
        double capB;
        double LFFA;
        double LFFB;

        bool DirichletCond;
        double Tsat;

        // Use Fluxes as in Bulk Convection
        HeatConvectionInBulk NegFlux;
        HeatConvectionInBulk PosFlux;

        

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            //if(this.MaterialInterface) {

                U_NegFict = U_Pos;
                U_PosFict = U_Neg;

            //} else {
            //    throw new NotImplementedException();
            //}
        }


        public double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;


            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);

            //Flux for negativ side
            double FlxNeg;
            if (DirichletCond) {

                double r = 0.0;

                // Calculate central part
                // ======================
                r += Tsat * (ParamsNeg[0] * cp.n[0] + ParamsNeg[1] * cp.n[1]);
                if (m_D == 3) {
                    r += Tsat * ParamsNeg[2] * cp.n[2];
                }

                // Calculate dissipative part
                // ==========================

                double[] VelocityMeanIn = new double[m_D];
                for (int d = 0; d < m_D; d++) {
                    VelocityMeanIn[d] = ParamsNeg[m_D + d];
                }

                double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.n, true);

                double uJump = U_Neg[0] - Tsat;

                r += LambdaIn * uJump * LFFA;

                FlxNeg = capA * r;


            } else { 

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;

                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
                //Console.WriteLine("FlxNeg = {0}", FlxNeg);
            }

            // Flux for positive side
            double FlxPos;
            if (DirichletCond) {

                double r = 0.0;

                // Calculate central part
                // ======================
                r += Tsat * (ParamsPos[0] * cp.n[0] + ParamsPos[1] * cp.n[1]);
                if (m_D == 3) {
                    r += Tsat * ParamsPos[2] * cp.n[2];
                }

                // Calculate dissipative part
                // ==========================

                double[] VelocityMeanOut = new double[m_D];
                for (int d = 0; d < m_D; d++) {
                    VelocityMeanOut[d] = ParamsPos[m_D + d];
                }

                double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.n, true);

                double uJump = Tsat - U_Pos[0];

                r += LambdaOut * uJump * LFFB;

                FlxPos = capB * r;

            } else {

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
                //Console.WriteLine("FlxPos = {0}", FlxPos);
            }

            if(movingmesh)
                return 0.0;
            else
                return FlxNeg * v_Neg - FlxPos * v_Pos;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
    }


    public class HeatConvectionAtLevelSet_Divergence : EvaporationAtLevelSet {


        public HeatConvectionAtLevelSet_Divergence(int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            ThermalParameters thermParams, double _Rint, double _sigma, bool _DiriCond = true) {
            //double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.capA = thermParams.c_A * _rhoA;
            this.capB = thermParams.c_B * _rhoB;

            this.m_LsTrk = lsTrk;


            this.kA = thermParams.k_A;
            this.kB = thermParams.k_B;
            this.hVapA = thermParams.hVap_A;
            this.Rint = _Rint;

            this.Tsat = thermParams.T_sat;
            this.sigma = _sigma;
            this.pc = thermParams.pc;

            this.DirichletCond = _DiriCond;
        }
 
        double rhoA;
        double rhoB;
        double capA;
        double capB;

        bool DirichletCond;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }


        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - GeneralizedDivergenceAtLevelSet: {0}", M);

            return M;

        }

        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.ParamsNeg, cp.ParamsPos, cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;


            double Tavg = (DirichletCond) ? Tsat : 0.5 * (U_Neg[0] + U_Pos[0]);
            //double T_sat = 0.0;

            double uAxN = -M * (1 / rhoA) * (Tavg);
            double uBxN = -M * (1 / rhoB) * (Tavg);

            // transform from species B to A: we call this the "A-fictitious" value
            double uAxN_fict;
            //uAxN_fict = (1 / rhoA) * (rhoB * uBxN);
            uAxN_fict = uBxN;

            // transform from species A to B: we call this the "B-fictitious" value
            double uBxN_fict;
            //uBxN_fict = (1 / rhoB) * (rhoA * uAxN);
            uBxN_fict = uAxN;


            // compute the fluxes: note that for the continuity equation, we use not a real flux,
            // but some kind of penalization, therefore the fluxes have opposite signs!
            double FlxNeg = -Flux(uAxN, uAxN_fict); // flux on A-side
            double FlxPos = +Flux(uBxN_fict, uBxN);  // flux on B-side

            FlxNeg *= capA;
            FlxPos *= capB;

            double Ret = FlxNeg * vA - FlxPos * vB;

            return -Ret;
        }


        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out) {
            return 0.5 * (UxN_in - UxN_out);
        }


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }


    }


}
