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
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using System.Collections;


namespace BoSSS.Solution.XheatCommon {


    public class ConvectionAtLevelSet_nonMaterialLLF : EvaporationAtLevelSet {


        public ConvectionAtLevelSet_nonMaterialLLF(int _d, int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            ThermalParameters thermParams, double _Rint, double _sigma) {
            //double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.m_d = _d;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_LsTrk = lsTrk;

            this.kA = thermParams.k_A;
            this.kB = thermParams.k_B;
            this.hVapA = thermParams.hVap_A;
            this.Rint = _Rint;

            this.Tsat = thermParams.T_sat;
            this.sigma = _sigma;
            this.pc = thermParams.pc;
        }

        int m_d;
        double rhoA;
        double rhoB;



        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - ConvectionAtLevelSet_nonMaterialLLF: {0}", M);

            return M;

        }


        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.ParamsNeg.GetSubVector(2 * D, D + 3), cp.ParamsPos.GetSubVector(2 * D, D + 3), cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;


            double[] VelocityMeanIn = new double[D];
            double[] VelocityMeanOut = new double[D];
            for (int d = 0; d < D; d++) {
                VelocityMeanIn[d] = cp.ParamsNeg[D + d];
                VelocityMeanOut[d] = cp.ParamsPos[D + d];
            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.n, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.n, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);

            double uJump = -M * ((1 / rhoA) - (1 / rhoB)) * cp.n[m_d];

            double flx = Lambda * uJump * 0.8;

            return -flx * (rhoA * vA - rhoB * vB);
        }



        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D),
                    VariableNames.HeatFlux0Vector(D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
            }
        }


    }


    public class ConvectionAtLevelSet_Consistency : EvaporationAtLevelSet {


        public ConvectionAtLevelSet_Consistency(int _d, int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            double vorZeichen, bool RescaleConti, ThermalParameters thermParams, double _Rint, double _sigma) {
            //double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.D = _D;
            this.m_d = _d;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_LsTrk = lsTrk;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= rhoA;
                scaleB /= rhoB;
            }

            this.kA = thermParams.k_A;
            this.kB = thermParams.k_B;
            this.hVapA = thermParams.hVap_A;
            this.Rint = _Rint;

            this.Tsat = thermParams.T_sat;
            this.sigma = _sigma;
            this.pc = thermParams.pc;
        }

        int m_d;
        double rhoA;
        double rhoB;

        double scaleA;
        double scaleB;



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

            //Console.WriteLine("mEvap - ConvectionAtLevelSet_Divergence: {0}", M);

            return M;

        }


        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.ParamsNeg.GetSubVector(2 * D, D + 3), cp.ParamsPos.GetSubVector(2 * D, D + 3), cp.n, cp.jCell);
            if (M == 0.0)
                return 0.0;


            double Ucentral = 0.0;
            for (int d = 0; d < D; d++) {
                Ucentral += 0.5 * (cp.ParamsNeg[d] + cp.ParamsPos[d]) * cp.n[d];
            }

            double uAxN = Ucentral * (-M * (1 / rhoA) * cp.n[m_d]);
            double uBxN = Ucentral * (-M * (1 / rhoB) * cp.n[m_d]);

            uAxN += -M * (1 / rhoA) * 0.5 * (U_Neg[0] + U_Pos[0]);
            uBxN += -M * (1 / rhoB) * 0.5 * (U_Neg[0] + U_Pos[0]);

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

            FlxNeg *= rhoA;
            FlxPos *= rhoB;

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
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }


        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D),
                    VariableNames.HeatFlux0Vector(D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
            }
        }


    }

    
}
