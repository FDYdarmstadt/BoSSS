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


        public ConvectionAtLevelSet_nonMaterialLLF(int _d, int _D, LevelSetTracker lsTrk, ThermalParameters thermParams, double _sigma) 
            : base(_D, thermParams, _sigma) {
                                 
            this.m_d = _d;
        }

        int m_d;


        public override double InnerEdgeForm(ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {



            double M = ComputeEvaporationMass(cp.Parameters_IN.GetSubVector(2 * m_D, m_D + 3), cp.Parameters_OUT.GetSubVector(2 * m_D, m_D + 3), cp.Normal, cp.jCellIn);
            if (M == 0.0)
                return 0.0;

            double[] VelocityMeanIn = new double[m_D];
            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = cp.Parameters_IN[m_D + d];
                VelocityMeanOut[d] = cp.Parameters_OUT[m_D + d];

            }

            double LambdaIn;
            double LambdaOut;

            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.Normal, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.Normal, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);


            double uJump = -M * ((1 / m_rhoA) - (1 / m_rhoB)) * cp.Normal[m_d];


            double flx = Lambda * uJump * 0.8;

            return -flx * (m_rhoA * vA - m_rhoB * vB);
        }



        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D),
                    VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.MassFluxExtension);
            }
        }


    }


    public class ConvectionAtLevelSet_Consistency : EvaporationAtLevelSet {


        public ConvectionAtLevelSet_Consistency(int _d, int _D, LevelSetTracker lsTrk,
            double vorZeichen, bool RescaleConti, ThermalParameters thermParams, double _sigma) 
            : base(_D, thermParams, _sigma) {

            this.m_d = _d;

            scaleA = vorZeichen;
            scaleB = vorZeichen;

            if (RescaleConti) {
                scaleA /= m_rhoA;
                scaleB /= m_rhoB;
            }

        }

        int m_d;

        double scaleA;
        double scaleB;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }



        public override double InnerEdgeForm(ref CommonParams cp,

            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {


            double M = ComputeEvaporationMass(cp.Parameters_IN.GetSubVector(2 * m_D, m_D + 3), cp.Parameters_OUT.GetSubVector(2 * m_D, m_D + 3), cp.Normal, cp.jCellIn);

            if (M == 0.0)
                return 0.0;

            double Ucentral = 0.0;
            for (int d = 0; d < m_D; d++) {
                Ucentral += 0.5 * (cp.Parameters_IN[d] + cp.Parameters_OUT[d]) * cp.Normal[d];
            }

            double uAxN = Ucentral * (-M * (1 / m_rhoA) * cp.Normal[m_d]);
            double uBxN = Ucentral * (-M * (1 / m_rhoB) * cp.Normal[m_d]);


            uAxN += -M * (1 / m_rhoA) * 0.5 * (U_Neg[0] + U_Pos[0]);
            uBxN += -M * (1 / m_rhoB) * 0.5 * (U_Neg[0] + U_Pos[0]);

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

            FlxNeg *= m_rhoA;
            FlxPos *= m_rhoB;

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
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D),
                    VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.MassFluxExtension);
            }
        }


    }


    public class ConvectionAtLevelSet_MovingMesh : EvaporationAtLevelSet {

        public ConvectionAtLevelSet_MovingMesh(int _d, int _D, ThermalParameters thermParams, double _sigma) 
            : base(_D, thermParams, _sigma) {

            this.m_d = _d;

        }

        int m_d;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }



        public override double InnerEdgeForm(ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            throw new NotImplementedException("TODO");

            //double M = ComputeEvaporationMass(cp.ParamsNeg, cp.ParamsPos, cp.n, cp.jCell);

            //double s;
            //if (hVapA > 0) {
            //    s = (M / rhoB) * cp.n[m_d] + U_Pos[0];
            //} else {
            //    s = (M / rhoA) * cp.n[m_d] + U_Neg[0];
            //}

            //double FlxNeg;
            //double FlxPos;
            //if (hVapA > 0) {
            //    FlxNeg = -s * U_Neg[0];
            //    FlxPos = -FlxNeg;
            //} else {
            //    FlxNeg = -s * U_Pos[0];
            //    FlxPos = -FlxNeg;
            //}

            //return FlxNeg * vA - FlxPos * vB;
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }
    }
}
