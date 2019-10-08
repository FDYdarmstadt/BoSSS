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


    public class HeatConvectionAtLevelSet_LLF : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public HeatConvectionAtLevelSet_LLF(int _D, LevelSetTracker LsTrk, double _capA, double _capB, double _LFFA, double _LFFB, 
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

                U_NegFict = U_Pos;
                U_PosFict = U_Neg;
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
                double Tavg = Tsat; // 0.5 * (U_Neg[0] + Tsat);
                r += Tavg * (ParamsNeg[0] * cp.n[0] + ParamsNeg[1] * cp.n[1]);
                if (m_D == 3) {
                    r += Tavg * ParamsNeg[2] * cp.n[2];
                }

                // Calculate dissipative part
                // ==========================

                double[] VelocityMeanIn = new double[m_D];
                for (int d = 0; d < m_D; d++) {
                    VelocityMeanIn[d] = ParamsNeg[m_D + d];
                }

                double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.n, false);

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
                double Tavg = Tsat; // 0.5 * (Tsat +  U_Pos[0]);
                r += Tavg * (ParamsPos[0] * cp.n[0] + ParamsPos[1] * cp.n[1]);
                if (m_D == 3) {
                    r += Tavg * ParamsPos[2] * cp.n[2];
                }

                // Calculate dissipative part
                // ==========================

                double[] VelocityMeanOut = new double[m_D];
                for (int d = 0; d < m_D; d++) {
                    VelocityMeanOut[d] = ParamsPos[m_D + d];
                }

                double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.n, false);

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


    public class HeatConvectionAtLevelSet_WithEvaporation : EvaporationAtLevelSet {


        public HeatConvectionAtLevelSet_WithEvaporation(int _D, LevelSetTracker lsTrk, double _LFFA, double _LFFB,
            ThermalParameters thermParams, double _sigma)
            : base(_D, lsTrk, thermParams, _sigma) {

            this.LFFA = _LFFA;
            this.LFFB = _LFFB;

            this.capA = thermParams.c_A * thermParams.rho_A;
            this.capB = thermParams.c_B * thermParams.rho_B;

        }
 
        double capA;
        double capB;

        double LFFA;
        double LFFB;


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }


        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {

            U_NegFict = U_Pos;
            U_PosFict = U_Neg;
        }

        public override double LevelSetForm(ref Foundation.XDG.CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] U_NegFict, U_PosFict;
            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);


            double s = ComputeInterfaceNormalVelocity(cp.ParamsNeg.GetSubVector(2* m_D, m_D + 3), cp.ParamsPos.GetSubVector(2*m_D, m_D + 3), cp.n, cp.jCell);
            //Console.WriteLine("interfaceNormalVelocity = {0}", s);


            // Flux for negative side
            // ======================
            double FlxNeg = 0.0;

            // Calculate central part
            FlxNeg += m_Tsat * (ParamsNeg[0] * cp.n[0] + ParamsNeg[1] * cp.n[1]);
            if (m_D == 3) {
                FlxNeg += m_Tsat * ParamsNeg[2] * cp.n[2];
            }
            //FlxNeg -= Tsat * s;

            // Calculate dissipative part
            double[] VelocityMeanIn = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanIn[d] = ParamsNeg[m_D + d];
            }

            //double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, cp.n, false);
            double VA_n = 0.0;
            for (int d = 0; d < VelocityMeanIn.Length; d++)
                VA_n += VelocityMeanIn[d] * cp.n[d];

            double LambdaIn = Math.Abs(VA_n - s);

            double uJumpA = U_Neg[0] - m_Tsat;

            FlxNeg += LambdaIn * uJumpA * LFFA;

            FlxNeg *= capA;


            // Flux for positive side
            // ======================
            double FlxPos = 0.0;

            // Calculate central part
            FlxPos += m_Tsat * (ParamsPos[0] * cp.n[0] + ParamsPos[1] * cp.n[1]);
            if (m_D == 3) {
                FlxPos += m_Tsat * ParamsPos[2] * cp.n[2];
            }
            //FlxPos -= Tsat * s;

            // Calculate dissipative part
            double[] VelocityMeanOut = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                VelocityMeanOut[d] = ParamsPos[m_D + d];
            }

            //double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, cp.n, false);
            double VB_n = 0.0;
            for (int d = 0; d < VelocityMeanOut.Length; d++)
                VB_n += VelocityMeanOut[d] * cp.n[d];

            double LambdaOut = Math.Abs(VB_n - s);

            double uJumpB = m_Tsat - U_Pos[0];

            FlxPos += LambdaOut * uJumpB * LFFB;

            FlxPos *= capB;


            return FlxNeg * vA - FlxPos * vB;
        }


        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D),
                    VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
            }
        }


    }


    public class HeatConvectionAtLevelSet_Upwind : EvaporationAtLevelSet {


        public HeatConvectionAtLevelSet_Upwind(int _D,  LevelSetTracker lsTrk, double _capA, double _capB, 
            ThermalParameters thermParams, bool _movingmesh, bool _DiriCond, double _Tsat, double _sigma)
            : base(_D, lsTrk, thermParams, _sigma) {


            this.capA = thermParams.c_A * thermParams.rho_A;
            this.capB = thermParams.c_B * thermParams.rho_B;

            movingmesh = _movingmesh;

            DirichletCond = _DiriCond;

        }

        double capA;
        double capB;

        bool movingmesh;

        bool DirichletCond;        



        public override double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            throw new NotImplementedException("check for consistency");

            //double cNeg = 0.0;
            //double cPos = 0.0;
            //for (int d = 0; d < D; d++) {
            //    cNeg += inp.ParamsNeg[d] * inp.n[d];
            //    cPos += inp.ParamsPos[d] * inp.n[d];
            //}


            //double s = ComputeInterfaceNormalVelocity(inp.ParamsNeg, inp.ParamsPos, inp.n, inp.jCell);
            //double RelSpeedNeg = cNeg - s;
            //double RelSpeedPos = cPos - s;

            //double FluxNeg;
            //if (RelSpeedNeg >= 0) { // UP-wind with respect to relative speed
            //    FluxNeg = RelSpeedNeg * m_Tsat;
            //} else {
            //    FluxNeg = RelSpeedNeg * m_Tsat;
            //}
            //FluxNeg *= capA;

            //double FluxPos;
            //if (RelSpeedPos >= 0) { // UP-wind with respect to relative speed
            //    FluxPos = RelSpeedPos * m_Tsat;
            //} else {
            //    FluxPos = RelSpeedPos * m_Tsat;
            //}
            //FluxPos *= capB;

            //return FluxNeg * vA - FluxPos * vB;
        }





        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Temperature };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
            }
        }


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }
}
