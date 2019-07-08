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
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System.Collections;

namespace BoSSS.Solution.XheatCommon {


    public abstract class EvaporationAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        protected double hVapA;   // for the identification of the liquid phase
        protected double Rint;
        protected double Tsat;
        protected double rho;     // density of liquid phase 
        protected double sigma;
        protected double pc;

        protected double kA;
        protected double kB;

        protected int D;

        private double ComputeHeatFlux_Macro(double[] GradT_A, double[] GradT_B, double[] n) {

            double hVap = 0.0;
            double qEvap = 0.0;
            if (hVapA > 0) {
                hVap = hVapA;
                for (int d = 0; d < D; d++)
                    qEvap += (kA * GradT_A[d] - kB * GradT_B[d]) * n[d];
            } else {
                hVap = -hVapA;
                for (int d = 0; d < D; d++)
                    qEvap += (kB * GradT_B[d] - kA * GradT_A[d]) * n[d];
            }

            return qEvap;
        }

        private double ComputeHeatFlux_Micro(double T_A, double T_B, double curv, double p_disp) {

            double pc0 = (pc < 0.0) ? sigma * curv + p_disp : pc;      // augmented capillary pressure (without nonlinear evaporative masss part)

            double TintMin = 0.0;
            double hVap = 0.0;
            double qEvap = 0.0;
            if (hVapA > 0) {
                hVap = hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rho)));
                if (T_A > TintMin)
                    qEvap = -(T_A - TintMin) / Rint;
            } else {
                hVap = -hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rho)));
                if (T_B > TintMin)
                    qEvap = (T_B - TintMin) / Rint;
            }

            return qEvap;
        }

        protected double ComputeHeatFlux(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            if (hVapA == 0.0)
                return 0.0;

            double qEvap = 0.0;
            if (evapMicroRegion[jCell]) {
                Debug.Assert(paramsPos[D + 1] == paramsNeg[D + 1], "curvature must be continuous across interface");
                Debug.Assert(paramsPos[D + 2] == paramsNeg[D + 2], "disjoining pressure must be continuous across interface");

                qEvap = ComputeHeatFlux_Micro(paramsNeg[D], paramsPos[D], paramsNeg[D + 1], paramsNeg[D + 2]);
            } else {
                qEvap = ComputeHeatFlux_Macro(paramsNeg.GetSubVector(0, D), paramsPos.GetSubVector(0, D), N);
            }

            return qEvap;

        }


        public abstract double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB);


        protected LevelSetTracker m_LsTrk;
        BitArray evapMicroRegion;

        public virtual void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }


        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Temperature0Gradient(D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
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

        public virtual TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }


    }


    public class HeatFluxAtLevelSet : EvaporationAtLevelSet {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public HeatFluxAtLevelSet(LevelSetTracker LsTrk, double _rho, ThermalParameters thermParams, double _Rint, double _sigma) {
            m_LsTrk = LsTrk;

            this.rho = _rho;

            this.kA = thermParams.k_A;
            this.kB = thermParams.k_B;
            this.hVapA = thermParams.hVap_A;
            this.Rint = _Rint;

            this.Tsat = thermParams.T_sat;
            this.sigma = _sigma;
            this.pc = thermParams.pc;

            this.D = LsTrk.GridDat.SpatialDimension;
        }


        public override double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double qEvap = ComputeHeatFlux(cp.ParamsNeg, cp.ParamsPos, cp.n, cp.jCell);
            //Console.WriteLine("qEvap - HeatFluxAtLevelSet: {0}", qEvap);
            if (qEvap == 0.0)
                return 0.0;


            double FlxNeg = -0.5 * qEvap;
            double FlxPos = +0.5 * qEvap;

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }


    }



}
