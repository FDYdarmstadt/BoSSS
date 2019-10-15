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

        protected double m_hVap;

        // for micro regions
        protected double m_sigma;
        protected double m_pc;
        protected double m_fc;
        protected double m_Rc;
        protected double m_Tsat;

        //protected double rho;     // density of liquid phase 
        protected double m_rhoA;
        protected double m_rhoB;

        protected int m_D;


        public EvaporationAtLevelSet(int _D, LevelSetTracker _LsTrk,  ThermalParameters thermParams, double _sigma) {

            this.m_D = _D;
            this.m_LsTrk = _LsTrk;

            this.m_hVap = thermParams.hVap;

            this.m_rhoA = thermParams.rho_A;
            this.m_rhoB = thermParams.rho_B;

            this.m_pc = thermParams.pc;
            this.m_fc = thermParams.fc;
            this.m_Rc = thermParams.Rc;
            this.m_Tsat = thermParams.T_sat;

            this.m_sigma = _sigma;
                
        }


        private double ComputeHeatFlux_Macro(double[] HeatFlux_A, double[] HeatFlux_B, double[] n) {

            double qEvap = 0.0;
            for (int d = 0; d < m_D; d++)
                qEvap += (HeatFlux_B[d] - HeatFlux_A[d]) * n[d];

            return qEvap;
        }

        private double ComputeHeatFlux_Micro(double T_A, double T_B, double curv, double p_disp) {

            double pc = m_sigma * curv + p_disp;      // augmented capillary pressure (without nonlinear evaporative masss part)

            double Rint = 0.0;
            double TintMin = 0.0;
            double qEvap = 0.0;
            if (m_rhoA > m_rhoB) {
                Rint = ((2.0 - m_fc) / (2 * m_fc)) * m_Tsat * Math.Sqrt(2 * Math.PI * m_Rc * m_Tsat) / (m_rhoB * m_hVap.Pow2());
                TintMin = m_Tsat * (1 + (pc / (m_hVap * m_rhoA)));
                //if (T_A > TintMin)
                    qEvap = (T_A - TintMin) / Rint;
            } else {
                Rint = ((2.0 - m_fc) / (2 * m_fc)) * m_Tsat * Math.Sqrt(2 * Math.PI * m_Rc * m_Tsat) / (m_rhoA * m_hVap.Pow2());
                TintMin = m_Tsat * (1 + (pc / (m_hVap * m_rhoB)));
                //if (T_B > TintMin)
                    qEvap = (T_B - TintMin) / Rint;
            }

            return qEvap;
        }

        protected double ComputeHeatFlux(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            if (m_hVap == 0.0)
                return 0.0;

            if (MEvapIsPrescribd)
                return prescrbMEvap * m_hVap;

            double qEvap = 0.0;
            if (evapMicroRegion[jCell]) {
                Debug.Assert(paramsPos[m_D + 1] == paramsNeg[m_D + 1], "curvature must be continuous across interface");
                Debug.Assert(paramsPos[m_D + 2] == paramsNeg[m_D + 2], "disjoining pressure must be continuous across interface");

                qEvap = ComputeHeatFlux_Micro(paramsNeg[m_D], paramsPos[m_D], paramsNeg[m_D + 1], paramsNeg[m_D + 2]);
            } else {
                qEvap = ComputeHeatFlux_Macro(paramsNeg.GetSubVector(0, m_D), paramsPos.GetSubVector(0, m_D), N);
            }

            return qEvap;

        }

        protected double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double M = qEvap / m_hVap;

            return M;
        }



        protected double ComputeInterfaceNormalVelocity(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            double M = qEvap / m_hVap;

            double sNeg = 0.0;
            for (int d = 0; d < m_D; d++)
                sNeg += (paramsNeg[d] + (M / m_rhoA)) * N[d];

            double sPos = 0.0;
            for (int d = 0; d < m_D; d++)
                sPos += (paramsPos[d] + (M / m_rhoB)) * N[d];

            double s = (m_rhoA * sNeg + m_rhoB * sPos) / (m_rhoA + m_rhoB);     // density averaged, corresponding to the mean evo velocity 

            return s;

        }


        public abstract double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB);


        protected LevelSetTracker m_LsTrk;
        protected BitArray evapMicroRegion;

        bool MEvapIsPrescribd = false;
        double prescrbMEvap;

        public virtual void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            if (csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

            if (csA.UserDefinedValues.Keys.Contains("prescribedMassflux")) {
                MEvapIsPrescribd = true;
                prescrbMEvap = (double)csA.UserDefinedValues["prescribedMassflux"];
            }

        }


        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
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

}
