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
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon.Operator.Convection {
    public class ConvectionAtIB : ILevelSetForm {
        public ConvectionAtIB(int _d, int _D, LevelSetTracker LsTrk, double _LFFA, IncompressibleBoundaryCondMap _bcmap, Func<double, double>[] _uLevSet, Func<double, double>[] _wLevSet, double particleRadius, double fluidDensity, bool UseMovingMesh){
            m_LsTrk = LsTrk;
            m_D = _D;
            m_d = _d;
            this.uLevSet = _uLevSet;
            this.wLevSet = _wLevSet;
            LFFA = _LFFA;
            this.pRadius = particleRadius;
            //varMode = _varMode;
            fDensity = fluidDensity;
            m_UseMovingMesh = UseMovingMesh;

            NegFlux = new LinearizedConvection(_D, _bcmap, _d);
            //NegFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, fluidDensity, 0, _LFFA, double.NaN, LsTrk);
            //NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"), null);
        }
        LevelSetTracker m_LsTrk;
        int m_D;
        int m_d;
        Func<double, double>[] uLevSet, wLevSet;
        double pRadius;
        double fDensity;
        double LFFA;
        bool m_UseMovingMesh;

        // Use Fluxes as in Bulk Convection
        LinearizedConvection NegFlux;


        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Velocity_d(m_d) }; }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        /*

        // Flux over interface
        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp,
            double[] U_Neg, double[] U_Pos, double[,] GradU_Neg, double[,] GradU_Pos) {

            double[] _uLevSet = new double[2];

            //_uLevSet[0] = (uLevSet[m_d])(cp.time, cp.x);

            for (int d = 0; d < m_D; d++) {
                _uLevSet[d] = (uLevSet[d])(cp.time);
            }

            double[] uLevSet_temp = new double[1];
            uLevSet_temp[0] = (uLevSet[m_d])(cp.time);

            BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
            inp.Parameters_IN = cp.ParamsNeg;
            inp.Normale = cp.n;
            inp.iEdge = int.MinValue;
            inp.GridDat = this.m_LsTrk.GridDat;
            inp.X = cp.x;
            inp.time = cp.time;
            //inp.jCellIn = cp.jCell;
            //inp.jCellOut = cp.jCell;

            inp.Parameters_OUT = new double[inp.Parameters_IN.Length];

            //Outer values for Velocity and VelocityMean
            for (int j = 0; j < m_D; j++) {
                inp.Parameters_OUT[j] = (uLevSet[j])(cp.time);
                // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                inp.Parameters_OUT[m_D + j] = 0;
            }

            //FlxNeg = -this.NegFlux.IEF(ref inp, U_Neg, uLevSet_temp);

            FlxNeg = 0;

            FlxPos = 0;


        }
        */

        public double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] _uLevSet = new double[2];



            BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
            inp.Parameters_IN = cp.ParamsNeg;
            inp.Normale = cp.n;
            inp.iEdge = int.MinValue;
            inp.GridDat = this.m_LsTrk.GridDat;
            inp.X = cp.x;
            inp.time = cp.time;
            inp.Parameters_OUT = new double[inp.Parameters_IN.Length];



            //for (int d = 0; d < m_D; d++)
            //{
            //    _uLevSet[d] = (uLevSet[d])(cp.time);
            //}

            double[] uLevSet_temp = new double[1];
            if (m_d == 0) {
                uLevSet_temp[0] = (uLevSet[0])(cp.time) + pRadius * wLevSet[0](cp.time) * -cp.n[1];
            } else { uLevSet_temp[0] = (uLevSet[1])(cp.time) + pRadius * wLevSet[0](cp.time) * cp.n[0]; }

            //Outer values for Velocity and VelocityMean
            inp.Parameters_OUT[0] = (uLevSet[0])(cp.time) + pRadius * wLevSet[0](cp.time) * -cp.n[1];
            inp.Parameters_OUT[1] = (uLevSet[1])(cp.time) + pRadius * wLevSet[0](cp.time) * cp.n[0];
            // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
            inp.Parameters_OUT[2] = 0;
            inp.Parameters_OUT[3] = 0;


            double FlxNeg;
            if (m_UseMovingMesh == true) {
                FlxNeg = 0;

                return FlxNeg;
            }

            FlxNeg = this.NegFlux.InnerEdgeForm(ref inp, U_Neg, uLevSet_temp, null, null, v_Neg, 0, null, null);

            //FlxNeg = this.NegFlux.InnerEdgeForm(ref inp, U_Neg, uLevSet_temp, Grad_uA, Grad_uB, v_Neg, v_Pos, Grad_vA, Grad_vB);

            return FlxNeg;

        }
    }
}
