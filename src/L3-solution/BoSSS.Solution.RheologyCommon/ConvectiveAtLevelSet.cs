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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Solution.RheologyCommon {
    /// <summary>
    /// Volume integral of viscosity part of constitutive equations.
    /// </summary>
    public class ConvectiveAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        protected double WeissenbergA;
        protected double WeissenbergB;
        protected double[] pen1;
        protected double alpha;

        /// <summary>
        /// Initialize viscosity part
        /// </summary>
        public ConvectiveAtLevelSet(LevelSetTracker lstrk, int Component, double Weissenberg_a, double Weissenberg_b, double _alpha) {
            this.m_LsTrk = lstrk;
            this.Component = Component;
            this.WeissenbergA = Weissenberg_a;
            this.WeissenbergB = Weissenberg_b;
            this.alpha = _alpha;
        }

        /// <summary>
        /// Ordering of the dependencies
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                switch (Component) {
                    case 0:
                        return new string[] { VariableNames.StressXX };
                    case 1:
                        return new string[] { VariableNames.StressXY };
                    case 2:
                        return new string[] { VariableNames.StressYY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// Ordering of the parameters
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return VariableNames.Velocity0Vector(2);
            }
        }


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParams inp,
            double[] TA, double[] TB, double[,] Grad_TA, double[,] Grad_TB,
            double VA, double VB, double[] Grad_vA, double[] Grad_vB) {
            double[] Normale = inp.Normal;

            //Flux In
            double res1 = 0;
            double flxIn = 0;
            double n_u1 = 0;

            for (int d = 0; d < 2; d++) {
                n_u1 += Normale[d] * 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]);
            }

            double factor;
            if (n_u1 < 0)
                factor = this.alpha;
            else
                factor = 1.0 - this.alpha;

            res1 += factor * n_u1 * (TA[0] - TB[0]);
            flxIn = WeissenbergA * res1 * VA;


            //Flux out
            double res2 = 0;
            double flxOut = 0;
            double n_u2 = 0;

            Normale.ScaleV(-1.0);
            for (int d = 0; d < 2; d++) {
                n_u2 += Normale[d] * 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]);
            }

            double factor2;
            if (n_u2 < 0)
                factor2 = this.alpha;
            else
                factor2 = 1.0 - this.alpha;

            res2 += factor2 * n_u2 * (TA[0] - TB[0]);
            flxOut = WeissenbergB * res2 * VB;

            return flxIn + flxOut;
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            if (csA.UserDefinedValues.Keys.Contains("Weissenbergnumber")) {
                WeissenbergA = (double)csA.UserDefinedValues["Weissenbergnumber"];
                //Console.WriteLine("Weissenbergnumber = {0}", m_Weissenberg);
            }

            if (csB.UserDefinedValues.Keys.Contains("Weissenbergnumber")) {
                WeissenbergB = (double)csB.UserDefinedValues["Weissenbergnumber"];
                //Console.WriteLine("Weissenbergnumber = {0}", m_Weissenberg);
            }
        }

        //private static bool rem = true;

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }
    }
}
