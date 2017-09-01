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
using System.Threading.Tasks;

namespace BoSSS.Solution.Timestepping {



    /// <summary>
    /// Coefficients for Backward-Differentiation-Formulas (BDF) and Crank-Nicolson.
    /// </summary>
    public class BDFSchemeCoeffs {
               

        /// <summary>
        /// Empty constructor.
        /// </summary>
        public BDFSchemeCoeffs() {
        }

        /// <summary>
        /// Returns a Backward-Differentiation-Formula of order <paramref name="i"/>.
        /// </summary>
        static public BDFSchemeCoeffs BDF(int i) {
            BDFSchemeCoeffs R = new BDFSchemeCoeffs();

            switch (i) {
                
                case 1:
                R.theta1 = 1.0;
                R.theta0 = 0.0;
                R.beta = new[] { 1.0 };
                break;

                case 2:
                R.theta1 = 2.0 / 3.0;
                R.theta0 = 0.0;
                R.beta = new[] { 4.0 / 3.0, -1.0 / 3.0 };
                break;

                case 3:
                R.theta1 = 6.0 / 11.0;
                R.theta0 = 0.0;
                R.beta = new[] { 18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0 };
                break;

                case 4:
                R.theta1 = 12.0 / 25.0;
                R.theta0 = 0.0;
                R.beta = new[] { 48.0 / 25.0, -36.0 / 25.0, 16.0 / 25, -3.0 / 25.0 };
                break;

                case 5:
                R.theta1 = 60.0 / 137.0;
                R.theta0 = 0.0;
                R.beta = new[] { 300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0 };
                break;

                case 6:
                R.theta1 = 60.0 / 147.0;
                R.theta0 = 0.0;
                R.beta = new[] { 360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0 };
                break;

                default:
                throw new NotImplementedException("BDF scheme of order " + i + " is not supported.");
            }

            if (Math.Abs(R.beta.Sum() - 1.0) > 1.0e-10)
                throw new ApplicationException();

            return R;
        }

        /// <summary>
        /// Returns an Explicit-Euler formula.
        /// </summary>
        static public BDFSchemeCoeffs ExplicitEuler() {
            BDFSchemeCoeffs R = new BDFSchemeCoeffs();

            R.theta1 = 0.0;
            R.theta0 = 1.0;
            R.beta = new[] { 1.0 };

            if (Math.Abs(R.beta.Sum() - 1.0) > 1.0e-10)
                throw new ApplicationException();

            return R;
        }

        /// <summary>
        /// Returns a Crank-Nicolson formula.
        /// </summary>
        static public BDFSchemeCoeffs CrankNicolson() {
            BDFSchemeCoeffs R = new BDFSchemeCoeffs();

            R.theta1 = 0.5;
            R.theta0 = 0.5;
            R.beta = new[] { 1.0 };

            if (Math.Abs(R.beta.Sum() - 1.0) > 1.0e-10)
                throw new ApplicationException();

            return R;
        }


        /// <summary>
        /// Factor for the operator at t^1, i.e. 
        /// 0: explicit, 0.5: Crank-Nicoloson, 1: implicit;
        /// </summary>
        public double theta1;

        /// <summary>
        /// Factor for the operator at t^0, i.e. 
        /// 1: explicit, 0.5: Crank-Nicoloson, 0: implicit;
        /// </summary>
        public double theta0;

        /// <summary>
        /// BDF-factors for u[0], u[-1], u[-2], ...
        /// </summary>
        public double[] beta;

        /// <summary>
        /// Number of BDF-stages
        /// </summary>
        public int S {
            get {
                return beta.Length;
            }
        }
    }
}
