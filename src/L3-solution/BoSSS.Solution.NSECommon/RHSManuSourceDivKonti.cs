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
using BoSSS.Solution.Utils;


namespace BoSSS.Solution.NSECommon {


    /// <summary>
    /// Auxillary class to compute a source term from a manufactured solution for the continuity equation.
    /// Current implementation only supports 2D flows.
    /// Current manufactured solutions used is T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
    /// See also ControlManuSol() control function.
    /// </summary>
    public class RHSManuSourceDivKonti : BoSSS.Solution.Utils.LinearSource {

        double ReynoldsNumber;
        double[] MolarMasses;
        bool energyOK;
        bool speciesOK;
        /// <summary>
        /// Ctor.
        /// </summary>
        public RHSManuSourceDivKonti(double Reynolds, double[] MolarMasses, bool energyOK, bool speciesOK) {
            this.ReynoldsNumber = Reynolds;
            this.MolarMasses = MolarMasses;
            this.energyOK = energyOK;
            this.speciesOK = speciesOK;
        }

        /// <summary>
        /// None
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return new string[0]; } 
        }

        /// <summary>
        /// None
        /// </summary>
        public override IList<string> ParameterOrdering {
            get { return null; } 
        }

        //Manufactured solution for T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
        protected override double Source(double[] x, double[] parameters, double[] U) {

            double x_ = x[0];
            double y_ = x[1];
            double p0 = 1.0;
            double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3];
            double alpha1 = 0.3;
            double alpha2 = 0.6;
            double alpha3 = 0.1;
            double[] Coefficients = new double[] { alpha1, alpha2, alpha3 };

            double man1 =  -1*( -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(x_) * Math.Sin(x_ * y_) * y_ + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(x_) * (-alpha1 * Math.Sin(x_ * y_) * y_ / M1 - alpha2 * Math.Sin(x_ * y_) * y_ / M2 - alpha3 * Math.Sin(x_ * y_) * y_ / M3 + (alpha1 * Math.Sin(x_ * y_) * y_ + alpha2 * Math.Sin(x_ * y_) * y_ + alpha3 * Math.Sin(x_ * y_) * y_) / M4) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(x_) - p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(y_) * Math.Sin(x_ * y_) * x_ + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(y_) * (-alpha1 * Math.Sin(x_ * y_) * x_ / M1 - alpha2 * Math.Sin(x_ * y_) * x_ / M2 - alpha3 * Math.Sin(x_ * y_) * x_ / M3 + (alpha1 * Math.Sin(x_ * y_) * x_ + alpha2 * Math.Sin(x_ * y_) * x_ + alpha3 * Math.Sin(x_ * y_) * x_) / M4) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(y_) );
           
        
            if (!speciesOK) { // conti, mom  and energy equations
                man1 = -1*(-p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(x_) - p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(y_));
            }
            if (!energyOK && !speciesOK) { // conti and mom equations
                man1 = -1 * (Math.Sin(x_) + Math.Sin(y_));
            }

            return man1;

        }
    }
}
