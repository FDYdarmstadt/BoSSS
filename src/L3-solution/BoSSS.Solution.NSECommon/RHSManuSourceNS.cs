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
    /// Auxillary class to compute a source term from a manufactured solution for the Navier-Stokes momentum equations.
    /// Current implementation only supports 2D flows.
    /// Current manufactured solutions used is T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
    /// See also ControlManuSol() control function.
    /// </summary>
    public class RHSManuSourceNS : BoSSS.Solution.Utils.LinearSource {

        double[] MolarMasses;
        string direction;
        double Reynolds;
        bool energyOK;
        bool speciesOK;
        /// <summary>
        /// Ctor.
        /// <param name="MolarMasses">Array of the molar masses of the fuel, oxidizer and products.</param>
        /// <param name="direction">Can be "x" or "y".</param>        
        /// </summary>
        public RHSManuSourceNS(double Reynolds, double[] MolarMasses, string direction, bool energyOK, bool speciesOK) {
            this.MolarMasses = MolarMasses;
            this.direction = direction;
            this.Reynolds = Reynolds;
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

        //Manufactured solution for T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = cos(x*y), v = cos(x*y), p = sin(x*y).
        protected override double Source(double[] x, double[] parameters, double[] U) {

            double p0 = 1.0;
            double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; 
            double x_ = x[0];
            double y_ = x[1];
            double alpha1 = 0.3;
            double alpha2 = 0.6;
            double alpha3 = 0.1;
            double[] Coefficients = new double[] { alpha1, alpha2, alpha3 };

            double ConvectionTerm;
            double ViscTerm;
            double PressureGradientTerm;

            if (direction == "x") {
                ConvectionTerm = p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Pow(Math.Cos(x_), 0.2e1) * Math.Sin(x_ * y_) * y_ - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Pow(Math.Cos(x_), 0.2e1) * (-alpha1 * Math.Sin(x_ * y_) * y_ / M1 - alpha2 * Math.Sin(x_ * y_) * y_ / M2 - alpha3 * Math.Sin(x_ * y_) * y_ / M3 + (alpha1 * Math.Sin(x_ * y_) * y_ + alpha2 * Math.Sin(x_ * y_) * y_ + alpha3 * Math.Sin(x_ * y_) * y_) / M4) - 0.2e1 * p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(x_) * Math.Sin(x_) + p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(x_) * Math.Cos(y_) * Math.Sin(x_ * y_) * x_ - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(x_) * Math.Cos(y_) * (-alpha1 * Math.Sin(x_ * y_) * x_ / M1 - alpha2 * Math.Sin(x_ * y_) * x_ / M2 - alpha3 * Math.Sin(x_ * y_) * x_ / M3 + (alpha1 * Math.Sin(x_ * y_) * x_ + alpha2 * Math.Sin(x_ * y_) * x_ + alpha3 * Math.Sin(x_ * y_) * x_) / M4) - p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(x_) * Math.Sin(y_); // conti, mom, energy and species
          
                ViscTerm = -0.4e1 / 0.3e1 * Math.Cos(x_) / Reynolds; // TODO?
                PressureGradientTerm = y_ * Math.Cos(x_ * y_); // OK

       
                if (!speciesOK) { 
                ConvectionTerm = p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Pow(Math.Cos(x_), 0.2e1) * y_ * Math.Sin(x_ * y_) - 0.2e1 * p0 / Math.Cos(x_ * y_) * Math.Cos(x_) * Math.Sin(x_) + p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(x_) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) - p0 / Math.Cos(x_ * y_) * Math.Cos(x_) * Math.Sin(y_); // conti, mom and energy
                }

                if (!energyOK && !speciesOK) {
                    ConvectionTerm = -0.20e1 * Math.Cos(x_) * Math.Sin(x_) - 0.10e1 * Math.Cos(x_) * Math.Sin(y_); // conti and mom
                }

            }
            else if (direction == "y") {
                ConvectionTerm = p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(x_) * Math.Cos(y_) * Math.Sin(x_ * y_) * y_ - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(x_) * Math.Cos(y_) * (-alpha1 * Math.Sin(x_ * y_) * y_ / M1 - alpha2 * Math.Sin(x_ * y_) * y_ / M2 - alpha3 * Math.Sin(x_ * y_) * y_ / M3 + (alpha1 * Math.Sin(x_ * y_) * y_ + alpha2 * Math.Sin(x_ * y_) * y_ + alpha3 * Math.Sin(x_ * y_) * y_) / M4) - p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(x_) * Math.Cos(y_) + p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Pow(Math.Cos(y_), 0.2e1) * Math.Sin(x_ * y_) * x_ - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Pow(Math.Cos(y_), 0.2e1) * (-alpha1 * Math.Sin(x_ * y_) * x_ / M1 - alpha2 * Math.Sin(x_ * y_) * x_ / M2 - alpha3 * Math.Sin(x_ * y_) * x_ / M3 + (alpha1 * Math.Sin(x_ * y_) * x_ + alpha2 * Math.Sin(x_ * y_) * x_ + alpha3 * Math.Sin(x_ * y_) * x_) / M4) - 0.2e1 * p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.1e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(y_) * Math.Sin(y_);

                ViscTerm = -0.4e1 / 0.3e1 * Math.Cos(y_) / Reynolds; // TODO? maple says its 8/3 and not 4/3. maybe the sign is wrong
                PressureGradientTerm = x_ * Math.Cos(x_ * y_);

           
                if (!speciesOK) {
                    ConvectionTerm = p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(x_) * Math.Cos(y_) * y_ * Math.Sin(x_ * y_) - p0 / Math.Cos(x_ * y_) * Math.Sin(x_) * Math.Cos(y_) + p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Pow(Math.Cos(y_), 0.2e1) * x_ * Math.Sin(x_ * y_) - 0.2e1 * p0 / Math.Cos(x_ * y_) * Math.Cos(y_) * Math.Sin(y_); // conti, mom and energy
                }

                if (!energyOK && !speciesOK) {
                    ConvectionTerm = -0.10e1 * Math.Cos(y_) * Math.Sin(x_) - 0.20e1 * Math.Cos(y_) * Math.Sin(y_);
                }

            }
            else
                throw new ArgumentException("Specified direction not supported");


            return -(ConvectionTerm + ViscTerm + PressureGradientTerm);
        }
    }
}
