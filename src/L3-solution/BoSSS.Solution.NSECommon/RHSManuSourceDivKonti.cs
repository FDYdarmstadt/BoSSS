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
    //public class RHSManuSourceDivKonti : BoSSS.Solution.Utils.LinearSource {
    public class RHSManuSourceDivKonti : IVolumeForm, ISupportsJacobianComponent {
        double ReynoldsNumber;
        double[] MolarMasses;
        PhysicsMode physicsMode;

        /// <summary>
        /// Ctor.
        /// </summary>
        public RHSManuSourceDivKonti(double Reynolds, double[] MolarMasses, PhysicsMode physicsMode) {
            this.ReynoldsNumber = Reynolds;
            this.MolarMasses = MolarMasses;
            this.physicsMode = physicsMode;
        }


        /// <summary>
        /// None
        /// </summary>
        public IList<string> ArgumentOrdering {
            get { return new string[0]; }
        }

        /// <summary>
        /// None
        /// </summary>
        public IList<string> ParameterOrdering {
            get { return null; }
        }


        /// <summary>
        /// None
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.AllOn;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            ////Manufactured solution for T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
            double[] x = cpv.Xglobal;
            double x_ = x[0];
            double y_ = x[1];
            double p0 = 1.0; 

            double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
            double alpha1 = 0.3;
            double alpha2 = 0.4;
            double alpha3 = 0.1;
            double alpha4 = 0.2;

            double man1;

            double dRhoUdx;
            double dRhoVdy;

            switch(physicsMode) {
                case PhysicsMode.Incompressible:
                    dRhoUdx = 0.10e1 * Math.Sin(x_);
                    dRhoVdy = 0.10e1 * Math.Sin(y_);
                    man1 = -1 * (dRhoUdx + dRhoVdy);
                    break;
                case PhysicsMode.LowMach:      
                    dRhoUdx = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(x_);
                    dRhoVdy = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(y_);
                    man1 = -1 * ( dRhoUdx + dRhoVdy);
                    break;
                case PhysicsMode.Combustion:
                    dRhoUdx = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_);

                    dRhoVdy = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                    man1 = -1 * ( dRhoUdx + dRhoVdy);
                    break;

                default:
                    throw new NotImplementedException("should not happen");
            }





            return man1 * V;
        }
    }
}
