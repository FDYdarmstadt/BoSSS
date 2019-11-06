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
    public class RHSManuSourceDivKonti : IVolumeForm
    {
        double ReynoldsNumber;
        double[] MolarMasses;
        PhysicsMode physicsMode;
        double phystime;
        bool unsteady;



        /// <summary>
        /// Ctor.
        /// </summary>
        public RHSManuSourceDivKonti(double Reynolds, double[] MolarMasses, PhysicsMode physicsMode, double phystime, bool unsteady) {
            this.ReynoldsNumber = Reynolds;
            this.MolarMasses = MolarMasses;
            this.physicsMode = physicsMode;
            this.phystime = phystime;
            this.unsteady = unsteady;

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
 
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            //    throw new NotImplementedException();
            //}

            ////Manufactured solution for T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
            //protected override double Source(double[] x, double[] parameters, double[] U) {
            double[] x = cpv.Xglobal;
            double x_ = x[0];
            double y_ = x[1];
            double t_ = cpv.time;
            double p0 = 1.0;// ThermodynamicPressure.GetMeanValue(3);

            double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3];
            //double alpha1 = 0.3;
            //double alpha2 = 0.6;
            //double alpha3 = 0.1;
            double man1;


            //double conti = -p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * Math.Cos(x_ * t_) * y_ * t_ * Math.Sin(x_ * y_ * t_) + p0 / Math.Cos(x_ * y_ * t_) * t_ * Math.Sin(x_ * t_) - p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * Math.Cos(y_ * t_) * x_ * t_ * Math.Sin(x_ * y_ * t_) + p0 / Math.Cos(x_ * y_ * t_) * t_ * Math.Sin(y_ * t_);

            //return -1 * conti;















            /////////////////////////////////////////////////////////////////
            if (unsteady) {
                switch (physicsMode) {
                    case PhysicsMode.LowMach:
                        double dRhodt = p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * x_ * y_ * Math.Sin(x_ * y_ * t_);
                        double dRhoUdx = -p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * Math.Cos(x_ * t_) * y_ * t_ * Math.Sin(x_ * y_ * t_) + p0 / Math.Cos(x_ * y_ * t_) * t_ * Math.Sin(x_ * t_);
                        double dRhoVdy = -p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * Math.Cos(y_ * t_) * x_ * t_ * Math.Sin(x_ * y_ * t_) + p0 / Math.Cos(x_ * y_ * t_) * t_ * Math.Sin(y_ * t_);
                        man1 = -1 * (dRhodt*0 + dRhoUdx + dRhoVdy);
                        break;
                    case PhysicsMode.Combustion:
                        throw new NotImplementedException("TODO");
                    default:
                        throw new NotImplementedException("should not happen");
                }
            }
            else {
                switch (physicsMode) {
                    case PhysicsMode.LowMach:                    

                        double dRhodt = 0;
                        double dRhoUdx = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(x_);
                        double dRhoVdy = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(y_);
                        man1 = -1 * (dRhodt + dRhoUdx + dRhoVdy);
                        break;
                    case PhysicsMode.Combustion:
                        throw new NotImplementedException("TODO");
                    default:
                        throw new NotImplementedException("should not happen");
                }
            }



            //if (unsteady) {
            //    switch (physicsMode) {
            //        case PhysicsMode.LowMach:
            //            man1 = -1 * (p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * x_ * y_ * Math.Sin(x_ * y_ * t_) - p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * Math.Cos(x_ * t_) * y_ * t_ * Math.Sin(x_ * y_ * t_) + p0 / Math.Cos(x_ * y_ * t_) * t_ * Math.Sin(x_ * t_) - p0 * Math.Pow(Math.Cos(x_ * y_ * t_), -0.2e1) * Math.Cos(y_ * t_) * x_ * t_ * Math.Sin(x_ * y_ * t_) + p0 / Math.Cos(x_ * y_ * t_) * t_ * Math.Sin(y_ * t_));
            //            break;
            //        case PhysicsMode.Combustion:
            //            man1 = 0.0; //TODO
            //            break;
            //        default:
            //            throw new NotImplementedException("should not happen");
            //    }
            //}
            //else {
            //    switch (physicsMode) {
            //        case PhysicsMode.LowMach:
            //            man1 = -1 * (-p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(x_) - p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Sin(y_));// conti, momentum and energy
            //            break;
            //        case PhysicsMode.Combustion:
            //            man1 = -1 * (-p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_)) / M4) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(x_) - p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_)) / M4) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_)) / M4) * Math.Sin(y_)); // conti, momentum, energy and species   
            //            break;
            //        default:
            //            throw new NotImplementedException("should not happen");
            //    }
            //}







            return man1*V;
        }
    }
}
