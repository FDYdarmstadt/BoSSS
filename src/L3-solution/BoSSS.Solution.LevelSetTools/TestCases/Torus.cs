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
using ilPSP;

namespace BoSSS.Solution.LevelSetTools.TestCases {
    public class Torus {

        double bigR;
        double smallR;

        public Torus(double bigR, double smallR) {
            this.bigR = bigR;
            this.smallR = smallR;
        }

        /// <summary>
        /// Signed-Distance Level-Set Function for a Torus
        /// </summary>
        /// <param name="X">Cartesian Coordinates</param>
        /// <returns>Level-Set Value at X</returns>
        public double SignedDistance(double[] X) {
            return Math.Sqrt(X[2] * X[2] + (Math.Sqrt(X[0].Pow2() + X[1].Pow2()) - bigR).Pow2()) - smallR;
        }

        public double ExtensionAnalytical(double[] X, Func<double[],double> initial) {
            double[] R = CartesianToTorus(X);
            // project onto Surface , i.e. set small Radius to the value on the surface
            R[0] = smallR;
            // evaluate the initial Value at the interface
            double[] XSurf = TorusToCartesian(R);
            return initial(XSurf);
        }


        /// <summary>
        /// Transforms Cartesian to Torus Coordinates
        /// Index 0: radius
        /// Index 1: phi - toroidal angle (along big R)
        /// INdex 2: azimuthal angle (along small r)
        /// </summary>
        /// <param name="X"></param>
        /// <returns></returns>
        public double[] CartesianToTorus(double[] X) {
            double phi, theta;
            phi = Math.Atan2(X[1], X[0]);
            double xi = Math.Sqrt(X[0] * X[0] + X[1] * X[1]) - bigR;
            theta = Math.Atan2(X[2], xi);
            double rho = Math.Sqrt(xi * xi + X[2] * X[2]);
            return new double[] { rho, phi, theta };
        }

        /// <summary>
        /// Transforms Torus to Cartesian Coordinates
        /// Index 0: radius
        /// Index 1: phi - toroidal angle (along big R)
        /// INdex 2: azimuthal angle (along small r)
        /// </summary>
        /// <param name="X"></param>
        /// <returns></returns>
        public double[] TorusToCartesian(double[] R) {
            double x, y, z;
            double r        = R[0];
            double phi      = R[1];
            double theta    = R[2];
            double rCosTheta = r*Math.Cos(theta);
            x = (bigR + rCosTheta) * Math.Cos(phi);
            y = (bigR + rCosTheta) * Math.Sin(phi);
            z = r * Math.Sin(theta);
            return new double[] { x,y,z };
        }

    }
}
