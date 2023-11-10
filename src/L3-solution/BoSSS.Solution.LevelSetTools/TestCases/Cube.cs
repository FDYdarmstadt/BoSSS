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
    public class Cube {

        double L;

        public Cube(double length) {
            this.L = length;
        }

        /// <summary>
        /// Signed-Distance Level-Set Function for a Torus
        /// </summary>
        /// <param name="X">Cartesian Coordinates</param>
        /// <returns>Level-Set Value at X</returns>
        public double SignedDistance(double[] X) {
            double dx = Math.Max(Math.Abs(X[0]) - L / 2, 0.0);
            double dy = Math.Max(Math.Abs(X[1]) - L / 2, 0.0);
            double dz = Math.Max(Math.Abs(X[2]) - L / 2, 0.0);

            double distanceInside = Math.Min(Math.Max(dx, Math.Max(dy, dz)), 0.0);
            double distanceOutside = Math.Sqrt(dx * dx + dy * dy + dz * dz);

            return distanceInside + distanceOutside;
        }


        /// <summary>
        /// Signed-Distance 2D Level-Set Function for a Torus
        /// </summary>
        /// <param name="X">Cartesian Coordinates</param>
        /// <returns>Level-Set Value at X</returns>
        public double SignedDistance2D(double[] X) {
            double dx = Math.Max(Math.Abs(X[0]) - L / 2, 0.0);
            double dy = Math.Max(Math.Abs(X[1]) - L / 2, 0.0);

            double distanceInside = Math.Min(dx, dy);
            double distanceOutside = Math.Sqrt(dx * dx + dy * dy);

            return distanceInside + distanceOutside;
        }
    }
}
