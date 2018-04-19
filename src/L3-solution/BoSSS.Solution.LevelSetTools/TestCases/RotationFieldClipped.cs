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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;


namespace BoSSS.Solution.LevelSetTools.TestCases{

    /// <summary>
    /// Velocity field for a square box, where the Velocity is a solid body rotation in the interior,
    /// Parallel to the walls and zero in the corners
    /// </summary>
    public class RotationFieldClipped{
        double boundary;
        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="boundary">
        /// size of the box
        /// </param>
        public RotationFieldClipped(double boundary) {
            this.boundary = boundary;
        }

        private bool Interior(double[] X){
            return Math.Sqrt(X[0] * X[0] + X[1] * X[1]) < boundary;
        }

        /// <summary>
        /// A Function, which is 0 at 0 and 1 at 1
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double Smoothingfunction(double x){
            return x;
        }

        /// <summary>
        /// Velocity-Component in X-Direction
        /// </summary>
        /// <param name="X">
        /// Cartesian Coordinates
        /// </param>
        /// <returns>
        /// Value at Coordinates
        /// </returns>
        public double VelocityX(double[] X){

            if (Interior(X)) {
                return -X[1];
            }
            else {
                double phi = Math.Atan2(X[1], X[0]);
                double x0 = boundary * Math.Cos(phi);
                return -X[1] - Smoothingfunction((X[0] - x0).Abs() / (boundary - Math.Abs(x0))) * (-X[1]);
            }

        }

        /// <summary>
        /// Velocity-Component in Y-Direction
        /// </summary>
        /// <param name="X">
        /// Cartesian Coordinates
        /// </param>
        /// <returns>
        /// Value at Coordinates
        /// </returns>
        public double VelocityY(double[] X) {

            if (Interior(X)) {
                return X[0];
            }
            else {
                double phi = Math.Atan2(X[1], X[0]);
                double y0 = boundary * Math.Sin(phi);
                return X[0] - Smoothingfunction((X[1] - y0).Abs() / (boundary - Math.Abs(y0)))  * (X[0]);
            }

        }
    }

    
}
