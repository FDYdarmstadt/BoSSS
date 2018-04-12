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


namespace BoSSS.Solution.LevelSetTools.TestCases
{
    public class Flower
    {
        public Flower(double Radius, int subdivision, double Velocity = 0.5){
            this.subdivision = subdivision;
            this.Velocity = Velocity;
            this.Radius = Radius;
        }

        int subdivision;
        double Velocity;
        double Radius;

        public double[] VelocityField(double[] X, double t)
        {
            //angular Position
            double phi = Math.Atan2(X[1], X[0]);
            // Outward pointing normal Vector
            double vectorSize = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            double nx = X[0] / vectorSize;
            double ny = X[1] / vectorSize;

            double AngularOscillation = Math.Cos(phi * subdivision);
            double TemporalOscillation = Math.Cos(t * Math.PI) ;

            return new double[] { nx * AngularOscillation * Velocity * TemporalOscillation, ny * AngularOscillation * Velocity * TemporalOscillation };
        }

        /// <summary>
        /// Gives a distance in R-Direction to the exact Position of the interface
        /// </summary>
        /// <param name="x">Cartesian Coordinates of the Point</param>
        /// <param name="t">Time</param>
        /// <returns>Distance</returns>
        public double ContourError(double[]x, double t){
            double rPosition = x.L2Norm();
            return rPosition - rExact(x,t);
        }

        public double rExact(double[]X, double t) {
            double phi = Math.Atan2(X[1], X[0]);
            return Radius + Math.Cos(subdivision * phi) * Velocity / Math.PI * Math.Sin(Math.PI * t);
        }

        /// <summary>
        /// Calculates the Area enclosed by the flower contour at the time t
        /// </summary>
        /// This is an automatically generated and optimzied routine from maple:
        /// 
        /// r := R+cos(n*phi)*u*sin(Pi*t)/Pi;
        /// myArea := simplify(int((1/2)*r^2, phi = 0 .. 2*Pi));
        /// with(CodeGeneration)
        /// CSharp(myArea, optimize, declare = [t::float])
        /// 
        /// <param name="t">Time</param>
        /// <returns>Area</returns>
        public double EnclosedArea(double t){
            double t1 = Math.PI * t;
            double t2 = Math.Cos(t1);
            double t3 = t2 * t2;
            double t4 = Velocity * Velocity;
            double t5 = t4 * t3;
            double t6 = Math.PI * subdivision;
            double t7 = Math.Cos(t6);
            double t8 = t7 * t7;
            double t9 = t8 * t7;
            double t10 = Math.Sin(t6);
            double t19 = Math.Sin(t1);
            double t27 = t4 * subdivision;
            double t29 = Math.PI * Math.PI;
            double t31 = Radius * Radius;
            double t44 = -0.1e1 / subdivision / t29 * (-0.4e1 * t7 * t10 * Math.PI * Radius * t19 * Velocity - 0.2e1 * subdivision * t31 * t29 * Math.PI + t27 * Math.PI * t3 + t7 * t10 * t4 - 0.2e1 * t10 * t9 * t4 - t7 * t10 * t5 + 0.2e1 * t10 * t9 * t5 - Math.PI * t27) / 0.2e1;

            return t44;

        }
    }
}
