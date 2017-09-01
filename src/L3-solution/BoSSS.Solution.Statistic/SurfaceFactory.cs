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
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Statistic.QuadRules {
    
    /// <summary>
    /// 
    /// </summary>
    public class SurfaceFactory {


        /// <summary>
        /// Class implements a factory for some stationary level sets: Sphere, Torus, Plane
        /// </summary>
        public SurfaceFactory() {

        }
        /// <summary>
        /// Factory for some stationary level sets
        /// </summary>
        /// <param name="surf">To be replaced by the level set class....</param>
        /// <param name="specification">Switch for the different types of surfaces: Sphere, Torus or Plane</param>
        /// <param name="b">Basis for the level set</param>
        /// <param name="context">The context object</param>
        /// <param name="deltax">Parameter for the signum approximation, should correspond to the maximal side length in case of a Cartesian grid </param>
        public static void ProduceSurface(out Surface surf, String specification, Basis b, GridData context, double deltax) {

            switch (specification) {
                /* By default: a unit sphere is created.*/
                case "Sphere":
                    surf = new Sphere(context, b, deltax, 1.0);
                    break;
                /* By default: a torus with major radius 1.0 and minor radius 0.25 is created.*/
                case "Torus":
                    surf = new Torus(context, b, deltax, 1.0, 0.25);
                    break;
                case "Plane":
                    surf = new Plane(context, b, deltax);
                    break;

                default: Console.WriteLine("This option is not implemented. Please select Sphere,Torus or Plane.");
                    surf = null;
                    break;
            }
        }

    }
}
