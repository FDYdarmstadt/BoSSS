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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Runtime.Serialization;

namespace IntersectingLevelSetTest {

    internal class IntegrationTests {

        //Only works for resolution = 1, scaling required in eval method 
        [Test]
        public static void SquareTest([Values(1, 2, 3)]int resolution) {
            BoSSS.Solution.Application.InitMPI();
            Func<Vector, double> alpha = x => x.x;
            Func<Vector, double> beta = x => -x.y;
            (double intersection, double edge, double volume, double surface) = Integrals.Evaluate2D(alpha, beta, resolution, 1, 1);
            Assert.That(edge, Is.EqualTo(1).Within(1e-8));
            Assert.That(volume, Is.EqualTo(0.25).Within(1e-8));
            Assert.That(surface, Is.EqualTo(1).Within(1e-8));
            Assert.That(intersection, Is.EqualTo(2).Within(1e-8));
        }

        //Only works for resolution = 1, scaling required in eval method 
        [Test]
        public static void CubeTest([Values(1,2,3)]int resolution) {
            BoSSS.Solution.Application.InitMPI();
            Func<Vector, double> alpha = x => x.x;
            Func<Vector, double> beta = x => -x.y;
            (double intersection, double edge, double volume, double surface) = Integrals.Evaluate3D(alpha, beta, resolution, 1, 1);
            Assert.That(edge, Is.EqualTo(1.5).Within(1e-8));
            Assert.That(volume, Is.EqualTo(0.25).Within(1e-8));
            Assert.That(surface, Is.EqualTo(1).Within(1e-8));
            Assert.That(intersection, Is.EqualTo(1).Within(1e-8));
        }
    }
}