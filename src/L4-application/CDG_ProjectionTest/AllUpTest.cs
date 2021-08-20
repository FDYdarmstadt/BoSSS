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

using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;
using BoSSS.Foundation.ConstrainedDGprojection;

namespace BoSSS.Application.CDG_ProjectionTest {

    /// <summary>
    /// Unit test for the continuous L2-projection
    /// </summary>
    [TestFixture]
    static public class AllUpTest {


        [Test]
        static public void AllUp(
            [Values(0, 1, 2)] int caseNo,
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(2, 4, 8)] int gridResolution,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
            ) {

            CDGprojectionMain p = null;
            if (caseNo == 0 && dimension == 3
                || dimension == 3 && degree > 2 && gridResolution == 8)
                return;

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new CDGprojectionMain();
                p.projectionCase = caseNo;
                p.dimension = dimension;
                p.degree = degree;
                p.gridResolution = gridResolution;
                p.projectOnSameBasis = projectOnSameBasis;
                p.projectStrategy = projectStrategy;
                return p;
            });

            Assert.IsTrue(p.passed);
        }

        [Test]
        static public void AllUp_LegendrePolynomial(
            [Values(3, 4, 5)] int caseNo,
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(4, 8)] int gridResolution,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
            ) {

            CDGprojectionMain p = null;
            if (dimension == 3 && degree > 2 && gridResolution == 8)
                return;

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new CDGprojectionMain();
                p.projectionCase = caseNo;
                p.dimension = dimension;
                p.degree = degree;
                p.gridResolution = gridResolution;
                p.projectOnSameBasis = projectOnSameBasis;
                p.projectStrategy = projectStrategy;
                return p;
            });

            Assert.IsTrue(p.passed);
        }

        [Test]
        static public void AllUp_Cube(
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(2, 4, 8)] int gridResolution,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
            ) {

            CDGprojectionMain p = null;
            if (dimension == 3 && degree > 2 && gridResolution == 8 ||
                dimension == 3 && degree == 4 && gridResolution == 4 && !projectOnSameBasis)
                return;

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new CDGprojectionMain();
                p.projectionCase = 6;
                p.dimension = dimension;
                p.degree = degree;
                p.gridResolution = gridResolution;
                p.projectOnSameBasis = projectOnSameBasis;
                p.projectStrategy = projectStrategy;
                return p;
            });

            Assert.IsTrue(p.passed);
        }



    }



}
