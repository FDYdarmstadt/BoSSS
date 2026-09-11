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

        /// <summary>
        /// case 0, 1, 2
        /// </summary>
        [Test]
        static public void AllUp_globalOnly(
#if DEBUG
            [Values(0, 1, 2)] int caseNo,
            [Values(2)] int dimension,
            [Values(2)] int degree,
            [Values(2, 4)] int gridResolution,
            [Values(0, 1)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
#else
            [Values(0, 1, 2)] int caseNo,
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
#endif
            ) {
            AllUp(caseNo, dimension, degree, gridResolution, AMRlevel, projectOnSameBasis, ProjectionStrategy.globalOnly);
        }

        /// <summary>
        /// case 0
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case0(
#if DEBUG
            [Values(2)] int dimension,
            [Values(2)] int degree,
            [Values(2, 4)] int gridResolution,
            [Values(true, false)] bool projectOnSameBasis
#else
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(2, 4, 8)] int gridResolution,
            //[Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
#endif
            ) {
            AllUp(0, dimension, degree, gridResolution, 0, projectOnSameBasis, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 1, degree 2;
        /// split off from the former monolithic `AllUp_patchwiseOnly_case1` (which took
        /// almost two hours in Release, running all degrees in a single job) so that the
        /// per-degree buckets can be run/scheduled as separate, parallel jobs.
        /// </summary>
#if DEBUG
        [Test]
        static public void AllUp_patchwiseOnly_case1_p2(
            [Values(2)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 1)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
            ) {
            AllUp(1, dimension, 2, gridResolution, AMRlevel, projectOnSameBasis, ProjectionStrategy.patchwiseOnly);
        }
#else
        /// <summary>
        /// case 1, degree 2, gridResolution 2 and 4;
        /// the skip condition in <see cref="AllUp"/> only excludes the expensive
        /// dimension==3/gridResolution==8 combination for degree &gt; 2, so at degree 2 that
        /// combination is NOT skipped and dominates the runtime of the old, unsplit
        /// AllUp_patchwiseOnly_case1_p2 (~89 min for 24 cases, vs. ~30s for the same 24 cases
        /// at degree 3). Isolating gridResolution==8 into its own job (mirroring the existing
        /// case2_..._p4_grd2n4 / _grd8 split below) keeps this bucket fast.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case1_p2_grd2n4(
            [Values(2, 3)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
            ) {
            AllUp(1, dimension, 2, gridResolution, AMRlevel, projectOnSameBasis, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 1, degree 2, gridResolution 8, dimension 2;
        /// see <see cref="AllUp_patchwiseOnly_case1_p2_grd2n4"/>. Unlike the dimension==3
        /// counterpart below, dimension==2 at gridResolution==8 is cheap (only a few hundred
        /// cells), so it is split off into its own fast job rather than sharing a job with the
        /// expensive dimension==3 cases.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case1_p2_grd8_dim2(
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
            ) {
            AllUp(1, 2, 2, 8, AMRlevel, projectOnSameBasis, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 1, degree 2, gridResolution 8, dimension 3;
        /// this is the actual bottleneck bucket. Unlike the other buckets here, the base grid
        /// alone already has ~4096 cells at gridResolution 8 in 3D (AMRlevel just adds a
        /// handful more), so all four combinations are similarly expensive and none of them
        /// piggy-backs cheaply on another — split down to one job per combination so they run
        /// fully in parallel instead of queued behind each other in a single job.
        /// See <see cref="AllUp_patchwiseOnly_case1_p2_grd2n4"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case1_p2_grd8_dim3_amr0_sameDegree() {
            AllUp(1, 3, 2, 8, 0, true, ProjectionStrategy.patchwiseOnly);
        }

        [Test]
        static public void AllUp_patchwiseOnly_case1_p2_grd8_dim3_amr0_higherDegree() {
            AllUp(1, 3, 2, 8, 0, false, ProjectionStrategy.patchwiseOnly);
        }

        [Test]
        static public void AllUp_patchwiseOnly_case1_p2_grd8_dim3_amr2_sameDegree() {
            AllUp(1, 3, 2, 8, 2, true, ProjectionStrategy.patchwiseOnly);
        }

        [Test]
        static public void AllUp_patchwiseOnly_case1_p2_grd8_dim3_amr2_higherDegree() {
            AllUp(1, 3, 2, 8, 2, false, ProjectionStrategy.patchwiseOnly);
        }
#endif

#if !DEBUG
        /// <summary>
        /// case 1, degree 3; see <see cref="AllUp_patchwiseOnly_case1_p2"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case1_p3(
            [Values(2, 3)] int dimension,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
            ) {
            AllUp(1, dimension, 3, gridResolution, AMRlevel, projectOnSameBasis, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 1, degree 4; see <see cref="AllUp_patchwiseOnly_case1_p2"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case1_p4(
            [Values(2, 3)] int dimension,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis
            ) {
            AllUp(1, dimension, 4, gridResolution, AMRlevel, projectOnSameBasis, ProjectionStrategy.patchwiseOnly);
        }
#endif

#if DEBUG
        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p2(
            [Values(2)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 1)] int AMRlevel
            ) {
            AllUp(2, dimension, 2, gridResolution, AMRlevel, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p2(
            [Values(2)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 1)] int AMRlevel
            ) {
            AllUp(2, dimension, 2, gridResolution, AMRlevel, false, ProjectionStrategy.patchwiseOnly);
        }
#else
        /// <summary>
        /// case 2, degree 2, gridResolution 2 and 4;
        /// the skip condition in <see cref="AllUp"/> does not apply to case 2 at all (it only
        /// triggers for caseNo==0 or degree&gt;2), so at degree 2 the expensive
        /// dimension==3/gridResolution==8 combination runs unconditionally and dominated the
        /// runtime of the old, unsplit case2_sameDegree_p2/case2_higherDegree_p2
        /// (~44 min / ~55 min for 12 cases each) -- the same issue found and fixed for
        /// <see cref="AllUp_patchwiseOnly_case1_p2_grd2n4"/>, split here the same way.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p2_grd2n4(
            [Values(2, 3)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 2, gridResolution, AMRlevel, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2, degree 2, gridResolution 8, dimension 2; see
        /// <see cref="AllUp_patchwiseOnly_case2_sameDegree_p2_grd2n4"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p2_grd8_dim2(
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, 2, 2, 8, AMRlevel, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2, degree 2, gridResolution 8, dimension 3 -- the actual bottleneck, split down
        /// to one job per AMR level as for <see cref="AllUp_patchwiseOnly_case1_p2_grd8_dim3_amr0_sameDegree"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p2_grd8_dim3_amr0() {
            AllUp(2, 3, 2, 8, 0, true, ProjectionStrategy.patchwiseOnly);
        }

        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p2_grd8_dim3_amr2() {
            AllUp(2, 3, 2, 8, 2, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2, degree 2, gridResolution 2 and 4; see
        /// <see cref="AllUp_patchwiseOnly_case2_sameDegree_p2_grd2n4"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p2_grd2n4(
            [Values(2, 3)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 2, gridResolution, AMRlevel, false, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2, degree 2, gridResolution 8, dimension 2; see
        /// <see cref="AllUp_patchwiseOnly_case2_sameDegree_p2_grd2n4"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p2_grd8_dim2(
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, 2, 2, 8, AMRlevel, false, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2, degree 2, gridResolution 8, dimension 3 -- the actual bottleneck; see
        /// <see cref="AllUp_patchwiseOnly_case2_sameDegree_p2_grd8_dim3_amr0"/>.
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p2_grd8_dim3_amr0() {
            AllUp(2, 3, 2, 8, 0, false, ProjectionStrategy.patchwiseOnly);
        }

        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p2_grd8_dim3_amr2() {
            AllUp(2, 3, 2, 8, 2, false, ProjectionStrategy.patchwiseOnly);
        }
#endif


        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p3(
            [Values(2, 3)] int dimension,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 3, gridResolution, AMRlevel, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p3(
            [Values(2, 3)] int dimension,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 3, gridResolution, AMRlevel, false, ProjectionStrategy.patchwiseOnly);
        }
#if! DEBUG
        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p4_grd2n4(
            [Values(2, 3)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 4, gridResolution, AMRlevel, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p4_grd2n4(
            [Values(2, 3)] int dimension,
            [Values(2, 4)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 4, gridResolution, AMRlevel, false, ProjectionStrategy.patchwiseOnly);
        }
        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_sameDegree_p4_grd8(
            [Values(2, 3)] int dimension,
            [Values(8)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 4, gridResolution, AMRlevel, true, ProjectionStrategy.patchwiseOnly);
        }

        /// <summary>
        /// case 2
        /// </summary>
        [Test]
        static public void AllUp_patchwiseOnly_case2_higherDegree_p4_grd8(
            [Values(2, 3)] int dimension,
            [Values(8)] int gridResolution,
            [Values(0, 2)] int AMRlevel
            ) {
            AllUp(2, dimension, 4, gridResolution, AMRlevel, false, ProjectionStrategy.patchwiseOnly);
        }
#endif



        static public void AllUp(
#if DEBUG
            [Values(0, 1, 2)] int caseNo,
            [Values(2)] int dimension,
            [Values(2)] int degree,
            [Values(2, 4)] int gridResolution,
            [Values(0, 1)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
#else
            [Values(0, 1, 2)] int caseNo,
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
#endif
            ) {

            CDGprojectionMain p = null;
            if (caseNo == 0 && dimension == 3
                || dimension == 3 && degree > 2 && gridResolution == 8
                || caseNo == 0 && AMRlevel > 0)
                return;

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new CDGprojectionMain();
                p.projectionCase = caseNo;
                p.dimension = dimension;
                p.degree = degree;
                p.gridResolution = gridResolution;
                p.AMRlevel = AMRlevel;
                p.projectOnSameBasis = projectOnSameBasis;
                p.projectStrategy = projectStrategy;
                return p;
            });

            Assert.IsTrue(p.passed);
        }

        /// <summary>
        /// LegendrePolynomial: case 3, 4, 5
        /// </summary>
        [Test]
        static public void AllUp_LegendrePolynomial(
#if DEBUG
            [Values(3, 4, 5)] int caseNo,
            [Values(2)] int dimension,
            [Values(2)] int degree,
            [Values(0, 1)] int AMRlevel,
            [Values(4)] int gridResolution,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
#else
            [Values(3, 4, 5)] int caseNo,
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(0, 2)] int AMRlevel,
            [Values(4, 8)] int gridResolution,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
#endif
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
                p.AMRlevel = AMRlevel;
                p.projectOnSameBasis = projectOnSameBasis;
                p.projectStrategy = projectStrategy;
                return p;
            });

            Assert.IsTrue(p.passed);
        }

        /// <summary>
        /// Cube: case 6
        /// </summary>
        [Test]
        static public void AllUp_Cube(
#if DEBUG
            [Values(2)] int dimension,
            [Values(2)] int degree,
            [Values(2, 4)] int gridResolution,
            [Values(0, 1)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy

#else
            [Values(2, 3)] int dimension,
            [Values(2, 3, 4)] int degree,
            [Values(2, 4, 8)] int gridResolution,
            [Values(0, 2)] int AMRlevel,
            [Values(true, false)] bool projectOnSameBasis,
            [Values(ProjectionStrategy.globalOnly, ProjectionStrategy.patchwiseOnly)] ProjectionStrategy projectStrategy
#endif
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
                p.AMRlevel = AMRlevel;
                p.projectOnSameBasis = projectOnSameBasis;
                p.projectStrategy = projectStrategy;
                return p;
            });

            Assert.IsTrue(p.passed);
        }



    }



}
