using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XDGTest {

    /// <summary>
    /// collection of tests
    /// </summary>
    [TestFixture]
    public static class VariousTests {


        /// <summary>
        /// Tests that multiple calls to <see cref="LevelSetTracker.UpdateTracker"/> do no harm.
        /// </summary>
        [Test]
        public static void MultipleTrackerUpdateCalls([Values(0, 1, 2)] int NearRegionWidht) {

            var Nodes = GenericBlas.Linspace(-2, 2, 20);
            GridCommons grd = Grid2D.Cartesian2DGrid(Nodes, Nodes);
            var gdat = grd.GridData;
            var phi = new LevelSet(new Basis(gdat, 2), "Phi");
            var trk = new LevelSetTracker(gdat, CutCellQuadratureMethod.Saye, NearRegionWidht, new[] { "A", "B" }, phi);

            Func<double[], double> UA = X => X[0] * X[1];
            Func<double[], double> UB = X => X[0].Pow2();

            XDGField u = new XDGField(new XDGBasis(trk, 2), "u");
            // check that no one changes the default value accidentally or intentionally.
            // If you really want to change this, talk to Florian.
            Assert.AreEqual(u.UpdateBehaviour, BehaveUnder_LevSetMoovement.PreserveMemory, "Default update behavior has changed");
            u.GetSpeciesShadowField("A").ProjectField(UA);
            u.GetSpeciesShadowField("B").ProjectField(UB);

            void CheckError() {
                var Amask = trk.Regions.GetSpeciesMask("A");
                var Bmask = trk.Regions.GetSpeciesMask("B");
                double uAerr = u.GetSpeciesShadowField("A").L2Error(UA.Vectorize(), new CellQuadratureScheme(true, Amask));
                double uBerr = u.GetSpeciesShadowField("B").L2Error(UA.Vectorize(), new CellQuadratureScheme(true, Bmask));

                Assert.LessOrEqual(uAerr, 1.0e-10, "error for A to high");
                Assert.LessOrEqual(uAerr, 1.0e-10, "error for B to high");
            }

            for(int i = 0; i < 3; i++) {

                // check before update
                CheckError();

                // after one update
                trk.UpdateTracker(0.0);
                CheckError();

                // after more updates
                trk.UpdateTracker(0.0);
                trk.UpdateTracker(0.0);
                CheckError();

                // Push the stacks
                trk.PushStacks();
                trk.IncreaseHistoryLength(i+1);
            }
        }


        /// <summary>
        /// Line integrals over the intersection of a flat interface with the cell faces
        /// (<see cref="XQuadSchemeHelper.Get_SurfaceElement_EdgeQuadScheme"/>, used e.g. by the edge terms of the Laplace-Beltrami surface tension)
        /// on a 3D grid with hanging faces, against the exact values, for each edge:
        /// on a hanging face, the rule comes from the fine cell and must also be scaled with it
        /// (before the fix, the scaling was taken from the coarse cell and the line integrals on hanging faces were too large by a factor of 2).
        /// </summary>
        [Test]
        public static void SurfaceElementEdgeRule_HangingFaces3D() {

            // 4 x 4 x 4 cells on [-1, 1]^3, the half x < 0 refined once: hanging faces on x = 0
            var Nodes = GenericBlas.Linspace(-1, 1, 5);
            GridCommons grd = Grid3D.Cartesian3DGrid(Nodes, Nodes, Nodes);
            var gdat0 = grd.GridData;
            int[] toRefine = Enumerable.Range(0, gdat0.iLogicalCells.NoOfLocalUpdatedCells).Where(j => gdat0.iGeomCells.GetCenter(j)[0] < 0).ToArray();
            var gdat = gdat0.Adapt(toRefine, new int[0][], out _).GridData;

            // interface y = y0 (not on a grid plane)
            const double y0 = 0.13;
            Func<double[], double> PhiFunc = X => X[1] - y0;
            var phi = new LevelSet(new Basis(gdat, 2), "Phi");
            phi.ProjectField(PhiFunc);
            var trk = new LevelSetTracker(gdat, CutCellQuadratureMethod.Saye, 1, new[] { "A", "B" }, phi);
            trk.UpdateTracker(0.0);

            int order = 4;
            var schH = trk.GetXDGSpaceMetrics(trk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
            var rule = schH.Get_SurfaceElement_EdgeQuadScheme(trk.GetSpeciesId("A"), trk.GetSpeciesId("B"), 0).Compile(gdat, order);

            // integrals of 1 and of x + z over interface \cap edge
            var E2C = gdat.iGeomEdges.CellIndices;
            var E2Ctrafo = gdat.iGeomEdges.Edge2CellTrafoIndex;
            var integrals = new Dictionary<int, double[]>();
            EdgeQuadrature.GetQuadrature(new int[] { 2 }, gdat, rule,
                delegate (int e0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    for(int i = 0; i < Length; i++) {
                        var X = gdat.TransformLocal2Global(QR.Nodes.GetVolumeNodeSet(gdat, E2Ctrafo[e0 + i, 0], false), E2C[e0 + i, 0]);
                        for(int k = 0; k < QR.NoOfNodes; k++) {
                            EvalResult[i, k, 0] = 1.0;
                            EvalResult[i, k, 1] = X[k, 0] + X[k, 2];
                        }
                    }
                },
                delegate (int e0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        integrals[e0 + i] = new[] { ResultsOfIntegration[i, 0], ResultsOfIntegration[i, 1] };
                }).Execute();

            // exact values: the interface intersects the faces normal to x and z in straight lines
            var Kedge = gdat.iGeomEdges.EdgeRefElements.Single();
            int NoOfHanging = 0;
            for(int e = 0; e < gdat.iGeomEdges.Count; e++) {
                var Vtx = gdat.TransformLocal2Global(Kedge.Vertices.GetVolumeNodeSet(gdat, E2Ctrafo[e, 0], false), E2C[e, 0]);
                var min = new double[] { double.MaxValue, double.MaxValue, double.MaxValue };
                var max = new double[] { double.MinValue, double.MinValue, double.MinValue };
                for(int iVtx = 0; iVtx < Vtx.GetLength(0); iVtx++) {
                    for(int d = 0; d < 3; d++) {
                        min[d] = Math.Min(min[d], Vtx[iVtx, d]);
                        max[d] = Math.Max(max[d], Vtx[iVtx, d]);
                    }
                }
                int dNormal = Enumerable.Range(0, 3).OrderBy(d => max[d] - min[d]).First();
                double length = 0, intXZ = 0;
                if(dNormal != 1 && min[1] < y0 && y0 < max[1]) {
                    int dLine = dNormal == 0 ? 2 : 0; // direction of the line
                    length = max[dLine] - min[dLine];
                    intXZ = min[dNormal] * length + 0.5 * (max[dLine].Pow2() - min[dLine].Pow2());
                }

                bool hanging = E2C[e, 1] >= 0 && !(gdat.iGeomEdges.IsEdgeConformalWithCell1(e) && gdat.iGeomEdges.IsEdgeConformalWithCell2(e));
                if(hanging && length > 0)
                    NoOfHanging++;

                double[] I = integrals.ContainsKey(e) ? integrals[e] : new double[2];
                Assert.AreEqual(length, I[0], 1.0e-10, $"length of interface \\cap edge {e} (hanging: {hanging})");
                Assert.AreEqual(intXZ, I[1], 1.0e-10, $"integral of x + z over interface \\cap edge {e} (hanging: {hanging})");
            }
            Assert.Greater(NoOfHanging, 0, "no hanging face cut by the interface");
        }


    }
}
