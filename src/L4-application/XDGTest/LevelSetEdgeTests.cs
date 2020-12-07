using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.XQuadFactoryHelper;

namespace BoSSS.Application.XDGTest
{

    /// <summary>
    /// Tests to verify correct quadrature rules when the Level set is or nearly is parallel to an edge and eventually lies on that edge
    /// </summary>
    [TestFixture]
    public static class LevelSetEdgeTests
    {
        /// <summary>
        /// LevelSet lies parallel to an inner edge and the distance goes to zero
        /// </summary>
        [Test]
        public static void LevelSetEdgeTest([Values(MomentFittingVariants.Saye, MomentFittingVariants.OneStepGaussAndStokes)] MomentFittingVariants QuadType)
        {
            // simple 2x1 grid
            double[] Xnodes = GenericBlas.Linspace(-1, 1, 3);
            double[] Ynodes = GenericBlas.Linspace(0, 1, 2);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
            grd.EdgeTagNames.Add(1, "BoundaryEdge");
            grd.DefineEdgeTags(delegate (double[] X) {
                if (Math.Abs(X[0] + 1) <= 1e-8 || Math.Abs(X[0] - 1) <= 1e-8 || Math.Abs(X[1]) <= 1e-8 || Math.Abs(X[1] - 1) <= 1e-8) return (byte)1;
                else return (byte)0;
            });
            var grddat = new GridData(grd);

            double eps = 0.01;
            Func<double[], double, double> PhiFunc = (X, t) => X[0] + eps * (1.0 - t);

            // Basis
            Basis basis = new Basis(grddat, 1);

            // Fields
            SinglePhaseField Phi = new SinglePhaseField(basis, "phi");
            LevelSet LevSet = new LevelSet(basis, "LevelSet");
            LevelSetTracker LsTrk = new LevelSetTracker(grddat, QuadType, 1, new string[] { "A", "B" }, LevSet);

            for (double dt = 0.0; dt <= 2.0; dt = dt + 0.1)
            {

                // exact solutions
                // ========================================
                #region exact solutions
                double area_A = 1 - eps * (1.0 - dt);
                double area_B = 1 + eps * (1.0 - dt);

                double surface_ex = 1.0;

                double edge_ex_boundary_A = 3 - 2 * eps * (1.0 - dt);
                double edge_ex_boundary_B = 3 + 2 * eps * (1.0 - dt);
                #endregion

                // Update LevelSet
                // ========================================
                Phi.ProjectField(X => PhiFunc(X, dt));
                LevSet.Clear();
                LevSet.Acc(1.0, Phi);
                LsTrk.UpdateTracker(dt);

                // validate correct LsTrk update ??
                // ========================================

                // validate correct species Volumes
                // ========================================
                SpeciesId spcAId = LsTrk.SpeciesIdS[0];
                SpeciesId spcBId = LsTrk.SpeciesIdS[1];
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), 1, 1).XQuadSchemeHelper;

                #region Volumes                
                double areaA = Volume(spcAId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(areaA - area_A) <= 1e-12);

                double areaB = Volume(spcBId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(areaB - area_B) <= 1e-12);
                #endregion

                // validate correct Interface quadrature
                // ========================================

                #region Surface                

                double surfaceA = LevelSetArea(spcAId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(surfaceA - surface_ex) <= 1e-12);

                double surfaceB = LevelSetArea(spcBId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(surfaceB - surface_ex) <= 1e-12);
                #endregion

                // validate correct Edge quadrature
                // ========================================
                #region Boundary Edges
                double edgeA = EdgeArea(spcAId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(edgeA - edge_ex_boundary_A) <= 1e-12);

                double edgeB = EdgeArea(spcBId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(edgeB - edge_ex_boundary_B) <= 1e-12);
                #endregion
            }
        }
        static double Volume(SpeciesId spcId, XQuadSchemeHelper SchemeHelper, GridData GridDat)
        {
            double area = 0.0;
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, GridDat,
                vqs.Compile(GridDat, 1),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            area = area.MPISum();
            return area;
        }

        static double LevelSetArea(SpeciesId spcId, XQuadSchemeHelper SchemeHelper, GridData GridDat)
        {
            double surface = 0.0;
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, GridDat,
                surfElemVol.Compile(GridDat, 1),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            surface = surface.MPISum();
            return surface;
        }

        static double EdgeArea(SpeciesId spcId, XQuadSchemeHelper SchemeHelper, GridData GridDat)
        {
            double edge = 0.0;
            var edgeElemVol = SchemeHelper.GetEdgeQuadScheme(spcId);
            EdgeQuadrature.GetQuadrature(new int[] { 1 }, GridDat,
                edgeElemVol.Compile(GridDat, 1),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    for (int k = 0; k < Length; k++)
                    {
                        double scale = GridDat.Edges.EdgeTags[i0 + k] == 1 ? 1.0 : 0.0;
                        EvalResult.Storage[k] = scale;
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        edge += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            edge = edge.MPISum();
            return edge;
        }

        /// <summary>
        /// The level Set lies at a small angle to an inner edge.
        /// It is then moved over the edge.
        /// </summary>
        [Test]
        public static void LevelSetCornerTest([Values(MomentFittingVariants.Saye, MomentFittingVariants.OneStepGaussAndStokes)] MomentFittingVariants QuadType,
                                        [Values(5.0, 3.0, 1.0, 0.573)] double raw_angle)
        {

            double[] Xnodes = GenericBlas.Linspace(-1, 1, 3);
            double[] Ynodes = GenericBlas.Linspace(-1, 1, 3);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
            grd.EdgeTagNames.Add(1, "BoundaryEdge");
            grd.DefineEdgeTags(delegate (double[] X) {
                if (Math.Abs(X[0] + 1) <= 1e-8 || Math.Abs(X[0] - 1) <= 1e-8 || Math.Abs(X[1] + 1) <= 1e-8 || Math.Abs(X[1] - 1) <= 1e-8) return (byte)1;
                else return (byte)0;
            });
            var grddat = new GridData(grd);

            double eps = 0.01;
            double angle = Math.PI * raw_angle / 180.0;
            Func<double[], double, double> PhiFunc = (X, t) => X[0] * Math.Tan(angle) - (X[1] - eps * (1.0 - t));

            // Basis
            Basis basis = new Basis(grddat, 1);


            // Fields
            SinglePhaseField Phi = new SinglePhaseField(basis, "phi");
            LevelSet LevSet = new LevelSet(basis, "LevelSet");
            LevelSetTracker LsTrk = new LevelSetTracker(grddat, QuadType, 1, new string[] { "A", "B" }, LevSet);

            for (double dt = 0.0; dt <= 2.0; dt = dt + 0.1)
            {

                // exact solutions
                // ========================================
                #region exact solutions
                double area_A = 2 - 2 * eps * (1.0 - dt);
                double area_B = 2 + 2 * eps * (1.0 - dt);

                double surface_ex = 2 / Math.Cos(angle);

                double edge_ex_boundary_A = 4 - 2 * eps * (1.0 - dt);
                double edge_ex_boundary_B = 4 + 2 * eps * (1.0 - dt);
                #endregion

                // Update LevelSet
                // ========================================
                Phi.ProjectField(X => PhiFunc(X, dt));
                LevSet.Clear();
                LevSet.Acc(1.0, Phi);
                LsTrk.UpdateTracker(dt);

                // validate correct LsTrk update ??
                // ========================================

                // validate correct species Volumes
                // ========================================
                SpeciesId spcAId = LsTrk.SpeciesIdS[0];
                SpeciesId spcBId = LsTrk.SpeciesIdS[1];
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), 1, 1).XQuadSchemeHelper;

                #region Volumes                
                double areaA = Volume(spcAId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(areaA - area_A) <= 1e-12);

                double areaB = Volume(spcBId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(areaB - area_B) <= 1e-12);
                #endregion

                // validate correct Interface quadrature
                // ========================================
                #region Surface                

                double surfaceA = LevelSetArea(spcAId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(surfaceA - surface_ex) <= 1e-12);

                double surfaceB = LevelSetArea(spcBId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(surfaceB - surface_ex) <= 1e-12);
                #endregion

                // validate correct Edge quadrature
                // ========================================

                #region Boundary Edges
                double edgeA = EdgeArea(spcAId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(edgeA - edge_ex_boundary_A) <= 1e-12);

                double edgeB = EdgeArea(spcBId, SchemeHelper, LsTrk.GridDat);
                Assert.IsTrue(Math.Abs(edgeB - edge_ex_boundary_B) <= 1e-12);
                #endregion
            }
        }
    }
}
