using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading.Tasks.Dataflow;

namespace IntersectingLevelSetTest {
    static class Integrals {

        /// <summary>
        /// Evaluate on domain definde by negative region of both level sets
        /// on grid [-1,1]^2
        /// </summary>
        /// <param name="alpha">level set</param>
        /// <param name="beta">level set</param>
        /// <param name="resolution">Cells per dimension</param>
        /// <param name="levelSetDegree"></param>
        /// <param name="quadOrder"></param>
        /// <returns></returns>
        public static (double intersection, double edge, double volume, double surface) Evaluate2D(Func<Vector, double> alpha, Func<Vector, double> beta, int resolution, int levelSetDegree, int quadOrder) {
            LevelSetTracker LsTrk = CreateTracker(alpha.Vectorize(), beta.Vectorize(), resolution, levelSetDegree, 2);
            return Evaluate(LsTrk, quadOrder);
        }

        /// <summary>
        /// Evaluate on domain definde by negative region of both level sets
        /// on grid [-1,1]^3
        /// </summary>
        /// <param name="alpha">level set</param>
        /// <param name="beta">level set</param>
        /// <param name="resolution">Cells per dimension</param>
        /// <param name="levelSetDegree"></param>
        /// <param name="quadOrder"></param>
        /// <returns></returns>
        public static (double intersection, double edge, double volume, double surface) Evaluate3D(Func<Vector, double> alpha, Func<Vector, double> beta, int resolution, int levelSetDegree, int quadOrder) {
            LevelSetTracker LsTrk = CreateTracker(alpha.Vectorize(), beta.Vectorize(), resolution, levelSetDegree, 3);
            return Evaluate(LsTrk, quadOrder);
        }

        static LevelSetTracker CreateTracker(ScalarFunction alpha, ScalarFunction beta, int resolution, int degree, int dim) {
            double[] tics = GenericBlas.Linspace(-0.5, 0.5, resolution + 1);

            GridData grid = null;
            if(dim == 2) {
                grid = Grid2D.Cartesian2DGrid(tics, tics).GridData;
            } else if (dim == 3){
                grid = Grid3D.Cartesian3DGrid(tics, tics, tics).GridData;
            } else {
                throw new NotSupportedException();
            }

            string[,] speciesTable = new string[2, 2];
            speciesTable[0, 0] = "B"; // Liquid
            speciesTable[0, 1] = "C"; // Solid
            speciesTable[1, 0] = "A"; // Gas
            speciesTable[1, 1] = "C"; // Solid

            Basis basis = new Basis(grid, degree);
            LevelSet Alpha = new LevelSet(basis, "alpha");
            LevelSet Beta = new LevelSet(basis, "beta");
            Alpha.ProjectField(alpha);
            Beta.ProjectField(beta);

            LevelSetTracker tracker = new LevelSetTracker(grid, XQuadFactoryHelper.MomentFittingVariants.Saye, 2,
                speciesTable, Alpha, Beta);
            return tracker;
        }

        static (double intersection, double edge, double volume, double surface) Evaluate(LevelSetTracker LsTrk, int quadOrder) {
            //var schemes = new XQuadSchemeHelper(LsTrk, this.momentFittingVariant, LsTrk.SpeciesIdS.ToArray());
            SpeciesId id = LsTrk.GetSpeciesId("B");
            var schemes = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), quadOrder, 1).XQuadSchemeHelper;

            double edge = EvaluateEdges(LsTrk.GridDat, quadOrder, schemes, id);
            double volume = EvaluateVolume(LsTrk.GridDat, quadOrder, schemes, id);
            double surface = EvaluateSurface(LsTrk, quadOrder, schemes, id);
            double intersection = EvaluateIntersection(LsTrk, quadOrder, schemes, id);
            double surfaceBoundary = EvaluateSurfaceBoundary(LsTrk, quadOrder, schemes, id);
            return (intersection, edge, volume, surface);
        }

        static double EvaluateEdges(GridData gridData, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            var edgeScheme = schemes.GetEdgeQuadScheme(id);
            var EdgeB = EdgeQuadrature(edgeScheme.Compile(gridData, quadOrder), gridData);
            double integral = EdgeB.Sum();
            return integral;
        }

        static MultidimensionalArray EdgeQuadrature(ICompositeQuadRule<QuadRule> rule, GridData gridData) {
            int E = gridData.iLogicalEdges.Count;
            var ret = MultidimensionalArray.Create(E);

            BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, gridData,
                rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    ret.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 })
                        .Set(ResultsOfIntegration.ExtractSubArrayShallow(-1, 0));
                }).Execute();

            return ret;
        }

        static double EvaluateVolume(GridData gridData, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            var volScheme = schemes.GetVolumeQuadScheme(id);
            var Vol = CellQuadrature(volScheme.Compile(gridData, quadOrder), gridData);
            double integral = Vol.Sum();
            return integral;
        }

        /// <summary>
        /// Computes the volume of the cut cells/surface of the level-set according to <paramref name="rule"/>
        /// </summary>
        static MultidimensionalArray CellQuadrature(ICompositeQuadRule<QuadRule> rule, GridData gridData) {
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { 1 }, gridData,
                rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    //Integrand.Evaluate(i0, Length, 0, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    var A = ret.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 });
                    var B = ResultsOfIntegration.ExtractSubArrayShallow(-1, 0);
                    A.Set(B);

                }).Execute();

            return ret;
        }

        static double EvaluateSurface(LevelSetTracker lsTrkr, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            double integral = 0;
            for (int i = 0; i < lsTrkr.NoOfLevelSets; ++i) {
                CellQuadratureScheme surfScheme = schemes.GetLevelSetquadScheme(i, id, lsTrkr.Regions.GetCutCellMask4LevSet(i));
                var surf = CellQuadrature(surfScheme.Compile(lsTrkr.GridDat, quadOrder), lsTrkr.GridDat);
                integral += surf.Sum();
            }
            return integral;
        }

        static double EvaluateIntersection(LevelSetTracker lsTrkr, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            double integral = 0;
            for (int i = 0; i < lsTrkr.NoOfLevelSets; ++i) {
                CellQuadratureScheme surfScheme = schemes.GetContactLineQuadScheme( id, i);
                var surf = CellQuadrature(surfScheme.Compile(lsTrkr.GridDat, quadOrder), lsTrkr.GridDat);
                integral += surf.Sum();
            }
            return integral;
        }

        static double EvaluateSurfaceBoundary(LevelSetTracker lsTrkr, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            double integral = 0;
            for (int i = 0; i < lsTrkr.NoOfLevelSets; ++i) {
                EdgeQuadratureScheme surfScheme = schemes.Get_SurfaceElement_EdgeQuadScheme(id, i);
                var surf = EdgeQuadrature(surfScheme.Compile(lsTrkr.GridDat, quadOrder), lsTrkr.GridDat);
                integral += surf.Sum();
            }
            return integral;
        }
    }
}
