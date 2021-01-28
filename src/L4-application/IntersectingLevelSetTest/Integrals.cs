using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IntersectingLevelSetTest {
    class Integrals {

        public void Evaluate(LevelSetTracker LsTrk, int quadOrder, SpeciesId id, SpeciesId id1) {
            //var schemes = new XQuadSchemeHelper(LsTrk, this.momentFittingVariant, LsTrk.SpeciesIdS.ToArray());
            var schemes = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), quadOrder, 1).XQuadSchemeHelper;

            var cutCells = LsTrk.Regions.GetCutCellSubGrid().VolumeMask;


            double edge = EvaluateEdges(LsTrk.GridDat, quadOrder, schemes, id);
            double volume = EvaluateVolume(LsTrk.GridDat, quadOrder, schemes, id);
            double surface = EvaluateSurface(LsTrk, quadOrder, schemes, id, id1);
        }

        double EvaluateEdges(GridData gridData, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            var edgeScheme = schemes.GetEdgeQuadScheme(id);
            var EdgeB = EdgeQuadrature(edgeScheme.Compile(gridData, quadOrder), gridData);
            double integral = EdgeB.Sum();
            return integral;
        }

        private MultidimensionalArray EdgeQuadrature(ICompositeQuadRule<QuadRule> surfRule, GridData gridData) {
            int E = gridData.iLogicalEdges.Count;
            var ret = MultidimensionalArray.Create(E);

            BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, gridData,
                surfRule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    ret.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 })
                        .Set(ResultsOfIntegration.ExtractSubArrayShallow(-1, 0));
                }).Execute();

            return ret;
        }

        double EvaluateVolume(GridData gridData, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id) {
            var volScheme = schemes.GetVolumeQuadScheme(id);
            var Vol = CellQuadrature(volScheme.Compile(gridData, quadOrder), gridData);
            double integral = Vol.Sum();
            return integral;
        }

        /// <summary>
        /// Computes the volume of the cut cells/surface of the level-set according to <paramref name="surfRule"/>
        /// </summary>
        private MultidimensionalArray CellQuadrature(ICompositeQuadRule<QuadRule> surfRule, GridData gridData) {
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { 1 }, gridData,
                surfRule,
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

        double EvaluateSurface(LevelSetTracker lsTrkr, int quadOrder, XQuadSchemeHelper schemes, SpeciesId id, SpeciesId id2) {
            double integral = 0;
            for (int i = 0; i < lsTrkr.NoOfLevelSets; ++i) {
                CellQuadratureScheme surfScheme = schemes.GetLevelSetquadScheme(i, id, id2, lsTrkr.Regions.GetCutCellMask4LevSet(i));
                var surf = CellQuadrature(surfScheme.Compile(lsTrkr.GridDat, quadOrder), lsTrkr.GridDat);
                integral += surf.Sum();
            }
            return integral;
        }
    }
}
