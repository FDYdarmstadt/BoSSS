using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Solution;
using System;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid.Classic;
using ApplicationWithIDT.OptiLevelSets;

namespace BUIDT {
    internal class BUIDTQueries {
        /// <summary>
        /// the L2error of the Concentration is computed against a referenceSol
        /// </summary>
        /// <param name="fieldName"></param>
        /// <param name="referenceSolution"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        public static Query L2Error(string fieldName, Func<double[], double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is BUIDTMain program)) {
                    throw new NotSupportedException("Only valid for applications of type IDTControl");
                }
                // XDG field
                DGField dgField = app.IOFields.Single(f => f.Identification == fieldName);
                return dgField.L2Error(referenceSolution.Vectorize(time));
            };
        }
        /// <summary>
        /// the L2error of the splineLevelSet iso contour as 1D function y ->f(y) is computed against a referenceSol
        /// </summary>
        /// <param name="referenceSolution"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
        public static Query L2Error1D(Func<double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
               
                if(!(app is BUIDTMain program)) {
                    throw new NotSupportedException("Only valid for applications of type IDTControl");
                }
                if(!(program.LevelSetOpti is SplineOptiLevelSet splineLS)) {
                    throw new NotSupportedException("Only valid for splineLevelSet");
                }
                var grid = Grid1D.LineGrid(splineLS.y);
                Basis basis1D = new Basis(grid, 6);
                SinglePhaseField exfield1D = new SinglePhaseField(basis1D, "1dfield");
                exfield1D.ProjectField(y => referenceSolution(y));
                return exfield1D.L2Error(y => splineLS.Spline.Interpolate(y));
            };
        }
    }
}
