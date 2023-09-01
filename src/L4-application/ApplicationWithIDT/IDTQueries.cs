
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.XDGField;



/// <summary>
/// This stuff doesn't work like intended as the types dont match.
/// </summary>
namespace ApplicationWithIDT {
    public static class IDTQueries {

       
        public static Query L2ErrorPerSpecies(string fieldName, string speciesName, Func<double[], double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is Application<IDTControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type IDTControl");
                }

                IDTControl control = program.Control;
                LevelSetTracker lsTrk = program.LsTrk;

                // Quadrature
                SpeciesId speciesId = lsTrk.GetSpeciesId(speciesName);
                int quadDegree = control.NonlinearQuadratureDegree;
                var SchemeHelper = lsTrk.GetXDGSpaceMetrics(lsTrk.SpeciesIdS.ToArray(), program.Control.NonlinearQuadratureDegree, HistoryIndex: 1).XQuadSchemeHelper;
                CellQuadratureScheme scheme = SchemeHelper.GetVolumeQuadScheme(speciesId);

                // XDG field
                DGField dgField = app.IOFields.Single(f => f.Identification == fieldName);
                SpeciesShadowField shadowField = ((XDGField)dgField).GetSpeciesShadowField(speciesId);

                return shadowField.L2Error(referenceSolution.Vectorize(time), scheme);
            };
        }
        public static Query L2Error(string fieldName, Func<double[], double, double> referenceSolution) {
            return delegate (IApplication<AppControl> app, double time) {
                if(!(app is Application<IDTControl> program)) {
                    throw new NotSupportedException("Only valid for applications of type IDTControl");
                }
                // XDG field
                DGField dgField = app.IOFields.Single(f => f.Identification == fieldName);
                return dgField.L2Error(referenceSolution.Vectorize(time));
            };
        }

    }
}
