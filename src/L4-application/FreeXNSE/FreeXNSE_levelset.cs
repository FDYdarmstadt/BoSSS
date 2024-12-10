using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.XNSECommon;
using ilPSP.Tracing;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FreeXNSE {
    /// <summary>
    /// Computation of the fluid interface velocity for material interfaces, in the free surface setting, that is only phase "A" is used
    /// </summary>
    public class FreeSurfaceVelocity : ILevelSetParameter {
        protected int D;

        protected IList<string> parameters;

        protected int degree;

        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// ctor.
        /// </summary>
        public FreeSurfaceVelocity(string levelSetName, int D, int degree) {
            this.D = D;
            this.degree = degree;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
        }

        /// <summary>
        /// averaging velocity at interface
        /// </summary>
        public virtual void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(new FuncTrace()) {
                int D = levelSet.Tracker.GridDat.SpatialDimension;

                //Mean Velocity
                XDGField[] EvoVelocity; // = new XDGField[]
                try {
                    EvoVelocity = D.ForLoop(
                        d => (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d)]
                        );
                } catch {
                    Console.Error.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                    EvoVelocity = D.ForLoop(
                        d => (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0_d(d)]
                        );
                }

                DGField[] meanVelocity = new ConventionalDGField[D];

                for(int d = 0; d < D; d++) {
                    meanVelocity[d] = ParameterVarFields[ParameterNames[d]];
                    meanVelocity[d].Clear();
                    meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(levelSet.Tracker.SpeciesNames.First()));                    
                }
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[D];
            for(int d = 0; d < D; ++d) {
                Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(basis, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }
}
