using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    public class SmoothSinglePhaseFieldVariableCopy : ParameterS, ILevelSetParameter {

        string[] parameterNames;
        string[] fieldNames;
        string species;
        static int IDs = 0;
        int ID = 0;
        public SmoothSinglePhaseFieldVariableCopy(string species, string[] sourceFields, string[] parameterNames) {
            int D = sourceFields.Length;
            this.fieldNames = sourceFields;
            this.species = species;
            this.parameterNames = parameterNames;
            ID = IDs;
            IDs++;
        }

        public override IList<string> ParameterNames => parameterNames;

        public override DelParameterFactory Factory => ParameterFactory;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = fieldNames.Length;
            for(int d = 0; d < D; ++d) {
                string fieldName = fieldNames[d];

                ConventionalDGField field = ((XDGField)DomainVarFields[fieldName]).GetSpeciesShadowField(species);

                string paramName = parameterNames[d];
                SinglePhaseField parameter = (SinglePhaseField)ParameterVarFields[paramName];


                CellMask cells = levelSet.Tracker.Regions.GetSpeciesMask(species);

                L2PatchRecovery smoother = new L2PatchRecovery(field.Basis, parameter.Basis, cells);
                smoother.Perform(parameter, field);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            int D = fieldNames.Length;
            var copy = new (string, DGField)[D];

            for(int d = 0; d < D; ++d) {
                string fieldName = fieldNames[d];
                int degree = DomainVarFields[fieldName].Basis.Degree;
                Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
                string paramName = parameterNames[d];
                SinglePhaseField param = new SinglePhaseField(basis, paramName);
                copy[d] = (paramName, param);
            }

            return copy;
        }
    }
}
