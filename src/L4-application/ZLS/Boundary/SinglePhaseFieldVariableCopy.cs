using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    public class SinglePhaseFieldVariableCopy : ILevelSetParameter {

        string[] parameterNames;
        string[] fieldNames;
        string species;

        public SinglePhaseFieldVariableCopy(string levelSetName, string species, string[] xdgFieldNames){
            int D = xdgFieldNames.Length;
            this.fieldNames = xdgFieldNames;
            this.species = species;
            parameterNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();    
        }
        
        public IList<string> ParameterNames => parameterNames;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = fieldNames.Length;
            for(int d = 0; d < D; ++d) {
                string fieldName = fieldNames[d];
                DGField field = ((XDGField)DomainVarFields[fieldName]).GetSpeciesShadowField(species);

                string paramName = parameterNames[d];
                SinglePhaseField parameter = (SinglePhaseField)ParameterVarFields[paramName];

                parameter.Clear();
                parameter.AccLaidBack(1.0, field);
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

    public class SinglePhaseFieldParameterCopy : ILevelSetParameter {

        string[] parameterNames;
        string[] fieldNames;
        string species;
        int degree;
        IGridData gridDat;

        public SinglePhaseFieldParameterCopy(string levelSetName, string species, string[] xdgParameterFieldNames, int degree, IGridData gridDat) {
            int D = xdgParameterFieldNames.Length;
            this.fieldNames = xdgParameterFieldNames;
            this.species = species;
            parameterNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            this.degree = degree;
            this.gridDat = gridDat;
        }

        public IList<string> ParameterNames => parameterNames;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = fieldNames.Length;
            for(int d = 0; d < D; ++d) {
                string fieldName = fieldNames[d];
                DGField field = ((XDGField)ParameterVarFields[fieldName]).GetSpeciesShadowField(species);

                string paramName = parameterNames[d];
                SinglePhaseField parameter = (SinglePhaseField)ParameterVarFields[paramName];

                parameter.Clear();
                parameter.AccLaidBack(1.0, field);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            int D = fieldNames.Length;
            var copy = new (string, DGField)[D];
            for(int d = 0; d < D; ++d) {;
                Basis basis = new Basis(gridDat, degree);
                string paramName = parameterNames[d];
                SinglePhaseField param = new SinglePhaseField(basis, paramName);
                copy[d] = (paramName, param);
            }

            return copy;
        }
    }
}
