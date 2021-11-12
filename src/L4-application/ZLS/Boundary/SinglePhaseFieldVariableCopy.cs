using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Tecplot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    public class SinglePhaseFieldVariableCopy : ParameterS, ILevelSetParameter {

        string[] parameterNames;
        string[] fieldNames;
        string species;
        static int IDs = 0;
        int ID = 0;
        public SinglePhaseFieldVariableCopy( string species, string[] sourceFields, string[] parameterNames){
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
                
                DGField field = ((XDGField)DomainVarFields[fieldName]).GetSpeciesShadowField(species);
                
                string paramName = parameterNames[d];
                SinglePhaseField parameter = (SinglePhaseField)ParameterVarFields[paramName];

                parameter.CopyFrom(field);
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

    public class XDGFieldVariableCopy : ParameterS, ILevelSetParameter
    {

        string[] parameterNames;
        string[] fieldNames;
        static int IDs = 0;
        int ID = 0;
        public XDGFieldVariableCopy(string[] sourceFields, string[] parameterNames)
        {
            int D = sourceFields.Length;
            this.fieldNames = sourceFields;
            this.parameterNames = parameterNames;
            ID = IDs;
            IDs++;
        }

        public override IList<string> ParameterNames => parameterNames;

        public override DelParameterFactory Factory => ParameterFactory;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            int D = fieldNames.Length;
            for (int d = 0; d < D; ++d)
            {
                string fieldName = fieldNames[d];

                XDGField field = (XDGField)DomainVarFields[fieldName];

                string paramName = parameterNames[d];
                XDGField parameter = (XDGField)ParameterVarFields[paramName];
                
                parameter.CopyFrom(field);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            int D = fieldNames.Length;
            var copy = new (string, DGField)[D];

            for (int d = 0; d < D; ++d)
            {
                string fieldName = fieldNames[d];
                XDGBasis basis = ((XDGField)DomainVarFields[fieldName]).Basis;
                string paramName = parameterNames[d];
                XDGField param = new XDGField(basis, paramName);
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
