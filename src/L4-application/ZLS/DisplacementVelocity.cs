using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver
{
    class DisplacementVelocity : ParameterS, ILevelSetParameter
    {
        int D;
        string[] parameterNames;

        public DisplacementVelocity(int D)
        {
            this.D = D;
            parameterNames = VariableNames.Displacement0Vector(D);
        }

        public override IList<string> ParameterNames => parameterNames;

        public override DelParameterFactory Factory => ParameterFactory;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)
        {
            int D = parameterNames.Length;
            for (int d = 0; d < D; ++d)
            {
                string fieldName = BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d);

                XDGField field = ((XDGField)DomainVarFields[fieldName]);

                Console.WriteLine("Field: " + field.Identification);

                string paramName = parameterNames[d];
                XDGField parameter = (XDGField)ParameterVarFields[paramName];
                parameter.CopyFrom(field);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            var copy = new (string, DGField)[D];

            for (int d = 0; d < D; ++d)
            {
                string fieldName = BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d);
                int degree = DomainVarFields[fieldName].Basis.Degree;
                XDGBasis basis = (XDGBasis)DomainVarFields[fieldName].Basis;
                string paramName = parameterNames[d];
                XDGField param = new XDGField(basis, paramName);
                copy[d] = (paramName, param);
            }
            return copy;
        }
    }
}
