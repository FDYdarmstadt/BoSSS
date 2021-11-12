using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase
{
    public class Laplace : ParameterS, ILevelSetParameter
    {

        string[] parameterNames;
        string[] fieldNames;
        string species;
        static int IDs = 0;
        int ID = 0;
        DGField[] temp;

        public Laplace(string species, string[] sourceFields, string[] parameterNames)
        {
            int D = sourceFields.Length;
            this.fieldNames = sourceFields;
            this.species = species;
            this.parameterNames = parameterNames;
            ID = IDs;
            temp = new DGField[D];
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

                DGField field = ((XDGField)DomainVarFields[fieldName]).GetSpeciesShadowField(species);

                Console.WriteLine("Bump");

                string paramName = parameterNames[d];
                DGField parameter = (DGField)ParameterVarFields[paramName];
                parameter.Clear();
                var subgrid = levelSet.Tracker.Regions.GetSpeciesMask(species);
                parameter.Laplacian(1.0, field);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        {
            int D = fieldNames.Length;
            var copy = new (string, DGField)[D];

            for (int d = 0; d < D; ++d)
            {
                string fieldName = fieldNames[d];
                int degree = DomainVarFields[fieldName].Basis.Degree;
                Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
                string paramName = parameterNames[d];
                SinglePhaseField param = new SinglePhaseField(basis, paramName);
                temp[d] = new SinglePhaseField(basis, paramName + "temp");
                copy[d] = (paramName, param);

            }
            return copy;
        }
    }
}
