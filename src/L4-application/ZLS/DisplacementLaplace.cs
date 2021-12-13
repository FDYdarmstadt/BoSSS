using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace ZwoLevelSetSolver {
    class DisplacementLaplace : ParameterS, ILevelSetParameter {
        int D;
        string[] parameterNames;

        public DisplacementLaplace(int D) {
            this.D = D;
            parameterNames = VariableNames.DisplacementLaplaceVector(D);
        }

        public override IList<string> ParameterNames => parameterNames;

        public override DelParameterFactory Factory => ParameterFactory;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = parameterNames.Length;
            for(int d = 0; d < D; ++d) {
                string fieldName = VariableNames.DisplacementComponent(d);

                //var field = ((XDGField)DomainVarFields[fieldName]).GetSpeciesShadowField(levelSet.Tracker.GetSpeciesId("C"));
                var field = ((XDGField)DomainVarFields[fieldName]);

                Console.WriteLine("Field: " + field.Identification);

                string paramName = parameterNames[d];
                var parameter = ((XDGField)ParameterVarFields[paramName]);
                parameter.Clear();
                //XDGFieldDifferentiator.DerivativeByFlux(1, parameter, field, field.Basis.Tracker, 0,"C");
                XDGFieldDifferentiator.LaplacianByFlux(1, parameter, field, "C");
                parameter.Laplacian(1.0, field);
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var copy = new (string, DGField)[D];

            for(int d = 0; d < D; ++d) {
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

