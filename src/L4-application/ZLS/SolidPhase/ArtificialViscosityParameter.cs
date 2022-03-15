using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class ArtificialViscosityParameter : ParameterS, ILevelSetParameter {

        string[] parameterNames;
        string[] fieldNames;
        string species;


        public ArtificialViscosityParameter(string species, int D) {
            this.fieldNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
            this.species = species;
            this.parameterNames = new string[] { ZwoLevelSetSolver.VariableNames.ArtificialViscosity};
        }

        public override IList<string> ParameterNames => parameterNames;

        public override DelParameterFactory Factory => ParameterFactory;

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = fieldNames.Length;
            XDGField[] velocity = new XDGField[D];
            for(int d = 0; d < D; ++d) {
                string fieldName = fieldNames[d];
                velocity[d] = ((XDGField)DomainVarFields[fieldName]);
            }
            string paramName = parameterNames[0];
            XDGField parameter = (XDGField)ParameterVarFields[paramName];
            Update(velocity, parameter);

        }

        void Update(XDGField[] velocity, XDGField parameter) {
            
            parameter.Clear();
            //parameter.ProjectAbs(1.0, null, velocity);
            parameter.Divergence(1.0, velocity);
            //XDGFieldDifferentiator.DivergenceByFlux(1.0, parameter, velocity, velocity[0].Basis.Tracker, "C");
            //parameter.AccConstant(0.001);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            
            var copy = new (string, DGField)[1];
            int d = 0;
            string fieldName = fieldNames[d];
            
            XDGBasis basis = ((XDGField)DomainVarFields[fieldName]).Basis;
            string paramName = parameterNames[d];
            XDGField param = new XDGField(basis, paramName);
            copy[d] = (paramName, param);
            
            int D = basis.GridDat.SpatialDimension;
            XDGField[] velocity = new XDGField[D];
            for(int i = 0; i < D; ++i) {
                fieldName = fieldNames[i];
                velocity[i] = ((XDGField)DomainVarFields[fieldName]);
            }
            Update(velocity, param);
            Console.WriteLine("Do me");
            return copy;
        }
    }
}