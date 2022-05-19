using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Level set velocity, i.e. parameters with name <see cref="BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(string, IList{string})"/>
    /// </summary>
    public class LevelSetVelocity : ParameterS, ILevelSetParameter {

        string[] m_ParameterNames;
        int D;

        public override IList<string> ParameterNames {
            get {
                return m_ParameterNames;
            }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public LevelSetVelocity(string LsName, int D) : base() {
            this.D = D;
            m_ParameterNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(LsName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
        }

        public override DelPartialParameterUpdate Update {
            get {
                return InternalParameterUpdate;
            }
        }

        void InternalParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                DGField[] velocities = new ConventionalDGField[D];
                for (int d = 0; d < D; d++) {
                    var velocity = ParameterVarFields[ParameterNames[d]];
                    velocity.Clear();

                    XDGField xVelocity = (XDGField)DomainVarFields["Var_" + ParameterNames[d]];
                    velocity.Acc(1.0, xVelocity.ProjectToSinglePhaseField(1));

                    velocities[d] = velocity;
                }
            } 
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            InternalParameterUpdate(time, DomainVarFields, ParameterVarFields);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocities = new (string, DGField)[D];
            for (int d = 0; d < D; ++d) {
                string paramName = ParameterNames[d];
                var bv = DomainVarFields["Var_" + paramName].Basis; // XDG Basis
                var b = new Basis(bv.GridDat, bv.Degree); // DG Basis
                DGField lsVelocity = new SinglePhaseField(b, paramName);
                velocities[d] = (paramName, lsVelocity);
            }
            return velocities;
        }
    }
}
