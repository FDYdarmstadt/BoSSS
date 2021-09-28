using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver {
    public class ZLS_Control : XNSE_Control {

        public Solid Material = new HardSiliconeRubber();

        public int Degree { get; private set; }

        public ZLS_Control() : base() { }

        public override Type GetSolverType() {
            return typeof(ZLS);
        }

        public ZLS_Control(int p) {
            UseImmersedBoundary = true;
            Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            Option_LevelSetEvolution2 = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            Degree = p;
            SetDGdegree(p);

            FieldOptions.Add(VariableNames.DisplacementX, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            FieldOptions.Add(VariableNames.DisplacementY, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }
    }
}
