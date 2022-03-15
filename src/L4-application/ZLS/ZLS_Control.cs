using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver {
    public class ZLS_Control : XNSE_Control {
        [DataMember]
        public Solid Material = new AllOne();

        [DataMember]
        public int Degree { get; private set; }

        public bool DisplacementExtension = false;

        public double ArtificialViscosity = 0.000;
        
        public double ExtensionArtificialViscosity = 0.000;

        public bool VelocityContinuity = true;

        public ZLS_Control() {
            NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            UseImmersedBoundary = true;
            Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            Option_LevelSetEvolution2 = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
        }

        public override Type GetSolverType() {
            return typeof(ZLS);
        }

        public ZLS_Control(int p) : this() {
            
            Degree = p;
            SetDGdegree(p);

            FieldOptions.Add(VariableNames.DisplacementX, new FieldOpts() {
                //Degree = p + DisplacementDegOffset,
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            FieldOptions.Add(VariableNames.DisplacementY, new FieldOpts() {
                //Degree = p + DisplacementDegOffset,
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

        }

        //public static int DisplacementDegOffset = 0;
    }


}
