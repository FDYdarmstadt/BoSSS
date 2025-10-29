using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSFE_Solver;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using HFSISolver.SolidPhase;
using Newtonsoft.Json;


namespace HFSISolver {
    public class HFSI_Control : XNSFE_Control {
        [DataMember]
        public Solid Material = new AllOne();

        [DataMember]
        public int Degree { get; private set; }

        public HFSI_Control() {
            NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            UseImmersedBoundary = true;
            Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            Option_LevelSetEvolution2 = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            AdvancedDiscretizationOptions.ViscosityMode = BoSSS.Solution.XNSECommon.ViscosityMode.FullySymmetric;
            AdvancedDiscretizationOptions.PenaltySafety = 1;
        }

        public override Type GetSolverType() {
            return typeof(HFSI);
        }

        public HFSI_Control(int p) : this() {

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

            FieldOptions.Add(VariableNames.DisplacementZ, new FieldOpts() {
                //Degree = p + DisplacementDegOffset,
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

        }

        /// <summary>
        /// Exact solution for Displacement, for each species (either A or B or C).
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public IDictionary<string, Func<double[], double, double>[]> ExactSolutionDisplacement;
    }


}
