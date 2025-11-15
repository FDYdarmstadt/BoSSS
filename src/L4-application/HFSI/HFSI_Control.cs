using BoSSS.Application.XNSFE_Solver;
using BoSSS.Solution.Control;
using System;
using System.Runtime.Serialization;
using HFSISolver.SolidPhase;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;


namespace HFSISolver {
    public class HFSI_Control : XNSFE_Control {
        [DataMember]
        public Solid Material = new AllOne();

        //[DataMember]
        //public int Degree { get; private set; }

        public HFSI_Control() {
            NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            UseImmersedBoundary = true;
            Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            Option_LevelSetEvolution2 = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            AdvancedDiscretizationOptions.ViscosityMode = BoSSS.Solution.XNSECommon.ViscosityMode.FullySymmetric;
            AdvancedDiscretizationOptions.PenaltySafety = 1;
            base.StokesExtentionUseBCmap = StokesExtentionBoundaryOption.ZeroFlux;
            base.CutCellQuadratureType = BoSSS.Foundation.XDG.CutCellQuadratureMethod.Saye;
        }

        public override Type GetSolverType() {
            return typeof(HFSI);
        }

        public HFSI_Control(int p) : this() {

            SetDGdegree(p);

    
        }

        public override void SetDGdegree(int p) {
            base.SetDGdegree(p);

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
       
    }


}
