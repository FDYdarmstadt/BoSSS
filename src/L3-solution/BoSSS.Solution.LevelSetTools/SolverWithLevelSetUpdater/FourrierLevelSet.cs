using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public class FourrierLevelSet : LevelSet {
        public FourierLevSetBase Fourier_LevSet;

        public FourrierLevelSet(FourierLevSetControl control, Basis b, string id) : base( b, id) {
            Fourier_LevSet = FourierLevelSetFactory.Build(control);
        }
    }

    public class FourierEvolver : Parameter, ILevelSetParameter, ILevelSetEvolver {
        /// <summary>
        /// specialized timestepper (Runge-Kutta-based) for the evoultion of the Fourier-LS
        /// </summary>
        FourierLevSetTimestepper Fourier_Timestepper;

        int curvatureDegree;

        string levelSetName;

        IList<string> parameters;

        public IList<string> VariableNames => null;

        public override IList<string> ParameterNames => parameters;

        IList<string> ILevelSetParameter.ParameterNames {
            get {
                return new[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };
            }
        }

        IList<string> ILevelSetEvolver.ParameterNames {
            get {
                return BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(2));
            }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public FourierEvolver(string interfaceName, FourrierLevelSet ls , FourierLevSetControl control, int curvatureDegree) {
            parameters = new[] { BoSSS.Solution.NSECommon.VariableNames.Curvature };
            this.levelSetName = interfaceName;

            if (control == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            //create specialized fourier timestepper
            Fourier_Timestepper = FourierLevelSetFactory.Build_Timestepper(control, ls.Fourier_LevSet.GetFLSproperty(),
                ls.Fourier_LevSet.ComputeChangerate, ls.Fourier_LevSet.EvolveFourierLS);

            this.curvatureDegree = curvatureDegree;
        }

        public void LevelSetParameterUpdate(
            DualLevelSet levelSet,
            double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            Fourier_Timestepper.updateFourierLevSet();
            VectorField<SinglePhaseField> filtLevSetGradient;
            SinglePhaseField Curvature = (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature];
            FourrierLevelSet ls = (FourrierLevelSet)levelSet.DGLevelSet;
            ls.Fourier_LevSet.ProjectToDGCurvature(Curvature, out filtLevSetGradient, levelSet.Tracker.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex));
        }

        public void MovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            SinglePhaseField[] meanVelocity = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityX)],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName,BoSSS.Solution.NSECommon.VariableNames.VelocityY)],
            };
            
            Fourier_Timestepper.moveLevelSet(dt, meanVelocity, levelSet.Tracker.Regions.GetNearFieldMask(1));
            if (incremental)
                Fourier_Timestepper.updateFourierLevSet();
            FourrierLevelSet ls = (FourrierLevelSet)levelSet.DGLevelSet;
            ls.Fourier_LevSet.ProjectToDGLevelSet(levelSet.DGLevelSet, levelSet.Tracker);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            //Curvature
            (string ParameterName, DGField ParamField)[] fields = new (string, DGField)[1];
            IGridData gridData = DomainVarFields.First().Value.GridDat;
            Basis curvatureBasis = new Basis(gridData, curvatureDegree);
            string curvatureName = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            fields[0] = (curvatureName, new SinglePhaseField(curvatureBasis, curvatureName));
            return fields;
        }
    }
}
