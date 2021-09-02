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
    
    /// <summary>
    /// Derivative of <see cref="LevelSet"/>, mainly for the aesthetics 
    /// of having a pair of evolver and level-set.
    /// </summary>
    public class FourierLevelSet : LevelSet {
        public FourierLevSetBase Fourier_LevSet;

        public FourierLevelSet(FourierLevSetControl control, Basis b, string id) : base( b, id) {
            Fourier_LevSet = FourierLevelSetFactory.Build(control);
        }
    }

    /// <summary>
    /// Evolution of the Fourier level-set by moving the Fourier nodes with the flow velocity.
    /// Obviously, this also involves some interpolation.
    /// 
    /// The idea is to describe an interface explicitly by means of Fourier series.
    /// This explicit representation is then converted into an implicit one and projected 
    /// onto a DG field -- from this point on, it is handled as any interface representation
    /// for the XDG method.
    /// Obviously, this is only possible for simple topologies, e.g. height in dependence of x-coordinate (e.g. layers of fluids)
    /// or radius in dependence of angel (2D bubbles and droplets).
    /// 
    /// The motivation behind the Fourier level-set is to be a very high accurate reference; 
    /// Since the representation is infinitely differentiable, the curvature can be computed up to machine accuracy,
    /// ruling out any artificial oscillations from in-precise derivatives.
    /// </summary>
    /// <remarks>
    /// - Mainly a driver around <see cref="FourierLevSetTimestepper"/>
    /// - Developed and maintained by Martin Smuda mainly through 2016 to 2020
    /// - Also used in Master/Bachelor thesis of J. Triebwasser and O. Yotov.
    /// </remarks>
    public class FourierEvolver : ParameterS, ILevelSetParameter, ILevelSetEvolver {
        
        /// <summary>
        /// specialized timestepper (Runge-Kutta-based) for the evolution of the Fourier-LS
        /// </summary>
        FourierLevSetTimestepper Fourier_Timestepper;

        int curvatureDegree;

        string levelSetName;

        IList<string> parameters;

        public IList<string> VariableNames => null;

        public override IList<string> ParameterNames => parameters;

        // nothing to do
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => null;

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

        /// <summary>
        /// <see cref="ILevelSetEvolver.InternalFields"/>; here, empty;
        /// </summary>
        public IDictionary<string, DGField> InternalFields { 
            get { return null; }
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public FourierEvolver(string interfaceName, FourierLevelSet ls , FourierLevSetControl control, int curvatureDegree) {
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
            FourierLevelSet ls = (FourierLevelSet)levelSet.DGLevelSet;
            ls.Fourier_LevSet.ProjectToDGCurvature(Curvature, out filtLevSetGradient, levelSet.Tracker.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex));
        }

        public void MovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if(levelSet.Tracker.GridDat.SpatialDimension != 2)
                throw new NotSupportedException("Only supported in 2D.");

            SinglePhaseField[] meanVelocity = new SinglePhaseField[]
            {
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityX)],
                (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName,BoSSS.Solution.NSECommon.VariableNames.VelocityY)],
            };

            Fourier_Timestepper.updateFourierLevSet();
            Fourier_Timestepper.moveLevelSet(dt, meanVelocity, levelSet.Tracker.Regions.GetNearFieldMask(1));
            if (incremental)
                Fourier_Timestepper.updateFourierLevSet();
            FourierLevelSet ls = (FourierLevelSet)levelSet.DGLevelSet;
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
