using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.ParameterizedLevelSet;
using ilPSP;
using ilPSP.Tracing;
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
    abstract public class ParameterizedLevelSet : LevelSet {
        

        public ParameterizedLevelSet(Basis b, string id) : base(b, id) {
        }

        protected abstract void ImplicitFromCurve(MultidimensionalArray X, MultidimensionalArray phi);

        public void Project() {
            this.ProjectField(1.0, ImplicitFromCurve);
        }

        abstract public double[] Parameters { get; set; }



        abstract public MyRealF_Base GetIntegrandForMinimi(CellQuadratureScheme _levSetQuadScheme, int _quadOrder);
    }


    /// <summary>
    /// Derivative of <see cref="LevelSet"/>, mainly for the aesthetics 
    /// of having a pair of evolver and level-set.
    /// </summary>
    public class ParameterizedLevelSetEllipse : ParameterizedLevelSet {
        internal double xSemiAxis;
        internal double ySemiAxis;
        internal double yCenter;

        public ParameterizedLevelSetEllipse(ParameterizedLevelSetControlEllipse control, Basis b, string id) : base(b, id) {
            this.xSemiAxis = control.xSemiAxis;
            this.ySemiAxis = control.ySemiAxis;
            this.yCenter = control.yCenter;

        }

        public override double[] Parameters {
            get {
                return new double[] { xSemiAxis, ySemiAxis, yCenter };
            }
            set {
                if (value.Length != 3)
                    throw new ArgumentException("expecting exactly 3 parameters for the ellipse");
                xSemiAxis = value[0];
                ySemiAxis = value[1];
                yCenter = value[2];
            }
        }

        public override MyRealF_Base GetIntegrandForMinimi(CellQuadratureScheme _levSetQuadScheme, int _quadOrder) {
            return new MyRealF_Ellipse(_levSetQuadScheme, this.GridDat, _quadOrder);
        }

        
        protected override void ImplicitFromCurve(MultidimensionalArray X, MultidimensionalArray phi) {
            int NoOfNodes = X.GetLength(0);
            for (int i = 0; i < NoOfNodes; i++) {
                double x = X[i, 0];
                double y = X[i, 1];
                phi[i] = y - yCenter + Math.Sqrt(ySemiAxis.Pow2() * (1 - x.Pow2() / xSemiAxis.Pow2()));
            }
        }
        
        
    }




    /// <summary>
    /// 
    /// </summary>
    /// <remarks>
    /// 
    /// </remarks>
    public class ParameterizedLevelSetEvolver : ILevelSetEvolver {

        /// <summary>
        /// specialized timestepper  for the evolution of the Parameterized-LS
        /// </summary>
        ParameterizedLevelSetTimeStepper Parameterized_TimeStepper;

        IList<string> parameters;

        int m_HMForder;

        int SpatialDimension;

        ParameterizedLevelSet m_ls;

        string levelSetName;
        public IList<string> ParameterNames => parameters;

        public IList<string> VariableNames => new string[] { };

        // nothing to do
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => null;

        /// <summary>
        /// <see cref="ILevelSetEvolver.InternalFields"/>; here, empty;
        /// </summary>
        public IDictionary<string, DGField> InternalFields {
            get { return null; }
        }


        public ParameterizedLevelSetEvolver(string interfaceName, ParameterizedLevelSet ls, ParameterizedLevelSetControl control, int hMForder, int D) {
            this.levelSetName = interfaceName;
            this.m_ls = ls;
            this.m_HMForder = hMForder;
            this.SpatialDimension = D;
           
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(interfaceName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));


            if (control == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of ParameterizedLevelSetControl!");

            //create specialized parameterized timestepper
            Parameterized_TimeStepper = ParameterizedLevelSetFactory.Build_Timestepper(control);
        }


        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;

                levelSet.DGLevelSet.Clear();
                levelSet.CGLevelSet.Clear();

                ParameterizedLevelSet ls = (ParameterizedLevelSet)levelSet.DGLevelSet;
                if (!object.ReferenceEquals(ls, m_ls)) {
                    throw new ApplicationException("level-set mismatch");
                }

                var LsTrk = levelSet.Tracker;
                int D = LsTrk.GridDat.SpatialDimension;
                SinglePhaseField[] meanVelocity = D.ForLoop(
                    d => (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))]
                    );

                var quadScheme = levelSet.Tracker.GetXDGSpaceMetrics(levelSet.Tracker.SpeciesIdS, this.m_HMForder).XQuadSchemeHelper.GetLevelSetquadScheme(levelSet.LevelSetIndex, levelSet.Tracker.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex));

                //update of elliptic parameters
                var myF = ls.GetIntegrandForMinimi(quadScheme, m_HMForder);
                var Param1 = Parameterized_TimeStepper.MoveLevelSet(dt, time, meanVelocity, ls.Parameters, levelSet.CGLevelSet.GridDat, myF);

                ls.Parameters = Param1;
                ((ParameterizedLevelSet)levelSet.DGLevelSet).Project();

            }
        }


            
    }
}
