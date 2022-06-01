using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Application.LsTest {


    /// <summary>
    /// This solver merely projects some prescribed velocity field onto a xdg field
    /// and advects the level set.
    /// This process is however embedded in the "normal" XdgTimestepper frame, thus allows to investigate overall stability and behavior of e.g. coupled level set handling.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public partial class SolverWithLevelSetUpdaterTestCenter : SolverWithLevelSetUpdater<SolverWithLevelSetUpdaterTestControl> {
        protected override int NoOfLevelSets {
            get {
                return Control.NoOfLevelSets;
            }
        }

        protected override Array SpeciesTable {
            get {
                int N = Control.NoOfLevelSets;
                switch (N) {
                    case 1: {
                            return new[] { "A", "B" };
                            break;
                        }
                    case 2: {
                            var r = new string[2, 2];
                            r[0, 0] = "A"; // row 1st level set
                            r[0, 1] = "C"; // column 2nd level set
                            r[1, 0] = "B";
                            r[1, 1] = "C";
                            return r;
                            break;
                        }
                    default: {
                            throw new NotImplementedException();
                        }
                }
            }
        }

        public override int QuadOrder() {
            if (Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye
                && Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes) {
                throw new ArgumentException($"The SolverWithLevelSetUpdater solver is only verified for cut-cell quadrature rules " +
                    $"{XQuadFactoryHelper.MomentFittingVariants.Saye} and {XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes}; " +
                    $"you have set {Control.CutCellQuadratureType}, so you are notified that you reach into unknown territory; " +
                    $"If you do not know how to remove this exception, you should better return now!");
            }

            //QuadOrder
            int degU = LevelSetDegree();
            int quadOrder = degU * 3; // Use high quadorder, to correctly capture the advection term
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                //See remarks in XNSE
                quadOrder *= 2;
                quadOrder += 1;
            }

            return quadOrder;
        }

        /// <summary>
        /// We only solve the level set advection, use the level set degree (of DG level set, which is set over CG Level set, this is bonkers)
        /// </summary>
        protected int LevelSetDegree() {
            int p = -1;
            for (int i = 0; i < this.NoOfLevelSets; i++) {
                if (this.Control.FieldOptions.TryGetValue(VariableNames.LevelSetCGidx(i), out FieldOpts v)) {
                    p = Math.Max(v.Degree, p);
                } 
            }
            if(p == -1) throw new Exception("No Level Set Degree specified");
            return p;
        }

        private int LevelSetDegree(string LsName) {
            int p = -1;
            for (int i = 0; i < this.NoOfLevelSets; i++) {
                if (this.Control.FieldOptions.TryGetValue(LsName, out FieldOpts v)) {
                    p = Math.Max(v.Degree, p);
                }
            }
            if (p == -1) throw new Exception("No Level Set Degree specified");
            return p;
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            int i = 0;
            foreach (var field in CurrentState.Fields) {
                var config = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { field.Basis.Degree },
                    mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { i }
                };
                configsLevel.Add(config);
                i++;
            }
        }

        private IncompressibleBoundaryCondMap m_boundaryMap;

        /// <summary>
        /// Relation between 
        /// - edge tags (<see cref="Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>, passed to equation components via <see cref="BoSSS.Foundation.CommonParams.EdgeTag"/>)
        /// - boundary conditions specified in the control object (<see cref="AppControl.BoundaryValues"/>)
        /// </summary>
        protected IncompressibleBoundaryCondMap boundaryMap {
            get {
                if (m_boundaryMap == null)
                    m_boundaryMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);
                return m_boundaryMap;
            }
        }

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            return boundaryMap;
        }

        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            int D = GridData.SpatialDimension;
            string LsName = VariableNames.LevelSetCGidx(iLevSet);
            ILevelSetParameter levelSetVelocity = new LevelSetVelocity(LsName, D);
            return levelSetVelocity;
        }

        /// <summary>
        /// Dummy operator
        /// </summary>
        /// <param name="D"></param>
        /// <param name="levelSetUpdater"></param>
        /// <returns></returns>
        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater) {

            // Collect equation components
            Dictionary<string, List<IEquationComponent>> components = new Dictionary<string, List<IEquationComponent>>();
            for (int iLevSet = 0; iLevSet < Control.NoOfLevelSets; iLevSet++) {
                for (int d = 0; d < D; d++) {
                    var variable = "Var_" + VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(iLevSet), VariableNames.VelocityVector(D)[d]);
                    if (Control.InitialValues_Evaluators_TimeDep.TryGetValue(variable, out var func)) {
                        if (!components.ContainsKey(variable)) components[variable] = new List<IEquationComponent>();                        
                        components[variable].Add(new MultiPhasePrescribedVariable(variable, func));                        
                    }
                }
            }

            // Get Spatial Operator
            XSpatialOperatorMk2 XOP = new XSpatialOperatorMk2(components.Keys.Select(k => k).ToList(), components.Keys.Select(k => k).ToList(), (x,y,z) => this.QuadOrder(), LsTrk.SpeciesNames); // these are dummy fields
            foreach(var kvp in components) {
                foreach(var comp in kvp.Value) {
                    XOP.EquationComponents[kvp.Key].Add(comp);
                }
            }
            XOP.AgglomerationThreshold = 0.1;
            AgglomerationAlgorithm.RecoverFromAgglomerationFail = true;
            if (AgglomerationAlgorithm.RecoverFromAgglomerationFail)
                Console.WriteLine("Careful activated experimental agglomeration fail recovery - i.e. do not agglomerate when no target is found");
            XOP.TemporalOperator = new ConstantXTemporalOperator(XOP, 0.0);

            // === level set related parameters === //
            for (int iLevSet = 0; iLevSet < Control.NoOfLevelSets; iLevSet++) {
                string LsName = VariableNames.LevelSetCGidx(iLevSet);
                Normals normalsParameter = new Normals(LsName, D, LevelSetDegree(LsName));
                LsUpdater.AddLevelSetParameter(LsName, normalsParameter);

                GradientAndCurvature lsBGradient = new GradientAndCurvature(LsName, LevelSetDegree(LsName), LevelSetDegree(LsName), QuadOrder(), D);
                LsUpdater.AddLevelSetParameter(LsName, lsBGradient);
            }
            XOP.IsLinear = true;

            //final settings
            XOP.Commit();
            return XOP;
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            if((int)this.Control.TimeSteppingScheme >= 100)
                dt = Math.Min(Control.dtFixed, Control.Endtime - phystime);
            else
                dt = this.GetTimestep();

            Console.WriteLine("Starting timestep {0}, time {2}, dt {1}", TimestepNo, dt, phystime);
            bool success = Timestepping.Solve(phystime, dt);
            return dt;
        }
    }

    /// <summary>
    /// Minimal volume source term
    /// </summary>
    public class PrescribedVariable : IVolumeForm, ISupportsJacobianComponent {
        string[] variableName;
        Func<double[], double, double> m_func;

        public PrescribedVariable(string variable, Func<double[], double, double> func) {
            variableName = new string[]
            {
                variable
            };
            m_func = func;
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V |TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering => variableName;

        public IList<string> ParameterOrdering => new string[0];

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return (U[0] - m_func(cpv.Xglobal, cpv.time)) * V ;
        }
    }

    public class MultiPhasePrescribedVariable : PrescribedVariable, ISpeciesFilter {

        public MultiPhasePrescribedVariable(string parameter, Func<double[], double, double> func) : base(parameter, func) {
        }

        public string ValidSpecies => null; // Valid in the whole domain
    }
}
