using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {
    abstract class XCommon<T> : XdgApplicationWithSolver<T> where T : XNSE_Control, new() {

        public abstract void SetOperatorEquations(int D, OperatorFactory opFactory);
        public abstract void SetOperatorParameter(int D, OperatorFactory opFactory);
        public abstract void SetSpatialOperator(out XSpatialOperatorMk2 XOP, int D, OperatorFactory opFactory);

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {

            OperatorFactory opFactory = new OperatorFactory();

            SetOperatorEquations(D, opFactory);

            SetOperatorParameter(D, opFactory);
           
            //Get Spatial Operator
            SetSpatialOperator(out XSpatialOperatorMk2 XOP, D, opFactory);

            return XOP;
        }

        public abstract void MultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel);

        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        protected override MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    var configsLevel = new List<MultigridOperator.ChangeOfBasisConfig>();

                    MultigridConfigLevel(configsLevel);

                    configs[iLevel] = configsLevel.ToArray();
                }
                return configs;
            }
        }

        protected abstract int QuadOrder();

        protected LevelSetUpdater lsUpdater;
        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;        

        protected override LevelSetTracker InstantiateTracker() {
            if (Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye
                && Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes) {
                throw new ArgumentException($"The XNSE solver is only verified for cut-cell quadrature rules " +
                    $"{XQuadFactoryHelper.MomentFittingVariants.Saye} and {XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes}; " +
                    $"you have set {Control.CutCellQuadratureType}, so you are notified that you reach into unknown territory; " +
                    $"If you do not know how to remove this exception, you should better return now!");
            }
            LevelSet levelSet = new LevelSet(new Basis(GridData, Control.FieldOptions["Phi"].Degree), "Phi");
            levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
            lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
            switch (Control.Option_LevelSetEvolution) {
                case LevelSetEvolution.Fourier:
                    var fourrier = new FourierEvolver(Control, QuadOrder());
                    lsUpdater.AddLevelSetParameter("Phi", fourrier);
                    lsUpdater.AddEvolver("Phi", fourrier);
                    break;
                case LevelSetEvolution.FastMarching:
                    var fastMarcher = new FastMarcher(Control, QuadOrder());
                    lsUpdater.AddLevelSetParameter("Phi", fastMarcher);
                    lsUpdater.AddEvolver("Phi", fastMarcher);
                    break;
                case LevelSetEvolution.None:
                    lsUpdater.AddLevelSetParameter("Phi", new CurvatureProvider(Control, QuadOrder()));
                    break;
                default:
                    throw new NotImplementedException();
            }

            return lsUpdater.Tracker;
        }

        public override double UpdateLevelset(DGField[] domainFields, double time, double dt, double UnderRelax, bool incremental) {
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Length);
            for (int iVar = 0; iVar < domainFields.Length; iVar++) {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            double residual = lsUpdater.UpdateLevelSets(DomainVarsDict, ParameterVarsDict, time, dt, UnderRelax, incremental);
            Console.WriteLine("Residual of level-set update: " + residual);
            return 0.0;
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            base.CreateEquationsAndSolvers(L);

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            lsUpdater.UpdateParameters(ParameterVarsDict, 0.0);

        }                    

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetFixedTimestep();
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"done with timestep {TimestepNo}");
            return dt;
        }
    }
}
