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
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Application.XNSE_Solver {
    abstract class XCommon<T> : XdgApplicationWithSolver<T> where T : XNSE_Control, new() {

        public abstract void SetOperatorEquations(int D, OperatorFactory opFactory);
        public abstract void SetOperatorParameter(int D, OperatorFactory opFactory);
        public virtual void SetLevelSetParameter(int D, OperatorFactory opFactory) {

            if (Control.AdvancedDiscretizationOptions.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                MaxSigma maxSigmaParameter = new MaxSigma(Control.PhysicalParameters, Control.AdvancedDiscretizationOptions, QuadOrder(), Control.dtFixed);
                opFactory.AddParameter(maxSigmaParameter);
                lsUpdater.AddLevelSetParameter("Phi", maxSigmaParameter);
            }

            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                    BeltramiGradient lsGradient = BeltramiGradient.CreateFrom(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter("Phi", lsGradient);
                    break;
                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                    BeltramiGradientAndCurvature lsGradientAndCurvature =
                        BeltramiGradientAndCurvature.CreateFrom(Control, "Phi", QuadOrder(), D);
                    opFactory.AddParameter(lsGradientAndCurvature);
                    lsUpdater.AddLevelSetParameter("Phi", lsGradientAndCurvature);
                    break;
                case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                    var fourrier = new FourierEvolver(
                        Control,
                        QuadOrder(),
                        Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Curvature].Degree);
                    lsUpdater.AddLevelSetParameter("Phi", fourrier);
                    lsUpdater.AddEvolver("Phi", fourrier);
                    opFactory.AddParameter(fourrier);
                    break;
                default:
                    break;
            }
        }

        public abstract void SetSpatialOperator(out XSpatialOperatorMk2 XOP, int D, OperatorFactory opFactory);

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {

            OperatorFactory opFactory = new OperatorFactory();

            SetOperatorEquations(D, opFactory);

            SetOperatorParameter(D, opFactory);

            SetLevelSetParameter(D, opFactory);
            
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
                    //Add this stuff later on
                    break;
                case LevelSetEvolution.FastMarching:
                    var fastMarcher = new FastMarcher(Control, QuadOrder(), levelSet.GridDat.SpatialDimension);
                    lsUpdater.AddEvolver("Phi", fastMarcher);
                    break;
                case LevelSetEvolution.None:
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

            var domainFields = CurrentState.Fields;
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Count);
            for (int iVar = 0; iVar < domainFields.Count; iVar++) {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            lsUpdater.InitializeParameters(DomainVarsDict, ParameterVarsDict);
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
