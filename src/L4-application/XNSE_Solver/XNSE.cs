using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Control;
using System.Runtime.CompilerServices;
using NUnit.Framework.Constraints;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Foundation.IO;
using BoSSS.Solution.AdvancedSolvers;

namespace BoSSS.Application.XNSE_Solver
{
    class XNSE : XdgApplicationWithSolver<XNSE_Control>
    {
        LevelSetUpdater lsUpdater;

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        protected override MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig
        {
            get
            {
                int pVel;
                if(this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v))
                {
                    pVel = v.Degree;
                }
                else if(this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1))
                {
                    pVel = v1.Degree;
                }
                else
                {
                    throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
                }
                int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
                int D = this.GridData.SpatialDimension;

                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++)
                {
                    var configsLevel = new List<MultigridOperator.ChangeOfBasisConfig>();
                    if (this.Control.UseSchurBlockPrec)
                    {
                        // using a Schur complement for velocity & pressure
                        var confMomConti = new MultigridOperator.ChangeOfBasisConfig();
                        for (int d = 0; d < D; d++)
                        {
                            d.AddToArray(ref confMomConti.VarIndex);
                            //Math.Max(1, pVel - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                            pVel.AddToArray(ref confMomConti.DegreeS);
                        }
                        D.AddToArray(ref confMomConti.VarIndex);
                        //Math.Max(0, pPrs - iLevel).AddToArray(ref confMomConti.DegreeS); // global p-multi-grid
                        pPrs.AddToArray(ref confMomConti.DegreeS);

                        confMomConti.mode = MultigridOperator.Mode.SchurComplement;

                        configsLevel.Add(confMomConti);
                    }
                    else
                    {
                        // configurations for velocity
                        for (int d = 0; d < D; d++)
                        {
                            var configVel_d = new MultigridOperator.ChangeOfBasisConfig()
                            {
                                DegreeS = new int[] { pVel },
                                //DegreeS = new int[] { Math.Max(1, pVel - iLevel) },
                                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                                VarIndex = new int[] { d }
                            };
                            configsLevel.Add(configVel_d);
                        }
                        // configuration for pressure
                        var configPres = new MultigridOperator.ChangeOfBasisConfig()
                        {
                            DegreeS = new int[] { pPrs },
                            //DegreeS = new int[] { Math.Max(0, pPrs - iLevel) },
                            mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                            VarIndex = new int[] { D }
                        };
                        configsLevel.Add(configPres);
                    }
                    configs[iLevel] = configsLevel.ToArray();
                }
                return configs;
            }
        }

        int QuadOrder()
        {
            //QuadOrder
            int degU;
            if (Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts field))
            {
                degU = field.Degree;
            }
            else if (Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts field1))
            {
                degU = field1.Degree;
            }
            else
            {
                throw new Exception("Velocity not found!");
            }
            int quadOrder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) 
            {
                quadOrder *= 2;
                quadOrder += 1;
            }
            return quadOrder;
        }

        protected override LevelSetTracker InstantiateTracker() 
        {
            if (Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye
                && Control.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)
            {
                throw new ArgumentException($"The XNSE solver is only verified for cut-cell quadrature rules " +
                    $"{XQuadFactoryHelper.MomentFittingVariants.Saye} and {XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes}; " +
                    $"you have set {Control.CutCellQuadratureType}, so you are notified that you reach into unknown territory; " +
                    $"If you do not know how to remove this exception, you should better return now!");
            }
            LevelSet levelSet = new LevelSet(new Basis(GridData, Control.FieldOptions["Phi"].Degree), "Phi");
            levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
            lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
            switch (Control.Option_LevelSetEvolution)
            {
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

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        {
            base.CreateEquationsAndSolvers(L);

            var domainFields = CurrentState.Fields;
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Count);
            for (int iVar = 0; iVar < domainFields.Count; iVar++)
            {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++)
            {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            lsUpdater.InitializeParameters(DomainVarsDict, ParameterVarsDict);
            lsUpdater.UpdateParameters(ParameterVarsDict, 0.0);
            
        }

        public override double UpdateLevelset(DGField[] domainFields, double time, double dt, double UnderRelax, bool incremental)
        {
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Length);
            for (int iVar = 0; iVar < domainFields.Length; iVar++)
            {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++)
            {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            double residual = lsUpdater.UpdateLevelSets(DomainVarsDict, ParameterVarsDict, time, dt, UnderRelax, incremental);
            Console.WriteLine(residual);
            return 0.0;
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {

            int quadOrder = QuadOrder();
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            IncompressibleMultiphaseBoundaryCondMap boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            //Build Equations
            OperatorFactory opFactory = new OperatorFactory();
            for (int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NavierStokes("A", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A));
                opFactory.AddEquation(new NavierStokes("B", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, LsTrk, config));
                opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, LsTrk, config));
            }
            opFactory.AddParameter(new Velocity0(D));
            Velocity0Mean v0Mean = new Velocity0Mean(D, LsTrk, quadOrder);
            opFactory.AddParameter(v0Mean);
            lsUpdater.AddLevelSetParameter("Phi", v0Mean);

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);
            lsUpdater.AddLevelSetParameter("Phi", normalsParameter);
            
            if (config.isContinuity)
            {
                opFactory.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
                opFactory.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
                opFactory.AddEquation(new InterfaceContinuity(config, D, LsTrk));
            }

            if (Control.AdvancedDiscretizationOptions.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine)
            {
                MaxSigma maxSigmaParameter = new MaxSigma(Control.PhysicalParameters, Control.AdvancedDiscretizationOptions, QuadOrder(), Control.dtFixed);
                opFactory.AddParameter(maxSigmaParameter);
                lsUpdater.AddLevelSetParameter("Phi", maxSigmaParameter);
            }
            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode)
            {
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
                        BeltramiGradientAndCurvature.CreateFrom(Control, "Phi", quadOrder, D);
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
            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(quadOrder);
            
            //final settings
            XOP.FreeMeanValue[VariableNames.Pressure] = !boundaryMap.DirichletPressureBoundary;
            XOP.LinearizationHint = LinearizationHint.AdHoc;
            XOP.Commit();

            // test the ordering: Should not make a difference anyways!
            Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityX) < XOP.DomainVar.IndexOf(VariableNames.VelocityY));
            Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityY) < XOP.DomainVar.IndexOf(VariableNames.Pressure));
            Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationX) < XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationY));
            Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationY) < XOP.CodomainVar.IndexOf(EquationNames.ContinuityEquation));
            if(D > 2) {
                Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityY) < XOP.DomainVar.IndexOf(VariableNames.VelocityZ));
                Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityZ) < XOP.DomainVar.IndexOf(VariableNames.Pressure));
                Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationY) < XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationZ));
                Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationZ) < XOP.CodomainVar.IndexOf(EquationNames.ContinuityEquation));
            }
            return XOP;
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt)
        {
            //Update Calls
            dt = GetFixedTimestep();
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            return dt;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1)
        {
            Tecplot.PlotFields(this.m_RegisteredFields, "XNSE_Solver" + timestepNo, physTime, superSampling);
            if(Timestepping?.Parameters != null)
            {
                Tecplot.PlotFields(Timestepping.Parameters, "XNSE_Solver_Params" + timestepNo, physTime, superSampling);
            }
        }
    }
}
