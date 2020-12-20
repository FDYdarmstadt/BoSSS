using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;

namespace BoSSS.Application.XNSE_Solver {
    public class XHeat : SolverWithLevelSetUpdater<XNSE_Control> 
    {
        ThermalMultiphaseBoundaryCondMap boundaryMap;

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel)
        {

            int D = this.GridData.SpatialDimension;

            int pTemp = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Temperature].Degree;
            // configuration for Temperature
            var confTemp = new MultigridOperator.ChangeOfBasisConfig()
            {
                DegreeS = new int[] { pTemp }, //Math.Max(1, pTemp - iLevel) },
                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Temperature) }
            };
            configsLevel.Add(confTemp);

            // configuration for auxiliary heat flux
            if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
            {
                int pFlux;
                if (this.Control.FieldOptions.TryGetValue("HeatFlux*", out FieldOpts f))
                {
                    pFlux = f.Degree;
                }
                else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.HeatFluxX, out FieldOpts f1))
                {
                    pFlux = f1.Degree;
                }
                else
                {
                    throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of HeatFlux not found");
                }
                for (int d = 0; d < D; d++)
                {
                    var confHeatFlux = new MultigridOperator.ChangeOfBasisConfig()
                    {
                        DegreeS = new int[] { pFlux }, // Math.Max(1, pFlux - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.HeatFluxVectorComponent(d)) }
                    };
                    configsLevel.Add(confHeatFlux);
                }
            }
        }

        int QuadOrder() {
            //QuadOrder
            int degT;
            if (Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.Temperature, out FieldOpts field)) {
                degT = field.Degree;
            } else {
                throw new Exception("Temperature not found!");
            }
            int quadOrder = degT * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if (this.Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                quadOrder *= 2;
                quadOrder += 1;
            }
            return quadOrder;
        }

        int VelocityDegree()
        {
            int pVel;
            if (this.Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts v))
            {
                pVel = v.Degree;
            }
            else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts v1))
            {
                pVel = v1.Degree;
            }
            else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.Temperature, out FieldOpts t1))
            {
                Console.WriteLine("Degree of Velocity not found, using Temperature Degree");
                pVel = t1.Degree;
            }
            else
            {
                throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of Velocity not found");
            }
            return pVel;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1) {
            Tecplot.PlotFields(this.m_RegisteredFields, "XHEAT_Solver" + timestepNo, physTime, superSampling);
            if (Timestepping?.Parameters != null) {
                Tecplot.PlotFields(Timestepping.Parameters, "XHEAT_Solver_Params" + timestepNo, physTime, superSampling);
            }
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater)
        {
            OperatorFactory opFactory = new OperatorFactory();
            boundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());
            DefineSystem(D, opFactory, levelSetUpdater);

            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(QuadOrder());

            //final settings
            XOP.AgglomerationThreshold = this.Control.AgglomerationThreshold;
            XOP.Commit();

            return XOP;
        }

        public void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater)
        {
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);

            int quadOrder = QuadOrder();
            // add Heat equation components
            // ============================
            opFactory.AddEquation(new Heat("A", lsUpdater.Tracker, D, boundaryMap, config));
            opFactory.AddEquation(new Heat("B", lsUpdater.Tracker, D, boundaryMap, config));

            if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
            {
                for (int d = 0; d < D; ++d)
                {
                    opFactory.AddEquation(new HeatFlux("A", d, lsUpdater.Tracker, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFlux("B", d, lsUpdater.Tracker, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, boundaryMap, lsUpdater.Tracker, config));
                }
            }
            opFactory.AddEquation(new HeatInterface("A", "B", D, boundaryMap, lsUpdater.Tracker, config));
            opFactory.AddCoefficient(new EvapMicroRegion());

            // Add Velocity parameters as prescribed variables
            for (int d = 0; d < D; d++)
                opFactory.AddParameter(Velocity0Prescribed.CreateFrom(lsUpdater.Tracker, d, D, Control));

            Velocity0MeanPrescribed v0Mean = new Velocity0MeanPrescribed(D, lsUpdater.Tracker, quadOrder);
            opFactory.AddParameter(v0Mean);
            lsUpdater.AddLevelSetParameter("Phi", v0Mean);

            Normals normalsParameter = new Normals(D, ((LevelSet)lsUpdater.Tracker.LevelSets[0]).Basis.Degree);
            opFactory.AddParameter(normalsParameter);
            lsUpdater.AddLevelSetParameter("Phi", normalsParameter);
            // even though we have no surface tension here, we still need this
            switch (Control.AdvancedDiscretizationOptions.SST_isotropicMode) {
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine:
                    MaxSigma maxSigmaParameter = new MaxSigma(Control.PhysicalParameters, Control.AdvancedDiscretizationOptions, QuadOrder(), Control.dtFixed);
                    opFactory.AddParameter(maxSigmaParameter);
                    lsUpdater.AddLevelSetParameter("Phi", maxSigmaParameter);
                    BeltramiGradient lsBGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter("Phi", lsBGradient);
                    break;
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux:
                case SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local:
                    BeltramiGradient lsGradient = FromControl.BeltramiGradient(Control, "Phi", D);
                    lsUpdater.AddLevelSetParameter("Phi", lsGradient);
                    break;
                case SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint:
                case SurfaceStressTensor_IsotropicMode.Curvature_Projected:
                case SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean:
                    BeltramiGradientAndCurvature lsGradientAndCurvature =
                        FromControl.BeltramiGradientAndCurvature(Control, "Phi", quadOrder, D);
                    opFactory.AddParameter(lsGradientAndCurvature);
                    lsUpdater.AddLevelSetParameter("Phi", lsGradientAndCurvature);
                    break;
                case SurfaceStressTensor_IsotropicMode.Curvature_Fourier:
                    var fourrier = new FourierEvolver(
                        "Phi",
                        Control.FourierLevSetControl,
                        Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Curvature].Degree);
                    lsUpdater.AddLevelSetParameter("Phi", fourrier);
                    lsUpdater.AddEvolver("Phi", fourrier);
                    opFactory.AddParameter(fourrier);
                    break;
                default:
                    break;
            }
        }

        protected override LevelSetUpdater InstantiateLevelSetUpdater()
        {
            
            int levelSetDegree = Control.FieldOptions["Phi"].Degree;
            LevelSet levelSet = new LevelSet(new Basis(GridData, levelSetDegree), "Phi");
            levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
            LevelSetUpdater lsUpdater;
            switch (Control.Option_LevelSetEvolution)
            {
                case LevelSetEvolution.Fourier:
                    if (Control.EnforceLevelSetConservation) {
                        throw new NotSupportedException("mass conservation correction currently not supported");
                    }
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet, Control.LSContiProjectionMethod);
                    lsUpdater.AddLevelSetParameter("Phi", new LevelSetVelocity("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters));
                    break;
                case LevelSetEvolution.FastMarching:
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
                    var fastMarcher = new FastMarchingEvolver("Phi", QuadOrder(), levelSet.GridDat.SpatialDimension);
                    lsUpdater.AddEvolver("Phi", fastMarcher);
                    lsUpdater.AddLevelSetParameter("Phi", new LevelSetVelocity("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters));
                    break;
                case LevelSetEvolution.SplineLS:
                    SplineLevelSet SplineLevelSet = new SplineLevelSet(Control.Phi0Initial, levelSet.Basis, "Phi", (int)Math.Sqrt(levelSet.DOFLocal));
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, SplineLevelSet);
                    var SplineEvolver = new SplineLevelSetEvolver("Phi", (GridData)SplineLevelSet.GridDat);
                    lsUpdater.AddEvolver("Phi", SplineEvolver);
                    lsUpdater.AddLevelSetParameter("Phi", new LevelSetVelocity("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters));
                    break;
                case LevelSetEvolution.None:
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
                    break;
                default:
                    throw new NotImplementedException();
            }
            return lsUpdater;
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            //Update Calls
            dt = GetFixedTimestep();
            Console.WriteLine($"Starting timesetp {TimestepNo}, dt={dt}");
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            Console.WriteLine($"done with timestep {TimestepNo}");
            return dt;
        }
    }
}
