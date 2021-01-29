using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Extension of the <see cref="XNSE"/>-solver for additional heat transfer.
    /// (The 'F' stands for Fourier equation, i.e. Heat equation.)
    /// </summary>
    public class XNSFE : XNSE {

        private void AddXHeatMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            int D = this.GridData.SpatialDimension;

            int pTemp = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Temperature].Degree;
            // configuration for Temperature
            var confTemp = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pTemp }, //Math.Max(1, pTemp - iLevel) },
                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Temperature) }
            };
            configsLevel.Add(confTemp);

            // configuration for auxiliary heat flux
            if (this.Control.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                int pFlux;
                if (this.Control.FieldOptions.TryGetValue("HeatFlux*", out FieldOpts f)) {
                    pFlux = f.Degree;
                } else if (this.Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.HeatFluxX, out FieldOpts f1)) {
                    pFlux = f1.Degree;
                } else {
                    throw new Exception("MultigridOperator.ChangeOfBasisConfig: Degree of HeatFlux not found");
                }
                for (int d = 0; d < D; d++) {
                    var confHeatFlux = new MultigridOperator.ChangeOfBasisConfig() {
                        DegreeS = new int[] { pFlux }, // Math.Max(1, pFlux - iLevel) },
                        mode = MultigridOperator.Mode.Eye,
                        VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.HeatFluxVectorComponent(d)) }
                    };
                    configsLevel.Add(confHeatFlux);
                }
            }
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            base.AddMultigridConfigLevel(configsLevel);
            AddXHeatMultigridConfigLevel(configsLevel);
        }

        protected override void DefineSystem(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            base.DefineSystem(D, opFactory, lsUpdater);
            AddXHeat(D, opFactory, lsUpdater);
        }

        private void AddXHeat(int D, OperatorFactory opFactory, LevelSetUpdater lsUpdater) {
            int quadOrder = QuadOrder();
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            ThermalMultiphaseBoundaryCondMap heatBoundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            // add Heat equation components
            // ============================
            opFactory.AddEquation(new Heat("A", lsUpdater.Tracker, D, heatBoundaryMap, config));
            opFactory.AddEquation(new Heat("B", lsUpdater.Tracker, D, heatBoundaryMap, config));

            if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                for (int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new HeatFlux("A", d, lsUpdater.Tracker, D, heatBoundaryMap, config));
                    opFactory.AddEquation(new HeatFlux("B", d, lsUpdater.Tracker, D, heatBoundaryMap, config));
                    opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, heatBoundaryMap, lsUpdater.Tracker, config));
                }
            }
            opFactory.AddEquation(new HeatInterface("A", "B", D, heatBoundaryMap, lsUpdater.Tracker, config));
            opFactory.AddCoefficient(new EvapMicroRegion());
            if (config.prescribedMassflux != null)
                opFactory.AddCoefficient(new PrescribedMassFlux(config));

            // add Evaporation at Interface components, heads-up depends only on parameters
            // ============================
            if (config.isEvaporation) {
                //opFactory.AddParameter(new Temperature0());
                //opFactory.AddParameter(new HeatFlux0(D, lsUpdater.Tracker, config));
                var MassFluxExt = new MassFluxExtension_Evaporation(config);
                opFactory.AddParameter(MassFluxExt);
                lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, MassFluxExt);

                opFactory.AddEquation(new HeatInterface_MassFlux("A", "B", D, heatBoundaryMap, lsUpdater.Tracker, config));

                for (int d = 0; d < D; ++d)
                    opFactory.AddEquation(new InterfaceNSE_MassFlux("A", "B", D, d, lsUpdater.Tracker, config));

                if (config.isContinuity)
                    opFactory.AddEquation(new InterfaceContinuity_MassFlux("A", "B", D, lsUpdater.Tracker, config));
            }
        }


        protected override ILevelSetParameter GetLevelSetVelocity(int iLevSet) {
            if(iLevSet == 0) {
                // Main Difference to base implementation:
                //var levelSetVelocity = new LevelSetVelocityEvaporative("Phi", GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters, new XNSFE_OperatorConfiguration(Control);
                var levelSetVelocity = new LevelSetVelocityGeneralNonMaterial(VariableNames.LevelSetCG, GridData.SpatialDimension, VelocityDegree(), Control.InterVelocAverage, Control.PhysicalParameters);
                return levelSetVelocity;
            } else {
                return base.GetLevelSetVelocity(iLevSet);
            }
        }

        /*

        protected override LevelSetUpdater InstantiateLevelSetUpdater() {
            int levelSetDegree = Control.FieldOptions["Phi"].Degree;
            LevelSetUpdater lsUpdater;

           
            
            switch(Control.Option_LevelSetEvolution) {
                case LevelSetEvolution.Fourier: {
                    if(Control.EnforceLevelSetConservation) {
                        throw new NotSupportedException("mass conservation correction currently not supported");
                    }
                    FourierLevelSet fourrierLevelSet = new FourierLevelSet(Control.FourierLevSetControl, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                    fourrierLevelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, fourrierLevelSet, VariableNames.LevelSetCG);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                    break;
                }
                case LevelSetEvolution.FastMarching: {
                    LevelSet levelSet = new LevelSet(new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                    levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                    var fastMarcher = new FastMarchingEvolver(VariableNames.LevelSetCG, QuadOrder(), levelSet.GridDat.SpatialDimension);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet, VariableNames.LevelSetCG);
                    lsUpdater.AddEvolver(VariableNames.LevelSetCG, fastMarcher);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                    break;
                }
                case LevelSetEvolution.StokesExtension: {

                    throw new NotImplementedException("todo");
                }
                case LevelSetEvolution.SplineLS: {
                    int nodeCount = 30;
                    Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                    SplineLevelSet SplineLevelSet = new SplineLevelSet(Control.Phi0Initial, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG, nodeCount);
                    var SplineEvolver = new SplineLevelSetEvolver(VariableNames.LevelSetCG, (GridData)SplineLevelSet.GridDat);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, SplineLevelSet, VariableNames.LevelSetCG);
                    lsUpdater.AddEvolver(VariableNames.LevelSetCG, SplineEvolver);
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                    break;
                }
                case LevelSetEvolution.None: {
                    LevelSet levelSet1 = new LevelSet(new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                    levelSet1.ProjectField(Control.InitialValues_Evaluators["Phi"]);
                    
                    lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet1, VariableNames.LevelSetCG);
                    break;
                }
                default:
                throw new NotImplementedException();
            }
            
            return lsUpdater;
        }

        */



    }
}