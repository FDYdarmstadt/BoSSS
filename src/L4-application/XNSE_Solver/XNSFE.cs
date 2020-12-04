using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver 
{
    class XNSFE : XNSE 
    {
        void AddXHeatMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel)
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

        public override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) 
        {
            base.AddMultigridConfigLevel(configsLevel);
            AddXHeatMultigridConfigLevel(configsLevel);
        }

        protected override void SetOperator(int D, OperatorFactory opFactory)
        {
            base.SetOperator(D, opFactory);
            AddXHeat(D, opFactory);
        }

        void AddXHeat(int D, OperatorFactory opFactory)
        {
            int quadOrder = QuadOrder();
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            ThermalMultiphaseBoundaryCondMap heatBoundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            // add Heat equation components
            // ============================
            opFactory.AddEquation(new Heat("A", lsUpdater.Tracker, D, heatBoundaryMap, config));
            opFactory.AddEquation(new Heat("B", lsUpdater.Tracker, D, heatBoundaryMap, config));

            if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP)
            {
                for (int d = 0; d < D; ++d)
                {
                    opFactory.AddEquation(new HeatFlux("A", d, lsUpdater.Tracker, D, heatBoundaryMap, config));
                    opFactory.AddEquation(new HeatFlux("B", d, lsUpdater.Tracker, D, heatBoundaryMap, config));
                    opFactory.AddEquation(new HeatFluxInterface("A", "B", D, d, heatBoundaryMap, lsUpdater.Tracker, config));
                }
            }
            opFactory.AddEquation(new HeatInterface("A", "B", D, heatBoundaryMap, lsUpdater.Tracker, config));
            opFactory.AddCoefficient(new EvapMicroRegion());

            // add Evaporation at Interface components, heads-up depends only on parameters
            // ============================
            if (config.isEvaporation) {

                opFactory.AddParameter(new Temperature0());
                opFactory.AddParameter(new HeatFlux0(D, lsUpdater.Tracker, config));
                var MassFluxExt = new MassFluxExtension(config);
                opFactory.AddParameter(MassFluxExt);
                lsUpdater.AddLevelSetParameter("Phi", MassFluxExt);

                for (int d = 0; d < D; ++d)
                    opFactory.AddEquation(new InterfaceNSE_Evaporation("A", "B", D, d, lsUpdater.Tracker, config));
                
                if (config.isContinuity)
                    opFactory.AddEquation(new InterfaceContinuity_Evaporation("A", "B", D, lsUpdater.Tracker, config));
            }
        }
        
    }
}
