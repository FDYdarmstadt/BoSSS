using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {
    class XNSFE : XNSE {


        protected override void MultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {

            base.MultigridConfigLevel(configsLevel);

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


        protected override void GetOperatorComponents(int D, OperatorFactory opFactory) {

            base.GetOperatorComponents(D, opFactory);

            int quadOrder = QuadOrder();
            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            ThermalMultiphaseBoundaryCondMap boundaryMap = new ThermalMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            // add Heat equation components
            // ============================
            opFactory.AddEquation(new Heat("A", LsTrk, D, boundaryMap, config));
            opFactory.AddEquation(new Heat("B", LsTrk, D, boundaryMap, config));
            if (config.conductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                for (int d = 0; d < D; ++d) {
                    opFactory.AddEquation(new HeatFlux("A", d, LsTrk, D, boundaryMap, config));
                    opFactory.AddEquation(new HeatFlux("B", d, LsTrk, D, boundaryMap, config));
                }
            }
            opFactory.AddEquation(new HeatInterface("A", "B", D, boundaryMap, LsTrk, config));


            // add Evaporation interface components
            // ====================================
            if (config.isEvaporation) {
                throw new NotImplementedException();
                //opFactory.AddEquation(new HeatInterfaceEvaporation("A", "B", d, D, boundaryMap, LsTrk, config));

                //if (config.isContinuity)
                //    opFactory.AddEquation(new HeatInterfaceContinuityEvaporation("A", "B", d, D, boundaryMap, LsTrk, config));

            }
        }

        
    }
}
