using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {
    class XNSFE : XCommon<XNSE_Control> {

        XNSE xnse = new XNSE();
        private void GetXNSEOperatorComponents(int D, OperatorFactory opFactory) {
            xnse.SetOperatorEquations(D, opFactory);
        }
        private void XNSEMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            xnse.MultigridConfigLevel(configsLevel);
        }

        XHeat xheat = new XHeat();
        private void GetXHEATOperatorComponents(int D, OperatorFactory opFactory) {
            xheat.SetOperatorEquations(D, opFactory);
        }
        private void XHEATMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {
            xheat.MultigridConfigLevel(configsLevel);
        }
        public override void MultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel) {

            XNSEMultigridConfigLevel(configsLevel);
            XHEATMultigridConfigLevel(configsLevel);
            
        }

        protected void GetXNSFEOperatorComponents(int D, OperatorFactory opFactory) {

            //base.GetOperatorComponents(D, opFactory);

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

        public override void SetOperatorEquations(int D, OperatorFactory opFactory) {
            GetXNSEOperatorComponents(D, opFactory);
            GetXHEATOperatorComponents(D, opFactory);
        }

        protected override int QuadOrder() {
            throw new NotImplementedException();
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {
            throw new NotImplementedException();
        }

        public override void SetOperatorParameter(int D, OperatorFactory opFactory) {
            throw new NotImplementedException();
        }

        public override void SetSpatialOperator(out XSpatialOperatorMk2 XOP, int D, OperatorFactory opFactory) {
            throw new NotImplementedException();
        }
    }
}
