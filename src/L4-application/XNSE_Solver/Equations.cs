using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.RheologyCommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSE_Solver {
    class InterfaceNSE_Evaporation : SurfaceEquation {

        string codomainName;
        string phaseA, phaseB;
        public InterfaceNSE_Evaporation(string phaseA,
            string phaseB,
            int dimension,
            int d,
            LevelSetTracker LsTrk,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.MomentumEquationComponent(d);
            AddInterfaceNSE_Evaporation(dimension, d, LsTrk, config);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.HeatFlux0Vector(dimension));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Temperature0);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension);
            AddCoefficient("EvapMicroRegion");
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        private void AddInterfaceNSE_Evaporation(int D, int d, LevelSetTracker lsTrk, XNSFE_OperatorConfiguration config) {

            PhysicalParameters physParams = config.getPhysParams;
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            double sigma = physParams.Sigma;

            // from XNSFE_OperatorComponents
            if (config.isTransport) {
                if (!config.isMovingMesh) {
                    AddComponent(new MassFluxAtInterface(d, D, lsTrk, thermParams, sigma, config.isMovingMesh));
                    AddComponent(new ConvectionAtLevelSet_nonMaterialLLF(d, D, lsTrk, thermParams, sigma));
                    AddComponent(new ConvectionAtLevelSet_Consistency(d, D, lsTrk, dntParams.ContiSign, dntParams.RescaleConti, thermParams, sigma));
                }
            } else {
                AddComponent(new MassFluxAtInterface(d, D, lsTrk, thermParams, sigma, config.isMovingMesh));
            }

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_withEvap(lsTrk, physParams.mu_A, physParams.mu_B, dntParams.PenaltySafety, d, thermParams, sigma));
            }

        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    class InterfaceNSE_MassFlux : SurfaceEquation {

        string codomainName;
        string phaseA, phaseB;
        public InterfaceNSE_MassFlux(string phaseA,
            string phaseB,
            int dimension,
            int d,
            LevelSetTracker LsTrk,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.MomentumEquationComponent(d);
            AddInterfaceNSE_MassFlux(dimension, d, LsTrk, config);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension);
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        private void AddInterfaceNSE_MassFlux(int D, int d, LevelSetTracker lsTrk, XNSFE_OperatorConfiguration config) {

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_withMassFlux(lsTrk, dntParams.PenaltySafety, d, physParams));
            }

            AddComponent(new MassFluxAtLevelSet_withMassFlux(d, D, lsTrk, physParams));
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    class InterfaceContinuity_Evaporation : SurfaceEquation {

        string codomainName;
        string phaseA, phaseB;
        public InterfaceContinuity_Evaporation(string phaseA,
            string phaseB,
            int dimension,
            LevelSetTracker LsTrk,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.ContinuityEquation;
            AddInterfaceContinuity_Evaporation(dimension, LsTrk, config);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.HeatFlux0Vector(dimension));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Temperature0);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Curvature);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension);
            AddCoefficient("EvapMicroRegion");
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        private void AddInterfaceContinuity_Evaporation(int D, LevelSetTracker lsTrk, XNSFE_OperatorConfiguration config) {

            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            var divEvap = new DivergenceAtLevelSet_withEvaporation(D, lsTrk, dntParams.ContiSign, dntParams.RescaleConti, thermParams, config.getPhysParams.Sigma);
            AddComponent(divEvap);
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    class InterfaceContinuity_MassFlux : SurfaceEquation {

        string codomainName;
        string phaseA, phaseB;
        public InterfaceContinuity_MassFlux(string phaseA,
            string phaseB,
            int dimension,
            LevelSetTracker LsTrk,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.ContinuityEquation;
            AddInterfaceContinuity_MassFlux(dimension, LsTrk, config);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension);
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        private void AddInterfaceContinuity_MassFlux(int D, LevelSetTracker lsTrk, XNSFE_OperatorConfiguration config) {

            PhysicalParameters physicalParameters = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            var divEvap = new DivergenceAtLevelSet_withMassFlux(D, lsTrk, dntParams.ContiSign, dntParams.RescaleConti, physicalParameters);
            AddComponent(divEvap);
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }
}
