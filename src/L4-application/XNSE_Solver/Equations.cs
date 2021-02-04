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
                    AddComponent(new ConvectionAtLevelSet_Consistency(d, D, lsTrk, -1, false, thermParams, sigma));
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
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension);
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        private void AddInterfaceNSE_MassFlux(int D, int d, LevelSetTracker lsTrk, XNSFE_OperatorConfiguration config) {

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // from XNSFE_OperatorComponents
            if (config.isTransport) {
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
                if (!config.isMovingMesh) {
                    // the following terms decode the condition at the interface (consider the similarity to the rankine hugoniot condition)
                    // for the moving mesh discretization this condition is already contained in the convective terms
                    // therefore we only need these terms when using splitting...
                    AddComponent(new MassFluxAtLevelSet_withMassFlux(d, D, lsTrk, physParams, config.isMovingMesh));
                    AddComponent(new ConvectionAtLevelSet_nonMaterialLLF_withMassFlux(d, D, lsTrk, physParams));
                    AddComponent(new ConvectionAtLevelSet_Consistency_withMassFlux(d, D, lsTrk, -1, false, physParams));
                } else {
                    AddComponent(new ConvectionAtLevelSet_MovingMesh_withMassFlux(d, D, lsTrk, physParams));
                }
            } else {
                //  ... and when the convective terms are turned off we still need the contribution below
                AddComponent(new MassFluxAtLevelSet_withMassFlux(d, D, lsTrk, physParams, config.isMovingMesh));
            }           

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_withMassFlux(lsTrk, dntParams.PenaltySafety, d, physParams));
            }            
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

            var divEvap = new DivergenceAtLevelSet_withEvaporation(D, lsTrk, -1, false, thermParams, config.getPhysParams.Sigma);
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

            var divEvap = new DivergenceAtLevelSet_withMassFlux(D, lsTrk, -1, false, physicalParameters);
            AddComponent(divEvap);
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    public class HeatInterface_MassFlux : SurfaceEquation {


        string codomainName;
        string phaseA, phaseB;
        public HeatInterface_MassFlux(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            LevelSetTracker LsTrk,
            IXHeat_Configuration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.HeatEquation;
            AddInterfaceHeatEq(dimension, boundaryMap, LsTrk, config);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);

            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;


        //Methode aus der XNSF_OperatorFactory
        void AddInterfaceHeatEq(
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            LevelSetTracker LsTrk,
            IXHeat_Configuration config) {

            PhysicalParameters physParams = config.getPhysParams;
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double capA = thermParams.rho_A * thermParams.c_A;
            double LFFA = dntParams.LFFA;
            double kA = thermParams.k_A;

            double capB = thermParams.rho_B * thermParams.c_B;
            double LFFB = dntParams.LFFB;
            double kB = thermParams.k_B;

            double Tsat = thermParams.T_sat;

            // convective part
            // ================
            if (thermParams.IncludeConvection) {
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(dimension));
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(dimension));
                if (config.isMovingMesh) {
                    AddComponent(new HeatConvectionAtLevelSet_MovingMesh_withMassflux(dimension, LsTrk, Tsat, config.getPhysParams, thermParams));
                } else {
                    AddComponent(new HeatConvectionAtLevelSet_LLF_withMassflux(dimension, LsTrk, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh, Tsat, physParams));
                }
            }

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {

                double penalty = dntParams.PenaltySafety;

                var Visc = new ConductivityAtLevelSet_withMassflux(LsTrk, kA, kB, penalty * 1.0, Tsat);
                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }

        }
    }
}
