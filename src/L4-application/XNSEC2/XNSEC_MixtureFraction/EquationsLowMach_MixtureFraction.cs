using BoSSS.Application.XNSEC;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Continuity equation for the Low-Mach equations, where density varies locally with temperature and concentrations in the bulk.
    /// This class is used for calculating the system using the mixture fractions instead of the temperature and concentrations
    /// </summary>
    public class BulkContinuity_MF : BulkEquation {
        private string speciesName;
        private string codomainName;

        public BulkContinuity_MF(
            int D,
            string spcName,
            XNSEC_OperatorConfiguration config,
            IncompressibleBoundaryCondMap BcMap,
            MaterialLaw EoS,
            double dt) {
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;

            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MixtureFraction);
            speciesName = spcName;

            for (int d = 0; d < D; ++d) {
                var conti = new Solution.XNSECommon.Operator.Continuity.DivergenceInSpeciesBulk_CentralDifferenceNewton(spcName, d, BcMap, D, EoS, NoOfChemicalSpecies);
                AddComponent(conti);
            }

            //Temporal term contribution:
            //Implicit Euler:  d(rho) / dt = (rho ^ n_t - rho_(t - 1)) / delta t, n: newton iteration counter
            if (!config.isSteady && config.timeDerivativeConti_OK) {
                var drho_dt = new BoSSS.Solution.XNSECommon.LowMach_TimeDerivativeConti(spcName, EoS, dt, NoOfChemicalSpecies);
                AddComponent(drho_dt);
                AddParameter("Density_t0");
                AddParameter("Density_t00");
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 0.0;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// Low-Mach momentum equations in the bulk phase
    /// This class is used for calculating the system using the mixture fractions instead of the temperature and concentrations
    /// </summary>
    public class BulkNavierStokes_MF : BulkEquation {
        private string speciesName;
        private string codomainName;

        /// <summary>
        ///
        /// </summary>
        /// <param name="spcName"></param>
        /// <param name="d">
        /// Momentum component index
        /// </param>
        /// <param name="LsTrk"></param>
        /// <param name="D">
        /// Spatial dimension
        /// </param>
        /// <param name="boundaryMap"></param>
        /// <param name="config"></param>
        public BulkNavierStokes_MF(
            string spcName,
            int d,
            int D,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS) {
            double Reynolds = config.Reynolds;
            double Froude = config.Froude;
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;

            Vector gravityDirection = config.gravityDirection;
            speciesName = spcName;
            codomainName = EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));
            AddVariableNames(NSECommon.VariableNames.MixtureFraction);

            DoNotTouchParameters dntParams = config.getDntParams;
            double penalty = dntParams.PenaltySafety;

            // Convective term
            // =================

            if (config.physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustionConvectionInSpeciesBulk_LLF_Newton(spcName, D, boundaryMap, d, EoS, NoOfChemicalSpecies);
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
                AddComponent(conv);
            }

            // pressure gradient
            // =================

            var pres = new Solution.XNSECommon.Operator.Pressure.PressureInSpeciesBulk(d, boundaryMap, spcName);
            AddComponent(pres);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure);

            // viscous operator
            // ================
            var viscOption = ViscosityOption.VariableViscosityDimensionless;
            var visc = new Solution.XNSECommon.Operator.Viscosity.LowMachViscosityInSpeciesBulk_AllTerms(spcName, penalty, d, D, boundaryMap, viscOption, 1, Reynolds, EoS);
            AddComponent(visc);

            // Gravity
            //==================
            var gravity = new BoSSS.Solution.XNSECommon.LowMach_Gravity(spcName, gravityDirection, d, Froude, boundaryMap.PhysMode, EoS, NoOfChemicalSpecies);
            AddComponent(gravity);

            if (config.PlotAdditionalParameters) {
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Rho);

                var conv = new BoSSS.Solution.NSECommon.DummyParameter(1);
                AddComponent(conv);
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// Low-Mach mixture fractions equations in the bulk phase
    /// </summary>
    public class BulkMixtureFraction_MF : BulkEquation {
        private string speciesName;
        private string codomainName;

        public BulkMixtureFraction_MF(
            string spcName,
            int D,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS) {
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;
            double Reynolds = config.Reynolds;
            double Prandtl = config.Prandtl;

            speciesName = spcName;
            codomainName = EquationNames.MixtureFractionEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(NSECommon.VariableNames.MixtureFraction);

            DoNotTouchParameters dntParams = config.getDntParams;

            // Convection
            if (config.physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF(spcName, D, NoOfChemicalSpecies, boundaryMap, EoS, 0);
                AddComponent(conv);
            }
            // Diffusion
            double penalty = dntParams.PenaltySafety;
            var heatConduction = new LowMachMassFractionDiffusionBulk(spcName, penalty, boundaryMap, EoS, Reynolds, Prandtl, false);
            AddComponent(heatConduction);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// Heat conduction term of the energy equation for LowMach solver.
    /// \/ *(k \/ T)
    /// </summary>
    public class LowMachMassFractionDiffusionBulk : SIPDiffusionMixtureFraction, ISpeciesFilter {

        public LowMachMassFractionDiffusionBulk(string spcName, double PenaltyBase,
                                     IncompressibleBoundaryCondMap BcMap,
                                     MaterialLaw EoS,
                                     double Reynolds,
                                     double Prandtl,
                                     bool prmsOK) : base(PenaltyBase, BcMap, EoS, Reynolds, Prandtl, prmsOK) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }

    public class MixtureFractionInterface_MF : SurfaceEquation {
        private string codomainName;

        public MixtureFractionInterface_MF(XNSEC_OperatorConfiguration config, int D, MaterialLaw EoS_A, MaterialLaw EoS_B) {
            codomainName = EquationNames.MixtureFractionEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MixtureFraction);
            DoNotTouchParameters dntParams = config.getDntParams;

            PhysicalParameters physParams = config.getPhysParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double penalty = dntParams.PenaltySafety;

            // set components
            var divPen = new Interface_MixtureFractionDiffusivity_LowMach(D, EoS_A, config.Reynolds, config.Prandtl, penalty, FirstSpeciesName, SecondSpeciesName);
            AddComponent(divPen);
        }

        public override string FirstSpeciesName => "A";

        public override string SecondSpeciesName => "B";

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// Newtonian momentum equation, (fluid/fluid) interface part;
    /// Using the mixture fraction formulation
    /// This provides coupling of two phases/components of a multiphase flow.
    /// </summary>
    public class NSEInterface_MF : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        //Methode aus der XNSF_OperatorFactory
        public NSEInterface_MF(
            string phaseA,
            string phaseB,
            int d,
            int dimension,
            IncompressibleBoundaryCondMap boundaryMap,
            INSE_Configuration config,
            MaterialLaw EoS_A,
            MaterialLaw EoS_B,
            int NoOfChemComp) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.MomentumEquationComponent(d);
            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double LFFA = dntParams.LFFA;
            double LFFB = dntParams.LFFB;
            double muA = physParams.mu_A;
            double muB = physParams.mu_B;


            // convective operator
            // ===================
            if (physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.ConvectionAtLevelSet_LLF_Newton_LowMach(d, dimension, rhoA, rhoB, LFFA, LFFB, physParams.Material, boundaryMap, FirstSpeciesName, SecondSpeciesName, EoS_A, EoS_B, -1);
                AddComponent(conv);
            }
            // pressure gradient
            // =================
            if (config.isPressureGradient) { //OK
                var presLs = new Solution.XNSECommon.Operator.Pressure.PressureFormAtLevelSet(d, dimension);
                AddComponent(presLs);
            }

            // viscous operator
            // ================
            if (config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {
                double penalty = dntParams.PenaltySafety;
                switch (dntParams.ViscosityMode) {
                    case ViscosityMode.FullySymmetric:
                        AddComponent(new Solution.XNSECommon.Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(dimension, muA, muB, penalty, d, false));
                        break;

                    default:
                        throw new NotImplementedException();
                }
            }
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(dimension).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;

       
    }


    public class InterfaceContinuity_MF: SurfaceEquation {
        private string codomainName;

        public InterfaceContinuity_MF(INSE_Configuration config, int D, LevelSetTracker LsTrk, bool isMaterialInterface) {
            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));

            PhysicalParameters physParams = config.getPhysParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;

            // set components
            var divPen = new Solution.XNSECommon.Operator.Continuity.DivergenceAtLevelSetLowMach(D, LsTrk, rhoA, rhoB, isMaterialInterface, -1, true);
            AddComponent(divPen);
        }

        public override string FirstSpeciesName => "A";

        public override string SecondSpeciesName => "B";

        public override string CodomainName => codomainName;
    }


}