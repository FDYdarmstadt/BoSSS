using BoSSS.Application.XNSEC;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Continuity equation for the Low-Mach equations, where density varies locally with temperature and concentrations in the bulk.
    /// This class is used for calculating the system using the mixture fractions instead of the temperature and concentrations
    /// </summary>
    public class LowMachContinuity_MixtureFractions : BulkEquation {
        private string speciesName;
        private string codomainName;

        public LowMachContinuity_MixtureFractions(
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

            for(int d = 0; d < D; ++d) {
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
    public class LowMachNavierStokes_MixtureFractions : BulkEquation {
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
        public LowMachNavierStokes_MixtureFractions(
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

            if(config.physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustionConvectionInSpeciesBulk_LLF_Newton(spcName, D, boundaryMap, d, EoS, NoOfChemicalSpecies);
                AddComponent(conv);
                //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d]);
                //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]);
            }

            // pressure gradient
            // =================

            var pres = new Solution.XNSECommon.Operator.Pressure.PressureInSpeciesBulk(d, boundaryMap, spcName);
            AddComponent(pres);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure);

            // viscous operator
            // ================
            var viscOption = ViscosityOption.VariableViscosityDimensionless;
            //AddCoefficient("SlipLengths");
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
    public class LowMachMixtureFraction : BulkEquation {
        private string speciesName;
        private string codomainName;

        public LowMachMixtureFraction(
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
            var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF(spcName, D, NoOfChemicalSpecies, boundaryMap, EoS, 0);
            AddComponent(conv);

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
}