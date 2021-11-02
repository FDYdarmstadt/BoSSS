using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using ilPSP;
using ilPSP.Utils;
using System;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Continuity equation for the Low-Mach equations, where density varies locally with temperature and concentrations in the bulk.
    /// </summary>
    public class LowMachContinuity : BulkEquation {
        private string speciesName;
        private string codomainName;

        public LowMachContinuity(
            int D,
            string spcName,
            XNSEC_OperatorConfiguration config,
            IncompressibleBoundaryCondMap BcMap,
            MaterialLaw EoS,
            double dt) {
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;

            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MassFractions(NoOfChemicalSpecies));

            speciesName = spcName;
            double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

            for(int d = 0; d < D; ++d) {
                var conti = new Solution.XNSECommon.Operator.Continuity.DivergenceInSpeciesBulk_CentralDifference(spcName, d, BcMap, D, EoS, NoOfChemicalSpecies);
                AddComponent(conti);
            }

            // manufactured solution
            if(config.manSolSource_OK) {
                var MS_Momentum = new BoSSS.Solution.XNSECommon.LowMach_ContiManSolution(spcName, -1, MolarMasses, PhysicsMode.Combustion, false, null);
                AddComponent(MS_Momentum);
            }

            //Temporal term contribution:
            //Implicit Euler:  d(rho) / dt = (rho ^ n_t - rho_(t - 1)) / delta t, n: newton iteration counter
            if(!config.isSteady && config.timeDerivativeConti_OK) {
                var drho_dt = new BoSSS.Solution.XNSECommon.LowMach_TimeDerivativeConti(spcName, EoS, dt, NoOfChemicalSpecies);
                AddComponent(drho_dt);
                AddParameter("Density_t0");
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 0.0;

        public override string CodomainName => codomainName;
    }

    public class InterfaceContinuityLowMach : SurfaceEquation {
        private string codomainName;
        public InterfaceContinuityLowMach(INSE_Configuration config, int D, LevelSetTracker LsTrk, bool isMaterialInterface) {
            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));

            PhysicalParameters physParams = config.getPhysParams;
            //DoNotTouchParameters dntParams = config.getDntParams;

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

    /// <summary>
    /// LowMach, Newtonian momentum equation, (fluid/fluid) interface part;
    /// This provides coupling of two phases/components of a multiphase flow, the bulk phases are defines through <see cref="NavierStokes"/>.
    /// </summary>
    public class NSEInterface_LowMach : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        //Methode aus der XNSF_OperatorFactory
        public NSEInterface_LowMach(
            string phaseA,
            string phaseB,
            int d,
            int dimension,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS_A,
            MaterialLaw EoS_B,
            bool isMovingMesh) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.MomentumEquationComponent(d);
            AddInterfaceNSE_LowMach(dimension, d, boundaryMap, config, isMovingMesh, EoS_A, EoS_B);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(dimension).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MassFractions(config.NoOfChemicalSpecies));
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;

        private void AddInterfaceNSE_LowMach(
            int dimension,
            int d,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            bool isMovingMesh,
            MaterialLaw EoS_A,
            MaterialLaw EoS_B
            ) {
            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;
            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double LFFA = dntParams.LFFA;
            double LFFB = dntParams.LFFB;
            double muA = physParams.mu_A;
            double muB = physParams.mu_B;

            // convective operator
            // ===================
            if(physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.ConvectionAtLevelSet_LLF_Newton_LowMach(d, dimension, rhoA, rhoB, LFFA, LFFB, physParams.Material, boundaryMap, isMovingMesh, FirstSpeciesName, SecondSpeciesName, EoS_A, EoS_B, NoOfChemicalSpecies);
                AddComponent(conv);
            }
            if(isMovingMesh && (physParams.IncludeConvection && config.isTransport == false)) {
                // if Moving mesh, we need the interface transport term somehow

                throw new NotImplementedException("Something missing here.");
            }

            // pressure gradient
            // =================
            if(config.isPressureGradient) {
                var presLs = new Solution.XNSECommon.Operator.Pressure.PressureFormAtLevelSet(d, dimension);
                AddComponent(presLs);
            }

            // viscous operator
            // ================
            if(config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {
                double penalty = dntParams.PenaltySafety;
                switch(dntParams.ViscosityMode) {
                    case ViscosityMode.Standard:
                    AddComponent(new Solution.XNSECommon.Operator.Viscosity.ViscosityAtLevelSet_Standard(muA, muB, penalty * 1.0, dimension, d, true));
                    break;

                    case ViscosityMode.TransposeTermMissing:
                    AddComponent(new Solution.XNSECommon.Operator.Viscosity.ViscosityAtLevelSet_Standard(muA, muB, penalty * 1.0, dimension, d, false));
                    break;

                    case ViscosityMode.FullySymmetric:
                    AddComponent(new Solution.XNSECommon.Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(dimension, muA, muB, penalty, d));
                    break;

                    default:
                    throw new NotImplementedException();
                }
            }
        }
    }

    /// <summary>
    /// Low-Mach momentum equations in the bulk phase
    /// </summary>
    public class LowMachNavierStokes : BulkEquation {
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
        public LowMachNavierStokes(
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
            AddVariableNames(NSECommon.VariableNames.Temperature);
            AddVariableNames(NSECommon.VariableNames.MassFractions(NoOfChemicalSpecies));
            //SpeciesId spcId = LsTrk.GetSpeciesId(spcName);
            DoNotTouchParameters dntParams = config.getDntParams;
            double penalty = dntParams.PenaltySafety;

            // Convective term
            // =================

            if(config.physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustionConvectionInSpeciesBulk_LLF(spcName, D, boundaryMap, d, EoS, NoOfChemicalSpecies);
                AddComponent(conv);

                //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d]);
                //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]);
            }

            AddParameter(BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Rho);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Mu);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.cp);


            // pressure gradient
            // =================
            if(config.isPressureGradient) {
                var pres = new Solution.XNSECommon.Operator.Pressure.PressureInSpeciesBulk(d, boundaryMap, spcName);
                AddComponent(pres);
            }

            // viscous operator
            // ================
            var viscOption = ViscosityOption.VariableViscosityDimensionless;
            AddCoefficient("SlipLengths");
            AddCoefficient("Reynolds");
            if(config.isViscous &&
                (spcName == "A" && !(config.physParams.mu_A == 0.0)) ||
                (spcName == "B" && !(config.physParams.mu_B == 0.0))) {
                var visc = new Solution.XNSECommon.Operator.Viscosity.LowMachViscosityInSpeciesBulk_AllTerms(spcName, penalty, d, D, boundaryMap, viscOption, 1, Reynolds, EoS);
                AddComponent(visc);
            }
            // Gravity
            //==================
            if(gravityDirection.L2Norm() >0.0) { 
            var gravityLowMach = new BoSSS.Solution.XNSECommon.LowMach_Gravity(spcName, gravityDirection, d, Froude, boundaryMap.PhysMode, EoS, NoOfChemicalSpecies);
            AddComponent(gravityLowMach);
            }
            // gravity & more general volume force
            // ================
            if(config.isGravity) {
                string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
                string gravityOfSpecies = gravity + "#" + SpeciesName;
                var gravityComponent = new Solution.XNSECommon.Operator.MultiPhaseSource(gravityOfSpecies, speciesName);
                AddComponent(gravityComponent);
                AddParameter(gravityOfSpecies);
            }
            if(config.isVolForce) {
                string volforce = BoSSS.Solution.NSECommon.VariableNames.VolumeForceVector(D)[d];
                string volforceOfSpecies = volforce + "#" + SpeciesName;
                var volforceComponent = new Solution.XNSECommon.Operator.MultiPhaseSource(volforceOfSpecies, speciesName);
                AddComponent(volforceComponent);
                AddParameter(volforceOfSpecies);
            }

            // Manufactured Solutions source
            //=========================================
            string direction = d == 0 ? "x" : "y";
            double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };
            bool rhoOne = false;
            if(config.manSolSource_OK) {
                var MS_Momentum = new BoSSS.Solution.XNSECommon.LowMach_MomentumManSolution(spcName, Reynolds, Froude, MolarMasses, direction, PhysicsMode.Combustion, rhoOne, null);
                AddComponent(MS_Momentum);
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }

    public class IdentityEquation : BulkEquation {
        private string speciesName;
        private string m_codomainName;

        public IdentityEquation(
            string spcName,
            string variableName,
            string codomainName
            ) {
            speciesName = spcName;
            m_codomainName = codomainName;
            AddVariableNames(variableName);

            // Add identity
            var identity = new BoSSS.Solution.XNSECommon.identityTerm(spcName, variableName);
            AddComponent(identity);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => m_codomainName;
    }

    /// <summary>
    /// Low-Mach momentum equations in the bulk phase
    /// </summary>
    public class LowMachEnergy : BulkEquation {
        private string speciesName;
        private string codomainName;

        public LowMachEnergy(
            string spcName,
            int D,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            MaterialLaw_MultipleSpecies EoS,
            double HeatReleaseFactor,
            double[] ReactionRateConstants,
            double[] MolarMasses,
            double TRef,
            double cpRef,
            double dt,
            BoSSS.Solution.NSECommon.SIPDiffusionTemperature.ThermalWallType myThermalWallType
            ) {
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;
            double Reynolds = config.Reynolds;
            double Prandtl = config.Prandtl;

            speciesName = spcName;
            codomainName = EquationNames.HeatEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));
            AddVariableNames(NSECommon.VariableNames.Temperature);
            AddVariableNames(NSECommon.VariableNames.MassFractions(NoOfChemicalSpecies));

            DoNotTouchParameters dntParams = config.getDntParams;

            // Convection
            if(config.getPhysParams.IncludeConvection) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF(spcName, D, NoOfChemicalSpecies, boundaryMap, EoS, 0);
                AddComponent(conv);
            }
            //Heat Conduction
            double penalty = dntParams.PenaltySafety;
            var heatConduction = new BoSSS.Solution.XNSECommon.Operator.Viscosity.LowMachEnergyConductionBulk(spcName, penalty, boundaryMap, EoS, Reynolds, Prandtl, false, myThermalWallType);
            AddComponent(heatConduction);

            if(config.includeReactionTerms) {
                var ReactionTerm = new BoSSS.Solution.XNSECommon.LowMach_HeatSource(spcName, HeatReleaseFactor, ReactionRateConstants, MolarMasses, EoS, TRef, cpRef, config.VariableReactionRateParameters);
                AddComponent(ReactionTerm);
                AddParameter("kReact");

            }

            //Temporal term contribution:
            //Implicit Euler:  -d(p0) / dt = -(p0^{n}_{t} - p0_(t - 1)) / delta t, n: newton iteration counter
            if(!config.isSteady && config.timeDerivativeEnergyp0_OK) {
                //-(p0^{n}_{t} - p0_(t - 1)) / delta t
                var dtp0_dt = new BoSSS.Solution.XNSECommon.LowMach_TimeDerivativep0(spcName, dt);
                AddComponent(dtp0_dt);
                //AddParameter(BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure);
                //AddParameter(BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure + "_t0");
                //AddParameter("dp0dt");
            }
            AddParameter("dp0dt");
            // Manufactured Solutions source
            //=========================================
            if(config.manSolSource_OK) {
                double Schmidt = Prandtl;
                double[] StoichiometricCoefficients = null;
                var MS_Energy = new BoSSS.Solution.XNSECommon.LowMach_ScalarManSolution(spcName, HeatReleaseFactor, Reynolds, Prandtl, Schmidt, StoichiometricCoefficients, ReactionRateConstants, MolarMasses, EoS, "Temperature", PhysicsMode.Combustion);
                AddComponent(MS_Energy);
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }

    public class HeatInterface_LowMach : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        //Methode aus der XNSF_OperatorFactory
        public HeatInterface_LowMach(
            string phaseA,
            string phaseB,
            int dimension,
            IncompressibleBoundaryCondMap boundaryMap,
            IXHeat_Configuration config) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.HeatEquation;
            AddInterfaceHeatEq(dimension, boundaryMap, config);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            //AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;

        //Methode aus der XNSF_OperatorFactory
        private void AddInterfaceHeatEq(
            int dimension,
            IncompressibleBoundaryCondMap boundaryMap,
            IXHeat_Configuration config) {
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
            if(thermParams.IncludeConvection) {
                Console.WriteLine("include heat convection");
                DefineConvective(dimension, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh);
            }

            // viscous operator (laplace)
            // ==========================

            double penalty = dntParams.PenaltySafety;

            //var Visc = new ConductivityAtLevelSet(LsTrk, kA, kB, penalty * 1.0, Tsat);
            //var Visc = new ConductivityAtLevelSet_material(dimension, kA, kB, penalty * 1.0, Tsat);
            //AddComponent(Visc);
        }

        protected virtual void DefineConvective(int dimension, double capA, double capB, double LFFA, double LFFB, IncompressibleBoundaryCondMap boundaryMap, bool isMovingMesh) {
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(dimension));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(dimension));
            //AddComponent(new HeatConvectionAtLevelSet_LLF(dimension, LsTrk, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh, Tsat));
            //AddComponent(new HeatConvectionAtLevelSet_LLF_material(dimension, capA, capB, LFFA, LFFB, boundaryMap, isMovingMesh));            
            //AddComponent(new HeatConvectionAtLevelSet_LLF_material_Newton_Hamiltonian(dimension, capA, capB, LFFA, LFFB, boundaryMap, isMovingMesh, FirstSpeciesName, SecondSpeciesName));
        }
    }

    /// <summary>
    /// Low-Mach Mass fraction in the bulk phase
    /// </summary>
    public class LowMachMassFraction : BulkEquation {
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
        public LowMachMassFraction(
            string spcName,
            int D,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            MaterialLaw_MultipleSpecies EoS,
            int ChemicalSpeciesCounter,
            double[] ReactionRateConstants,
            double[] StoichiometricCoefficients,
            double[] MolarMasses) {
            double Reynolds = config.Reynolds;
            double Prandtl = config.Prandtl;
            double[] Lewis = config.Lewis;
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;

            speciesName = spcName;
            codomainName = EquationNames.SpeciesMassBalanceName(ChemicalSpeciesCounter);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));
            AddVariableNames(NSECommon.VariableNames.Temperature);
            AddVariableNames(NSECommon.VariableNames.MassFractions(NoOfChemicalSpecies));

            DoNotTouchParameters dntParams = config.getDntParams;

            // Convection of species i
            //===============
            if(config.getPhysParams.IncludeConvection) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF(spcName, D, NoOfChemicalSpecies, boundaryMap, EoS, ChemicalSpeciesCounter + 1);
                AddComponent(conv);
            }
            // Mass diffusion of species i
            //===============
            double penalty = dntParams.PenaltySafety;
            var massDiffusion_i = new BoSSS.Solution.XNSECommon.Operator.Viscosity.LowMachSpeciesBalanceDiffusionBulk(spcName, penalty, boundaryMap, EoS, Reynolds, Prandtl, Lewis, ChemicalSpeciesCounter, NoOfChemicalSpecies);
            AddComponent(massDiffusion_i);

            // Mass source/sink term of species i
            //===================================
            if(config.includeReactionTerms) {
                var massReaction_i = new BoSSS.Solution.XNSECommon.LowMach_MassFractionSource(spcName, ReactionRateConstants, StoichiometricCoefficients, MolarMasses, EoS, NoOfChemicalSpecies, ChemicalSpeciesCounter, 300, 1.0, config.VariableReactionRateParameters);
                AddComponent(massReaction_i);
            }

            if(config.manSolSource_OK) {
                double Schmidt = Prandtl;
                double HeatReleaseFactor = -10000; // not needed
                var MS_Energy = new BoSSS.Solution.XNSECommon.LowMach_ScalarManSolution(spcName, HeatReleaseFactor, Reynolds, Prandtl, Schmidt, StoichiometricCoefficients, ReactionRateConstants, MolarMasses, EoS, "MassFraction", PhysicsMode.Combustion, ChemicalSpeciesCounter, true);
                AddComponent(MS_Energy);
            }
        }

        public override string SpeciesName => speciesName;
        public override double MassScale => 1.0;
        public override string CodomainName => codomainName;
    }
}