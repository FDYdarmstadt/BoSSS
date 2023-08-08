using BoSSS.Application.XNSEC;
using BoSSS.Application.XNSFE_Solver;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Linq;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Continuity equation for the Low-Mach equations, where density varies locally with temperature and mass fractions in the bulk.
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
            double dt,
            Func<double[], double, double> ManSol,
            NonLinearSolverCode NonLinSolverCode) {
            int NoOfChemicalSpecies = config.NoOfChemicalSpecies;

            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MassFractions(NoOfChemicalSpecies));

            speciesName = spcName;

            //Temporal term contribution:
            //Implicit Euler:  d(rho) / dt = (rho ^ n_t - rho_(t - 1)) / delta t, n: newton iteration counter
            if (!config.isSteady && config.timeDerivativeConti_OK) {
                var drho_dt = new BoSSS.Solution.XNSECommon.LowMach_TimeDerivativeConti(spcName, EoS, dt, NoOfChemicalSpecies);
                AddComponent(drho_dt);
                AddParameter("Density_t0");
                AddParameter("Density_t00");
            }

            // Divergence term
            // \/ . (\rho \vec{u})
            for (int d = 0; d < D; ++d) {
                var contiNewton = new Solution.XNSECommon.Operator.Continuity.DivergenceInSpeciesBulk_CentralDifferenceNewton(spcName, d, BcMap, D, EoS, NoOfChemicalSpecies);
                AddComponent(contiNewton);
            }

            // manufactured solution
            if (config.manSolSource_OK) {
                var MS_conti = new BoSSS.Solution.XNSECommon.LowMach_ManSolution(spcName, ManSol);
                AddComponent(MS_conti);
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
    /// same as <see cref="InterfaceContinuity_Evaporation"/> but using Newton solver compatible components
    /// </summary>
    public class InterfaceContinuity_Evaporation_Newton_LowMach : InterfaceContinuity_Evaporation {

        public InterfaceContinuity_Evaporation_Newton_LowMach(string phaseA,
            string phaseB,
            int dimension,
            XNSFE_OperatorConfiguration config) : base(phaseA, phaseB, dimension, config) {
        }

        protected override void AddInterfaceContinuity_Evaporation(int D, XNSFE_OperatorConfiguration config) {
            AddComponent(new DivergenceAtLevelSet_Evaporation_StrongCoupling_LowMach(D, -1, false, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
        }
    }

    /// <summary>
    /// same as <see cref="InterfaceContinuity_Evaporation"/> but using Newton solver compatible components
    /// </summary>
    public class InterfaceContinuity_Evaporation_Newton_LowMach_MF : InterfaceContinuity_Evaporation {

        public InterfaceContinuity_Evaporation_Newton_LowMach_MF(string phaseA,
            string phaseB,
            int dimension,
            XNSFE_OperatorConfiguration config) : base(phaseA, phaseB, dimension, config) {
        }

        protected override void AddInterfaceContinuity_Evaporation(int D, XNSFE_OperatorConfiguration config) {
            AddComponent(new DivergenceAtLevelSet_Evaporation_StrongCoupling_LowMach_MF(D, -1, false, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
        }
    }

    public class InterfaceNSE_Evaporation_MF : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        public InterfaceNSE_Evaporation_MF(string phaseA,
            string phaseB,
            int dimension,
            int d,
            XNSFE_OperatorConfiguration config) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;
            int D = dimension;
            codomainName = EquationNames.MomentumEquationComponent(d);

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            if (config.isTransport) {
                // the following terms decode the condition at the interface (consider the similarity to the rankine hugoniot condition)
                AddComponent(new ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton_MF(d, D, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
                AddComponent(new ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton_MF(d, D, -1, false, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
            }

            if (config.isRecoilPressure) {
                AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling_MF(d, D, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
            }

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling_MF(dntParams.PenaltySafety, d, D, config.getThermParams, physParams, FirstSpeciesName, SecondSpeciesName));
            }

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d));
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// LowMach, Newtonian momentum equation, (fluid/fluid) interface part;
    /// This provides coupling of two phases/components of a multiphase flow, the bulk phases are defined through <see cref="NavierStokes"/>.
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
            if (physParams.IncludeConvection && config.isTransport) {
                var conv = new Solution.XNSECommon.Operator.Convection.ConvectionAtLevelSet_LLF_Newton_LowMach(d, dimension, rhoA, rhoB, LFFA, LFFB, physParams.Material, boundaryMap, FirstSpeciesName, SecondSpeciesName, EoS_A, EoS_B, NoOfChemicalSpecies);
                AddComponent(conv);
            }

            // pressure gradient
            // =================
            if (config.isPressureGradient) {
                var presLs = new Solution.XNSECommon.Operator.Pressure.PressureFormAtLevelSet(d, dimension);
                AddComponent(presLs);
            }

            // viscous operator
            // ================
            if (config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {
                double penalty = dntParams.PenaltySafety;
                switch (dntParams.ViscosityMode) {
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
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(dimension).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MassFractions(config.NoOfChemicalSpecies));
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    public class InterfaceNSE_Evaporation_LowMach : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        public InterfaceNSE_Evaporation_LowMach(string phaseA,
            string phaseB,
            int dimension,
            int d,
            XNSFE_OperatorConfiguration config) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.MomentumEquationComponent(d);
            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;
            int D = dimension;

            if (config.isTransport) {
                    // the following terms decode the condition at the interface (consider the similarity to the rankine hugoniot condition)
                    
                    if (config.isRecoilPressure) {
                        AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling(d, D, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    }
                    AddComponent(new ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton(d, D, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
                    AddComponent(new ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton(d, D, -1, false, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
                 
            } else {
                //  ... and when the convective terms are turned off we still need the contribution below
                if (config.isRecoilPressure) {
                    AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling(d, D, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                }
            }

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling_LowMach(dntParams.PenaltySafety, d, D, config.getThermParams, physParams, FirstSpeciesName, SecondSpeciesName));
            } 

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d));
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");
        }

      

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// Low-Mach momentum equations in the bulk phase
    /// </summary>
    public class LowMachMomentumEquations : BulkEquation {
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
        public LowMachMomentumEquations(
            string spcName,
            int d,
            int D,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS,
            Func<double[], double, double> ManSol,
            NonLinearSolverCode NonLinSolverCode) {
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

            if (config.physParams.IncludeConvection && config.isTransport) {
                if (NonLinSolverCode == NonLinearSolverCode.Newton) {
                    var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustionConvectionInSpeciesBulk_LLF_Newton(spcName, D, boundaryMap, d, EoS, NoOfChemicalSpecies);
                    AddComponent(conv);
                    AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]); // we still use this for penalty computation
                } else if (NonLinSolverCode == NonLinearSolverCode.Picard) {
                    var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustionConvectionInSpeciesBulk_LLF(spcName, D, boundaryMap, d, EoS, NoOfChemicalSpecies);
                    AddComponent(conv);
                    AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d]);
                    AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]);
                    AddParameter(BoSSS.Solution.NSECommon.VariableNames.Temperature0);
                    AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFraction0_0);
                }
            }

            AddParameter(BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure);

            // pressure gradient
            // =================
            if (config.isPressureGradient) {
                var pres = new Solution.XNSECommon.Operator.Pressure.PressureInSpeciesBulk(d, boundaryMap, spcName);
                AddComponent(pres);
            }

            // viscous operator
            // ================
            var viscOption = ViscosityOption.VariableViscosityDimensionless;
            AddCoefficient("SlipLengths");
            AddCoefficient("Reynolds");
            if (config.isViscous &&
                (spcName == "A" && !(config.physParams.mu_A == 0.0)) ||
                (spcName == "B" && !(config.physParams.mu_B == 0.0))) {
                var visc = new Solution.XNSECommon.Operator.Viscosity.LowMachViscosityInSpeciesBulk_AllTerms(spcName, penalty, d, D, boundaryMap, viscOption, 1, Reynolds, EoS);
                AddComponent(visc);
            }
            // Gravity
            //==================
            if (gravityDirection.L2Norm() > 0.0) {
                var gravityLowMach = new BoSSS.Solution.XNSECommon.LowMach_Gravity(spcName, gravityDirection, d, Froude, boundaryMap.PhysMode, EoS, NoOfChemicalSpecies);
                AddComponent(gravityLowMach);
            }
            // gravity & more general volume force
            // ================
            if (config.isGravity) {
                string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
                string gravityOfSpecies = gravity + "#" + SpeciesName;
                var gravityComponent = new Solution.XNSECommon.Operator.MultiPhaseSource(gravityOfSpecies, speciesName);
                AddComponent(gravityComponent);
                AddParameter(gravityOfSpecies);
            }
            if (config.isVolForce) {
                string volforce = BoSSS.Solution.NSECommon.VariableNames.VolumeForceVector(D)[d];
                string volforceOfSpecies = volforce + "#" + SpeciesName;
                var volforceComponent = new Solution.XNSECommon.Operator.MultiPhaseSource(volforceOfSpecies, speciesName);
                AddComponent(volforceComponent);
                AddParameter(volforceOfSpecies);
            }

            // Manufactured Solutions source
            //=========================================
            string direction = d == 0 ? "x" : "y";
            bool rhoOne = false;

            if (config.manSolSource_OK) {
                var MS_Momentum = new BoSSS.Solution.XNSECommon.LowMach_ManSolution(spcName, ManSol);

                AddComponent(MS_Momentum);
            }

            if (config.PlotAdditionalParameters) {
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Rho);
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.Mu);
                AddParameter(BoSSS.Solution.NSECommon.VariableNames.cp);
                var conv = new BoSSS.Solution.NSECommon.DummyParameter();
                AddComponent(conv);
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
            BoSSS.Solution.NSECommon.SIPDiffusionTemperature.ThermalWallType myThermalWallType,
            Func<double[], double, double> ManSol
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
            if (config.physParams.IncludeConvection) {
                var conv = new Solution.XNSECommon.Operator.Convection.LowMachCombustion_ScalarConvectionInSpeciesBulk_LLF(spcName, D, NoOfChemicalSpecies, boundaryMap, EoS, 0);
                AddComponent(conv);
            }
            //Heat Conduction
            double penalty = dntParams.PenaltySafety;
            var heatConduction = new BoSSS.Solution.XNSECommon.Operator.Viscosity.LowMachEnergyConductionBulk(spcName, penalty, boundaryMap, EoS, Reynolds, Prandtl, false, myThermalWallType);
            AddComponent(heatConduction);

            if (config.includeReactionTerms) {
                Console.WriteLine("including reactive terms!!!!");
                var ReactionTerm = new BoSSS.Solution.XNSECommon.LowMach_HeatSource(spcName, HeatReleaseFactor, ReactionRateConstants, MolarMasses, EoS, TRef, cpRef, config.VariableReactionRateParameters, config.NoOfChemicalSpecies);
                AddComponent(ReactionTerm);
                AddParameter("kReact");
            }

            //Temporal term contribution:
            //Implicit Euler:  -d(p0) / dt = -(p0^{n}_{t} - p0_(t - 1)) / delta t, n: newton iteration counter
            if (!config.isSteady && config.timeDerivativeEnergyp0_OK) {
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
            if (config.manSolSource_OK) {
                double Schmidt = Prandtl;
                double[] StoichiometricCoefficients = null;
                var MS_Energy = new BoSSS.Solution.XNSECommon.LowMach_ManSolution(spcName, ManSol);
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
            IXHeat_Configuration config,
            MaterialLaw EoS_A,
            MaterialLaw EoS_B, double Reynolds, double Prandtl) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.HeatEquation;
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
                Console.WriteLine("include heat convection");
                DefineConvective(dimension, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh);
            }

            // viscous operator (laplace)
            // ==========================

            double penalty = dntParams.PenaltySafety;

            var Visc = new ConductivityAtLevelSet_material(dimension, kA, kB, penalty * 1.0, Tsat, FirstSpeciesName, SecondSpeciesName);
                //var Visc = new Interface_EnergyConduction_LowMach(dimension, EoS_A, EoS_B, Reynolds, Prandtl, penalty, Tsat, FirstSpeciesName, SecondSpeciesName);
                AddComponent(Visc); AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            //AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;

        protected virtual void DefineConvective(int dimension, double capA, double capB, double LFFA, double LFFB, IncompressibleBoundaryCondMap boundaryMap, bool isMovingMesh) {
            //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(dimension));
            //AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(dimension));
            //AddComponent(new HeatConvectionAtLevelSet_LLF(dimension, LsTrk, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh, Tsat));
            //AddComponent(new HeatConvectionAtLevelSet_LLF_material(dimension, capA, capB, LFFA, LFFB, boundaryMap, isMovingMesh));
            //AddComponent(new HeatConvectionAtLevelSet_LLF_material_Newton_Hamiltonian(dimension, capA, capB, LFFA, LFFB, boundaryMap, isMovingMesh, FirstSpeciesName, SecondSpeciesName));
        }
    }

    public class HeatInterface_Evaporation_LowMach : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        public HeatInterface_Evaporation_LowMach(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            XNSFE_OperatorConfiguration config, MaterialLaw EoS_A, MaterialLaw EoS_B, double Reynolds, double Prandtl, double ValueAtInterface) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.HeatEquation;
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

            // convective part
            // ================
            if (thermParams.IncludeConvection) {
                AddComponent(new HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian_LowMach(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, ValueAtInterface, thermParams, FirstSpeciesName, SecondSpeciesName));
            }
            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                double penalty = dntParams.PenaltySafety;
                //var Visc = new ConductivityAtLevelSet_withMassflux(dimension, kA, kB, penalty * 1.0, Tsat, config.isMaterialAtContactLine);
                var Visc = new Interface_EnergyConduction_LowMach(dimension, EoS_A, EoS_B, Reynolds, Prandtl, penalty, ValueAtInterface, FirstSpeciesName, SecondSpeciesName);

                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);

            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    public class SpeciesMassTransferInterface_Evaporation_LowMach : SurfaceEquation {
        private string codomainName;
        private string phaseA, phaseB;

        public SpeciesMassTransferInterface_Evaporation_LowMach(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            XNSFE_OperatorConfiguration config,
            MaterialLaw EoS_A,
            MaterialLaw EoS_B,
            double Reynolds,
            double Prandtl
            , double ValueAtInterface,
            int Y_component) : base() {
            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.SpeciesMassBalanceName(Y_component);
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

            // convective part
            // ================
            if (thermParams.IncludeConvection) {
                AddComponent(new SpeciesConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, ValueAtInterface, thermParams, FirstSpeciesName, SecondSpeciesName, Y_component));
            }

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                double penalty = dntParams.PenaltySafety;
                var Visc = new Interface_MassDiffusion_LowMach(dimension, EoS_A, EoS_B, Reynolds, Prandtl, penalty, ValueAtInterface, FirstSpeciesName, SecondSpeciesName, Y_component);
                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MassFraction_n(Y_component));

            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
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
            double[] MolarMasses,
            XNSEC_Control control) {
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
            if (config.getPhysParams.IncludeConvection) {
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
            if (config.includeReactionTerms) {

                var massReaction_i = new BoSSS.Solution.XNSECommon.LowMach_MassFractionSource(spcName, ReactionRateConstants, StoichiometricCoefficients, MolarMasses, EoS, NoOfChemicalSpecies, ChemicalSpeciesCounter, 300, 1.0, config.VariableReactionRateParameters);
                AddComponent(massReaction_i);
            }

            if (config.manSolSource_OK) {
                Func<double[], double, double> MF;
                switch (ChemicalSpeciesCounter) {
                    case 0:
                        MF = control.ManufacturedSolution_Species0;
                        break;

                    case 1:
                        MF = control.ManufacturedSolution_Species1;
                        break;

                    case 2:
                        MF = control.ManufacturedSolution_Species2;
                        break;

                    case 3:
                        MF = control.ManufacturedSolution_Species3;
                        break;

                    default:
                        MF = control.ManufacturedSolution_Species3;
                        break;
                }
                var MS = new BoSSS.Solution.XNSECommon.LowMach_ManSolution(spcName, MF);
                AddComponent(MS);
            }
        }

        public override string SpeciesName => speciesName;
        public override double MassScale => 1.0;
        public override string CodomainName => codomainName;
    }

    //public class MassFractionInterface_Evaporation : SurfaceEquation {
    //    private string codomainName;
    //    private string phaseA, phaseB;

    //    public MassFractionInterface_Evaporation(
    //        string phaseA,
    //        string phaseB,
    //        int dimension,
    //        int chemicalComponentIndex,
    //        XNSFE_OperatorConfiguration config) : base() {
    //        this.phaseA = phaseA;
    //        this.phaseB = phaseB;

    //        codomainName = EquationNames.SpeciesMassBalanceName(chemicalComponentIndex);
    //        PhysicalParameters physParams = config.getPhysParams;
    //        ThermalParameters thermParams = config.getThermParams;
    //        DoNotTouchParameters dntParams = config.getDntParams;

    //        // set species arguments
    //        double capA = thermParams.rho_A;
    //        double LFFA = dntParams.LFFA;
    //        double rhoD_A = thermParams.k_A;

    //        double capB = thermParams.rho_B;
    //        double LFFB = dntParams.LFFB;
    //        double rhoD_B = thermParams.k_B;

    //        double InterfaceMassFraction = chemicalComponentIndex == 0 ? 1.0 : 0.0;  //TODO!!!!

    //        // convective part
    //        // ================
    //        if (thermParams.IncludeConvection) {
    //            AddComponent(new SpeciesConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, InterfaceMassFraction, thermParams, FirstSpeciesName, SecondSpeciesName, chemicalComponentIndex));
    //        }

    //        // viscous operator (laplace)
    //        // ==========================
    //        if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
    //            double penalty = dntParams.PenaltySafety;
    //            var Visc = new Interface_MassDiffusion_LowMach(dimension, rhoD_A, rhoD_B, penalty, InterfaceMassFraction, chemicalComponentIndex, config.isMaterialAtContactLine);
    //            AddComponent(Visc);
    //        } else {
    //            throw new NotImplementedException();
    //        }

    //        AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MassFraction_n(chemicalComponentIndex));

    //        AddCoefficient("EvapMicroRegion");
    //    }

    //    public override string FirstSpeciesName => phaseA;

    //    public override string SecondSpeciesName => phaseB;

    //    public override string CodomainName => codomainName;
    //}

    public class NSEimmersedBoundary_Newton_LowMach : SurfaceEquation {
        protected string m_codomainName;
        protected string m_fluidPhase;
        protected string m_solidPhase;
        protected int m_iLevSet;

        //Methode aus der XNSF_OperatorFactory
        public NSEimmersedBoundary_Newton_LowMach(
            string fluidPhase,
            string solidPhase,
            int iLevSet,
            int d,
            int D,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
                        MaterialLaw EoS,
            bool isMovingMesh,
            PhysicsMode physicsMode) : base() //
        {
            m_fluidPhase = fluidPhase;
            m_solidPhase = solidPhase;
            m_iLevSet = iLevSet;
            m_codomainName = EquationNames.MomentumEquationComponent(d);
            AddInterfaceNSE(D, d, boundaryMap, config, EoS, isMovingMesh, physicsMode);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D).Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure));

            AddParameter(NSECommon.VariableNames.AsLevelSetVariable(NSECommon.VariableNames.LevelSetCGidx(m_iLevSet), NSECommon.VariableNames.VelocityVector(D)).ToArray());
        }

        private void AddInterfaceNSE(
            int D,
            int d,
            IncompressibleBoundaryCondMap boundaryMap,
            XNSEC_OperatorConfiguration config,
                        MaterialLaw EoS,
            bool isMovingMesh,
            PhysicsMode physicsMode) {
            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rho, LFF, mu;
            switch (this.m_fluidPhase) {
                case "A":
                    rho = physParams.rho_A;
                    LFF = dntParams.LFFA;
                    mu = physParams.mu_A;
                    break;

                case "B":
                    rho = physParams.rho_B;
                    LFF = dntParams.LFFB;
                    mu = physParams.mu_B;
                    break;

                default: throw new NotSupportedException($"Unknown fluid species: {this.m_fluidPhase}");
            }

            // convective operator
            // ===================
            if (physParams.IncludeConvection && config.isTransport) {
                var ConvIB = new BoSSS.Solution.NSECommon.Operator.Convection.ConvectionAtIB_Newton(
                           d, D, LFF, boundaryMap, rho, isMovingMesh,
                           m_iLevSet, m_fluidPhase, m_solidPhase, true);

                AddComponent(ConvIB);
            }

            // pressure gradient
            // =================
            if (config.isPressureGradient) {
                var presLs = new BoSSS.Solution.NSECommon.Operator.Pressure.PressureFormAtIB(d, D, m_iLevSet, m_fluidPhase, m_solidPhase);
                AddComponent(presLs);
            }

            // viscous operator
            // ================
            if (config.isViscous && (mu != 0.0)) {
                double penalty = dntParams.PenaltySafety;
                if (dntParams.IBM_BoundaryType == IBM_BoundaryType.NoSlip) {
                    switch (dntParams.ViscosityMode) {
                        case ViscosityMode.Standard:
                        case ViscosityMode.TransposeTermMissing:
                            //AddComponent(new BoSSS.Solution.NSECommon.Operator.Viscosity.ViscosityAtIB(d, D, penalty, mu, m_iLevSet, m_fluidPhase, m_solidPhase, true));
                            break;

                        case ViscosityMode.FullySymmetric:
                            //throw new NotImplementedException("todo");
                            var viscOption = ViscosityOption.VariableViscosityDimensionless;

                            AddComponent(new BoSSS.Solution.NSECommon.Operator.Viscosity.ViscosityAtIB_FullySymmetric_LowMach(d, D, penalty, viscOption, physicsMode, mu, config.Reynolds, EoS, m_iLevSet, m_fluidPhase, m_solidPhase, true));
                            break;

                        case ViscosityMode.Viscoelastic:
                            throw new NotImplementedException("todo");

                        default:
                            throw new NotImplementedException();
                    }
                } else {
                    throw new NotImplementedException();
                }
            }
        }

        public override string FirstSpeciesName {
            get { return m_fluidPhase; }
        }

        public override string SecondSpeciesName {
            get { return m_solidPhase; }
        }

        public override string CodomainName {
            get {
                return m_codomainName;
            }
        }
    }

    public class ImmersedBoundaryContinuity_LowMach : SurfaceEquation {
        private string codomainName;

        /// <summary>
        ///
        /// </summary>
        public ImmersedBoundaryContinuity_LowMach(string fluidPhase, string solidPhase, int iLevSet, INSE_Configuration config, int D) {
            codomainName = EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));

            m_FirstSpeciesName = fluidPhase;
            m_SecondSpeciesName = solidPhase;

            // set components
            var divPen = new BoSSS.Solution.NSECommon.Operator.Continuity.DivergenceAtIB_LowMach(D, iLevSet, FirstSpeciesName, SecondSpeciesName, true, -1);

            AddComponent(divPen);
            AddParameter(NSECommon.VariableNames.AsLevelSetVariable(NSECommon.VariableNames.LevelSetCGidx(iLevSet), NSECommon.VariableNames.VelocityVector(D)).ToArray());
        }

        private string m_FirstSpeciesName;
        private string m_SecondSpeciesName;

        public override string FirstSpeciesName {
            get {
                return m_FirstSpeciesName;
            }
        }

        public override string SecondSpeciesName {
            get {
                return m_SecondSpeciesName;
            }
        }

        public override string CodomainName => codomainName;
    }

    public class ImmersedBoundaryHeat_LowMach : SurfaceEquation {
        private string codomainName;
        private string fluidPhase, solidPhase;

        //Methode aus der XNSF_OperatorFactory
        public ImmersedBoundaryHeat_LowMach(
            string phaseA,
            string phaseB,
            int iLevSet,
            int dimension,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS) : base() {
            this.fluidPhase = phaseA;
            this.solidPhase = phaseB;

            codomainName = EquationNames.HeatEquation;
            AddInterfaceHeatEq(iLevSet, dimension, config, EoS);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => fluidPhase;

        public override string SecondSpeciesName => solidPhase;

        public override string CodomainName => codomainName;

        //Methode aus der XNSF_OperatorFactory
        private void AddInterfaceHeatEq(
            int iLevSet,
            int dimension,
            XNSEC_OperatorConfiguration config,
               MaterialLaw EoS) {
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double kF;
            switch (fluidPhase) {
                case "A": { kF = thermParams.k_A; break; }
                case "B": { kF = thermParams.k_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            double kS = thermParams.k_C;
            double Tboundary = thermParams.T_sat;

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                double penalty = dntParams.PenaltySafety;

                var Visc = new ConductivityAt_ImmersedBoundary_LowMach(dimension, EoS, config.Reynolds, config.Prandtl, penalty * 1.0, Tboundary, FirstSpeciesName, SecondSpeciesName, iLevSet: iLevSet);
                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }
        }
    }

    public class ImmersedBoundaryMF_LowMach : SurfaceEquation {
        private string codomainName;
        private string fluidPhase, solidPhase;

        public ImmersedBoundaryMF_LowMach(
            string phaseA,
            string phaseB,
            int iLevSet,
            int dimension,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS,
            int mfComp,
            int noOfSpecies) : base() {
            this.fluidPhase = phaseA;
            this.solidPhase = phaseB;

            codomainName = EquationNames.SpeciesMassBalanceName(mfComp);
            AddInterfaceMFEq(iLevSet, dimension, config, EoS, mfComp, noOfSpecies);
            AddVariableNames(NSECommon.VariableNames.Temperature);
            AddVariableNames(NSECommon.VariableNames.MassFractions(noOfSpecies));

            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => fluidPhase;

        public override string SecondSpeciesName => solidPhase;

        public override string CodomainName => codomainName;

        //Methode aus der XNSF_OperatorFactory
        private void AddInterfaceMFEq(
            int iLevSet,
            int dimension,
            XNSEC_OperatorConfiguration config,
               MaterialLaw EoS,
               int mfComp,
              int noOfSpecies) {
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double kF;
            switch (fluidPhase) {
                case "A": { kF = thermParams.k_A; break; }
                case "B": { kF = thermParams.k_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            double kS = thermParams.k_C;
            double Tboundary = thermParams.T_sat;

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                double penalty = dntParams.PenaltySafety;

                double[] Le = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };
                var Visc = new SpeciesMassDiffusivity_AtLevelSet_material_LowMach(dimension, EoS, config.Reynolds, config.Prandtl, Le, mfComp, noOfSpecies, penalty, 1.0, FirstSpeciesName, SecondSpeciesName, iLevSet: iLevSet);
                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }
        }
    }

    public class ImmersedBoundaryMixtureFraction_LowMach : SurfaceEquation {
        private string codomainName;
        private string fluidPhase, solidPhase;

        //Methode aus der XNSF_OperatorFactory
        public ImmersedBoundaryMixtureFraction_LowMach(
            string phaseA,
            string phaseB,
            int iLevSet,
            int dimension,
            XNSEC_OperatorConfiguration config,
            MaterialLaw EoS) : base() {
            this.fluidPhase = phaseA;
            this.solidPhase = phaseB;

            codomainName = EquationNames.MixtureFractionEquation;
            AddMixtureFractionImmersedBoundaryEq(iLevSet, dimension, config, EoS);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.MixtureFraction);
            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => fluidPhase;

        public override string SecondSpeciesName => solidPhase;

        public override string CodomainName => codomainName;

        //Methode aus der XNSF_OperatorFactory
        private void AddMixtureFractionImmersedBoundaryEq(
            int iLevSet,
            int dimension,
            XNSEC_OperatorConfiguration config,
               MaterialLaw EoS) {
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double kF;
            switch (fluidPhase) {
                case "A": { kF = thermParams.k_A; break; }
                case "B": { kF = thermParams.k_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                double penalty = dntParams.PenaltySafety;

                var Visc = new MixtureFractionDiffusivity_AtLevelSet_material_LowMach(dimension, EoS, config.Reynolds, config.Prandtl, penalty, FirstSpeciesName, SecondSpeciesName, iLevSet: iLevSet);
                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }
        }
    }
}