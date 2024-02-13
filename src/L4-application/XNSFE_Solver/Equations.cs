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

namespace BoSSS.Application.XNSFE_Solver {

    // Outdated forms from old XNSFE (pre 2021), disabled 4/2021 - MR
    /*
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
                    AddComponent(new MassFluxAtInterface(d, D, thermParams, sigma, config.isMovingMesh));
                    AddComponent(new ConvectionAtLevelSet_nonMaterialLLF(d, D, lsTrk, thermParams, sigma));
                    AddComponent(new ConvectionAtLevelSet_Consistency(d, D, lsTrk, -1, false, thermParams, sigma));
                }
            } else {
                AddComponent(new MassFluxAtInterface(d, D, thermParams, sigma, config.isMovingMesh));
            }

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_withEvap(lsTrk.GridDat.SpatialDimension, physParams.mu_A, physParams.mu_B, dntParams.PenaltySafety, d, thermParams, sigma));
            }

        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }
    */

    public class InterfaceNSE_Evaporation : SurfaceEquation {

        string codomainName;
        string phaseA, phaseB;
        public InterfaceNSE_Evaporation(string phaseA,
            string phaseB,
            int dimension,
            int d,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.MomentumEquationComponent(d);
            AddInterfaceNSE_Evaporation(dimension, d, config);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d));
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        private void AddInterfaceNSE_Evaporation(int D, int d, XNSFE_OperatorConfiguration config) {

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // from XNSFE_OperatorComponents            
            if (config.isTransport) {
               DefineConvective(d, D, config);
            } else {
                //  ... and when the convective terms are turned off we still need the contribution below
                if (config.isRecoilPressure) {
                    AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling(d, D, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    if (physParams.slipI != 0) {
                        AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling_Tangential(d, D, config.physParams, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    }
                }
            }           

            if (config.isViscous) {
                AddComponent(new ViscosityAtLevelSet_FullySymmetric_Evaporation_StrongCoupling(dntParams.PenaltySafety, d, D, config.getThermParams, physParams, FirstSpeciesName, SecondSpeciesName));                
            }            
        }

        protected virtual void DefineConvective(int d, int D, XNSFE_OperatorConfiguration config) {
            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
            if (!config.isMovingMesh) {
                // the following terms decode the condition at the interface (consider the similarity to the rankine hugoniot condition)
                // for the moving mesh discretization this condition is already contained in the convective terms
                // therefore we only need these terms when using splitting...                
                AddComponent(new MassFluxAtLevelSet_withMassFlux(d, D, physParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                AddComponent(new ConvectionAtLevelSet_nonMaterialLLF_withMassFlux(d, D, physParams, FirstSpeciesName, SecondSpeciesName));
                AddComponent(new ConvectionAtLevelSet_Consistency_withMassFlux(d, D, -1, false, physParams, FirstSpeciesName, SecondSpeciesName));
                
            } else {
                AddComponent(new ConvectionAtLevelSet_MovingMesh_withMassFlux(d, D, physParams, FirstSpeciesName, SecondSpeciesName));
            }
        }        
    

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

   public class InterfaceNSE_Evaporation_Newton : InterfaceNSE_Evaporation {

        public InterfaceNSE_Evaporation_Newton(string phaseA,
            string phaseB,
            int dimension,
            int d,
            XNSFE_OperatorConfiguration config) : base(phaseA, phaseB, dimension, d, config) {
        }

        protected override void DefineConvective(int d, int D, XNSFE_OperatorConfiguration config) {
            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;
            if (!config.isMovingMesh) {
                // the following terms decode the condition at the interface (consider the similarity to the rankine hugoniot condition)
                // for the moving mesh discretization this condition is already contained in the convective terms
                // therefore we only need these terms when using splitting...
                if (config.isRecoilPressure) {
                    AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling(d, D, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    if (physParams.slipI != 0) {
                        AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling_Tangential(d, D, config.physParams, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    }
                }
                AddComponent(new ConvectionAtLevelSet_nonMaterialLLF_Evaporation_StrongCoupling_Newton(d, D, config.getPhysParams, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
                AddComponent(new ConvectionAtLevelSet_Consistency_Evaporation_StrongCoupling_Newton(d, D, -1, false, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
            } else {
                if (config.isRecoilPressure) {
                    AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling(d, D, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    if(physParams.slipI != 0) {
                        AddComponent(new MassFluxAtLevelSet_Evaporation_StrongCoupling_Tangential(d, D, config.physParams, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
                    }
                }
            }
        }
    }

    

    public class InterfaceContinuity_Evaporation : SurfaceEquation {

        string codomainName;
        string phaseA, phaseB;
        public InterfaceContinuity_Evaporation(string phaseA,
            string phaseB,
            int dimension,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.ContinuityEquation;
            AddInterfaceContinuity_Evaporation(dimension, config);
            if (config.prescribedMassflux != null)
                AddCoefficient("PrescribedMassFlux");

        }

        protected virtual void AddInterfaceContinuity_Evaporation(int D, XNSFE_OperatorConfiguration config) {
            PhysicalParameters physicalParameters = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension);
            AddComponent(new DivergenceAtLevelSet_withMassFlux(D, -1, false, physicalParameters, FirstSpeciesName, SecondSpeciesName));                     
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// same as <see cref="InterfaceContinuity_Evaporation"/> but using Newton solver compatible components
    /// </summary>
    public class InterfaceContinuity_Evaporation_Newton : InterfaceContinuity_Evaporation {

        public InterfaceContinuity_Evaporation_Newton(string phaseA,
            string phaseB,
            int dimension,
            XNSFE_OperatorConfiguration config) : base(phaseA, phaseB, dimension, config) {
        }

        protected override void AddInterfaceContinuity_Evaporation(int D, XNSFE_OperatorConfiguration config) {
            AddComponent(new DivergenceAtLevelSet_Evaporation_StrongCoupling(D, -1, false, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
        }
    }

    public class HeatInterface_Evaporation : SurfaceEquation {


        string codomainName;
        string phaseA, phaseB;
        public HeatInterface_Evaporation(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            XNSFE_OperatorConfiguration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.HeatEquation;
            AddInterfaceHeatEq(dimension, boundaryMap, config);
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
            XNSFE_OperatorConfiguration config) {

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
                DefineConvective(dimension, Tsat, config);                
            }

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {

                double penalty = dntParams.PenaltySafety;

                var Visc = new ConductivityAtLevelSet_withMassflux(dimension, kA, kB, penalty * 1.0, Tsat, config.isMaterialAtContactLine);
                AddComponent(Visc);
            } else {
                throw new NotImplementedException();
            }

        }

        protected virtual void DefineConvective(int dimension, double Tsat, XNSFE_OperatorConfiguration config) {
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(dimension));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(dimension));
            if (config.isMovingMesh) {
                AddComponent(new HeatConvectionAtLevelSet_MovingMesh_withMassflux(dimension, Tsat, config.getPhysParams, config.getThermParams, FirstSpeciesName, SecondSpeciesName));
            } else {
                throw new NotImplementedException("Evaporation only implemented with use of Newton-solver!");
                //AddComponent(new HeatConvectionAtLevelSet_LLF_withMassflux(dimension, LsTrk, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh, Tsat, physParams));
                //AddComponent(new HeatConvectionAtLevelSet_Direct(dimension, LsTrk, capA, capB, Tsat, physParams, LsTrk.GetSpeciesId(phaseA)));
                //AddComponent(new HeatConvectionAtLevelSet_Direct(dimension, LsTrk, capA, capB, Tsat, physParams, LsTrk.GetSpeciesId(phaseB)));                
            }
        }
    }

    /// <summary>
    /// same as <see cref="HeatInterface_Evaporation"/> but using Newton solver compatible components
    /// </summary>
    public class HeatInterface_Evaporation_Newton : HeatInterface_Evaporation {
        public HeatInterface_Evaporation_Newton(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            XNSFE_OperatorConfiguration config) : base(phaseA, phaseB, dimension, boundaryMap, config) {
        }

        protected override void DefineConvective(int dimension, double Tsat, XNSFE_OperatorConfiguration config) {
            if (config.isMovingMesh) {
                AddComponent(new HeatFluxAtLevelSet_Evaporation_StrongCoupling(dimension, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
            } else {
                ThermalParameters thermParams = config.getThermParams;
                DoNotTouchParameters dntParams = config.getDntParams;

                // set species arguments
                double capA = thermParams.rho_A * thermParams.c_A;
                double LFFA = dntParams.LFFA;
                double kA = thermParams.k_A;

                double capB = thermParams.rho_B * thermParams.c_B;
                double LFFB = dntParams.LFFB;
                double kB = thermParams.k_B;
               // AddComponent(new HeatConvectionAtLevelSet_LLF_withMassflux_StrongCoupling(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, Tsat, thermParams, FirstSpeciesName, SecondSpeciesName));
                AddComponent(new HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, Tsat, thermParams, FirstSpeciesName, SecondSpeciesName));
            }
        }
    }




    /// <summary>
    /// same as <see cref="HeatInterface_Evaporation"/> but using Newton solver compatible components
    /// </summary>
    public class HeatInterface_Evaporation_Newton_LowMach : HeatInterface_Evaporation {
        public HeatInterface_Evaporation_Newton_LowMach(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            XNSFE_OperatorConfiguration config) : base(phaseA, phaseB, dimension, boundaryMap, config) {
        }

        protected override void DefineConvective(int dimension, double Tsat, XNSFE_OperatorConfiguration config) {
            if (config.isMovingMesh) {
                AddComponent(new HeatFluxAtLevelSet_Evaporation_StrongCoupling(dimension, config.getThermParams, config.isMovingMesh, FirstSpeciesName, SecondSpeciesName));
            } else {
                ThermalParameters thermParams = config.getThermParams;
                DoNotTouchParameters dntParams = config.getDntParams;

                // set species arguments
                double capA = thermParams.rho_A * thermParams.c_A;
                double LFFA = dntParams.LFFA;
                double kA = thermParams.k_A;

                double capB = thermParams.rho_B * thermParams.c_B;
                double LFFB = dntParams.LFFB;
                double kB = thermParams.k_B;
                // AddComponent(new HeatConvectionAtLevelSet_LLF_withMassflux_StrongCoupling(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, Tsat, thermParams, FirstSpeciesName, SecondSpeciesName));
                AddComponent(new HeatConvectionAtLevelSet_LLF_Evaporation_StrongCoupling_Hamiltonian_LowMach(dimension, capA, capB, LFFA, LFFB, config.isMovingMesh, Tsat, thermParams, FirstSpeciesName, SecondSpeciesName));
            }
        }
    }


    class NavierStokesBuoyancy : SpatialEquation {

        string codomainName;
        string speciesName;
        public NavierStokesBuoyancy(string spc,
            int dimension,
            int d,
            XNSFE_OperatorConfiguration config) : base() {

            this.speciesName = spc;

            codomainName = EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);

            var thermparams = config.getThermParams;

            string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(dimension)[d];
            string gravityOfSpecies = gravity + "#" + speciesName;
            // keep in mind the special treatment of gravity in XNSE G = -rho * g
            var buoyancyComponent = new Solution.XheatCommon.BoussinesqApproximation_Buoyancy(speciesName, gravityOfSpecies, d, thermparams);
            AddComponent(buoyancyComponent);
            AddParameter(gravityOfSpecies);
        }        

        public override string CodomainName => codomainName;
    }
}
