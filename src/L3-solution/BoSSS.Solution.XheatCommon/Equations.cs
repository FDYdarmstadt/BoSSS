using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Solution.XheatCommon {
    public class Heat : BulkEquation {
        string speciesName;

        string codomainName;

        int D;

        double cap;

        public Heat(
            string spcName,
            int D,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IHeat_Configuration config) {

            speciesName = spcName;
            codomainName = EquationNames.HeatEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            if (config.getConductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.HeatFluxVector(D));

            this.D = D;

            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double capSpc, LFFSpc;
            switch (spcName) {
                case "A": { capSpc = thermParams.rho_A * thermParams.c_A; LFFSpc = dntParams.LFFA; break; }
                case "B": { capSpc = thermParams.rho_B * thermParams.c_B; LFFSpc = dntParams.LFFB; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            cap = capSpc;

            // convective part
            // ================
            if (thermParams.IncludeConvection) {
                DefineConvective(D, boundaryMap, spcName, capSpc, LFFSpc, config.useUpwind);               
            }


            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                double penalty = dntParams.PenaltySafety;

                var Visc = new ConductivityInSpeciesBulk(
                    penalty, //dntParams.UseGhostPenalties ? 0.0 : penalty, 
                    1.0,
                    boundaryMap, D, spcName, thermParams.k_A, thermParams.k_B);

                AddComponent(Visc);

                //if (dntParams.UseGhostPenalties) {
                //    var ViscPenalty = new ConductivityInSpeciesBulk(penalty * 1.0, 0.0, boundaryMap, D,
                //        spcName, spcId, thermParams.k_A, thermParams.k_B);
                //    AddGhostComponent(ViscPenalty);
                //}
            } else {
                // Local DG add divergence term
                AddComponent(new HeatFluxDivergenceInSpeciesBulk(D, boundaryMap, spcName));
            }
        }

        protected virtual void DefineConvective(int D, ThermalMultiphaseBoundaryCondMap boundaryMap, string spcName, double capSpc, double LFFSpc, bool useUpwind) {
            IEquationComponent conv;            
            if (useUpwind) {
                conv = new HeatConvectionInSpeciesBulk_Upwind(D, boundaryMap, spcName, capSpc);
            } else {
                conv = new HeatConvectionInSpeciesBulk_LLF(D, boundaryMap, spcName, capSpc, LFFSpc);
            }
            AddComponent(conv);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => cap;

        public override string CodomainName => codomainName;
    }

    /// <summary>
    /// Same as <see cref="Heat"/>, but with <see cref="HeatConvectionInSpeciesBulk_LLF_Newton"/> to work with Newton Solver
    /// </summary>
    public class Heat_Newton : Heat {
        
        public Heat_Newton(
            string spcName,
            int D,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IHeat_Configuration config) : base(spcName, D, boundaryMap, config){
           
        }

        protected override void DefineConvective(int D, ThermalMultiphaseBoundaryCondMap boundaryMap, string spcName, double capSpc, double LFFSpc, bool useUpwind) {
            IEquationComponent conv;
            if (useUpwind) {
                throw new NotImplementedException();
            } else {
                //conv = new HeatConvectionInSpeciesBulk_LLF_Newton(D, boundaryMap, spcName, spcId, capSpc, LFFSpc, LsTrk);
                conv = new HeatConvectionInSpeciesBulk_Hamiltonian_Newton(D, boundaryMap, spcName, capSpc, LFFSpc);
            }
            AddComponent(conv);
        }
    }

    public class HeatFlux : BulkEquation {
        string speciesName;

        string codomainName;

        int D;

        double cap = 0.0;

        public HeatFlux(
            string spcName,
            int d,
            int D,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IHeat_Configuration config) {

            speciesName = spcName;
            codomainName = EquationNames.AuxHeatFluxComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.HeatFluxVector(D)[d]);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            this.D = D;

            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double kSpc;
            switch (spcName) {
                case "A": { kSpc = thermParams.k_A; break; }
                case "B": { kSpc = thermParams.k_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // viscous operator (laplace)
            // ==========================
            AddComponent(new AuxiliaryHeatFlux_Identity(d, spcName));   // cell local
            AddComponent(new TemperatureGradientInSpeciesBulk(D, d, boundaryMap, spcName, kSpc));
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => cap;

        public override string CodomainName => codomainName;
    }

    public class HeatInterface : SurfaceEquation {


        string codomainName;
        string phaseA, phaseB;
        //Methode aus der XNSF_OperatorFactory
        public HeatInterface(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IXHeat_Configuration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.HeatEquation;
            AddInterfaceHeatEq(dimension, boundaryMap, config);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            if (config.getConductMode != ConductivityInSpeciesBulk.ConductivityMode.SIP) AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.HeatFluxVector(dimension));
            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;


        //Methode aus der XNSF_OperatorFactory
        void AddInterfaceHeatEq(
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
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
            if (thermParams.IncludeConvection) {
                Console.WriteLine("include heat convection");
                DefineConvective(dimension, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh);                
            }

            // viscous operator (laplace)
            // ==========================
            if (config.getConductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) {

                double penalty = dntParams.PenaltySafety;

                //var Visc = new ConductivityAtLevelSet(LsTrk, kA, kB, penalty * 1.0, Tsat);
                var Visc = new ConductivityAtLevelSet_material(dimension, kA, kB, penalty * 1.0, Tsat);
                AddComponent(Visc);
            } else {
                AddComponent(new HeatFluxDivergencetAtLevelSet(dimension));
            }

        }

        protected virtual void DefineConvective(int dimension, double capA, double capB, double LFFA, double LFFB, ThermalMultiphaseBoundaryCondMap boundaryMap, bool isMovingMesh) {
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(dimension));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(dimension));
            //AddComponent(new HeatConvectionAtLevelSet_LLF(dimension, LsTrk, capA, capB, LFFA, LFFB, boundaryMap, config.isMovingMesh, Tsat));
            AddComponent(new HeatConvectionAtLevelSet_LLF_material(dimension, capA, capB, LFFA, LFFB, boundaryMap, isMovingMesh));
        }
    }

    /// <summary>
    /// Same as <see cref="HeatInterface"/>, but with <see cref="HeatConvectionAtLevelSet_LLF_material_Newton"/> to work with Newton Solver
    /// </summary>
    public class HeatInterface_Newton : HeatInterface {

        public HeatInterface_Newton(
            string phaseA,
            string phaseB,
            int dimension,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IXHeat_Configuration config) : base(phaseA, phaseB, dimension, boundaryMap, config) {           
        }        

        protected override void DefineConvective(int dimension, double capA, double capB, double LFFA, double LFFB, ThermalMultiphaseBoundaryCondMap boundaryMap, bool isMovingMesh) {
            AddComponent(new HeatConvectionAtLevelSet_LLF_material_Newton_Hamiltonian(dimension, capA, capB, LFFA, LFFB, boundaryMap, isMovingMesh, FirstSpeciesName, SecondSpeciesName));
        }
    }

    public class HeatFluxInterface : SurfaceEquation {


        string codomainName;
        string phaseA, phaseB;
        //Methode aus der XNSF_OperatorFactory
        public HeatFluxInterface(
            string phaseA,
            string phaseB,
            int dimension,
            int d,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IXHeat_Configuration config) : base() {

            this.phaseA = phaseA;
            this.phaseB = phaseB;

            codomainName = EquationNames.AuxHeatFlux(dimension)[d];
            AddInterfaceHeatEq(dimension, d, boundaryMap, config);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Temperature);
            AddCoefficient("EvapMicroRegion");
        }

        public override string FirstSpeciesName => phaseA;

        public override string SecondSpeciesName => phaseB;

        public override string CodomainName => codomainName;


        //Methode aus der XNSF_OperatorFactory
        void AddInterfaceHeatEq(
            int dimension,
            int d,
            ThermalMultiphaseBoundaryCondMap boundaryMap,
            IXHeat_Configuration config) {

            ThermalParameters thermParams = config.getThermParams;

            // set species arguments
            double kA = thermParams.k_A;
            double kB = thermParams.k_B;

            double Tsat = thermParams.T_sat;

            AddComponent(new TemperatureGradientAtLevelSet(d, kA, kB, Tsat));

        }
    }
}
