using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Low-Mach Mass fraction in the bulk phase
    /// </summary>
    public class LowMachUnsteadyEquationPart : BulkEquation {
        private string speciesName;
        private string codomainName;
        private double m_massScale;

        public LowMachUnsteadyEquationPart(
            string spcName,
            int D,
            string varname,
            string _codomainName,
            int NoOfChemicalSpecies,
            MaterialLaw EoS,
            double massScale = 1.0,
            double heatCapacityRatio = 1.0) {
            m_massScale = massScale;
            speciesName = spcName;
            codomainName = _codomainName;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(NSECommon.VariableNames.Pressure);
            if (EoS is MaterialLawMixtureFractionNew) {
                AddVariableNames(NSECommon.VariableNames.MixtureFraction);
            } else {
                AddVariableNames(NSECommon.VariableNames.Temperature);
                AddVariableNames(NSECommon.VariableNames.MassFractions(NoOfChemicalSpecies));
            }
            var timeDerivative = new LowMachUnsteadyTerm(spcName, EoS, varname, D, NoOfChemicalSpecies);
            AddComponent(timeDerivative);
        }

        public override string SpeciesName => speciesName;
        public override double MassScale => m_massScale;
        public override string CodomainName => codomainName;
    }



    /// <summary>
    /// Low-Mach Mass fraction in the bulk phase
    /// </summary>
    public class LowMachUnsteadyEquationPart_MF : BulkEquation {
        private string speciesName;
        private string codomainName;
        private double m_massScale;

        public LowMachUnsteadyEquationPart_MF(
            string spcName,
            int D,
            string varname,
            string _codomainName,
            int NoOfChemicalSpecies,
            MaterialLaw EoS,
            double massScale = 1.0,
            double heatCapacityRatio = 1.0) {
            m_massScale = massScale;
            speciesName = spcName;
            codomainName = _codomainName;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(NSECommon.VariableNames.Pressure);
            AddVariableNames(NSECommon.VariableNames.MixtureFraction);
          
            var timeDerivative = new LowMachUnsteadyTerm(spcName, EoS, varname, D);
            AddComponent(timeDerivative);
        }

        public override string SpeciesName => speciesName;
        public override double MassScale => m_massScale;
        public override string CodomainName => codomainName;
    }



    public class LowMachUnsteadyTerm : MassMatrixLowMachComponent, ISpeciesFilter {
        // Constructor for Combustionsimulation
        public LowMachUnsteadyTerm(string spcName, MaterialLaw EoS, string varname, int spatDim, int NumberOfReactants) : base(EoS, varname, spatDim, NumberOfReactants) {
            ValidSpecies = spcName;
        }
        // Constructor for MixtureFraction  simulation
        public LowMachUnsteadyTerm(string spcName, MaterialLaw EoS, string varname, int spatDim) : base(EoS, varname, spatDim) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }
}