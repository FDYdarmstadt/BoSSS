using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Residual;

namespace BoSSS.Solution.CompressibleFlowCommon {

    public interface ICompressibleConfiguration {
        IEquationOfState EquationOfState { get; set; }

        IViscosityLaw ViscosityLaw { get; set; }

        double MachNumber { get; set; }

        double ReynoldsNumber { get; set; }

        double PrandtlNumber { get; set; }

        double FroudeNumber { get; set; }

        double ViscosityRatio { get; set; }

        int PrintInterval { get; set; }

        int ResidualInterval { get; set; }

        ResidualLoggerTypes ResidualLoggerType { get; set; }

        IDictionary<string, double> ResidualBasedTerminationCriteria { get; set; }

        IReadOnlyDictionary<Variable, int> VariableToDegreeMap { get; }

        Material GetMaterial();

        void SetVariableDegree(Variable variable, int degree);

        int DensityDegree { get; }

        int MomentumDegree { get; }

        int EnergyDegree { get; }
    }

    public interface ICompressibleControl {
        ICompressibleConfiguration CompressibleConfiguration { get; }
    }
}
