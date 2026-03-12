using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
using BoSSS.Foundation;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Residual;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.CompressibleFlowCommon {

    [Serializable]
    [DataContract]
    public class CompressibleConfiguration : ICompressibleConfiguration {

        [DataMember]
        //[NotNull]
        public IEquationOfState EquationOfState { get; set; } = IdealGas.Air;

        [DataMember]
        public IViscosityLaw ViscosityLaw { get; set; } = new ConstantViscosity();

        [DataMember]
        //[ExclusiveLowerBound(0.0)]
        public double MachNumber { get; set; } = 1.0;

        [DataMember]
        public double ReynoldsNumber { get; set; } = 1.0;

        [DataMember]
        public double PrandtlNumber { get; set; } = 0.71;

        [DataMember]
        //[InclusiveLowerBound(0.0)]
        public double FroudeNumber { get; set; }

        [DataMember]
        //[InclusiveLowerBound(0.0)]
        public double ViscosityRatio { get; set; } = 0.0;

        [DataMember]
        //[InclusiveLowerBound(0.0)]
        public int PrintInterval { get; set; } = 1;

        [DataMember]
        //[InclusiveLowerBound(0)]
        public int ResidualInterval { get; set; } = 0;

        [DataMember]
        public ResidualLoggerTypes ResidualLoggerType { get; set; } = ResidualLoggerTypes.None;

        [DataMember]
        public IDictionary<string, double> ResidualBasedTerminationCriteria { get; set; } = new Dictionary<string, double>();

        [DataMember]
        protected Dictionary<Variable, int> variableFields = new Dictionary<Variable, int>();

        public IReadOnlyDictionary<Variable, int> VariableToDegreeMap => variableFields;

        public virtual Material GetMaterial() {
            return new Material(EquationOfState, ViscosityLaw, MachNumber, ReynoldsNumber, PrandtlNumber, FroudeNumber, ViscosityRatio);
        }

        public void SetVariableDegree(Variable variable, int degree) {
            variableFields[variable] = degree;
        }

        public int DensityDegree => variableFields[CompressibleVariables.Density];

        public int MomentumDegree => variableFields[CompressibleVariables.Momentum.xComponent];

        public int EnergyDegree => variableFields[CompressibleVariables.Energy];
    }
}
