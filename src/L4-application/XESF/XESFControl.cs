
using BoSSS.Foundation.IO;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
//using XDGShock.Fluxes;
using XESF.Fluxes;
using ApplicationWithIDT;
using BoSSS.Foundation;
using System.Linq;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Residual;

namespace XESF {
    public class XESFControl : IDTControl, ICompressibleControl {
        public XESFControl():base() {
            base.NoOfMultigridLevels = 1;
            base.quadOrderFunc = (int[] A, int[] B, int[] C) =>  Math.Abs(2*A.Max()) + Math.Abs(C.Max()) + Math.Max(this.LevelSetDegree,this.LevelSetTwoDegree);
            base.CutCellQuadratureType = BoSSS.Foundation.XDG.CutCellQuadratureMethod.Saye;
        }

        [DataMember]
        public CompressibleConfiguration CompressibleConfiguration { get; private set; } = new CompressibleConfiguration();

        ICompressibleConfiguration ICompressibleControl.CompressibleConfiguration => CompressibleConfiguration;

        //[NotNull]
        public IEquationOfState EquationOfState {
            get => CompressibleConfiguration.EquationOfState;
            set => CompressibleConfiguration.EquationOfState = value;
        }

        public IViscosityLaw ViscosityLaw {
            get => CompressibleConfiguration.ViscosityLaw;
            set => CompressibleConfiguration.ViscosityLaw = value;
        }

        //[ExclusiveLowerBound(0.0)]
        public double MachNumber {
            get => CompressibleConfiguration.MachNumber;
            set => CompressibleConfiguration.MachNumber = value;
        }

        public double ReynoldsNumber {
            get => CompressibleConfiguration.ReynoldsNumber;
            set => CompressibleConfiguration.ReynoldsNumber = value;
        }

        public double PrandtlNumber {
            get => CompressibleConfiguration.PrandtlNumber;
            set => CompressibleConfiguration.PrandtlNumber = value;
        }

        //[InclusiveLowerBound(0.0)]
        public double FroudeNumber {
            get => CompressibleConfiguration.FroudeNumber;
            set => CompressibleConfiguration.FroudeNumber = value;
        }

        //[InclusiveLowerBound(0.0)]
        public double ViscosityRatio {
            get => CompressibleConfiguration.ViscosityRatio;
            set => CompressibleConfiguration.ViscosityRatio = value;
        }

        //[InclusiveLowerBound(0.0)]
        public int PrintInterval {
            get => CompressibleConfiguration.PrintInterval;
            set => CompressibleConfiguration.PrintInterval = value;
        }

        //[InclusiveLowerBound(0)]
        public int ResidualInterval {
            get => CompressibleConfiguration.ResidualInterval;
            set => CompressibleConfiguration.ResidualInterval = value;
        }

        public ResidualLoggerTypes ResidualLoggerType {
            get => CompressibleConfiguration.ResidualLoggerType;
            set => CompressibleConfiguration.ResidualLoggerType = value;
        }

        public IDictionary<string, double> ResidualBasedTerminationCriteria {
            get => CompressibleConfiguration.ResidualBasedTerminationCriteria;
            set => CompressibleConfiguration.ResidualBasedTerminationCriteria = value;
        }

        public IReadOnlyDictionary<Variable, int> VariableToDegreeMap => CompressibleConfiguration.VariableToDegreeMap;

        public int DensityDegree => CompressibleConfiguration.DensityDegree;

        public int MomentumDegree => CompressibleConfiguration.MomentumDegree;

        public int EnergyDegree => CompressibleConfiguration.EnergyDegree;

        public Material GetMaterial() {
            return CompressibleConfiguration.GetMaterial();
        }

        public void AddVariable(Variable variable, int degree, bool saveToDB = true) {
            CompressibleConfiguration.SetVariableDegree(variable, degree);

            FieldOpts.SaveToDBOpt option = saveToDB ? FieldOpts.SaveToDBOpt.TRUE : FieldOpts.SaveToDBOpt.FALSE;
            var fieldOpts = new FieldOpts() {
                Degree = degree,
                SaveToDB = option
            };

            if (FieldOptions.ContainsKey(variable)) {
                FieldOptions[variable] = fieldOpts;
            } else {
                FieldOptions.Add(variable, fieldOpts);
            }
        }

        public override Type GetSolverType() {
            return typeof(XESFMain);
        }
      
        public string PointPath { get; set; }
        
        public ConvectiveBulkFluxes ConvectiveBulkFlux { get; set; } = ConvectiveBulkFluxes.OptimizedHLLC;

        public FluxVersion FluxVersion { get; set; } = FluxVersion.Optimized;

        public ConvectiveInterfaceFluxes ConvectiveInterfaceFlux_LsOne { get; set; } = ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var;

        public ConvectiveInterfaceFluxes ConvectiveInterfaceFlux_LsTwo { get; set; } = ConvectiveInterfaceFluxes.OptimizedHLLCInterface;
        public int IVTimestepNumber { get; set; } = 0;
        public int StartDegree { get; set; } = 0;
        public double ExactEnthalpy { get; internal set; }


        //public SensorTypes SensorType { get; internal set; };
        //public string SensorVariable { get; internal set; } = null;
        //public double SensorLimit { get; internal set; } = double.MinValue;
        //public ArtificialViscosityLawTypes ArtificialViscosityLawType { get; internal set; }
        //public DiffusiveBulkFluxes DiffusiveBulkFlux { get; internal set; }
        //public DiffusiveInterfaceFluxes DiffusiveInterfaceFlux { get; internal set; }

    }

}
