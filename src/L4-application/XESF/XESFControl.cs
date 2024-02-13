
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

namespace XESF {
    public class XESFControl : IDTControl {
        public XESFControl():base() {
            base.quadOrderFunc = (int[] A, int[] B, int[] C) =>  Math.Abs(2*A.Max()) + Math.Abs(C.Max()) + Math.Max(this.LevelSetDegree,this.LevelSetTwoDegree);
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