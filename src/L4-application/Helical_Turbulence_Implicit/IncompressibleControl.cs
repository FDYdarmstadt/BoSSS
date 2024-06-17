using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE {
    
    /// <summary>
    /// PDE-solver-control object which defines configuration options for nonlinear and linear equation solvers
    /// </summary>
    [Serializable]
    [DataContract]
    public class IncompressibleControl : HelicalControl {

        /// <summary>
        /// Constructor
        /// </summary>
        public IncompressibleControl() {
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
        }


        //
        // add configuration options specific to your solver here
        //


        /// <summary>
        /// fluid density, unit is [mass / volume ]
        /// </summary>
        [DataMember]
        public double Density = 1.0;

        /// <summary>
        /// dynamic viscosity, unit is [mass / (length*time)]
        /// </summary>
        [DataMember]
        public double Viscosity = 1.0;

        /// <summary>
        /// Penalty factor, must be positive, usually in the order of 1 to 10; 
        /// higher values improve stability but also increase numerical dissipation.
        /// </summary>
        [DataMember]
        public double PenaltySafety = 2.0;


        /// <summary>
        /// required for BoSSSpad workflow management
        /// </summary>
        public override Type GetSolverType() {
            return typeof(Helical_Turbulence_Implicit_Main);
        }

        /// <summary>
        /// Setting DG degree for all fields at once
        /// </summary>
        public override void SetDGdegree(int k) {
            if (k < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");

            base.FieldOptions.Clear();
            this.AddFieldOption("Velocity*", k);
            this.AddFieldOption("Pressure", k - 1);
        }

        /// <summary>
        /// checking values
        /// </summary>
        public override void Verify() {
            base.Verify();

            if(PenaltySafety <= 0)
                throw new ArgumentException("penalty must be positive");
            if(Density <= 0)
                throw new ArgumentException("density must be positive");
            if(Viscosity <= 0)
                throw new ArgumentException("viscosity must be positive");
        }
    }
}
