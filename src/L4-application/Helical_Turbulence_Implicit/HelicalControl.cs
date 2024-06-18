using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Text.Json.Serialization;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE {
    [Flags]
    public enum TermSwitch {
        // z momentum
        Pressure_MomXI = 1,  //   0001
        Viscosity_1stOrder_MomXI = 2, // 0010
        Viscosity_2ndOrder_MomXI = 4, // 0100

        // eta momentum
        Viscosity_1stOrder_MomETA = 8, // 1000
        Viscosity_2ndOrder_MomETA = 16, // 10000

        //r momentum
        Pressure_MomR = 32,  //   100000
        Viscosity_1stOrder_MomR = 64, // 1000000
        Viscosity_2ndOrder_MomR = 128, // 10000000
        AllOn = 0x0FFFFFFF


        //NotUsed = 256, // 512 1024
    }


    public enum BoundaryTypeE {
        Dirichlet,
        Neumann
    }

    [Serializable]
    [DataContract]
    public class HelicalControl : AppControlSolver {

        /// <summary>
        /// experimental, used if <see cref="Config.UseDiagonalPmg"/> is not set.
        /// Then low order and high order blocks are both solved by direct solver.
        /// </summary>
        private ISparseSolver hiSolver;

        [DataMember]
        public int Resolution_R;
        [DataMember]
        public int Resolution_Xi;

        // -------------------------------------------- //
        /* fk, 08dec23: maybe these options are to complex for the user /end-user
        // Block-Preconditiond for the velocity/momentum-block of the saddle-point system
        [DataMember]
        //public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.Eye;
        public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
        //public MultigridOperator.Mode VelocityBlockPrecondMode = MultigridOperator.Mode.LeftInverse_DiagBlock;

        // Block-Preconditiond for the pressure-block of the saddle-point system
        [DataMember]
        public MultigridOperator.Mode PressureBlockPrecondMode = MultigridOperator.Mode.Eye;
        */
        //---------------------------------------------- //

        /// <summary>
        /// Factor for penalty
        /// </summary>
        [DataMember]
        public double penaltySafety = 4;

        [DataMember]
        public bool PressureReferencePoint = true;

        [DataMember]
        public TermSwitch TermSwitch =
                TermSwitch.Viscosity_1stOrder_MomXI
                | TermSwitch.Pressure_MomXI
                | TermSwitch.Viscosity_2ndOrder_MomXI
                | TermSwitch.Viscosity_1stOrder_MomETA
                | TermSwitch.Viscosity_2ndOrder_MomETA
                | TermSwitch.Viscosity_1stOrder_MomR
                | TermSwitch.Pressure_MomR
                | TermSwitch.Viscosity_2ndOrder_MomR;

        [DataMember]
        public bool ExactResidual = true;


        [DataMember]
        public bool HagenPoisseulle = false;

        /// <summary>
        /// solving steady-state or transient
        /// </summary>
        [JsonIgnore]
        public bool steady {
            get {
                switch(base.TimesteppingMode) {
                    case _TimesteppingMode.Steady: return true;
                    case _TimesteppingMode.Transient: return false;
                    default: throw new NotImplementedException();
                }
            }
            set {
                if(value)
                    base.TimesteppingMode = _TimesteppingMode.Steady;
                else
                    base.TimesteppingMode = _TimesteppingMode.Transient;
            }
        }

        /// <summary>
        /// - true: solve Navier-Stokes
        /// - false: solve Stokes, convective terms deactivated.
        /// </summary>
        [DataMember]
        public bool NavierStokes = true;




        /// <summary>
        /// DG-Degree
        /// </summary>
        [DataMember]
        public int dg_degree = -1;

        //public override void SetDGdegree(int degree) {

        //    base.FieldOptions.Clear();
        //    base.FieldOptions.Add("Pressure", new FieldOpts() {
        //        Degree = degree - 1,
        //        SaveToDB = FieldOpts.SaveToDBOpt.TRUE
        //    });
        //    base.FieldOptions.Add("ur", new FieldOpts() {
        //        Degree = degree,
        //        SaveToDB = FieldOpts.SaveToDBOpt.TRUE
        //    });
        //    base.FieldOptions.Add("uxi", new FieldOpts() {
        //        Degree = degree,
        //        SaveToDB = FieldOpts.SaveToDBOpt.TRUE
        //    });
        //    base.AddFieldOption("ueta", degree);

        //    base.FieldOptions.Add("PhiDG", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //    base.FieldOptions.Add("Phi", new FieldOpts() { Degree = degree, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
        //}

        [DataMember]
        public bool R0fixOn = false;

        [DataMember]
        public double rMin = 0.1;

        [DataMember]
        public double maxAmpli = 0.01;

        [DataMember]
        public double rMax = 1;

        [DataMember]
        public string grid;

        [DataMember]
        public int restartTimeStep;
    }

}
