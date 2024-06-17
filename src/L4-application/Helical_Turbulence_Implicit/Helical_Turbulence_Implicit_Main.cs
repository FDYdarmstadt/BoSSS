using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE {

    /// <summary>
    /// The solver for the Helically Symetric Incompressible Navier-Stokes equation.
    /// </summary>
    public class Helical_Turbulence_Implicit_Main : BoSSS.Solution.XdgTimestepping.DgApplicationWithSolver<IncompressibleControl> {

        static void Main(string[] args) {

            InitMPI();


            var c = BoSSS.Application.IncompressibleNSE.ControlExamples.ChannelFlow();
            var solver = new Helical_Turbulence_Implicit_Main();
            solver.Init(c);
            solver.RunSolverMode();
            Process.Start("mpiexec");


            _Main(args, false, delegate () {
                var p = new Helical_Turbulence_Implicit_Main();
                return p;
            });

        }

#pragma warning disable 649

        //        //#############################################################################################################################
        //        //#############################################################################################################################
        //        //#############################################################################################################################
        //        //#############################################################################################################################
        //        //#############################################################################################################################

        //        /// <summary>
        //        /// velocity
        //        /// </summary>
        //        [InstantiateFromControlFile(new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" },
        //            null,
        //            true, true,
        //            IOListOption.ControlFileDetermined)]
        //        public VectorField<SinglePhaseField> Velocity;

        //        /// <summary>
        //        /// Volume Force, dimension is acceleration, i.e. length per time-square.
        //        /// </summary>
        //        [InstantiateFromControlFile(
        //            new string[] { "Gravity_R", "Gravity_ETA", "Gravity_XI" },
        //            new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" },
        //            true, true,
        //            IOListOption.ControlFileDetermined)]
        //        public VectorField<SinglePhaseField> Gravity;

        //        /// <summary>
        //        /// Residual in the momentum equation.
        //        /// </summary>
        //        [InstantiateFromControlFile(new string[] { "ResidualMomentum_R", "ResidualMomentum_ETA", "ResidualMomentum_XI" },
        //            new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" },
        //            true, true,
        //            IOListOption.ControlFileDetermined)]
        //        public VectorField<SinglePhaseField> ResidualMomentum;

        //        /// <summary>
        //        /// Pressure
        //        /// </summary>
        //        [InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        //        SinglePhaseField Pressure;

        //        /// <summary>
        //        /// Residual of the continuity equation
        //        /// </summary>
        //        [InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        //        SinglePhaseField ResidualContinuity;
        //#pragma warning restore 649

        //        /// <summary>
        //        /// Links edge tags (<see cref="BoSSS.Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>) and
        //        /// boundary conditions in the control object (<see cref="BoSSS.Solution.Control.AppControl.BoundaryValues"/>).
        //        /// </summary>
        //        protected IncompressibleBoundaryCondMap boundaryCondMap;


        //        /// <summary>
        //        /// vector of helical velocity names
        //        /// </summary>
        //        /// <param name="D">
        //        /// spatial dimension
        //        /// </param>
        //        public static string[] VelocityVector_Helical(int D) {
        //            if (D == 2)
        //                return new string[] { "Velocity_R", "Velocity_XI" };
        //            else if (D == 3)
        //                return new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" };
        //            else
        //                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        //        }

        //        /// <summary>
        //        /// vector of helical gravity/volume force names
        //        /// </summary>
        //        /// <param name="D">
        //        /// spatial dimension
        //        /// </param>
        //        public static string[] GravityVector_Helical(int D) {
        //            if (D == 2)
        //                return new string[] { "Gravity_R", "Gravity_XI" };
        //            else if (D == 3)
        //                return new string[] { "Gravity_R", "Gravity_ETA", "Gravity_XI" };
        //            else
        //                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        //        }

        //        /// <summary>
        //        /// Declaration of the spatial HELICAL operator
        //        /// </summary>
        //        protected override DifferentialOperator GetOperatorInstance(int D) {

        //            // instantiate boundary condition mapping
        //            // ======================================
        //            boundaryCondMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

        //            // instantiate operator
        //            // ====================
        //            string[] CodName = (new[] { "ResidualMomentum_R", "ResidualMomentum_ETA", "ResidualMomentum_XI" }).GetSubVector(0, D).Cat("ResidualConti");

        //            var op = new DifferentialOperator(
        //                __DomainVar: VelocityVector_Helical(D).Cat(VariableNames.Pressure),
        //                __ParameterVar: GravityVector_Helical(D),
        //                __CoDomainVar: CodName,
        //                QuadOrderFunc: QuadOrderFunc.NonLinear(2));

        //            op.LinearizationHint = LinearizationHint.GetJacobiOperator;

        //            // Temporal Operator
        //            // =================

        //            var TempOp = new ConstantTemporalOperator(op, 0.0); // init with entire diagonal set to 0.0
        //            op.TemporalOperator = TempOp;

        //            for (int d = 0; d < D; d++)
        //                TempOp.SetDiagonal(CodName[d], Control.Density); // set momentum equation entries to density

        //            // Pressure Reference
        //            // ==================

        //            // if there is no Dirichlet boundary condition,
        //            // the mean value of the pressure is free:
        //            op.FreeMeanValue[VariableNames.Pressure] = !boundaryCondMap.DirichletPressureBoundary;

        //            // Momentum Equation
        //            // =================



        //            //diffOp_explicit.EquationComponents["rmom"].Add(new convectiveRmom());
        //            //diffOp_explicit.EquationComponents["etamom"].Add(new convectiveETAmom());
        //            //diffOp_explicit.EquationComponents["zmom"].Add(new convectiveXImom());

        //            // convective part:
        //            {
        //                for (int d = 0; d < D; d++) {

        //                    var comps = op.EquationComponents[CodName[d]];

        //                    var ConvBulk = new LocalLaxFriedrichsConvection(D, boundaryCondMap, d, Control.Density, null);
        //                    comps.Add(ConvBulk); // bulk component
        //                }
        //            }


        //            // pressure part:
        //            {
        //                for (int d = 0; d < D; d++) {
        //                    var comps = op.EquationComponents[CodName[d]];
        //                    var pres = new PressureGradientLin_d(d, boundaryCondMap);
        //                    comps.Add(pres); // bulk component

        //                }
        //            }

        //            // viscous part:
        //            {
        //                for (int d = 0; d < D; d++) {
        //                    var comps = op.EquationComponents[CodName[d]];

        //                    double penalty_bulk = this.Control.PenaltySafety;

        //                    var Visc = new SipViscosity_GradU(penalty_bulk, d, D, boundaryCondMap,
        //                        ViscosityOption.ConstantViscosity,
        //                        constantViscosityValue: Control.Viscosity);
        //                    comps.Add(Visc); // bulk component GradUTerm
        //                }
        //            }


        //            // Continuity equation
        //            // ===================
        //            {
        //                for (int d = 0; d < D; d++) {
        //                    var src = new Divergence_DerivativeSource(d, D);
        //                    var flx = new Divergence_DerivativeSource_Flux(d, boundaryCondMap);
        //                    op.EquationComponents[CodName[D]].Add(src);
        //                    op.EquationComponents[CodName[D]].Add(flx);
        //                }


        //                //IBM_Op.EquationComponents["div"].Add(new PressureStabilization(1, 1.0 / this.Control.PhysicalParameters.mu_A));
        //            }

        //            // Gravity parameter
        //            // =================

        //            op.ParameterFactories.Add(delegate (IReadOnlyDictionary<string, DGField> DomainVarFields) {
        //                return D.ForLoop(d => (VariableNames.Gravity_d(d), this.Gravity[d] as DGField));
        //            });


        //            // commit & return
        //            // ===============
        //            op.Commit();
        //            return op;
        //        }

        //#############################################################################################################################
        //#############################################################################################################################
        //#############################################################################################################################
        //#############################################################################################################################
        //#############################################################################################################################

        /// <summary>
        /// velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" },
            null,
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Velocity;

        /// <summary>
        /// Volume Force, dimension is acceleration, i.e. length per time-square.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { "Gravity_R", "Gravity_ETA", "Gravity_XI" },
            new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Gravity;

        /// <summary>
        /// Residual in the momentum equation.
        /// </summary>
        [InstantiateFromControlFile(new string[] { "ResidualMomentum_R", "ResidualMomentum_ETA", "ResidualMomentum_XI" },
            new string[] { "Velocity_R", "Velocity_ETA", "Velocity_XI" },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> ResidualMomentum;

        /// <summary>
        /// Pressure
        /// </summary>
        [InstantiateFromControlFile(VariableNames.Pressure, null, IOListOption.ControlFileDetermined)]
        SinglePhaseField Pressure;

        /// <summary>
        /// Residual of the continuity equation
        /// </summary>
        [InstantiateFromControlFile("ResidualConti", VariableNames.Pressure, IOListOption.ControlFileDetermined)]
        SinglePhaseField ResidualContinuity;
#pragma warning restore 649

        /// <summary>
        /// Links edge tags (<see cref="BoSSS.Foundation.Grid.IGeometricalEdgeData.EdgeTags"/>) and
        /// boundary conditions in the control object (<see cref="BoSSS.Solution.Control.AppControl.BoundaryValues"/>).
        /// </summary>
        protected IncompressibleBoundaryCondMap boundaryCondMap;


        /// <summary>
        /// Declaration of the spatial operator
        /// </summary>
        protected override DifferentialOperator GetOperatorInstance(int D) {

            // instantiate boundary condition mapping
            // ======================================
            boundaryCondMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

            // instantiate operator
            // ====================
            string[] CodName = (new[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" }).GetSubVector(0, D).Cat("ResidualConti");

            var op = new DifferentialOperator(
                __DomainVar: VariableNames.VelocityVector(D).Cat(VariableNames.Pressure),
                __ParameterVar: VariableNames.GravityVector(D),
                __CoDomainVar: CodName,
                QuadOrderFunc: QuadOrderFunc.NonLinear(2));

            op.LinearizationHint = LinearizationHint.GetJacobiOperator;

            // Temporal Operator
            // =================

            var TempOp = new ConstantTemporalOperator(op, 0.0); // init with entire diagonal set to 0.0
            op.TemporalOperator = TempOp;

            for (int d = 0; d < D; d++)
                TempOp.SetDiagonal(CodName[d], Control.Density); // set momentum equation entries to density

            // Pressure Reference
            // ==================

            // if there is no Dirichlet boundary condition,
            // the mean value of the pressure is free:
            op.FreeMeanValue[VariableNames.Pressure] = !boundaryCondMap.DirichletPressureBoundary;

            // Momentum Equation
            // =================

            // convective part:
            {
                for (int d = 0; d < D; d++) {

                    var comps = op.EquationComponents[CodName[d]];

                    var ConvBulk = new LocalLaxFriedrichsConvection(D, boundaryCondMap, d, Control.Density, null);
                    comps.Add(ConvBulk); // bulk component
                }
            }

            // pressure part:
            {
                for (int d = 0; d < D; d++) {
                    var comps = op.EquationComponents[CodName[d]];
                    var pres = new PressureGradientLin_d(d, boundaryCondMap);
                    comps.Add(pres); // bulk component

                }
            }

            // viscous part:
            {
                for (int d = 0; d < D; d++) {
                    var comps = op.EquationComponents[CodName[d]];

                    double penalty_bulk = this.Control.PenaltySafety;

                    var Visc = new SipViscosity_GradU(penalty_bulk, d, D, boundaryCondMap,
                        ViscosityOption.ConstantViscosity,
                        constantViscosityValue: Control.Viscosity);
                    comps.Add(Visc); // bulk component GradUTerm
                }
            }


            // Continuity equation
            // ===================
            {
                for (int d = 0; d < D; d++) {
                    var src = new Divergence_DerivativeSource(d, D);
                    var flx = new Divergence_DerivativeSource_Flux(d, boundaryCondMap);
                    op.EquationComponents[CodName[D]].Add(src);
                    op.EquationComponents[CodName[D]].Add(flx);
                }


                //IBM_Op.EquationComponents["div"].Add(new PressureStabilization(1, 1.0 / this.Control.PhysicalParameters.mu_A));
            }

            // Gravity parameter
            // =================

            op.ParameterFactories.Add(delegate (IReadOnlyDictionary<string, DGField> DomainVarFields) {
                return D.ForLoop(d => (VariableNames.Gravity_d(d), this.Gravity[d] as DGField));
            });


            // commit & return
            // ===============
            op.Commit();
            return op;
        }


        /// <summary>
        /// Declares Velocity and Pressure as those variables that we want to solve for
        /// </summary>
        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            return Velocity.Cat(Pressure);
        }



        /// <summary>
        /// Returns the fields where we want to store our residuals
        /// </summary>
        public override IEnumerable<DGField> InstantiateResidualFields() {
            return ResidualMomentum.Cat(ResidualContinuity);
        }


        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            dt = Control.dtFixed;

            Console.WriteLine(" ");
            Console.WriteLine($"Starting time-step #{TimestepNo}, dt = {dt}...");

            Solve(phystime, dt);

            Console.WriteLine("done.");
            return dt;
        }
    }
}
