using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE {

    /// <summary>
    /// A minimal solver for the incompressible Navier-Stokes equation.
    /// </summary>
    public class IncompressibleNSEMain : BoSSS.Solution.XdgTimestepping.DgApplicationWithSolver<IncompressibleControl> {

        static void Main(string[] args) {
            _Main(args, false, delegate () {
                var p = new IncompressibleNSEMain();
                return p;
            });
        }

#pragma warning disable 649
        /// <summary>
        /// velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            null,
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Velocity;

        /// <summary>
        /// Volume Force, dimension is acceleration, i.e. length per time-square.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { VariableNames.GravityX, VariableNames.GravityY, VariableNames.GravityZ },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Gravity;

        /// <summary>
        /// Residual in the momentum equation.
        /// </summary>
        [InstantiateFromControlFile(new string[] { "ResidualMomentumX", "ResidualMomentumY", "ResidualMomentumZ" },
            new string[] { VariableNames.VelocityX, VariableNames.VelocityY, VariableNames.VelocityZ },
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
