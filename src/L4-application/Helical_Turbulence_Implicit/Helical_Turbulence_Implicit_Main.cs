using BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations;
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


            var c = BoSSS.Application.IncompressibleNSE.DNS_Hagen_Poiseulle.HagenPoiseulle();
            c.ImmediatePlotPeriod = 1;
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

        /// <summary>
        /// velocity
        /// </summary>
        [InstantiateFromControlFile(new string[] { "Velocity_R", "Velocity_XI", "Velocity_ETA" },
            null,
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Velocity;

        /// <summary>
        /// Volume Force, dimension is acceleration, i.e. length per time-square.
        /// </summary>
        [InstantiateFromControlFile(
            new string[] { "Gravity_R", "Gravity_XI", "Gravity_ETA"},
            new string[] { "Velocity_R", "Velocity_XI", "Velocity_ETA" },
            true, true,
            IOListOption.ControlFileDetermined)]
        public VectorField<SinglePhaseField> Gravity;

        /// <summary>
        /// Residual in the momentum equation.
        /// </summary>
        [InstantiateFromControlFile(new string[] { "ResidualMomentum_R", "ResidualMomentum_XI", "ResidualMomentum_ETA" },
            new string[] { "Velocity_R", "Velocity_XI", "Velocity_ETA" },
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
        /// vector of helical velocity names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] VelocityVector_Helical(int D) {
            if (D == 2)
                return new string[] { "Velocity_R", "Velocity_XI" };
            else if (D == 3)
                return new string[] { "Velocity_R",  "Velocity_XI","Velocity_ETA" };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// vector of helical gravity/volume force names
        /// </summary>
        /// <param name="D">
        /// spatial dimension
        /// </param>
        public static string[] GravityVector_Helical(int D) {
            if (D == 2)
                return new string[] { "Gravity_R", "Gravity_XI" };
            else if (D == 3)
                return new string[] { "Gravity_R", "Gravity_XI" ,"Gravity_ETA" };
            else
                throw new NotSupportedException("unsupported spatial dimension: D = " + D + ".");
        }

        /// <summary>
        /// Declaration of the spatial HELICAL operator
        /// </summary>
        protected override DifferentialOperator GetOperatorInstance(int D_) {
            int D = D_ +1;

            // instantiate boundary condition mapping
            // ======================================
            boundaryCondMap = new IncompressibleBoundaryCondMap(this.GridData, this.Control.BoundaryValues, PhysicsMode.Incompressible);

            // instantiate operator
            // ====================
            string[] CodName = (new[] { "ResidualMomentum_R",  "ResidualMomentum_XI", "ResidualMomentum_ETA" }).GetSubVector(0, D).Cat("ResidualConti");

            // instantiate Values 
            // ====================
            int order = this.Control.dg_degree;
            double penalty = this.Control.penaltySafety * order * order;
            int noOfCells = this.Control.Resolution_Xi;
            Globals.MaxAmp = this.Control.maxAmpli;

            var op = new DifferentialOperator(
                __DomainVar: VelocityVector_Helical(D).Cat(VariableNames.Pressure),
                __ParameterVar: GravityVector_Helical(D),
                __CoDomainVar: CodName,
                QuadOrderFunc: QuadOrderFunc.NonLinear(2));

            op.LinearizationHint = LinearizationHint.GetJacobiOperator;

            // Temporal Operator
            // =================
            var tmpOp = new DependentTemporalOperator(op);

            for (int d = 0; d < D; d++) {
                var comps = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.TransientTerm(1, d);
                tmpOp.EquationComponents[CodName[d]].Add(comps); // set momentum equation entries to density
            }
            op.TemporalOperator = tmpOp;

            // Pressure Reference
            // ==================

            // if there is no Dirichlet boundary condition,
            // the mean value of the pressure is free:
            op.FreeMeanValue[VariableNames.Pressure] = !boundaryCondMap.DirichletPressureBoundary;

            // Momentum Equation
            // =================
            // convective part:
            {
                var comps_conv_R = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations.convectiveRmom();
                op.EquationComponents["ResidualMomentum_R"].Add(comps_conv_R);       // bulk component

                var comps_conv_XI = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations.convectiveXImom();
                op.EquationComponents["ResidualMomentum_XI"].Add(comps_conv_XI);     // bulk component

                var comps_conv_ETA = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations.convectiveETAmom();
                op.EquationComponents["ResidualMomentum_ETA"].Add(comps_conv_ETA);   // bulk component
            }
            // pressure part:
            // ===================
            {
                var comps_pressure_R = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations.Gradient_pressure_rMomentum(0);
                op.EquationComponents["ResidualMomentum_R"].Add(comps_pressure_R);       // bulk component
                var comps_pressure_XI = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations.Gradient_pressure_xiMomentum(1);
                op.EquationComponents["ResidualMomentum_XI"].Add(comps_pressure_XI);       // bulk component
            }

            // viscous part:
            // ===================
            {
                var comps_visc_R = new rMomentum(this.Control.TermSwitch, penalty, PenaltyFactor);
                op.EquationComponents["ResidualMomentum_R"].Add(comps_visc_R);      // bulk component
                var comps_visc_XI = new xiMomentum(this.Control.TermSwitch, penalty, PenaltyFactor);
                op.EquationComponents["ResidualMomentum_XI"].Add(comps_visc_XI);   // bulk component
                var comps_visc_ETA = new etaMomentum(noOfCells, this.Control.TermSwitch, penalty, PenaltyFactor);
                op.EquationComponents["ResidualMomentum_ETA"].Add(comps_visc_ETA);  // bulk component
            }

            // Continuity equation
            // ===================
            {
                var myContiNew = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ContinuityEquation.Conti(noOfCells);
                op.EquationComponents["ResidualConti"].Add(myContiNew);
            }
            // Forcing Terms
            // ===================
            {
                if (Control.HagenPoisseulle) {
                    var comps_forcing_ETA = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ForcingTerms.Hagen_Poiseulle.ForcingTermEta_Hagen();
                    op.EquationComponents["ResidualMomentum_ETA"].Add(comps_forcing_ETA);
                    var comps_forcing_XI = new BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ForcingTerms.Hagen_Poiseulle.ForcingTermXi_Hagen();
                    op.EquationComponents["ResidualMomentum_XI"].Add(comps_forcing_XI);
                }
            }
            Console.WriteLine("Wo soll ich da zuweisen?!");
            Globals.DirichletValue_uR = (X => this.Control.BoundaryValues["Dirichlet_inner_wall"].Values["Velocity_R"].Evaluate(X, 0));
            Globals.DirichletValue_uXi = (X => this.Control.BoundaryValues["Dirichlet_inner_wall"].Values["Velocity_XI"].Evaluate(X, 0));
            Globals.DirichletValue_uEta = (X => this.Control.BoundaryValues["Dirichlet_inner_wall"].Values["Velocity_ETA"].Evaluate(X, 0));
            Globals.BoundaryType = X => BoundaryTypeE.Dirichlet;
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

        public double PenaltyFactor(double penaltyFactor, int jCellIn, int jCellOut, MultidimensionalArray cj) {
            double cj_in = cj[jCellIn];
            double eta = penaltyFactor * cj_in;
            if (jCellOut >= 0) {
                double cj_out = cj[jCellOut];
                eta = Math.Max(eta, penaltyFactor * cj_out);
            }
            return eta;
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
