using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.TimeStepping;
using ilPSP;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Level-set evolution using a div-free velocity extension computed from Stokes equation.
    /// Due to the div-free velocity extension, the level-set evolution can be implemented as scalar convection,
    /// which has proven to be very stable in DG.
    /// </summary>
    /// <remarks>
    /// - original implemented by Fk, jan21, <see cref="StokesExtensionEvolver"/>
    /// - cloned and modified by Lauritz, april21
    /// </remarks>
    public class ImplicitStokesExtensionEvolver : ILevelSetEvolver {


        /// <summary>
        /// ctor
        /// </summary>
        public ImplicitStokesExtensionEvolver(string levelSetName, int hMForder, int D, IncompressibleBoundaryCondMap bcMap, double AgglomThreshold, IGridData grd, bool fullStokes = true) {
            for (int d = 0; d < D; d++) {
                if (!bcMap.bndFunction.ContainsKey(NSECommon.VariableNames.Velocity_d(d)))
                    throw new ArgumentException($"Missing boundary condition for variable {NSECommon.VariableNames.Velocity_d(d)}.");
            }
            this.fullStokes = fullStokes;
            this.SpatialDimension = D;
            this.AgglomThreshold = AgglomThreshold;
            this.m_HMForder = hMForder;
            this.levelSetName = levelSetName;
            this.bcmap = bcMap;
            parameters = NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            this.m_grd = grd;
        }

        IGridData m_grd;
        int SpatialDimension;
        double AgglomThreshold;
        int m_HMForder;
        string levelSetName;
        string[] parameters;
        IncompressibleBoundaryCondMap bcmap;
        bool fullStokes;

        /// <summary>
        /// should only be the interface velocity vector; typically, a phase-averaged velocity.
        /// </summary>
        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// nix
        /// </summary>
        public IList<string> VariableNames => null;

        
        /// <summary>
        /// Currently empty
        /// </summary>
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => null;


        /// <summary>
        /// Provides access to the internally constructed extension velocity.
        /// <see cref="ILevelSetEvolver.InternalFields"/>
        /// </summary>
        public IDictionary<string, DGField> InternalFields {
            get {
                var Ret = new Dictionary<string, DGField>();

                if (extensionVelocity != null) {
                    foreach (var f in extensionVelocity)
                        Ret.Add(f.Identification, f);
                }

                return Ret;
            }
        }


        BDFTimestepper implicitTimeStepper;

        RungeKutta explicitTimeStepper;

        private (BDFTimestepper, RungeKutta) InitializeTimeStepper(SinglePhaseField levelSet, SinglePhaseField[] Velocity) {
            var diffOp = new DifferentialOperator(new string[] { "Phi" },
                Solution.NSECommon.VariableNames.VelocityVector(this.SpatialDimension),
                new string[] { "codom1" },
                QuadOrderFunc.NonLinear(1));
            diffOp.EquationComponents["codom1"].Add(new StokesExtension.ScalarTransportFlux(this.bcmap, this.SpatialDimension));
            diffOp.TemporalOperator = new ConstantTemporalOperator(diffOp, 1);
            diffOp.IsLinear = true;
            diffOp.Commit();

            ISparseSolver Solver() {
                return new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
            }

            BDFTimestepper bdf = new BDFTimestepper(diffOp, levelSet.Mapping, Velocity, 3, Solver, false);
            RungeKutta rk = new RungeKutta(RungeKuttaScheme.TVD3, diffOp, levelSet.Mapping, new CoordinateMapping(Velocity));
            return (bdf, rk);
        }


        SinglePhaseField[] extensionVelocity;

        /// <summary>
        /// 
        /// </summary>
        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            int D = levelSet.Tracker.GridDat.SpatialDimension;

            SinglePhaseField[] meanVelocity = D.ForLoop(
                d => (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))]
                );

            if (extensionVelocity == null) {
                extensionVelocity = D.ForLoop(d => new SinglePhaseField(meanVelocity[d].Basis, $"ExtensionVelocity[{d}]"));
            } else {
                foreach (var f in extensionVelocity)
                    f.Clear();
            }

            var ExtVelBuilder = new StokesExtension.StokesExtension(D, this.bcmap, this.m_HMForder, this.AgglomThreshold, fullStokes);
            ExtVelBuilder.SolveExtension(levelSet.LevelSetIndex, levelSet.Tracker, meanVelocity, extensionVelocity);

            if (implicitTimeStepper == null) {
                (implicitTimeStepper, explicitTimeStepper)  = InitializeTimeStepper(levelSet.DGLevelSet, extensionVelocity);
            }
            if (!ReferenceEquals(implicitTimeStepper.Mapping.Fields[0], levelSet.DGLevelSet) 
                || ! ReferenceEquals(explicitTimeStepper.Mapping.Fields[0], levelSet.DGLevelSet)) {
                throw new Exception("Something went wrong with the internal pointer magic of the levelSetTracker. Definitely a weakness of ObjectOrientation.");
            }
            if(Math.Abs(time - internalTime) < dt * 1e-10 || first) {
                implicitTimeStepper.FinishTimeStep();
                internalTime = time + dt;
                first = false;
                Console.WriteLine("pushed LevelSet! Created initial guess!");
                //explicitTimeStepper.Perform(dt);
            }
            Console.WriteLine("iterating levelSetposition");
            implicitTimeStepper.Perform(dt);
        }
        bool first = true;
        double internalTime = 0.0;
    }
}
