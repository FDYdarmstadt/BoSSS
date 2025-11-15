using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.Smoothing;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.TimeStepping;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel.Design;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Level-set evolution using a div-free velocity extension computed from Stokes equation.
    /// Due to the div-free velocity extension, the level-set evolution can be implemented as scalar convection,
    /// which has proven to be very stable in DG.
    /// </summary>
    /// <remarks>
    /// - implemented by Fk, jan21
    /// - added AdamsBashforth by Lb, april21
    /// </remarks>
    public class StokesExtensionEvolver : ILevelSetEvolver {


        /// <summary>
        /// ctor
        /// </summary>
        public StokesExtensionEvolver(string levelSetName, int hMForder, int D, IncompressibleBoundaryCondMap bcMap, double AgglomThreshold, IGridData grd, bool fullStokes = true, StokesExtentionBoundaryOption useBCmap = StokesExtentionBoundaryOption.FreeSlipAtWall) {
            for(int d = 0; d < D; d++) {
                if(!bcMap.bndFunction.ContainsKey(NSECommon.VariableNames.Velocity_d(d)))
                    throw new ArgumentException($"Missing boundary condition for variable {NSECommon.VariableNames.Velocity_d(d)}.");
            }
            this.fullStokes = fullStokes;
            this.useBCmap = useBCmap;
            this.SpatialDimension = D;
            this.AgglomThreshold = AgglomThreshold;
            this.m_HMForder = hMForder;
            this.levelSetName = levelSetName;
            this.bcmap = bcMap;
            parameters = NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            this.m_grd = grd;
            timeStepOrder = 2;

            //this.ReInit_Period = ReInitPeriod;
            //ReInit_Control = new EllipticReInitAlgoControl();
        }

        IGridData m_grd;
        int SpatialDimension;
        double AgglomThreshold;
        int m_HMForder;
        string levelSetName;
        string[] parameters;
        int timeStepOrder;
        bool fullStokes;
        StokesExtentionBoundaryOption useBCmap;
        IncompressibleBoundaryCondMap bcmap;

        /// <summary>
        /// should only be the interface velocity vector; typically, a phase-averaged velocity.
        /// </summary>
        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// nix
        /// </summary>
        public IList<string> VariableNames => null;

        // nothing to do
        //public Func<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>, bool> AfterMovePhaseInterface => Reinitialize;

        /// <summary>
        /// re-initialization
        /// </summary>
        public bool AfterMovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            return Reinitialize(levelSet, time, dt, incremental, DomainVarFields, ParameterVarFields);
        }


        /// <summary>
        /// Provides access to the internally constructed extension velocity.
        /// <see cref="ILevelSetEvolver.InternalFields"/>
        /// </summary>
        public IDictionary<string, DGField> InternalFields { 
            get {
                var Ret = new Dictionary<string, DGField>();

                if(extensionVelocity != null) {
                    foreach(var f in extensionVelocity)
                        Ret.Add(f.Identification, f);
                }

                return Ret;
            }
        }

        ITimeStepper timeStepper;

        private ITimeStepper InitializeAdamsBashforth(SinglePhaseField levelSet, SinglePhaseField[] Velocity) {
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

            AdamsBashforth abf = new AdamsBashforth(diffOp, levelSet.Mapping, new CoordinateMapping(Velocity), timeStepOrder);
            return abf;
        }



        private ITimeStepper InitializeRungeKutta(SinglePhaseField levelSet, SinglePhaseField[] Velocity) {
            var diffOp = new DifferentialOperator(new string[] { "Phi" },
                Solution.NSECommon.VariableNames.VelocityVector(this.SpatialDimension),
                new string[] { "codom1" },
                QuadOrderFunc.NonLinear(1));
            diffOp.EquationComponents["codom1"].Add(new StokesExtension.ScalarTransportFlux(this.bcmap, this.SpatialDimension));
            diffOp.TemporalOperator = new ConstantTemporalOperator(diffOp, 1);
            diffOp.IsLinear = true;
            diffOp.Commit();

            var Timestepper = new RungeKutta(RungeKuttaScheme.TVD3, diffOp, levelSet.Mapping, new CoordinateMapping(Velocity));
            return Timestepper;

        }

        SinglePhaseField[] extensionVelocity;

        /// <summary>
        /// 
        /// </summary>
        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(var tr = new FuncTrace()) {
                var g = levelSet.Tracker.GridDat;
                int D = g.SpatialDimension;

                SinglePhaseField[] meanVelocity = D.ForLoop(
                    d => (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))]
                    );

                if(extensionVelocity == null) {
                    extensionVelocity = D.ForLoop(d => new SinglePhaseField(meanVelocity[d].Basis, $"ExtensionVelocity[{d}]"));
                } else {
                    foreach(var f in extensionVelocity)
                        f.Clear();
                }

                var ExtVelBuilder = new StokesExtension.StokesExtension(D, this.bcmap, this.m_HMForder, this.AgglomThreshold, fullStokes, useBCMap: useBCmap);
                ExtVelBuilder.SolveExtension(levelSet.LevelSetIndex, levelSet.Tracker, meanVelocity, extensionVelocity);

                if(timeStepper == null) {
                    //timeStepper = InitializeAdamsBashforth(levelSet.DGLevelSet, extensionVelocity);
                    timeStepper = InitializeRungeKutta(levelSet.DGLevelSet, extensionVelocity);
                }
                if(!ReferenceEquals(timeStepper.Mapping.Fields[0], levelSet.DGLevelSet)) {
                    throw new Exception("Something went wrong with the internal pointer magic of the levelSetTracker. Definitely a weakness of ObjectOrientation.");
                }

                double degree = timeStepper.Mapping.Fields[0].Basis.Degree;
                var dtCFL = g.ComputeCFLTime(extensionVelocity, dt * 1000);
                dtCFL *= 1.0 / (2 * degree + 1);

                if(dt > dtCFL) {

                    int N = (int)Math.Ceiling(dt/dtCFL);
                    double dtMod = dt / N;
                    tr.Warning($"CFL violation in level-set evolution: original dt = {dt}, max. stable CFL dt = {dtCFL}; performing {N} sub-steps with dt = {dtMod}");
                    for(int n  = 0; n < N; n++)
                        timeStepper.Perform(dtMod);

                } else {
                    timeStepper.Perform(dt);
                }

                tr.Info("time in LS evolver: " + timeStepper.Time);
            }
        }


        private EllipticReInitAlgoControl ReInit_Control = new EllipticReInitAlgoControl();
        private int ReInit_TimestepIndex = 0;
        private int ReInit_Period = 0;

        public void InitializeReInit(EllipticReInitAlgoControl RI_ctrl, int RI_period, int RI_tsI) {
            ReInit_Control = RI_ctrl;
            ReInit_Period = RI_period;
            ReInit_TimestepIndex = RI_tsI + 1;
        }

        /// <summary>
        /// 
        /// </summary>
        public bool Reinitialize(
            DualLevelSet phaseInterface,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {


            bool changed = false;
            ReInit_TimestepIndex++; // increment first
            if (phaseInterface.LevelSetIndex > 0) 
                return false; //skip the second level set

            // after level-set evolution and for initializing non-signed-distance level set fields
            if (ReInit_Period > 0 && ReInit_TimestepIndex % ReInit_Period == 0) {

                Console.WriteLine("Performing ReInit");
                ReInit_Control.Potential = ReInitPotential.BastingSingleWell;
                EllipticReInit.EllipticReInit ReInitPDE = new EllipticReInit.EllipticReInit(phaseInterface.Tracker, ReInit_Control, phaseInterface.DGLevelSet);
                ReInitPDE.ReInitialize(); // Restriction: phaseInterface.Tracker.Regions.GetCutCellSubGrid());

                //FastMarchReinit FastMarchReinitSolver = new FastMarchReinit(phaseInterface.DGLevelSet.Basis);
                //CellMask Accepted = phaseInterface.Tracker.Regions.GetCutCellMask();
                ////CellMask ActiveField = phaseInterface.Tracker.Regions.GetNearFieldMask(1);
                //CellMask NegativeField = phaseInterface.Tracker.Regions.GetSpeciesMask("A");
                //FastMarchReinitSolver.FirstOrderReinit(phaseInterface.DGLevelSet, Accepted, NegativeField, null);
                changed = true;
            }

            
            return changed;
        }
    }
}
