using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
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
        public StokesExtensionEvolver(string levelSetName, int hMForder, int D, IncompressibleBoundaryCondMap bcMap, double AgglomThreshold, IGridData grd, bool fullStokes = true) {
            for(int d = 0; d < D; d++) {
                if(!bcMap.bndFunction.ContainsKey(NSECommon.VariableNames.Velocity_d(d)))
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
            timeStepOrder = 2;
        }

        IGridData m_grd;
        int SpatialDimension;
        double AgglomThreshold;
        int m_HMForder;
        string levelSetName;
        string[] parameters;
        int timeStepOrder;
        bool fullStokes;
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
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => Reinitialize;


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
            var diffOp = new SpatialOperator(new string[] { "Phi" },
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
            var diffOp = new SpatialOperator(new string[] { "Phi" },
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
        static int cnt = 0;
        /// <summary>
        /// 
        /// </summary>
        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(var tr = new FuncTrace()) {

                int D = levelSet.Tracker.GridDat.SpatialDimension;

                SinglePhaseField[] meanVelocity = D.ForLoop(
                    d => (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))]
                    );

                if(extensionVelocity == null) {
                    extensionVelocity = D.ForLoop(d => new SinglePhaseField(meanVelocity[d].Basis, $"ExtensionVelocity[{d}]"));
                } else {
                    foreach(var f in extensionVelocity)
                        f.Clear();
                }

                var ExtVelBuilder = new StokesExtension.StokesExtension(D, this.bcmap, this.m_HMForder, this.AgglomThreshold, fullStokes);
                ExtVelBuilder.SolveExtension(levelSet.LevelSetIndex, levelSet.Tracker, meanVelocity, extensionVelocity);

                if(timeStepper == null) {
                    //timeStepper = InitializeAdamsBashforth(levelSet.DGLevelSet, extensionVelocity);
                    timeStepper = InitializeRungeKutta(levelSet.DGLevelSet, extensionVelocity);
                }
                if(!ReferenceEquals(timeStepper.Mapping.Fields[0], levelSet.DGLevelSet)) {
                    throw new Exception("Something went wrong with the internal pointer magic of the levelSetTracker. Definitely a weakness of ObjectOrientation.");
                }

                timeStepper.Perform(dt);

                tr.Info("time in LS evolver: " + timeStepper.Time );
            }
        }

        static int lastReinit = 0;
        public void Reinitialize(
            DualLevelSet phaseInterface,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            int rcnt = 0;

            // Not implemented
            REINIT:

            double Jump_NORM = 0.0;
            double GradJump_NORM = 0.0;

            EdgeMask Edges = phaseInterface.Tracker.Regions.GetCutCellMask().GetAllInnerEdgesMask();
            EdgeQuadratureScheme cqs = new EdgeQuadratureScheme(true, Edges);

            EdgeQuadrature.GetQuadrature(new int[] { 2 }, phaseInterface.DGLevelSet.GridDat,
                cqs.Compile(phaseInterface.DGLevelSet.GridDat, 2 * phaseInterface.DGLevelSet.Basis.Degree),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray PhiIN = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray PhiOUT = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray PhiGradIN = MultidimensionalArray.Create(Length, QR.NoOfNodes, phaseInterface.DGLevelSet.GridDat.SpatialDimension);
                    MultidimensionalArray PhiGradOUT = MultidimensionalArray.Create(Length, QR.NoOfNodes, phaseInterface.DGLevelSet.GridDat.SpatialDimension);

                    phaseInterface.DGLevelSet.EvaluateEdge(i0, Length, QR.Nodes, PhiIN, PhiOUT, null, null, PhiGradIN, PhiGradOUT);
                    MultidimensionalArray PhiGradJump = PhiGradIN - PhiGradOUT;
                    MultidimensionalArray PhiJump = PhiIN - PhiOUT;


                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Multiply(1.0, PhiJump, PhiJump, 0.0, "ik", "ik", "ik");
                    EvalResult.ExtractSubArrayShallow(-1, -1, 1).Multiply(1.0, PhiGradJump, PhiGradJump, 0.0, "ik", "ikj", "ikj");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        Jump_NORM += ResultsOfIntegration[i, 0];
                        GradJump_NORM += ResultsOfIntegration[i, 1];
                    }
                }
            ).Execute();

            Jump_NORM = Jump_NORM.MPISum().Sqrt();
            GradJump_NORM = GradJump_NORM.MPISum().Sqrt();

            Console.WriteLine("DG-LevelSet Jump : {0} || Gradient Jump : {1}", Jump_NORM, GradJump_NORM);

            Jump_NORM = 0.0;
            GradJump_NORM = 0.0;

            // doesnt make sense, CG levelset isnt updated yet ...
            EdgeQuadrature.GetQuadrature(new int[] { 2 }, phaseInterface.DGLevelSet.GridDat,
                cqs.Compile(phaseInterface.DGLevelSet.GridDat, 2 * phaseInterface.DGLevelSet.Basis.Degree),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray PhiIN = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray PhiOUT = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray PhiGradIN = MultidimensionalArray.Create(Length, QR.NoOfNodes, phaseInterface.DGLevelSet.GridDat.SpatialDimension);
                    MultidimensionalArray PhiGradOUT = MultidimensionalArray.Create(Length, QR.NoOfNodes, phaseInterface.DGLevelSet.GridDat.SpatialDimension);

                    phaseInterface.CGLevelSet.EvaluateEdge(i0, Length, QR.Nodes, PhiIN, PhiOUT, null, null, PhiGradIN, PhiGradOUT);
                    MultidimensionalArray PhiGradJump = PhiGradIN - PhiGradOUT;
                    MultidimensionalArray PhiJump = PhiIN - PhiOUT;


                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Multiply(1.0, PhiJump, PhiJump, 0.0, "ik", "ik", "ik");
                    EvalResult.ExtractSubArrayShallow(-1, -1, 1).Multiply(1.0, PhiGradJump, PhiGradJump, 0.0, "ik", "ikj", "ikj");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        Jump_NORM += ResultsOfIntegration[i, 0];
                        GradJump_NORM += ResultsOfIntegration[i, 1];
                    }
                }
            ).Execute();

            Jump_NORM = Jump_NORM.MPISum().Sqrt();
            GradJump_NORM = GradJump_NORM.MPISum().Sqrt();

            Console.WriteLine("CG-LevelSet Jump : {0} || Gradient Jump : {1}", Jump_NORM, GradJump_NORM);

            double Cont_DIFF = 0.0;
            double Cont_CTRL = 0.0;
            double Cont_LEN = 0.0;

            int id = phaseInterface.LevelSetIndex;
            CellQuadratureScheme lqs = phaseInterface.Tracker.GetXDGSpaceMetrics(phaseInterface.Tracker.SpeciesIdS, 2 * phaseInterface.DGLevelSet.Basis.Degree).XQuadSchemeHelper.GetLevelSetquadScheme(id, phaseInterface.Tracker.Regions.GetCutCellMask4LevSet(id));

            CellQuadrature.GetQuadrature(new int[] { 3 }, phaseInterface.DGLevelSet.GridDat,
                lqs.Compile(phaseInterface.DGLevelSet.GridDat, 2 * phaseInterface.DGLevelSet.Basis.Degree),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray PhiDG = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray PhiCG = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray PhiCGGrad = MultidimensionalArray.Create(Length, QR.NoOfNodes, phaseInterface.DGLevelSet.GridDat.SpatialDimension);
                    MultidimensionalArray PhiCGGradNorm = MultidimensionalArray.Create(Length, QR.NoOfNodes);

                    phaseInterface.DGLevelSet.Evaluate(i0, Length, QR.Nodes, PhiDG);
                    phaseInterface.CGLevelSet.Evaluate(i0, Length, QR.Nodes, PhiCG);
                    phaseInterface.CGLevelSet.EvaluateGradient(i0, Length, QR.Nodes, PhiCGGrad);

                    // difference
                    PhiDG.Acc(-1.0, PhiCG);
                    PhiDG.ApplyAll(x => Math.Abs(x));

                    // control, should be near zero
                    PhiCG.ApplyAll(x => Math.Abs(x));

                    // normalize by cg levelset gradient norm
                    PhiCGGradNorm.Multiply(1.0, PhiCGGrad, PhiCGGrad, 0.0, "ik", "ikj", "ikj");
                    PhiCGGradNorm.ApplyAll(x => 1.0/Math.Sqrt(x));

                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Multiply(1.0, PhiDG, PhiCGGradNorm, 0.0, "ik", "ik", "ik");
                    EvalResult.ExtractSubArrayShallow(-1, -1, 1).Multiply(1.0, PhiCG, PhiCGGradNorm, 0.0, "ik", "ik", "ik");
                    EvalResult.ExtractSubArrayShallow(-1, -1, 2).SetAll(1.0); // length of interface

                }, delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        Cont_DIFF += ResultsOfIntegration[i, 0];
                        Cont_CTRL += ResultsOfIntegration[i, 1];
                        Cont_LEN += ResultsOfIntegration[i, 2];
                    }
                }).Execute();

            Cont_DIFF *= 1.0 / Cont_LEN;
            Cont_CTRL *= 1.0 / Cont_LEN;

            Console.WriteLine("CG-DG Contour Difference : {0} || Contour Control : {1} || Contour Length : {2}", Cont_DIFF, Cont_CTRL, Cont_LEN);

            if(GradJump_NORM > 1e0 && lastReinit < 0 && false) {
                Console.WriteLine("Performing Reinit");
                var ReInit_Control = new EllipticReInitAlgoControl();
                ReInit_Control.Potential = ReInitPotential.BastingSingleWell;
                EllipticReInit.EllipticReInit ReInitPDE = new EllipticReInit.EllipticReInit(phaseInterface.Tracker, ReInit_Control, phaseInterface.DGLevelSet);
                ReInitPDE.ReInitialize(Restriction: phaseInterface.Tracker.Regions.GetCutCellSubGrid());

                BoSSS.Solution.LevelSetTools.Reinit.FastMarch.FastMarchReinit FastMarchReinitSolver = new BoSSS.Solution.LevelSetTools.Reinit.FastMarch.FastMarchReinit(phaseInterface.DGLevelSet.Basis);
                CellMask Accepted = phaseInterface.Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex);
                CellMask ActiveField = phaseInterface.Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex);
                CellMask NegativeField = phaseInterface.Tracker.Regions.GetSpeciesMask("A");
                do {
                    ActiveField = ActiveField.Union(ActiveField.AllNeighbourCells());
                    FastMarchReinitSolver.FirstOrderReinit(phaseInterface.DGLevelSet, Accepted, NegativeField, ActiveField);
                    Accepted = Accepted.Union(ActiveField);
                } while (ActiveField != CellMask.GetFullMask(m_grd));
                lastReinit = 100;
            }
            lastReinit--;

            double jTrsh = 0.1;
            double gTrsh = 10.0;

            //if (rcnt == 0) {
            //    var jp = new JumpPenalization(JumpPenalization.jumpPenalizationTerms.JumpGradJump2, Jump_NORM * GradJump_NORM);
            //    jp.ImplicitEuler(dt * 0.001, new SubGrid(CellMask.GetFullMask(phaseInterface.DGLevelSet.GridDat)), phaseInterface.DGLevelSet);
            //    rcnt++;
            //    if (rcnt < 100) goto REINIT;
            //}

            //if ((GradJump_NORM > gTrsh) || (Jump_NORM > jTrsh)) {
            //    var jp = new JumpPenalization(JumpPenalization.jumpPenalizationTerms.GradJump2, 10.0);
            //    jp.ImplicitEuler(dt * 0.001, new SubGrid(CellMask.GetFullMask(phaseInterface.DGLevelSet.GridDat)), phaseInterface.DGLevelSet);
            //    rcnt++;
            //    if (rcnt < 100) goto REINIT;
            //}

            //if ((GradJump_NORM > gTrsh) && (Jump_NORM > jTrsh)) {
            //    var jp = new JumpPenalization(JumpPenalization.jumpPenalizationTerms.GradJump2, 1.0);
            //    jp.ImplicitEuler(dt * 0.001, new SubGrid(CellMask.GetFullMask(phaseInterface.DGLevelSet.GridDat)), phaseInterface.DGLevelSet);
            //    rcnt++;
            //    if (rcnt < 100) goto REINIT;
            //} else if (Jump_NORM > jTrsh) {
            //    var jp = new JumpPenalization(JumpPenalization.jumpPenalizationTerms.Jump, 1.0);
            //    jp.ImplicitEuler(dt * 0.001, new SubGrid(CellMask.GetFullMask(phaseInterface.DGLevelSet.GridDat)), phaseInterface.DGLevelSet);
            //    rcnt++;
            //    if (rcnt < 100) goto REINIT;
            //} else if (GradJump_NORM > gTrsh) {
            //    var jp = new JumpPenalization(JumpPenalization.jumpPenalizationTerms.GradJump2, 1.0);
            //    jp.ImplicitEuler(dt * 0.001, new SubGrid(CellMask.GetFullMask(phaseInterface.DGLevelSet.GridDat)), phaseInterface.DGLevelSet);
            //    rcnt++;
            //    if (rcnt < 100) goto REINIT;
            //}
        }
        }
}
