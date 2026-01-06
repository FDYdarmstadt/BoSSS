using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
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
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using MathNet.Numerics.Providers.LinearAlgebra;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
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
        public StokesExtensionEvolver(string levelSetName, int hMForder, int D, IncompressibleBoundaryCondMap bcMap, double AgglomThreshold, IGridData grd, 
            bool fullStokes = true, StokesExtentionBoundaryOption useBCmap = StokesExtentionBoundaryOption.FreeSlipAtWall, EllipticReInitAlgoControl ReInitControl = null, int ReInitPeriod = 0) {
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

            //SetupLogging(m_grd.MpiRank);
            LogValues = new double[LogValueNames.Length];

            ReInit_Control = ReInitControl != null ? ReInitControl : new EllipticReInitAlgoControl();
            this.ReInit_Period = ReInitPeriod;

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

                //var dgPltList = new List<DGField>();
                //SinglePhaseField oldDG = new SinglePhaseField(levelSet.DGLevelSet.Basis, "oldDGLevelSet");
                //oldDG.Acc(1.0, levelSet.DGLevelSet);
                //dgPltList.Add(oldDG);

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
                //dgPltList.AddRange(meanVelocity);

                var ExtVelBuilder = new StokesExtension.StokesExtension(D, this.bcmap, this.m_HMForder, this.AgglomThreshold, fullStokes, useBCMap: useBCmap);
                ExtVelBuilder.SolveExtension(levelSet.LevelSetIndex, levelSet.Tracker, meanVelocity, extensionVelocity);
                //dgPltList.AddRange(extensionVelocity);

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

                //SinglePhaseField newDG = new SinglePhaseField(levelSet.DGLevelSet.Basis, "newDGLevelSet");
                //newDG.Acc(1.0, levelSet.DGLevelSet);
                //dgPltList.Add(newDG);

                //int timestep = (int)(time / dt);
                //Tecplot.Tecplot.PlotFields(dgPltList, $"MovePhaseInterfaceStokesExt-{timestep}", time, 2);

                tr.Info("time in LS evolver: " + timeStepper.Time);
            }
        }


        private EllipticReInitAlgoControl ReInit_Control;
        private int ReInit_Period;
        private int ReInit_TimestepIndex = 0;

        public void SetReInitTimestepNumber(int RI_Timestep) {
            ReInit_TimestepIndex = RI_Timestep;
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
            using(var tr = new FuncTrace()) {
                //tr.InfoToConsole = true;

                bool changed = false;
                ReInit_TimestepIndex++;

                if(phaseInterface.LevelSetIndex > 0)
                    return false; //skip the second level set

                // after level-set evolution and for initializing non-signed-distance level set fields
                //if (ReInit_Period > 0 && ReInit_TimestepIndex % ReInit_Period == 0) {
                if(CheckForReInitialization(phaseInterface, ReInit_Control, tr, time, dt)) {

                    int timestep = (int)(time / dt);
                    Console.WriteLine("===============================================");
                    Console.WriteLine($"Performing ReInit in time step {timestep} ... ");
                    Console.WriteLine("===============================================");
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

        protected bool CheckForReInitialization(DualLevelSet phaseInterface, EllipticReInitAlgoControl ReInit_Control, 
            FuncTrace tr, double time, double dt) {

            // check for the Eikonal properties
            if(ReInit_Control.UseAdaptiveReInit) {

                // absolute value of the LS-Gradient
                var LsTrk = phaseInterface.Tracker;

                //var grdDat = phaseInterface.DGLevelSet.GridDat;
                //int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
                //BitArray fullBA = new BitArray(J);
                //fullBA.SetAll(true);
                //CellMask fullMask = new CellMask(grdDat, fullBA);

                var CCmask = LsTrk.Regions.GetCutCellMask();
                //phaseInterface.DGLevelSet.GradientNorm(out SinglePhaseField GradNormCC, CCmask);

                var NBmask = LsTrk.Regions.GetNearFieldMask(1);
                phaseInterface.DGLevelSet.GradientNorm(out SinglePhaseField AbsGradNB, NBmask);
                AbsGradNB.Identification = "AbsoluteLevelSetGradient";
                //Tecplot.Tecplot.PlotFields(new DGField[] { AbsGradNB }, $"Eikonal-{ReInit_TimestepIndex}", ReInit_TimestepIndex, 2);

                double NBarea = 0.0;
                double CCarea = 0.0;
                foreach(int j in NBmask.ItemEnum) {
                    NBarea += LsTrk.GridDat.iGeomCells.GetCellVolume(j);
                    if(CCmask.Contains(j))
                        CCarea += LsTrk.GridDat.iGeomCells.GetCellVolume(j);
                }

                Func<double[], double> oneFunc = X => 1.0;
                int order = phaseInterface.DGLevelSet.Basis.Degree * 2;
                double GradNormNB = AbsGradNB.L2Error(oneFunc.Vectorize(), order, new CellQuadratureScheme(UseDefaultFactories: true, domain: NBmask));
                double GradNormCC = AbsGradNB.L2Error(oneFunc.Vectorize(), order, new CellQuadratureScheme(UseDefaultFactories: true, domain: CCmask));

                //AbsGradNB.GetExtremalValues(out double minVal, out double maxVal);
                CellQuadratureScheme scheme = new CellQuadratureScheme(domain: NBmask);
                var rule = scheme.SaveCompile(LsTrk.GridDat, order);
                double[] localError = AbsGradNB.LocalLxError(
                    (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) { 
                        results.SetAll(1.0);
                    }, null, rule);
                for(int j = 0; j < localError.Length; j++) {
                    localError[j] /= LsTrk.GridDat.iGeomCells.GetCellVolume(j);
                }
                double minVal = localError.Min();
                double maxVal = localError.Max();

                double MeanAbsGradNB = AbsGradNB.GetMeanValueTotal(NBmask);
                double MeanAbsGradCC = AbsGradNB.GetMeanValueTotal(CCmask);

                //tr.Info($"L2-Norm (NB/CC) = {GradNormNB} / {GradNormCC}; minVal/maxVal = {minVal} / {maxVal}; MeanTotalValue (NB/CC) = {MeanAbsGradNB} / {MeanAbsGradCC}");

                Basis basis = new Basis(LsTrk.GridDat, 0);
                SinglePhaseField L2error = new SinglePhaseField(basis, "L2error-LSGgradient");
                CellQuadrature.GetQuadrature(
                    new int[] { 1 }, LsTrk.GridDat, (new CellQuadratureScheme(true, NBmask)).Compile(LsTrk.GridDat, order),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                        phaseInterface.DGLevelSet.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                        for(int i = 0; i < Length; i++) {
                            for(int j = 0; j < QR.NoOfNodes; j++) {
                                EvalResult[i, j, 0] = (EvalResult[i, j, 0] - 1.0).Pow2();
                            }
                        }
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                        for(int i = 0; i < Length; i++) {
                            L2error.SetMeanValue(i, ResultsOfIntegration[i, 0].Sqrt());
                        }
                    }
                ).Execute();
                //Tecplot.Tecplot.PlotFields(new DGField[] { AbsGradNB, L2error }, $"Eikonal-{ReInit_TimestepIndex}", ReInit_TimestepIndex, 2);

                // distance to CG-LevelSet
                double IntDist = 0.0;
                double IntArea = 0.0;

                //Basis basis = new Basis(LsTrk.GridDat, order);
                //SinglePhaseField DGLevelSetPow2 = new SinglePhaseField(basis);
                //DGLevelSetPow2.Clear();
                //DGLevelSetPow2.ProjectPow(1.0, phaseInterface.DGLevelSet, 2);

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
                CellQuadratureScheme LSQuad = SchemeHelper.GetLevelSetQuadScheme(phaseInterface.LevelSetIndex, CCmask, order);
                CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                    LSQuad.Compile(LsTrk.GridDat, order),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        phaseInterface.DGLevelSet.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                        for(int i = 0; i < Length; i++) {
                            for(int j = 0; j < QR.NoOfNodes; j++) {
                                EvalResult[i, j, 0] = EvalResult[i, j, 0].Abs();
                            }
                        }
                        EvalResult.ExtractSubArrayShallow(-1, -1, 1).SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            IntDist += ResultsOfIntegration[i, 0];
                            IntArea += ResultsOfIntegration[i, 1];
                        }
                    }
                ).Execute();

                //L2IntDist = L2IntDist.Sqrt();

                //tr.Info($"Distance (relative to interface area/length) from DG Interface to CG Interface: {IntDist/IntArea}");
                //if(this.Log != null) {
                //if(this.m_grd.MpiRank == 0) {
                //    this.Log = new StreamWriter(
                //        new FileStream(LogName, FileMode.Append, FileAccess.Write, FileShare.Read));

                //    int timestep = (int)(time / dt);
                //    string line = String.Format(LogFormat, timestep, time, GradNormNB, GradNormCC, MeanAbsGradNB, MeanAbsGradCC, minVal, maxVal, IntDist / IntArea);
                //    Log.WriteLine(line);
                //    Log.Flush();
                //    Log.Close();
                //}

                // store values for logging
                LogValues = new double[] { GradNormNB / NBarea, GradNormCC / CCarea, minVal, maxVal, MeanAbsGradNB, MeanAbsGradCC, IntDist / IntArea };

                if(GradNormNB / NBarea > 50.0) {
                    tr.Info($"L2-norm per narrowband area over thershhold");
                    return true;
                }

                if(GradNormCC / CCarea > 100.0) {
                    tr.Info($"L2-norm per cut-cell area over thershhold");
                    return true;
                }

                if(maxVal > 0.2) {
                    tr.Info($"cell-local L2-norm max value over thershhold");
                    return true;
                }


            } else {
                if(ReInit_Period > 0 && ReInit_TimestepIndex % ReInit_Period == 0) {
                    return true;
                }
            }
            return false;
        }


        #region logging

        string[] LogValueNames = new string[] { "L2-Norm_NB (rel)", "L2-Norm_CC (rel)", "L2minVal(cell,rel)", "L2maxVal(cell,rel)", "MeanTotalValue_NB", "MeanTotalValue_CC", "InterfaceDist(rel)" };

        double[] LogValues;

        public double[] GetLogValues() { 
            return LogValues;
        }

        //TextWriter Log;

        //string LogName = "StokesExtensionEvolver_Logging.txt";

        //string LogFormat = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}";

        //void SetupLogging(int MPIrank) {

        //    if(MPIrank == 0) {
        //        this.Log = new StreamWriter(
        //                        new FileStream(LogName, FileMode.Append, FileAccess.Write, FileShare.Read));

        //        string header = String.Format(LogFormat, "#timestep", "time", "L2-Norm_NB", "L2-Norm_CC", "MeanTotalValue_NB", "MeanTotalValue_CC", "minVal", "maxVal", "InterfaceDist");
        //        Log.WriteLine(header);
        //        Log.Flush();
        //        Log.Close();
        //    }
        //}

        #endregion
    }
}
