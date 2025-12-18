using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.TimeStepping;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using MPI.Wrappers;
using NUnit.Framework.Constraints;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Solution.Control.AppControl;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Level-set evolution, based on the Phasefield approach.
    /// For the evolution of the Phasefield a div-free velocity extension computed from Stokes equation is used <see cref="StokesExtension"/>.
    /// Due to the div-free velocity extension, the level-set evolution can be implemented as scalar convection,
    /// which has proven to be very stable in DG.
    /// </summary>
    /// <remarks>
    /// - implemented by MR, jan21
    /// </remarks>
    public class PhasefieldEvolver : ILevelSetEvolver {

        // <summary>
        /// ctor
        /// </summary>
        public PhasefieldEvolver(string levelSetName, int hMForder, int D, IncompressibleBoundaryCondMap bcMap, SolverWithLevelSetUpdaterControl control, double AgglomThreshold, IGridData grd) {
            for (int d = 0; d < D; d++) {
                if (!bcMap.bndFunction.ContainsKey(NSECommon.VariableNames.Velocity_d(d)))
                    throw new ArgumentException($"Missing boundary condition for variable {NSECommon.VariableNames.Velocity_d(d)}.");
            }

            this.SpatialDimension = D;
            this.AgglomThreshold = AgglomThreshold;
            this.m_HMForder = hMForder;
            this.levelSetName = levelSetName;
            this.bcmap = bcMap;
            parameters = NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            this.m_grd = grd;
            this.m_bndVals = control.BoundaryValues;
            m_control = control;

            using (FuncTrace ft = new FuncTrace()) {
                ft.Info(String.Format("mobility : {0}; cahn : {1}; theta : {2}", control.PhasefieldControl.diff, control.PhasefieldControl.cahn, control.PhasefieldControl.theta));
            }
        }

        SolverWithLevelSetUpdaterControl m_control;
        IGridData m_grd;
        int SpatialDimension;
        double AgglomThreshold;
        int m_HMForder;
        string levelSetName;
        string[] parameters;
        IncompressibleBoundaryCondMap bcmap;

        /// <summary>
        /// should only be the interface velocity vector; typically, a phase-averaged velocity.
        /// </summary>
        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// nix
        /// </summary>
        public IList<string> VariableNames => null;

        
        /// <summary>
        /// <see cref="MassCorrection(DualLevelSet, double, double, bool, IReadOnlyDictionary{string, DGField}, IReadOnlyDictionary{string, DGField})"/>
        /// </summary>
        public bool AfterMovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            return MassCorrection(levelSet, time, dt, incremental, DomainVarFields, ParameterVarFields);
        }

        static Dictionary<string,double> mass;
        /// <summary>
        /// The order parameter is conserved, we can try and "lift" the interface to conserve the mass, i.e. the area of one phase.
        /// Only alid when the mass of the phases doesnt change, not when evaporation or inflows of both phases are present.
        /// </summary>
        /// <param name="phaseInterface"></param>
        /// <param name="time"></param>
        /// <param name="dt"></param>
        /// <param name="incremental"></param>
        /// <param name="DomainVarFields"></param>
        /// <param name="ParameterVarFields"></param>
        /// <returns>
        /// - true: the level-set-field has been changed
        /// - false: the level-set-field remains unchanged
        /// </returns>
        public bool MassCorrection(DualLevelSet phaseInterface,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            bool changed = false;
            if (m_control.PhasefieldControl.CorrectionType != Correction.Mass)
                return changed;

            using (FuncTrace ft = new FuncTrace()) {
                if (mass == null) {
                    mass = new Dictionary<string, double>();
                    mass[phaseInterface.Identification] = ComputeMass(phaseInterface.Tracker); // based on cg level set, which is not yet updated!
                }

                var CorrectedDGLevelSet = phaseInterface.DGLevelSet.CloneAs();
                var CorrectedLsTrk = new LevelSetTracker(phaseInterface.Tracker.GridDat, phaseInterface.Tracker.CutCellQuadratureType, 1, new string[] { "A", "B" }, CorrectedDGLevelSet);

                Queue<double> massNew = new Queue<double>();
                massNew.Enqueue(ComputeMass(CorrectedLsTrk));
                double massUncorrected = massNew.Peek();

                double massDiff = massNew.Peek() - mass[phaseInterface.Identification];

                int i = 0;
                while (massDiff.Abs() > 1e-6) {
                    changed = true;
                    // FD sensitivity of the mass/area
                    double correction = Math.Sign(massDiff) * 1e-10;

                    // take the correction guess and calculate a forward difference to approximate the derivative
                    ProjectCorrection(CorrectedDGLevelSet, phaseInterface.DGLevelSet, correction);

                    // update LsTracker
                    CorrectedLsTrk.UpdateTracker(0.0);
                    massNew.Enqueue(ComputeMass(CorrectedLsTrk));

                    correction = -(massDiff) / ((-massNew.Dequeue() + massNew.Dequeue()) / (correction));

                    double initial = massDiff;
                    bool finished = false;
                    int k = 0;
                    //while (massDiff.Abs() - initial.Abs() >= 0.0 && step > 1e-12)
                    while (!finished) {
                        double step = Math.Pow(0.5, k);
                        // compute and project 
                        // step one calculate distance field phiDist = 0.5 * log(Max(1+c, eps)/Max(1-c, eps)) * sqrt(2) * Cahn
                        // step two project the new phasefield phiNew = tanh((cDist + correction)/(sqrt(2) * Cahn))
                        // ===================
                        ProjectCorrection(CorrectedDGLevelSet, phaseInterface.DGLevelSet, correction * step);

                        // update LsTracker
                        CorrectedLsTrk.UpdateTracker(0.0);
                        massNew.Enqueue(ComputeMass(CorrectedLsTrk));
                        massDiff = massNew.Peek() - mass[phaseInterface.Identification];

                        if (massDiff.Abs() < (1 - 1e-4 * step) * initial.Abs()) {
                            finished = true;

                            // update field
                            phaseInterface.DGLevelSet.Clear();
                            phaseInterface.DGLevelSet.Acc(1.0, CorrectedDGLevelSet);

                            ft.Info($"" +
                                $"converged with stepsize:  {step}, correction: {correction}\n" +
                                $"                     dM:  {massDiff}");
                            if (k > 0)
                                ft.Info($"Finished Linesearch in {k} iterations");
                        } else if (Math.Abs(correction * step) < 1e-15) {
                            ft.Info($"Linesearch failed after {k} iterations");
                            finished = true;
                            massDiff = 0.0;
                        } else {
                            massNew.Dequeue();
                            k++;
                        }
                    }
                    i++;
                }

                ft.Info($"Performed Mass Correction in {i} iterations: \n" +
                    $"\toriginal mass:      {mass[phaseInterface.Identification]:N6}\n" +
                    $"\tuncorrected mass:   {massUncorrected:N6}\n" +
                    $"\tcorrected mass:     {massNew.Dequeue():N6}");

            }
            return changed;
        }
        
        void ProjectCorrection(LevelSet PhiNew, LevelSet PhiOld, double correction) {
            PhiNew.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Phi
                    // -----------------------------
                    PhiOld.Evaluate(j0, Len, NS, result);

                    // compute the pointwise values of the new level set
                    // -----------------------------

                    result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * Math.Sqrt(2) * m_control.PhasefieldControl.cahn);
                    result.ApplyAll(x => Math.Tanh((x + correction) / (Math.Sqrt(2) * m_control.PhasefieldControl.cahn)));
                }
            );

            //Tecplot.Tecplot.PlotFields(new DGField[] { PhiNew , PhiOld}, "Phasefield", 0.0, 2);
        }

        double ComputeMass(LevelSetTracker levelSetTracker) {
            int order = 0;
            // area of bubble
            double area = 0.0;
            SpeciesId spcId = levelSetTracker.SpeciesIdS[1];
            var SchemeHelper = levelSetTracker.GetXDGSpaceMetrics(levelSetTracker.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, levelSetTracker.GridDat,
                vqs.Compile(levelSetTracker.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            area = area.MPISum();

            return area;
        }

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

        SinglePhaseField potential;
        SinglePhaseField[] extensionVelocity;
        XdgTimestepping.XdgTimestepping timestepper;
        private XdgTimestepping.XdgTimestepping GetTimestepper(SinglePhaseField levelSet, SinglePhaseField[] Velocity) {

            int D = this.m_grd.SpatialDimension;
            var m_bcMap = new BoundaryCondMap<BoundaryType>(this.m_grd, BoundaryTranslator(this.m_bndVals), "phi");

            potential = new SinglePhaseField(levelSet.Basis, "mu");
            var res_phi = new SinglePhaseField(levelSet.Basis, "Phasefield");
            var res_mu = new SinglePhaseField(levelSet.Basis, "Potential");

            var diffOp = new DifferentialOperator(new string[] { "phi", "mu" },
                Solution.NSECommon.VariableNames.VelocityVector(this.SpatialDimension),
                new string[] { "Phasefield", "Potential" },
                QuadOrderFunc.NonLinear(3));

            diffOp.EquationComponents["Phasefield"].Add(new phi_Flux(D, () => Velocity, m_bcMap));
            diffOp.EquationComponents["Phasefield"].Add(new phi_Diffusion(D, 2.6, m_control.PhasefieldControl.diff, 0.0, m_bcMap));

            diffOp.EquationComponents["Potential"].Add(new mu_Diffusion(D, 2.6, m_control.PhasefieldControl.cahn, m_bcMap, contactAngle: new Formula(m_control.PhasefieldControl.theta)));
            diffOp.EquationComponents["Potential"].Add(new mu_Source());

            double[] MassScales = { 1.0, 0.0 };
            diffOp.TemporalOperator = new ConstantTemporalOperator(diffOp, MassScales);
            diffOp.LinearizationHint = LinearizationHint.GetJacobiOperator;
            diffOp.IsLinear = false;
            diffOp.Commit();

            //var Timestepper = new BDFTimestepper(diffOp, new SinglePhaseField[] { levelSet, potential }, Velocity, 1, () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver(), false);

            XdgTimestepping.XdgTimestepping Timestepper = new XdgTimestepping.XdgTimestepping(
                diffOp,
                new List<SinglePhaseField> { levelSet, potential },
                new List<SinglePhaseField> { res_phi, res_mu },
                m_control.PhasefieldControl.TimeSteppingScheme);

            Timestepper.RegisterResidualLogger(new ResidualLogger(levelSet.GridDat.MpiRank, new DatabaseDriver(NullFileSystemDriver.Instance), Guid.Empty));

            return Timestepper;

        }

        IDictionary<string, BoundaryValueCollection> m_bndVals;
        /// <summary>
        /// Translate the Boundary Values of the Base Control to the appropriate Cahn-Hilliard Boundary Values
        /// </summary>
        /// <param name="_BndIn"></param>
        /// <returns></returns>
        private IDictionary<string, BoundaryValueCollection> BoundaryTranslator(in IDictionary<string, BoundaryValueCollection> _BndIn) {
            // automatische generierung passender Cahn Hilliard Bedingungen, wie funktioniert das mit den EdgeTagNames und Boundaryfunctions etc.?
            IDictionary<string, BoundaryValueCollection> BndOut = new Dictionary<string, BoundaryValueCollection>();
            BndOut = _BndIn;
            return BndOut;
        }

        int iTimestep;
        /// <summary>
        /// 
        /// </summary>
        /// <param name="levelSet"></param>
        /// <param name="time"></param>
        /// <param name="dt"></param>
        /// <param name="incremental"></param>
        /// <param name="DomainVarFields"></param>
        /// <param name="ParameterVarFields"></param>
        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (var tr = new FuncTrace()) {                

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

                var ExtVelBuilder = new StokesExtension.StokesExtension(D, this.bcmap, this.m_HMForder, this.AgglomThreshold, true, StokesExtentionBoundaryOption.useBcMap);
                ExtVelBuilder.SolveExtension(levelSet.LevelSetIndex, levelSet.Tracker, meanVelocity, extensionVelocity);

                if (timestepper == null) {
                    timestepper = GetTimestepper(levelSet.DGLevelSet, extensionVelocity);
                    //Phasefield = new Phasefield(null, levelSet.C0LevelSet, levelSet.DGLevelSet, levelSet.Tracker, extensionVelocity, m_grd, m_control, null);
                }
                timestepper.LsTrk.UpdateTracker(time);

                //Phasefield.UpdateFields(levelSet.C0LevelSet, levelSet.DGLevelSet, levelSet.Tracker, extensionVelocity, m_grd, m_control, null);
                if (!ReferenceEquals(timestepper.CurrentState.Fields[0], levelSet.DGLevelSet)) {
                    throw new Exception("Something went wrong with the internal pointer magic of the levelSetTracker. Definitely a weakness of ObjectOrientation.");
                }

                bool success = timestepper.Solve(time, dt);
                //Tecplot.Tecplot.PlotFields(timestepper.CurrentState.Cat(extensionVelocity), "Phasefield", time, 2);

                tr.Info("Phasefield evolver, total phi: " + LevelSetIntegral(levelSet.DGLevelSet));
            }

            double LevelSetIntegral(DGField phi) {
                // total concentration, careful above values are "old" when CorrectionTracker is not updated, this value is always "new"
                double concentration = 0.0;
                var tqs = new CellQuadratureScheme();
                CellQuadrature.GetQuadrature(new int[] { 1 }, phi.GridDat,
                    tqs.Compile(phi.GridDat, phi.Basis.Degree * 2 + 2),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        phi.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            concentration += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
                concentration = concentration.MPISum();
                return concentration;
            }
        }
    }
}
