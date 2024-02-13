using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.TimeStepping;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
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
    public  class PhasefieldEvolver : ILevelSetEvolver {

        // <summary>
        /// ctor
        /// </summary>
        public PhasefieldEvolver(string levelSetName, int hMForder, int D, IncompressibleBoundaryCondMap bcMap, AppControl control, double AgglomThreshold, IGridData grd) {
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
        }
        AppControl m_control;
        IGridData m_grd;
        int SpatialDimension;
        double AgglomThreshold;
        int m_HMForder;
        string levelSetName;
        string[] parameters;
        IncompressibleBoundaryCondMap bcmap;
        Phasefield Phasefield;

        /// <summary>
        /// should only be the interface velocity vector; typically, a phase-averaged velocity.
        /// </summary>
        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// nix
        /// </summary>
        public IList<string> VariableNames => null;

        // nothing to do
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

        SinglePhaseField potential;
        SinglePhaseField[] extensionVelocity;
        BDFTimestepper bdf;
        private BDFTimestepper GetTimestepper(SinglePhaseField levelSet, SinglePhaseField[] Velocity) {

            int D = this.m_grd.SpatialDimension;
            var m_bcMap = new BoundaryCondMap<BoundaryType>(this.m_grd, BoundaryTranslator(this.m_bndVals), "phi");

            var diffOp = new DifferentialOperator(new string[] { "phi", "mu" },
                Solution.NSECommon.VariableNames.VelocityVector(this.SpatialDimension),
                new string[] { "Phasefield", "Potential" },
                QuadOrderFunc.NonLinear(2));


            //diffOp.EquationComponents["Phasefield"].Add(new phi_Flux(D, m_bcMap));
            diffOp.EquationComponents["Phasefield"].Add(new phi_Diffusion(D, 2.6, -0.1, 0.0, m_bcMap));

            diffOp.EquationComponents["Potential"].Add(new mu_Diffusion(D, 2.6, m_cahn, m_bcMap));
            diffOp.EquationComponents["Potential"].Add(new mu_Source());


            double[] MassScales = { 1.0, 0.0 };
            diffOp.TemporalOperator = new ConstantTemporalOperator(diffOp, MassScales);
            diffOp.LinearizationHint = LinearizationHint.GetJacobiOperator;
            diffOp.IsLinear = false;
            diffOp.Commit();

            var Timestepper = new BDFTimestepper(diffOp, new SinglePhaseField[] { levelSet, potential }, Velocity, 1, () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver(), false);
            
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

            var ExtVelBuilder = new StokesExtension.StokesExtension(D, this.bcmap, this.m_HMForder, this.AgglomThreshold, true);
            ExtVelBuilder.SolveExtension(levelSet.LevelSetIndex, levelSet.Tracker, meanVelocity, extensionVelocity);

            // Here Phasefield Movement            
            if (Phasefield == null) {
                iTimestep = 0;
                Phasefield = new Phasefield(null, levelSet.CGLevelSet, levelSet.DGLevelSet, levelSet.Tracker, extensionVelocity, m_grd, m_control, null);
                Phasefield.InitCH();
            }

            Phasefield.UpdateFields(levelSet.CGLevelSet, levelSet.DGLevelSet, levelSet.Tracker, extensionVelocity, m_grd, m_control, null);
            Phasefield.MovePhasefield(iTimestep, dt, time);
            iTimestep++;
        }

        double m_cahn;

        /// <summary>
        /// initialize the $tanh$ profile of the Phasefield from a signed distance level set
        /// </summary>
        /// <param name="levelSet"></param>
        private void InitializePhasefield(LevelSet levelSet) {

            SinglePhaseField phasefield = new SinglePhaseField(levelSet.Basis);

            // interface thickness = f(pOrder, D), WIP
            double dInterface = 2.0 * Math.Pow(2.0, this.m_grd.SpatialDimension) / levelSet.Basis.Degree;
            // dInterface * 1/4.164 * hmin
            double hmin;
            // set the interface width based on the largest CutCell
            // for now just a rough estimate, based on global grid size
            hmin = this.m_grd.iGeomCells.h_max.Max();
            m_cahn = dInterface * (1.0 / 4.164) * hmin / Math.Sqrt(2);

            // compute and project 
            // step one calculate distance field phiDist = 0.5 * log(Max(1+phi, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn_old
            // step two project the new phasefield phiNew = tanh(phiDist/(sqrt(2) * Cahn_new))
            // here done in one step, with default quadscheme
            // ===================
            phasefield.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Phi
                    // -----------------------------
                    levelSet.Evaluate(j0, Len, NS, result);

                    // compute the pointwise values of the new level set
                    // -----------------------------

                    result.ApplyAll(x => Math.Tanh(x/(Math.Sqrt(2) * m_cahn)));
                }
            );

            levelSet.Clear();
            levelSet.Acc(1.0, phasefield);

            potential = new SinglePhaseField(levelSet.Basis);
        }
    }
}
