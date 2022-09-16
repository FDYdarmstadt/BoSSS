using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
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
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Level-set evolution using a div-free prescribed velocity field
    /// </summary>
    /// <remarks>
    /// - implemented by Mr, jun22
    /// </remarks>
    public class PrescribedVelocityEvolver : ILevelSetEvolver {


        /// <summary>
        /// ctor
        /// </summary>
        public PrescribedVelocityEvolver(string levelSetName, int D ) {
            this.SpatialDimension = D;
            this.levelSetName = levelSetName;
            parameters = NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
        }

        int SpatialDimension;
        string levelSetName;
        string[] parameters;

        /// <summary>
        /// should only be the interface velocity vector; typically, a phase-averaged velocity.
        /// </summary>
        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// nix
        /// </summary>
        public IList<string> VariableNames => null;

        // nothing to do
        public Func<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>, bool> AfterMovePhaseInterface => null;


        /// <summary>
        /// 
        /// <see cref="ILevelSetEvolver.InternalFields"/>
        /// </summary>
        public IDictionary<string, DGField> InternalFields { 
            get {                
                return null;
            }
        }

        ITimeStepper timeStepper;

        /// <summary>
        /// Flux for a scalar transport equation
        /// </summary>
        internal class UpwindFlux : LinearFlux {

            public UpwindFlux(int D) {
                m_spatDim = D;
            }

            int m_spatDim;

            protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
                Vector n = inp.Normal;

                var vel = ((Vector)inp.Parameters_IN);

                if (n * vel >= 0) {
                    // flow from inside 
                    return (vel * Uin[0]) * n;
                } else {
                    // flow from outside into the domain
                    return (vel * Uin[0]) * n;
                }                
            }

            /// <summary>
            /// An upwind flux
            /// </summary>
            protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
                Vector n = inp.Normal;

                var vel = 0.5 * (((Vector)inp.Parameters_IN) + ((Vector)inp.Parameters_OUT));

                if (vel * n > 0)
                    return (vel * Uin[0]) * n;
                else
                    return (vel * Uout[0]) * n;

            }

            /// <summary>
            /// `$ \vec{u} \varphi `$
            /// </summary>
            protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
                for (int d = 0; d < inp.D; d++) {
                    output[d] = inp.Parameters[d] * U[0];
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public override IList<string> ArgumentOrdering {
                get { return new string[] { "Phi" }; }
            }

            /// <summary>
            /// the transport velocity
            /// </summary>
            public override IList<string> ParameterOrdering {
                get {
                    return Solution.NSECommon.VariableNames.VelocityVector(this.m_spatDim);
                }
            }
        }

        private ITimeStepper InitializeRungeKutta(SinglePhaseField levelSet, SinglePhaseField[] Velocity) {
            var diffOp = new SpatialOperator(new string[] { "Phi" },
                Solution.NSECommon.VariableNames.VelocityVector(this.SpatialDimension),
                new string[] { "codom1" },
                QuadOrderFunc.NonLinear(1));
            diffOp.EquationComponents["codom1"].Add(new UpwindFlux(this.SpatialDimension));
            diffOp.TemporalOperator = new ConstantTemporalOperator(diffOp, 1);
            diffOp.IsLinear = true;
            diffOp.Commit();

            var Timestepper = new RungeKutta(RungeKuttaScheme.TVD3, diffOp, levelSet.Mapping, new CoordinateMapping(Velocity));
            return Timestepper;
        }

        
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

                if(timeStepper == null) {
                    //timeStepper = InitializeAdamsBashforth(levelSet.DGLevelSet, extensionVelocity);
                    timeStepper = InitializeRungeKutta(levelSet.DGLevelSet, meanVelocity);
                }
                if(!ReferenceEquals(timeStepper.Mapping.Fields[0], levelSet.DGLevelSet)) {
                    throw new Exception("Something went wrong with the internal pointer magic of the levelSetTracker. Definitely a weakness of ObjectOrientation.");
                }

                timeStepper.Perform(dt);

                tr.Info("time in LS evolver: " + timeStepper.Time );
            }
        }
        
    }
}
