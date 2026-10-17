using System;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using NUnit.Framework;
using BoSSS.Solution.XdgTimestepping;
using MPI.Wrappers;

namespace BoSSS.Application.TraceDGtest {


    public class TraceDGtestMain : XdgApplicationWithSolver<TraceDGtestControl> {
        static void Main(string[] args) {
            //InitMPI(args);
            //ilPSP.Environment.NumThreads = 1;
            //BoSSS.Application.TraceDGtest.UnitTests.CartesianWithLevelSetAMR(2, true);
            //FinalizeMPI();
            //Assert.IsTrue(false, "remove me");

            _Main(args, false, () => new TraceDGtestMain());
        }

       
        protected override XDifferentialOperatorMk2 GetOperatorInstance(int D) {
            // Define domain and codomain variable names
            var domainVars = new[] { "SurfaceConcentration" };
            var parameters = new[] { "RHSconcentration" };
            var codomainVars = new[] { "SurfaceBalance" };

            // Use species names from the tracker
            var species = _LsTrk.SpeciesNames;

            // Create the operator
            var op = new XDifferentialOperatorMk2(
                domainVars,
                parameters,
                codomainVars,
                QuadOrderFunc.MaxDegTimesTwo(),
                species
            );
            op.AgglomerationThreshold = this.Control.AgglomerationThreshold;


            // Add a simple identity flux for each species (replace with your physics as needed)
            op.SurfaceElementOperator_Ls0.EquationComponents["SurfaceBalance"].Add(
                new IdentitityOperator("SurfaceConcentration", "RHSconcentration", () => this.RHS));
            //op.EquationComponents["SurfaceBalance"].Add(
            //                new IdentitityOperator("SurfaceConcentration", "RHS"));

            op.Commit();
            return op;

        }

       

#pragma warning disable 649, 8618, 169
        /// <summary>
        /// 
        /// </summary>
        [LevelSetTracker("-:A +:B", 1)]
        LevelSetTracker _LsTrk;

        /// <summary>
        /// continuous projection of <see cref="PhiDG"/>
        /// </summary>
        [InstantiateFromControlFile("Phi", null, IOListOption.Always)]
        LevelSet Phi;

        // <summary>
        // the discontinuous level-set
        // </summary>
        //[InstantiateFromControlFile("PhiDG", "VelocityX", IOListOption.Always)]
        //SinglePhaseField PhiDG;


        [InstantiateFromControlFile(["VelocityX", "VelocityY", "VelocityZ"],
            ioListOpt: IOListOption.Always)]
        public VectorField<XDGField> Velocity;


        [InstantiateFromControlFile(FieldIdentification: "SurfaceConcentration", DGDegreeImpliedBy: "Velocity*", ioListOpt: IOListOption.Always)]
        public TraceDGField SurfaceConcentration;
        
        [InstantiateFromControlFile(FieldIdentification: "RHS", DGDegreeImpliedBy: "Velocity*", ioListOpt: IOListOption.Always)]
        public TraceDGField RHS;

        protected override LevelSetHandling LevelSetHandling => LevelSetHandling.Coupled_Once;
#pragma warning restore 649, 8618, 169


        protected override IEnumerable<DGField> InstantiateSolutionFields() {
            return [SurfaceConcentration];
        }
    }

    class IdentitityOperator : IVolumeForm, IParameterHandling {

        public IdentitityOperator(string _DomVer, string _RHSparam, Func<DGField> _GetRHS) {
            m_DomVer = _DomVer;
            m_RHSparam = _RHSparam;
            m_GetRHS = _GetRHS;
        }

        Func<DGField> m_GetRHS;
        string m_DomVer;
        string m_RHSparam;

        public TermActivationFlags VolTerms => TermActivationFlags.V | TermActivationFlags.UxV;

        public IList<string> ArgumentOrdering => [m_DomVer];

        public IList<string> ParameterOrdering => [m_RHSparam];

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return U[0] * V - cpv.Parameters[0] * V;
        }

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var rhsDG = m_GetRHS();
            if(!object.ReferenceEquals(Parameters[0], rhsDG)) {
                Parameters[0].Clear();
                Parameters[0].AccLaidBack(1.0, rhsDG);
            }
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var b0 = Arguments[0].Basis;
            return [ new SinglePhaseField(new Basis(b0.GridDat, b0.Degree)) ];
        }
    }

}
