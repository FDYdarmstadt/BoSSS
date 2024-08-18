using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Errors against exact solution (<see cref="XNSE_Control.ExactSolutionVelocity"/>, <see cref="XNSE_Control.ExactSolutionPressure"/>)
    /// </summary>
    [Serializable]
    public class L2ErrorLogger: L2ErrorLogger<XNSE_Control> { }

    /// <summary>
    /// Errors against exact solution (<see cref="XNSE_Control.ExactSolutionVelocity"/>, <see cref="XNSE_Control.ExactSolutionPressure"/>)
    /// </summary>
    [Serializable]
    public class L2ErrorLogger<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {
        
        /// <summary>
        /// Null: no log-file will be created; only queries will be saved
        /// </summary>
        protected override string LogFileName => null;

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            ComputeL2Error(PhysTime);
        }

        public override void Setup(IApplication solverMain) {
            base.Setup(solverMain);
        }

        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
        internal double[] ComputeL2Error(double time) {
            int D = this.SolverMainOverride.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            string[] fluidSpecies;

            if (this.Control.ExactSolutionVelocity == null && this.Control.ExactSolutionPressure == null)
                // nothing to do
                return Ret;

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


            // Velocity error
            // ==============
            if (this.Control.ExactSolutionVelocity != null) {
                fluidSpecies = this.Control.ExactSolutionVelocity.Keys.ToArray();

                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                foreach (var spc in fluidSpecies) {
                    L2Error_Species.Add(spc, new double[D]);

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                    for (int d = 0; d < D; d++) {
                        ConventionalDGField Vel_d = this.CurrentVel[d].GetSpeciesShadowField(spc);

                        L2Error_Species[spc][d] = Vel_d.L2Error(this.Control.ExactSolutionVelocity[spc][d].Vectorize(time), order, scheme);
                        L2Error[d] += L2Error_Species[spc][d].Pow2();

                        base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
                    }
                }
                L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

                for (int d = 0; d < D; d++) {
                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                    Ret[d] = L2Error[d];
                }

            }


            // pressure error
            // ==============
            if (this.Control.ExactSolutionPressure != null) {
                fluidSpecies = this.Control.ExactSolutionPressure.Keys.ToArray();

                // pass 1: mean value of pressure difference
                double DiffInt = 0, Volume = 0;
                foreach (var spc in fluidSpecies) {

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.SolverMainOverride.GridData, order);

                    DiffInt += this.CurrentPressure.GetSpeciesShadowField(spc).LxError(this.Control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b), rule);
                    Volume +=  this.CurrentPressure.GetSpeciesShadowField(spc).LxError(null, (X, a, b) => (1.0), rule); // exact volume of the 
                }
                double PressureDiffMean = DiffInt / Volume;


                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach (var spc in fluidSpecies) {

                    SpeciesId spId = this.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(this.GridData, order);

                    double IdV = this.CurrentPressure.GetSpeciesShadowField(spc).LxError(this.Control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, IdV.Sqrt());

                    base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
                }


                L2Error = L2Error.Sqrt();
                base.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
                Ret[D] = L2Error;

            } //*/

            return Ret;
        }
    }
}
