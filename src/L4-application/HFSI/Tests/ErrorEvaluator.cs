using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using HFSISolver;
using BoSSS.Solution.Tecplot;

namespace HFSISolver.Tests
{

    /// <summary>
    /// error evaluator for ZLS-based tests 
    /// computes error for: velocity, pressure 
    /// </summary>
    /// <remarks>
    /// Seems redundant with <see cref="L2ErrorLogger"/>
    /// </remarks>
    public class ZLSErrorEvaluator<T> : BoSSS.Application.XNSE_Solver.Tests.XNSEErrorEvaluator<T> where T: HFSI_Control, new() {
        
        public ZLSErrorEvaluator(IApplication solver) : base(solver){

        }

        new protected ApplicationWithSolver<T> solver {
            get {
                return (ApplicationWithSolver<T>)(base.solver);
            }
        }

        public double ComputeTemperatureError(IDictionary<string, Func<double[], double, double>> exactTemperature, double time) {
            int D = solver.GridData.SpatialDimension;
            var FluidSpecies = exactTemperature.Keys.ToArray();

            int order = 0;
            if(solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            //order = Math.Max(temperature.Basis.Degree * 2, order);
            Console.WriteLine("temperature error with order " + order);

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            double L2Error = 0;
            Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

            foreach(var spc in FluidSpecies) {

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                string temperatureName = BoSSS.Solution.NSECommon.VariableNames.Temperature;
                ConventionalDGField temperature = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == temperatureName)).GetSpeciesShadowField(spc);
                var rule = scheme.Compile(solver.GridData, order);

                double IdV = temperature.LxError(exactTemperature[spc].Vectorize(time), (X, a, b) => (a - b).Pow2(), rule);
                L2Error += IdV;
                L2Error_Species.Add(spc, L2Error.Sqrt());

                solver.QueryHandler.ValueQuery("L2err_" + BoSSS.Solution.NSECommon.VariableNames.Temperature + "#" + spc, L2Error_Species[spc], true);
            }


            L2Error = L2Error.Sqrt();
            solver.QueryHandler.ValueQuery("L2err_" + BoSSS.Solution.NSECommon.VariableNames.Temperature, L2Error, true);
            return L2Error;
        }

        public double[] ComputeDisplacementError(IDictionary<string, Func<double[], double, double>[]> exactDisplacement, double time) {
            int D = solver.GridData.SpatialDimension;
            var FluidSpecies = exactDisplacement.Keys.ToArray();
            //var FluidSpecies = new[] { "C" };
            double[] Ret = new double[D];

            int order = 0;
            if(solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
            double[] L2Error = new double[D];

            foreach(var spc in FluidSpecies) {
                L2Error_Species.Add(spc, new double[D]);

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                for(int d = 0; d < D; d++) {
                    string displacement = HFSISolver.VariableNames.DisplacementVector(D)[d];
                    //string displacement = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d];
                    ConventionalDGField Dis_d = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == displacement)).GetSpeciesShadowField(spc);

                    double LxError = Dis_d.LxError(exactDisplacement[spc][d].Vectorize(time), null, scheme.SaveCompile(Dis_d.GridDat, order));
                    LxError = (LxError > -1e-12 && LxError < 0) ? 0.0 : LxError; // Avoid nans if solution is too close to the analytic one
                    L2Error_Species[spc][d] = LxError.Sqrt();

                    L2Error[d] += L2Error_Species[spc][d].Pow2();

                    solver.QueryHandler.ValueQuery("L2err_" + HFSISolver.VariableNames.DisplacementComponent(d) + "#" + spc, L2Error_Species[spc][d], true);
                }
            }
            L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

            for(int d = 0; d < D; d++) {
                solver.QueryHandler.ValueQuery("L2err_" + HFSISolver.VariableNames.DisplacementComponent(d), L2Error[d], true);
                Ret[d] = L2Error[d];
            }
            return L2Error;
        }


        IDictionary<string, Func<double[], double, double>[]> GetExactSolutionDisplacement(AppControl control) {
            int D = solver.GridData.SpatialDimension;
            return base.GetExactSolution(control, VariableNames.DisplacementVector(D));
        }

        IDictionary<string, Func<double[], double, double>> GetExactSolutionTemperature(AppControl control) {
            return base.GetExactSolution(control, BoSSS.Solution.NSECommon.VariableNames.Temperature);
        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="HFSI_Control.ExactSolutionVelocity"/> and <see cref="HFSI_Control.ExactSolutionPressure"/>).
        /// </summary>
        public double[] ComputeL2Error(double time, HFSI_Control control) {
            int D = solver.GridData.SpatialDimension;
           
            IDictionary<string, Func<double[], double, double>[]> ExactSolutionDisplacement = this.GetExactSolutionDisplacement(control);
            IDictionary<string, Func<double[], double, double>> ExactSolutionTemperature = this.GetExactSolutionTemperature(control);

            var Ret0 = base.ComputeL2Error(time, control);

            double[] Ret1 = new double[1], Ret2 = new double[D];
            if(ExactSolutionTemperature != null) {
                Ret1 = [ ComputeTemperatureError(ExactSolutionTemperature, time) ];
            }
            if(ExactSolutionDisplacement != null) {
                Ret2 = ComputeDisplacementError(ExactSolutionDisplacement, time);
            }
            return Ret0.Cat(Ret1, Ret2);
        }
    }

   
}
