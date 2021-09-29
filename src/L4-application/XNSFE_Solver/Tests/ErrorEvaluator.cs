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
using BoSSS.Application.XNSE_Solver.Tests;

namespace BoSSS.Application.XNSFE_Solver.Tests
{

    class XHeatErrorEvaluator<T> : ErrorEvaluator where T: XNSE_Control, new() {

        public XHeatErrorEvaluator(IApplication solver) : base(solver) {

        }
        new protected SolverWithLevelSetUpdater<T> solver {
            get {
                return (SolverWithLevelSetUpdater<T>)(base.solver);
            }
        }


        public double ComputeTemperatureError(IDictionary<string, Func<double[], double, double>> exactTemperature, double time) {
            int D = solver.GridData.SpatialDimension;

            int order = 0;
            if (solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            double L2Error = 0;
            Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

            foreach (var spc in solver.LsTrk.SpeciesNames) {

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                string temperatureName = VariableNames.Temperature;
                ConventionalDGField temperature = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == temperatureName)).GetSpeciesShadowField(spc);
                var rule = scheme.Compile(solver.GridData, order);

                double IdV = temperature.LxError(exactTemperature[spc].Vectorize(time), (X, a, b) => (a - b).Pow2(), rule);
                L2Error += IdV;
                L2Error_Species.Add(spc, L2Error.Sqrt());

                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Temperature + "#" + spc, L2Error_Species[spc], true);
            }


            L2Error = L2Error.Sqrt();
            solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Temperature, L2Error, true);
            return L2Error;            
        }

        public double ComputeEnergyError(Func<double, double> exactEnergy, double time) {
            int D = solver.GridData.SpatialDimension;

            int order = 0;
            if (solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            double TotalEnergy = 0;

            foreach (var spc in solver.LsTrk.SpeciesNames) {

                double c, rho;
                switch (spc) {
                    case "A": { c = solver.Control.ThermalParameters.c_A; rho = solver.Control.ThermalParameters.rho_A; break; }
                    case "B": { c = solver.Control.ThermalParameters.c_B; rho = solver.Control.ThermalParameters.rho_B; break; }
                    default: { throw new ArgumentException(); }
                }

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                string temperatureName = VariableNames.Temperature;
                ConventionalDGField temperature = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == temperatureName)).GetSpeciesShadowField(spc);
                var rule = scheme.Compile(solver.GridData, order);

                double E = 0.0;
                CellQuadrature.GetQuadrature(new int[] { 1 }, solver.LsTrk.GridDat, rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    temperature.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1,-1,0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        E += c * rho * ResultsOfIntegration[i, 0];
                }).Execute();
                TotalEnergy += E;
            }

            double EnergyError = (exactEnergy(time) - TotalEnergy).Abs() / exactEnergy(time).Abs(); // relative Error in Energy, scaled by some factor
            return EnergyError;
        }

        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionTemperature"/>.
        /// </summary>
        public override double[] ComputeL2Error(double time, XNSE_Control control) {
            double[] Ret = new double[1];

            if (control.ExactSolutionTemperature != null) {
                double error = ComputeTemperatureError(control.ExactSolutionTemperature, time);
                Ret[0] = error;
            }
            
            return Ret;
        }
    }
}
