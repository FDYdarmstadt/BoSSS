/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSEC {

    static public partial class NUnitTest {
        private class CombustionErrorEvaluator<T> : ErrorEvaluator where T : XNSEC_Control, new() {

            public CombustionErrorEvaluator(IApplication solver) : base(solver) {
            }

            new protected Solution.LevelSetTools.SolverWithLevelSetUpdater.SolverWithLevelSetUpdater<T> solver {
                get {
                    return (Solution.LevelSetTools.SolverWithLevelSetUpdater.SolverWithLevelSetUpdater<T>)(base.solver);
                }
            }

            /// <summary>
            /// Computes the L2 Error of the actual solution against the exact solution in the control object
            /// (<see cref="XNSE_Control.ExactSolutionTemperature"/>.
            /// </summary>
            public double ComputeTemperatureError(IDictionary<string, Func<double[], double, double>> exactTemperature, double time) {
                int D = solver.GridData.SpatialDimension;

                int order = 0;
                if(solver.LsTrk.GetCachedOrders().Count > 0) {
                    order = solver.LsTrk.GetCachedOrders().Max();
                } else {
                    order = 1;
                }

                var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach(var spc in solver.LsTrk.SpeciesNames) {
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

            /// <summary>
            /// Computes the L2 Error of the actual solution against the exact solution in the control object
            /// (<see cref="XNSE_Control.ExactSolutionTemperature"/>.
            /// </summary>
            public double ComputeMassFractionError(IDictionary<string, Func<double[], double, double>[]> exactTemperature, double time, int y) {
                int D = solver.GridData.SpatialDimension;

                int order = 0;
                if(solver.LsTrk.GetCachedOrders().Count > 0) {
                    order = solver.LsTrk.GetCachedOrders().Max();
                } else {
                    order = 1;
                }

                var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach(var spc in solver.LsTrk.SpeciesNames) {
                    SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                    string variablename = VariableNames.MassFraction_n(y);
                    ConventionalDGField massFraction_y = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == variablename)).GetSpeciesShadowField(spc);
                    var rule = scheme.Compile(solver.GridData, order);

                    double IdV = massFraction_y.LxError(exactTemperature[spc][y].Vectorize(time), (X, a, b) => (a - b).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, L2Error.Sqrt());

                    solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Temperature + "#" + spc, L2Error_Species[spc], true);
                }

                L2Error = L2Error.Sqrt();
                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Temperature, L2Error, true);
                return L2Error;
            }

            public override double[] ComputeL2Error(double time, XNSE_Control control) {
                return ComputeL2Error(time, (XNSEC_Control)control); ;
            }

            public double[] ComputeL2Error(double time, XNSEC_Control control) {
                int NoOfSpcs = control.NumberOfChemicalSpecies;
                double[] Ret = new double[1 + NoOfSpcs]; // temperature and all mass fractions

                if(control.ExactSolutionTemperature != null) {
                    double error = ComputeTemperatureError(control.ExactSolutionTemperature, time);
                    Ret[0] = error;
                }

                if(control.ExactSolutionMassFractions != null) {
                    for(int y = 0; y < NoOfSpcs; y++) {
                        double error = ComputeMassFractionError(control.ExactSolutionMassFractions, time, y);
                        Ret[y + 1] = error;
                    }
                }

                return Ret;
            }
        }

        private class MixtureFractionErrorEvaluator<T> : ErrorEvaluator where T : XNSEC_Control, new() {

            public MixtureFractionErrorEvaluator(IApplication solver) : base(solver) {
            }

            new protected Solution.LevelSetTools.SolverWithLevelSetUpdater.SolverWithLevelSetUpdater<T> solver {
                get {
                    return (Solution.LevelSetTools.SolverWithLevelSetUpdater.SolverWithLevelSetUpdater<T>)(base.solver);
                }
            }

            /// <summary>
            /// Computes the L2 Error of the actual solution against the exact solution in the control object
            /// </summary>
            public double ComputeMixtureFractionError(IDictionary<string, Func<double[], double, double>> exactMixtureFraction, double time) {
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

                    string MF_Field = VariableNames.MixtureFraction;
                    ConventionalDGField temperature = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == MF_Field)).GetSpeciesShadowField(spc);
                    var rule = scheme.Compile(solver.GridData, order);

                    double IdV = temperature.LxError(exactMixtureFraction[spc].Vectorize(time), (X, a, b) => (a - b).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, L2Error.Sqrt());

                    solver.QueryHandler.ValueQuery("L2err_" + VariableNames.MixtureFraction + "#" + spc, L2Error_Species[spc], true);
                }

                L2Error = L2Error.Sqrt();
                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.MixtureFraction, L2Error, true);
                return L2Error;
            }

         
            public override double[] ComputeL2Error(double time, XNSE_Control control) {
                return ComputeL2Error(time, (XNSEC_Control)control); ;
            }

            public double[] ComputeL2Error(double time, XNSEC_Control control) {
                
                double[] Ret = new double[1]; // MixtureFraction field
                if (control.ExactSolutionMixtureFraction!= null) {
                    double error = ComputeMixtureFractionError(control.ExactSolutionMixtureFraction, time);
                    Ret[0] = error;
                }

  

                return Ret;
            }
        }

    }
}