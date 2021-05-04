using BoSSS.Application.XNSE_Solver.Legacy;
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {


    /// <summary>
    /// Post-processing specific to <see cref="RisingBubble"/>
    /// </summary>
    [Serializable]
    public class MassfluxLogging : XNSEinSituPostProcessingModule {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "Massflux";

        /// <summary>
        /// hard-coded name for the Rising bubble
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// Header for the rising bubble log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "#timestep", "time", "mass-liq", "mass-vap", "mass-total", "masschange-evap", "masschange-vapor", "masschange-liquid", "masschange-total", "interface length");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using (new FuncTrace()) {
                var R = ComputeBenchmarkQuantities_Massflux();
                AppendToLog(iTimestep);
                AppendToLog(PhysTime);
                AppendToLog(R);
            }
        }

        ThermalParameters m_thermParams;
        ThermalParameters ThermParams {
            get {
                if(m_thermParams == null) {
                    m_thermParams = this.Control.ThermalParameters;
                }
                return m_thermParams;
            }
        }
        (double mass_liq, double mass_vap, double mass_evap, double mass_total, double interface_length) LastState;
        internal (double mass_liq, double mass_vap, double mass_total, double mass_evap, double mass_evap_v, double mass_evap_l, double mass_evap_t, double interface_length) ComputeBenchmarkQuantities_Massflux() {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // mass of liquid
            double mass_liq = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        mass_liq += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            mass_liq *= ThermParams.rho_A;

            // mass of vapor
            double mass_vap = 0.0;
            spcId = LsTrk.SpeciesIdS[1];
            vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        mass_vap += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            mass_vap *= m_thermParams.rho_B;           

            // evaporated/condensed mass
            double mass_evap = 0.0;
            var lqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            if (this.ThermParams.hVap != 0.0) {                
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    lqs.Compile(LsTrk.GridDat, order),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        CurrentMassFlux.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            mass_evap += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
                double dt = this.Control.GetFixedTimestep();
                mass_evap *= dt;
            }

            // total mass
            double mass_total = mass_liq + mass_vap;

            // interface length
            double interface_length = 0.0;
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                lqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.AccConstant(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        interface_length += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            // calculate differences
            double mass_evap_v = mass_vap - LastState.mass_vap;
            double mass_evap_l = mass_liq - LastState.mass_liq;
            double mass_evap_t = mass_total - LastState.mass_total;

            // save current state
            LastState = (mass_liq, mass_vap, mass_evap, mass_total, interface_length);

            return (mass_liq, mass_vap, mass_total, mass_evap, mass_evap_v, mass_evap_l, mass_evap_t, interface_length);
        }

        /// <summary>
        /// current massflux parameter
        /// </summary>
        protected DGField CurrentMassFlux {
            get {
                if (this.SolverMain is XNSE_SolverMain oldSolver) {
                    throw new NotImplementedException();
                } else if (this.SolverMain is XNSFE newSolver) {
                    int D = this.SolverMain.GridData.SpatialDimension;
                    IReadOnlyDictionary<string, DGField> parameters = newSolver.LsUpdater.Parameters;

                    DGField ret = null;
                    for (int i = 0; i < 3; ++i) {
                        if (parameters.TryGetValue(VariableNames.MassFluxExtension, out DGField velocityField)) {
                            ret = velocityField;
                        } else {
                            throw new ApplicationException("Unable to identify mass flux extension field.");
                        }
                    }

                    return ret;
                } else {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// current ls velocity y parameter
        /// </summary>
        protected DGField CurrentVelocityYLevelSet {
            get {
                if (this.SolverMain is XNSE_SolverMain oldSolver) {
                    throw new NotImplementedException();
                } else if (this.SolverMain is XNSFE newSolver) {
                    int D = this.SolverMain.GridData.SpatialDimension;
                    IReadOnlyDictionary<string, DGField> parameters = newSolver.LsUpdater.Parameters;

                    DGField ret = null;
                    for (int i = 0; i < 3; ++i) {
                        if (parameters.TryGetValue(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCG, VariableNames.Velocity_d(1)), out DGField velocityField)) {
                            ret = velocityField;
                        } else {
                            throw new ApplicationException("Unable to identify level set velocity y field.");
                        }
                    }

                    return ret;
                } else {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// current temperature solution
        /// </summary>
        protected XDGField CurrentTemperature {
            get {
                if (this.SolverMain is XNSE_SolverMain oldSolver) {
                    throw new NotImplementedException();
                } else if (this.SolverMain is XNSE newSolver) {
                    int D = this.SolverMain.GridData.SpatialDimension;

                    var ret = newSolver.CurrentState.Fields.ElementAt(D+1) as XDGField;
                    
                    if (ret.Identification != VariableNames.Temperature)
                        throw new ApplicationException("Unable to identify temperature field.");
                    

                    return ret;
                } else {
                    throw new NotImplementedException();
                }
            }
        }
    }    
}
