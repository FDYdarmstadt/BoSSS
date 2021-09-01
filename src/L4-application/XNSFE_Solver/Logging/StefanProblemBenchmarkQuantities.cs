using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases {


    /// <summary>
    /// Post-processing specific to <see cref="RisingBubble"/>
    /// </summary>
    [Serializable]
    public class StefanProblemBenchmarkQuantities : XNSFEinSituPostProcessingModule {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "StefanProblem";

        /// <summary>
        /// hard-coded name for the StefanProblem
        /// </summary>
        protected override string LogFileName => LogfileName;


        private byte OutletEdge;
        public StefanProblemBenchmarkQuantities(byte OutletTag) {
            OutletEdge = OutletTag;
        }

        /// <summary>
        /// Header for the rising bubble log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "#timestep", "time", "interface-x-pos-min", "interface-x-pos-max", "mass-vapor", "mass-liquid", "massflux-interface", "massflux-outlet");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using (new FuncTrace()) {
                var R = ComputeBenchmarkQuantities_Massflux();
                ITuple RT = R;
                double[] RR = new double[RT.Length];
                for (int i = 0; i < RT.Length; i++)
                    RR[i] = (double)RT[i];

                AppendToLog(iTimestep);
                AppendToLog(PhysTime);
                AppendToLog(RR);

                base.QueryResultTable.LogValue("xPos", (R.interface_xpos_max + R.interface_xpos_min)/2); // log mean interface position
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
        internal (double interface_xpos_min, double interface_xpos_max, double mass_vapor, double mass_liquid, double mass_flux_interface, double mass_flux_outlet) ComputeBenchmarkQuantities_Massflux() {

            int D = SolverMainOverride.Grid.SpatialDimension;

            double mass_flux_interface = 0.0;
            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


            double interface_xpos_min = double.MaxValue;
            double interface_xpos_max = double.MinValue;

            MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
            for (int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                if (InterfacePoints[i, 0] < interface_xpos_min) {
                    interface_xpos_min = InterfacePoints[i, 0];
                }
                if (InterfacePoints[i, 0] > interface_xpos_max) {
                    interface_xpos_max = InterfacePoints[i, 0];
                }
            }
           

            // mass of liquid
            double mass_liquid = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        mass_liquid += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            mass_liquid *= ThermParams.rho_A;

            // mass of vapor
            double mass_vapor = 0.0;
            spcId = LsTrk.SpeciesIdS[1];
            vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        mass_vapor += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            mass_vapor *= m_thermParams.rho_B;
            
            var lqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                lqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    CurrentMassFlux.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        mass_flux_interface += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            double dt = this.Control.GetFixedTimestep();
            mass_flux_interface *= dt;


            double mass_flux_outlet = 0.0;
            var em = new EdgeMask(this.GridData, this.GridData.EdgeTagNames[OutletEdge]);
            EdgeQuadratureScheme eqs = SchemeHelper.GetEdgeQuadScheme(LsTrk.GetSpeciesId("A"), IntegrationDomain: em);
            var normals = this.GridData.Edges.NormalsForAffine;
            EdgeQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                eqs.Compile(LsTrk.GridDat, order), //  agg.HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray DummyArray = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                    for (int d = 0; d < D; d++) {
                        CurrentVel[d].GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d), DummyArray.ExtractSubArrayShallow(-1, -1, d), null, null, null, null, 0, 0.0);
                    }                    
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            mass_flux_outlet += ResultsOfIntegration[i, d] * normals[i0 + i, d];
                        }
                    }
                }
            ).Execute();
            mass_flux_outlet *= ThermParams.rho_A * dt;

            return (interface_xpos_min, interface_xpos_max, mass_vapor, mass_liquid, mass_flux_interface, mass_flux_outlet);
        }
    }    
}
