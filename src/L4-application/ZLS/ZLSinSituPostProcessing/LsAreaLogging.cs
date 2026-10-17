using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ZLSinSituPostProcessing
{
    public class LsAreaLogging : ZLSinSituPostProcessingModule
    {

        public LsAreaLogging()
        {

        }
        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "BenchmarkQuantities_LsAreaLogging";

        /// <summary>
        /// hard-coded name for the ConvergenceTest_ZLS
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// Header for the ConvergenceTest_ZLS log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter)
        {
            string header = String.Format("{0}\t{1}\t{2}", "#timestep", "phystime","area");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime)
        {
            using (new FuncTrace())
            {
                double R = ComputeBenchmarkQuantities_area();

                AppendToLog(TimestepNo);
                AppendToLog(phystime);
                AppendToLog(R);

                base.QueryResultTable.LogValue("area", R);

            }
        }
        internal double ComputeBenchmarkQuantities_area()
        {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0)
            {
                order = LsTrk.GetCachedOrders().Max();
            }
            else
            {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // area
            double area = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[1];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return area;

        }

    }
}
