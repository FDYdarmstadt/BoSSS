using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    public class Dropletlike : Dropletlike<XNSE_Control> { }

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    public class Dropletlike<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// SemiAxis
        /// </summary>
        public const string LogfileName = "SemiAxis";


        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => LogfileName;

        private int InnerSpecies;
        public Dropletlike(int InnerSpecies) {
            this.InnerSpecies = InnerSpecies;
        }
        public Dropletlike() : this(0) {
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", "#timestep", "time", "semi axis x", "semi axis y", "area", "perimeter");
            textWriter.WriteLine(header);

        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {

            MultidimensionalArray interP = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
            int numP = interP.Lengths[0];
            double[] Xcoord = new double[numP];
            double[] Ycoord = new double[numP];
            for (int i = 0; i < numP; i++) {
                Xcoord[i] = interP[i, 0];
                Ycoord[i] = interP[i, 1];
            }
            double semiAxisX = Xcoord.Max() - Xcoord.Min();
            double semiAxisY = Ycoord.Max() - Ycoord.Min();

            double[] sphereProps = this.ComputeSphericalPorperties();

            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", TimestepNo, phystime, semiAxisX, semiAxisY, sphereProps[0], sphereProps[1]);
            Log.WriteLine(line);
            Log.Flush();

            return;
        }



        double[] ComputeSphericalPorperties() {

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), this.m_HMForder, 1).XQuadSchemeHelper;

            // area/volume
            double volume = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[InnerSpecies];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        volume += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            // surface
            double surface = 0.0;
            //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId, 0);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                surfElemVol.Compile(LsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return new double[] { volume, surface };

        }
    }    
}
