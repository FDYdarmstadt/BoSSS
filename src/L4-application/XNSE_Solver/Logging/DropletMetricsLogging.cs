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

using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;
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
    /// <see cref="DropletMetricsLogging<T>"/>
    /// </summary>
    [Serializable]
    public class DropletMetricsLogging : DropletMetricsLogging<XNSE_Control> {
    }

    /// <summary>
    /// For single droplet/bubble computations with origin in (0,0,0):
    /// extraction of the major axis lengths (theta = 0 and theta = 90, where theta describes the inclination angle) and volume  
    /// </summary>
    [Serializable]
    public class DropletMetricsLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {


        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => "DropletMetrics";


        /// <summary>
        /// CSV header
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {

            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", "#timestep", "time", "theta0", "theta90x", "theta90y", "volume");
            Log.WriteLine(header);

        }


        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double physTime) {
            using(new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                double[] majorAxis = new double[3]; 
                switch (D) {
                    case 2: {
                        majorAxis = getMajorAxis2D();
                        break;
                    }
                    case 3: {
                        majorAxis = getMajorAxis3D();
                        break;
                    }
                    default:
                        throw new NotSupportedException("Not supported spatial dimension");
                }

                // volume(area)
                double volume = 0.0;
                int InnerSpecies = 0;
                SpeciesId spcId = LsTrk.SpeciesIdS[InnerSpecies];
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), this.m_HMForder, 1).XQuadSchemeHelper;
                var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    vqs.Compile(LsTrk.GridDat, this.m_HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++)
                            volume += ResultsOfIntegration[i, 0];
                    }
                ).Execute();


                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", TimestepNo, physTime, majorAxis[0], majorAxis[1], majorAxis[2], volume);
                Log.WriteLine(line);
                Log.Flush();
            }

        }

        /// <summary>
        /// implementation for 3D settings
        /// </summary>
        /// <returns></returns>
        private double[] getMajorAxis3D() {

            MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

            double r_theta0 = 0.0;
            double dist_zAxis_min = double.MaxValue;
            double r_theta90x = 0.0;
            double dist_xAxis_min = double.MaxValue;
            double r_theta90y = 0.0;
            double dist_yAxis_min = double.MaxValue;
            for(int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                double xCoord = InterfacePoints[i, 0];
                double yCoord = InterfacePoints[i, 1];
                double zCoord = InterfacePoints[i, 2];

                double dist_xAxis = Math.Sqrt(yCoord.Pow2() + zCoord.Pow2());
                if(dist_xAxis < dist_xAxis_min && xCoord > 0.0) {
                    dist_xAxis_min = dist_xAxis;
                    r_theta90x = xCoord;
                }
                double dist_yAxis = Math.Sqrt(xCoord.Pow2() + zCoord.Pow2());
                if(dist_yAxis < dist_yAxis_min && yCoord > 0.0) {
                    dist_yAxis_min = dist_yAxis;
                    r_theta90y = yCoord;
                }
                double dist_zAxis = Math.Sqrt(xCoord.Pow2() + yCoord.Pow2());
                if(dist_zAxis < dist_zAxis_min && zCoord > 0.0) {
                    dist_zAxis_min = dist_zAxis;
                    r_theta0 = zCoord;
                }
            }

            return new double[] { r_theta0, r_theta90x, r_theta90y };
        }

        /// <summary>
        /// implementation for 2D settings
        /// </summary>
        /// <returns></returns>
        private double[] getMajorAxis2D() {

            MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

            double r_theta0 = 0.0;
            double dist_zAxis_min = double.MaxValue;
            double r_theta90x = 0.0;
            double dist_xAxis_min = double.MaxValue;
            double r_theta90y = 0.0;
            double dist_yAxis_min = double.MaxValue;
            for(int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                double xCoord = InterfacePoints[i, 0];
                double yCoord = InterfacePoints[i, 1];
                double zCoord = InterfacePoints[i, 2];

                double dist_xAxis = Math.Sqrt(yCoord.Pow2() + zCoord.Pow2());
                if(dist_xAxis < dist_xAxis_min && xCoord > 0.0) {
                    dist_xAxis_min = dist_xAxis;
                    r_theta90x = xCoord;
                }
                double dist_yAxis = Math.Sqrt(xCoord.Pow2() + zCoord.Pow2());
                if(dist_yAxis < dist_yAxis_min && yCoord > 0.0) {
                    dist_yAxis_min = dist_yAxis;
                    r_theta0 = yCoord;
                }
            }

            return new double[] { r_theta0, r_theta90x, 0.0};
        }

    }

}
