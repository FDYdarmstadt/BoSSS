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
using System.Diagnostics;
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
        /// true if the droplet is axis symmetric
        /// </summary>
        public bool AxisSymmetric = false;


        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => "DropletMetrics";


        /// <summary>
        /// CSV header
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {

            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "northPole", "L", "Wx", "Wy", "volume");
            Log.WriteLine(header);

        }


        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double physTime) {
            using(new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                double[] majorAxis = new double[4]; 
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


                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, physTime, majorAxis[0], majorAxis[1], majorAxis[2], majorAxis[3], volume);
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

            // droplet width W along x-Axis
            double r_theta90x1 = 0.0;
            double dist_xAxis_min1 = double.MaxValue;
            double r_theta90x2 = 0.0;
            double dist_xAxis_min2 = double.MaxValue;
            // droplet width W along y-Axis
            double r_theta90y1 = 0.0;
            double dist_yAxis_min1 = double.MaxValue;
            double r_theta90y2 = 0.0;
            double dist_yAxis_min2 = double.MaxValue;
            // droplet length L along z-Axis 
            double r_theta0 = 0.0;
            double dist_zAxis_min1 = double.MaxValue;
            double r_theta180 = 0.0;
            double dist_zAxis_min2 = double.MaxValue;

            for(int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                double xCoord = InterfacePoints[i, 0];
                double yCoord = InterfacePoints[i, 1];
                double zCoord = InterfacePoints[i, 2];

                // droplet width W along x-Axis
                double dist_xAxis = Math.Sqrt(yCoord.Pow2() + zCoord.Pow2());
                if(dist_xAxis < dist_xAxis_min1 && xCoord > 0.0) {
                    dist_xAxis_min1 = dist_xAxis;
                    r_theta90x1 = xCoord;
                }
                if(!AxisSymmetric && dist_xAxis < dist_xAxis_min2 && xCoord < 0.0) {
                    dist_xAxis_min2 = dist_xAxis;
                    r_theta90x2 = xCoord;
                }


                // droplet width W along y-Axis
                double dist_yAxis = Math.Sqrt(xCoord.Pow2() + zCoord.Pow2());
                if(dist_yAxis < dist_yAxis_min1 && yCoord > 0.0) {
                    dist_yAxis_min1 = dist_yAxis;
                    r_theta90y1 = yCoord;
                }
                if(!AxisSymmetric && dist_yAxis < dist_yAxis_min2 && yCoord < 0.0) {
                    dist_yAxis_min2 = dist_yAxis;
                    r_theta90y2 = yCoord;
                }

                // droplet length L along z-Axis 
                double dist_zAxis = Math.Sqrt(xCoord.Pow2() + yCoord.Pow2());
                if(dist_zAxis < dist_zAxis_min1 && zCoord > 0.0) {
                    dist_zAxis_min1 = dist_zAxis;
                    r_theta0 = zCoord;
                }
                if(dist_zAxis < dist_zAxis_min2 && zCoord < 0.0) {
                    dist_zAxis_min2 = dist_zAxis;
                    r_theta180 = zCoord;
                }  
                
            }

            double Wx = r_theta90x1 - r_theta90x2;
            if(Wx < 0)
                throw new ArgumentOutOfRangeException("Droplet length Wx is negative");

            double Wy = r_theta90y1 - r_theta90y2;
            if(Wy < 0)
                throw new ArgumentOutOfRangeException("Droplet length Wy is negative");

            double L = r_theta0 - r_theta180;
            if(L < 0)
                throw new ArgumentOutOfRangeException("Droplet length L is negative");

            return new double[] { r_theta0, L, Wx, Wy };
        }

        /// <summary>
        /// implementation for 2D settings
        /// </summary>
        /// <returns></returns>
        private double[] getMajorAxis2D() {

            MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);

            double r_theta0 = 0.0;
            double r_theta180 = 0.0;
            double r_theta90x1 = 0.0;
            double r_theta90x2 = 0.0;
            double dist_xAxis_min1 = double.MaxValue;
            double dist_xAxis_min2 = double.MaxValue;
            double dist_yAxis_min1 = double.MaxValue;
            double dist_yAxis_min2 = double.MaxValue;
            for(int i = 0; i < InterfacePoints.Lengths[0]; i++) {
                double xCoord = InterfacePoints[i, 0];
                double yCoord = InterfacePoints[i, 1];
                double zCoord = InterfacePoints[i, 2];

                double dist_xAxis = Math.Sqrt(yCoord.Pow2() + zCoord.Pow2());
                if(dist_xAxis < dist_xAxis_min1 && xCoord > r_theta90x2) {
                    dist_xAxis_min1 = dist_xAxis;
                    r_theta90x1 = xCoord;
                }
                if(dist_xAxis < dist_xAxis_min2 && xCoord < r_theta90x1) {
                    dist_xAxis_min2 = dist_xAxis;
                    r_theta90x2 = xCoord;
                }

                double dist_yAxis = Math.Sqrt(xCoord.Pow2() + zCoord.Pow2());
                if(dist_yAxis < dist_yAxis_min1 && yCoord > r_theta180) {
                    dist_yAxis_min1 = dist_yAxis;
                    r_theta0 = yCoord;
                }
                if(dist_yAxis < dist_yAxis_min2 && yCoord < r_theta0) {
                    dist_yAxis_min2 = dist_yAxis;
                    r_theta180 = yCoord;
                }
            }

            double W = r_theta90x1 - r_theta90x2;
            if(W < 0)
                throw new ArgumentOutOfRangeException("Droplet length W is negative");

            double L = r_theta0 - r_theta180;
            if(L < 0)
                throw new ArgumentOutOfRangeException("Droplet length L is negative");

            return new double[] { r_theta0, L, W, 0.0};
        }

    }

}
