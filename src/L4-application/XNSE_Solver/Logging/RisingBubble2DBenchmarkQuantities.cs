using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Post-processing specific to <see cref="RisingBubble"/>
    /// </summary>
    [Serializable]
    public class RisingBubble2DBenchmarkQuantities : RisingBubble2DBenchmarkQuantities<XNSE_Control> { }

    /// <summary>
    /// Post-processing specific to <see cref="RisingBubble"/>
    /// </summary>
    [Serializable]
    public class RisingBubble2DBenchmarkQuantities<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "BenchmarkQuantities_RisingBubble";

        /// <summary>
        /// hard-coded name for the Rising bubble
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// Header for the rising bubble log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "#timestep", "time", "area", "center of mass - x", "center of mass - y", "circularity", "mean velocity x", "mean velocity y");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using (new FuncTrace()) {
                //var R = ComputeBenchmarkQuantities_RisingBubble();
                var R = RisingBubbleBenchmarkQuantitites.ComputeBenchmarkQuantities_RisingBubble2D(LsTrk, CurrentVel, m_HMForder);
                ITuple RT = R;
                double[] RR = new double[RT.Length];
                for (int i = 0; i < RT.Length; i++)
                    RR[i] = (double)RT[i];

                AppendToLog(iTimestep);
                AppendToLog(PhysTime);
                AppendToLog(RR);


                base.QueryResultTable.LogValue("area", R.area);
                base.QueryResultTable.LogValue("yCM", R.centerY);
                base.QueryResultTable.LogValue("circ", R.circularity);
                base.QueryResultTable.LogValue("riseV", R.MeanVelocityY);

            }
        }


        internal (double area, double centerX, double centerY, double circularity, double VelocityAtCenterX, double VelocityAtCenterY) ComputeBenchmarkQuantities_RisingBubble() {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = m_HMForder;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // area of bubble
            double area = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
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

            // center of mass/geometric center (for incompressible fluid)
            int D = SolverMainOverride.Grid.SpatialDimension;
            MultidimensionalArray center = MultidimensionalArray.Create(1, D);
            CellQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    var nodes_global = LsTrk.GridDat.GlobalNodes.GetValue_Cell(QR.Nodes, i0, Length);
                    EvalResult.Acc(1.0, nodes_global);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            center[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();

            center.Scale(1.0 / area);

            // rise velocity
            MultidimensionalArray VelocityAtCenter = MultidimensionalArray.Create(1, D);

            // integral computation
            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    for (int d = 0; d < D; d++) {
                        CurrentVel[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            VelocityAtCenter[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();
            VelocityAtCenter.Scale(1.0 / area);

            double v_rise = VelocityAtCenter[0, 1];

            //Console.WriteLine("rise velocity = " + v_rise);


            // circularity
            double diamtr_c = Math.Sqrt(4 * area / Math.PI);
            double perimtr_b = 0.0;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetQuadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        perimtr_b += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            double circ = Math.PI * diamtr_c / perimtr_b;

            return (area, center[0, 0], center[0, 1], circ, VelocityAtCenter[0, 0], VelocityAtCenter[0, 1]);
        }
    }

    /// <summary>
    /// actual computation of the rising bubble benchmark quantitites
    /// may be used within the <see cref="InSituPostProcessingModule"> and for postprocessing within notebooks
    /// </summary>
    public static class RisingBubbleBenchmarkQuantitites {


        public static (double area, double centerX, double centerY, double circularity, double MeanVelocityX, double MeanVelocityY) ComputeBenchmarkQuantities_RisingBubble2D(LevelSetTracker LsTrk, XDGField[] VelocityField, int HMForder = 0) {

            int order = GetOrder(LsTrk, VelocityField, HMForder);
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);

            // area of bubble
            double area = ComputeBubbleVolume(LsTrk, vqs, order);

            // center of mass/geometric center (for incompressible fluid)
            double[] center = ComputeCenterOfMass(LsTrk, vqs, order, area);

            // rise velocity
            double[] meanVel = ComputeMeanVelocity(LsTrk, vqs, order, VelocityField, area);
            //double v_rise = meanVel[1];     // upwards velocity in y-Direction

            // circularity
            double circ = ComputeCircularity(LsTrk, SchemeHelper, order, area);


            return (area, center[0], center[1], circ, meanVel[0], meanVel[1]);
        }


        public static (double area, double centerX, double centerY, double centerZ, double sphericity, double MeanVelocityX, double MeanVelocityY, double MeanVelocityZ) ComputeBenchmarkQuantities_RisingBubble3D(LevelSetTracker LsTrk, XDGField[] VelocityField, int HMForder = 0) {

            int order = GetOrder(LsTrk, VelocityField, HMForder);
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);

            // area of bubble
            double area = ComputeBubbleVolume(LsTrk, vqs, order);

            // center of mass/geometric center (for incompressible fluid)
            double[] center = ComputeCenterOfMass(LsTrk, vqs, order, area);

            // rise velocity
            double[] meanVel = ComputeMeanVelocity(LsTrk, vqs, order, VelocityField, area);
            double v_rise = meanVel[1];     // upwards velocity in y-Direction

            // circularity
            double spher = ComputeSphericity(LsTrk, SchemeHelper, order, area);


            return (area, center[0], center[1], center[2], spher, meanVel[0], meanVel[1], meanVel[2]);
        }



        internal static int GetOrder(LevelSetTracker LsTrk, XDGField[] VelocityField, int HMForder = 0) {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                if (HMForder > 0) {
                    order = HMForder;
                } else {
                    int degU = VelocityField[0].Basis.Degree;
                    int quadOrder = 2 * degU;
                    if (LsTrk.CutCellQuadratureType == CutCellQuadratureMethod.Saye) {
                        quadOrder *= 2;
                        quadOrder += 1;
                    }
                    order = quadOrder;
                }
            }

            return order;
        }

        /// <summary>
        /// in 2D case computing the bubble area
        /// </summary>
        /// <returns></returns>
        internal static double ComputeBubbleVolume(LevelSetTracker LsTrk, CellQuadratureScheme vqs, int order) {

            double volume = 0.0;

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        volume += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return volume;
        }

        /// <summary>
        /// Computing the center of mass
        /// </summary>
        /// <returns></returns>
        internal static double[] ComputeCenterOfMass(LevelSetTracker LsTrk, CellQuadratureScheme vqs, int order, double volume) {

            int D = LsTrk.GridDat.SpatialDimension;
            MultidimensionalArray center = MultidimensionalArray.Create(1, D);
            CellQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    NodeSet nodes_global = QR.Nodes.CloneAs();
                    for (int i = i0; i < i0 + Length; i++) {
                        LsTrk.GridDat.TransformLocal2Global(QR.Nodes, i, 1, nodes_global, i);
                        EvalResult.AccSubArray(1.0, nodes_global, new int[] { i - i0, -1, -1 });
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            center[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();

            center.Scale(1.0 / volume);

            return center.ExtractSubArrayShallow(new int[] {0, -1}).To1DArray();
        }

        /// <summary>
        /// Mean bubble velocity
        /// </summary>
        /// <returns></returns>
        internal static double[] ComputeMeanVelocity(LevelSetTracker LsTrk, CellQuadratureScheme vqs, int order, XDGField[] VelocityField, double volume) {

            int D = LsTrk.GridDat.SpatialDimension;
            MultidimensionalArray MeanVelocity = MultidimensionalArray.Create(1, D);
            CellQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    for (int d = 0; d < D; d++) {
                        VelocityField[d].Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            MeanVelocity[0, d] += ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();
            MeanVelocity.Scale(1.0 / volume);

            return MeanVelocity.ExtractSubArrayShallow(new int[] { 0, -1 }).To1DArray();
        }

        /// <summary>
        /// Only for 2D benchmarking
        /// </summary>
        /// <returns></returns>
        internal static double ComputeCircularity(LevelSetTracker LsTrk, XQuadSchemeHelper SchemeHelper, int order, double area) {

            double diamtr_c = Math.Sqrt(4 * area / Math.PI);
            double perimtr_b = 0.0;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetQuadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        perimtr_b += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return Math.PI * diamtr_c / perimtr_b;
        }

        /// <summary>
        /// Only for 3D benchmarking
        /// </summary>
        /// <returns></returns>
        internal static double ComputeSphericity(LevelSetTracker LsTrk, XQuadSchemeHelper SchemeHelper, int order, double area) {

            return 0;
        }


    }
}
