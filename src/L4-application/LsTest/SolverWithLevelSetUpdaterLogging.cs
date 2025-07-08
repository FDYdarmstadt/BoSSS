using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.LsTest {

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    public class SolverWithLevelSetUpdaterLogging : SolverWithLevelSetUpdaterLogging<SolverWithLevelSetUpdaterTestControl> { }

    /// <summary>
    /// Benchmark quantities for droplet-testcases, <see cref="Droplet"/>
    /// </summary>
    [Serializable]
    public class SolverWithLevelSetUpdaterLogging<T> : InSituPostProcessingModule where T : SolverWithLevelSetUpdaterControl, new() {

        /// <summary>
        /// SemiAxis
        /// </summary>
        public const string LogfileName = "LevelSetMeasurements";


        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => LogfileName;

        public bool exactSolGiven = false;
        public int iLevSet = 0;
        public double vol_ana;
        public double surf_ana;
        public Func<double[], double, double> contour;

        public SolverWithLevelSetUpdaterLogging() {
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "volume", "volume-err","surface", "surface-err", "contour-err");
            textWriter.WriteLine(header);

        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {

            double[] Props = this.ComputePorperties(phystime);

            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime, Props[0], Props[1], Props[2], Props[3], Props[4]);
            Log.WriteLine(line);
            Log.Flush();

            return;
        }



        double[] ComputePorperties(double phystime) {

            int order = 0;
            if(LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // area/volume
            double volume = 0.0;
            SpeciesId spcId = LsTrk.SpeciesIdS[0]; // inner species always!
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
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

            // surface
            double surface = 0.0;
            var surfElemVol = SchemeHelper.GetLevelSetquadScheme(iLevSet, LsTrk.Regions.GetCutCellMask4LevSet(iLevSet));
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                surfElemVol.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            double[] Props = new double[5];
            Props[0] = volume;
            Props[2] = surface;

            if (exactSolGiven) {
                // contour- error
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(iLevSet, LsTrk.Regions.GetCutCellMask4LevSet(iLevSet));
                int Norm = 2;
                double conterr = DGField.IntegralOverEx(cqs, ((X, U, j) => contour(X, phystime).Pow(Norm)), 2, new DGField[] { (DGField)LsTrk.LevelSets[iLevSet] });
                Props[1] = Math.Abs(volume - vol_ana);
                Props[3] = Math.Abs(surface - surf_ana);
                Props[4] = conterr;
            } else {
                Props[1] = -1; // clearly mark, that these are not evaluated
                Props[3] = -1;
                Props[4] = -1;
            }

            return Props;
        }
    }

    /// <summary>
    /// Benchmark quantities for ls quality
    /// </summary>
    [Serializable]
    public class SolverWithLevelSetUpdaterQualityLogging : SolverWithLevelSetUpdaterQualityLogging<SolverWithLevelSetUpdaterTestControl> { }

    /// <summary>
    /// Benchmark quantities for ls quality
    /// </summary>
    [Serializable]
    public class SolverWithLevelSetUpdaterQualityLogging<T> : InSituPostProcessingModule where T : SolverWithLevelSetUpdaterControl, new() {

        /// <summary>
        /// SemiAxis
        /// </summary>
        public const string LogfileName = "LevelSetQuality";


        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => LogfileName;

        public int iLevSet = 0;

        public SolverWithLevelSetUpdaterQualityLogging() {
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}", "#timestep", "time", "DG-CG-err", "interface-err");
            textWriter.WriteLine(header);

        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {

            double[] Props = this.ComputePorperties(phystime);

            string line = String.Format("{0}\t{1}\t{2}\t{3}", TimestepNo, phystime, Props[0], Props[1]);
            Log.WriteLine(line);
            Log.Flush();

            return;
        }



        double[] ComputePorperties(double phystime) {

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            var PhiCG = SolverMain.IOFields.Single(f => f.Identification == VariableNames.LevelSetCGidx(iLevSet));
            var PhiDG = SolverMain.IOFields.Single(f => f.Identification == VariableNames.LevelSetDGidx(iLevSet));

            double[] Props = new double[2];

            // L2 error between CG and DG levelSet
            Props[0] = PhiCG.L2Error(PhiDG, LsTrk.Regions.GetNearMask4LevSet(iLevSet, LsTrk.NearRegionWidth));           

            // interface error (L2 - norm of PhiDG along interface computed by PhiCG)            
            var surfElemVol = SchemeHelper.GetLevelSetquadScheme(iLevSet, LsTrk.Regions.GetCutCellMask4LevSet(iLevSet));
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                surfElemVol.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    var Result = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    PhiDG.Evaluate(i0, Length, QR.Nodes, Result);
                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Multiply(1.0, Result, Result, 0.0, "ij", "ij", "ij");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        Props[1] += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            Props[1] = Props[1].MPISum().Sqrt();

            return Props;
        }
    }

    /// <summary>
    /// Benchmark quantities for level set test cases see: https://doi.org/10.1016/j.cam.2017.02.016 (36) + (37)
    /// </summary>
    [Serializable]
    public class LevelSetTestLogging : LevelSetTestLogging<SolverWithLevelSetUpdaterTestControl> { }

    /// <summary>
    /// <see cref="LevelSetTestLogging"/>
    /// </summary>
    [Serializable]
    public class LevelSetTestLogging<T> : InSituPostProcessingModule where T : SolverWithLevelSetUpdaterControl, new() {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "LevelSetTestLogging";


        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => LogfileName;

        public int iLevSet = 0;
        public double interface_measure; // Length in 2D, area in 3D
        public Func<double[], double, double> exact_levelset; // often we only know the cyclical "exact" solution, set Logging period accordingly

        public LevelSetTestLogging() {
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "#timestep", "time", "CG-L1-err", "CG-volume-err", "CG-surface-err", "DG-L1-err", "DG-volume-err", "DG-surface-err");
            textWriter.WriteLine(header);

        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {

            double[] Props = this.ComputeProperties(phystime);

            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", TimestepNo, phystime, Props[0], Props[1], Props[2], Props[3], Props[4], Props[5]);
            Log.WriteLine(line);
            Log.Flush();

            return;
        }



        double[] ComputeProperties(double phystime) {

            double[] Props = new double[7];

            var PhiCG = SolverMain.IOFields.Single(f => f.Identification == VariableNames.LevelSetCGidx(iLevSet));
            var PhiDG = SolverMain.IOFields.Single(f => f.Identification == VariableNames.LevelSetDGidx(iLevSet));

            // follow the exact definition of the paper, do not use XDG Quadrature
            CellQuadrature.GetQuadrature(new int[] { 7 }, LsTrk.GridDat, (new CellQuadratureScheme()).Compile(LsTrk.GridDat, PhiCG.Basis.Degree * 2), 
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                    var CGVals = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    var DGVals = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    var ExactVals = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                                        
                    // Evaluation
                    PhiCG.Evaluate(i0, Length, QR.Nodes, CGVals);
                    PhiDG.Evaluate(i0, Length, QR.Nodes, DGVals);
                    for (int i = 0; i < Length; i++) {
                        var GlobalNodes = MultidimensionalArray.Create(QR.NoOfNodes, QR.Nodes.SpatialDimension);
                        LsTrk.GridDat.TransformLocal2Global(QR.Nodes, GlobalNodes, i0 + i);
                        var _ExactVals = QR.NoOfNodes.ForLoop(j => exact_levelset(GlobalNodes.GetRow(j), phystime)).ToArray();
                        ExactVals.ExtractSubArrayShallow(i, -1).SetVector(_ExactVals);
                    }

                    // L1-norm
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 0).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, CGVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(-1.0, ExactVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 0).ApplyAll(x => Math.Abs(x));

                    _EvalResult.ExtractSubArrayShallow(-1, -1, 3).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 3).Acc(1.0, DGVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 3).Acc(-1.0, ExactVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 3).ApplyAll(x => Math.Abs(x));

                    // Apply Heaviside
                    CGVals.ApplyAll(x => x.Heaviside());
                    DGVals.ApplyAll(x => x.Heaviside());
                    ExactVals.ApplyAll(x => x.Heaviside());

                    // volume
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 1).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 1).Acc(1.0, CGVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 1).Acc(-1.0, ExactVals);

                    _EvalResult.ExtractSubArrayShallow(-1, -1, 4).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 4).Acc(1.0, CGVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 4).Acc(-1.0, ExactVals);

                    // exact volume
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 6).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 6).Acc(1.0, ExactVals);

                    // interface
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 2).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 2).Acc(1.0, CGVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 2).Acc(-1.0, ExactVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 2).ApplyAll(x => Math.Abs(x));

                    _EvalResult.ExtractSubArrayShallow(-1, -1, 5).Clear();
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 5).Acc(1.0, CGVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 5).Acc(-1.0, ExactVals);
                    _EvalResult.ExtractSubArrayShallow(-1, -1, 5).ApplyAll(x => Math.Abs(x));

                }, 
                delegate (int i0, int Length, MultidimensionalArray _IntegrationResult) { 
                    for(int i = 0; i < Length; i++) {
                        for (int j = 0; j < Props.Length; j++) {
                            Props[j] = _IntegrationResult[i, j];
                        }
                    }
                }).Execute();

            for (int j = 0; j < Props.Length; j++) {
                Props[j] = Props[j].MPISum();
            }

            // L1-norm is finished, CG Value makes no sense here, due to farfield being +/- 1

            // volume
            Props[1] = 1.0 / Props[6] * Math.Abs(Props[1]);
            Props[4] = 1.0 / Props[6] * Math.Abs(Props[4]);

            // interface
            Props[2] = 1.0 / interface_measure * Props[2];
            Props[5] = 1.0 / interface_measure * Props[5];

            return Props.Take(6).ToArray();
        }
    }
}
