/* =======================================================================
Copyright 2026 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using ilPSP.Tracing;

namespace ilPSP.Utils {
    public static partial class BLAS {
        // Numerical startup diagnostics live in their own type, apart from the driver and loading policy.
        static class ProviderCalibration {
            internal static ProviderCalibrationResult Run() {
                if(ilPSP.Environment.InParallelSection)
                    throw new NotSupportedException("calibration must not be called within a parallel section");
                if(!OpenBLASAvailable)
                    throw new NotSupportedException("OpenBLAS unavailable: " + OpenBLASUnavailableReason);

                var ret = new ProviderCalibrationResult();
                var total = Stopwatch.StartNew();
                var ob = s_OpenBLAS.Value.Blas;

                // (M, N, K, transA, transB), Fortran order
                var gemmCases = new (int M, int N, int K, char ta, char tb)[] {
                    (35, 35, 49, 'n', 'n'),    // LECEdgeQuadrature2, batch of cells
                    (35, 35, 49, 't', 'n'),
                    (35, 700, 1029, 'n', 'n'), // LECVolumeQuadrature2, 3D p4
                    (56, 20, 1536, 'n', 't'),  // LECVolumeQuadrature2/NECQuadratureVolume
                    (60, 35, 3, 'n', 'n'),     // short reduction
                    (10, 10, 10, 'n', 'n'),    // tiny
                    (20, 20, 20, 'n', 't'),
                    (128, 128, 128, 'n', 'n'), // medium
                    (300, 300, 300, 't', 'n'),
                };
                var gemvCases = new (int M, int N, char t)[] {
                    (1029, 35, 't'),
                    (35, 1029, 'n'),
                    (200, 200, 't'),
                    (16, 16, 't'),
                };

                int maxA = 0, maxB = 0, maxC = 0;
                foreach(var c in gemmCases) {
                    maxA = Math.Max(maxA, c.M * c.K);
                    maxB = Math.Max(maxB, c.K * c.N);
                    maxC = Math.Max(maxC, c.M * c.N);
                }
                const int Nthr = 512;
                maxA = Math.Max(maxA, Nthr * Nthr); maxB = Math.Max(maxB, Nthr * Nthr); maxC = Math.Max(maxC, Nthr * Nthr);

                var rnd = new Random(4711);
                double[] A = new double[maxA], B = new double[maxB], C0 = new double[maxC], C1 = new double[maxC], C2 = new double[maxC];
                for(int i = 0; i < A.Length; i++) A[i] = rnd.NextDouble() - 0.5;
                for(int i = 0; i < B.Length; i++) B[i] = rnd.NextDouble() - 0.5;
                for(int i = 0; i < C0.Length; i++) C0[i] = rnd.NextDouble() - 0.5;

                bool budgetExceeded = false;
                string correctnessFailure = null;
                try {
                    // sequential: OpenBLAS count 1, sequential MKL; 
                    // the active provider (m_BLAS) is not changed, the OpenBLAS instance is called directly.
                    lock(s_ProviderLock) {
                        SetOpenBLASThreads(1);
                    }

                    unsafe {
                        fixed(double* pA = A, pB = B, pC0 = C0, pC1 = C1, pC2 = C2) {
                            double* _pA = pA, _pB = pB, _pC0 = pC0, _pC1 = pC1, _pC2 = pC2;

                            foreach(var c in gemmCases) {
                                if(total.Elapsed.TotalSeconds > CalibrationBudget) {
                                    budgetExceeded = true;
                                    break;
                                }
                                int M = c.M, N = c.N, K = c.K, ta = c.ta, tb = c.tb;
                                int lda = c.ta == 'n' ? M : K, ldb = c.tb == 'n' ? K : N, ldc = M;
                                double alpha = 1.25, beta = 0.5;
                                int len = M * N;

                                void Run(UnsafeDBLAS lib, double* pC) {
                                    lib.DGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha, _pA, ref lda, _pB, ref ldb, ref beta, pC, ref ldc);
                                }

                                // correctness: both from the same initial C
                                Buffer.MemoryCopy(_pC0, _pC1, len * sizeof(double), len * sizeof(double));
                                Buffer.MemoryCopy(_pC0, _pC2, len * sizeof(double), len * sizeof(double));
                                Run(m_seq_BLAS, _pC1);
                                Run(ob, _pC2);
                                double err = ErrorRatio(_pC1, _pC2, len, K);

                                var cs = new ProviderCalibrationResult.Case() {
                                    Name = $"dgemm {c.ta}{c.tb} {M}x{N}x{K}",
                                    ErrorRatio = err,
                                };
                                ret.Cases.Add(cs);
                                if(!(err <= 1.0)) {
                                    correctnessFailure = $"{cs.Name}: deviation {err:0.0e0} x round-off bound";
                                    break;
                                }
                                // beta = 0.5 repeatedly applied keeps the values bounded
                                (cs.MKL, cs.OpenBLAS) = TimePair(() => Run(m_seq_BLAS, _pC1), () => Run(ob, _pC2), 2.0 * M * N * K);
                            }

                            if(!budgetExceeded && correctnessFailure == null) {
                                foreach(var c in gemvCases) {
                                    if(total.Elapsed.TotalSeconds > CalibrationBudget) {
                                        budgetExceeded = true;
                                        break;
                                    }
                                    int M = c.M, N = c.N, t = c.t, lda = M, inc = 1;
                                    double alpha = 1.25, beta = 0.5;
                                    int lenY = c.t == 'n' ? M : N, K = c.t == 'n' ? N : M;
                                    void Run(UnsafeDBLAS lib, double* pY) {
                                        lib.DGEMV(ref t, ref M, ref N, ref alpha, _pA, ref lda, _pB, ref inc, ref beta, pY, ref inc);
                                    }
                                    Buffer.MemoryCopy(_pC0, _pC1, lenY * sizeof(double), lenY * sizeof(double));
                                    Buffer.MemoryCopy(_pC0, _pC2, lenY * sizeof(double), lenY * sizeof(double));
                                    Run(m_seq_BLAS, _pC1);
                                    Run(ob, _pC2);
                                    double err = ErrorRatio(_pC1, _pC2, lenY, K);
                                    var cs = new ProviderCalibrationResult.Case() {
                                        Name = $"dgemv {c.t} {M}x{N}",
                                        ErrorRatio = err,
                                    };
                                    ret.Cases.Add(cs);
                                    if(!(err <= 1.0)) {
                                        correctnessFailure = $"{cs.Name}: deviation {err:0.0e0} x round-off bound";
                                        break;
                                    }
                                    (cs.MKL, cs.OpenBLAS) = TimePair(() => Run(m_seq_BLAS, _pC1), () => Run(ob, _pC2), 2.0 * M * N);
                                }
                            }

                            // threaded: only relevant if ordinary BLAS runs threaded
                            int budget = OpenBLASThreadBudget;
                            if(!budgetExceeded && correctnessFailure == null && budget > 1 && ilPSP.Environment.OpenMPenabled) {
                                lock(s_ProviderLock) {
                                    SetOpenBLASThreads(budget);
                                }
                                int M = Nthr, N = Nthr, K = Nthr, ta = 'n', tb = 'n';
                                double alpha = 1.0, beta = 0.5;
                                void Run(UnsafeDBLAS lib, double* pC) {
                                    lib.DGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha, _pA, ref M, _pB, ref K, ref beta, pC, ref M);
                                }
                                var cs = new ProviderCalibrationResult.Case() {
                                    Name = $"dgemm nn {M}x{N}x{K} ({budget} thr)",
                                    Threaded = true,
                                    ErrorRatio = double.NaN
                                };
                                (cs.MKL, cs.OpenBLAS) = TimePair(() => Run(m_omp_BLAS, _pC1), () => Run(ob, _pC2), 2.0 * M * N * K);
                                ret.Cases.Add(cs);
                                ret.ThreadedRatio = cs.Ratio;
                            }
                        }
                    }
                } finally {
                    lock(s_ProviderLock) {
                        if(s_Provider == BLASProvider.OpenBLAS)
                            ApplyProviderAndMode();
                        else
                            SetOpenBLASThreads(OpenBLASThreadBudget); // as if loaded with OPENBLAS_NUM_THREADS, e.g. for PETSc
                    }
                }

                var seq = ret.Cases.Where(c => !c.Threaded && c.MKL > 0 && c.OpenBLAS > 0).ToArray();
                ret.SeqGeoMeanRatio = seq.Length > 0 ? Math.Exp(seq.Average(c => Math.Log(c.Ratio))) : double.NaN;
                ret.SeqWorstRatio = seq.Length > 0 ? seq.Max(c => c.Ratio) : double.NaN;
                ret.WallTime = total.Elapsed.TotalSeconds;

                string stats = string.Format(CultureInfo.InvariantCulture,
                    "OpenBLAS/MKL time ratio: sequential geo. mean {0:0.000}, worst {1:0.000}, threaded {2:0.000}; {3} cases in {4:0.000} s",
                    ret.SeqGeoMeanRatio, ret.SeqWorstRatio, ret.ThreadedRatio, ret.Cases.Count, ret.WallTime);
                if(correctnessFailure != null) {
                    ret.OpenBLASWins = false;
                    ret.Summary = "MKL kept, OpenBLAS failed correctness check (" + correctnessFailure + ")";
                } else if(budgetExceeded) {
                    ret.OpenBLASWins = false;
                    ret.Summary = $"MKL kept, calibration exceeded budget of {CalibrationBudget} s (incomplete); " + stats;
                } else if(!(ret.SeqGeoMeanRatio <= 1.0 / CalibrationWinMargin)) {
                    ret.OpenBLASWins = false;
                    ret.Summary = $"MKL kept, OpenBLAS not faster by margin {CalibrationWinMargin}; " + stats;
                } else if(!(ret.SeqWorstRatio <= CalibrationWorstCaseLimit)) {
                    ret.OpenBLASWins = false;
                    ret.Summary = $"MKL kept, OpenBLAS slower than limit {CalibrationWorstCaseLimit} on some case; " + stats;
                } else if(!double.IsNaN(ret.ThreadedRatio) && !(ret.ThreadedRatio <= CalibrationThreadedLimit)) {
                    ret.OpenBLASWins = false;
                    ret.Summary = $"MKL kept, threaded OpenBLAS slower than limit {CalibrationThreadedLimit}; " + stats;
                } else {
                    ret.OpenBLASWins = true;
                    ret.Summary = "OpenBLAS selected; " + stats;
                }
                return ret;
            }

            /// <summary>
            /// max. deviation between <paramref name="X"/> and <paramref name="Y"/>, relative to a round-off bound for a reduction of length <paramref name="K"/>:
            /// inputs are in [-0.5, 0.5] and alpha = 1.25, so the sum of absolute products is bounded by 0.3125*K, |beta*C| by 0.5;
            /// both results may deviate from the exact one by (K + 2)*eps times that.
            /// </summary>
            unsafe static double ErrorRatio(double* X, double* Y, int L, int K) {
                double maxDiff = 0;
                for(int i = 0; i < L; i++) {
                    double d = Math.Abs(X[i] - Y[i]);
                    if(double.IsNaN(d))
                        return double.PositiveInfinity;
                    maxDiff = Math.Max(maxDiff, d);
                }
                double bound = 2.0 * (K + 2) * MachineEps * (0.3125 * K + 0.5);
                return maxDiff / bound;
            }

            /// <summary>
            /// best time per call (seconds) of each action; batches of both are alternated to cancel drift (e.g. CPU frequency)
            /// </summary>
            static (double tA, double tB) TimePair(Action a, Action b, double flops) {
                // warm-up
                a(); b(); a(); b();

                // batch size: about 0.2 ms of work, assuming 5 GFlop/s
                int reps = (int)Math.Max(1, Math.Min(1000, 0.2e-3 * 5e9 / flops));
                double bestA = double.MaxValue, bestB = double.MaxValue;
                var sw = new Stopwatch();
                for(int r = 0; r < 5; r++) {
                    sw.Restart();
                    for(int i = 0; i < reps; i++) a();
                    bestA = Math.Min(bestA, sw.Elapsed.TotalSeconds / reps);

                    sw.Restart();
                    for(int i = 0; i < reps; i++) b();
                    bestB = Math.Min(bestB, sw.Elapsed.TotalSeconds / reps);
                }
                return (bestA, bestB);
            }
        }
    }
}
