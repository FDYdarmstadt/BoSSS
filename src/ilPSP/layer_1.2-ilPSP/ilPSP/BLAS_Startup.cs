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
        // ==========================================================================================
        // startup selection
        // ==========================================================================================

        /// <summary>
        /// Environment variable for the choice of the BLAS provider: `mkl` (default), `auto` (calibration at startup) or `openblas`
        /// </summary>
        public const string ProviderEnvVar = "BOSSS_BLAS_PROVIDER";

        static bool s_StartupSelectionDone = false;

        /// <summary>
        /// timings of the last startup calibration; null if not executed
        /// </summary>
        public static ProviderCalibrationResult LastCalibration {
            get;
            private set;
        }

        /// <summary>
        /// Selects the BLAS provider, once per process:
        /// must be called after native library search paths and the CPU allocation (<see cref="Environment.InitThreading"/>) are set up,
        /// and before parallel numerical work.
        /// According to <see cref="ProviderEnvVar"/>: MKL by default (on AMD EPYC Zen 4, MKL with the fake-Intel probe outperformed OpenBLAS 0.3.34);
        /// for `auto`, a short calibration on representative shapes, which keeps MKL for ties, incomplete or inconclusive measurements and if OpenBLAS is unavailable.
        /// In the default case, the OpenBLAS library is not even loaded.
        /// </summary>
        public static void SelectProviderAtStartup() {
            if(s_StartupSelectionDone)
                return;
            using(var tr = new FuncTrace()) {
                s_StartupSelectionDone = true;
                string envVal = System.Environment.GetEnvironmentVariable(ProviderEnvVar);
                string req = (envVal ?? "mkl").Trim().ToLowerInvariant();
                if(req.Length == 0) {
                    req = "mkl";
                    envVal = null;
                }
                RequestedProvider = req;

                switch(req) {
                    case "mkl":
                    ProviderSelectionReason = envVal == null ? $"default ({ProviderEnvVar} not set)" : $"explicit request ({ProviderEnvVar}=mkl)";
                    break;

                    case "openblas":
                    if(TryActivateProvider(BLASProvider.OpenBLAS, out string reason)) {
                        ProviderSelectionReason = $"explicit request ({ProviderEnvVar}=openblas)";
                    } else {
                        ProviderSelectionReason = $"explicit request ({ProviderEnvVar}=openblas) failed, using MKL: {reason}";
                        tr.Warning(ProviderSelectionReason);
                    }
                    break;

                    case "auto":
                    if(!OpenBLASAvailable) {
                        ProviderSelectionReason = "auto: OpenBLAS unavailable, using MKL (" + OpenBLASUnavailableReason + ")";
                        break;
                    }
                    try {
                        var cal = CalibrateProviders();
                        LastCalibration = cal;
                        if(cal.OpenBLASWins) {
                            if(TryActivateProvider(BLASProvider.OpenBLAS, out string r2))
                                ProviderSelectionReason = "auto: " + cal.Summary;
                            else
                                ProviderSelectionReason = "auto: activation failed, using MKL: " + r2;
                        } else {
                            ProviderSelectionReason = "auto: " + cal.Summary;
                        }
                    } catch(Exception e) when(e is not ArithmeticException) {
                        tr.Warning("BLAS provider calibration failed: " + e);
                        ProviderSelectionReason = "auto: calibration failed, using MKL: " + e.Message;
                        lock(s_ProviderLock) {
                            s_Provider = BLASProvider.MKL;
                            ApplyProviderAndMode();
                        }
                    }
                    break;

                    default:
                    ProviderSelectionReason = $"unknown value {ProviderEnvVar}={req} (expected mkl, auto or openblas), using MKL";
                    tr.Warning(ProviderSelectionReason);
                    Console.Error.WriteLine("Warning: " + ProviderSelectionReason);
                    break;
                }
                tr.Info($"BLAS provider: {Provider}; {ProviderSelectionReason}");
            }
        }

        /// <summary>
        /// Result of <see cref="CalibrateProviders"/>
        /// </summary>
        public sealed class ProviderCalibrationResult {

            /// <summary>
            /// one measured case
            /// </summary>
            public sealed class Case {
                /// <summary>e.g. `dgemm nn 35x35x49`</summary>
                public string Name;
                /// <summary>seconds per call, best of several batches</summary>
                public double MKL, OpenBLAS;
                /// <summary>OpenBLAS/MKL; less than 1 means OpenBLAS is faster</summary>
                public double Ratio => OpenBLAS / MKL;
                /// <summary>maximum deviation, relative to the expected round-off</summary>
                public double ErrorRatio;
                /// <summary>threaded (true) or sequential (false) measurement</summary>
                public bool Threaded;
            }

            /// <summary> measured cases </summary>
            public List<Case> Cases = new List<Case>();

            /// <summary> geometric mean of the sequential <see cref="Case.Ratio"/> </summary>
            public double SeqGeoMeanRatio;

            /// <summary> worst sequential <see cref="Case.Ratio"/> </summary>
            public double SeqWorstRatio;

            /// <summary> threaded ratio, NaN if not measured </summary>
            public double ThreadedRatio = double.NaN;

            /// <summary> wall time of the calibration in seconds, including the first calls </summary>
            public double WallTime;

            /// <summary> decision </summary>
            public bool OpenBLASWins;

            /// <summary> decision and its reason </summary>
            public string Summary;

            /// <summary>
            /// table of all cases
            /// </summary>
            public override string ToString() {
                var stw = new StringBuilder();
                stw.AppendLine(Summary);
                foreach(var c in Cases)
                    stw.AppendLine(string.Format(CultureInfo.InvariantCulture, "  {0,-28} {1}  MKL {2,10:0.000} us  OpenBLAS {3,10:0.000} us  ratio {4:0.000}  err {5:0.00}",
                        c.Name, c.Threaded ? "thr" : "seq", c.MKL * 1e6, c.OpenBLAS * 1e6, c.Ratio, c.ErrorRatio));
                return stw.ToString();
            }
        }

        /// <summary>
        /// OpenBLAS is selected only if it is faster by at least this factor in the geometric mean over the sequential cases
        /// </summary>
        public const double CalibrationWinMargin = 1.10;

        /// <summary>
        /// OpenBLAS is rejected if any sequential case is slower than MKL by this factor (protects small shapes)
        /// </summary>
        public const double CalibrationWorstCaseLimit = 1.25;

        /// <summary>
        /// OpenBLAS is rejected if threaded OpenBLAS is slower than threaded MKL by this factor
        /// </summary>
        public const double CalibrationThreadedLimit = 1.10;

        /// <summary>
        /// Calibration is aborted (MKL is kept) if it exceeds this wall time (seconds)
        /// </summary>
        public const double CalibrationBudget = 0.5;

        /// <summary>
        /// Correctness-checks and times OpenBLAS against MKL on representative DGEMM/DGEMV shapes (from the Multiply census, `CENSUS-4.md`,
        /// plus some medium square cases), sequential and, if the thread budget is larger than one, threaded.
        /// Uses private buffers; the provider and threading mode are restored on exit.
        /// Must not be called within a parallel section.
        /// </summary>
        public static ProviderCalibrationResult CalibrateProviders() => ProviderCalibration.Run();
    }
}
