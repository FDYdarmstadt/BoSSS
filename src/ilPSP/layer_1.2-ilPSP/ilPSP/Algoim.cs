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


using System;
using System.Collections.Generic;
using MPI.Wrappers.Utils;
using MPI.Wrappers;
using ilPSP;
using ilPSP.Utils;
using System.IO;
using System.Globalization;
using System.Diagnostics;
using System.Linq;
using static MPI.Wrappers.Utils.DynLibLoader;
using System.Threading;
using System.Collections.Concurrent;
using System.Runtime.InteropServices;
using static ilPSP.Utils.UnsafeDBLAS;
using static ilPSP.Utils.UnsafeAlgoim;

namespace ilPSP.Utils {


    // Configurations for the DynLibLoader
    internal class Algoim_Libstuff {
        // workaround for .NET bug:
        // https://connect.microsoft.com/VisualStudio/feedback/details/635365/runtimehelpers-initializearray-fails-on-64b-framework
        public static PlatformID[] GetPlatformID(Parallelism par) {
            switch (par) {
                case Parallelism.SEQ: return new PlatformID[] { PlatformID.Win32NT, PlatformID.Win32NT, PlatformID.Win32NT, PlatformID.Unix };
                default: throw new ArgumentOutOfRangeException();
            }
        }

        public static string[] GetLibname(Parallelism par) {
            switch (par) {
                case Parallelism.SEQ: return new string[] { "CppAlgoim.dll", "CppAlgoim.dll", "CppAlgoim.dll", "CppAlgoim.so" };
                default: throw new ArgumentOutOfRangeException();
            }
        }

        public static GetNameMangling[] GetGetNameMangling(Parallelism par) {
            switch (par) {
                case Parallelism.SEQ: return new GetNameMangling[] { s => s, s => s, s => s, DynLibLoader.BoSSS_Prefix };
                default: throw new ArgumentOutOfRangeException();
            }
        }

        public static string[][][] GetPrequesiteLibraries(Parallelism par) {
            return new string[GetLibname(par).Length][][];
        }

        public static int[] GetPointerSizeFilter(Parallelism par) {
            var r = new int[GetLibname(par).Length];
            r.SetAll(-1);
            return r;
        }
    }


    
    /// <summary>
    /// subset of Algoim
    /// </summary>
    public sealed class UnsafeAlgoim : DynLibLoader {


        [StructLayout(LayoutKind.Sequential)]
        public struct QuadSchemeUnmanaged {
            public int dimension;
            public int size;
            public IntPtr nodes;
            public IntPtr weights;

            // Method to free memory allocated for nodes and weights
            public void FreeMemory() {
                Marshal.FreeHGlobal(nodes);
                Marshal.FreeHGlobal(weights);
            }
        }

        // Define the QuadScheme struct in C# (memory management is handled automatically by the garbage collector)
        public struct QuadScheme {
            public int dimension;
            public int length;
            public double[] nodes;
            public double[] weights;

            // Constructor to create QuadScheme from QuadSchemeUnmanaged
            public QuadScheme(QuadSchemeUnmanaged unmanagedQuadScheme) {
                dimension = unmanagedQuadScheme.dimension;
                length = unmanagedQuadScheme.size;

                // Copy data from IntPtr to managed double[] arrays
                nodes = new double[length * dimension];
                Marshal.Copy(unmanagedQuadScheme.nodes, nodes, 0, length * dimension);

                weights = new double[length];
                Marshal.Copy(unmanagedQuadScheme.weights, weights, 0, length);
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        public UnsafeAlgoim(Parallelism par) :
            base(Algoim_Libstuff.GetLibname(par),
                 Algoim_Libstuff.GetPrequesiteLibraries(par),
                 Algoim_Libstuff.GetGetNameMangling(par),
                 Algoim_Libstuff.GetPlatformID(par),
                 Algoim_Libstuff.GetPointerSizeFilter(par)) { }

#pragma warning disable 649
        _GetVolumeScheme GetVolumeScheme;
        _GetSurfaceScheme GetSurfaceScheme;
#pragma warning restore 649

        // Defines a delegate that can point to the method matching its signature.
        public unsafe delegate QuadSchemeUnmanaged _GetVolumeScheme(int dim, int q, int[] sizes, double[] coordinates, double[] LSvalues);

        public unsafe delegate QuadSchemeUnmanaged _GetSurfaceScheme(int dim, int q, int[] sizes, double[] coordinates, double[] LSvalues);

        public unsafe _GetVolumeScheme getUnmanagedVolumeScheme {
            get { return GetVolumeScheme; }
        }

        public unsafe _GetSurfaceScheme getUnmanagedSurfaceScheme {
            get { return GetSurfaceScheme; }
        }


    }



    /// <summary>
    /// some parts of the Algoim interface, which are used by BoSSS;
    /// </summary>
    static public class Algoim {

        /// <summary>
        /// the machine double accuracy: the smallest x, so that 1 + x > 1
        /// </summary>
        public static double MachineEps {
            get {
                double machEps = 1.0d;

                do {
                    machEps /= 2.0d;
                } while((1.0 + machEps) != 1.0);

                return 2 * machEps;
            }
        }

        public readonly static UnsafeAlgoim m_seq_Algoim;
        //public readonly static UnsafeAlgoim m_omp_Algoim;

        static UnsafeAlgoim m_Algoim;

        public static QuadScheme GetSurfaceQuadratureRules() {

            // Hardcoded example values
            // Define points_1dy array
            double[] points_1dy = { 4.0, 3.0, 4.0, 0.0, -1.0, 0.0, 4.0, 3.0, 4.0 };

            // Define points_1dx array
            double[] points_1dx = { -1.0, 0.0, 1.0 };
            double[] l = new double[points_1dx.Length * 2];
            Array.Copy(points_1dx, 0, l, 0, points_1dx.Length);
            Array.Copy(points_1dx, 0, l, points_1dx.Length, points_1dx.Length);

            int[] s = { 3, 3 };

            QuadSchemeUnmanaged retC = m_Algoim.getUnmanagedSurfaceScheme(2, 5, s, l, points_1dy);
            QuadScheme ret = new QuadScheme(retC);
            retC.FreeMemory();
            return ret;
        }

        public static QuadScheme GetVolumeQuadratureRules() {

            // Hardcoded example values
            // Define points_1dy array
            double[] points_1dy = { 4.0, 3.0, 4.0, 0.0, -1.0, 0.0, 4.0, 3.0, 4.0 };

            // Define points_1dx array
            double[] points_1dx = { -1.0, 0.0, 1.0 };
            double[] l = new double[points_1dx.Length * 2];
            Array.Copy(points_1dx, 0, l, 0, points_1dx.Length);
            Array.Copy(points_1dx, 0, l, points_1dx.Length, points_1dx.Length);

            int[] s = { 3, 3 };

            QuadSchemeUnmanaged retC = m_Algoim.getUnmanagedVolumeScheme(2, 5, s, l, points_1dy);
            QuadScheme ret = new QuadScheme(retC);
            retC.FreeMemory();
            return ret;
        }

        internal static void ActivateSEQ() {
            m_Algoim = m_seq_Algoim;
        }


        /// <summary>
        /// static ctor
        /// </summary>
        static Algoim() {
            m_seq_Algoim = new UnsafeAlgoim(Parallelism.SEQ);
            m_Algoim = m_seq_Algoim;
        }


    }
}
