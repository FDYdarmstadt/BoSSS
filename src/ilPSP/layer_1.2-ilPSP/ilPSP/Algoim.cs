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
using static ilPSP.Utils.UnsafeAlgoim;
using System.Xml;
using System.Runtime.CompilerServices;

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
                case Parallelism.SEQ: return new string[] { "Algoimwrapper.dll", "Algoimwrapper.dll", "Algoimwrapper.dll", "libBoSSSnative_seq.so" };
                default: throw new ArgumentOutOfRangeException();
            }
        }

        public static GetNameMangling[] GetGetNameMangling(Parallelism par) {
            switch (par) {
                case Parallelism.SEQ: return new GetNameMangling[] { DynLibLoader.BoSSS_Prefix, DynLibLoader.BoSSS_Prefix, DynLibLoader.BoSSS_Prefix, DynLibLoader.BoSSS_Prefix };
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

		[StructLayout(LayoutKind.Sequential)]
        public struct QuadSchemeComboUnmanaged {
            public int dimension;
            public int sizeSurf;   //size of surface quadrature rule
            public int sizeVol;    //size of voluume quadrature rule
            public IntPtr nodes;   //nodes of surface(first) and volume(second) quadrature
            public IntPtr weights; //weights of surface(first) and volume(second) quadrature

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

            public void OutputQuadratureRuleAsVtpXML(string filePath) {
                var q = this;
                int dim = q.dimension;

                if (dim != 2 && dim != 3) {
                    Console.Error.WriteLine("XML output is supported only for 2D and 3D schemes.");
                }

                try {
                    using (XmlWriter writer = XmlWriter.Create(filePath, new XmlWriterSettings { Indent = true })) {
                        writer.WriteStartDocument();
                        writer.WriteStartElement("VTKFile");
                        writer.WriteAttributeString("type", "PolyData");
                        writer.WriteAttributeString("version", "0.1");
                        writer.WriteAttributeString("byte_order", "LittleEndian");

                        writer.WriteStartElement("PolyData");
                        writer.WriteStartElement("Piece");
                        writer.WriteAttributeString("NumberOfPoints", q.length.ToString());
                        writer.WriteAttributeString("NumberOfVerts", q.length.ToString());
                        writer.WriteAttributeString("NumberOfLines", "0");
                        writer.WriteAttributeString("NumberOfStrips", "0");
                        writer.WriteAttributeString("NumberOfPolys", "0");

                        // Points
                        writer.WriteStartElement("Points");
                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Float32");
                        writer.WriteAttributeString("Name", "Points");
                        writer.WriteAttributeString("NumberOfComponents", "3");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.length; i++) {
                            writer.WriteString($"{q.nodes[i * dim]} {q.nodes[i * dim + 1]} {(dim == 3 ? q.nodes[i * dim + 2] : 0.0)}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // Points

                        // Verts
                        writer.WriteStartElement("Verts");
                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Int32");
                        writer.WriteAttributeString("Name", "connectivity");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.length; i++) {
                            writer.WriteString($"{i}\n");
                        }

                        writer.WriteEndElement(); // DataArray

                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Int32");
                        writer.WriteAttributeString("Name", "offsets");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 1; i <= q.length; i++) {
                            writer.WriteString($"{i}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // Verts

                        // PointData
                        writer.WriteStartElement("PointData");
                        writer.WriteAttributeString("Scalars", "w");

                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Float32");
                        writer.WriteAttributeString("Name", "w");
                        writer.WriteAttributeString("NumberOfComponents", "1");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.length; i++) {
                            writer.WriteString($"{q.weights[i]}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // PointData

                        writer.WriteEndElement(); // Piece
                        writer.WriteEndElement(); // PolyData
                        writer.WriteEndElement(); // VTKFile

                        writer.WriteEndDocument();
                    }
                } catch (Exception ex) {
                    Console.WriteLine("Error opening file: " + ex.Message);
                }
            }
        }

        public struct QuadSchemeCombo {
            public int dimension;
            public int lengthSurf;   //length of surface quadrature rule
            public int lengthVol;    //length of volume quadrature rule
            public double[] nodes;
            public double[] weights;

            // Constructor to create QuadScheme from QuadSchemeUnmanaged
            public QuadSchemeCombo(QuadSchemeComboUnmanaged unmanagedQuadScheme) {
                dimension = unmanagedQuadScheme.dimension;
                lengthSurf = unmanagedQuadScheme.sizeSurf;
                lengthVol = unmanagedQuadScheme.sizeVol;
                int length = lengthSurf + lengthVol;    //total length of nodes and weights (surface + volume)

                // Copy data from IntPtr to managed double[] arrays
                nodes = new double[length * dimension];
                Marshal.Copy(unmanagedQuadScheme.nodes, nodes, 0, length * dimension);

                weights = new double[length];
                Marshal.Copy(unmanagedQuadScheme.weights, weights, 0, length);
            }

            public void OutputQuadratureRuleAsVtpXML(string filePath) {
                var q = this;
                int dim = q.dimension;
                string filenameSurf = filePath + "_surf.vtp";
                string filenameVol = filePath + "_vol.vtp";

                if (dim != 2 && dim != 3) {
                    Console.Error.WriteLine("XML output is supported only for 2D and 3D schemes.");
                }

                try {
                    using (XmlWriter writer = XmlWriter.Create(filenameSurf, new XmlWriterSettings { Indent = true })) {
                        writer.WriteStartDocument();
                        writer.WriteStartElement("VTKFile");
                        writer.WriteAttributeString("type", "PolyData");
                        writer.WriteAttributeString("version", "0.1");
                        writer.WriteAttributeString("byte_order", "LittleEndian");

                        writer.WriteStartElement("PolyData");
                        writer.WriteStartElement("Piece");
                        writer.WriteAttributeString("NumberOfPoints", q.lengthSurf.ToString());
                        writer.WriteAttributeString("NumberOfVerts", q.lengthSurf.ToString());
                        writer.WriteAttributeString("NumberOfLines", "0");
                        writer.WriteAttributeString("NumberOfStrips", "0");
                        writer.WriteAttributeString("NumberOfPolys", "0");

                        // Points
                        writer.WriteStartElement("Points");
                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Float32");
                        writer.WriteAttributeString("Name", "Points");
                        writer.WriteAttributeString("NumberOfComponents", "3");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.lengthSurf; i++) {
                            writer.WriteString($"{q.nodes[i * dim]} {q.nodes[i * dim + 1]} {(dim == 3 ? q.nodes[i * dim + 2] : 0.0)}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // Points

                        // Verts
                        writer.WriteStartElement("Verts");
                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Int32");
                        writer.WriteAttributeString("Name", "connectivity");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.lengthSurf; i++) {
                            writer.WriteString($"{i}\n");
                        }

                        writer.WriteEndElement(); // DataArray

                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Int32");
                        writer.WriteAttributeString("Name", "offsets");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 1; i <= q.lengthSurf; i++) {
                            writer.WriteString($"{i}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // Verts

                        // PointData
                        writer.WriteStartElement("PointData");
                        writer.WriteAttributeString("Scalars", "w");

                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Float32");
                        writer.WriteAttributeString("Name", "w");
                        writer.WriteAttributeString("NumberOfComponents", "1");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.lengthSurf; i++) {
                            writer.WriteString($"{q.weights[i]}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // PointData

                        writer.WriteEndElement(); // Piece
                        writer.WriteEndElement(); // PolyData
                        writer.WriteEndElement(); // VTKFile

                        writer.WriteEndDocument();
                    }
                    using (XmlWriter writer = XmlWriter.Create(filenameVol, new XmlWriterSettings { Indent = true })) {
                        writer.WriteStartDocument();
                        writer.WriteStartElement("VTKFile");
                        writer.WriteAttributeString("type", "PolyData");
                        writer.WriteAttributeString("version", "0.1");
                        writer.WriteAttributeString("byte_order", "LittleEndian");

                        writer.WriteStartElement("PolyData");
                        writer.WriteStartElement("Piece");
                        writer.WriteAttributeString("NumberOfPoints", q.lengthVol.ToString());
                        writer.WriteAttributeString("NumberOfVerts", q.lengthVol.ToString());
                        writer.WriteAttributeString("NumberOfLines", "0");
                        writer.WriteAttributeString("NumberOfStrips", "0");
                        writer.WriteAttributeString("NumberOfPolys", "0");

                        // Points
                        writer.WriteStartElement("Points");
                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Float32");
                        writer.WriteAttributeString("Name", "Points");
                        writer.WriteAttributeString("NumberOfComponents", "3");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = q.lengthSurf; i < q.lengthSurf + q.lengthVol; i++) {
                            writer.WriteString($"{q.nodes[i * dim]} {q.nodes[i * dim + 1]} {(dim == 3 ? q.nodes[i * dim + 2] : 0.0)}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // Points

                        // Verts
                        writer.WriteStartElement("Verts");
                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Int32");
                        writer.WriteAttributeString("Name", "connectivity");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 0; i < q.lengthVol; i++) {
                            writer.WriteString($"{i}\n");
                        }

                        writer.WriteEndElement(); // DataArray

                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Int32");
                        writer.WriteAttributeString("Name", "offsets");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = 1; i <= q.lengthVol; i++) {
                            writer.WriteString($"{i}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // Verts

                        // PointData
                        writer.WriteStartElement("PointData");
                        writer.WriteAttributeString("Scalars", "w");

                        writer.WriteStartElement("DataArray");
                        writer.WriteAttributeString("type", "Float32");
                        writer.WriteAttributeString("Name", "w");
                        writer.WriteAttributeString("NumberOfComponents", "1");
                        writer.WriteAttributeString("format", "ascii");

                        for (int i = q.lengthSurf; i < q.lengthSurf + q.lengthVol; i++) {
                            writer.WriteString($"{q.weights[i]}\n");
                        }

                        writer.WriteEndElement(); // DataArray
                        writer.WriteEndElement(); // PointData

                        writer.WriteEndElement(); // Piece
                        writer.WriteEndElement(); // PolyData
                        writer.WriteEndElement(); // VTKFile

                        writer.WriteEndDocument();
                    }

                } catch (Exception ex) {
                    Console.WriteLine("Error opening file: " + ex.Message);
                }
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
        _GetComboScheme GetComboScheme;
        _GetVolumeSchemeTwoLS GetVolumeSchemeTwoLS;
        _GetSurfaceSchemeTwoLS GetSurfaceSchemeTwoLS;
		_GetComboSchemeTwoLS GetComboSchemeTwoLS;
#pragma warning restore 649

        // Defines a delegate that can point to the method matching its signature.
        // The number of parameters can be seen too many but utilizing int/double arrays instead of a user defined struct,
        // facilities memory management and handover it to garbage collector.
        public unsafe delegate QuadSchemeUnmanaged _GetVolumeScheme(int dim, int p, int q, int[] sizes, double[] coordinates, double[] LSvalues);

		public unsafe delegate QuadSchemeUnmanaged _GetSurfaceScheme(int dim, int p, int q, int[] sizes, double[] coordinates, double[] LSvalues);

        public unsafe delegate QuadSchemeComboUnmanaged _GetComboScheme(int dim, int p, int q, int[] sizes, double[] coordinates, double[] LSvalues);

        public unsafe delegate QuadSchemeUnmanaged* _GetVolumeSchemeTwoLS(int dim, int p1, int p2, int q, int[] sizes1, int[] sizes2, double[] coordinates1, double[] coordinates2, double[] LSvalues1, double[] LSvalues2);

        public unsafe delegate QuadSchemeUnmanaged* _GetSurfaceSchemeTwoLS(int dim, int p1, int p2, int q, int[] sizes1, int[] sizes2, double[] coordinates1, double[] coordinates2, double[] LSvalues1, double[] LSvalues2);

		public unsafe delegate QuadSchemeComboUnmanaged _GetComboSchemeTwoLS(int dim, int p1, int p2, int q, int[] sizes1, int[] sizes2, double[] coordinates1, double[] coordinates2, double[] LSvalues1, double[] LSvalues2);


        public unsafe _GetVolumeScheme getUnmanagedVolumeScheme {
            get { return GetVolumeScheme; }
        }

        public unsafe _GetSurfaceScheme getUnmanagedSurfaceScheme {
            get { return GetSurfaceScheme; }
        }

        public unsafe _GetComboScheme getUnmanagedComboScheme {
            get { return GetComboScheme; }
        }

        public unsafe _GetVolumeSchemeTwoLS getUnmanagedVolumeSchemeTwoLS {
            get { return GetVolumeSchemeTwoLS; }
        }

        public unsafe _GetSurfaceSchemeTwoLS getUnmanagedSurfaceSchemeTwoLS {
            get { 
                return GetSurfaceSchemeTwoLS; 
            }
        }

		public unsafe _GetComboSchemeTwoLS getUnmanagedComboSchemeTwoLS {
            get { return GetComboSchemeTwoLS; }
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

        /// <summary>
        /// Returns the volume quadrature rules for the given parameters
        /// </summary>
        /// <param name="dim">dimension of space</param>
        /// <param name="p">degree of level set (will be used for Berstein pol. interpolation)</param>
        /// <param name="q">quadrature order</param>
        /// <param name="lengths"> array for the lengths in each axis</param>
        /// <param name="x">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths)</param>
        /// <param name="y">concatenated array for the level set values at nodes (its length = multiplication of lengths)</param>
        /// <returns></returns>
        public static QuadScheme GetSurfaceQuadratureRules(int dim, int p, int q, int[] lengths, double[] x, double[] y) {
            QuadSchemeUnmanaged retC;
            QuadScheme ret;
            if(dim > 1) {
                retC = m_Algoim.getUnmanagedSurfaceScheme(dim, p, q, lengths, x, y);
                ret = new QuadScheme(retC);

            } else if(dim == 1) {
                //Create a fake plane with repeating the same elements
                dim = 2; //upgrade dimension
                double[] pseudoPoints = new double[3] { -1.0, 0.0, 1.0 }; //3 points should suffice, = x 
                int[] pseudoLengths = { lengths[0], pseudoPoints.Length };

                double[] pseudoX = new double[pseudoLengths[0] + pseudoLengths[1]];
                Array.Copy(x, pseudoX, pseudoLengths[0]);
                Array.Copy(pseudoPoints, 0, pseudoX, pseudoLengths[0], pseudoLengths[1]);

                double[] pseudoY = new double[pseudoLengths[0] * pseudoLengths[1]];
                for(int i = 0; i < pseudoLengths[1]; i++) { 
                    Array.Copy(y, 0, pseudoY, i * pseudoLengths[0], pseudoLengths[0]);
                }

                // call for pseudo plane
                retC = m_Algoim.getUnmanagedSurfaceScheme(dim, p, q, pseudoLengths, pseudoX, pseudoY);
                var retPseudo = new QuadScheme(retC); 

                if (retPseudo.length > 0) { // if not empty
                    Debug.Assert(Math.Abs(retPseudo.weights.Sum() - 2) < 1e-6); //single point quadrature rule should have a weight=2
                     Debug.Assert(
                             retPseudo.nodes.Where((_, index) => index % 2 == 0).All(n => Math.Abs(n - retPseudo.nodes[0]) < 1e-6),
                              "All non-pseudo nodes should be identical.");

                    // Create a new quadscheme for the real case (1D)
                    ret = new QuadScheme {
                            dimension = dim - 1,
                            length = 1,
                            weights = new[] { 1.0 }, //all should be same except the pseudo axis
                            nodes = new[] { retPseudo.nodes[0] } //use the first as all should be same except the pseudo axis
                        };
                } else {
                    retPseudo.dimension = dim - 1;
                    ret = retPseudo;
                }
            } else {
                throw new ArgumentOutOfRangeException("Surface quadrules should be at least 1-dimensional");
            }

            retC.FreeMemory();
            return ret;
        }

		/// <summary>
		/// Returns the surface quadrature rules for two level sets
		/// First element is surface rules for level set 1 set and the second element is for level set 2
		/// </summary>
		/// <param name="dim">dimension of space</param>
		/// <param name="p1">degree of level set 1 (will be used for Berstein pol. interpolation)</param>
		/// <param name="p2">degree of level set 2 (will be used for Berstein pol. interpolation)</param>
		/// <param name="q">quadrature order</param>
		/// <param name="lengths1">array for the lengths in each axis for level set 1</param>
		/// <param name="lengths2">array for the lengths in each axis for level set 2</param>
		/// <param name="x1">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths), for level set 1</param>
		/// <param name="x2">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths), for level set 2</param>
		/// <param name="y1">concatenated array for the level set values at nodes (its length = multiplication of lengths), for level set 1</param>
		/// <param name="y2">concatenated array for the level set values at nodes (its length = multiplication of lengths), for level set 2</param>
		/// <returns>Quadrature scheme with nodes and weights (size determined by Algoim)</returns>
		public static QuadScheme[] GetSurfaceQuadratureRules(int dim, int p1, int p2, int q, int[] lengths1, int[] lengths2, double[] x1, double[] x2, double[] y1, double[] y2) {
			QuadScheme[] retArray = new QuadScheme[2]; //first one for level set1 and the second one for level set 2
			unsafe {
				QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedSurfaceSchemeTwoLS(dim, p1, p2, q, lengths1, lengths2, x1, x2, y1, y2);

				for (int k = 0; k < 2; ++k) {
					var scheme = retC[k];
					QuadScheme ret = new QuadScheme(scheme);
					scheme.FreeMemory();
					retArray[k] = ret;
				}
			}
			return retArray;
		}

		/// <summary>
		/// Returns the volume quadrature rules for the given parameters
		/// </summary>
		/// <param name="dim">dimension of space</param>
		/// <param name="p">degree of level set (will be used for Berstein pol. interpolation)</param>
		/// <param name="q">quadrature order</param>
		/// <param name="lengths"> array for the lengths in each axis</param>
		/// <param name="x">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths)</param>
		/// <param name="y">concatenated array for the level set values at nodes (its length = multiplication of lengths)</param>
		/// <returns></returns>
		public static QuadScheme GetVolumeQuadratureRules(int dim, int p, int q, int[] lengths, double[] x, double[] y) {

            QuadSchemeUnmanaged retC = m_Algoim.getUnmanagedVolumeScheme(dim, p, q, lengths, x, y);
			QuadScheme ret = new QuadScheme(retC);
            retC.FreeMemory();
            return ret;
        }

		/// <summary>
		/// Returns the volume quadrature rules for the given parameters
		/// - First element is for level set 1 < 0 and level set 2 < 0
		/// - Second element is for level set 1 < 0 and level set 2 > 0
		/// - Third element is for level set 1 > 0 and level set 2 < 0
        /// - Fourth element is for level set 1 > 0 and level set 2 > 0
		/// </summary>
		/// <param name="dim">dimension of space</param>
		/// <param name="dim">dimension of space</param>
		/// <param name="p1">degree of level set 1 (will be used for Bernstein pol. interpolation)</param>
		/// <param name="p2">degree of level set 2 (will be used for Bernstein pol. interpolation)</param>
		/// <param name="q">quadrature order</param>
		/// <param name="lengths1">array for the lengths in each axis for level set 1</param>
		/// <param name="lengths2">array for the lengths in each axis for level set 2</param>
		/// <param name="x1">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths), for level set 1</param>
		/// <param name="x2">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths), for level set 2</param>
		/// <param name="y1">concatenated array for the level set values at nodes (its length = multiplication of lengths), for level set 1</param>
		/// <param name="y2">concatenated array for the level set values at nodes (its length = multiplication of lengths), for level set 2</param>
		/// <returns>Quadrature scheme with nodes and weights (size determined by Algoim)</returns>
		public static QuadScheme[] GetVolumeQuadratureRules(int dim, int p1, int p2, int q, int[] lengths1, int[] lengths2, double[] x1, double[] x2, double[] y1, double[] y2) {
			QuadScheme[] retArray = new QuadScheme[4];
			unsafe {
				QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedVolumeSchemeTwoLS(dim, p1, p2, q, lengths1, lengths2, x1, x2, y1, y2);
				for (int k = 0; k < 4; ++k) {
					var scheme = retC[k];
					QuadScheme ret = new QuadScheme(scheme);
					scheme.FreeMemory();
					retArray[k] = ret;
				}
			}
			return retArray;
		}

		/// <summary>
		/// Returns the surface + volume quadrature rules for the given parameters
		/// </summary>
		/// <param name="dim">dimension of space</param>
		/// <param name="p">degree of level set (will be used for Berstein pol. interpolation)</param>
		/// <param name="q">quadrature order</param>
		/// <param name="lengths"> array for the lengths in each axis</param>
		/// <param name="x">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths)</param>
		/// <param name="y">concatenated array for the level set values at nodes (its length = multiplication of lengths)</param>
		/// <returns></returns>
		public static QuadSchemeCombo GetComboQuadratureRules(int dim, int p, int q, int[] lengths, double[] x, double[] y) {

            QuadSchemeComboUnmanaged retC = m_Algoim.getUnmanagedComboScheme(dim, p, q, lengths, x, y);
            QuadSchemeCombo ret = new QuadSchemeCombo(retC);
            retC.FreeMemory();
            return ret;
        }

		/// <summary>
		/// Returns the volume quadrature rules for the given parameters
		/// </summary>
		/// <param name="dim">dimension of space</param>
		/// <param name="p1">degree of level set 1 (will be used for Berstein pol. interpolation)</param>
		/// <param name="p2">degree of level set 2 (will be used for Berstein pol. interpolation)</param>
		/// <param name="q">quadrature order</param>
		/// <param name="lengths1">array for the lengths in each axis for level set 1</param>
		/// <param name="lengths2">array for the lengths in each axis for level set 2</param>
		/// <param name="x1">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths), for level set 1</param>
		/// <param name="x2">concatenated array for nodes in each axis (not repeated) (its length = sum of lengths), for level set 2</param>
		/// <param name="y1">concatenated array for the level set values at nodes (its length = multiplication of lengths), for level set 1</param>
		/// <param name="y2">concatenated array for the level set values at nodes (its length = multiplication of lengths), for level set 2</param>
		/// <returns>Quadrature scheme with nodes and weights (size determined by Algoim)</returns>
		public static QuadSchemeCombo GetComboQuadratureRules(int dim, int p1, int p2, int q, int[] lengths1, int[] lengths2, double[] x1, double[] x2, double[] y1, double[] y2) {

			QuadSchemeComboUnmanaged retC = m_Algoim.getUnmanagedComboSchemeTwoLS(dim, p1, p2, q, lengths1, lengths2, x1, x2, y1, y2);
			QuadSchemeCombo ret = new QuadSchemeCombo(retC);
			retC.FreeMemory();
			return ret;
		}

		public static QuadScheme GetSurfaceQuadratureRulesTest() {

            // Hardcoded example values
            // Define points_1dy array
            double[] points_1dy = { 4.0, 3.0, 4.0, 0.0, -1.0, 0.0, 4.0, 3.0, 4.0 };

            // Define points_1dx array
            double[] points_1dx = { -1.0, 0.0, 1.0 };
            double[] l = new double[points_1dx.Length * 2];
            Array.Copy(points_1dx, 0, l, 0, points_1dx.Length);
            Array.Copy(points_1dx, 0, l, points_1dx.Length, points_1dx.Length);

            int[] s = { 3, 3 };

            QuadSchemeUnmanaged retC = m_Algoim.getUnmanagedSurfaceScheme(2, 3, 5, s, l, points_1dy);
            QuadScheme ret = new QuadScheme(retC);
            retC.FreeMemory();
            ret.OutputQuadratureRuleAsVtpXML("AlgoimSurfaceTest.vtp");
            return ret;
        }

        public static QuadScheme GetVolumeQuadratureRulesTest() {

            // Hardcoded example values
            // Define points_1dy array
            double[] points_1dy = { 4.0, 3.0, 4.0, 0.0, -1.0, 0.0, 4.0, 3.0, 4.0 };

            // Define points_1dx array
            double[] points_1dx = { -1.0, 0.0, 1.0 };
            double[] l = new double[points_1dx.Length * 2];
            Array.Copy(points_1dx, 0, l, 0, points_1dx.Length);
            Array.Copy(points_1dx, 0, l, points_1dx.Length, points_1dx.Length);

            int[] s = { 3, 3 };

            QuadSchemeUnmanaged retC = m_Algoim.getUnmanagedVolumeScheme(2, 3, 5, s, l, points_1dy);
            QuadScheme ret = new QuadScheme(retC);
            retC.FreeMemory();
            ret.OutputQuadratureRuleAsVtpXML("AlgoimVolTest.vtp");
            return ret;
        }

        public static QuadSchemeCombo GetComboQuadratureRulesTest() {

            // Hardcoded example values
            // Define points_1dy array
            double[] points_1dy = { 4.0, 3.0, 4.0, 0.0, -1.0, 0.0, 4.0, 3.0, 4.0 };

            // Define points_1dx array
            double[] points_1dx = { -1.0, 0.0, 1.0 };
            double[] l = new double[points_1dx.Length * 2];
            Array.Copy(points_1dx, 0, l, 0, points_1dx.Length);
            Array.Copy(points_1dx, 0, l, points_1dx.Length, points_1dx.Length);

            int[] s = { 3, 3 };

            QuadSchemeComboUnmanaged retC = m_Algoim.getUnmanagedComboScheme(2, 3, 5, s, l, points_1dy);
            QuadSchemeCombo ret = new QuadSchemeCombo(retC);
            retC.FreeMemory();
            ret.OutputQuadratureRuleAsVtpXML("AlgoimComboTest");
            return ret;
        }

		public static QuadScheme[] GetSurfaceQuadratureRulesTwoLSTest(int quadOrder = 5) {

            Func<double, double, double> LS1 = (x1, x2) => (-x2); // a linear line at y=0, positive under it: R^2 -> R, i.e., (x1,x2) -> y
            Func<double, double, double> LS2 = (x1, x2) => (x1 - x2); // a linear line at x=y, positive under it: R^2 -> R, i.e., (x1,x2) -> y

			int n = 2; // number of points in one axis (degree + 1) 
			int spaceDim = 2;

			//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
			double[] pointsA = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);
			double[] pointsB = pointsA; //same for B

			//Cartesian pair product for the points
			int numberOfCombinations = (int)Math.Pow(pointsA.Length, spaceDim);
			MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(pointsA, spaceDim);

			// initiate arrays and assign values
			double[] yA = new double[numberOfCombinations];
			double[] yB = new double[numberOfCombinations];

			for (int i = 0; i < numberOfCombinations; i++) {
				yA[i] = LS1(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
				yB[i] = LS2(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
			}


			double[] xA = new double[pointsA.Length * 2]; //same for xB
			Array.Copy(pointsA, 0, xA, 0, pointsA.Length);
			Array.Copy(pointsA, 0, xA, pointsA.Length, pointsA.Length);
            double[] xB = xA;

			int[] sA = { n, n };
			int[] sB = sA;   //same for sB
			
            QuadScheme[] retArray = new QuadScheme[2];
			unsafe { 
			    QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedSurfaceSchemeTwoLS(spaceDim, n, n, quadOrder, sA, sB, xA, xB, yA, yB);

				for (int k = 0; k < 2; ++k) {
					var scheme = retC[k];
					QuadScheme ret = new QuadScheme(scheme);
					scheme.FreeMemory();
                    retArray[k] = ret;
					ret.OutputQuadratureRuleAsVtpXML("AlgoimSurfTwoLSTest" + k + ".vtp");
				}
			}

            return retArray;
		}

		public static QuadScheme[] GetVolumeQuadratureRulesTwoLSTest() {

            Func<double, double, double> LS1 = (x1, x2) => (-x2); // a linear line at y=0, positive under it: R^2 -> R, i.e., (x1,x2) -> y
            Func<double, double, double> LS2 = (x1, x2) => (x1 - x2); // a linear line at x=y, positive under it: R^2 -> R, i.e., (x1,x2) -> y

			int n = 2; // number of points in one axis (degree + 1) 
			int spaceDim = 2;

			//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
			double[] pointsA = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);
			double[] pointsB = pointsA; //same for B

			//Cartesian pair product for the points
			int numberOfCombinations = (int)Math.Pow(pointsA.Length, spaceDim);
			MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(pointsA, spaceDim);

			// initiate arrays and assign values
			double[] yA = new double[numberOfCombinations];
			double[] yB = new double[numberOfCombinations];

			for (int i = 0; i < numberOfCombinations; i++) {
				yA[i] = LS1(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
				yB[i] = LS2(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
			}


			double[] xA = new double[pointsA.Length * spaceDim]; //same for xB
			Array.Copy(pointsA, 0, xA, 0, pointsA.Length);
			Array.Copy(pointsA, 0, xA, pointsA.Length, pointsA.Length);
			double[] xB = xA;

			int[] sA = { n, n };
			int[] sB = sA;   //same for sB

			QuadScheme[] retArray = new QuadScheme[4];
			unsafe {
                QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedVolumeSchemeTwoLS(spaceDim, n, n, 5, sA, sB, xA, xB, yA, yB);

                for (int k = 0; k < 4; ++k) {
                    var scheme = retC[k];
                    QuadScheme ret = new QuadScheme(scheme);
                    scheme.FreeMemory();
                    retArray[k] = ret;
                    ret.OutputQuadratureRuleAsVtpXML("AlgoimVolTwoLSTest" + k + ".vtp");
                }
            }

			double area = 0;

			foreach (var w in retArray[3].weights)
				area += w;

			Console.WriteLine($"Calculated area {area}, theoretical 3/2");
			return retArray;
		}

		public static QuadScheme[] GetVolumeQuadratureRulesTwoLSTest3D(int quadOrder = 5) {
			Func<double, double, double, double> LS1 = (x1, x2, x3) => (-x2); // a linear plane at y=0, positive under it: R^3 -> R, i.e., (x1,x2) -> y
			Func<double, double, double, double> LS2 = (x1, x2, x3) => (x1 - x2); // a linear plane at x=y

			int n = 2; // number of points in one axis (degree + 1) 
			int spaceDim = 3;

			//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
			double[] pointsA = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);
			double[] pointsB = pointsA; //same for B

			//Cartesian pair product for the points
			int numberOfCombinations = (int)Math.Pow(pointsA.Length, spaceDim);
			MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(pointsA, spaceDim);

			// initiate arrays and assign values
			double[] yA = new double[numberOfCombinations];
			double[] yB = new double[numberOfCombinations];

			for (int i = 0; i < numberOfCombinations; i++) {
				yA[i] = LS1(combinations[i, 0], combinations[i, 1], combinations[i, 2]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
				yB[i] = LS2(combinations[i, 0], combinations[i, 1], combinations[i, 2]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
			}


			double[] xA = new double[pointsA.Length * spaceDim]; //same for xB
			Array.Copy(pointsA, 0, xA, 0, pointsA.Length);
			Array.Copy(pointsA, 0, xA, pointsA.Length, pointsA.Length);
			Array.Copy(pointsA, 0, xA, 2*pointsA.Length, pointsA.Length);

			double[] xB = xA;

			int[] sA = { n, n, n };
			int[] sB = sA;   //same for sB

			QuadScheme[] retArray = new QuadScheme[4];
			unsafe {
				QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedVolumeSchemeTwoLS(spaceDim, n, n, quadOrder, sA, sB, xA, xB, yA, yB);

				for (int k = 0; k < 4; ++k) {
					var scheme = retC[k];
					QuadScheme ret = new QuadScheme(scheme);
					scheme.FreeMemory();
					retArray[k] = ret;
					ret.OutputQuadratureRuleAsVtpXML("AlgoimVolTwoLSTest3D" + k + ".vtp");
				}
			}

			double area = 0;

			foreach (var w in retArray[3].weights)
				area += w;

			Console.WriteLine($"Calculated volume {area}, theoretical 3");
			return retArray;
		}

		public static QuadScheme[] GetSurfaceQuadratureRulesTwoLSpConvergenceTest(int degree) {
			Func<double, double, double> LS1 = (x1, x2) => (-x2); // a linear line at y=0, positive under it: R^2 -> R, i.e., (x1,x2) -> y
            Func<double, double, double> LS2 = (x1, x2) => (x2 - Math.Cos(x1 * 0.5 * Math.PI)); //(x2 - Math.Cos(x1 * 0.5 * Math.PI)); 

			int n = 5; // number of points in one axis (degree + 1) 
			int spaceDim = 2;

			//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
			double[] pointsA = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);
			double[] pointsB = pointsA; //same for B

			//Cartesian pair product for the points
			int numberOfCombinations = (int)Math.Pow(pointsA.Length, spaceDim);
			MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(pointsA, spaceDim);

			// initiate arrays and assign values
			double[] yA = new double[numberOfCombinations];
			double[] yB = new double[numberOfCombinations];

			for (int i = 0; i < numberOfCombinations; i++) {
				yA[i] = LS1(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
				yB[i] = LS2(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
			}


			double[] xA = new double[pointsA.Length * 2]; //same for xB
			Array.Copy(pointsA, 0, xA, 0, pointsA.Length);
			Array.Copy(pointsA, 0, xA, pointsA.Length, pointsA.Length);
			double[] xB = xA;

			int[] sA = { n, n };
			int[] sB = sA;   //same for sB

			QuadScheme[] retArray = new QuadScheme[2];
			unsafe {
				QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedSurfaceSchemeTwoLS(2, n, n, degree, sA, sB, xA, xB, yA, yB);

				for (int k = 0; k < 2; ++k) {
					var scheme = retC[k];
					QuadScheme ret = new QuadScheme(scheme);
					scheme.FreeMemory();
					retArray[k] = ret;
					ret.OutputQuadratureRuleAsVtpXML("AlgoimSurfTwoLSTestP" + degree + "_" + k + ".vtp");
				}
			}

			return retArray;
		}

		public static QuadScheme[] GetVolumeQuadratureRulesTwoLSpConvergenceTest(int degree) {
			Func<double, double, double> LS1 = (x1, x2) => (-x2); // a linear line at y=0, positive above it: R^2 -> R, i.e., (x1,x2) -> -y
			Func<double, double, double> LS2 = (x1, x2) => (x2 - Math.Cos(x1 * 0.5 * Math.PI)); 

			int n = 5; // number of points in one axis (degree + 1) 
			int spaceDim = 2;

			//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
			double[] pointsA = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);
			double[] pointsB = pointsA; //same for B

			//Cartesian pair product for the points
			int numberOfCombinations = (int)Math.Pow(pointsA.Length, spaceDim);
			MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(pointsA, spaceDim);

			// initiate arrays and assign values
			double[] yA = new double[numberOfCombinations];
			double[] yB = new double[numberOfCombinations];

			for (int i = 0; i < numberOfCombinations; i++) {
				yA[i] = LS1(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
				yB[i] = LS2(combinations[i, 0], combinations[i, 1]); // calculate the level set function y=LS(x) where x coordinate vector x=(x1, x2)
			}


			double[] xA = new double[pointsA.Length * 2]; //same for xB
			Array.Copy(pointsA, 0, xA, 0, pointsA.Length);
			Array.Copy(pointsA, 0, xA, pointsA.Length, pointsA.Length);
			double[] xB = xA;

			int[] sA = { n, n };
			int[] sB = sA;   //same for sB

			QuadScheme[] retArray = new QuadScheme[4];
			unsafe {
				QuadSchemeUnmanaged* retC = m_Algoim.getUnmanagedVolumeSchemeTwoLS(2, n, n, degree, sA, sB, xA, xB, yA, yB);

				for (int k = 0; k < 4; ++k) {
					var scheme = retC[k];
					QuadScheme ret = new QuadScheme(scheme);
					scheme.FreeMemory();
					retArray[k] = ret;
					ret.OutputQuadratureRuleAsVtpXML("AlgoimVolTwoLSTestP" + degree + "_"+ k + ".vtp");
				}
			}

			double area = 0;

			foreach (var w in retArray[0].weights)
				area += w;

			Console.WriteLine($"P={degree} - calculated area={area}, theoretical area={4 / Math.PI}, abs. error={area - 4 / Math.PI}");
			return retArray;
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
