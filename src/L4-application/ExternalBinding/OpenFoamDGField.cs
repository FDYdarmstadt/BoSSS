using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.Connectors;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Runtime.InteropServices;
using System.Text.Json;
using BoSSS.Solution.Tecplot;

namespace BoSSS.Application.ExternalBinding {


    /// <summary>
    /// Wrapper around one or more DG fields to be used in the OpenFOAM binding foam-dg
    /// </summary>
    public class OpenFoamDGField : CoordinateVector, IForeignLanguageProxy {

        int m_degree;
        int m_noOfCoefficients;
        int m_noOfComponents;

        static DGField[] CtorHelper(OpenFOAMGrid g, int degree, int NoOfComponents, string fieldName) {
            if (NoOfComponents <= 0)
                throw new ArgumentException("number of field components must be positive.");
            int managedDegree = ReadManagedDegree(degree, fieldName);
            var basis = new Basis(g.GridData, managedDegree);
            return NoOfComponents.ForLoop(i => new SinglePhaseField(basis,
                String.IsNullOrEmpty(fieldName) ? "Field" + i : (NoOfComponents == 1 ? fieldName : fieldName + i)));
        }

        // Read on the managed side only. In particular, do not change OpenFOAM's
        // compile-time storage degree when the case requests more managed modes.
        static int ReadManagedDegree(int nativeDegree, string fieldName) {
            string path = Path.GetFullPath(Path.Combine("bosss", "config.json"));
            if (String.IsNullOrEmpty(fieldName) || !File.Exists(path))
                return nativeDegree;
            try {
                using (var document = JsonDocument.Parse(File.ReadAllText(path))) {
                    var root = document.RootElement;
                    if (root.ValueKind != JsonValueKind.Object)
                        throw new FormatException("Expected a JSON object.");
                    // Preserve the historical case-file spelling.
                    if (!root.TryGetProperty("DgDgrees", out var degrees))
                        return nativeDegree;
                    if (degrees.ValueKind != JsonValueKind.Object)
                        throw new FormatException("DgDgrees must be an object mapping field names to degrees.");
                    if (!degrees.TryGetProperty(fieldName, out var value))
                        return nativeDegree;
                    if (value.ValueKind != JsonValueKind.Number || !value.TryGetInt32(out int degree) || degree < 0)
                        throw new FormatException("The field degree must be a non-negative integer.");
                    return degree;
                }
            } catch (Exception error) when (error is JsonException || error is FormatException || error is IOException) {
                throw new InvalidDataException("Cannot read DG degree for field '" + fieldName + "' from '" + path + "'.", error);
            }
        }

        /// <summary>Creates an unnamed field using the supplied degree.</summary>
        [CodeGenExport]
        public OpenFoamDGField(OpenFOAMGrid g, int degree, int NoOfComponents)
            : this(g, degree, NoOfComponents, (string)null) {
        }

        /// <summary>
        /// Creates a named field; the UTF-8 name selects an optional managed degree
        /// in the current case's bosss/config.json. The pointer is only read during construction.
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField(OpenFOAMGrid g, int degree, int NoOfComponents, IntPtr fieldNameInC)
            : this(g, degree, NoOfComponents,
                fieldNameInC == IntPtr.Zero ? null : Marshal.PtrToStringUTF8(fieldNameInC)) {
        }

        /// <summary>Managed equivalent of the named native constructor.</summary>
        public OpenFoamDGField(OpenFOAMGrid g, int degree, int NoOfComponents, string fieldName)
            : base(CtorHelper(g, degree, NoOfComponents, fieldName)) {
            m_degree = this.Fields[0].Basis.Degree;
            m_noOfCoefficients = this.Fields[0].Basis.Length;
            m_noOfComponents = NoOfComponents;
            if (!String.IsNullOrEmpty(fieldName)) {
                Console.WriteLine("BoSSS field '{0}': native degree {1}, managed degree {2}, {3} coefficients per cell.",
                    fieldName, degree, m_degree, m_noOfCoefficients);
            }
        }

        /// <summary>
        /// Setter for DG coordinates
        /// </summary>
        /// <param name="f">field index</param>
        /// <param name="j">cell index</param>
        /// <param name="n">DG coordinate/mode index</param>
        /// <param name="val">new value</param>
        [CodeGenExport]
        public void SetDGcoordinate(int f, int j, int n, double val) {
            this[Mapping.LocalUniqueCoordinateIndex(f, j, n)] = val;
        }

        /// <summary>
        /// Getter for DG coordinates
        /// </summary>
        /// <param name="f">field index</param>
        /// <param name="j">cell index</param>
        /// <param name="n">DG coordinate/mode index</param>
        /// <returns>value of respective dg coordinate</returns>
        [CodeGenExport]
        public double GetDGcoordinate(int f, int j, int n) {
            if (n >= this.m_noOfCoefficients) {
                return 0;
            }
            try{
                return this[Mapping.LocalUniqueCoordinateIndex(f, j, n)];
            } catch (Exception e) {
                Console.WriteLine(e);
                Console.WriteLine("f: " + f);
                Console.WriteLine("j: " + j);
                Console.WriteLine("n: " + n);
                Console.WriteLine("length: " + this.Count);
                Console.WriteLine("local unique index: " + Mapping.LocalUniqueCoordinateIndex(f, j, n));
                throw e;
            }
        }

        /// <summary>
        /// Set the mean value of a cell
        /// </summary>
        /// <param name="f">field index</param>
        /// <param name="j">cell index</param>
        /// <returns>value of respective dg coordinate</returns>
        [CodeGenExport]
        public void SetMean(int f, int j, double val) {
            SinglePhaseField spf = this.Fields[f] as SinglePhaseField;
            spf.SetMeanValue(j, val);
        }


        /// <summary>
        /// Obtain the mean value of a cell (useful for projection down to DG0 for use in FVM schemes)
        /// </summary>
        /// <param name="f">field index</param>
        /// <param name="j">cell index</param>
        /// <returns>value of respective dg coordinate</returns>
        [CodeGenExport]
        public double GetMean(int f, int j) {
            SinglePhaseField sff = this.Fields[f] as SinglePhaseField;
            return sff.GetMeanValue(j);
        }

        IntPtr m_ForeignPtr;

        /// <summary>
        /// %
        /// </summary>
        public void _SetForeignPointer(IntPtr ptr) {
            if(ptr == IntPtr.Zero) {
                m_ForeignPtr = IntPtr.Zero;
            } else {

                if(m_ForeignPtr != IntPtr.Zero) {
                    throw new ApplicationException("already registered");
                }
                m_ForeignPtr = ptr;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IntPtr _GetForeignPointer() {
            return m_ForeignPtr;
        }
    }
}
