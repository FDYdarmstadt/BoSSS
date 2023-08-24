
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.Connectors;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Tecplot;
using System.IO;

namespace BoSSS.Application.ExternalBinding {


    /// <summary>
    /// Wrapper around one surfaceScalarField to be used in the OpenFOAM binding foam-dg
    /// </summary>
    public class OpenFoamSurfaceField : List<double>, IForeignLanguageProxy {

        OpenFOAMGrid m_grid;

        /// <summary>
        /// Ctor
        /// /summary>
        [CodeGenExport]
        unsafe public OpenFoamSurfaceField(OpenFOAMGrid g, double* values, int nFaces)
            : base(ConstructorHelper(values, nFaces)) //
        {
            m_grid = g;
        }

        /// <summary>
        /// Ctor
        /// /summary>
        unsafe static public double[] ConstructorHelper(double* values, int nFaces)
        {
            double[] safeVals = new double[nFaces];
            for (int i = 0; i < nFaces; i++) {
                safeVals[i] = values[i];
            }
            return safeVals;
        }


        public double GetFlux(int cellIndex, int surfaceIndex){
            // Console.WriteLine();
            // Console.WriteLine(surfaceIndex);
            // Console.WriteLine(m_grid.BoSSSiEdgeToOpenFOAMFace.Count);
            try {
                int faceIndex = m_grid.BoSSSiEdgeToOpenFOAMFace[surfaceIndex];
                // Console.WriteLine(faceIndex);
                if (Double.IsNaN(this[faceIndex]))
                {
                    Console.WriteLine();
                    Console.WriteLine(surfaceIndex);
                    Console.WriteLine(m_grid.BoSSSiEdgeToOpenFOAMFace.Length);
                    Console.WriteLine(faceIndex);
                }
                double faceArea = m_grid.GridData.Edges.GetEdgeArea(surfaceIndex);

                // using (var fs = new FileStream("./out2.txt", FileMode.Append))
                // using (var sw = new StreamWriter(fs))
                // {
                //     sw.WriteLine("faceArea: " + faceArea);
                // }
                return this[faceIndex]/faceArea;
            } catch (Exception e) {
                Console.WriteLine("surfaceIndex: " + surfaceIndex);
                Console.WriteLine("length " + m_grid.BoSSSiEdgeToOpenFOAMFace.Length);
                if (surfaceIndex < m_grid.BoSSSiEdgeToOpenFOAMFace.Length && surfaceIndex >= 0){
                    Console.WriteLine(m_grid.BoSSSiEdgeToOpenFOAMFace[surfaceIndex]);
                }
                throw;
            }
        }

        public double GetFluxIN(ref CommonParams inp){
            return this.GetFlux(inp.jCellIn, inp.iEdge);
        }

        public double GetFluxOUT(ref CommonParams inp){
            return this.GetFlux(inp.jCellOut, inp.iEdge);
        }

        public double GetFlux(ref CommonParamsBnd bnp){
            return this.GetFlux(bnp.jCellIn, bnp.iEdge);
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
