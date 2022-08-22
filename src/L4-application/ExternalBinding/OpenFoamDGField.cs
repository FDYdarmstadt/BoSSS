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

namespace BoSSS.Application.ExternalBinding {

    
    /// <summary>
    /// Wrapper around one or more DG fields to be used in the OpenFOAM binding foam-dg
    /// </summary>
    public class OpenFoamDGField : CoordinateVector, IForeignLanguageProxy {

        static DGField[] CtorHelper(OpenFOAMGrid g, int degree, int NoOfComponents) {
            if(NoOfComponents <= 0)
                throw new ArgumentException("number of field components must be positive.");
            var b = new Basis(g.GridData, degree);

            return NoOfComponents.ForLoop(i => new SinglePhaseField(b, "Field" + i));
        }

        /// <summary>
        /// Ctor
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField(OpenFOAMGrid g, int degree, int NoOfComponents) 
            : base(CtorHelper(g, degree, NoOfComponents)) //
        {
        }

        /// <summary>
        /// Getter for DG coordinates
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
            return this[Mapping.LocalUniqueCoordinateIndex(f, j, n)];
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
