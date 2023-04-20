
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
    /// Wrapper around one surfaceScalarField to be used in the OpenFOAM binding foam-dg
    /// </summary>
    public class OpenFoamSurfaceField : List<double>, IForeignLanguageProxy {

        OpenFOAMGrid m_grid;

        /// <summary>
        /// Ctor
        /// /summary>
        [CodeGenExport]
        public OpenFoamSurfaceField(OpenFOAMGrid g, double[] values)
            : base(values) //
        {
            m_grid = g;
        }


        public double GetFlux(int cellIndex, int surfaceIndex){
            int faceIndex = m_grid.Cells2Faces[cellIndex][surfaceIndex];
            return this[faceIndex];
        }

        public double GetFluxIN(CommonParams inp){
            return this.GetFlux(inp.jCellIn, inp.iEdge);
        }

        public double GetFluxOUT(CommonParams inp){
            return this.GetFlux(inp.jCellOut, inp.iEdge);
        }

        public double GetFlux(CommonParamsBnd bnp){
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
