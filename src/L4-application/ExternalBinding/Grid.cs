using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.ExternalBinding {


    /// <summary>
    /// Instantiation of BoSSS grids
    /// </summary>
    public static class Grid_ {

        /// <summary>
        /// Create a BoSSS grid from an OpenFOAM mesh
        /// </summary>
        /// <param name="_ref"></param>
        /// <param name="Dim"></param>
        /// <param name="NoOfPoints"></param>
        /// <param name="NoOfFaces"></param>
        /// <param name="coordinates"></param>
        /// <param name="FaceIndices"></param>
        /// <param name="NoOfFaceVertices"></param>
        /// <param name="ierr"></param>
        unsafe public static void CreateGrid(out int _ref, 
            ref int Dim,
            ref int NoOfPoints,
            ref int NoOfFaces,
            double* coordinates,
            int** FaceIndices,
            int* NoOfFaceVertices,
            out int ierr) {



            _ref = -1;
            ierr = -1;

        }
    }

    
}
