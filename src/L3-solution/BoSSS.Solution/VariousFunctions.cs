using BoSSS.Foundation;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Utils
{
    /// <summary>
    /// various extension functions
    /// </summary>
    static public class VariousFunctions
    {

        /// <summary>
        /// Accumulates an identity matrix, but only for selected variabels
        /// </summary>
        /// <param name="M">matrix to be modified (input/output)</param>
        /// <param name="factor"></param>
        /// <param name="iVar"></param>
        static public void AccEyeSp(this ISparseMatrix M, UnsetteledCoordinateMapping map, int[] iVar, double factor = 1.0) {
            using (new FuncTrace()) {
                if (M.RowPartitioning.LocalLength != M.ColPartition.LocalLength)
                    throw new ArgumentException("supported only for quadratic matrices");
                if (map.LocalLength != M.RowPartitioning.LocalLength)
                    throw new ArgumentException();

                int J = map.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                for(int j = 0; j < J; j++) {
                    foreach(int f in iVar) {
                        int Np = map.BasisS[f].GetLength(j);

                        for(int n = 0; n < Np; n++) {
                            int idx = map.GlobalUniqueCoordinateIndex(f, j, n);

                            M.SetDiagonalElement(idx, M.GetDiagonalElement(idx) + factor);
                        }
                    }
                }
            }
        }
    }
}
