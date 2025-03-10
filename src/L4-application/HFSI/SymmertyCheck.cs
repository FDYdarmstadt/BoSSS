using System;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections.Generic;
using Microsoft.CodeAnalysis.CSharp.Syntax;

using MPI.Wrappers;
using ilPSP;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Statistic;

namespace HFSISolver {

    /// <summary>
    /// 
    /// </summary>
    public static class SymmertyCheck {


        public static MultidimensionalArray GetCellCenters(this IGridData gridData) {

            NodeSet[] centers = gridData.iGeomCells.RefElements.Select(Kref => Kref.Center).ToArray();
            int J = gridData.iGeomCells.NoOfLocalUpdatedCells;
            int D = gridData.SpatialDimension;

            MultidimensionalArray Centers = MultidimensionalArray.Create(J, 1, D);

            int j = 0;
            while(j < J) {
                int iKref = gridData.iGeomCells.GetRefElementIndex(j);
                int Jx = gridData.iGeomCells.GetNoOfSimilarConsecutiveCells(CellInfo.CellType_Mask, j, int.MaxValue);
                var Centers_j = Centers.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j + Jx - 1, 0, D - 1 });

                gridData.TransformLocal2Global(centers[iKref], j, Jx, Centers_j);
                j += Jx;
            }

            return Centers.ResizeShallow(J, D);
        }


        public static MultidimensionalArray Allgather(this MultidimensionalArray local) {

            int D = local.GetLength(1);
            double[] global = local.Storage.MPIAllGatherv();

            var ret = MultidimensionalArray.CreateWrapper(global, global.Length / D, D);
            return ret;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="axis">
        /// 0, 1, 2 for assessing symmetry along the x-axis, y-axis and z-axis, respectively.
        /// </param>
        /// <returns></returns>
        public static double GridAsymmetry(this IGridData gridData, int axis) {

            MultidimensionalArray Centers, MirrorCenters;
            {
                // create cell centers
                var CentersLocal = gridData.GetCellCenters();

                // create mirror centers
                var MirrorCenterLocal = CentersLocal.CloneAs();
                MirrorCenterLocal.ExtractSubArrayShallow(-1, axis).Scale(-1);

                Centers = CentersLocal.Allgather();
                MirrorCenters = MirrorCenterLocal.Allgather();
            }

            var pl = new PointLocalization(Centers);
            double eps = gridData.iGeomCells.h_min.To1DArray().Min().MPIMin() * 0.0001;

            double Asymm = 0.0;
            List<int> foundPoints = new List<int>();
            for(int iPt = 0; iPt < MirrorCenters.GetLength(0); iPt++) {
                var pt = MirrorCenters.GetRowPt(iPt);

                pl.FindNearPoints(foundPoints, eps, pt);
                if(foundPoints.Count <= 0) {
                    Asymm += pl.PointsBB.Diameter * 10; // no symmetric point found -- worst case 

                } else {
                    foreach(int iPtF in foundPoints) {
                        var foundPoint = pl.Points.GetRowPt(iPtF);
                        Asymm += foundPoint.Dist(pt);
                    }
                }
            }


            return Asymm;

        }



        public static double[] FieldAsymmetry(this IEnumerable<DGField> f, int axis, bool[] antiSymmetric) {
            if(antiSymmetric.Count() != f.Count())
                throw new ArgumentException("correlating arrays are supposed to have the same length.");
            int NoOfFields = f.Count();

            MultidimensionalArray Centers, MirrorCenters;
            {
                // create cell centers
                Centers = f.First().Basis.GridDat.GetCellCenters();

                // create mirror centers
                MirrorCenters = Centers.CloneAs();
                MirrorCenters.ExtractSubArrayShallow(-1, axis).Scale(-1);
            }

            var EvenResult = f.EvaluateParallel(Centers);
            var MirrorResult = f.EvaluateParallel(MirrorCenters);

            for(int iF = 0; iF < NoOfFields; iF++) {
                if(antiSymmetric[iF])
                    MirrorResult.ExtractSubArrayShallow(-1, iF).Scale(-1);
            }

            var Error = EvenResult - MirrorResult;

            double[] LocalErrors = NoOfFields.ForLoop(ifld => Error.ExtractSubArrayShallow(ifld, -1).L2Norm().Pow2());

            return LocalErrors.MPISum().Select(x => x.Sqrt()).ToArray();
        }


        public static double FieldAsymmetry(this DGField f, int axis, bool antiSymmetric) {
            return (new DGField[] { f }).FieldAsymmetry(axis, new bool[] { antiSymmetric })[0];
        }
    }
}
