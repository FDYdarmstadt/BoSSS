using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation {
    

    /// <summary>
    /// various tests to validate the grid adaptation
    /// </summary>
    public static class AMRtests {


        /// <summary>
        /// Verifies whether a gird is symmetric;
        /// 
        /// Background: For initially symmetric meshes, given that the refinement is symmetric, also the result must be symmetric.
        /// One can employ this test to verify this condition.
        /// </summary>
        public static void MeshSymmetryTest(this IGrid g, bool SymmX = true, bool SymmY = true, bool SymmZ = true) {
            var loclCenters = GetCellCenters(g);
            var globCenters = GatherOnRank0(loclCenters);

            int D = g.SpatialDimension;
            if(D <= 1) {
                SymmY = false;
                SymmZ = false;
            }
            if(D <= 2) {
                SymmZ = false;
            }


            var eps = g.iGridData.iGeomCells.h_min.Min().MPIMin() * 1.0e-8;

            if(g.CellPartitioning.MpiRank == 0) {

                var reps = ExtractSymmetricRepresentatives(globCenters, SymmX, SymmY, SymmZ);
                var mirrors = Mirror(reps, SymmX, SymmY, SymmZ);


                bool mirrors_to_glop = PointsAreSubsets(mirrors, globCenters, eps);
                bool glob_to_mirrors = PointsAreSubsets(globCenters, mirrors, eps);

                Assert.IsTrue(glob_to_mirrors && mirrors_to_glop, "Grid Symmetry test failed");

                if(glob_to_mirrors && mirrors_to_glop) {
                    Console.WriteLine("================================================");
                    Console.WriteLine("Grid seems to be symmetric.");
                    Console.WriteLine("================================================");
                }
            }
        }

        static bool PointsAreSubsets(MultidimensionalArray A, MultidimensionalArray B, double eps) {
            var PL = new PointLocalization(B);

            int K = A.NoOfRows;
            for(int k = 0; k < K; k++) {
                var pt = A.GetRowPt(k);

                var nearest = PL.FindNextPoint(pt);
                double dist = (nearest.nearest_pt - pt).L2Norm();
                if(dist > eps)
                    return false;
            }
            return true;
        }




        public static MultidimensionalArray ExtractSymmetricRepresentatives(MultidimensionalArray All, bool SymmX, bool SymmY, bool SymmZ) {
            int D = All.NoOfCols;
            int K = All.NoOfRows;
            var ret = new List<Vector>();
            bool[] Symm = [SymmX, SymmY, SymmZ];
            for(int k = 0; k < K; k++) {
                var pt = All.GetRowPt(k);

                bool isRep = true;
                for(int d = 0; d < D; d++) {
                    if(!Symm[d])
                        continue; // we are not interested in symmetry w.r.t. dimension d
                    
                    if(pt[d] < 0) 
                        // point is on the wrong side of the symmetry plane
                        isRep = false;
                }



                if(isRep)
                    ret.Add(pt);

            }

            var rret = MultidimensionalArray.Create(ret.Count, D);
            for(int k = 0; k < ret.Count; k++) {
                rret.SetRowPt(k, ret[k]);
            }
            return rret;
        }

        public static MultidimensionalArray Mirror(MultidimensionalArray Representatives, bool SymmX, bool SymmY, bool SymmZ) {
            int Lout = Representatives.NoOfRows;
            int D = Representatives.NoOfCols;
            if(SymmX)
                Lout *= 2;
            if(SymmY)
                Lout *= 2;
            if(SymmZ)
                Lout *= 2;

            var ret = MultidimensionalArray.Create(Lout, D);
            int lins = Representatives.NoOfRows;
            ret.ExtractSubArrayShallow([0, 0], [lins - 1, D - 1]).Acc(1.0, Representatives);

            bool[] Symm = [SymmX, SymmY, SymmZ];
            for(int d = 0; d < D; d++) {
                if(Symm[d]) {
                    var beforSymPlane = ret.ExtractSubArrayShallow([0, 0], [lins - 1, D - 1]);
                    beforSymPlane.ExtractSubArrayShallow(-1, d).Scale(-1);
                    var afterSymPlane = ret.ExtractSubArrayShallow([lins, 0], [2 * lins - 1, D - 1]);
                    afterSymPlane.Acc(1.0, beforSymPlane);
                    beforSymPlane.ExtractSubArrayShallow(-1, d).Scale(-1); // undo the sign change
                    lins *= 2;
                }
            }


            return ret;
        }



        public static MultidimensionalArray GetCellCenters(IGrid grid) {
            if(grid == null)
                throw new ArgumentNullException(nameof(grid));
            
            var gridData = grid.iGridData;
            int Jloc = gridData.iGeomCells.NoOfLocalUpdatedCells;
            int D = grid.SpatialDimension;

            var centers = MultidimensionalArray.Create(Jloc, D);

            for(int j = 0; j < Jloc; j++) {
                Vector center = gridData.iGeomCells.GetCenter(j);
                centers.SetRowPt(j, center);
            }

            return centers;
        }


        public static MultidimensionalArray GatherOnRank0(MultidimensionalArray M) {
            int D = M.NoOfCols;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            if(!D.MPIEquals())
                throw new ArgumentException();

            int L = M.NoOfRows;
            int[] Ls = L.MPIGather(0);
            if(M.Storage.Length != M.Length || M.IsContinuous == false)
                M = M.CloneAs(); 

            double[] gatehred = M.Storage.MPIGatherv(Ls?.Select(l => l*D).ToArray());
            if(rank == 0) {
                return MultidimensionalArray.CreateWrapper(gatehred, Ls.Sum(), D);
            } else {
                return null;
            }
        }


    }
}
