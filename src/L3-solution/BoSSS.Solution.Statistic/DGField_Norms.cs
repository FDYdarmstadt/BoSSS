using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Statistic;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Statistic {

    /// <summary>
    /// L2 and H1 distances for DG fields on different meshes
    /// </summary>
    public static class DGField_Norms {

        /// <summary>
        /// Projects a DG field <paramref name="source"/>, which may be defined on some different mesh,
        /// onto the DG field <paramref name="target"/>.
        /// </summary>
        static public void ProjectFromForeignGrid(this ConventionalDGField target, double alpha, ConventionalDGField source, CellQuadratureScheme scheme = null) {
            using (new FuncTrace()) {
                //Console.WriteLine(string.Format("Projecting {0} onto {1}... ", source.Identification, target.Identification));
                int maxDeg = Math.Max(target.Basis.Degree, source.Basis.Degree);
                var CompQuadRule = scheme.SaveCompile(target.GridDat, maxDeg * 3 + 3); // use over-integration
                int D = target.GridDat.SpatialDimension;



                if (object.ReferenceEquals(source.GridDat, target.GridDat)) {
                    // +++++++++++++++++
                    // equal grid branch
                    // +++++++++++++++++

                    target.ProjectField(alpha, source.Evaluate, CompQuadRule);

                } else {
                    // +++++++++++++++++++++
                    // different grid branch
                    // +++++++++++++++++++++

                    if (source.GridDat.SpatialDimension != D) {
                        throw new ArgumentException("Spatial Dimension Mismatch.");
                    }

                    var eval = new FieldEvaluation(GridHelper.ExtractGridData(source.GridDat));


                    // 
                    // Difficulty with MPI parallelism:
                    // While the different meshes may represent the same geometrical domain, 
                    // their MPI partitioning may be different.
                    //
                    // Solution Approach: we circulate unlocated points among all MPI processes.
                    //

                    int MPIsize = target.GridDat.MpiSize;


                    // pass 1: collect all points
                    // ==========================
                    MultidimensionalArray allNodes = null;
                    void CollectPoints(MultidimensionalArray input, MultidimensionalArray output) {

                        Debug.Assert(input.Dimension == 2);
                        Debug.Assert(input.NoOfCols == D);

                        if (allNodes == null) {
                            allNodes = input.CloneAs();
                        } else {
                            /*
                            var newNodes = MultidimensionalArray.Create(allNodes.NoOfRows + input.NoOfRows, D);
                            newNodes.ExtractSubArrayShallow(new[] { 0, 0 }, new[] { allNodes.NoOfRows - 1, D - 1 })
                                .Acc(1.0, allNodes);
                            newNodes.ExtractSubArrayShallow(new[] {  allNodes.NoOfRows, 0 }, new[] { newNodes.NoOfRows - 1, D - 1 })
                                .Acc(1.0, allNodes);
                            */
                            allNodes = allNodes.CatVert(input);
                        }
                        Debug.Assert(allNodes.NoOfCols == D);

                    }
                    target.ProjectField(alpha, CollectPoints, CompQuadRule);

                    int L = allNodes != null ? allNodes.GetLength(0) : 0;

                    // evaluate 
                    // ========

                    //allNodes = MultidimensionalArray.Create(2, 2);
                    //allNodes[0, 0] = -1.5;
                    //allNodes[0, 1] = 2.0;
                    //allNodes[1, 0] = -0.5;
                    //allNodes[1, 1] = 2.0;
                    //L = 2;
                    var Res = L > 0 ? MultidimensionalArray.Create(L, 1) : default(MultidimensionalArray);
                    int NoOfUnlocated = eval.EvaluateParallel(1.0, new DGField[] { source }, allNodes, 0.0, Res);

                    int TotalNumberOfUnlocated = NoOfUnlocated.MPISum();
                    if (TotalNumberOfUnlocated > 0) {
                        Console.Error.WriteLine($"WARNING: {TotalNumberOfUnlocated} unlocalized points in 'ProjectFromForeignGrid(...)'");
                    }

                    // perform the real projection
                    // ===========================


                    int lc = 0;
                    void ProjectionIntegrand(MultidimensionalArray input, MultidimensionalArray output) {

                        Debug.Assert(input.Dimension == 2);
                        Debug.Assert(input.NoOfCols == D);

                        int LL = input.GetLength(0);
                        output.Set(Res.ExtractSubArrayShallow(new int[] { lc, 0 }, new int[] { lc + LL - 1, -1 }));
                        lc += LL;
                    }
                    target.ProjectField(alpha, ProjectionIntegrand, CompQuadRule);

                }
            }
        }


        static public double L2Distance(this XDGField A, XDGField B, string[] speciesNames, bool IgnoreMeanValue = false) {
<<<<<<< HEAD
=======

>>>>>>> b7f11eee5f (before adding new level-set evolver)
            XDGField fine, coarse;
            if (A.GridDat.CellPartitioning.TotalLength > B.GridDat.CellPartitioning.TotalLength) {
                fine = A;
                coarse = B;
            } else {
                fine = B;
                coarse = A;
            }

<<<<<<< HEAD

=======
            //Debugger.Launch();
>>>>>>> b7f11eee5f (before adding new level-set evolver)

            var trackerB = coarse.Basis.Tracker;
            int maxDeg = Math.Max(A.Basis.Degree, B.Basis.Degree);
            int quadOrder = maxDeg * 3 + 3;

            var schemeFactory = trackerB.GetXDGSpaceMetrics(speciesNames.Select(spc => trackerB.GetSpeciesId(spc)), quadOrder).XQuadSchemeHelper;

            double totNorm = 0;
<<<<<<< HEAD
            foreach (string spc in speciesNames) {

=======
            foreach(string spc in speciesNames) {
            
>>>>>>> b7f11eee5f (before adding new level-set evolver)
                var spc_id = trackerB.GetSpeciesId(spc);

                if (IgnoreMeanValue == true)
                    throw new NotImplementedException();

<<<<<<< HEAD

=======
                
>>>>>>> b7f11eee5f (before adding new level-set evolver)
                var A_spc = fine.GetSpeciesShadowField(spc);
                var B_spc = coarse.GetSpeciesShadowField(spc);

                totNorm += A_spc.L2Distance(B_spc, scheme: schemeFactory.GetVolumeQuadScheme(spc_id)).Pow2();
            }

            totNorm = totNorm.Sqrt();
            return totNorm;

        }



<<<<<<< HEAD
=======

>>>>>>> b7f11eee5f (before adding new level-set evolver)

        /// <summary>
        /// Approximate L2 distance between two DG fields; this also supports DG fields on different meshes, 
        /// it could be used for convergence studies.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="IgnoreMeanValue">
        /// if true, the mean value (mean over entire domain) will be subtracted - this mainly useful for comparing pressures 
        /// </param>
        /// <param name="scheme">
        /// a cell quadrature scheme on the coarse of the two meshes
        /// </param>
        /// <returns></returns>
        static public double L2Distance(this ConventionalDGField A, ConventionalDGField B, bool IgnoreMeanValue = false, CellQuadratureScheme scheme = null) {
            using (var tr = new FuncTrace()) {
                int maxDeg = Math.Max(A.Basis.Degree, B.Basis.Degree);
                int quadOrder = maxDeg * 3 + 3;
                tr.Info($"L2Distance {A.Identification} -- {B.Identification}...");
                tr.Info($"Quad order degree: {quadOrder}, for field {A.Identification} deg = {A.Basis.Degree}, field {B.Identification} deg = {B.Basis.Degree}");


                if (A.GridDat.SpatialDimension != B.GridDat.SpatialDimension)
                    throw new ArgumentException("Both fields must have the same spatial dimension.");

                if (object.ReferenceEquals(A.GridDat, B.GridDat) && false) {
                    // ++++++++++++++
                    // equal meshes
                    // ++++++++++++++
                    CellMask domain = scheme != null ? scheme.Domain : null;

                    tr.Info("equal meshes");
                    double errPow2 = A.L2Error(B, domain).Pow2();
                    tr.Info("error^2 = " + errPow2);
                    if (IgnoreMeanValue) {
                        // domain volume
                        double Vol = 0;
                        int J = A.GridDat.iGeomCells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++) {
                            Vol += A.GridDat.iGeomCells.GetCellVolume(j);
                        }
                        Vol = Vol.MPISum();

                        // mean value
                        double mean = A.GetMeanValueTotal(domain) - B.GetMeanValueTotal(domain);
                        tr.Info("mean value = " + mean);

                        // Note: for a field p, we have 
                        // || p - <p> ||^2 = ||p||^2 - <p>^2*vol
                        return Math.Sqrt((errPow2 - mean * mean * Vol).Abs());
                    } else {
                        tr.Info("norm WITH mean value");
                        return Math.Sqrt(errPow2);
                    }

                } else {
                    // ++++++++++++++++++
                    // different meshes
                    // ++++++++++++++++++

                    tr.Info("Different meshes...");
                    DGField fine, coarse;
                    if (A.GridDat.CellPartitioning.TotalLength > B.GridDat.CellPartitioning.TotalLength) {
                        fine = A;
                        coarse = B;
                    } else {
                        fine = B;
                        coarse = A;
                    }

<<<<<<< HEAD

=======
                    
>>>>>>> b7f11eee5f (before adding new level-set evolver)
                    var CompQuadRule = scheme.SaveCompile(coarse.GridDat, quadOrder); // use over-integration
                    var eval = new FieldEvaluation(GridHelper.ExtractGridData(fine.GridDat));

                    void FineEval(MultidimensionalArray input, MultidimensionalArray output) {
                        int L = input.GetLength(0);
                        Debug.Assert(output.GetLength(0) == L);

                        eval.Evaluate(1.0, new DGField[] { fine }, input, 0.0, output.ResizeShallow(L, 1));

                    }


                    double errPow2 = coarse.LxError(FineEval, (double[] X, double fC, double fF) => (fC - fF).Pow2(), CompQuadRule, Quadrature_ChunkDataLimitOverride: int.MaxValue);
                    tr.Info("error^2 = " + errPow2);

                    if (IgnoreMeanValue == true) {

                        // domain volume
                        double Vol = 0;
                        int J = coarse.GridDat.iGeomCells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++) {
                            Vol += coarse.GridDat.iGeomCells.GetCellVolume(j);
                        }
                        Vol = Vol.MPISum();

                        // mean value times domain volume 
                        double meanVol = coarse.LxError(FineEval, (double[] X, double fC, double fF) => fC - fF, CompQuadRule, Quadrature_ChunkDataLimitOverride: int.MaxValue);
                        tr.Info("mean value * domain volume = " + meanVol + ", domain volume = " + Vol);

                        // Note: for a field p, we have 
                        // || p - <p> ||^2 = ||p||^2 - <p>^2*vol
                        return Math.Sqrt(Math.Max(0, errPow2 - meanVol * meanVol / Vol));
                    } else {
                        tr.Info("norm WITH mean value");
                        return Math.Sqrt(errPow2);
                    }
                }
            }
        }

        /// <summary>
        /// Approximate H1 distance (difference in the H1 norm) between two DG fields; this also supports DG fields on different meshes, 
        /// it could be used for convergence studies.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="scheme">
        /// a cell quadrature scheme on the coarse of the two meshes
        /// </param>
        /// <returns></returns>
        static public double H1Distance(this ConventionalDGField A, ConventionalDGField B, CellQuadratureScheme scheme = null) {
            if (A.GridDat.SpatialDimension != B.GridDat.SpatialDimension)
                throw new ArgumentException("Both fields must have the same spatial dimension.");

            int D = A.GridDat.SpatialDimension;

            double Acc = 0.0;
            Acc += L2Distance(A, B, false, scheme).Pow2();
            for (int d = 0; d < D; d++) {
                ConventionalDGField dA_dd = new SinglePhaseField(A.Basis);
                dA_dd.Derivative(1.0, A, d, scheme != null ? scheme.Domain : null);

                ConventionalDGField dB_dd = new SinglePhaseField(B.Basis);
                dB_dd.Derivative(1.0, B, d, scheme != null ? scheme.Domain : null);

                Acc += L2Distance(dA_dd, dB_dd, false, scheme).Pow2();
            }

            return Acc.Sqrt();
        }


        /// <summary>
        /// Approximate H1 norm (aka. Sobolev norm) of a DG field;
        /// </summary>
        /// <param name="A"></param>
        /// <param name="mask">
        /// a cell quadrature scheme on the coarse of the two meshes
        /// </param>
        /// <returns></returns>
        static public double H1Norm(this DGField A, CellMask mask = null) {

            int D = A.GridDat.SpatialDimension;

            double Acc = 0.0;
            Acc += A.L2Norm();
            for (int d = 0; d < D; d++) {
                DGField dA_dd = A.CloneAs();
                dA_dd.Clear();
                dA_dd.Derivative(1.0, A, d, mask);

                Acc += dA_dd.L2Norm(mask).Pow2();
            }

            return Acc.Sqrt();
        }


        /// <summary>
        /// L2 norm of a DG field, without the mean value (i.e. typically a norm used for hydrodynamic pressure, where the mean value is not relevant).
        /// I.e., for some field $`f `$, we compute $ \left| f -  \langle f \rangle \right| $, where $` \langle f \rangle `$ is the average value;
        /// Therefore, we have the relation
        /// ```math
        ///   \left| f -  \langle f \rangle \right|^2 = \left| f \right|^2 - \langle f \rangle^2 | \Omega |
        /// ```
        /// </summary>
        static public double L2Norm_IgnoreMean(this DGField A, CellMask domain = null) {

            if (domain == null)
                domain = CellMask.GetFullMask(A.GridDat);

            // L2 norm
            double errPow2 = A.L2Norm(domain).Pow2();

            // domain volume
            double Vol = 0;
            foreach (int j in domain.ItemEnum) {
                Vol += A.GridDat.iGeomCells.GetCellVolume(j);
            }
            Vol = Vol.MPISum();

            // mean value
            double mean = A.GetMeanValueTotal(domain);

            // Note: for a field p, we have 
            // || p - <p> ||^2 = ||p||^2 - <p>^2*vol
            return Math.Sqrt(Math.Max(0, errPow2 - mean * mean * Vol));

        }
    }
}
