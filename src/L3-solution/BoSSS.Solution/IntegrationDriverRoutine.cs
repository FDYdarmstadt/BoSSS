using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.Utils {
    
    /// <summary>
    /// Driver routines for integral evaluation
    /// </summary>
    public static class IntegrationDriverRoutine {


        /// <summary>
        /// Resolves the quadrature infrastructure (grid data, scheme helper) and
        /// determines the number of components of the integrand <paramref name="func"/>.
        /// </summary>
        /// <param name="lsTrk">the level-set tracker</param>
        /// <param name="quadOrder">quadrature order</param>
        /// <param name="U">DG fields the integrand may depend on</param>
        /// <param name="func">the integrand (normal-dependent for surface integrals, otherwise null)</param>
        /// <param name="funcSurf">the surface integrand; exactly one of <paramref name="func"/>/<paramref name="funcSurf"/> must be non-null</param>
        static int GetComponentCount(LevelSetTracker lsTrk, int quadOrder, DGField[] U,
                                     Func<Vector, double[], double[]> func,
                                     Func<Vector, Vector, double[], double[]> funcSurf) {
            int D = lsTrk.GridDat.SpatialDimension;
            var x = new Vector(D);
            var u = new double[U?.Length ?? 0];

            double[] probe;
            if(func != null)
                probe = func(x, u);
            else if(funcSurf != null)
                probe = funcSurf(x, x, u);
            else
                throw new ArgumentException("no integrand given");

            if(probe == null || probe.Length < 1)
                throw new ArgumentException("integrand must return at least one component");
            return probe.Length;
        }


        /// <summary>
        /// scalar version of <see cref="VolumeIntegral(LevelSetTracker, SpeciesId, Func{Vector, double[], double[]}, CellMask, int, int, DGField[])"/>
        /// </summary>
        public static double[] VolumeIntegral(this LevelSetTracker lsTrk, SpeciesId Spc, Func<Vector, double[], double> func, CellMask cm = null, int HistoryIndex = 1, int quadOrder = -1, DGField[] U = null) {
            var res = VolumeIntegral(lsTrk, Spc,
                (Vector X, double[] U) => new double[] { func(X, U) },
                cm:cm, HistoryIndex:HistoryIndex, quadOrder:quadOrder, U:U);
            return res.ExtractSubArrayShallow(-1, 0).To1DArray();
        }



        /// <summary>
        /// Driver interface for numerical integration;
        /// 
        /// For each cell $ K_j $ in <paramref name="cm"/>, this computes
        /// \[
        ///   I_j = \int_{K_j \cap S} f(\vec{x}, u_1(\vec{x}), ... , u_n(\vec{x}) ) dV 
        /// \]
        /// where:
        /// - $ S $ is the domain of the species <paramref name="Spc"/>
        /// - the position $ \vec{x} $
        /// - $ u_i $ are the DG fields provided by <paramref name="U"/>
        /// - $ f(\vec{x}, u_1(\vec{x}), ... , u_n(\vec{x}) ) $ is the function <paramref name="func"/>
        /// <param name="lsTrk">
        /// required XDG tracker <see cref="LevelSetTracker.GetXDGSpaceMetrics(int, int)"/>
        /// </param>
        /// </summary>
        /// <param name="Spc">
        /// species index
        /// </param>
        /// <param name="cm">
        /// optional cell-mask; if null, the entire domain is taken
        /// </param>
        /// <param name="func">
        /// the integrand
        /// </param>
        /// <param name="quadOrder">
        /// quadrature order
        /// </param>
        /// <param name="U">
        /// DG fields on which <paramref name="func"/> might depend on
        /// </param>
        /// <param name="HistoryIndex">
        /// time level, 1 is current
        /// </param>
        /// <returns>
        /// the integrand value for each cell
        /// - 1st index: corresponds to the cell index in <paramref name="cm"/>
        /// - 2nd index: integrand component, i.e., correlates with the return values of <paramref name="func"/>
        /// </returns>
        public static MultidimensionalArray VolumeIntegral(this LevelSetTracker lsTrk, SpeciesId Spc, Func<Vector, double[], double[]> func, CellMask cm = null, int HistoryIndex = 1, int quadOrder = -1, DGField[] U = null) {
            var gd = lsTrk.GridDat;
            XDGSpaceMetrics metrics = GetMetrics(lsTrk, HistoryIndex, ref quadOrder, ref U);
            var schH = metrics.XQuadSchemeHelper;
            int J = gd.iLogicalCells.Count;
            int nComp = GetComponentCount(lsTrk, quadOrder, U, func, null);
            int nFields = U?.Length ?? 0;

            // integration domain: species domain, optionally restricted by 'cm'
            // (GetVolumeQuadScheme intersects the species mask with 'cm' internally)
            var volRule = schH.GetVolumeQuadScheme(Spc, true, cm).Compile(gd, quadOrder);

            var result = MultidimensionalArray.Create(J, nComp);

            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { nComp }, gd,
                volRule,
                _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    int K = QR.NoOfNodes;
                    var physNodes = gd.GlobalNodes.GetValue_Cell(QR.Nodes, i0, Length);

                    var fieldVals = new MultidimensionalArray[nFields];
                    for(int f = 0; f < nFields; f++) {
                        fieldVals[f] = MultidimensionalArray.Create(Length, K);
                        U[f].Evaluate(i0, Length, QR.Nodes, fieldVals[f], 0, 0.0);
                    }

                    var x = new Vector(gd.SpatialDimension);
                    var u = new double[nFields];
                    for(int i = 0; i < Length; i++) {
                        for(int k = 0; k < K; k++) {
                            for(int d = 0; d < x.Dim; d++)
                                x[d] = physNodes[i, k, d];
                            for(int f = 0; f < nFields; f++)
                                u[f] = fieldVals[f][i, k];

                            double[] fv = func(x, u);
                            for(int c = 0; c < nComp; c++)
                                EvalResult[i, k, c] = fv[c];
                        }
                    }
                },
                _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        int jCell = i + i0;
                        for(int c = 0; c < nComp; c++)
                            result[jCell, c] += ResultsOfIntegration[i, c];
                    }
                }).Execute();

            return result;
        }

        private static XDGSpaceMetrics GetMetrics(LevelSetTracker lsTrk, int HistoryIndex, ref int quadOrder, ref DGField[] U) {
            if(U == null)
                U = [];
            if(quadOrder < 0) {
                quadOrder = 2; // some reasonable minimum
                foreach(var u in U) {
                    quadOrder = Math.Max(u.Basis.Degree * 2, quadOrder);
                }

                var _q = quadOrder;
                var cached = lsTrk.GetCachedOrders(HistoryIndex).Where(p => p >= _q);
                if(cached.Any()) {
                    quadOrder = cached.Min();
                }
            }
            var metrics = lsTrk.GetXDGSpaceMetrics(CutCellsQuadOrder: quadOrder, HistoryIndex: HistoryIndex);
            return metrics;
        }

        /// <summary>
        /// scalar version of <see cref="CellSurfaceIntegral(LevelSetTracker, SpeciesId, Func{Vector, Vector, double[], double[]}, CellMask, int, int, DGField[])"/>
        /// </summary>
        public static double[] CellSurfaceIntegral(this LevelSetTracker lsTrk, SpeciesId Spc, Func<Vector, Vector, double[], double> func, CellMask cm = null, int HistoryIndex = 1, int quadOrder = -1, DGField[] U = null) {
            var res = CellSurfaceIntegral(lsTrk, Spc, 
                (Vector X, Vector N, double[] U) => new double[] { func(X, N, U) },
                cm: cm, HistoryIndex: HistoryIndex, quadOrder: quadOrder, U: U);
            return res.ExtractSubArrayShallow(-1, 0).To1DArray();
        }


        /// <summary>
        /// Driver interface for numerical integration;
        /// 
        /// For each cell $ K_j $ in <paramref name="cm"/>, this computes
        /// \[
        ///   I_j = \oint_{\partial K_j \cap S} f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) dV 
        /// \]
        /// where:
        /// - $ S $ is the domain of the species <paramref name="Spc"/>
        /// - the position $ \vec{x} $
        /// - the outer cell normal $ \vec{n} $
        /// - $ u_i $ are the DG fields provided by <paramref name="U"/>
        /// - $ f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) $ is the function <paramref name="func"/>, which might depend on the position $ \vec{x} $ and the values of the DG fields
        /// - Note: the dependence of $ f $ on $ \vec{n} $ must be linear, i.e., $ f(\vec{x}, -\vec{n}, \ldots ) =  -f(\vec{x}, \vec{n}, \ldots ) $
        /// </summary>
        /// <param name="lsTrk">
        /// required XDG tracker <see cref="LevelSetTracker.GetXDGSpaceMetrics(int, int)"/>
        /// </param>
        /// <param name="Spc">
        /// species index
        /// </param>
        /// <param name="cm">
        /// optional cell-mask; if null, the entire domain is taken
        /// </param>
        /// <param name="func">
        /// the integrand
        /// </param>
        /// <param name="quadOrder">
        /// quadrature order
        /// </param>
        /// <param name="U">
        /// DG fields on which <paramref name="func"/> might depend on
        /// </param>
        /// <param name="HistoryIndex">
        /// time level, 1 is current
        /// </param>
        /// <returns>
        /// the integrand value for each cell
        /// - 1st index: corresponds to the cell index in <paramref name="cm"/>
        /// - 2nd index: integrand component, i.e., correlates with the return values of <paramref name="func"/>
        /// </returns>
        public static MultidimensionalArray CellSurfaceIntegral(this LevelSetTracker lsTrk, SpeciesId Spc, Func<Vector, Vector, double[], double[]> func, CellMask cm = null, int HistoryIndex = 1, int quadOrder = -1, DGField[] U = null) {
            var gd = lsTrk.GridDat;
            XDGSpaceMetrics metrics = GetMetrics(lsTrk, HistoryIndex, ref quadOrder, ref U);
            var schH = metrics.XQuadSchemeHelper;
            int J = gd.iLogicalCells.Count;
            int nComp = GetComponentCount(lsTrk, quadOrder, U, null, func);
            int nFields = U?.Length ?? 0;

            int[,] E2C = gd.iLogicalEdges.CellIndices;
            var Edge2CellTrafoIndex = gd.iGeomEdges.Edge2CellTrafoIndex;

            // optional cell-mask: only accumulate contributions for cells in 'cm'
            System.Collections.BitArray cmMask = (cm != null) ? cm.GetBitMaskWithExternal() : null;

            // edge domain: all inner + domain-boundary edges of the species subgrid
            var edgeScheme = schH.GetEdgeQuadScheme(Spc);
            var edgeRule = edgeScheme.Compile(gd, quadOrder);

            var result = MultidimensionalArray.Create(J, nComp);

            BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                new int[] { nComp }, gd,
                edgeRule,
                _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    int K = QR.NoOfNodes;
                    int D = gd.SpatialDimension;

                    // outward normals w.r.t. the IN-cell of each edge
                    var normals = gd.iGeomEdges.NormalsCache.GetNormals_Edge(QR.Nodes, i0, Length);
                    var physNodes = gd.GlobalNodes.GetValue_EdgeSV(QR.Nodes, i0, Length);

                    var x = new Vector(D);
                    var n = new Vector(D);
                    var u = new double[nFields];

                    for(int i = 0; i < Length; i++) {
                        int iEdge = i0 + i;
                        int jCellIn = E2C[iEdge, 0];
                        int iTrf = Edge2CellTrafoIndex[iEdge, 0];

                        // reference nodes within the IN-cell, to evaluate DG fields
                        NodeSet cellNodes = QR.Nodes.GetVolumeNodeSet(gd, iTrf, true);

                        var fieldVals = new MultidimensionalArray[nFields];
                        for(int f = 0; f < nFields; f++) {
                            fieldVals[f] = MultidimensionalArray.Create(1, K);
                            U[f].Evaluate(jCellIn, 1, cellNodes, fieldVals[f], 0, 0.0);
                        }

                        for(int k = 0; k < K; k++) {
                            for(int d = 0; d < D; d++) {
                                x[d] = physNodes[i, k, d];
                                n[d] = normals[i, k, d];
                            }
                            for(int f = 0; f < nFields; f++)
                                u[f] = fieldVals[f][0, k];

                            double[] fv = func(x, n, u);
                            for(int c = 0; c < nComp; c++)
                                EvalResult[i, k, c] = fv[c];
                        }
                    }
                },
                _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        int iEdge = i0 + i;
                        int jCell0 = E2C[iEdge, 0]; // IN-cell:  normal is the outward normal
                        int jCell1 = E2C[iEdge, 1]; // OUT-cell: outward normal is negated

                        for(int c = 0; c < nComp; c++) {
                            double val = ResultsOfIntegration[i, c];
                            if(cmMask == null || cmMask[jCell0])
                                result[jCell0, c] += val;
                            if(jCell1 >= 0 && (cmMask == null || cmMask[jCell1]))
                                result[jCell1, c] -= val;
                        }
                    }
                }).Execute();

            return result;
        }

        /// <summary>
        /// scalar version of <see cref="InterfaceIntegral(LevelSetTracker, int, SpeciesId, SpeciesId, Func{Vector, Vector, double[], double[]}, CellMask, int, int, DGField[])"/>
        /// </summary>
        public static double[] InterfaceIntegral(this LevelSetTracker lsTrk, int iLevSet, SpeciesId SpcA, SpeciesId SpcB, Func<Vector, Vector, double[], double> func, CellMask cm = null, int HistoryIndex = 1, int quadOrder = -1, DGField[] U = null) {
            var res = InterfaceIntegral(lsTrk, iLevSet, SpcA, SpcB, 
                (Vector X, Vector N, double[] U) => new double[] { func(X, N, U) },
                cm: cm, HistoryIndex: HistoryIndex, quadOrder: quadOrder, U: U);
            return res.ExtractSubArrayShallow(-1, 0).To1DArray();
        }


        /// <summary>
        /// Driver interface for numerical integration;
        /// 
        /// For each cell $ K_j $ in <paramref name="cm"/>, this computes
        /// \[
        ///   I_j = \oint_{\mathfrak{I} \cap \overline{S_1} \cap \overline{S_2}} f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) dV 
        /// \]
        /// where:
        /// - $ S_1 $ is the domain of the negative/first species <paramref name="SpcA"/>
        /// - $ S_2 $ is the domain of the positive/second species <paramref name="SpcB"/>
        /// - the position $ \vec{x} $
        /// - the level-set/interface normal $ \vec{n} $, pointing from species <paramref name="SpcA"/> to <paramref name="SpcB"/>
        /// - $ u_i $ are the DG fields provided by <paramref name="U"/>
        /// - $ f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) $ is the function <paramref name="func"/>, which might depend on the position $ \vec{x} $ and the values of the DG fields
        /// - Note: the dependence of $ f $ on $ \vec{n} $ must be linear, i.e., $ f(\vec{x}, -\vec{n}, \ldots ) =  -f(\vec{x}, \vec{n}, \ldots ) $
        /// </summary>
        /// <param name="lsTrk">
        /// required XDG tracker <see cref="LevelSetTracker.GetXDGSpaceMetrics(int, int)"/>
        /// </param>
        /// <param name="iLevSet">
        /// level-set/interface index
        /// </param>
        /// <param name="SpcA">
        /// first species index
        /// </param>
        /// <param name="SpcB">
        /// second species index
        /// </param>
        /// <param name="cm">
        /// optional cell-mask; if null, the entire domain is taken
        /// </param>
        /// <param name="func">
        /// the integrand
        /// </param>
        /// <param name="quadOrder">
        /// quadrature order
        /// </param>
        /// <param name="U">
        /// DG fields on which <paramref name="func"/> might depend on
        /// </param>
        /// <param name="HistoryIndex">
        /// time level, 1 is current
        /// </param>
        /// <returns>
        /// the integrand value for each cell
        /// - 1st index: corresponds to the cell index in <paramref name="cm"/>
        /// - 2nd index: integrand component, i.e., correlates with the return values of <paramref name="func"/>
        /// </returns>
        public static MultidimensionalArray InterfaceIntegral(this LevelSetTracker lsTrk, int iLevSet, SpeciesId SpcA, SpeciesId SpcB, Func<Vector, Vector, double[], double[]> func, CellMask cm = null, int HistoryIndex = 1, int quadOrder = -1, DGField[] U = null) {
            var gd = lsTrk.GridDat;
            XDGSpaceMetrics metrics = GetMetrics(lsTrk, HistoryIndex, ref quadOrder, ref U);
            var schH = metrics.XQuadSchemeHelper;
            int J = gd.iLogicalCells.Count;
            int nComp = GetComponentCount(lsTrk, quadOrder, U, null, func);
            int nFields = U?.Length ?? 0;

            // integration domain: interface of level-set 'iLevSet' between species A and B,
            // optionally restricted by 'cm'
            CellMask integrationDom = metrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet)
                .Intersect(metrics.LevelSetRegions.GetSpeciesMask(SpcA))
                .Intersect(metrics.LevelSetRegions.GetSpeciesMask(SpcB));
            if(cm != null)
                integrationDom = integrationDom.Intersect(cm);

            var surfScheme = schH.GetLevelSetQuadScheme(iLevSet, SpcA, SpcB, integrationDom);
            var surfRule = surfScheme.Compile(gd, quadOrder);

            var result = MultidimensionalArray.Create(J, nComp);

            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { nComp }, gd,
                surfRule,
                _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    int K = QR.NoOfNodes;
                    int D = gd.SpatialDimension;

                    // interface normal, pointing from species A (level-set < 0) to species B (level-set > 0)
                    var normals = lsTrk.DataHistories[iLevSet].Current.GetLevelSetNormals(QR.Nodes, i0, Length);
                    var physNodes = gd.GlobalNodes.GetValue_Cell(QR.Nodes, i0, Length);

                    var fieldVals = new MultidimensionalArray[nFields];
                    for(int f = 0; f < nFields; f++) {
                        fieldVals[f] = MultidimensionalArray.Create(Length, K);
                        U[f].Evaluate(i0, Length, QR.Nodes, fieldVals[f], 0, 0.0);
                    }

                    var x = new Vector(D);
                    var n = new Vector(D);
                    var u = new double[nFields];
                    for(int i = 0; i < Length; i++) {
                        for(int k = 0; k < K; k++) {
                            for(int d = 0; d < D; d++) {
                                x[d] = physNodes[i, k, d];
                                n[d] = normals[i, k, d];
                            }
                            for(int f = 0; f < nFields; f++)
                                u[f] = fieldVals[f][i, k];

                            double[] fv = func(x, n, u);
                            for(int c = 0; c < nComp; c++)
                                EvalResult[i, k, c] = fv[c];
                        }
                    }
                },
                _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++) {
                        int jCell = i + i0;
                        for(int c = 0; c < nComp; c++)
                            result[jCell, c] += ResultsOfIntegration[i, c];
                    }
                }).Execute();

            return result;
        }
    }
}
