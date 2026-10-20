using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.CurvedElementTests {

    /// <summary>
    /// Shared helpers for tests on curved elements with rotationally periodic boundaries (see <see cref="TaylorCouetteTests"/>):
    /// grid, exact Taylor-Couette solution, and the velocity jump across (rotated) periodic and interior edges.
    /// </summary>
    public static class CurvedWedgeHelpers {

        /// <summary>
        /// Curved 2D annulus grid, 1 &lt;= r &lt;= 2:
        /// either a wedge which is rotationally periodic in the azimuthal direction,
        /// or the full ring, without any periodic boundary (the last cell in azimuthal direction is a conformal neighbor of the first one).
        /// </summary>
        /// <param name="resolution">
        /// refinement level: (4 x 10) * <paramref name="resolution"/> cells in (r, theta) direction for the wedge;
        /// the full ring has the same cell size, i.e. (4 x 54) * <paramref name="resolution"/> cells
        /// </param>
        /// <param name="fullRing">false: periodic wedge, true: full ring</param>
        internal static GridCommons CreateGrid(int resolution = 1, bool fullRing = false) {
            double[] rNodes = GenericBlas.Linspace(1, 2, 4 * resolution + 1);
            double[] sNodes = fullRing
                ? GenericBlas.Linspace(0, 1, 54 * resolution + 1)
                : GenericBlas.Linspace(0, 0.187, 10 * resolution + 1);
            return Grid2D.CurvedSquareGrid(rNodes, sNodes, CellType.Square_16, PeriodicS: !fullRing);
        }

        /// <summary>
        /// Exact solution of laminar Taylor-Couette flow between two concentric cylinders:
        /// azimuthal velocity \f$ u_\theta(r) = A r + B / r \f$, with
        /// \f$ A = (\Omega_o r_o^2 - \Omega_i r_i^2)/(r_o^2 - r_i^2) \f$ and
        /// \f$ B = (\Omega_i - \Omega_o) r_i^2 r_o^2 / (r_o^2 - r_i^2) \f$.
        /// </summary>
        /// <param name="d">Cartesian velocity component (0: x, 1: y)</param>
        /// <param name="r_i">inner radius</param>
        /// <param name="r_o">outer radius</param>
        /// <param name="Omega_i">angular velocity of the inner cylinder</param>
        /// <param name="Omega_o">angular velocity of the outer cylinder</param>
        internal static ScalarFunction TaylorCouetteVelocity(int d, double r_i, double r_o, double Omega_i, double Omega_o) {
            var (A, B) = TaylorCouetteCoefficients(r_i, r_o, Omega_i, Omega_o);

            return delegate (MultidimensionalArray X, MultidimensionalArray U) {
                int K = X.GetLength(0);
                for(int k = 0; k < K; k++) {
                    double x = X[k, 0], y = X[k, 1];
                    double r = Math.Sqrt(x * x + y * y);
                    double uTheta = A * r + B / r;
                    // e_theta = (-y/r, x/r)
                    U[k] = d == 0 ? -uTheta * y / r : uTheta * x / r;
                }
            };
        }

        /// <summary>
        /// Coefficients \f$ A, B \f$ of the Taylor-Couette velocity \f$ u_\theta(r) = A r + B / r \f$, see <see cref="TaylorCouetteVelocity"/>.
        /// </summary>
        static (double A, double B) TaylorCouetteCoefficients(double r_i, double r_o, double Omega_i, double Omega_o) {
            double A = (Omega_o * r_o * r_o - Omega_i * r_i * r_i) / (r_o * r_o - r_i * r_i);
            double B = (Omega_i - Omega_o) * r_i * r_i * r_o * r_o / (r_o * r_o - r_i * r_i);
            return (A, B);
        }

        /// <summary>
        /// Exact pressure of laminar Taylor-Couette flow (up to a constant), from \f$ \partial_r p = \rho u_\theta^2 / r \f$:
        /// \f$ p(r) = \rho ( A^2 r^2 / 2 + 2 A B \ln r - B^2 / (2 r^2) ) \f$, see <see cref="TaylorCouetteVelocity"/> for \f$ A, B \f$;
        /// for Stokes flow (no convection), use \f$ \rho = 0 \f$, i.e. constant pressure.
        /// </summary>
        /// <param name="rho">density</param>
        /// <param name="r_i">inner radius</param>
        /// <param name="r_o">outer radius</param>
        /// <param name="Omega_i">angular velocity of the inner cylinder</param>
        /// <param name="Omega_o">angular velocity of the outer cylinder</param>
        internal static ScalarFunction TaylorCouettePressure(double rho, double r_i, double r_o, double Omega_i, double Omega_o) {
            var (A, B) = TaylorCouetteCoefficients(r_i, r_o, Omega_i, Omega_o);

            return delegate (MultidimensionalArray X, MultidimensionalArray P) {
                int K = X.GetLength(0);
                for(int k = 0; k < K; k++) {
                    double x = X[k, 0], y = X[k, 1];
                    double r2 = x * x + y * y;
                    P[k] = rho * (0.5 * A * A * r2 + A * B * Math.Log(r2) - 0.5 * B * B / r2); // note: 2 ln(r) = ln(r^2)
                }
            };
        }

        /// <summary>
        /// L2 norms of the velocity jump over periodic and interior edges
        /// </summary>
        internal struct JumpNorms {
            /// <summary> jump across periodic edges, OUT-velocity rotated into the IN-frame </summary>
            public double Periodic;
            /// <summary> jump across periodic edges, without rotation </summary>
            public double PeriodicUnrotated;
            /// <summary> jump across interior (non-periodic) edges </summary>
            public double Interior;
            /// <summary> total length (area in 3D) of periodic edges </summary>
            public double LengthPeriodic;
            /// <summary> total length (area in 3D) of interior edges </summary>
            public double LengthInterior;
        }

        /// <summary>
        /// Computes the L2 norm of the jump \f$ \vec{u}_{IN} - \mathbf{M} \vec{u}_{OUT} \f$ over all inner edges,
        /// where \f$ \mathbf{M} \f$ is the matrix of the periodic transformation which maps the OUT-side of an edge onto the IN-side
        /// (identity for non-periodic edges).
        /// </summary>
        /// <param name="u">velocity components</param>
        /// <param name="quadOrder">quadrature order</param>
        internal static JumpNorms VelocityJumps(SinglePhaseField[] u, int quadOrder) {
            var gdat = u[0].GridDat;
            int D = gdat.SpatialDimension;
            if(u.Length != D)
                throw new ArgumentException("expecting one field per spatial direction");

            int[,] CellIndices = gdat.iGeomEdges.CellIndices;
            int[,] E2C = gdat.iGeomEdges.Edge2CellTrafoIndex;
            byte[] EdgeTags = gdat.iGeomEdges.EdgeTags;
            int Jup = gdat.iGeomCells.NoOfLocalUpdatedCells;

            // split the inner edges into periodic and interior (non-periodic) ones
            int E = gdat.iGeomEdges.Count;
            var periodicEdges = new BitArray(E);
            var interiorEdges = new BitArray(E);
            for(int e = 0; e < E; e++) {
                if(CellIndices[e, 1] < 0)
                    continue; // boundary edge
                if(EdgeTags[e] >= GridCommons.FIRST_PERIODIC_BC_TAG)
                    periodicEdges[e] = true;
                else
                    interiorEdges[e] = true;
            }
            var periodicMask = new EdgeMask(gdat, periodicEdges);
            var interiorMask = new EdgeMask(gdat, interiorEdges);

            // same convention as in DGField.JumpNorm: count an inter-process edge only on one process
            bool CountEdge(int iEdge) {
                int jCell1 = CellIndices[iEdge, 1];
                if(jCell1 < Jup)
                    return true;
                long Gidx0 = gdat.GetLogicalCellIndex(CellIndices[iEdge, 0]) + gdat.CellPartitioning.i0;
                long Gidx1 = gdat.iParallel.GlobalIndicesExternalCells[gdat.GetLogicalCellIndex(jCell1) - Jup];
                return Gidx0 < Gidx1;
            }

            EdgeQuadratureScheme Scheme(EdgeMask m) => new EdgeQuadratureScheme(false, m).AddFixedOrderRules(gdat, quadOrder);

            // plain (un-rotated) jump of the velocity vector, from the jumps of the individual components
            double PlainJump(EdgeMask m) => Math.Sqrt(u.Sum(f => f.JumpNorm(Scheme(m)).Pow2()));

            double TotalLength(EdgeMask m) => m.ItemEnum.Where(CountEdge).Sum(e => gdat.iGeomEdges.GetEdgeArea(e)).MPISum();

            // periodic transformations and their inverses: one of them maps the OUT-side of a periodic edge onto the IN-side
            var Candidates = new List<AffineTrafo>();
            foreach(var T in gdat.Grid.PeriodicTrafo) {
                Candidates.Add(T);
                Candidates.Add(T.Invert());
            }

            foreach(var f in u)
                f.MPIExchange();

            // rotated jump across periodic edges: no built-in function for that
            double accRot = 0.0;
            EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, gdat,
                Scheme(periodicMask).Compile(gdat, quadOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    int K = NS.NoOfNodes;

                    var uIN = new MultidimensionalArray[D];
                    var uOT = new MultidimensionalArray[D];
                    for(int d = 0; d < D; d++) {
                        uIN[d] = MultidimensionalArray.Create(Length, K);
                        uOT[d] = MultidimensionalArray.Create(Length, K);
                        u[d].EvaluateEdge(i0, Length, NS, uIN[d], uOT[d], null, null, null, null, 0, 0.0);
                    }

                    var X0 = MultidimensionalArray.Create(K, D);
                    var X1 = MultidimensionalArray.Create(K, D);
                    var X1mapped = MultidimensionalArray.Create(K, D);
                    var jump = new double[D];

                    for(int i = 0; i < Length; i++) {
                        int iEdge = i + i0;

                        // find the transformation which maps the OUT-side quadrature nodes onto the IN-side ones
                        gdat.TransformLocal2Global(NS.GetVolumeNodeSet(gdat, E2C[iEdge, 0], false), X0, CellIndices[iEdge, 0]);
                        gdat.TransformLocal2Global(NS.GetVolumeNodeSet(gdat, E2C[iEdge, 1], false), X1, CellIndices[iEdge, 1]);

                        MultidimensionalArray M = null;
                        double minDist = double.MaxValue;
                        foreach(var T in Candidates) {
                            T.Transform(X1, X1mapped);
                            X1mapped.Acc(-1.0, X0);
                            double dist = X1mapped.InfNorm();
                            if(dist < minDist) {
                                minDist = dist;
                                M = T.Matrix;
                            }
                        }
                        Assert.Less(minDist, 1.0e-8, $"no periodic transformation maps the OUT-side of edge {iEdge} onto its IN-side.");

                        // jump u_IN - M*u_OUT
                        for(int k = 0; k < K; k++) {
                            for(int d = 0; d < D; d++) {
                                jump[d] = uIN[d][i, k];
                                for(int dd = 0; dd < D; dd++)
                                    jump[d] -= M[d, dd] * uOT[dd][i, k];
                            }
                            EvalResult[i, k, 0] = jump.L2NormPow2();
                        }
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for(int i = 0; i < Length; i++) {
                        if(CountEdge(i + i0))
                            accRot += ResultsOfIntegration[i, 0];
                    }
                }).Execute();

            return new JumpNorms() {
                Periodic = Math.Sqrt(accRot.MPISum()),
                PeriodicUnrotated = PlainJump(periodicMask),
                LengthPeriodic = TotalLength(periodicMask),
                Interior = PlainJump(interiorMask),
                LengthInterior = TotalLength(interiorMask)
            };
        }
    }
}
