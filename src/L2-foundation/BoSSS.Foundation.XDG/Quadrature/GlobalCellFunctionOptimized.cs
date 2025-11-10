using System;
using System.Collections.Generic;
using System.Text;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using ilPSP;
using IntersectingQuadrature.Tensor;
using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {



    /// <summary>
    /// Fast in-cell function using our own Chebyshev tensor polynomial.
    /// One-time: sample BoSSS at a tensor Chebyshev–Lobatto grid, compute coefficients by separable DCT-I.
    /// Runtime: evaluate by separable Clenshaw-style contractions, incl. gradient and Hessian.
    /// </summary>
    class GlobalCellFunctionOptimized : IScalarFunction {

        readonly LevelSetTracker.LevelSetData levelSet;
        readonly RefElement Kref;
        readonly int jCell;

        // Polynomial sizes per axis (n = degree+1)
        readonly int[] n;
        // Chebyshev coefficients (tensor): c[k1], c[k1,k2], or c[k1,k2,k3]
        // Store dense double[] with row-major layout for cache-friendly contractions.
        readonly double[] coeff;    // flattened
        readonly int[] stride;      // stride for index mapping

        // perf: precompute Chebyshev-Lobatto nodes per axis
        readonly double[][] nodesCL;
        // dimension
        public int M { get; private set; }

        // timers/counters (kept from your API)
        public int EvalCounterV = 0;
        public int EvalCounterG = 0;
        public int EvalCounterH = 0;
        public readonly Stopwatch stwV = new Stopwatch();
        public readonly Stopwatch stwG = new Stopwatch();
        public readonly Stopwatch stwH = new Stopwatch();


        //readonly GlobalCellFunction reference; 


        /// <summary></summary>
        /// <param name="nPerAxis">
        /// Polynomial grid sizes per axis (n = degree+1).
        /// If null, you must provide it elsewhere; setup cost is paid once.
        /// </param>
        /// <param name="jCell"></param>
        /// <param name="levelSet"></param>
        public GlobalCellFunctionOptimized(LevelSetTracker.LevelSetData levelSet, int jCell, int[] nPerAxis) {
            this.levelSet = levelSet;
            this.jCell = jCell;
            this.M = levelSet.GridDat.SpatialDimension;
            this.Kref = levelSet.GridDat.iGeomCells.GetRefElement(jCell);
            //reference = new GlobalCellFunction(levelSet, jCell);


            if(nPerAxis == null || nPerAxis.Length != M)
                throw new ArgumentException("nPerAxis must be provided with length equal to spatial dimension M.");

            // Sanity: tensor-product reference only
            if(!IsTensorProductRefElement(Kref))
                throw new NotSupportedException("This fast path assumes a tensor-product reference element ([-1,1]^d).");

            n = (int[])nPerAxis.Clone();
            stride = new int[M];
            int total = 1;
            for(int d = M - 1; d >= 0; --d) {
                stride[d] = total;
                total *= n[d];
            }
            coeff = new double[total];

            // Build per-axis Chebyshev–Lobatto nodes
            nodesCL = new double[M][];
            for(int d = 0; d < M; ++d)
                nodesCL[d] = ChebLobatto(n[d]);

            // 1) Build tensor NodeSet of all CL nodes
            NodeSet X = BuildTensorNodeSet(nodesCL);

            // 2) Sample BoSSS once (batched)
            //    Shape from BoSSS is typically [K, 1]; read as V[i]
            MultidimensionalArray V = levelSet.GetLevSetValues(X, jCell, 1);

            // 3) Move samples into a dense tensor 'F' (flattened)
            double[] F = new double[total];
            for(int k = 0; k < total; ++k)
                F[k] = V[0, k];

            // 4) Separable DCT-I along each axis to produce Chebyshev coefficients
            //    In-place transforms on 'F' into 'coeff'.
            Array.Copy(F, coeff, total);
            for(int d = 0; d < M; ++d)
                DCT1_AlongAxis_InPlace(coeff, n, stride, d);
            NormalizeChebyshevCoeffsInPlace();
        }

        // --------------------------
        // Public evaluation routines
        // --------------------------

        public double Evaluate(Tensor1 x) {
            stwV.Start();
            // If x is in global coords: map to reference here (affine) before use.
            double val = (M == 1) ? Eval1D(x) : (M == 2) ? Eval2D(x) : Eval3D(x);
            EvalCounterV++;
            stwV.Stop();

            /*
            var valRef = reference.Evaluate(x);
            if(Math.Abs(valRef - val) > 1.0e-10)
                throw new ArithmeticException("optimized version sux");
            */
            return val;
        }

        public (double evaluation, Tensor1 gradient) EvaluateAndGradient(Tensor1 x) {
            stwG.Start();
            double val;
            var grad = Tensor1.Zeros(M);

            if(M == 1) {
                EvalGrad1D(x, out val, out double gx);
                grad[0] = gx;
            } else if(M == 2) {
                EvalGrad2D(x, out val, out double gx, out double gy);
                grad[0] = gx; grad[1] = gy;
            } else {
                EvalGrad3D(x, out val, out double gx, out double gy, out double gz);
                grad[0] = gx; grad[1] = gy; grad[2] = gz;
            }

            EvalCounterG++;
            stwG.Stop();

            /*
            var valgradRef = reference.EvaluateAndGradient(x);
            if(Math.Abs(valgradRef.evaluation - val) > 1.0e-10) {
                
                Tensor1 x2 = Tensor1.Zeros(M);
                x2[0] = 0.1234;
                x2[1] = -0.3467;

                EvalGrad2D(x2, out double val2, out double gx2, out double gy2);

                var ref2 = reference.EvaluateAndGradient(x2);




                throw new ArithmeticException("optimized version sux 1");
            }
            var gradDist = (valgradRef.gradient - grad);
            double abs_gradDist = 0;
            for(int i = 0; i < M; i++) {
                abs_gradDist += gradDist[i].Abs();
            }
            if(abs_gradDist > 1.0e-10) {

                Tensor1 x2 = Tensor1.Zeros(M);
                x2[0] = 0.1234;
                x2[1] = -0.3467;

                EvalGrad2D(x2, out double val2, out double gx2, out double gy2);

                var ref2 = reference.EvaluateAndGradient(x2);




                throw new ArithmeticException("optimized version sux 2");
            }


            */
            return (val, grad);
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            stwH.Start();
            double val;
            var grad = Tensor1.Zeros(M);
            var hess = Tensor2.Zeros(M);

            if(M == 1) {
                EvalGradHess1D(x, out val, out double gx, out double hxx);
                grad[0] = gx; hess[0, 0] = hxx;
            } else if(M == 2) {
                EvalGradHess2D(x, out val,
                               out double gx, out double gy,
                               out double hxx, out double hxy, out double hyy);
                grad[0] = gx; grad[1] = gy;
                hess[0, 0] = hxx; hess[0, 1] = hxy;
                hess[1, 0] = hxy; hess[1, 1] = hyy;
            } else {
                EvalGradHess3D(x, out val,
                               out double gx, out double gy, out double gz,
                               out double hxx, out double hxy, out double hxz,
                               out double hyy, out double hyz, out double hzz);
                grad[0] = gx; grad[1] = gy; grad[2] = gz;
                hess[0, 0] = hxx; hess[0, 1] = hxy; hess[0, 2] = hxz;
                hess[1, 0] = hxy; hess[1, 1] = hyy; hess[1, 2] = hyz;
                hess[2, 0] = hxz; hess[2, 1] = hyz; hess[2, 2] = hzz;
            }

            EvalCounterH++;
            stwH.Stop();

            /*
            var valgradRef = reference.EvaluateAndGradientAndHessian(x);
            if(Math.Abs(valgradRef.evaluation - val) > 1.0e-10)
                throw new ArithmeticException("optimized version sux 1b");
            var gradDist = (valgradRef.gradient - grad);
            double abs_gradDist = 0;
            for(int i = 0; i < M; i++) {
                abs_gradDist += gradDist[i].Abs();
            }
            if(abs_gradDist > 1.0e-10)
                throw new ArithmeticException("optimized version sux 2b");
            var hessDist = (valgradRef.hessian - hess);
            double abs_hessDist = 0;
            for(int i = 0; i < M; i++) {
                for(int j = 0; j < M; j++) {
                    abs_hessDist += hessDist[i, j].Abs();
                }
            }
            if(abs_hessDist > 1.0e-10)
                throw new ArithmeticException("optimized version sux 3b");


            */
            return (val, grad, hess);
        }

        // ---------------
        // 1D evaluation
        // ---------------

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        double Eval1D(Tensor1 x) {
            double xx = x[0]; // assumed reference coord in [-1,1]
            int n0 = n[0];

            // T_k and U_k-1 at x
            var T = new double[n0];
            var U = new double[Math.Max(1, n0 - 1)];
            Build_T_and_U(xx, n0, T, U);

            // f = sum_k c_k T_k
            double s = 0.0;
            int off = 0;
            for(int k = 0; k < n0; ++k, off += 1)
                s += coeff[off] * T[k];
            return s;
        }

        void EvalGrad1D(Tensor1 x, out double f, out double fx) {
            double xx = x[0];
            int n0 = n[0];

            var T = new double[n0];
            var U = new double[Math.Max(1, n0 - 1)];
            var dU = new double[Math.Max(1, n0 - 1)];
            Build_T_U_dU(xx, n0, T, U, dU);

            double s = 0.0, sx = 0.0;
            int off = 0;
            for(int k = 0; k < n0; ++k, off += 1) {
                double ck = coeff[off];
                s += ck * T[k];
                if(k > 0) sx += ck * k * U[k - 1];
            }
            f = s;
            fx = sx;
        }

        void EvalGradHess1D(Tensor1 x, out double f, out double fx, out double fxx) {
            double xx = x[0];
            int n0 = n[0];

            var T = new double[n0];
            var U = new double[Math.Max(1, n0 - 1)];
            var dU = new double[Math.Max(1, n0 - 1)];
            Build_T_U_dU(xx, n0, T, U, dU);

            double s = 0.0, sx = 0.0, sxx = 0.0;
            int off = 0;
            for(int k = 0; k < n0; ++k, off += 1) {
                double ck = coeff[off];
                s += ck * T[k];
                if(k > 0) {
                    double ku = k * U[k - 1];
                    sx += ck * ku;
                    sxx += ck * k * dU[k - 1];
                }
            }
            f = s; fx = sx; fxx = sxx;
        }

        // ---------------
        // 2D evaluation
        // ---------------

        /*
        double Eval2D(Tensor1 x) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            Build_T_and_U(x0, n0, T0, U0);
            Build_T_and_U(x1, n1, T1, U1);

            // Contract: s = sum_{k0,k1} c[k0,k1] T0[k0] T1[k1]
            double s = 0.0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double a = 0.0;
                for(int k1 = 0; k1 < n1; ++k1, ++idx)
                    a += coeff[idx] * T1[k1];
                s += a * T0[k0];
            }
            return s;
        }*/

        double Eval2D(Tensor1 x) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            Build_T_and_U(x0, n0, T0, U0);
            Build_T_and_U(x1, n1, T1, U1);

            double s = 0.0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double a = 0.0;
                for(int k1 = 0; k1 < n1; ++k1, ++idx)
                    a += coeff[idx] * T1[k1];
                s += a * T0[k0];
            }
            return s;
        }

        /*
        void EvalGrad2D(Tensor1 x, out double f, out double fx, out double fy) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            Build_T_and_U(x0, n0, T0, U0);
            Build_T_and_U(x1, n1, T1, U1);

            double s = 0.0, sx = 0.0, sy = 0.0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                // along axis-1
                double a = 0.0, ay = 0.0;
                for(int k1 = 0; k1 < n1; ++k1, ++idx) {
                    double c = coeff[idx];
                    a += c * T1[k1];
                    if(k1 > 0) ay += c * (k1 * U1[k1 - 1]);
                }
                s += a * T0[k0];
                if(k0 > 0) sx += (k0 * U0[k0 - 1]) * a;
                sy += T0[k0] * ay;
            }
            f = s; fx = sx; fy = sy;
        }*/

        void EvalGrad2D(Tensor1 x, out double f, out double fx, out double fy) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            Build_T_and_U(x0, n0, T0, U0);
            Build_T_and_U(x1, n1, T1, U1);

            double s = 0.0, sx = 0.0, sy = 0.0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double a = 0.0, ay = 0.0;
                for(int k1 = 0; k1 < n1; ++k1, ++idx) {
                    double c = coeff[idx];
                    a += c * T1[k1];                        // undiff axis 1: already has 1/2 on k1=0
                    if(k1 > 0) ay += c * (k1 * U1[k1 - 1]); // diff axis 1: no 1/2 (k1>0 anyway)
                }
                s += a * T0[k0];                              // undiff axis 0: has 1/2 on k0=0
                if(k0 > 0) sx += (k0 * U0[k0 - 1]) * a;      // diff axis 0: no 1/2 (k0>0)
                sy += T0[k0] * ay;                             // undiff axis 0: already 1/2 on k0=0
            }
            f = s; fx = sx; fy = sy;
        }




        /*
        void EvalGradHess2D(Tensor1 x, out double f,
                            out double fx, out double fy,
                            out double fxx, out double fxy, out double fyy) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var dU0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            var dU1 = new double[Math.Max(1, n1 - 1)];
            Build_T_U_dU(x0, n0, T0, U0, dU0);
            Build_T_U_dU(x1, n1, T1, U1, dU1);

            double s = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, syy = 0.0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double a = 0.0, ay = 0.0, ayy = 0.0;
                for(int k1 = 0; k1 < n1; ++k1, ++idx) {
                    double c = coeff[idx];
                    a += c * T1[k1];
                    if(k1 > 0) {
                        double ky = k1 * U1[k1 - 1];
                        ay += c * ky;
                        ayy += c * k1 * dU1[k1 - 1];
                    }
                }
                s += a * T0[k0];
                sy += T0[k0] * ay;
                syy += T0[k0] * ayy;

                if(k0 > 0) {
                    double kxU = k0 * U0[k0 - 1];
                    sx += kxU * a;
                    sxy += kxU * ay;
                    sxx += (k0 * dU0[k0 - 1]) * a;
                }
            }
            f = s; fx = sx; fy = sy; fxx = sxx; fxy = sxy; fyy = syy;
        }*/

        // ---------------
        // 3D evaluation
        // ---------------

        double Eval3D(Tensor1 x) {
            double x0 = x[0], x1 = x[1], x2 = x[2];
            int n0 = n[0], n1 = n[1], n2 = n[2];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            var T2 = new double[n2]; var U2 = new double[Math.Max(1, n2 - 1)];
            Build_T_and_U(x0, n0, T0, U0);
            Build_T_and_U(x1, n1, T1, U1);
            Build_T_and_U(x2, n2, T2, U2);

            // (((coeff • T2) • T1) • T0)
            double s = 0.0;
            int idx012 = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double a0 = 0.0;
                for(int k1 = 0; k1 < n1; ++k1) {
                    double a1 = 0.0;
                    for(int k2 = 0; k2 < n2; ++k2, ++idx012)
                        a1 += coeff[idx012] * T2[k2];
                    a0 += a1 * T1[k1];
                }
                s += a0 * T0[k0];
            }
            return s;
        }

        void EvalGradHess2D(Tensor1 x, out double f,
                    out double fx, out double fy,
                    out double fxx, out double fxy, out double fyy) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)]; var dU0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)]; var dU1 = new double[Math.Max(1, n1 - 1)];
            Build_T_U_dU(x0, n0, T0, U0, dU0);
            Build_T_U_dU(x1, n1, T1, U1, dU1);

            // undifferentiated axes must carry 1/2 on zero modes
            T0[0] *= 0.5;
            T1[0] *= 0.5;

            double s = 0, sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double a = 0, ay = 0, ayy = 0;
                for(int k1 = 0; k1 < n1; ++k1, ++idx) {
                    double c = coeff[idx];
                    a += c * T1[k1];                            // undiff axis 1 (has 1/2 if k1=0)
                    if(k1 > 0) {
                        ay += c * (k1 * U1[k1 - 1]);           // diff axis 1
                        ayy += c * (k1 * dU1[k1 - 1]);          // diff^2 axis 1
                    }
                }
                s += a * T0[k0];                               // undiff axis 0 (has 1/2 if k0=0)
                sy += T0[k0] * ay;                              // undiff axis 0
                syy += T0[k0] * ayy;                             // undiff axis 0

                if(k0 > 0) {
                    double kxU = k0 * U0[k0 - 1];
                    double kxDU = k0 * dU0[k0 - 1];
                    sx += kxU * a;                             // diff axis 0
                    sxy += kxU * ay;                            // diff axis 0, diff axis 1
                    sxx += kxDU * a;                             // diff^2 axis 0
                }
            }
            f = s; fx = sx; fy = sy; fxx = sxx; fxy = sxy; fyy = syy;
        }


        void EvalGrad3D(Tensor1 x, out double f,
                         out double fx, out double fy, out double fz) {
            double x0 = x[0], x1 = x[1], x2 = x[2];
            int n0 = n[0], n1 = n[1], n2 = n[2];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            var T2 = new double[n2]; var U2 = new double[Math.Max(1, n2 - 1)];
            Build_T_and_U(x0, n0, T0, U0);
            Build_T_and_U(x1, n1, T1, U1);
            Build_T_and_U(x2, n2, T2, U2);

            double s = 0.0, sx = 0.0, sy = 0.0, sz = 0.0;
            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double A = 0.0, Ay = 0.0;// Az = 0.0;
                for(int k1 = 0; k1 < n1; ++k1) {
                    double B = 0.0, Bz = 0.0;
                    for(int k2 = 0; k2 < n2; ++k2, ++idx) {
                        double c = coeff[idx];
                        B += c * T2[k2];
                        if(k2 > 0) Bz += c * (k2 * U2[k2 - 1]);
                    }
                    A += B * T1[k1];
                    Ay += Bz * T1[k1];
                    if(k1 > 0) Ay += (k1 * U1[k1 - 1]) * B; // accumulate dy from axis-1 too
                }
                s += A * T0[k0];
                sy += Ay * T0[k0];
                if(k0 > 0) sx += (k0 * U0[k0 - 1]) * A;
                sz += A * 0.0; // already included in Ay via Bz; add direct z-term below
            }

            // We missed direct z contribution multiplied with T0*T1; fold it properly:
            // More directly, recompute sz with a second pass limited to z-terms:
            idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double Az = 0.0;
                for(int k1 = 0; k1 < n1; ++k1) {
                    double Bz = 0.0;
                    for(int k2 = 0; k2 < n2; ++k2, ++idx)
                        if(k2 > 0) Bz += coeff[idx] * (k2 * U2[k2 - 1]);
                    Az += Bz * T1[k1];
                }
                sz += Az * T0[k0];
            }

            f = s; fx = sx; fy = sy; fz = sz;
        }

        void EvalGradHess3D(Tensor1 x, out double f,
                            out double fx, out double fy, out double fz,
                            out double fxx, out double fxy, out double fxz,
                            out double fyy, out double fyz, out double fzz) {
            double x0 = x[0], x1 = x[1], x2 = x[2];
            int n0 = n[0], n1 = n[1], n2 = n[2];

            var T0 = new double[n0]; var U0 = new double[Math.Max(1, n0 - 1)];
            var dU0 = new double[Math.Max(1, n0 - 1)];
            var T1 = new double[n1]; var U1 = new double[Math.Max(1, n1 - 1)];
            var dU1 = new double[Math.Max(1, n1 - 1)];
            var T2 = new double[n2]; var U2 = new double[Math.Max(1, n2 - 1)];
            var dU2 = new double[Math.Max(1, n2 - 1)];
            Build_T_U_dU(x0, n0, T0, U0, dU0);
            Build_T_U_dU(x1, n1, T1, U1, dU1);
            Build_T_U_dU(x2, n2, T2, U2, dU2);

            double S = 0.0, Sx = 0.0, Sy = 0.0, Sz = 0.0;
            double Sxx = 0.0, Sxy = 0.0, Sxz = 0.0, Syy = 0.0, Syz = 0.0, Szz = 0.0;

            int idx = 0;
            for(int k0 = 0; k0 < n0; ++k0) {
                double A = 0.0, Ay = 0.0, Ayy = 0.0, Az = 0.0, Azz = 0.0, Ayz = 0.0;
                for(int k1 = 0; k1 < n1; ++k1) {
                    double B = 0.0, Bz = 0.0, Bzz = 0.0;
                    for(int k2 = 0; k2 < n2; ++k2, ++idx) {
                        double c = coeff[idx];
                        B += c * T2[k2];
                        if(k2 > 0) {
                            double kzU = k2 * U2[k2 - 1];
                            Bz += c * kzU;
                            Bzz += c * (k2 * dU2[k2 - 1]);
                        }
                    }
                    // along axis-1 derivatives
                    A += B * T1[k1];
                    Az += Bz * T1[k1];
                    Azz += Bzz * T1[k1];

                    if(k1 > 0) {
                        double kyU = k1 * U1[k1 - 1];
                        Ay += kyU * B;
                        Ayz += kyU * Bz;
                        Ayy += (k1 * dU1[k1 - 1]) * B;
                    }
                }

                S += A * T0[k0];
                Sy += Ay * T0[k0];
                Sz += Az * T0[k0];
                Syy += Ayy * T0[k0];
                Syz += Ayz * T0[k0];
                Szz += Azz * T0[k0];

                if(k0 > 0) {
                    double kxU = k0 * U0[k0 - 1];
                    Sx += kxU * A;
                    Sxy += kxU * Ay;
                    Sxz += kxU * Az;
                    Sxx += (k0 * dU0[k0 - 1]) * A;
                }
            }

            f = S;
            fx = Sx; fy = Sy; fz = Sz;
            fxx = Sxx; fxy = Sxy; fxz = Sxz;
            fyy = Syy; fyz = Syz; fzz = Szz;
        }

        // -------------------------
        // Chebyshev helpers
        // -------------------------

        // Chebyshev–Lobatto nodes j=0..n-1: cos(pi*j/(n-1))
        static double[] ChebLobatto(int n) {
            var x = new double[n];
            if(n == 1) { x[0] = 1.0; return x; }
            double denom = n - 1;
            for(int j = 0; j < n; ++j)
                x[j] = Math.Cos(Math.PI * j / denom);
            return x;
        }

        // Build T_k(x) for k=0..n-1 and U_k(x) for k=0..n-2
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void Build_T_and_U(double x, int n, double[] T, double[] U) {
            if(n == 1) { T[0] = 1.0; return; }
            T[0] = 1.0; T[1] = x;
            for(int k = 2; k < n; ++k)
                T[k] = 2.0 * x * T[k - 1] - T[k - 2];

            int m = n - 1;
            if(m <= 0) return;
            U[0] = 1.0; if(m == 1) return;
            U[1] = 2.0 * x;
            for(int k = 2; k < m; ++k)
                U[k] = 2.0 * x * U[k - 1] - U[k - 2];
        }

        // Also build dU_k/dx via stable three-term recurrence:
        // dU_0=0, dU_1=2, dU_{k+1} = 2*U_k + 2*x*dU_k - dU_{k-1}
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void Build_T_U_dU(double x, int n, double[] T, double[] U, double[] dU) {
            Build_T_and_U(x, n, T, U);
            int m = n - 1;
            if(m <= 0) return;
            dU[0] = 0.0;
            if(m == 1) return;
            dU[1] = 2.0;
            for(int k = 2; k < m; ++k)
                dU[k] = 2.0 * U[k - 1] + 2.0 * x * dU[k - 1] - dU[k - 2];
        }

        // -------------------------
        // Separable DCT-I (O(n^2))
        // -------------------------

        // In-place DCT-I along axis 'ax' of a flattened tensor with extents n[] and given strides.
        // Scaling: classic Clenshaw–Curtis interpolation:
        //   c_k = 2/(N-1) * sum_{j=0}^{N-1} s_j * f_j * cos(pi*k*j/(N-1)),  s_0=s_{N-1}=1/2 else 1
        // After this, 'coeff' holds Chebyshev T_k coefficients suitable for direct evaluation via T_k.
        static void DCT1_AlongAxis_InPlace(double[] a, int[] ext, int[] stride, int ax) {
            int N = ext[ax];
            if(N == 1) return;

            // Precompute cos table
            double denom = N - 1;
            var COS = new double[N * N];
            for(int k = 0; k < N; ++k) {
                for(int j = 0; j < N; ++j) {
                    COS[k * N + j] = Math.Cos(Math.PI * k * j / denom);
                }
            }

            // Iterate all lines along axis ax
            int lines = a.Length / N;
            int step = stride[ax];
            int block = N * step;

            // We iterate lexicographically over all indices except 'ax'
            // Using a mixed-radix counter
            int[] idx = new int[ext.Length];
            while(true) {
                // Base offset for this line
                int baseOff = 0;
                for(int d = 0; d < ext.Length; ++d)
                    baseOff += idx[d] * stride[d];

                // Gather f_j along this line
                // (indices along ax vary, others fixed)
                var fj = new double[N];
                for(int j = 0, off = baseOff; j < N; ++j, off += step)
                    fj[j] = a[off];

                // Apply Clenshaw–Curtis DCT-I scaling
                double s0 = 0.5 * fj[0];
                double sN = 0.5 * fj[N - 1];
                var outk = new double[N];
                for(int k = 0; k < N; ++k) {
                    double acc = s0 * COS[k * N + 0] + sN * COS[k * N + (N - 1)];
                    for(int j = 1; j < N - 1; ++j)
                        acc += fj[j] * COS[k * N + j];
                    acc *= (2.0 / (N - 1));
                    // No extra halving here; coefficients are ready for direct T_k evaluation with our recurrences.
                    outk[k] = acc;
                }

                // Scatter back along the line
                for(int k = 0, off = baseOff; k < N; ++k, off += step)
                    a[off] = outk[k];

                // Advance mixed-radix counter (skip axis 'ax')
                int d2;
                for(d2 = ext.Length - 1; d2 >= 0; --d2) {
                    if(d2 == ax) continue;
                    idx[d2]++;
                    if(idx[d2] < ext[d2]) break;
                    idx[d2] = 0;
                }
                if(d2 < 0) break; // done
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        void NormalizeChebyshevCoeffsInPlace() {
            int total = coeff.Length;
            if(M == 1) {
                // 1D: halve the k=0 term
                if(n[0] > 0) coeff[0] *= 0.5;
                return;
            }

            // General d-D (up to 3D here)
            for(int lin = 0; lin < total; ++lin) {
                int t = lin;
                int mZeros = 0;
                // decode indices in the same (last-axis-fastest) order as 'stride'
                for(int d = M - 1; d >= 0; --d) {
                    int kd = t % n[d];
                    if(kd == 0) mZeros++;
                    t /= n[d];
                }
                if(mZeros == 1) coeff[lin] *= 0.5;
                else if(mZeros == 2) coeff[lin] *= 0.25;
                else if(mZeros == 3) coeff[lin] *= 0.125; // 3D
                                                          // else mZeros==0 => no scaling
            }
        }

        // -------------------------
        // NodeSet assembly (tensor)
        // -------------------------

        NodeSet BuildTensorNodeSet(double[][] perAxisNodes) {
            // total nodes
            int K = 1;
            for(int d = 0; d < M; ++d) K *= n[d];

            var X = new NodeSet(Kref, K, M, false);

            // Fill tensor grid in row-major order: k = k0*n1*n2 + k1*n2 + k2 ...
            int[] idx = new int[M];
            for(int k = 0; k < K; ++k) {
                // decode multi-index
                int t = k;
                for(int d = M - 1; d >= 0; --d) {
                    idx[d] = t % n[d];
                    t /= n[d];
                }
                for(int d = 0; d < M; ++d) {
                    // Reference element assumed [-1,1]^d
                    X[k, d] = perAxisNodes[d][idx[d]];
                }
            }
            X.LockForever();
            return X;
        }

        static bool IsTensorProductRefElement(RefElement Kref) {
           if(Kref.GetType() == typeof(Line))
                return true;
            if(Kref.GetType() == typeof(Square))
                return true;
            if(Kref.GetType() == typeof(Cube))
                return true;

            return false;
        }
    }
}


