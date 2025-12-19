using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using IntersectingQuadrature.Tensor;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {



    /// <summary>
    /// Fast in-cell function using our own Chebyshev tensor polynomial.
    /// One-time: sample BoSSS at a tensor Chebyshev–Lobatto grid, compute coefficients by separable DCT-I.
    /// Runtime: evaluate by separable Clenshaw-style contractions, incl. gradient and Hessian.
    /// </summary>
    class GlobalCellFunctionOptimized : IScalarFunction {

        readonly RefElement Kref;
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


#if DEBUG
        readonly IScalarFunction reference;
#endif



        /// <summary></summary>
        public GlobalCellFunctionOptimized(RefElement Kref, IScalarFunction other, int degree) : this(Kref, Kref.SpatialDimension.ForLoop(d => degree + 1)) {
#if DEBUG
            reference = other;
#endif

            // 1) Build tensor NodeSet of all CL nodes & Sample BoSSS
            Tensor1[] xk = BuildTensors();
            MultidimensionalArray V = MultidimensionalArray.Create(1, xk.Length);
            for(int k = 0; k < xk.Length; k++) {
                V[0, k] = other.Evaluate(xk[k]);
            }

            // 2) compute internal coefficients
            ConstructorCommon(V);

#if DEBUG
            // check if we recover the nodal values
            {

                MultidimensionalArray val_this = MultidimensionalArray.Create(1, xk.Length);
                for(int k = 0; k < xk.Length; k++) {
                    

                    val_this[0, k] = Evaluate(xk[k]);

                }

                if(val_this.L2Dist(V) > 1.0e-10)
                    throw new ArithmeticException("mismatch between nodal values and interpolation");

            }
#endif
        }


        /// <summary></summary>
        /// <param name="nPerAxis">
        /// Polynomial grid sizes per axis (n = degree+1).
        /// If null, you must provide it elsewhere; setup cost is paid once.
        /// </param>
        /// <param name="jCell"></param>
        /// <param name="levelSet"></param>
        public GlobalCellFunctionOptimized(LevelSetTracker.LevelSetData levelSet, int jCell, int[] nPerAxis) : this(levelSet.GridDat.iGeomCells.GetRefElement(jCell), nPerAxis) {
#if DEBUG
            reference = new GlobalCellFunction(levelSet, jCell);
#endif

            // 1) Build tensor NodeSet of all CL nodes & Sample BoSSS
            NodeSet X = BuildTensorNodeSet();
            MultidimensionalArray V = levelSet.GetLevSetValues(X, jCell, 1);

            // 2) compute internal coefficients
            ConstructorCommon(V);

            // check if we recover the nodal values
#if DEBUG
            {

                MultidimensionalArray val_this = MultidimensionalArray.Create(1, X.NoOfNodes);
                for(int k = 0; k < X.NoOfNodes; k++) {
                    Tensor1 xk = Tensor1.Zeros(M);
                    for(int d = 0; d < M; d++) {
                        xk[d] = X[k, d];
                    }

                    val_this[0, k] = Evaluate(xk);

                }

                if(val_this.L2Dist(V) > 1.0e-10)
                    throw new ArithmeticException("mismatch between nodal values and interpolation");

                for(int k = 0; k < X.NoOfNodes; k++) {
                    Tensor1 xk = Tensor1.Zeros(M);
                    for(int d = 0; d < M; d++) {
                        xk[d] = X[k, d];
                    }

                    EvaluateAndGradient(xk);
                    EvaluateAndGradientAndHessian(xk);
                }


            }
#endif
        }


        GlobalCellFunctionOptimized(RefElement Kref, int[] nPerAxis) {
            this.Kref = Kref;
            this.M = Kref.SpatialDimension;

            if(nPerAxis == null || nPerAxis.Length != M)
                throw new ArgumentException("nPerAxis must be provided with length equal to spatial dimension M.");
            if(!IsTensorProductRefElement(Kref))
                throw new NotSupportedException("This fast path assumes a tensor-product reference element ([-1,1]^D).");

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
        }



        private void ConstructorCommon(MultidimensionalArray V) {
            // 1) Move samples into a dense tensor 'F' (flattened)
            int total = coeff.Length;
            double[] F = new double[total];
            for(int k = 0; k < total; ++k)
                F[k] = V[0, k];

            // 2) Separable DCT-I along each axis to produce Chebyshev coefficients
            //    In-place transforms on 'F' into 'coeff'.
            Array.Copy(F, coeff, total);
            for(int ax = 0; ax < M; ++ax) {
                DCT1_AlongAxis_InPlace(coeff, n, stride, ax);
                HalfEndpointsAlongAxisInPlace(coeff, n, stride, ax); // <-- add this
            }
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

#if DEBUG
            var valgradRef = reference.EvaluateAndGradient(x);
            if(Math.Abs(valgradRef.evaluation - val) > 1.0e-10) {
                throw new ArithmeticException("optimized version sux 1");
            }
            var gradDist = (valgradRef.gradient - grad);
            double abs_gradDist = 0;
            for(int i = 0; i < M; i++) {
                abs_gradDist += gradDist[i].Abs();
            }
            if(abs_gradDist > 1.0e-10) {

                double[] dxi = new double[M];
                double eps = GenericBlas.MachineEps.Sqrt();
                double[] vp = new double[M];
                double[] vm = new double[M];
                double[] vp1 = new double[M];
                double[] vm1 = new double[M];
                for(int i = 0; i < M; i++) {
                    var xp = (Tensor1) x.Clone();
                    xp[i] += eps;
                    vp[i] = reference.Evaluate(xp);
                    vp1[i] = this.Evaluate(xp);

                    var xm = (Tensor1) x.Clone();
                    xm[i] -= eps;
                    vm[i] = reference.Evaluate(xm);
                    vm1[i] = this.Evaluate(xm);

                    dxi[i] = (vp[i] - vm[i])/(2*eps);
                }
                

                throw new ArithmeticException("optimized version sux 2: " + abs_gradDist);
            }
#endif


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

#if DEBUG
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
#endif

            return (val, grad, hess);
        }

        // ---------------
        // 1D evaluation
        // ---------------

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        double Eval1D(Tensor1 X) {
            double x = X[0]; // assumed reference coord in [-1,1]
            int n0 = n[0];

            unsafe {
                // T_k 
                double* T = stackalloc double[n0];
                Build_T(x, n0, T);

                // f = sum_k c_k T_k
                double s = 0.0;
                int idx = 0;
                for(int k = 0; k < n0; ++k, idx += 1)
                    s += coeff[idx] * T[k];
                return s;
            }
        }

        void EvalGrad1D(Tensor1 X, out double f, out double fx) {
            double xx = X[0];
            int n0 = n[0];

            unsafe {
                double* T =  stackalloc double[n0];
                double* dT = stackalloc double[n0];
                Build_T_dT(xx, n0, T, dT);

                double s = 0.0, sx = 0.0;
                int idx = 0;
                for(int k = 0; k < n0; ++k, ++idx) {
                    double ck = coeff[idx];
                    s += ck * T[k];
                    sx += ck * dT[k];
                }
                f = s;
                fx = sx;
            }
        }

        void EvalGradHess1D(Tensor1 x, out double f, out double fx, out double fxx) {
            double xx = x[0];
            int n0 = n[0];

            

            unsafe {
                double* T =  stackalloc double[n0];
                double* dT = stackalloc double[n0];
                double* ddT = stackalloc double[n0];
                Build_T_dT_ddT(xx, n0, T, dT, ddT);

                double s = 0.0, sx = 0.0, sxx = 0.0;
                int idx = 0;
                for(int k = 0; k < n0; ++k, ++idx) {
                    double ck = coeff[idx];
                    s += ck * T[k];
                    sx += ck * dT[k];
                    sxx += ck * ddT[k];
                }
                f = s;
                fx = sx;
                fxx = sxx;
            }
        }

        // ---------------
        // 2D evaluation
        // ---------------

       

        double Eval2D(Tensor1 X) {
            double x0 = X[0], x1 = X[1];
            int n0 = n[0], n1 = n[1];

            unsafe {
                double* T0 = stackalloc double[n0];
                double* T1 = stackalloc double[n1];
                Build_T(x0, n0, T0);
                Build_T(x1, n1, T1);

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
        }

    

        void EvalGrad2D(Tensor1 X, out double f, out double fx, out double fy) {
            double x0 = X[0], x1 = X[1];
            int n0 = n[0], n1 = n[1];

            unsafe {
                double* T0 = stackalloc double[n0];
                double* dT0 = stackalloc double[n0];
                double* T1 = stackalloc double[n1];
                double* dT1 = stackalloc double[n0];
                Build_T_dT(x0, n0, T0, dT0);
                Build_T_dT(x1, n1, T1, dT1);

                double s = 0.0, sx = 0.0, sy = 0.0;
                int idx = 0;
                for(int k0 = 0; k0 < n0; ++k0) {
                    double a = 0.0, ax = 0.0, ay = 0.0;
                    for(int k1 = 0; k1 < n1; ++k1, ++idx) {
                        double ck = coeff[idx];
                        a += ck * T1[k1];
                        //ax += ck * T1[k1];
                        ay += ck * dT1[k1];
                    }
                    ax = a;
                    s += a * T0[k0];
                    sx += ax * dT0[k0];
                    sy += ay * T0[k0];
                }
                f = s; fx = sx; fy = sy;
            }
        }

        void EvalGradHess2D(Tensor1 x, out double f,
                    out double fx, out double fy,
                    out double fxx, out double fxy, out double fyy) {
            double x0 = x[0], x1 = x[1];
            int n0 = n[0], n1 = n[1];

           
            unsafe {
                double* T0 = stackalloc double[n0];
                double* dT0 = stackalloc double[n0];
                double* ddT0 = stackalloc double[n0];
                double* T1 = stackalloc double[n1];
                double* dT1 = stackalloc double[n0];
                double* ddT1 = stackalloc double[n0];
                Build_T_dT_ddT(x0, n0, T0, dT0, ddT0);
                Build_T_dT_ddT(x1, n1, T1, dT1, ddT1);

                double s = 0.0, sx = 0.0, sy = 0.0, sxx = 0, sxy = 0, syy = 0;
                int idx = 0;
                for(int k0 = 0; k0 < n0; ++k0) {
                    double a = 0.0, ax = 0.0, ay = 0.0, axx = 0, axy = 0, ayy = 0;
                    for(int k1 = 0; k1 < n1; ++k1, ++idx) {
                        double ck = coeff[idx];
                        a += ck * T1[k1];
                        //ax += ck * T1[k1];
                        ay += ck * dT1[k1];

                        //axx += ck*T1[k1];
                        axy += ck*dT1[k1];
                        ayy += ck*ddT1[k1];
                    }
                    ax = a;
                    axx = a;
                    s += a * T0[k0];
                    sx += ax * dT0[k0];
                    sy += ay * T0[k0];
                    sxx += axx * ddT0[k0];
                    sxy += axy * dT0[k0];
                    syy += ayy * T0[k0];

                }
                f = s; fx = sx; fy = sy; fxx = sxx; fxy = sxy; fyy = syy;
            }
        }

        // ---------------
        // 3D evaluation
        // ---------------

        double Eval3D(Tensor1 X) {
            double x0 = X[0], x1 = X[1], x2 = X[2];
            int n0 = n[0], n1 = n[1], n2 = n[2];

            unsafe {
                double* T0 = stackalloc double[n0];
                double* T1 = stackalloc double[n1];
                double* T2 = stackalloc double[n2];
                Build_T(x0, n0, T0);
                Build_T(x1, n1, T1);
                Build_T(x2, n2, T2);

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
        }

        void EvalGrad3D(Tensor1 x, out double f,
                         out double fx, out double fy, out double fz) {
            double x0 = x[0], x1 = x[1], x2 = x[2];
            int n0 = n[0], n1 = n[1], n2 = n[2];

            /*{
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
            }*/

            unsafe {
                double* T0 = stackalloc double[n0];
                double* dT0 = stackalloc double[n0];
                double* T1 = stackalloc double[n1];
                double* dT1 = stackalloc double[n0];
                double* T2 = stackalloc double[n1];
                double* dT2 = stackalloc double[n0];
                Build_T_dT(x0, n0, T0, dT0);
                Build_T_dT(x1, n1, T1, dT1);
                Build_T_dT(x2, n2, T2, dT2);

                // (((coeff • T2) • T1) • T0)
                double s = 0.0, sx = 0.0, sy = 0.0, sz = 0.0;
                int idx012 = 0;
                for(int k0 = 0; k0 < n0; ++k0) {
                    double a0 = 0.0, a0x = 0.0, a0y = 0.0, a0z = 0.0;
                    for(int k1 = 0; k1 < n1; ++k1) {
                        double a1 = 0.0, a1x = 0.0, a1y = 0.0, a1z = 0.0;
                        for(int k2 = 0; k2 < n2; ++k2, ++idx012) {
                            double ck = coeff[idx012];
                            a1 += ck * T2[k2];
                            //a1x += ck * T2[k2];
                            //a1y += ck * T2[k2];
                            a1z += ck * dT2[k2];
                        }
                        a1x = a1;
                        a1y = a1;
                        
                        a0 += a1 * T1[k1];
                        a0x += a1x *  T1[k1];
                        a0y += a1y * dT1[k1];
                        a0z += a1z *  T1[k1];
                    }
                    s += a0 * T0[k0];
                    sx += a0x * dT0[k0];
                    sy += a0y *  T0[k0];
                    sz += a0z *  T0[k0];

                }
                f = s; fx = sx; fy = sy; fz = sz;
            }
        }

        void EvalGradHess3D(Tensor1 x, out double f,
                            out double fx, out double fy, out double fz,
                            out double fxx, out double fxy, out double fxz,
                            out double fyy, out double fyz, out double fzz) {
            double x0 = x[0], x1 = x[1], x2 = x[2];
            int n0 = n[0], n1 = n[1], n2 = n[2];

            
            unsafe {
                double* T0 = stackalloc double[n0];
                double* dT0 = stackalloc double[n0];
                double* ddT0 = stackalloc double[n0];
                double* T1 = stackalloc double[n1];
                double* dT1 = stackalloc double[n0];
                double* ddT1 = stackalloc double[n0];
                double* T2 = stackalloc double[n1];
                double* dT2 = stackalloc double[n0];
                double* ddT2 = stackalloc double[n0];
                Build_T_dT_ddT(x0, n0, T0, dT0, ddT0);
                Build_T_dT_ddT(x1, n1, T1, dT1, ddT1);
                Build_T_dT_ddT(x2, n2, T2, dT2, ddT2);

                // (((coeff • T2) • T1) • T0)
                double s = 0.0, sx = 0.0, sy = 0.0, sz = 0.0;
                double sxx = 0.0, sxy = 0.0, sxz = 0.0, syy = 0.0, syz = 0.0, szz = 0.0;

                int idx012 = 0;
                for(int k0 = 0; k0 < n0; ++k0) {
                    double a0 = 0.0, a0x = 0.0, a0y = 0.0, a0z = 0.0;
                    double a0xx = 0.0, a0xy = 0.0, a0xz = 0.0, a0yy = 0.0, a0yz = 0.0, a0zz = 0.0;

                    for(int k1 = 0; k1 < n1; ++k1) {
                        double a1 = 0.0, a1x = 0.0, a1y = 0.0, a1z = 0.0;
                        double a1xx = 0.0, a1xy = 0.0, a1xz = 0.0, a1yy = 0.0, a1yz = 0.0, a1zz = 0.0;

                        for(int k2 = 0; k2 < n2; ++k2, ++idx012) {
                            double ck = coeff[idx012];
                            a1 += ck * T2[k2];
                            //a1x += ck * T2[k2];
                            //a1y += ck * T2[k2];
                            a1z += ck * dT2[k2];

                            //axx += ck*T2[k2];
                            //a1xy += ck*T2[k2];
                            //a1xz += ck*dT2[k2];
                            //a1yy += ck*T2[k2];
                            //a1yz += ck*dT2[k2];
                            a1zz += ck*ddT2[k2];
                        }
                        a1x = a1;
                        a1y = a1;
                        a1xx = a1;
                        a1xy = a1;
                        a1xz = a1z;
                        a1yy = a1;
                        a1yz = a1z;
                        
                        a0 += a1 * T1[k1];
                        a0x += a1x *  T1[k1];
                        a0y += a1y * dT1[k1];
                        a0z += a1z *  T1[k1];
                        a0xx += a1xx * T1[k1];
                        a0xy += a1xy * dT1[k1];
                        a0xz += a1xz * T1[k1];
                        a0yy += a1yy * ddT1[k1];
                        a0yz += a1yz * dT1[k1];
                        a0zz += a1zz * T1[k1];
                    }
                    s += a0 * T0[k0];
                    sx += a0x * dT0[k0];
                    sy += a0y *  T0[k0];
                    sz += a0z *  T0[k0];
                    sxx += a0xx * ddT0[k0];
                    sxy += a0xy * dT0[k0];
                    sxz += a0xz * dT0[k0];
                    syy += a0yy * T0[k0];
                    syz += a0yz * T0[k0];
                    szz += a0zz * T0[k0];

                }
                f = s; fx = sx; fy = sy; fz = sz;
                fxx = sxx; fxy = sxy; fxz = sxz;
                fyy = syy; fyz = syz; fzz = szz;
            }
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

        
        // Build T_k(x) for k=0..n-1
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static unsafe void Build_T(double x, int n, double* T) {
            if(n == 1) {
                T[0] = 1.0;
                return;
            }
            T[0] = 1.0;
            T[1] = x;
            for(int k = 2; k < n; ++k)
                T[k] = 2.0 * x * T[k - 1] - T[k - 2];
        }

        // Build T_k(x) and dT_k(x)/dx for k=0..n-1
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static unsafe void Build_T_dT(double x, int n, double* T, double* dT) {
            Build_T(x, n, T);
            if(n == 1) {
                dT[0] = 0.0;
                return;
            }
            dT[0] = 0.0;
            dT[1] = 1;
            for(int k = 2; k < n; ++k)
                dT[k] = 2*T[k-1] + 2*x*dT[k-1] - dT[k-2];
        }

        // Build T_k(x) and dT_k(x)/dx and ddT_k(x)/dxdx for k=0..n-1
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static unsafe void Build_T_dT_ddT(double x, int n, double* T, double* dT, double* ddT) {
            Build_T_dT(x, n, T, dT);
            if(n == 1) {
                ddT[0] = 0.0;
                return;
            }
            ddT[0] = 0.0;
            ddT[1] = 0.0;
            for(int k = 2; k < n; ++k)
                ddT[k] = 4*dT[k-1] + 2*x*ddT[k-1] - ddT[k-2];
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
        static void HalfEndpointsAlongAxisInPlace(double[] a, int[] ext, int[] stride, int ax) {
            int N = ext[ax];
            if(N <= 1) return;

            int step = stride[ax];

            // Mixed-radix counter over all dims except 'ax'
            int D = ext.Length;
            var idx = new int[D];

            while(true) {
                // base offset for this line (k_ax = 0)
                int baseOff = 0;
                for(int d = 0; d < D; ++d)
                    baseOff += idx[d] * stride[d];

                // scale k_ax = 0
                a[baseOff] *= 0.5;

                // scale k_ax = N-1
                a[baseOff + (N - 1) * step] *= 0.5;

                // advance idx (skip axis 'ax')
                int d2;
                for(d2 = D - 1; d2 >= 0; --d2) {
                    if(d2 == ax) continue;
                    if(++idx[d2] < ext[d2]) break;
                    idx[d2] = 0;
                }
                if(d2 < 0) break;
            }
        }



        // -------------------------
        // NodeSet assembly (tensor)
        // -------------------------

        NodeSet BuildTensorNodeSet() {
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
                    X[k, d] = this.nodesCL[d][idx[d]];
                }
            }
            X.LockForever();
            return X;
        }

        Tensor1[] BuildTensors() {
            // total nodes
            int K = 1;
            for(int d = 0; d < M; ++d) 
                K *= n[d];

            var ret = new Tensor1[K];

            // Fill tensor grid in row-major order: k = k0*n1*n2 + k1*n2 + k2 ...
            int[] idx = new int[M];
            for(int k = 0; k < K; ++k) {
                // decode multi-index
                int t = k;
                for(int d = M - 1; d >= 0; --d) {
                    idx[d] = t % n[d];
                    t /= n[d];
                }
                ret[k] = Tensor1.Zeros(M);
                for(int d = 0; d < M; ++d) {
                    // Reference element assumed [-1,1]^d
                    ret[k][d] = nodesCL[d][idx[d]];
                }
            }

            return ret;
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


