using System;
using ilPSP.Utils;

namespace ilPSP {
    public partial class MultidimensionalArray {

        // ==========================================================================================
        // Multiply paths via optional native BLAS extensions (PLAN-5 Phase I)
        // ------------------------------------------------------------------------------------------
        // Candidates are only enabled in the selection (TryMultiplyBlas) where they were measured to
        // win, see .agent/kummer/tensorproductoptimization/RESULTS-5.md; the switches allow benchmarks
        // (GettF replay) to evaluate them on census shapes.
        // ==========================================================================================

        /// <summary>
        /// Switch for <see cref="TryBlisGett"/> in the selection
        /// </summary>
        internal static bool s_BlisGett = false;

        /// <summary>
        /// minimal GEMM size (m*n*k per BLIS call) for <see cref="TryBlisGett"/>
        /// </summary>
        internal static long s_BlisGettMinMnk = 64 * 64 * 64;

        /// <summary>
        /// Switch for one MKL strided-batch call instead of the loop over batch entries in <see cref="TryDirectGemm"/>;
        /// OFF: no measurable gain over the loop of DGEMM calls (AMD EPYC 9654, RESULTS-5.md).
        /// </summary>
        internal static bool s_MklBatchDirect = false;

        /// <summary>
        /// Switch for <see cref="TryJitGemm"/>: MKL JIT kernels for small GEMMs with BLAS-compatible layout.
        /// ON: on AMD EPYC 9654 1.5x-9x faster than the blocked SIMD path and 1.1x-5x faster than DGEMM for n = 4...35 (RESULTS-5.md);
        /// results differ from the generated kernels at round-off level (like all BLAS paths).
        /// </summary>
        internal static bool s_MklJitGemm = true;

        /// <summary>
        /// maximal m*n*k for <see cref="TryJitGemm"/> (above, DGEMM is as fast)
        /// </summary>
        internal const long JIT_MAX_MNK = 64 * 64 * 64;

        /// <summary>
        /// Small GEMM (batch) with the layout conditions of <see cref="TryDirectGemm"/>, but without its size limits,
        /// via a cached MKL JIT kernel (shape, leading dimensions, alpha, beta). False, without touching T, if not eligible,
        /// JIT is unavailable, MKL declined the shape or the kernel cache is full.
        /// </summary>
        static unsafe bool TryJitGemm(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta) {
            if(runDim != 2 && runDim != 3) return false;
            int row = runDim - 2;
            if(!JitGemmKernel(lengths + row, ct + row, ca + row, cb + row, k, sa, sb, alpha, beta, out var kern))
                return false;
            int m = lengths[row], n = lengths[row + 1], ldc = ct[row];
            int batches = row == 0 ? 1 : lengths[0];
            if(beta == 0.0) {
                // all batch entries first, see JitGemmScale
                double* tq = t;
                for(int q = 0; q < batches; q++, tq += (row != 0 ? ct[0] : 0))
                    JitGemmScale(tq, m, n, ldc, beta);
            }
            for(int q = 0; q < batches; q++) {
                // column-major view: T^T (n x m) = op(B^T) * op(A^T)
                kern.Invoke(b, a, t);
                if(row != 0) { t += ct[0]; a += ca[0]; b += cb[0]; }
            }
            return true;
        }

        /// <summary>
        /// JIT kernel for one row-major GEMM T[i,j] = alpha*sum_l A[i,l] B[l,j] + beta*T[i,j] (two running loops i, j, one summation loop),
        /// if the layout is BLAS-compatible (as in <see cref="TryDirectGemm"/>, without its size limits)
        /// </summary>
        static unsafe bool JitGemmKernel(int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta, out MKLJit.Kernel kern) {
            kern = default;
            int m = lengths[0], n = lengths[1];
            if(m < 1 || n < 1 || k < 1 || (long)m * n * k > JIT_MAX_MNK) return false;
            if(ca[1] != 0 || cb[0] != 0 || ct[1] != 1 || ct[0] < n) return false;
            int transB, ldb;
            if(sa[0] == 1 && ca[0] >= k) { transB = MKLBatch.NoTrans; ldb = ca[0]; }
            else if(ca[0] == 1 && sa[0] >= m) { transB = MKLBatch.Trans; ldb = sa[0]; }
            else return false;
            int transA, lda;
            if(cb[1] == 1 && sb[0] >= n) { transA = MKLBatch.NoTrans; lda = sb[0]; }
            else if(sb[0] == 1 && cb[1] >= k) { transA = MKLBatch.Trans; lda = cb[1]; }
            else return false;
            // beta = 0: T is scaled first (keeps NaN/Inf propagation of the generated kernels), kernel with beta = 1
            return MKLJit.TryGetKernel(MKLBatch.ColMajor, transA, transB, n, m, k, alpha, lda, ldb, beta == 0.0 ? 1.0 : beta, ct[0], out kern);
        }

        /// <summary>
        /// T *= beta for one m x n block (for beta = 0: keeps NaN/Inf propagation of the generated kernels; the kernel then uses beta = 1).
        /// Must be done for all blocks before the first kernel call: scaling directly before the kernel of the same block stalls the
        /// kernel's vector loads (store-to-load forwarding; 13x slower on AMD EPYC 9654, RESULTS-5.md).
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveOptimization)] // tier-0/OSR code of this loop is ~10x slower (RESULTS-5.md)
        static unsafe void JitGemmScale(double* t, int m, int n, int ldc, double beta) {
            int W = System.Numerics.Vector<double>.Count;
            var vb = new System.Numerics.Vector<double>(beta);
            for(int i = 0; i < m; i++) {
                double* ti = t + (long)i * ldc;
                int j = 0;
                for(; j + W <= n; j += W)
                    System.Runtime.CompilerServices.Unsafe.Write(ti + j, System.Runtime.CompilerServices.Unsafe.Read<System.Numerics.Vector<double>>(ti + j) * vb);
                for(; j < n; j++)
                    ti[j] *= beta;
            }
        }


        /// <summary>
        /// General contraction as a sequence of BLIS GEMMs with arbitrary row/column strides, directly on the
        /// operand storage (no packing): the running loops are classified as M (varies in A only), N (B only) and
        /// batch (both); the longest M loop, the longest N loop and the longest summation loop form the GEMM
        /// dimensions, all other loops (including a second summation loop) are executed outside.
        /// The caller excludes aliasing; T entries must be pairwise distinct.
        /// </summary>
        static unsafe bool TryBlisGett(int runDim, int sumDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int* sums, int* sa, int* sb, double alpha, double beta, long minMnk) {
            const int MAXD = 8;
            if(runDim < 2 || runDim > MAXD || sumDim < 1 || sumDim > 2) return false;
            int iM = -1, iN = -1;
            for(int d = 0; d < runDim; d++) {
                if(lengths[d] <= 0) return false;
                if(lengths[d] == 1) continue;
                if(ca[d] != 0 && cb[d] == 0) { if(iM < 0 || lengths[d] > lengths[iM]) iM = d; }
                else if(ca[d] == 0 && cb[d] != 0) { if(iN < 0 || lengths[d] > lengths[iN]) iN = d; }
                else if(ca[d] == 0 && cb[d] == 0) return false; // T-only (broadcast) loop
            }
            if(iM < 0 || iN < 0) return false;
            int kS = 0;
            if(sumDim > 1 && sums[1] > sums[0]) kS = 1;
            long m = lengths[iM], n = lengths[iN], k = sums[kS];
            if(m < 4 || n < 4 || k < 4 || m * n * k < minMnk) return false;
            if(ct[iM] == 0 || ct[iN] == 0 || sa[kS] == 0 || sb[kS] == 0) return false;
            if(!BLIS.Available) return false;
            if(!GettTNonOverlapping(runDim, lengths, ct)) return false;

            // outer loops: running loops other than iM, iN (with length > 1), and the other summation loop
            int* oDim = stackalloc int[MAXD + 1];
            int* oLen = stackalloc int[MAXD + 1];
            int* oT = stackalloc int[MAXD + 1];
            int* oA = stackalloc int[MAXD + 1];
            int* oB = stackalloc int[MAXD + 1];
            int nO = 0;
            for(int d = 0; d < runDim; d++) {
                if(d == iM || d == iN || lengths[d] == 1) continue;
                oLen[nO] = lengths[d]; oT[nO] = ct[d]; oA[nO] = ca[d]; oB[nO] = cb[d]; nO++;
            }
            if(sumDim > 1) {
                int s = 1 - kS;
                if(sums[s] <= 0) return false;
                oLen[nO] = sums[s]; oT[nO] = 0; oA[nO] = sa[s]; oB[nO] = sb[s]; nO++;
            }

            // beta first (keeps NaN/Inf propagation of beta = 0 as in the generated kernels), then accumulate with beta = 1
            if(beta != 1.0)
                ScaleStrided(runDim, t, lengths, ct, beta);

            int* idx = stackalloc int[MAXD + 1];
            for(int j = 0; j < nO; j++) idx[j] = 0;
            while(true) {
                long ot = 0, oa = 0, ob = 0;
                for(int j = 0; j < nO; j++) { ot += (long)idx[j] * oT[j]; oa += (long)idx[j] * oA[j]; ob += (long)idx[j] * oB[j]; }
                BLIS.Dgemm(false, false, m, n, k, alpha,
                    a + oa, ca[iM], sa[kS],
                    b + ob, sb[kS], cb[iN],
                    1.0, t + ot, ct[iM], ct[iN]);
                int jj = nO - 1;
                for(; jj >= 0; jj--) { if(++idx[jj] < oLen[jj]) break; idx[jj] = 0; }
                if(jj < 0) break;
            }
            return true;
        }

        /// <summary> T *= beta over all running loops </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveOptimization)] // tier-0/OSR code of this loop is ~10x slower (RESULTS-5.md)
        static unsafe void ScaleStrided(int runDim, double* t, int* lengths, int* ct, double beta) {
            // innermost loop: the one with the smallest |stride|
            int inner = runDim - 1;
            for(int d = 0; d < runDim; d++)
                if(Math.Abs(ct[d]) < Math.Abs(ct[inner])) inner = d;
            int L = lengths[inner], c = ct[inner];
            int* idx = stackalloc int[runDim];
            for(int d = 0; d < runDim; d++) idx[d] = 0;
            while(true) {
                long o = 0;
                for(int d = 0; d < runDim; d++) if(d != inner) o += (long)idx[d] * ct[d];
                double* tp = t + o;
                for(int i = 0; i < L; i++) tp[i * c] *= beta;
                int jj = runDim - 1;
                for(; jj >= 0; jj--) {
                    if(jj == inner) continue;
                    if(++idx[jj] < lengths[jj]) break;
                    idx[jj] = 0;
                }
                if(jj < 0) break;
            }
        }

        /// <summary>
        /// Batched variant of the loop in <see cref="TryDirectGemm"/>: one MKL strided-batch call (column-major view, as there).
        /// </summary>
        static unsafe bool DirectGemmMklBatch(int batches, int transA, int transB, int n, int m, int k, double alpha,
            double* b, int lda, int strideB, double* a, int ldb, int strideA, double beta, double* t, int ldc, int strideT,
            int rowLen, int colLen) {
            if(batches < 2 || !MKLBatch.Available) return false;
            // MKL_INT is 32 bit: the last batch entry must be addressable
            if((long)strideT * (batches - 1) > int.MaxValue || (long)strideA * (batches - 1) > int.MaxValue || (long)strideB * (batches - 1) > int.MaxValue)
                return false;
            // batch entries are independent in MKL (no order): the T blocks must not overlap
            if(strideT <= 0 || strideT < (long)(rowLen - 1) * ldc + colLen || strideA < 0 || strideB < 0) return false;
            double nativeBeta = beta;
            if(beta == 0.0) {
                for(int q = 0; q < batches; q++) {
                    double* tq = t + (long)q * strideT;
                    for(int i = 0; i < rowLen; i++) for(int j = 0; j < colLen; j++) tq[i * ldc + j] *= beta;
                }
                nativeBeta = 1.0;
            }
            MKLBatch.DgemmBatchStrided(MKLBatch.ColMajor,
                transA == 'n' ? MKLBatch.NoTrans : MKLBatch.Trans, transB == 'n' ? MKLBatch.NoTrans : MKLBatch.Trans,
                n, m, k, alpha, b, lda, strideB, a, ldb, strideA, nativeBeta, t, ldc, strideT, batches);
            return true;
        }
    }
}
