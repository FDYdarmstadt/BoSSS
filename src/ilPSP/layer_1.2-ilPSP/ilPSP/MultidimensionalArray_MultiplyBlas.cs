using System;
using ilPSP.Utils;

namespace ilPSP {
    public partial class MultidimensionalArray {
        /// <summary>
        /// Recognizes a matrix product (optionally batched over the first running
        /// index) from actual strides, not from symbolic index names. No packing.
        /// The caller excludes transformed indices and shared output/input storage.
        /// </summary>
        static unsafe bool TryMultiplyBlas(int runDim, int sumDim,
            double* t, double* a, double* b, int* lengths,
            int* ct, int* ca, int* cb, int* sums, int* sa, int* sb,
            double alpha, double beta) {
            if(sumDim != 1 || (runDim != 2 && runDim != 3))
                return false;

            int row = runDim - 2, col = row + 1;
            int m = lengths[row], n = lengths[col], k = sums[0];
            // Single-thread measurements: native-call overhead dominates tiny
            // products. Keep short DG coordinate reductions on the unrolled path.
            // Require at least a 16^3 product's work, with no extent below eight.
            if(m < 8 || n < 8 || k < 8 || (long)m * n * k < 4096 || alpha == 0.0)
                return false;
            if(ca[col] != 0 || cb[row] != 0 || ct[col] != 1 || ct[row] < n)
                return false;

            // Row-major C is column-major C^T = B^T A^T. Either input may
            // already be transposed; padded leading dimensions are allowed.
            int transB, ldb;
            if(sa[0] == 1 && ca[row] >= k) {
                transB = 'n'; ldb = ca[row];
            } else if(ca[row] == 1 && sa[0] >= m) {
                transB = 't'; ldb = sa[0];
            } else {
                return false;
            }
            int transA, lda;
            if(cb[col] == 1 && sb[0] >= n) {
                transA = 'n'; lda = sb[0];
            } else if(sb[0] == 1 && cb[col] >= k) {
                transA = 't'; lda = cb[col];
            } else {
                return false;
            }

            int batches = row == 0 ? 1 : lengths[0];
            if(row != 0 && (ct[0] < (long)(m - 1) * ct[row] + n || ca[0] < 0 || cb[0] < 0))
                return false;
            int ldc = ct[row];
            for(int batch = 0; batch < batches; batch++) {
                double nativeBeta = beta;
                if(beta == 0.0) {
                    // Unlike DGEMM, Multiply evaluates beta*C even for beta=0.
                    // Preserve NaN/Infinity propagation in the old destination.
                    for(int i = 0; i < m; i++)
                        for(int j = 0; j < n; j++)
                            t[i * ldc + j] *= beta;
                    nativeBeta = 1.0;
                }
                BLAS.m_seq_BLAS.DGEMM(ref transA, ref transB,
                    ref n, ref m, ref k, ref alpha, b, ref lda, a, ref ldb,
                    ref nativeBeta, t, ref ldc);
                if(row != 0) {
                    t += ct[0]; a += ca[0]; b += cb[0];
                }
            }
            return true;
        }
    }
}
