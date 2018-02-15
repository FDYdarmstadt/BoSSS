/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ilPSP.Utils;

namespace ilPSP {

    partial class MultidimensionalArray {

        /// <summary>
        /// adds <paramref name="other"/>*<paramref name="scl"/> to this array,
        /// overwriting this array; this array must be 1-dimensional and the
        /// length of <paramref name="other"/> must be equal to the length of
        /// the 0-th dimension of this array;
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scl">
        /// scaling of <paramref name="other"/> during accumulation
        /// </param>
        public void AccVector<T>(double scl, T other) where T : IList<double> {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");
            if (this.Dimension != 1)
                throw new NotSupportedException("This array must be 1-dimensional to support the called method");
            if (other.Count != m_StorageLayout.m_Length0)
                throw new ArgumentException("mismatch in length of 0-th dimension.", "other");

            int I = m_StorageLayout.m_Length0;
            int ind;
            for (int i = 0; i < I; i++) {
                ind = Index(i);
                m_Storage[ind] += other[i] * scl;
            }
        }

        /// <summary>
        /// adds <paramref name="other"/>*<paramref name="scl"/> to this array,
        /// overwriting this array; this array must be 2-dimensional and the
        /// length of each dimension must match this array;
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scl">
        /// scaling of <paramref name="other"/> during accumulation
        /// </param>
        public void Acc2DArray(double scl, double[,] other) {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");
            if (this.Dimension != 2)
                throw new NotSupportedException("This array must be 2-dimensional to support the called method");
            if (other.GetLength(0) != m_StorageLayout.m_Length0)
                throw new ArgumentException("mismatch in length of 0-th dimension.", "other");
            if (other.GetLength(1) != m_StorageLayout.m_Length1)
                throw new ArgumentException("mismatch in length of 1-st dimension.", "other");

            int I = m_StorageLayout.m_Length0;
            int J = m_StorageLayout.m_Length1;
            int ind;
            for (int i = 0; i < I; i++) {
                for (int j = 0; j < J; j++) {
                    ind = Index(i, j);

                    m_Storage[ind] += other[i, j] * scl;
                }
            }
        }

        /// <summary>
        /// see <see cref="AccMatrix"/>
        /// </summary>
        void IMatrix.Acc(double scl, IMatrix other) {
             this.AccMatrix(scl, other);
        }

        /// <summary>
        /// adds <paramref name="other"/>*<paramref name="scl"/> to this array,
        /// overwriting this array; this array must be 2-dimensional and the
        /// length of each dimension must match this array;
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scl">
        /// scaling of <paramref name="other"/> during accumulation
        /// </param>
        public void AccMatrix(double scl, IMatrix other) {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");
            if (this.Dimension != 2)
                throw new NotSupportedException("This array must be 2-dimensional to support the called method");
            if (other.NoOfRows != m_StorageLayout.m_Length0)
                throw new ArgumentException("mismatch in length of 0-th dimension.", "other");
            if (other.NoOfCols != m_StorageLayout.m_Length1)
                throw new ArgumentException("mismatch in length of 1-st dimension.", "other");

            int I = m_StorageLayout.m_Length0;
            int J = m_StorageLayout.m_Length1;
            int ind;
            for (int i = 0; i < I; i++) {
                for (int j = 0; j < J; j++) {
                    ind = Index(i, j);

                    m_Storage[ind] += other[i, j] * scl;
                }
            }
        }

        /// <summary>
        /// adds <paramref name="other"/>*<paramref name="scl"/> to this array
        /// overwriting this array; this array must have the same dimension as
        /// <paramref name="other"/> and the length of each dimension must match.
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scl">
        /// scaling of <paramref name="other"/> during accumulation
        /// </param>
        public void Acc(double scl, MultidimensionalArray other) {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");
            if (this.Dimension != other.Dimension)
                throw new ArgumentException("mismatch in number of dimensions.", "other");

            if (other.m_Dimension != this.m_Dimension)
                throw new ArgumentException("Number of dimensions must be equal.");
            for (int i = this.m_Dimension - 1; i >= 0; i--) {
                if (other.GetLength(i) != this.GetLength(i))
                    throw new ArgumentException("Mismatch in dimension " + i);
            }

            // Optimized accumulation for continuous pieces of memory
            if (this.IsContinious && other.IsContinious) {
                int[] index = new int[Dimension];
                int thisOffset = this.Index(index);
                int otherOffset = other.Index(index);

                unsafe {
                    fixed (double* pThis = &this.Storage[thisOffset], pOther = &other.Storage[otherOffset]) {
                        //for (int i = 0; i < this.Length; i++) {
                        //    *(pThis + i) += scl * *(pOther + i);
                        //}
                        BLAS.daxpy(this.Length, scl, pOther, 1, pThis, 1);
                    }
                }

                return;
            }

            // Standard versions
            switch(this.Dimension) {
                case 2: {
                    int L0 = this.GetLength(0);
                    int L1 = this.GetLength(1);

                    for(int i0 = 0; i0 < L0; i0++)
                        for(int i1 = 0; i1 < L1; i1++) {
                            int ind_this = this.Index(i0, i1);
                            int ind_othr = other.Index(i0, i1);

                            this.m_Storage[ind_this] += other.m_Storage[ind_othr] * scl;
                        }
                    return;
                }

                default: {
                    double[] other_stor = other.m_Storage;

                    ApplyAll(delegate(int[] idx, ref double entry) {
                        int ind_other = other.Index(idx);
                        entry += scl * other_stor[ind_other];
                    });
                    return;
                }
            }
        }


        /// <summary>
        /// Adds the constant <paramref name="a"/> to all entries.
        /// </summary>
        public void AccConstant(double a) {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");
            

            // Standard versions
            switch(this.Dimension) {
                case 2: {
                    int L0 = this.GetLength(0);
                    int L1 = this.GetLength(1);

                    for(int i0 = 0; i0 < L0; i0++)
                        for(int i1 = 0; i1 < L1; i1++) {
                            int ind_this = this.Index(i0, i1);
                            
                            this.m_Storage[ind_this] += a;
                        }
                    return;
                }

                default: {
                    
                    ApplyAll(delegate(int[] idx, ref double entry) {
                        
                        entry += a;
                    });
                    return;
                }
            }
        }

        /// <summary>
        /// accumulates <paramref name="x"/>*<paramref name="alpha"/>
        /// to a sub-section of this array.
        /// </summary>
        /// <param name="x">values to accumulate</param>
        /// <param name="alpha">
        /// scaling for the values <paramref name="x"/>
        /// </param>
        /// <param name="Istart">start indices (including)</param>
        /// <param name="Iend">end indices (including)</param>
        public void AccSubArray(double alpha, MultidimensionalArray x, int[] Istart, int[] Iend) {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);
            sub.Acc(alpha, x);
        }

        /// <summary>
        /// accumulates <paramref name="x"/>*<paramref name="alpha"/>
        /// to a sub-section of this array.
        /// </summary>
        /// <param name="x">values to accumulate</param>
        /// <param name="alpha">
        /// scaling for the values <paramref name="x"/>
        /// </param>
        /// <param name="SubArrayIdx"></param>
        public void AccSubArray(double alpha, MultidimensionalArray x, params int[] SubArrayIdx) {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);
            sub.Acc(alpha, x);
        }

        /// <summary>
        /// used by <see cref="ApplyAll(ApplyAllOperation)"/>.
        /// </summary>
        /// <param name="entry">
        /// the value of some entry of this array (input/output)
        /// </param>
        /// <param name="index">
        /// index of <paramref name="entry"/>
        /// </param>
        public delegate void ApplyAllOperation(int[] index, ref double entry);

        /// <summary>
        /// applies the operation <paramref name="op"/> to all entries of this array
        /// </summary>
        public void ApplyAll(ApplyAllOperation op) {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");

            // special treatment for empty arrays
            for (int i = this.Dimension - 1; i >= 0; i--) {
                if (this.GetLength(i) <= 0)
                    return;
            }

            // do all...
            unsafe {
                int[] _index = new int[m_Dimension];
                int index_lin = this.Index(_index);
                fixed (int* index = _index) {
                    fixed (double* pStorage = this.m_Storage) {
                        ApplyAllRec(0, index, _index, index_lin, pStorage, op);
                    }
                }
            }
        }

        /// <summary>
        /// applies the operation <paramref name="op"/> to all entries of this array
        /// </summary>
        public void ApplyAll(Func<double, double> op) {
            if(m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");

            if (this.IsContinious) {
                // optimized branch 
                int L = this.Length;
                int i0 = this.m_Offset;
                L += i0;
                double[] st = this.m_Storage;
                for (int i = i0; i < L; i++) {
                    st[i] = op(st[i]);
                }
            } else {
                ApplyAll(delegate(int[] i, ref double d) {
                    d = op(d);
                });
            }
        }


        /// <summary>
        /// applies the action <paramref name="op"/> to all entries of this array
        /// </summary>
        public void ApplyAll(Action<int[], double> op) {
            bool oldState = this.m_LockedForever;
            this.m_LockedForever = false;
            try {
                ApplyAll(delegate(int[] i, ref double d) {
                    op(i, d);
                });
            } catch(Exception e) {
                throw e;
            } finally {
                this.m_LockedForever = oldState;
            }
        }


        /// <summary>
        /// applies the action <paramref name="op"/> to all entries of this array
        /// </summary>
        public void ApplyAll(Action<double> op) {
            if (this.IsContinious) {
                // optimized branch
                int L = this.Length;
                int i0 = this.m_Offset;
                L += i0;
                double[] st = this.m_Storage;
                for (int i = i0; i < L; i++) {
                    op(st[i]);
                }
            } else {
                ApplyAll(delegate(int[] i, ref double d) {
                    op(d);
                });
            }
        }

        /// <summary>
        /// applies the operation <paramref name="op"/> to all entries of this
        /// array
        /// </summary>
        unsafe private void ApplyAllRec(int dim, int* index, int[] _index, int index_lin, double* pStorage, ApplyAllOperation op) {
            Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
            int L = this.GetLength(dim);
            int C = this.GetCycle(dim);
            if (dim == m_Dimension - 1) {
                // end of recursion

                for (index[dim] = 0; index[dim] < L; index[dim]++) {
                    op(_index, ref pStorage[index_lin]);
                    index_lin += C;
                }
            } else {
                for (index[dim] = 0; index[dim] < L; index[dim]++) {
                    ApplyAllRec(dim + 1, index, _index, index_lin, pStorage, op);
                    index_lin += C;
                }
            }
        }

        /// <summary>
        /// multiplies all entries of this array by <paramref name="alpha"/>
        /// </summary>
        /// <param name="alpha"></param>
        public void Scale(double alpha) {
            if(this.IsContinious) {
                unsafe {
                    fixed(double* pStorage = this.m_Storage) {
                        BLAS.dscal(this.Length, alpha, pStorage + this.m_Offset, 1);
                    }
                }
            } else {
                ApplyAll(d => d * alpha);
            }
        }

        /// <summary>
        /// computes the inner product of this array and another array <paramref name="other"/>.
        /// </summary>
        /// <param name="other">
        /// must be of the same size a this array
        /// </param>
        /// <returns></returns>
        public double InnerProduct(MultidimensionalArray other) {
            if(this.Dimension != other.Dimension)
                throw new ArgumentException("mismatch in number of dimensions");
            int D = this.Dimension;
            for(int k = 0; k < D; k++) {
                if(this.GetLength(k) != other.GetLength(k))
                    throw new ArgumentException("mismatch in dimension " + k);
            }

            if(this.IsContinious && other.IsContinious) {
                unsafe {
                    fixed(double* pStorageThis = this.m_Storage, pStorageOther = other.m_Storage) {
                        return BLAS.ddot(this.Length, pStorageThis, 1, pStorageOther, 1);
                    }
                }

            } else {
                double acc = 0;
                this.ApplyAll(delegate(int[] index, ref double entry) {
                    acc += entry * other[index];
                });
                return acc;
            }
        }


        /// <summary>
        /// sets all entries in this array to <paramref name="x"/>
        /// </summary>
        /// <param name="x"></param>
        public void SetAll(double x) {
            ApplyAll(d => x);
        }

        /// <summary>
        /// generalized multiplication, triple product.
        /// </summary>
        public void Multiply(double scale, MultidimensionalArray A, MultidimensionalArray B, MultidimensionalArray C, double thisscale, string Tindex, string Aindex, string Bindex, string Cindex) {
            // check arguments
            // ===============
            if (Tindex.Length != this.Dimension)
                throw new ArgumentException();
            if (Aindex.Length != A.Dimension)
                throw new ArgumentException();
            if (Bindex.Length != B.Dimension)
                throw new ArgumentException();
            if (Cindex.Length != C.Dimension)
                throw new ArgumentException();

            // determine indices for intermediate result
            // =========================================

            var _Aindex = Aindex.ToCharArray();
            var _Bindex = Bindex.ToCharArray();
            var _Cindex = Cindex.ToCharArray();
            var _Tindex = Tindex.ToCharArray();

            List<char> _Qindex = new List<char>();
            List<int> _Qlength = new List<int>();


            for (int i = 0; i < _Aindex.Length; i++) {
                char g = _Aindex[i];
                int L = A.GetLength(i);

                if (Array.IndexOf<char>(_Tindex, g) >= 0 || Array.IndexOf<char>(_Cindex, g) >= 0) {
                    _Qindex.Add(g);
                    _Qlength.Add(L);
                }
            }

            for (int i = 0; i < _Bindex.Length; i++) {
                char g = _Bindex[i];
                int L = B.GetLength(i);

                if (!_Qindex.Contains(g)) {
                    if (Array.IndexOf<char>(_Tindex, g) >= 0 || Array.IndexOf<char>(_Cindex, g) >= 0) {
                        _Qindex.Add(g);
                        _Qlength.Add(L);
                    }
                }
            }

            // execute multiplication
            // ======================

            // Q = A*B
            //var Q = MultidimensionalArray.Create(_Qlength.ToArray());
            int iBuf;
            var Q = TempBuffer.GetTempMultidimensionalarray(out iBuf, _Qlength.ToArray());
            string Qindex = new string(_Qindex.ToArray());
            Q.Multiply(1.0, A, B, 0.0, Qindex, Aindex, Bindex);

            // this = scale*Q*C + this*thisscale
            this.Multiply(scale, Q, C, thisscale, Tindex, Qindex, Cindex);
            TempBuffer.FreeTempBuffer(iBuf);
        }
        

        /// <summary>
        /// Caches the infromation extracted form the string arguments of 
        /// the tensor-multiplication routine <see cref="Multiply(double,MultidimensionalArray,MultidimensionalArray, double, string, string, string)"/>.
        /// </summary>
        public struct MultiplyProgram {
            internal const int MAX_SUM_LOOPS = 4;

            internal int DT;
            internal int DA;
            internal int DB;

            internal int NoOfSumCycles;



            internal void CheckArgs(MultidimensionalArray T, MultidimensionalArray A, MultidimensionalArray B) {
                if(T.Dimension != DT)
                    throw new ArgumentException("Wrong dimension of result array.");
                if(A.Dimension != DA)
                    throw new ArgumentException("Wrong dimension of array A.");
                if(B.Dimension != DB)
                    throw new ArgumentException("Wrong dimension of array B.");

              
                for(int i = 0; i < RunACount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkA = (int)((RunRankAStore & mask) >> shift);
                    int sumA = (int)((RunLoopAStore & mask) >> shift);
                    Debug.Assert(sumA >= 0 && sumA < DT);

                    if((T.GetLength(sumA) != A.GetLength(rnkA)) && (iTrafoIdx != sumA))
                        throw new ArgumentException(string.Format("Wrong length of {0}-th rank of array A.", rnkA));
                } 
                
                for(int i = 0; i < RunBCount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkB = (int)((RunRankBStore & mask) >> shift);
                    int sumB = (int)((RunLoopBStore & mask) >> shift);
                    Debug.Assert(sumB >= 0 && sumB < DT);

                    if((T.GetLength(sumB) != B.GetLength(rnkB)) && (iTrafoIdx != sumB))
                        throw new ArgumentException(string.Format("Wrong length of {0}-th rank of array B.", rnkB));
                }
                
                int[] lenSum = new int[NoOfSumCycles];

                for(int i = 0; i < SumACount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkA = (int)((SumRankAStore & mask) >> shift);
                    int sumA = (int)((SumLoopAStore & mask) >> shift);
                    Debug.Assert(sumA >= 0 && sumA < NoOfSumCycles);
                    lenSum[sumA] = A.GetLength(rnkA);
                }
                
                for(int i = 0; i < SumBCount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkB = (int)((SumRankBStore & mask) >> shift);
                    int sumB = (int)((SumLoopBStore & mask) >> shift);
                    Debug.Assert(sumB >= 0 && sumB < NoOfSumCycles);
                    if(B.GetLength(rnkB) != lenSum[sumB])
                        throw new ArgumentException(string.Format("Wrong length of {0}-th rank of array B.", rnkB));
                }
            }

            unsafe internal void GetCyclesAndRanks(int* cycRunA, int* cycRunB, int* lenSum, int* cycSumA, int* cycSumB, MultidimensionalArray A, MultidimensionalArray B) {

                
#if DEBUG
                for(int i = 0; i < NoOfSumCycles; i++) {
                    Debug.Assert(cycSumA[i] == 0);
                    Debug.Assert(cycSumB[i] == 0);
                }

                for(int i = 0; i < DT; i++) {
                    Debug.Assert(cycRunA[i] == 0);
                    Debug.Assert(cycRunB[i] == 0);
                }
#endif 
                for(int i = 0; i < RunACount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkA = (int)((RunRankAStore & mask) >> shift);
                    int runA = (int)((RunLoopAStore & mask) >> shift);
                    Debug.Assert(runA >= 0 && runA < DT);
                    cycRunA[runA] += A.GetCycle(rnkA);
                }
                for(int i = 0; i < RunBCount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkB = (int)((RunRankBStore & mask) >> shift);
                    int runB = (int)((RunLoopBStore & mask) >> shift);
                    Debug.Assert(runB >= 0 && runB < DT);
                    cycRunB[runB] += B.GetCycle(rnkB);
                }
                
                
                for(int i = 0; i < SumACount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkA = (int)((SumRankAStore & mask) >> shift);
                    int sumA = (int)((SumLoopAStore & mask) >> shift);
                    Debug.Assert(sumA >= 0 && sumA < NoOfSumCycles);
                    cycSumA[sumA] += A.GetCycle(rnkA);
                    lenSum[sumA] = A.GetLength(rnkA);
                }
                for(int i = 0; i < SumBCount; i++) {
                    int shift = i << 2;
                    uint mask = ((uint)(0xf)) << shift;

                    int rnkB = (int)((SumRankBStore & mask) >> shift);
                    int sumB = (int)((SumLoopBStore & mask) >> shift);
                    Debug.Assert(sumB >= 0 && sumB < NoOfSumCycles);
                    cycSumB[sumB] += B.GetCycle(rnkB);
                }
            }

            int SumACount;
            uint SumRankAStore;
            uint SumLoopAStore;

            int SumBCount;
            uint SumRankBStore;
            uint SumLoopBStore;

            int RunACount;
            uint RunRankAStore; // 4 bytes == 0x8 nibbles
            uint RunLoopAStore; // 4 bytes == 0x8 nibbles

            int RunBCount;
            uint RunRankBStore; // 4 bytes == 0x8 nibbles
            uint RunLoopBStore;

            internal int iTrafoIdx;

            /// <summary>
            /// Defines howto run across array 'A'.
            /// </summary>
            /// <param name="iLoop">Running loop index.</param>
            /// <param name="rankA">Rank index of array 'A' that corresponds with running loop <paramref name="iLoop"/>.</param>
            void AddRunRankA(int iLoop, int rankA) {
                if(rankA < 0 || rankA >= 0x8 || rankA >= DA)
                    throw new ArgumentOutOfRangeException("Rank index out of range, supporting at max. 8 ranks.");
                if(iLoop < 0 || iLoop >= DT)
                    throw new ArgumentOutOfRangeException();

                int shift = RunACount * 4;
                RunRankAStore |= (((uint)rankA) << shift);
                RunLoopAStore |= (((uint)iLoop) << shift);
                RunACount++;
            }

            /// <summary>
            /// Defines howto run across array 'B'.
            /// </summary>
            /// <param name="iLoop">Running loop index.</param>
            /// <param name="rankB">Rank index of array 'B' that corresponds with running loop <paramref name="iLoop"/>.</param>
            void AddRunRankB(int iLoop, int rankB) {
                if(rankB < 0 || rankB >= 0x10 || rankB >= DB)
                    throw new ArgumentOutOfRangeException("Rank index out of range, supporting at max. 8 ranks.");
                if(iLoop < 0 || iLoop >= DT)
                    throw new ArgumentOutOfRangeException();

                int shift = RunBCount * 4;
                RunRankBStore |= (((uint)rankB) << shift);
                RunLoopBStore |= (((uint)iLoop) << shift);
                RunBCount++;
            }

            void AddSumRankA(int iLoop, int rankA) {
                if(rankA < 0 || rankA >= 0x10)
                    throw new ArgumentOutOfRangeException();
                if(iLoop < 0 || iLoop >= MAX_SUM_LOOPS)
                    throw new ArgumentOutOfRangeException();

                int shift = SumACount * 4;
                SumRankAStore |= (((uint)rankA) << shift);
                SumLoopAStore |= (((uint)iLoop) << shift);
                SumACount++;
            }

            void AddSumRankB(int iLoop, int rankB) {
                if(rankB < 0 || rankB >= 0x10)
                    throw new ArgumentOutOfRangeException();
                if(iLoop < 0 || iLoop >= MAX_SUM_LOOPS)
                    throw new ArgumentOutOfRangeException();

                int shift = SumBCount * 4;
                SumRankBStore |= (((uint)rankB) << shift);
                SumLoopBStore |= (((uint)iLoop) << shift);
                SumBCount++;
            }

            unsafe static bool IdentifyTrafoIndex(ref char[] IdxNames, bool* useTrafo, char* TrafoName, char*TrafoVarName, bool* trafoIndexMarker) {
                bool r = false;
                int Len = IdxNames.Length;
                for(int i = 1; i < Len; i++) {
                    if(IdxNames[i] == '(') {
                        if(IdxNames.Length < (i + 3) || IdxNames[i + 2] != ')') {
                            throw new ArgumentException("something wrong with index transformation notation.");
                        }

                        if(!*useTrafo) {
                            *TrafoName = IdxNames[i - 1];
                            *TrafoVarName = IdxNames[i + 1];
                        } else {
                            if(*TrafoName != IdxNames[i - 1])
                                throw new ArgumentException("Only one index trafo is supported.");
                            if(*TrafoVarName != IdxNames[i + 1])
                                throw new ArgumentException("Index trafo can only be applied to one index.");
                        }
                        *useTrafo = true;
                        r = true;

                        Debug.Assert(trafoIndexMarker[i - 1] == false);
                        trafoIndexMarker[i - 1] = true;
                        IdxNames[i - 1] = IdxNames[i + 1]; // now, we have identified that the (i - 1) --th index is transformed; 
                        //                  since the pattern is something like "T(i)" from (i - 1) to (i + 2). 
                        //                  This info is stored now, we reduce the pattern to soley "i".

                        for(int k = i + 3; k < Len; k++) {
                            IdxNames[k - 3] = IdxNames[k];
                        }
                        Len -= 3;
                    }
                }
                if(Len != IdxNames.Length) {
                    Array.Resize(ref IdxNames, Len);
                }

                return r;
            }

            internal int TrfT0Sw;
            internal int TrfA0Sw;
            internal int TrfB0Sw;

            /// <summary>
            /// Compiles a 'program' to compute e.g. a tensor product
            /// like 
            /// \f[
            ///    T_{i j n} = \sum{k m} A_{i k m} B_{k j n m}.
            /// \f]
            /// In this example, <paramref name="Tindex"/> = "ijn", <paramref name="Aindex"/> = "ikm" and <paramref name="Bindex"/> = "kjnm".
            /// </summary>
            /// <param name="Tindex">Indices into result array.</param>
            /// <param name="Aindex">Indices into first tensor.</param>
            /// <param name="Bindex">Indices into second tensor.</param>
            /// <param name="detectTrafo">
            /// Enables index-transformation.
            /// </param>
            /// <returns></returns>
            public static MultiplyProgram Compile(string Tindex, string Aindex, string Bindex, bool detectTrafo = false) {
                MultiplyProgram p = default(MultiplyProgram);
                                

                unsafe {
                    
                    char[] TIdxNames = Tindex.ToCharArray();
                    char[] AIdxNames = Aindex.ToCharArray();
                    char[] BIdxNames = Bindex.ToCharArray();

                    char TrafoName = (char)0, TrafoVarName = (char)0;
                    bool useTrafo = false;
                    bool* TtrafoIndexMarker  = stackalloc bool[TIdxNames.Length];
                    bool* AtrafoIndexMarker  = stackalloc bool[Aindex.Length];
                    bool* BtrafoIndexMarker  = stackalloc bool[Bindex.Length];
                    if(detectTrafo) {
                        p.TrfT0Sw = IdentifyTrafoIndex(ref TIdxNames, &useTrafo, &TrafoName, &TrafoVarName, TtrafoIndexMarker) ? 1 : 0;
                        if(p.TrfT0Sw != 0)
                            throw new ArgumentException("Index transformation is not supported for the result array.");

                        p.TrfA0Sw = IdentifyTrafoIndex(ref AIdxNames, &useTrafo, &TrafoName, &TrafoVarName, AtrafoIndexMarker) ? 1 : 0;
                        p.TrfB0Sw = IdentifyTrafoIndex(ref BIdxNames, &useTrafo, &TrafoName, &TrafoVarName, BtrafoIndexMarker) ? 1 : 0;

                        if((p.TrfT0Sw * p.TrfA0Sw * p.TrfB0Sw) != 0)
                            throw new ArgumentException("Index-transformation can only be applied to one or two operands.");
                    }

                    p.DT = TIdxNames.Length;
                    p.DA = AIdxNames.Length;
                    p.DB = BIdxNames.Length;
                    
                    // check
                    p.iTrafoIdx = int.MinValue;
                    for(int i = 0; i < p.DT; i++) {
                        if(useTrafo && TIdxNames[i] == TrafoVarName) {
                            if(p.iTrafoIdx >= 0) {
                                throw new ArgumentException("Only one index of can be transformed.");
                            } else {
                                p.iTrafoIdx = i;
                            }
                        }

                        for(int j = i + 1; j < p.DT; j++) {
                            if(TIdxNames[i] == TIdxNames[j])
                                throw new ArgumentException("found non-unique running index name.", "ThisIndex");
                        }
                    }
                    
                    // Identify running and summation indices, setup running cycles (array A)
                    // ======================================================================
                    bool* bSumMarkerA = stackalloc bool[p.DA];
                    for(int j = 0; j < p.DA; j++) {
                        bSumMarkerA[j] = true;
                        for(int i = 0; i < p.DT; i++) {
                            if(AIdxNames[j] == TIdxNames[i]) {
                                //if(A.GetLength(j) != this.GetLength(i))
                                //    throw new ArgumentException(string.Format("length mismatch in running index '{0}', '{1}'-th index in this, '{2}'-th index in 'A'", AIdxNames[j], i, j), "A");
                                //cycRunA[i] += A.GetCycle(j);
                                
                                p.AddRunRankA(i, j);
                                
                                bSumMarkerA[j] = false;
                                break;
                            }
                        }

                        if(bSumMarkerA[j] && AtrafoIndexMarker[j])
                            throw new NotSupportedException("Index transformation is not supported for summation indices.");
                    }

                    // Identify running and summation indices, setup running cycles (array B)
                    // ======================================================================
                    bool* bSumMarkerB = stackalloc bool[p.DB];
                    for(int j = 0; j < p.DB; j++) { // loop over B-ranks
                        bSumMarkerB[j] = true;
                        for(int i = 0; i < p.DT; i++) { // loop over running indices (T-ranks)
                            if(BIdxNames[j] == TIdxNames[i]) {
                                //if(B.GetLength(j) != this.GetLength(i))
                                //    throw new ArgumentException(string.Format("length mismatch in running index '{0}', '{1}'-th index in this, '{2}'-th index in 'B'", BIdxNames[j], i, j), "B");
                                //cycRunB[i] += B.GetCycle(j);
                                
                                p.AddRunRankB(i, j);
                                
                                bSumMarkerB[j] = false;
                                break;
                            }
                        }

                        if(bSumMarkerB[j] && BtrafoIndexMarker[j])
                            throw new NotSupportedException("Index transformation is not supported for summation indices.");
                    }

                    // collect & check summation cycles
                    // ================================
                    p.NoOfSumCycles = 0;
                    char* SumIdxNames = stackalloc char[MAX_SUM_LOOPS];
                    
                    for(int j1 = 0; j1 < p.DA; j1++) {

                        if(bSumMarkerA[j1]) {

                            bool twiceOrMore = false;
                            int k;
                            for(k = 0; k < p.NoOfSumCycles; k++) {
                                if(SumIdxNames[k] == AIdxNames[j1]) {
                                    // summation index occurs more than once!
                                    //cycSumA[k] += A.GetCycle(j1);
                                    p.AddSumRankA(k, j1);
                                    twiceOrMore = true;
                                    break;
                                }
                            }

                            if(!twiceOrMore) {
                                if(k >= MAX_SUM_LOOPS)
                                    throw new NotSupportedException("More than '" + MAX_SUM_LOOPS + "' summation loops are currently not supported.");
                                SumIdxNames[k] = AIdxNames[j1];
                                //cycSumA[k] += A.GetCycle(j1);
                                //lenSum[k] = A.GetLength(j1);
                                p.AddSumRankA(k, j1);
                                p.NoOfSumCycles++;
                            }

                            bool bfound = false;
                            for(int j2 = 0; j2 < p.DB; j2++) {
                                if(!bSumMarkerB[j2])
                                    continue;

                                if(BIdxNames[j2] == AIdxNames[j1]) {
                                    bfound = true;

                                    //if(A.GetLength(j1) != B.GetLength(j2))
                                    //    throw new ArgumentException(string.Format("length mismatch for summation index '{0}'", AIdxNames[j1]));

                                    //cycSumB[k] += B.GetCycle(j2);
                                    p.AddSumRankB(k, j2);
                                    bSumMarkerB[j2] = false;
                                }

                            }

                            if(bfound == false && twiceOrMore == false)
                                throw new ArgumentException(string.Format("summation index '{0}' present in 'A', but not in 'B'", AIdxNames[j1]));
                        }
                    }

                    for(int j2 = 0; j2 < p.DB; j2++) {
                        if(bSumMarkerB[j2])
                            throw new ArgumentException(string.Format("summation index '{0}' present in 'B', but not in 'A'", BIdxNames[j2]));
                    }

                    return p;
                }
            }

        }

        /// <summary>
        /// Generalized tensor multiplication with index-transformation. 
        /// </summary>
        public void Multiply(double scale, MultidimensionalArray A, MultidimensionalArray B, double thisscale, ref MultiplyProgram mp, int[] IndexTrafo = null) {

            unsafe {
                fixed(int* pIndexTrafo = IndexTrafo) {
                    Multiply(scale, A, B, thisscale, ref mp, pIndexTrafo, IndexTrafo != null ? IndexTrafo.Length : 0);

                }
            }

        }

        /// <summary>
        /// Generalized tensor multiplication with index-transformation. 
        /// </summary>
        unsafe public void Multiply(double scale, MultidimensionalArray A, MultidimensionalArray B, double thisscale, 
            ref MultiplyProgram mp, int* IndexTrafo, int IndexTrafo_Length,
            int trfPreOffset_A = 0, int trfCycle_A = 1, int trfPostOffset_A = 0, int trfPreOffset_B = 0, int trfCycle_B = 1, int trfPostOffset_B = 0) {
            if(mp.DT != this.Dimension)
                throw new ArgumentException();
            if(mp.DA != A.Dimension)
                throw new ArgumentException();
            if(mp.DB != B.Dimension)
                throw new ArgumentException();

            unsafe {
                int DT = mp.DT;
                //int DA = mp.DA;
                //int DB = mp.DB;

                // running cycles:
                int* cycRunT = stackalloc int[3 * DT];
                int* cycRunA = cycRunT + DT;
                int* cycRunB = cycRunA + DT;

                int* lenRun = stackalloc int[DT];

                for(int i = 0; i < DT; i++) {
                    cycRunT[i] = this.GetCycle(i);

                    lenRun[i] = this.GetLength(i);
                }

                int* cycSumA = stackalloc int[MultiplyProgram.MAX_SUM_LOOPS];
                int* cycSumB = stackalloc int[MultiplyProgram.MAX_SUM_LOOPS];
                int* lenSum = stackalloc int[MultiplyProgram.MAX_SUM_LOOPS];

#if DEBUG
                mp.CheckArgs(this, A, B);
#endif
                
                mp.GetCyclesAndRanks(cycRunA, cycRunB, lenSum, cycSumA, cycSumB, A, B);
                
                if(mp.NoOfSumCycles == 2) {
                    // for better loop unrolling, make sure the inner loop is the smaller one
                    if(lenSum[1] > lenSum[0]) {
                        SwapInt(cycSumA + 0, cycSumA + 1);
                        SwapInt(cycSumB + 0, cycSumB + 1);
                        SwapInt(lenSum + 0, lenSum + 1);
                    }
                }


                // Execute Tensor Multiplication
                // =============================
                fixed(double* pTstor = this.m_Storage, pAstor = A.m_Storage, pBstor = B.m_Storage) {
                    double* pT = pTstor + this.m_Offset;
                    double* pA = pAstor + A.m_Offset;
                    double* pB = pBstor + B.m_Offset;

                    if(mp.iTrafoIdx >= 0) {
                        
                        if(mp.iTrafoIdx > 0) {
                            // transformed cycle MUST be the outer-most, i.e. the first one.
                            // => need to shift some cycles.

                            int kk = mp.iTrafoIdx;
                            SwapInt(lenRun + 0, lenRun + kk);
                            SwapInt(cycRunT + 0, cycRunT + kk);
                            SwapInt(cycRunA + 0, cycRunA + kk);
                            SwapInt(cycRunB + 0, cycRunB + kk);
                        }

                        
                        if(IndexTrafo == null)
                            throw new ArgumentException("Index transformation required.");

                        int* pIndexTrafo = IndexTrafo;
                        {
                            //Debug.Assert(mp.NoOfSumCycles == 1); Debug.Assert(DT == 2); __MultiplyWTrafo_Sum1_FOR2( 
                            MultiplyWTrafo_Dispatch(
                                DT, mp.NoOfSumCycles, pT, pA, pB, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scale, thisscale,
                                pIndexTrafo, IndexTrafo_Length, mp.TrfT0Sw, mp.TrfA0Sw, mp.TrfB0Sw,
                                0, 1, 0, trfPreOffset_A, trfCycle_A, trfPostOffset_A, trfPreOffset_B, trfCycle_B, trfPostOffset_B);
                        }



                    } else {


                        

                        Multiply_Dispatch(DT, mp.NoOfSumCycles, pT, pA, pB, lenRun, cycRunT, cycRunA, cycRunB, lenSum, cycSumA, cycSumB, scale, thisscale);
                    }
                }
            }
        }

        static unsafe void SwapInt(int* x, int* y) {
            if(x != y) {
                // XOR swap: bringt nichts, ist aber cool.
                *x ^= *y;
                *y ^= *x;
                *x ^= *y;
            }
        }


        /// <summary>
        /// Generalized tensor multiplication; 
        /// </summary>
        public void Multiply(double scale, MultidimensionalArray A, MultidimensionalArray B, double thisscale, string ThisIndex, string Aindex, string Bindex) {
            MultiplyProgram mp = MultiplyProgram.Compile(ThisIndex, Aindex, Bindex, false);
            mp.CheckArgs(this, A, B);
            this.Multiply(scale, A, B, thisscale, ref mp, null);
        }

        /*
        unsafe static private void __MultiplyWTrafo_Sum1_FOR2(int RunDim, int SumDim,
            double* pT_org, double* pA_org, double* pB_org, 
            int* lenRun, int* cycRunT, int* cycRunA, int* cycRunB, 
            int* lenSum, int* cycSumA, int* cycSumB, 
            double scl, double Tscl,
            int* Trf, int TrfEnd, int trfT0sw, int trfA0sw, int trfB0sw
            ) {
            Debug.Assert(RunDim == 2);
            int I0 = lenRun[0];
            int cTi0 = cycRunT[0];
            int cAi0 = cycRunA[0];
            int cBi0 = cycRunB[0];
            int I1 = lenRun[1];
            int cTi1 = cycRunT[1];
            int cAi1 = cycRunA[1];
            int cBi1 = cycRunB[1];
            Debug.Assert(SumDim == 1);
            int K0 = lenSum[0];
            int cAk0 = cycSumA[0];
            int cBk0 = cycSumB[0];
            //double* pT_org = pT;
            //double* pA_org = pA;
            //double* pB_org = pB;
            double* pT = pT_org, pA = pA_org, pB = pB_org;
				
            for(int i0 = 0; i0 < I0; i0++) {
                //double* pT_Old_i0 = pT;
                //double* pA_Old_i0 = pA;
                //double* pB_Old_i0 = pB;
                Debug.Assert(i0 < TrfEnd);
                int Trf_i0 = Trf[i0];

                //if(pT_org != pT || pA_org != pA || pB_org != pB) {
                //    throw new ApplicationException("____ i0=" + i0 + " T-diff" + (pT_org - pT) + " A-diff" + (pA_org - pA) + " B-diff" + (pB_org - pB));
                //}


                //pT += cTi0 * (i0 * (1 - trfT0sw) + Trf_i0 * trfT0sw);
                //pA += cAi0 * (i0 * (1 - trfA0sw) + Trf_i0 * trfA0sw);
                //pB += cBi0 * (i0 * (1 - trfB0sw) + Trf_i0 * trfB0sw);

                pT = pT_org + cTi0 * (i0 * (1 - trfT0sw) + Trf_i0 * trfT0sw);
                pA = pA_org + cAi0 * (i0 * (1 - trfA0sw) + Trf_i0 * trfA0sw);
                pB = pB_org + cBi0 * (i0 * (1 - trfB0sw) + Trf_i0 * trfB0sw);
                for(int i1 = 0; i1 < I1; i1++) {
                    {
                        double acc = 0.0;
                        double* pA_Old_k0 = pA;
                        double* pB_Old_k0 = pB;
                        for(int k0 = 0; k0 < K0; k0++) {
                            acc += (*pB) * (*pA);
                            pA += cAk0;
                            pB += cBk0;
                        }
                        *pT = acc * scl + (*pT) * Tscl;
                        pA = pA_Old_k0;
                        pB = pB_Old_k0;
                    }
                    pT += cTi1;
                    pA += cAi1;
                    pB += cBi1;
                }
                
                //pT = pT_Old_i0;
                //pA = pA_Old_i0;
                //pB = pB_Old_i0;
            }
        }*/
    
    
    }
}

