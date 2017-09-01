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
using System.Text;
using System.Runtime.InteropServices;
using MPI.Wrappers;
using ilPSP.Utils;

namespace ilPSP.LinSolvers.monkey.mtCPU {
    
    /// <summary>
    /// Multi-threaded implementation of <see cref="MatrixBase"/>, that works completely on
    /// main processor/main memory.
    /// </summary>
    public class MtMatrix : MatrixBase {

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="M"></param>
        public MtMatrix(MsrMatrix M)
            : base(M) {
            base.PackMatrix(M);
        }
        
        override internal void SpMV_Local_Start(double alpha, VectorBase a, double beta, VectorBase acc) {
            // return immediately;
        }

        override internal void SpMV_Local_End(double alpha, VectorBase a, double beta, VectorBase acc) {
            // return immediately;
            // not needed in CPU implementation
        }

        override internal void SpMV_Local_Middle(double alpha, VectorBase a, double beta, VectorBase acc) {
            {
                MtVector __a = a as MtVector;
                MtVector __acc = acc as MtVector;
                if (__a == null)
                    throw new ArgumentException("a must be of type RefVector.", "a");
                if (__acc == null)
                    throw new ArgumentException("acc must be of type RefVector.", "acc");
            }


            //int[] TouchCount = new int[NoOfRows];

            // reference version:
            //{
            //    double[] _a_stor = _a.Storage;
            //    double[] _acc_stor = _acc.Storage;

            //    for (int Row = 0; Row < NoOfRows; Row++) {
            //        int RowStart = loaclPart.RowStart[Row];
            //        int RowEnd = loaclPart.RowStart[Row + 1];

            //        double rowacc = 0;
            //        for (int i = RowStart; i < RowEnd; i++) {
            //            int Col = loaclPart.ColInd[i];
            //            double MtxEntry = loaclPart.Val[i];

            //            rowacc += MtxEntry * _a_stor[Col];
            //        }

            //        _acc_stor[Row] = _acc_stor[Row] * beta + rowacc * alpha;
            //    }
            //}

            //// test code
            //double[] a_clone = (double[])_a.Storage.Clone();
            //double[] acc_clone = (double[])_acc.Storage.Clone();
            // test code end

                           // unsafe optimized version
            MtVector _a = a as MtVector;
            MtVector _acc = acc as MtVector;
            //MatrixBase.CSR loaclPart = (MatrixBase.CSR)this.m_LocalMtx;
            unsafe {
                double* _a_stor = _a.StorageAddr;
                double* _acc_stor = _acc.StorageAddr;
                int* _rowStart = this.LocalMatrixPin.pRowStart;
                int* ColInd = this.LocalMatrixPin.pColInd;
                double* _val = this.LocalMatrixPin.pVal;


                int NoOfRows = m_RowPart.LocalLength;
                ilPSP.Threading.Paralleism.For(0, NoOfRows, delegate(int i0, int iE) {

                    
                    double* MtxEntry = _val + _rowStart[i0];
                    int* __ColInd = ColInd + _rowStart[i0];
                    for (int Row = i0; Row < iE; Row++) {
                        int RowStart = _rowStart[Row];
                        int RowEnd = _rowStart[Row + 1];

                        double rowacc = 0;
                        for (int i = RowStart; i < RowEnd; i++) {
                            rowacc += *MtxEntry * _a_stor[*__ColInd];
                            __ColInd++;
                            MtxEntry++;
                        }

                        _acc_stor[Row] = _acc_stor[Row] * beta + rowacc * alpha;
                    }
                    

                });
            }

            //// Test code: test ok 05aug10, 16:36
            //base.m_LocalMtx.RefSpMv(alpha, a_clone, beta, acc_clone);
            //double err = 0;
            //for (int i = 0; i < this.RowPart.LocalLength; i++) {
            //    err += Math.Abs(acc_clone[i] - _acc.Storage[i]);
            //}
            //Console.WriteLine("err = " + err);
            //Console.WriteLine();
            //// test code end;
        }

        

        double[] m_acc_stor;
        double m_alpha = double.NaN;
        //double m_beta = double.NaN;

        override internal void SpMV_External_RecvCallBack(int procRank, IntPtr values) {
            External externalPart = ExtMatrix[procRank];

            // Reference version
            //{
            //    int K = externalPart.rowInd.Length;
            //    for (int k = 0; k < K; k++) {
            //        int Row = externalPart.rowInd[k];

            //        int RowStart = externalPart.RowStart[k];
            //        int RowEnd = externalPart.RowStart[k + 1];

            //        double rowacc = 0;
            //        for (int i = RowStart; i < RowEnd; i++) {
            //            int Col = externalPart.ColInd[i];
            //            double MtxEntry = externalPart.Val[i];

            //            rowacc += MtxEntry * values[Col];
            //        }

            //        m_acc_stor[Row] = m_acc_stor[Row] + rowacc * m_alpha;
            //    }
            //}

            // optimized version
            unsafe {
                double* _values = (double*)values;
                fixed (double* _Val = &externalPart.Val[0],
                               _acc_stor = &m_acc_stor[0]) {
                    fixed (int* _rowSt = &externalPart.RowStart[0],
                                _colInd = &externalPart.ColInd[0],
                                _rowInd = &externalPart.rowInd[0]) {

                        int K = externalPart.rowInd.Length;
                        double alpha = m_alpha;

                        int* Row = _rowInd;
                        int* Col = _colInd;
                        double* MtxEntry = _Val;

                        for (int k = 0; k < K; k++) {

                            int RowStart = _rowSt[k];
                            int RowEnd = _rowSt[k + 1];

                            double rowacc = 0;
                            for (int i = RowStart; i < RowEnd; i++) {
                                //int Col = externalPart.ColInd[i];
                                //double MtxEntry = externalPart.Val[i];

                                rowacc += *MtxEntry * _values[*Col];
                                Col++;
                                MtxEntry++;
                            }

                            //if (m_acc_stor[*Row] != 0.0)
                            //    throw new ApplicationException();

                            _acc_stor[*Row] += rowacc * alpha;
                            Row++;
                        }
                    }
                }
            }
        }

        

        /// <summary>
        /// 
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <param name="acc"></param>
        internal override void SpMV_External_Begin(double alpha, double beta, VectorBase acc) {
            m_alpha = alpha;
            //m_beta = beta;
            m_acc_stor = (acc as MtVector).Storage;
        }

        /// <summary>
        /// 
        /// </summary>
        internal override void SpMV_External_Finalize() {
            m_alpha = double.NaN;
            //m_beta = double.NaN;
            m_acc_stor = null;
        }

        /// <summary>
        /// see <see cref="MatrixBase.CreateVec"/>
        /// </summary>
        protected override VectorBase CreateVec<T>(T a, IPartitioning len, out bool CopyIsShallow) {
            if (a.Count != len.LocalLength)
                throw new ArgumentException("count of a must be at least 'len'!", "len");
            if (len.MPI_Comm != this.ColPartition.MPI_Comm)
                throw new ArgumentException();


            IPartitioning part_a = len;
            double[] _a_stor = a as double[]; // try to do a shallow copy and save mem. and time
            if (_a_stor == null) {
                // shallow copy didn't worked
                _a_stor = ArrayTools.List2Array<double>(a, 0, len.LocalLength);
                CopyIsShallow = false;
            } else {
                CopyIsShallow = true;
            }
            return new MtVector(part_a, _a_stor);
        }


        /// <summary>
        /// performs a GC pinning of matrix content
        /// </summary>
        public override void Lock() {
            base.Lock();
            LocalMatrixPin.Lock((MatrixBase.CSR)m_LocalMtx);
        }

        /// <summary>
        /// <see cref="LocalPin"/>
        /// </summary>
        LocalPin LocalMatrixPin;


        /// <summary>
        /// pinning for the local part of the matrix, within a lock/unlock section
        /// </summary>
        struct LocalPin {
            // GC pinning for Local part of the matrix
            GCHandle m_ValHandle;
            GCHandle m_ColIndHandle;
            GCHandle m_RowStartHandle;

            // pinned pointers to local part of the matrix, only valid in a look/unlock - section
            public unsafe double* pVal;
            public unsafe int* pColInd;
            public unsafe int* pRowStart;

            /// <summary>
            /// performs pinning and initializes the pointers
            /// </summary>
            /// <param name="LocalMatrix"></param>
            public void Lock(CSR LocalMatrix) {
                unsafe {
                    m_ValHandle = GCHandle.Alloc(LocalMatrix.Val, GCHandleType.Pinned);
                    pVal = (double*)Marshal.UnsafeAddrOfPinnedArrayElement(LocalMatrix.Val, 0);

                    m_ColIndHandle = GCHandle.Alloc(LocalMatrix.ColInd, GCHandleType.Pinned);
                    pColInd = (int*)Marshal.UnsafeAddrOfPinnedArrayElement(LocalMatrix.ColInd, 0);

                    m_RowStartHandle = GCHandle.Alloc(LocalMatrix.RowStart, GCHandleType.Pinned);
                    pRowStart = (int*)Marshal.UnsafeAddrOfPinnedArrayElement(LocalMatrix.RowStart, 0);
                }

            }

            /// <summary>
            /// releases the pinning, unlocks the pointers.
            /// </summary>
            public void Unlock() {
                m_ValHandle.Free();
                m_ColIndHandle.Free();
                m_RowStartHandle.Free();
                unsafe {
                    pVal = (double*)IntPtr.Zero;
                    pColInd = (int*)IntPtr.Zero;
                    pRowStart = (int*)IntPtr.Zero;
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override void Unlock() {
            base.Unlock();
            LocalMatrixPin.Unlock();
        }

        /// <summary>
        /// returns a <see cref="MatrixBase.CSR"/>-object.
        /// </summary>
        protected override MatrixBase.FormatBase AssembleFinalFormat(MatrixBase.TempCSR tmp) {
            //TestFormats(tmp);
            return new MatrixBase.CSR(tmp);
        }

        

        /// <summary>
        /// the poor Man's test unit
        /// </summary>
        void TestFormats(MatrixBase.TempCSR _tmp) {
            MatrixBase.CSR comp1 = new CSR(_tmp);
            //MatrixBase.CCBCSR comp2 = new CCBCSR(_tmp, 128, CellSize);
            //MatrixBase.ELLPACKmod comp2 = new ELLPACKmod(_tmp, 40, 4);
            MatrixBase.ManualCacheELLPACK comp2 = new ManualCacheELLPACK(_tmp, 1, 4, 2);
            
            double[] x = new double[this.ColPartition.LocalLength];
            double[] y = new double[this.ColPartition.LocalLength];

            double alpha = 1.234;
            double beta = 432.1;

            Random r = new Random(0);
            for (int i = 0; i < x.Length; i++)
                x[i] = r.NextDouble();
            for (int i = 0; i < y.Length; i++)
                y[i] = r.NextDouble();


            double[] res1 = (double[])y.Clone();
            comp1.RefSpMv(alpha, (double[])x.Clone(), beta, res1);
            
            double[] res2 = (double[])y.Clone();
            comp2.RefSpMv(alpha, (double[])x.Clone(), beta, res2);

            double err = 0;
            for( int i = 0; i < y.Length; i++) {
                err += Math.Abs(res1[i] - res2[i]);
            }

            Console.WriteLine("err = " + err);
            Console.Write("");
        }
         
    }
}
