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

namespace ilPSP.LinSolvers.monkey.CL
{
    /// <summary>
    /// Base class for all OpenCL matrices
    /// </summary>
    abstract public class clMatrix : MatrixBase, IDisposable
    {
        /// <summary>
        /// used by <see cref="SpMV_External_RecvCallBack"/>;
        /// key: MPI processor rank; <br/>
        /// content: indices into <see cref="h_ElementsToAcc"/>;
        /// </summary>
        SortedDictionary<int, int[]> acc_ind = new SortedDictionary<int, int[]>();

        /// <summary>
        /// <see cref="h_ElementsToAcc"/>
        /// </summary>
        private int[] h_IndicesToAccumulate;

        /// <summary>
        /// Device memory pointer for <see cref="h_IndicesToAccumulate"/>
        /// </summary>
        private cl_mem d_IndicesToAccumulate;

        /// <summary>
        /// This are the elemnts of the 
        /// page-locked memory, allocated by 'cuMemAllocHost';<br/>
        /// Let be y = bata*y + alpha*M*x the SpMV, whith scalars alpha and beta, vectors x and y and matrix M;<br/>
        /// Then, the 
        /// </summary>
        private IntPtr h_ElementsToAcc;

        /// <summary>
        /// Device memory pointer for <see cref="h_ElementsToAcc"/>
        /// </summary>
        private cl_mem d_ElementsToAcc;

        /// <summary>
        /// Device memory pointer for memory of result vector
        /// </summary>
        private cl_mem d_acc;

        private double m_alpha = double.NaN;

        /// <summary>
        /// Number of external elements needed by this part
        /// </summary>
        private int extSize;

        private int extlocalsize = 128;

        private int extglobalsize;

        /// <summary>
        /// Size of matrix kernel workgroup
        /// </summary>
        protected int localsize = 256;

        /// <summary>
        /// Global size of matrix kernel
        /// </summary>
        protected int globalsize;

        /// <summary>
        /// Device
        /// </summary>
        protected clDevice device;

        /// <summary>
        /// Matrix multiply kernel
        /// </summary>
        protected cl_kernel clmultiply;

        private cl_kernel claccext;

        private cl_event localEvent;

        private bool disposed = false;

        private void LMAA()
        {
            List<int> ins_ind = new List<int>();
            foreach (var kv in ExtMatrix)
            {
                int proc = kv.Key;
                External ext = kv.Value;
                ins_ind.AddRange(ext.rowInd);
            }

            // remove duplicate entries
            ins_ind.Sort();
            if (ins_ind.Count > 1)
            {
                for (int i = 1; i < ins_ind.Count; i++)
                {
                    while (i < ins_ind.Count && ins_ind[i] == ins_ind[i - 1])
                        ins_ind.RemoveAt(i);
                }
            }

            //build new lists
            foreach (int proc in ExtMatrix.Keys)
            {
                External ext = ExtMatrix[proc];
                int[] accIndProc = new int[ext.rowInd.Length];
                for (int i = 0; i < accIndProc.Length; i++)
                {
                    accIndProc[i] = ins_ind.BinarySearch(ext.rowInd[i]);
                    if (accIndProc[i] < 0)
                        // should never happen
                        throw new ApplicationException("error in alg");
                }
                acc_ind.Add(proc, accIndProc);
            }

            h_IndicesToAccumulate = ins_ind.ToArray();
            extSize = h_IndicesToAccumulate.Length;
        }

        /// <summary>
        /// Create matrix
        /// </summary>
        /// <param name="M">Original matrix</param>
        /// <param name="device">Corresponding OpenCL device</param>
        /// <param name="kernelName">Name of the kernel function</param>
        public clMatrix(MsrMatrix M, clDevice device, string kernelName)
            : base(M)
        {
            this.device = device;
            base.PackMatrix(M);
            this.clmultiply = cl.CreateKernel(device.matrixProgram, kernelName);
            this.claccext = cl.CreateKernel(device.matrixProgram, "accumulateExternal");
            disposed = false;
            
            LMAA();

            if (extSize > 0)
            {
                extglobalsize = extSize;
                int m = extSize % extlocalsize;
                if (m > 0)
                {
                    extglobalsize += extlocalsize - m;
                }

                h_ElementsToAcc = Marshal.AllocHGlobal(extSize * sizeof(double));
                d_ElementsToAcc = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_READ_ONLY, (uint)extSize * sizeof(double));

                d_IndicesToAccumulate = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_READ_ONLY, (uint)extSize * sizeof(int));
                cl.EnqueueWriteBuffer(device.cq, d_IndicesToAccumulate, true, 0, (uint)extSize * sizeof(int), h_IndicesToAccumulate);
            }
        }

        /// <summary>
        /// Destructor
        /// </summary>
        ~clMatrix()
        {
            Dispose();
        }

        /// <summary>
        /// Free resources
        /// </summary>
        public override void Dispose()
        {
            if (!disposed)
            {
                base.Dispose();

                if (m_IsLocked)
                    Unlock();

                if (extSize > 0)
                {
                    cl.ReleaseMemObject(d_IndicesToAccumulate);
                    cl.ReleaseMemObject(d_ElementsToAcc);
                    cl.ReleaseKernel(claccext);
                    cl.ReleaseKernel(clmultiply);
                    Marshal.FreeHGlobal(h_ElementsToAcc);
                }

                disposed = true;
            }
        }

        /// <summary>
        /// see <see cref="MatrixBase.CreateVec{T}(T,IPartitioning,out bool)"/>
        /// </summary>
        protected override VectorBase CreateVec<T>(T a, IPartitioning len, out bool CopyIsShallow)
        {
            if (a.Count != len.LocalLength)
                throw new ArgumentException("count of a must be at least 'len'!", "len");
            if (len.MPI_Comm != this.ColPartition.MPI_Comm)
                throw new ArgumentException();


            IPartitioning part_a = len;
            double[] _a_stor = ArrayTools.List2Array<double>(a, 0, len.LocalLength);
            CopyIsShallow = false;
            return new clVector(part_a, _a_stor, device);
        }

        internal override void SpMV_Local_Start(double alpha, VectorBase a, double beta, VectorBase acc)
        {
            if (!m_IsLocked)
                throw new ApplicationException("object must be locked.");

            clVector _a = a as clVector;
            clVector _acc = acc as clVector;

            if (_a == null)
                throw new ArgumentException("a must be of type clVector.", "a");
            if (_acc == null)
                throw new ArgumentException("acc must be of type clVector.", "acc");

            SetArguments(alpha, _a, beta, _acc);

            int[] local = { localsize };
            int[] global = { globalsize };
            localEvent = cl.EnqueueNDRangeKernel(device.cq, clmultiply, 1, global, local);
        }

        abstract internal void SetArguments(double alpha, clVector a, double beta, clVector acc);

        /// <summary>
        /// no operation
        /// </summary>
        internal override void SpMV_Local_Middle(double alpha, VectorBase a, double beta, VectorBase acc)
        {
            if (!m_IsLocked)
                throw new ApplicationException("object must be locked.");
        }

        /// <summary>
        /// blocks until the computation of the internal part is finished.
        /// </summary>
        internal override void SpMV_Local_End(double alpha, VectorBase a, double beta, VectorBase acc)
        {
            if (!m_IsLocked)
                throw new ApplicationException("object must be locked.");

            //cl_event_info_return info = cl.GetEventInfo(localEvent);
            cl.WaitForEvent(localEvent);
        }

        internal override void SpMV_External_Begin(double alpha, double beta, VectorBase acc)
        {
            m_alpha = alpha;
            clVector _acc = (clVector)acc;
            d_acc = _acc.GetDevicePointer();

            unsafe
            {
                double* _acc_stor = (double*)h_ElementsToAcc;
                for (int i = (int)extSize - 1; i >= 0; i--)
                {
                    *_acc_stor = 0;
                    _acc_stor++;
                }
            }
        }

        /// <summary>
        /// multiplies the external vector (received from MPI-process <paramref name="procRank"/>)
        /// with the corresponding part of the matrix.
        /// This is always done on the CPU.
        /// The result is stored in the buffer <see cref="h_ElementsToAcc"/>
        /// </summary>
        internal override void SpMV_External_RecvCallBack(int procRank, IntPtr values)
        {
            External externalPart = ExtMatrix[procRank];

            int[] accIndProc = acc_ind[procRank];

            unsafe
            {
                double* _acc_stor = (double*)h_ElementsToAcc;
                double* _values = (double*)values;
                fixed (double* _Val = &externalPart.Val[0])
                {
                    fixed (int* _rowSt = &externalPart.RowStart[0],
                                _colInd = &externalPart.ColInd[0])
                    {
                        int K = externalPart.rowInd.Length;

                        int* Col = _colInd;
                        double* MtxEntry = _Val;

                        for (int k = 0; k < K; k++)
                        {
                            int RowStart = _rowSt[k];
                            int RowEnd = _rowSt[k + 1];

                            double rowacc = 0;
                            for (int i = RowStart; i < RowEnd; i++)
                            {
                                rowacc += (*MtxEntry) * _values[*Col];
                                Col++;
                                MtxEntry++;
                            }

                            _acc_stor[accIndProc[k]] += rowacc;
                        }
                    }
                }
            }
        }
        
        internal override void SpMV_External_Finalize()
        {
            if (extSize > 0)
            {
                cl.EnqueueWriteBuffer(device.cq, d_ElementsToAcc, true, 0, (uint)extSize * sizeof(double), h_ElementsToAcc);

                cl.SetKernelArg(claccext, 0, d_acc);
                cl.SetKernelArg(claccext, 1, d_IndicesToAccumulate);
                cl.SetKernelArg(claccext, 2, d_ElementsToAcc);
                cl.SetKernelArg(claccext, 3, m_alpha);
                cl.SetKernelArg(claccext, 4, extSize);

                int[] local = { extlocalsize };
                int[] global = { extglobalsize };

                cl.EnqueueNDRangeKernel(device.cq, claccext, 1, global, local);
            }
        }
    }
}
