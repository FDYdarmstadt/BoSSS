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

namespace ilPSP.LinSolvers.monkey.CL
{
    /// <summary>
    /// OpenCL vector class
    /// </summary>
    public class clVector : VectorBase
    {
        private cl_kernel clscale, clacc, clmew, cldnrm2, clinnerprod;
        private clDevice device;

        private int size;
        private int groups;
        private int localsize = 256;
        private int globalsize;
        private int globalsizehalf;

        private IntPtr h_result;
        private double[] h_data;
        private cl_mem d_result;
        private cl_mem d_data;

        /// <summary>
        /// Create OpenCL vector
        /// </summary>
        /// <param name="p">Parition</param>
        /// <param name="device">Device</param>
        public clVector(IPartitioning p, clDevice device)
            : base(p)
        {
            h_data = new double[p.LocalLength];
            this.device = device;
            init(p, device.vectorProgram);
        }

        /// <summary>
        /// Create OpenCL vector with external memory
        /// </summary>
        /// <param name="p">Parition</param>
        /// <param name="content">Memory for this vector</param>
        /// <param name="device">Device</param>
        public clVector(IPartitioning p, double[] content, clDevice device)
            : base(p)
        {
            h_data = content;
            this.device = device;
            init(p, device.vectorProgram);
        }

        /// <summary>
        /// Destructor
        /// </summary>
        ~clVector()
        {
            if (this.IsLocked)
            {
                Unlock();
            }

            cl.ReleaseKernel(clscale);
            cl.ReleaseKernel(clacc);
            cl.ReleaseKernel(clmew);
            cl.ReleaseKernel(cldnrm2);
            cl.ReleaseKernel(clinnerprod);
        }

        private void init(IPartitioning p, cl_program program)
        {
            clscale = cl.CreateKernel(program, "scale");
            clacc = cl.CreateKernel(program, "acc");
            clmew = cl.CreateKernel(program, "mew");
            cldnrm2 = cl.CreateKernel(program, "dnrm2");
            clinnerprod = cl.CreateKernel(program, "innerprod");

            size = p.LocalLength;
            globalsize = size;
            int m = globalsize % localsize;
            if (m > 0)
            {
                globalsize += localsize - m;
            }
            
            globalsizehalf = globalsize / 2;
            m = globalsizehalf % localsize;
            if (m > 0)
            {
                globalsizehalf += localsize - m;
            }

            groups = globalsizehalf / localsize;
        }

        /// <summary>
        /// Returns the device pointer of this vector.
        /// </summary>
        /// <returns>Returns the device pointer.</returns>
        public cl_mem GetDevicePointer()
        {
            if (!this.IsLocked)
                throw new ApplicationException("works only in locked mode");

            return d_data;
        }

        /// <summary>
        /// Set some values in the vector.
        /// </summary>
        /// Vector must not be locked.
        /// <typeparam name="T">Array or kind of IList of doubles</typeparam>
        /// <param name="vals">Values</param>
        /// <param name="arrayIndex">Offset for values</param>
        /// <param name="insertAt">Offset for this vector</param>
        /// <param name="Length">Number of values</param>
        public override void SetValues<T>(T vals, int arrayIndex, int insertAt, int Length)
        {
            if (this.IsLocked)
                throw new ApplicationException("object is locked.");

            if (typeof(T) == typeof(double[]))
            {
                // optimized version
                Array.Copy(vals as double[], arrayIndex, h_data, insertAt, Length);
            }
            else
            {
                // version which works whith all kinds of IList<double>
                for (int i = 0; i < Length; i++)
                    h_data[i + insertAt] = vals[i + arrayIndex];
            }
        }

        /// <summary>
        /// Get some values from the vector.
        /// </summary>
        /// Vector must not be locked.
        /// <typeparam name="T">Array or kind of IList of doubles</typeparam>
        /// <param name="vals">Values</param>
        /// <param name="arrayIndex">Offset for values</param>
        /// <param name="readAt">Offset for this vector</param>
        /// <param name="Length">Number of values</param>
        public override void GetValues<T>(T vals, int arrayIndex, int readAt, int Length)
        {
            if (this.IsLocked)
                throw new ApplicationException("object is locked.");

            if (typeof(T) == typeof(double[]))
            {
                // optimized version
                Array.Copy(h_data, readAt, vals as double[], arrayIndex, Length);
            }
            else
            {
                // version which works whith all kinds of IList<double>
                for (int i = 0; i < Length; i++)
                    vals[i + arrayIndex] = h_data[i + readAt];
            }
        }

        /// <summary>
        /// Get vector element at given index.
        /// </summary>
        /// Vector must not be locked.
        /// <param name="index">Index</param>
        /// <returns>The element</returns>
        public override double this[int index]
        {
            get
            {
                if (this.IsLocked)
                    throw new ApplicationException("object is locked.");

                return h_data[index];
            }
            set
            {
                if (this.IsLocked)
                    throw new ApplicationException("object is locked.");

                h_data[index] = value;
            }
        }

        /// <summary>
        /// Set vector values to zero
        /// </summary>
        public override void Clear()
        {
            
            Array.Clear(h_data, 0, h_data.Length);
            
            if (this.IsLocked)
            {
                cl.EnqueueWriteBuffer(device.cq, d_data, true, 0, (uint)h_data.Length * sizeof(double), h_data);
            }
        }

        /// <summary>
        /// Scale vector with scalar alpha
        /// </summary>
        /// <param name="alpha">Scaling factor</param>
        public override void Scale(double alpha)
        {
            if (!this.IsLocked)
                throw new ApplicationException("works only in locked mode");

            cl.SetKernelArg(clscale, 0, d_data);
            cl.SetKernelArg(clscale, 1, alpha);
            cl.SetKernelArg(clscale, 2, size);

            int[] global = { globalsize };
            int[] local = { localsize };
            cl.EnqueueNDRangeKernel(device.cq, clscale, 1, global, local);
        }

        /// <summary>
        /// Accumulate vector element-wise with other vector scaled by alpha
        /// </summary>
        /// <param name="alpha">Scaling for other</param>
        /// <param name="other">Other vector</param>
        public override void Acc(double alpha, VectorBase other)
        {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            clVector _other = other as clVector;
            if (_other == null)
                throw new ArgumentException("other must be of type clVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            cl.SetKernelArg(clacc, 0, d_data);
            cl.SetKernelArg(clacc, 1, _other.GetDevicePointer());
            cl.SetKernelArg(clacc, 2, alpha);
            cl.SetKernelArg(clacc, 3, size);

            int[] global = { globalsize };
            int[] local = { localsize };
            cl.EnqueueNDRangeKernel(device.cq, clacc, 1, global, local);
        }

        /// <summary>
        /// Multiply with other vector element-wise
        /// </summary>
        /// <param name="other">Other vector</param>
        public override void MultiplyElementWise(VectorBase other)
        {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            clVector _other = other as clVector;
            if (_other == null)
                throw new ArgumentException("other must be of type clVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            cl.SetKernelArg(clmew, 0, d_data);
            cl.SetKernelArg(clmew, 1, _other.GetDevicePointer());
            cl.SetKernelArg(clmew, 2, size);

            int[] global = { globalsize };
            int[] local = { localsize };
            cl.EnqueueNDRangeKernel(device.cq, clmew, 1, global, local);
        }

        /// <summary>
        /// Get the squared length of this vector
        /// </summary>
        /// <returns>Returns the squared length</returns>
        public override double TwoNormSquare()
        {
            if (!this.IsLocked)
                throw new ApplicationException("works only in locked mode");

            double TwoNormSquareLocal = 0.0;

            cl.SetKernelArg(cldnrm2, 0, d_data);
            cl.SetKernelArg(cldnrm2, 1, d_result);
            cl.SetKernelArgLocalSize(cldnrm2, 2, (uint)localsize * sizeof(double));
            cl.SetKernelArg(cldnrm2, 3, size);

            int[] global = { globalsizehalf };
            int[] local = { localsize };
            cl.EnqueueNDRangeKernel(device.cq, cldnrm2, 1, global, local);

            IntPtr h_result;
            cl.EnqueueMapBuffer(device.cq, d_result, out h_result, true, cl_map_flags.CL_MAP_READ, 0, (uint)groups * sizeof(double));

            unsafe
            {
                double* ptr = (double*)h_result;
                for (int i = 0; i < groups; i++)
                {
                    TwoNormSquareLocal += ptr[i];
                }
            }

            cl.EnqueueUnmapMemObject(device.cq, d_result, h_result);

            double TwoNormSquareGlobal = double.NaN;
            unsafe
            {
                csMPI.Raw.Allreduce((IntPtr)(&TwoNormSquareLocal), (IntPtr)(&TwoNormSquareGlobal), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            }

            return TwoNormSquareGlobal;
        }

        /// <summary>
        /// Get the inner product with other vector
        /// </summary>
        /// <param name="other">Other vector</param>
        /// <returns>Returns the inner product</returns>
        public override double InnerProd(VectorBase other)
        {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            clVector _other = other as clVector;
            if (_other == null)
                throw new ArgumentException("other must be of type clVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            double InnerProdLocal = 0.0;

            cl.SetKernelArg(clinnerprod, 0, d_data);
            cl.SetKernelArg(clinnerprod, 1, _other.GetDevicePointer());
            cl.SetKernelArg(clinnerprod, 2, d_result);
            cl.SetKernelArgLocalSize(clinnerprod, 3, (uint)localsize * sizeof(double));
            cl.SetKernelArg(clinnerprod, 4, size);

            int[] global = { globalsizehalf };
            int[] local = { localsize };
            cl.EnqueueNDRangeKernel(device.cq, clinnerprod, 1, global, local);

            IntPtr h_result;
            cl.EnqueueMapBuffer(device.cq, d_result, out h_result, true, cl_map_flags.CL_MAP_READ, 0, (uint)groups * sizeof(double));
            
            unsafe
            {
                double* ptr = (double*)h_result;
                for (int i = 0; i < groups; i++)
                {
                    InnerProdLocal += ptr[i];
                }
            }

            cl.EnqueueUnmapMemObject(device.cq, d_result, h_result);

            double InnerProdGlobal = double.NaN;
            unsafe
            {
                csMPI.Raw.Allreduce((IntPtr)(&InnerProdLocal), (IntPtr)(&InnerProdGlobal), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            }

            return InnerProdGlobal;
        }

        /// <summary>
        /// Copy content from other vector
        /// </summary>
        /// <param name="other">Other vector</param>
        public override void CopyFrom(VectorBase other)
        {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            clVector _other = other as clVector;
            if (_other == null)
                throw new ArgumentException("other must be of type clVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            cl.EnqueueCopyBuffer(device.cq, _other.GetDevicePointer(), d_data, 0, 0, (uint)(size * sizeof(double)));
        }

        /// <summary>
        /// Swap contents with other vector
        /// </summary>
        /// <param name="other">Other vector</param>
        public override void Swap(VectorBase other)
        {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            clVector _other = other as clVector;
            if (_other == null)
                throw new ArgumentException("other must be of type clVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            cl_mem temp = _other.d_data;
            _other.d_data = this.d_data;
            this.d_data = temp;
        }

        /// <summary>
        /// Copy vector to device
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            h_result = Marshal.AllocHGlobal(groups * sizeof(double));
            d_data = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR, (uint)h_data.Length * sizeof(double), h_data);
            d_result = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_USE_HOST_PTR | cl_mem_flags.CL_MEM_WRITE_ONLY, (uint)groups * sizeof(double), h_result);
        }

        /// <summary>
        /// Copy back from device
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cl.EnqueueReadBuffer(device.cq, d_data, true, 0, (uint)h_data.Length * sizeof(double), h_data);
            cl.ReleaseMemObject(d_data);
            cl.ReleaseMemObject(d_result);
            Marshal.FreeHGlobal(h_result);
        }

        /// <summary>
        /// see <see cref="VectorBase.CopyPartlyFrom"/>
        /// </summary>
        public override void CopyPartlyFrom(VectorBase _src, int[] IdxThis, int PerThis, int[] IdxSrc, int PerSrc)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Create communication vector class
        /// </summary>
        /// <param name="M">Matrix for communication</param>
        /// <returns>Returns the commvector class</returns>
        public override CommVector CreateCommVector(MatrixBase M)
        {
            return new clCommVector(M, this);
        }

        class clCommVector : CommVector
        {
            private cl_kernel clfill;

            /// <summary>
            /// send buffer
            /// </summary>
            private IntPtr h_SendBuffer;

            /// <summary>
            /// Device memory pointer for <see cref="h_SendBuffer"/>
            /// </summary>
            private cl_mem d_SendBuffer;

            /// <summary>
            /// for the SpMv y=beta*y + alpha*M*x, the indices of those elements of x which must 
            /// be send to other processors.
            /// </summary>
            private int[] h_IndicesToSend;

            /// <summary>
            /// Device memory pointer for <see cref="h_IndicesToSend"/>
            /// </summary>
            private cl_mem d_IndicesToSend;

            /// <summary>
            /// clVector that this CommVector belongs to
            /// </summary>
            private clVector owner;

            private cl_event sendbufferevent;

            private int size;
            private int localsize = 256;
            private int globalsize;
            private bool disposed = false;

            internal clCommVector(MatrixBase M, clVector v)
                : base(M, v)
            {
                this.owner = v;

                clfill = cl.CreateKernel(owner.device.vectorProgram, "fillSendBuffer");

                IDictionary<int, int[]> comLists = M._SpmvCommPattern.ComLists;
                //int[] procranks = new int[comLists.Count]; // put all proccessor ranks in one list to have a unique ordering

                int totLen = 0;
                foreach (int procRnk in comLists.Keys)
                {
                    int l = comLists[procRnk].Length;
                    base.SendBuffersLengths[procRnk] = l;
                    totLen += l;
                }

                size = totLen;
                globalsize = size;
                int m = size % localsize;
                if (m > 0)
                {
                    globalsize += localsize - m;
                }
                
                if (size > 0)
                {
                    // alloc
                    h_IndicesToSend = new int[size];
                    d_IndicesToSend = cl.CreateBuffer(owner.device.env.context, cl_mem_flags.CL_MEM_READ_ONLY, (uint)size * sizeof(int));

                    h_SendBuffer = Marshal.AllocHGlobal(size * sizeof(double));
                    d_SendBuffer = cl.CreateBuffer(owner.device.env.context, cl_mem_flags.CL_MEM_WRITE_ONLY, (uint)size * sizeof(double));

                    // concat lists: 
                    int i0 = 0;
                    unsafe
                    {
                        double* P0 = (double*)h_SendBuffer;

                        foreach (int procRnk in comLists.Keys)
                        {
                            base.SendBuffers[procRnk] = (IntPtr)P0;  // startaddres for sending to process 'procRnk'

                            int l = base.SendBuffersLengths[procRnk];
                            P0 += l;
                            Array.Copy(comLists[procRnk], 0, h_IndicesToSend, i0, l); // concat comm list
                            i0 += l;
                        }
                    }

                    cl.EnqueueWriteBuffer(owner.device.cq, d_IndicesToSend, true, 0, (uint)size * sizeof(int), h_IndicesToSend);
                }
            }

            ~clCommVector()
            {
                Dispose();
            }

            public override void Dispose()
            {
                if (!disposed)
                {
                    base.Dispose();

                    if (size > 0)
                    {
                        cl.ReleaseMemObject(d_IndicesToSend);
                        cl.ReleaseMemObject(d_SendBuffer);
                        Marshal.FreeHGlobal(h_SendBuffer);
                    }

                    disposed = true;
                }
            }

            public override void FillSendBuffer()
            {
                if (!owner.IsLocked)
                    throw new ApplicationException("works only in locked mode");

                base.FillSendBuffer();

                if (size > 0)
                {
                    cl.SetKernelArg(clfill, 0, d_SendBuffer);
                    cl.SetKernelArg(clfill, 1, d_IndicesToSend);
                    cl.SetKernelArg(clfill, 2, owner.GetDevicePointer());
                    cl.SetKernelArg(clfill, 3, size);

                    int[] local = { localsize };
                    int[] global = { globalsize };
                    cl.EnqueueNDRangeKernel(owner.device.cq, clfill, 1, global, local);
                    
                    sendbufferevent = cl.EnqueueReadBuffer(owner.device.cq, d_SendBuffer, false, 0, (uint)size * sizeof(double), h_SendBuffer);
                }
            }

            public override void StartTransmissionImReturn()
            {
                if (size > 0)
                {
                    cl.WaitForEvent(sendbufferevent);
                }
                base.StartTransmissionImReturn();
            }
        }
    }
}
