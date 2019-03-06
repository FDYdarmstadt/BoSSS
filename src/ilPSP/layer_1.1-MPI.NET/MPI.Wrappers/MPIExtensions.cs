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
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;

namespace MPI.Wrappers {

    /// <summary>
    /// extension methods which perform MPI tasks
    /// </summary>
    public static class MPIExtensions {

        /// <summary>
        /// equal to <see cref="MPIBroadcast{T}(T, int, MPI_Comm)"/>, acting on
        /// the WORLD-communicator
        /// </summary>
        static public T MPIBroadcast<T>(this T o, int root) {
            return MPIBroadcast(o, root, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// broadcasts an object to all other processes in the 
        /// MPI communicator <paramref name="comm"/>.
        /// </summary>
        /// <param name="o">
        /// an arbitrary, serialize able object;
        /// ignored on all processes which are not <paramref name="root"/>
        /// </param>
        /// <param name="root">
        /// rank of the sender process
        /// </param>
        /// <param name="comm">MPI communicator</param>
        /// <returns>
        /// on the sender process, the input <paramref name="o"/> is returned;
        /// </returns>
        /// <typeparam name="T">
        /// type of the object to transmit
        /// </typeparam>
        static public T MPIBroadcast<T>(this T o, int root, MPI_Comm comm) {

            int MyRank;
            csMPI.Raw.Comm_Rank(comm, out MyRank);

            IFormatter _Formatter = new BinaryFormatter();

            // -----------------------------------------------------
            // 1st phase: serialize object and broadcast object size 
            // -----------------------------------------------------

            byte[] buffer = null;
            int Size = -1;
            MemoryStream ms = null;
            if (root == MyRank) {

                ms = new MemoryStream();
                _Formatter.Serialize(ms, o);
                Size = (int)ms.Position;
                buffer = ms.GetBuffer();
            }

            unsafe {
                csMPI.Raw.Bcast((IntPtr)(&Size), 4, csMPI.Raw._DATATYPE.BYTE, root, comm);
            }

            // ---------------------------
            // 2nd phase: broadcast object
            // ---------------------------

            if (buffer == null) {
                buffer = new byte[Size];
            }

            unsafe {
                fixed (byte* pBuffer = buffer) {
                    csMPI.Raw.Bcast((IntPtr)pBuffer, Size, csMPI.Raw._DATATYPE.BYTE, root, comm);
                }
            }

            if (MyRank == root) {
                ms.Dispose();
                return o;
            } else {
                T r;
                ms = new MemoryStream(buffer);
                r = (T)_Formatter.Deserialize(ms);
                ms.Dispose();
                return r;
            }
        }

        /// <summary>
        /// if <paramref name="e"/> is unequal to null on any of the calling
        /// MPI processes, this method broadcasts it to all other processes; 
        /// </summary>
        /// <param name="e">
        /// an exception, or null
        /// </param>
        /// <param name="comm">
        /// %
        /// </param>
        /// <remarks>
        /// The following code may cause a deadlock:
        /// <code>
        /// try {
        ///     // some code that may cause an exception on some MPI process,
        ///     // but not on all of them
        /// } catch(Exception e) {
        ///     // some statement that have an effect on control flow, e.g. 'continue;' or 
        ///     return false;
        /// }
        /// // some collective call, e.g. 
        /// csMPI.Raw.Allreduce( ... )
        /// </code>
        /// The collective call <c>csMPI.Raw.Allreduce( ... )</c> is never executed on some 
        /// MPI processes, therefore the application will ang in a deadlock.
        /// <code>
        /// try {
        ///     Exception e = null;
        ///     try { 
        ///         // some code that may cause an exception on some MPI process,
        ///         // but not on all of them
        ///     } catch(Exception ee) {
        ///         e = ee;
        ///     }
        ///     ExceptionBcast(e,csMPI.Raw._COMM.WORLD);
        /// } catch(Exception e2) {
        ///     // some statement that have an effect on control flow, e.g. 'continue;' or 
        ///     return false;
        /// }
        /// // some collective call, e.g. 
        /// csMPI.Raw.Allreduce( ... )
        /// </code>
        /// Note that, if <paramref name="e"/> is not serializable,
        /// only the message of <paramref name="e"/> will be broad-casted and
        /// wrapped into a new exception.
        /// </remarks>
        static public void ExceptionBcast(this Exception e, MPI_Comm comm) {
            int ExcSrc = int.MaxValue;
            ////Console.WriteLine("Reminder: test code in exception broadcast");
            //if (e != null)
            //    throw e;
            //else
            //    return;

            if (e != null)
                csMPI.Raw.Comm_Rank(comm, out ExcSrc);

            unsafe {
                int res = 0;
                csMPI.Raw.Allreduce((IntPtr)(&ExcSrc), (IntPtr)(&res), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, comm);
                ExcSrc = res;
            }


            if (ExcSrc < int.MaxValue) {
                // an exception occured on some process

                int myRank;
                csMPI.Raw.Comm_Rank(comm, out myRank);

                object reduced;

                if (myRank == ExcSrc) {
                    // sender branch


                    if (e.GetType().GetCustomAttributes(typeof(SerializableAttribute), false).Length > 0)
                        // exception is serializeable -> bcast exception itself
                        reduced = MPIBroadcast<object>(e, ExcSrc, comm);
                    else
                        // bcast exception message
                        reduced = MPIBroadcast<object>(e.GetType().Name + ": '" + e.Message + "'", ExcSrc, comm);
                } else {
                    // receiver branch

                    reduced = MPIBroadcast<object>(null, ExcSrc, comm);
                }

                if (reduced is string) {
                    throw new ApplicationException("On MPI Process #" + myRank + ": " + ((string)reduced));
                } else if (reduced is Exception) {
                    throw ((Exception)reduced);
                } else {
                    throw new ApplicationException("should never occur");
                }
            }
        }

        /// <summary>
        /// equal to <see cref="ExceptionBcast(Exception,MPI_Comm)"/>, acting
        /// on the WORLD-communicator
        /// </summary>
        /// <param name="e"></param>
        static public void ExceptionBcast(this Exception e) {
            ExceptionBcast(e, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// equal to <see cref="MPISum(int,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public int MPISum(this int i) {
            return MPISum(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the sum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public int MPISum(this int i, MPI_Comm comm) {
            int loc = i;
            unsafe {
                int glob = int.MinValue;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
                return glob;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIMax(double[],MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public double[] MPIMax(this double[] i) {
            return MPIMax(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the maximum of each entry of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public double[] MPIMax(this double[] i, MPI_Comm comm) {
            double[] R = new double[i.Length];
            unsafe {
                fixed (double* loc = i, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), i.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, comm);
                }
            }
            return R;
        }

        /// <summary>
        /// equal to <see cref="MPIMin(double[],MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public double[] MPIMin(this double[] i) {
            return MPIMin(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the minimum of each entry of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public double[] MPIMin(this double[] i, MPI_Comm comm) {
            double[] R = new double[i.Length];
            unsafe {
                fixed (double* loc = i, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), i.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, comm);
                }
            }
            return R;
        }




        /// <summary>
        /// equal to <see cref="MPIMax(int[],MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public int[] MPIMax(this int[] i) {
            return MPIMax(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the maximum of each entry of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public int[] MPIMax(this int[] i, MPI_Comm comm) {
            int[] R = new int[i.Length];
            unsafe {
                fixed (int* loc = i, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), i.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MAX, comm);
                }
            }
            return R;
        }

        /// <summary>
        /// equal to <see cref="MPIMin(int[],MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public int[] MPIMin(this int[] i) {
            return MPIMin(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the minimum of each entry of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public int[] MPIMin(this int[] i, MPI_Comm comm) {
            int[] R = new int[i.Length];
            unsafe {
                fixed (int* loc = i, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), i.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, comm);
                }
            }
            return R;
        }


        /// <summary>
        /// Equal to <see cref="MPISum(int[],MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public int[] MPISum(this int[] A) {
            return MPISum(A, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Returns, for each entry, the sum <paramref name="A"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public int[] MPISum(this int[] A, MPI_Comm comm) {
            int[] S = new int[A.Length];
            unsafe {
                fixed (int* pA = A, pS = S) {
                    csMPI.Raw.Allreduce(((IntPtr)(pA)), ((IntPtr)(pS)), A.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
                }
                return S;
            }
        }


        /// <summary>
        /// Equal to <see cref="MPISum(double[],MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public double[] MPISum(this double[] A) {
            return MPISum(A, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Returns, for each entry, the sum <paramref name="A"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public double[] MPISum(this double[] A, MPI_Comm comm) {
            double[] S = new double[A.Length];
            unsafe {
                fixed (double* pA = A, pS = S) {
                    csMPI.Raw.Allreduce(((IntPtr)(pA)), ((IntPtr)(pS)), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, comm);
                }
                return S;
            }
        }


        /// <summary>
        /// equal to <see cref="MPISum(double,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public double MPISum(this double i) {
            return MPISum(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the sum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public double MPISum(this double i, MPI_Comm comm) {
            double loc = i;
            unsafe {
                double glob = double.NaN;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, comm);
                return glob;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIMin(double,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public double MPIMin(this double i) {
            return MPIMin(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the minimum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public double MPIMin(this double i, MPI_Comm comm) {
            double loc = i;
            unsafe {
                double glob = double.NaN;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, comm);
                return glob;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIMax(double,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public double MPIMax(this double i) {
            return MPIMax(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the maximum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public double MPIMax(this double i, MPI_Comm comm) {
            double loc = i;
            unsafe {
                double glob = double.NaN;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, comm);
                return glob;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIMin(int,MPI_Comm)"/>, acting on the WORLD-communicator
        /// </summary>
        static public int MPIMin(this int i) {
            return MPIMin(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the minimum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public int MPIMin(this int i, MPI_Comm comm) {
            int loc = i;
            unsafe {
                int glob = int.MaxValue;
                csMPI.Raw.Allreduce(
                    (IntPtr)(&loc),
                    (IntPtr)(&glob),
                    1,
                    csMPI.Raw._DATATYPE.INT,
                    csMPI.Raw._OP.MIN,
                    comm);
                return glob;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIMax(int,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public int MPIMax(this int i) {
            return MPIMax(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the maximum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public int MPIMax(this int i, MPI_Comm comm) {
            int loc = i;
            unsafe {
                int glob = int.MinValue;
                csMPI.Raw.Allreduce(
                    (IntPtr)(&loc),
                    (IntPtr)(&glob),
                    1,
                    csMPI.Raw._DATATYPE.INT,
                    csMPI.Raw._OP.MAX,
                    comm);
                return glob;
            }
        }

        static public int[] MPIAllGather(this int i) {
            return i.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        static public int[] MPIAllGather(this int i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);

            int[] result = new int[size];
            unsafe {
                int sendBuffer = i;
                fixed (int* pResult = &result[0]) {
                    csMPI.Raw.Allgather(
                        (IntPtr)(&i),
                        1,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pResult,
                        1,
                        csMPI.Raw._DATATYPE.INT,
                        comm);
                }
            }

            return result;
        }

        static public double[] MPIAllGather(this double d) {
            return d.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        static public double[] MPIAllGather(this double d, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);

            double[] result = new double[size];
            unsafe {
                double sendBuffer = d;
                fixed (double* pResult = &result[0]) {
                    csMPI.Raw.Allgather(
                        (IntPtr)(&d),
                        1,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        (IntPtr)pResult,
                        1,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        comm);
                }
            }

            return result;
        }

        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static public int[] MPIAllGatherv(this int[] send, int[] recvcounts) {
            return send.Int_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static private int[] Int_MPIAllGatherv(this int[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            int[] result = new int[m_recvcounts.Sum()];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (int* pResult = &result[0], pSend = &send[0]) {
                    fixed (int* pRcvcounts = &m_recvcounts[0]) {
                        csMPI.Raw.Allgatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.INT,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.INT,
                            comm);
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Gathers all ulong[] send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static public ulong[] MPIAllGatherv(this ulong[] send, int[] recvcounts) {
            return send.Long_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }
        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static private ulong[] Long_MPIAllGatherv(this ulong[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            ulong[] result = new ulong[m_recvcounts.Sum()];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (ulong* pResult = &result[0], pSend = &send[0]) {
                    fixed (int* pRcvcounts = &m_recvcounts[0]) {
                        csMPI.Raw.Allgatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.UNSIGNED_LONG_LONG,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.UNSIGNED_LONG_LONG,
                            comm);
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// MPI-process with rank 0 gathers this int[] of all MPI-processes in the
        /// <paramref name="comm"/>-communicator with variable length. The length of the gathered int[] is specified by <paramref name="recvcount"/>
        /// </summary>
        /// <param name="recvcount">
        /// Length of the receive buffer
        /// </param>
        static public int[] MPIGatherv(this int[] send, int[] recvcount) {
            return send.MPIGatherv(
                recvcount,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }


        /// <summary>
        /// MPI-process with rank <paramref name="root"/> gathers this int[] of all MPI-processes in the
        /// <paramref name="comm"/>-communicator with variable length. The length of the gathered int[] is specified by <paramref name="recvcount"/>
        /// </summary>
        static public int[] MPIGatherv(this int[] send, int[] recvcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int[] result = rank == root ? new int[recvcounts.Sum()] : null;

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + recvcounts[i - 1];
                }

                fixed (int* pSend = &send[0], pRcvcounts = &recvcounts[0], pResult = result) {
                    Debug.Assert((rank == root) != (pResult == null));

                    csMPI.Raw.Gatherv(
                        (IntPtr)pSend,
                        send.Length,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pResult,
                        (IntPtr)pRcvcounts,
                        (IntPtr)displs,
                        csMPI.Raw._DATATYPE.INT,
                        root,
                        comm);
                }
            }

            return result;
        }
        /// <summary>
        /// MPI-process with rank 0 gathers this ulong[] of all MPI-processes in the
        /// <paramref name="comm"/>-communicator with variable length. The length of the gathered long[] is specified by <paramref name="recvcount"/>
        /// </summary>
        /// <param name="recvcount">
        /// Length of the receive buffer
        /// </param>
        static public ulong[] MPIGatherv(this ulong[] send, int[] recvcount) {
            return send.MPIGatherv(
                recvcount,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }
        /// <summary>
        /// MPI-process with rank <paramref name="root"/> gathers this ulong[] of all MPI-processes in the
        /// <paramref name="comm"/>-communicator with variable length. The length of the gathered long[] is specified by <paramref name="recvcount"/>
        /// </summary>
        /// <param name="recvcount">
        /// Length of the receive buffer
        /// </param>
        static public ulong[] MPIGatherv(this ulong[] send, int[] recvcount, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            ulong[] result = new ulong[recvcount.Sum()];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + recvcount[i - 1];
                }
                //LONG_LONG for long of 64 bits in size
                fixed (ulong* pSend = &send[0], pResult = &result[0]) {
                    fixed (int* pRcvcounts = &recvcount[0]) {
                        csMPI.Raw.Gatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.LONG_LONG,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.LONG_LONG,
                            root,
                            comm);
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Gatherv(IntPtr, int, MPI_Datatype, IntPtr, IntPtr, IntPtr, MPI_Datatype, int, MPI_Comm)"/>
        /// </summary>
        static public double[] MPIGatherv(this double[] send, int[] recvcounts) {
            return send.MPIGatherv(
                recvcounts,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Gatherv(IntPtr, int, MPI_Datatype, IntPtr, IntPtr, IntPtr, MPI_Datatype, int, MPI_Comm)"/>
        /// </summary>
        static public double[] MPIGatherv(this double[] send, int[] recvcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            double[] result = rank == root ? new double[recvcounts.Sum()] : null;

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + recvcounts[i - 1];
                }

                fixed (int*  pRcvcounts = &recvcounts[0]) {
                    fixed (double* pSend = &send[0], pResult = result) {
                        Debug.Assert((rank == root) != (pResult == null));

                        csMPI.Raw.Gatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.DOUBLE,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.DOUBLE,
                            root,
                            comm);
                    }
                }
            }

           
            return result;
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Scatterv(IntPtr, IntPtr, IntPtr, MPI_Datatype, IntPtr, int, MPI_Datatype, int, MPI_Comm)"/>.
        /// </summary>
        static public int[] MPIScatterv(this int[] send, int[] sendcounts) {
            return send.MPIScatterv(
                sendcounts,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Scatterv(IntPtr, IntPtr, IntPtr, MPI_Datatype, IntPtr, int, MPI_Datatype, int, MPI_Comm)"/>.
        /// </summary>
        static public int[] MPIScatterv(this int[] send, int[] sendcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);
            int[] result = new int[Math.Max(1, sendcounts[rank])];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + sendcounts[i - 1];
                }
                if(rank == root) {
                    if(send.Length < displs[size - 1] + sendcounts[size - 1])
                        throw new ArgumentException("Mismatch between send counts and send buffer size.");
                }


                //if (send == null || send.Length == 0) {
                //    // Dummy to avoid null pointer exception
                //    send = new int[1];
                //}

                fixed (int* pSend = send, pSendcounts = &sendcounts[0], pResult = &result[0]) {
                    csMPI.Raw.Scatterv(
                        (IntPtr)pSend,
                        (IntPtr)pSendcounts,
                        (IntPtr)displs,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pResult,
                        sendcounts[rank],
                        csMPI.Raw._DATATYPE.INT,
                        root,
                        comm);
                }
            }

            if (result.Length != sendcounts[rank]) {
                Debug.Assert(result.Length == 1);
                Debug.Assert(sendcounts[rank] == 0);
                Array.Resize(ref result, 0);
            }

            return result;
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Scatterv(IntPtr, IntPtr, IntPtr, MPI_Datatype, IntPtr, int, MPI_Datatype, int, MPI_Comm)"/>.
        /// </summary>
        static public double[] MPIScatterv(this double[] send, int[] sendcounts) {
            return send.MPIScatterv(
                sendcounts,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Scatterv(IntPtr, IntPtr, IntPtr, MPI_Datatype, IntPtr, int, MPI_Datatype, int, MPI_Comm)"/>.
        /// </summary>
        static public double[] MPIScatterv(this double[] send, int[] sendcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);
            double[] result = new double[sendcounts[rank]];

            unsafe
            {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + sendcounts[i - 1];
                    //sum += sendcounts[i];
                }
                if(rank == root) {
                    if(send.Length < displs[size - 1] + sendcounts[size - 1])
                        throw new ArgumentException("Mismatch between send counts and send buffer size.");
                }

                //if (send == null || send.Length == 0) {
                //    // Dummy to avoid null pointer exception
                //    send = new double[1];
                //}

                fixed (int* pSendcounts = &sendcounts[0]) {
                    fixed (double* pSend = send, pResult = &result[0]) {
                        csMPI.Raw.Scatterv(
                            (IntPtr)pSend,
                            (IntPtr)pSendcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.DOUBLE,
                            (IntPtr)pResult,
                            sendcounts[rank],
                            csMPI.Raw._DATATYPE.DOUBLE,
                            root,
                            comm);
                    }
                }
            }

            return result;
        }
    }
}
