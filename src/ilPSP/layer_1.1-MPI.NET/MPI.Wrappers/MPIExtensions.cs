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

using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using NUnit.Framework;
using System;
using System.Collections;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
//using System.Runtime.Serialization.Formatters.Binary;

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

        [Serializable]
        class JsonContainer<T> {
            public T PayLoad;
        }

        static JsonSerializer jsonFormatter {
            get {
                return new JsonSerializer() {
                    NullValueHandling = NullValueHandling.Include,
                    TypeNameHandling = TypeNameHandling.All,
                    ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                    ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
                    //TypeNameAssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Full
                    TypeNameAssemblyFormatHandling = TypeNameAssemblyFormatHandling.Full
                };
            }
        }

        static byte[] SerializeObject<T>(T o) {
            int Size;
            byte[] buffer;
            using(var ms = new MemoryStream()) {
                //_Formatter.Serialize(ms, o);
                using(var w = new BsonDataWriter(ms)) {
                    jsonFormatter.Serialize(w, new JsonContainer<T>() { PayLoad = o });
                    Size = (int)ms.Position;
                }
                buffer = ms.GetBuffer();
            }

            if(Size <= 0)
                throw new IOException("Error serializing object for MPI broadcast - size is 0");
            Array.Resize(ref buffer, Size);
            return buffer;
        }

        static T DeserializeObject<T>(byte[] buffer) {
            using(var ms = new MemoryStream(buffer)) {
                using(var w = new BsonDataReader(ms)) {
                    var containerObj = (JsonContainer<T>)jsonFormatter.Deserialize(w, typeof(JsonContainer<T>));
                    return containerObj.PayLoad;
                }
                
            }
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

            //IFormatter _Formatter = new BinaryFormatter();
            

            // -----------------------------------------------------
            // 1st phase: serialize object and broadcast object size 
            // -----------------------------------------------------

            byte[] buffer = null;
            int Size = -1;

            if(root == MyRank) {
                buffer = SerializeObject<T>(o);
                Size = buffer.Length;
            }

            unsafe {
                csMPI.Raw.Bcast((IntPtr)(&Size), 4, csMPI.Raw._DATATYPE.BYTE, root, comm);
            }

            // ---------------------------
            // 2nd phase: broadcast object
            // ---------------------------
            Debug.Assert(Size > 0);
            if(buffer == null) {
                buffer = new byte[Size];
            }

            unsafe {
                fixed(byte* pBuffer = buffer) {
                    csMPI.Raw.Bcast((IntPtr)pBuffer, Size, csMPI.Raw._DATATYPE.BYTE, root, comm);
                }
            }

            if(MyRank == root) {
                return o;
            } else {
                return DeserializeObject<T>(buffer);
            }
            
        }

        /// <summary>
        /// Gathers (serializeable) objects <paramref name="o"/> from each rank on rank <paramref name="root"/>.
        /// </summary>
        /// <param name="o"></param>
        /// <param name="root">
        /// MPI rank where the objects should be collected
        /// </param>
        /// <returns>
        /// - null on all ranks except <paramref name="root"/>
        /// - on rank <paramref name="root"/>, an array containing the objects from all ranks
        /// </returns>
        static public T[] MPIGatherO<T>(this T o, int root) {
            return o.MPIGatherO(root, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers (serializeable) objects <paramref name="o"/> from each rank on rank <paramref name="root"/>.
        /// </summary>
        /// <param name="o"></param>
        /// <param name="root">
        /// MPI rank where the objects should be collected
        /// </param>
        /// <param name="comm"></param>
        /// <returns>
        /// - null on all ranks except <paramref name="root"/>
        /// - on rank <paramref name="root"/>, an array containing the objects from all ranks
        /// </returns>
        static public T[] MPIGatherO<T>(this T o, int root, MPI_Comm comm) {

            csMPI.Raw.Comm_Rank(comm, out int MyRank);
            csMPI.Raw.Comm_Size(comm, out int MpiSize);

           
            // -----------------------------------------------------
            // 1st phase: serialize object and gather object size 
            // -----------------------------------------------------

            byte[] buffer = null;
            int Size;
            if(root != MyRank) {
                buffer = SerializeObject<T>(o);
                Size = buffer.Length;
            } else {
                buffer = new byte[0];
                Size = 0;
            }

            int[] Sizes = Size.MPIGather(root, comm);

            // -----------------------------------------------------
            // 2nd phase: gather data 
            // -----------------------------------------------------
            byte[] rcvBuffer = buffer.MPIGatherv(Sizes);

            // -----------------------------------------------------
            // 3rd phase: de-serialize
            // -----------------------------------------------------


            if(MyRank == root) {
                T[] ret = new T[MpiSize];

                int i0 = 0;
                for(int r = 0; r < MpiSize; r++) {
                    int sz = Sizes[r];

                    if(r == MyRank) {
                        ret[r] = o;
                    } else {
                        byte[] subBuffer = new byte[sz];
                        Array.Copy(rcvBuffer, i0, subBuffer, 0, sz);
                        ret[r] = DeserializeObject<T>(subBuffer);
                    }
                    i0 += sz;
                }

                return ret;
            } else {
                return null;
            }

        }


        /// <summary>
        /// Gathers (serializeable) objects <paramref name="o"/> from each rank on all processors.
        /// </summary>
        /// <param name="o"></param>
        /// <returns>
        /// an array of all objects on all ranks
        /// </returns>
        static public T[] MPIAllGatherO<T>(this T o) {
            return o.MPIAllGatherO(csMPI.Raw._COMM.WORLD);
        }
        /// <summary>
        /// Gathers (serializeable) objects <paramref name="o"/> from each rank on all processors.
        /// </summary>
        /// <param name="o"></param>
        /// <param name="comm"></param>
        /// <returns>
        /// an array of all objects on all ranks
        /// </returns>
        static public T[] MPIAllGatherO<T>(this T o, MPI_Comm comm) {

            csMPI.Raw.Comm_Rank(comm, out int MyRank);
            csMPI.Raw.Comm_Size(comm, out int MpiSize);

            //IFormatter _Formatter = new BinaryFormatter();

            // -----------------------------------------------------
            // 1st phase: serialize object and gather object size 
            // -----------------------------------------------------

            byte[] buffer = null;
            int Size;

            using(var ms = new MemoryStream()) {
                buffer = SerializeObject<T>(o);
                Size = buffer.Length;
            }

            int[] Sizes = Size.MPIAllGather(comm);

            // -----------------------------------------------------
            // 2nd phase: gather data 
            // -----------------------------------------------------
            byte[] rcvBuffer = buffer.MPIAllGatherv(Sizes);

            // -----------------------------------------------------
            // 3rd phase: de-serialize
            // -----------------------------------------------------
            T[] ret = new T[MpiSize];
            using(var ms = new MemoryStream(rcvBuffer)) {
                int i0 = 0;
                for(int r = 0; r < MpiSize; r++) {
                    int sz = Sizes[r];

                    byte[] subBuffer = new byte[sz];
                    Array.Copy(rcvBuffer, i0, subBuffer, 0, sz);
                    ret[r] = DeserializeObject<T>(subBuffer);

                    i0 += sz;
                }
            }
            return ret;
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

            if(e != null)
                csMPI.Raw.Comm_Rank(comm, out ExcSrc);

            unsafe {
                int res = 0;
                csMPI.Raw.Allreduce((IntPtr)(&ExcSrc), (IntPtr)(&res), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, comm);
                ExcSrc = res;
            }


            if(ExcSrc < int.MaxValue) {
                // an exception occured on some process

                int myRank;
                csMPI.Raw.Comm_Rank(comm, out myRank);

                object reduced;

                if(myRank == ExcSrc) {
                    // sender branch


                    if(e.GetType().GetCustomAttributes(typeof(SerializableAttribute), false).Length > 0)
                        // exception is serializeable -> bcast exception itself
                        reduced = MPIBroadcast<object>(e, ExcSrc, comm);
                    else
                        // bcast exception message
                        reduced = MPIBroadcast<object>(e.GetType().Name + ": '" + e.Message + "'", ExcSrc, comm);
                } else {
                    // receiver branch

                    reduced = MPIBroadcast<object>(null, ExcSrc, comm);
                }

                if(reduced is string) {
                    throw new ApplicationException("On MPI Process #" + myRank + ": " + ((string)reduced));
                } else if(reduced is Exception) {
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
        /// equal to <see cref="MPISum(long,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public long MPISum(this long i) {
            return MPISum(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the sum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public long MPISum(this long i, MPI_Comm comm) {
            long loc = i;
            unsafe {
                long glob = long.MinValue;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.LONG_LONG_INT, csMPI.Raw._OP.SUM, comm);
                return glob;
            }
        }


        /// <summary>
        /// equal to <see cref="MPIOr(bool,MPI_Comm)"/>, acting on the WORLD-communicator
        /// </summary>
        static public bool MPIOr(this bool i) {
            return MPIOr(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the logical or of <paramref name="b"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public bool MPIOr(this bool b, MPI_Comm comm) {
            int loc = b ? 1: 0;
            unsafe {
                int glob = 0;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
                return glob > 0;
            }
        }

        /// <summary>
        /// computes the logical or of <paramref name="b"/> on all MPI-processes in the
        /// WORLD--communicator and stores the result in-place
        /// </summary>
        static public void MPIOr(this BitArray b) {
            MPIOr(b, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// computes the logical or of <paramref name="b"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator and stores the result in-place
        /// </summary>
        static public void MPIOr(this BitArray b, MPI_Comm comm) {
            int L = b.Length;
            int Lx = L / 32 + 1;
            int[] s_buf = new int[Lx];
            int[] r_buf = new int[Lx];

            int k = 0;
            for(int lx = 0; lx < Lx; lx++) {
                int a = 0;
                for(int i = 0; i < 32; i++) {
                    bool b_k = b[k];
                    if(b_k) {
                        a |= (1 << i);
                    }
                    k++;
                    if(k >= L)
                        break;
                }
                s_buf[lx] = a;

                if(k >= L)
                    break;
            }


            unsafe {
                fixed(int* p_s_buf = s_buf, p_r_buf = r_buf) {
                    csMPI.Raw.Allreduce(((IntPtr)(p_s_buf)), ((IntPtr)(p_r_buf)), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.BOR, comm);
                }
            }

            k = 0;
            for(int lx = 0; lx < Lx; lx++) {
                int a = r_buf[lx];
                for(int i = 0; i < 32; i++) {
                    bool b_k = (a  & (1 << i)) != 0;
                    b[k] = b_k;
                    k++;
                    if(k >= L)
                        break;
                }

                if(k >= L)
                    break;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIAnd(System.Collections.Generic.IEnumerable{bool},MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public bool[] MPIAnd(this System.Collections.Generic.IEnumerable<bool> i) {
            return MPIAnd(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the logical and of <paramref name="bb"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public bool[] MPIAnd(this System.Collections.Generic.IEnumerable<bool> bb, MPI_Comm comm) {
            int[] loc = bb.Select(b => b ? 1 : 0).ToArray();
            int[] glb = new int[loc.Length];
            unsafe {
                fixed (int* ploc = loc, pglb = glb) {
                    csMPI.Raw.Allreduce(((IntPtr)ploc), ((IntPtr)pglb), loc.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.PROD, comm);
                }
            }
            return glb.Select(g => g != 0).ToArray();
        }

        /// <summary>
        /// equal to <see cref="MPIOr(System.Collections.Generic.IEnumerable{bool},MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public bool[] MPIOr(this System.Collections.Generic.IEnumerable<bool> i) {
            return MPIOr(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the logical or of <paramref name="bb"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public bool[] MPIOr(this System.Collections.Generic.IEnumerable<bool> bb, MPI_Comm comm) {
            int[] loc = bb.Select(b => b ? 1 : 0).ToArray();
            int[] glb = new int[loc.Length];
            unsafe {
                fixed (int* ploc = loc, pglb = glb) {
                    csMPI.Raw.Allreduce(((IntPtr)ploc), ((IntPtr)pglb), loc.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
                }
            }
            return glb.Select(g => g != 0).ToArray();
        }


        /// <summary>
        /// equal to <see cref="MPIAnd(int,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public bool MPIAnd(this bool i) {
            return MPIAnd(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the logical and of <paramref name="b"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public bool MPIAnd(this bool b, MPI_Comm comm) {
            int loc = b ? 1 : 0;
            unsafe {
                int glob = 0;
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.PROD, comm);
                return glob > 0;
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
                fixed(double* loc = i, glob = R) {
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
                fixed(double* loc = i, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), i.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, comm);
                }
            }
            return R;
        }




        /// <summary>
        /// equal to <see cref="MPIEqual(bool,MPI_Comm)"/> acting on WORLD-communicator
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        static public bool MPIEquals(this bool b) {
            return MPIEquals(b, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// 
        /// </summary>
        static public bool MPIEquals(this bool b, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int sz);
            if(sz <= 1)
                return true;

            int loc = (b ? 1 : 0);
            int glob = 0;
            unsafe {
                csMPI.Raw.Allreduce(((IntPtr)(&loc)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
            }
            return glob == loc*sz ? true : false;
        }

        /// <summary>
        /// equal to <see cref="MPIEquals(double,MPI_Comm)"/> acting on WORLD-communicator
        /// </summary>
        static public bool MPIEquals(this double i) {
            return MPIEquals(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns true on every process, if <paramref name="i"/> is equal at every process
        /// </summary>
        /// <param name="i"></param>
        /// <param name="comm"></param>
        /// <returns></returns>
        static public bool MPIEquals(this double i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int sz);
            if(sz <= 1)
                return true;

            unsafe {

                //
                // compute bit-wise equality of (a1, a2, ... , aN) as
                //   ( a1 & a2 & ... & aN ) | ( !a1 & !a2 & ... & !aN )
                // (the sum-aproach which was here before can have round-of errors);
                //

                ulong* pInp = (ulong*)(&i);

                ulong* inp2 = stackalloc ulong[2];
                inp2[0] = *pInp;    
                inp2[1] = ~(*pInp);
                
                ulong* glob1 = stackalloc ulong[2];

                csMPI.Raw.Allreduce(((IntPtr)(inp2)), ((IntPtr)(glob1)), 2, csMPI.Raw._DATATYPE.UNSIGNED_LONG_LONG, csMPI.Raw._OP.BAND, comm);

                ulong bitWiseEq = glob1[0] | glob1[1]; // now, all bits must be true;


                return bitWiseEq == ulong.MaxValue;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIEquals(int,MPI_Comm)"/> acting on WORLD-communicator
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        static public bool MPIEquals(this int i) {
            return MPIEquals(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns true on every process, if <paramref name="i"/> is equal at every process
        /// </summary>
        /// <param name="i"></param>
        /// <param name="comm"></param>
        /// <returns></returns>
        static public bool MPIEquals(this int i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int sz);
            if(sz <= 1)
                return true;

            int glob = int.MaxValue;
            unsafe {
                csMPI.Raw.Allreduce(((IntPtr)(&i)), ((IntPtr)(&glob)), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
            }
            return glob == i*sz ? true : false;
        }

        /// <summary>
        /// equal to <see cref="MPIEquals(double[],MPI_Comm)"/> acting on WORLD-communicator
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        static public bool[] MPIEquals(this double[] i) {
            return MPIEquals(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns true on every process for every entry of <paramref name="i"/> if they are equal at every process
        /// </summary>
        /// <param name="bAry"></param>
        /// <param name="comm"></param>
        /// <returns></returns>
        static public bool[] MPIEquals(this double[] bAry, MPI_Comm comm) {
            bool[] check = new bool[bAry.Length];

            csMPI.Raw.Comm_Size(comm, out int sz);
            if(sz <= 1) {
                for(int i = 0; i < check.Length; i++)
                    check[i] = true;
            }

            double[] Rmin = new double[bAry.Length];
            double[] Rmax = new double[bAry.Length];
            unsafe {
                fixed (void* loc = bAry, globMin = Rmin, globMax = Rmax) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(globMin)), bAry.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, comm);
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(globMax)), bAry.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, comm);
                }
            }
            for (int k = 0; k < bAry.Length; k++) {
                check[k] = (Rmin[k] == bAry[k]) && (Rmax[k] == bAry[k]);
            }
            return check;
        }

        /// <summary>
        /// equal to <see cref="MPIEquals(int[],MPI_Comm)"/> acting on WORLD-communicator
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        static public bool[] MPIEquals(this int[] i)
        {
            return MPIEquals(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns true on every process for every entry of <paramref name="iAry"/> if they are equal at every process
        /// </summary>
        /// <param name="iAry"></param>
        /// <param name="comm"></param>
        /// <returns></returns>
        static public bool[] MPIEquals(this int[] iAry, MPI_Comm comm) {
            bool[] check = new bool[iAry.Length];

            csMPI.Raw.Comm_Size(comm, out int sz);
            if(sz <= 1) {
                for(int i = 0; i < check.Length; i++)
                    check[i] = true;
            }

            int[] R = new int[iAry.Length];
            unsafe {
                fixed(int* loc = iAry, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), iAry.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, comm);
                }
            }
            for(int k = 0; k < iAry.Length; k++) {
                check[k] = R[k] == iAry[k]*sz ? true : false;
            }
            return check;
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
        /*
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
        static public int[] MPIAnd(this int[] i, MPI_Comm comm) {
            int[] R = new int[i.Length];
            unsafe {
                fixed (int* loc = i, glob = R) {
                    csMPI.Raw.Allreduce(((IntPtr)(loc)), ((IntPtr)(glob)), i.Length, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.A, comm);
                }
            }
            return R;
        }
        */

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
                    csMPI.Raw.Allreduce(((IntPtr)(pA)), ((IntPtr)(pS)), A.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, comm);
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

        /// <summary>
        /// equal to <see cref="MPIMax(int,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public long MPIMax(this long i) {
            return MPIMax(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the maximum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public long MPIMax(this long i, MPI_Comm comm) {
            long loc = i;
            unsafe {
                long glob = long.MinValue;
                csMPI.Raw.Allreduce(
                    (IntPtr)(&loc),
                    (IntPtr)(&glob),
                    1,
                    csMPI.Raw._DATATYPE.LONG_LONG,
                    csMPI.Raw._OP.MAX,
                    comm);
                return glob;
            }
        }

        /// <summary>
        /// equal to <see cref="MPIMin(int,MPI_Comm)"/>, acting on the
        /// WORLD-communicator
        /// </summary>
        static public long MPIMin(this long i) {
            return MPIMin(i, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// returns the maximum of <paramref name="i"/> on all MPI-processes in the
        /// <paramref name="comm"/>--communicator.
        /// </summary>
        static public long MPIMin(this long i, MPI_Comm comm) {
            long loc = i;
            unsafe {
                long glob = long.MaxValue;
                csMPI.Raw.Allreduce(
                    (IntPtr)(&loc),
                    (IntPtr)(&glob),
                    1,
                    csMPI.Raw._DATATYPE.LONG_LONG,
                    csMPI.Raw._OP.MIN,
                    comm);
                return glob;
            }
        }


        /// <summary>
        /// Gathers single numbers form each MPI rank in an array
        /// </summary>
        static public int[] MPIAllGather(this int i) {
            return i.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers single numbers form each MPI rank in an array
        /// </summary>
        static public int[] MPIAllGather(this int i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);

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

        /// <summary>
        /// Gathers single numbers form each MPI rank in an array at the <paramref name="root"/> rank
        /// </summary>
        static public int[] MPIGather(this int i, int root) {
            return i.MPIGather(root, csMPI.Raw._COMM.WORLD);
        }
        
        /// <summary>
        /// Gathers single numbers form each MPI rank in an array at the <paramref name="root"/> rank
        /// </summary>
        static public int[] MPIGather(this int i, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int[] result;
            if (rank == root)
                result = new int[size];
            else
                result = null;

            unsafe {
                int sendBuffer = i;
                fixed (int* pResult = result) {
                    csMPI.Raw.Gather(
                        (IntPtr)(&i),
                        1,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pResult,
                        1,
                        csMPI.Raw._DATATYPE.INT,
                        root,
                        comm);
                }
            }

            return result;
        }

        /// <summary>
        /// Scatters form <paramref name="root"/> rank to all other ranks
        /// </summary>
        static public int[] MPIScatter(this int[] L, int recvCount, int root) {
            return L.MPIScatter(recvCount, root, csMPI.Raw._COMM.WORLD);
        }
        
        /// <summary>
        /// Scatters form <paramref name="root"/> rank to all other ranks
        /// </summary>
        static public int[] MPIScatter(this int[] L, int recvCount, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            if(recvCount < 0)
                throw new ArgumentOutOfRangeException("Receive size must be greater than 0.");


            int SndCount;
            if(rank == root) {
                SndCount = recvCount;
                if(L.Length != recvCount * size)
                    throw new ArgumentOutOfRangeException("Send buffer length must be == receive count * number of processors.");
            } else {
                SndCount = 0;
            }

            int[] result = new int[recvCount];

            if(recvCount == 0)
                return result;

            unsafe {
                fixed(int* pL = L, pResult = result) {
                    csMPI.Raw.Scatter(
                        (IntPtr)(pL),
                        SndCount,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pResult,
                        recvCount,
                        csMPI.Raw._DATATYPE.INT,
                        root,
                        comm);

                }
            }

            return result;
        }



        /// <summary>
        /// Gathers boolean arrays form each MPI rank in an array
        /// </summary>
        static public bool[] MPIAllGather(this bool[] i) {
            return i.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers boolean arrays form each MPI rank in an array
        /// </summary>
        static public bool[] MPIAllGather(this bool[] i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int L = i.Length;
            bool[] result = new bool[L * size];


            unsafe {
                fixed (void* sendbuf = i, recvbuf = result) {
                    
                    csMPI.Raw.Allgather(
                        (IntPtr)sendbuf,
                        1,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)recvbuf,
                        1,
                        csMPI.Raw._DATATYPE.INT,
                        comm);
                }
            }
            return result;
        }

        /// <summary>
        /// Gathers integer arrays form each MPI rank in an array
        /// </summary>
        static public int[] MPIAllGather(this int[] i) {
            return i.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers integer arrays form each MPI rank in an array
        /// </summary>
        static public int[] MPIAllGather(this int[] i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int L = i.Length;
            int[] result = new int[size*L];


            unsafe {
                fixed (void* sendbuf = i, recvbuf = result) {
                    
                    csMPI.Raw.Allgather(
                        (IntPtr)sendbuf,
                        L,
                        csMPI.Raw._DATATYPE.INT,
                        (IntPtr)recvbuf,
                        L,
                        csMPI.Raw._DATATYPE.INT,
                        comm);
                }
            }
            return result;
        }


        /// <summary>
        /// Gathers double arrays form each MPI rank in an array
        /// </summary>
        static public double[] MPIAllGather(this double[] i) {
            return i.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers double arrays form each MPI rank in an array
        /// </summary>
        static public double[] MPIAllGather(this double[] i, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            //csMPI.Raw.Comm_Rank(comm, out int rank);

            int L = i.Length;
            double[] result = new double[size*L];


            unsafe {
                fixed (void* sendbuf = i, recvbuf = result) {

                    csMPI.Raw.Allgather(
                        (IntPtr)sendbuf,
                        L,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        (IntPtr)recvbuf,
                        L,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        comm);
                }
            }
            return result;
        }

        /// <summary>
        /// Gathers single numbers form each MPI rank in an array
        /// </summary>
        static public double[] MPIAllGather(this double d) {
            return d.MPIAllGather(csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers single numbers form each MPI rank in an array
        /// </summary>
        static public double[] MPIAllGather(this double d, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);

            double[] result = new double[size];
            unsafe {
                double sendBuffer = d;
                fixed (double* pResult = result) {
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
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public int[] MPIAllGatherv(this int[] send, MPI_Comm comm) {
            int sz = send.Length;
            int[] AllSz = sz.MPIAllGather(comm);
            return MPIAllGatherv(send, AllSz, comm);
        }

        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public int[] MPIAllGatherv(this int[] send) {
            return MPIAllGatherv(send, csMPI.Raw._COMM.WORLD);
        }


        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public double[] MPIAllGatherv(this double[] send, MPI_Comm comm) {
            int sz = send.Length;
            int[] AllSz = sz.MPIAllGather(comm);
            return MPIAllGatherv(send, AllSz, comm);
        }

        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public double[] MPIAllGatherv(this double[] send) {
            return MPIAllGatherv(send, csMPI.Raw._COMM.WORLD);
        }


        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public int[] MPIAllGatherv(this int[] send, int[] recvcounts) {
            return send.Int_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public int[] MPIAllGatherv(this int[] send, int[] recvcounts, MPI_Comm comm) {
            return send.Int_MPIAllGatherv(recvcounts, comm);
        }

        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static private int[] Int_MPIAllGatherv(this int[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            int rcs = m_recvcounts.Sum();
            if (rcs == 0)
                return new int[0];


            int[] result = new int[rcs];
            if (send.Length == 0)
                send = new int[1];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (int* pResult = result, pSend = send) {
                    fixed (int* pRcvcounts = m_recvcounts) {
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
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public double[] MPIAllGatherv(this double[] send, int[] recvcounts) {
            return send.Double_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers all int[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public double[] MPIAllGatherv(this double[] send, int[] recvcounts, MPI_Comm comm) {
            return send.Double_MPIAllGatherv(recvcounts, comm);
        }

        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static private double[] Double_MPIAllGatherv(this double[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            int rcs = m_recvcounts.Sum();
            if (rcs == 0)
                return new double[0];


            double[] result = new double[rcs];
            if (send.Length == 0)
                send = new double[1];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (double* pResult = result, pSend = send) {
                    fixed (int* pRcvcounts = m_recvcounts) {
                        csMPI.Raw.Allgatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.DOUBLE,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.DOUBLE,
                            comm);
                    }
                }
            }

            return result;
        }


        /// <summary>
        /// Gathers all byte[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public byte[] MPIAllGatherv(this byte[] send, int[] recvcounts) {
            return send.Byte_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every jth block of data is from the jth process.
        /// </summary>
        static private byte[] Byte_MPIAllGatherv(this byte[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            int rcs = m_recvcounts.Sum();
            if (rcs == 0)
                return new byte[0];


            byte[] result = new byte[rcs];
            if (send.Length == 0)
                send = new byte[1];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (byte* pResult = result, pSend = send) {
                    fixed (int* pRcvcounts = m_recvcounts) {
                        csMPI.Raw.Allgatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.BYTE,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.BYTE,
                            comm);
                    }
                }
            }

            return result;
        }
        

        /// <summary>
        /// Gathers all ulong[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public ulong[] MPIAllGatherv(this ulong[] send, int[] recvcounts) {
            return send.Long_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static private ulong[] Long_MPIAllGatherv(this ulong[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            int rcs = m_recvcounts.Sum();
            if (rcs == 0)
                return new ulong[0];

            ulong[] result = new ulong[rcs];

            if (send.Length == 0)
                send = new ulong[1];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (ulong* pResult = result, pSend = send) {
                    fixed (int* pRcvcounts = m_recvcounts) {
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
        /// Gathers all long[] send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static public long[] MPIAllGatherv(this long[] send, int[] recvcounts) {
            return send.Long_MPIAllGatherv(recvcounts, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathers all send Arrays on all MPI-processes, at which every j-th block of data is from the j-th process.
        /// </summary>
        static private long[] Long_MPIAllGatherv(this long[] send, int[] m_recvcounts, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            int rcs = m_recvcounts.Sum();
            if (rcs == 0)
                return new long[0];

            long[] result = new long[rcs];

            if (send.Length == 0)
                send = new long[1];

            unsafe {
                int* displs = stackalloc int[size];
                for (int i = 1; i < size; i++) {
                    displs[i] = displs[i - 1] + m_recvcounts[i - 1];
                }
                fixed (long* pResult = result, pSend = send) {
                    fixed (int* pRcvcounts = m_recvcounts) {
                        csMPI.Raw.Allgatherv(
                            (IntPtr)pSend,
                            send.Length,
                            csMPI.Raw._DATATYPE.LONG_LONG,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.LONG_LONG,
                            comm);
                    }
                }
            }
            
            return result;
        }


        /// <summary>
        /// Gathering of a jagged array on all processors
        /// </summary>
        static public long[][] MPI_AllGaterv(this long[][] send) {
            return MPI_AllGaterv(send, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Gathering of a jagged array on all processors
        /// </summary>
        static public long[][] MPI_AllGaterv(this long[][] send, MPI_Comm comm) {
            int Jloc = send.Length;
            int NoOfItemsToSend = 0;
            for(int j = 0; j < Jloc; j++) {
                if(send[j] != null)
                    NoOfItemsToSend += send[j].Length;
            }
            NoOfItemsToSend += Jloc;


            // compress the jagged array into a linear one which is efficient to send.
            long[] SendBuffer = new long[NoOfItemsToSend];
            int k = 0;
            for(int j = 0; j < Jloc; j++) {
                var Nj = send[j];
                int Lj = Nj != null ? Nj.Length : 0;
                SendBuffer[k] = -Lj - 1; k++; // the negative minus one encoding is not necessary, but it will raise an exception if the sequence gets messed up
                for(int i = 0; i < Lj; i++) {
                    SendBuffer[k + i] = Nj[i];
                }
                k += Lj;
            }
            Assert.AreEqual(k, NoOfItemsToSend);

            int[] NoOfItemsPerProc = NoOfItemsToSend.MPIAllGather(comm);
            int Jglob = Jloc.MPISum(comm);
            long[] RcvBuffer = SendBuffer.Long_MPIAllGatherv(NoOfItemsPerProc, comm);

            // unpack data to jagged array
            long[][] globalCellNeigbourship = new long[Jglob][];
            k = 0;
            long jG = 0;
            while(k < RcvBuffer.Length) {
                int L_jG = checked((int)(-(RcvBuffer[k] + 1))); k++;
                var gN_jG = new long[L_jG];
                globalCellNeigbourship[jG] = gN_jG;
                for(int i = 0; i < L_jG; i++) {
                    gN_jG[i] = RcvBuffer[k + i];
                }
                k += L_jG;
                jG++;
            }
            Assert.AreEqual(jG, Jglob);

            // finally;
            return globalCellNeigbourship;
        }


        /// <summary>
        /// MPI-process with rank 0 gathers this int[] over all MPI-processes in the world-communicator.
        /// The length of the gathered int[] is specified by <paramref name="recvcount"/>
        /// </summary>
        /// <param name="recvcount">
        /// number of items to receive from each rank
        /// </param>
        /// <param name="send">
        /// data to send from current process to rank 0
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
        /// <param name="recvcounts">
        /// Significant only at <paramref name="root"/> process. 
        /// number of items to receive from each process in the communicator <paramref name="comm"/>;
        /// </param>
        /// <param name="root">
        /// rank of receiving process
        /// </param>
        /// <param name="comm">
        /// communicator
        /// </param>
        /// <param name="send">
        /// data to send at each process.
        /// </param>
        static public int[] MPIGatherv(this int[] send, int[] recvcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int rcs = rank == root ? recvcounts.Sum() : 0;
            int[] result = rank == root ? new int[Math.Max(1, rcs)] : null;
            

            unsafe {
                int* displs = stackalloc int[size];
                if (rank == root) {
                    for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcounts[i - 1];
                    }
                }

                fixed (int* pSend = send, pRcvcounts = recvcounts, pResult = result) {
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

            if(result != null && result.Length > rcs) {
                Debug.Assert(rcs == 0);
                result = new int[0];
            }

            return result;
        }
        /// <summary>
        /// MPI-process with rank 0 gathers this ulong[] of all MPI-processes in the
        /// </summary>
        /// <param name="recvcount">
        /// number of items to receive from each sender
        /// </param>
        /// <param name="send">
        /// data to send
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
        /// number of items to receive from each sender
        /// </param>
        /// <param name="send">
        /// data to send
        /// </param>
        /// <param name="comm"></param>
        /// <param name="root">rank of receiver process</param>
        static public ulong[] MPIGatherv(this ulong[] send, int[] recvcount, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int rcs = rank == root ? recvcount.Sum() : 0;
            ulong[] result = rank == root ? new ulong[Math.Max(1, rcs)] : null;

            unsafe {
                int* displs = stackalloc int[size];
                if(rank == root)
                    for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcount[i - 1];
                    }
                //LONG_LONG for long of 64 bits in size
                fixed (ulong* pSend = send, pResult = result) {
                    fixed (int* pRcvcounts = recvcount) {
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

            if (result != null && result.Length > rcs) {
                Debug.Assert(rcs == 0);
                result = new ulong[0];
            }

            return result;
        }

        /// <summary>
        /// MPI-process with rank 0 gathers this ulong[] of all MPI-processes in the
        /// </summary>
        /// <param name="recvcount">
        /// number of items to receive from each sender
        /// </param>
        /// <param name="send">
        /// data to send
        /// </param>
        static public long[] MPIGatherv(this long[] send, int[] recvcount) {
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
        /// number of items to receive from each sender
        /// </param>
        /// <param name="send">
        /// data to send
        /// </param>
        /// <param name="comm"></param>
        /// <param name="root">rank of receiver process</param>
        static public long[] MPIGatherv(this long[] send, int[] recvcount, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int rcs = rank == root ? recvcount.Sum() : 0;
            long[] result = rank == root ? new long[Math.Max(1, rcs)] : null;

            unsafe {
                int* displs = stackalloc int[size];
                if(rank == root)
                    for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcount[i - 1];
                    }
                //LONG_LONG for long of 64 bits in size
                fixed (long* pSend = send, pResult = result) {
                    fixed (int* pRcvcounts = recvcount) {
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

            if (result != null && result.Length > rcs) {
                Debug.Assert(rcs == 0);
                result = new long[0];
            }

            return result;
        }


        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Gatherv"/>
        /// </summary>
        static public double[] MPIGatherv(this double[] send, int[] recvcounts) {
            return send.MPIGatherv(
                recvcounts,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Gatherv"/>
        /// </summary>
        static public double[] MPIGatherv(this double[] send, int[] recvcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int rcs = rank == root ? recvcounts.Sum() : 0;
            double[] result = rank == root ? new double[Math.Max(1, rcs)] : null;


            unsafe {
                int* displs = stackalloc int[size];
                if(rank == root)
                    for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcounts[i - 1];
                    }

                fixed (int*  pRcvcounts = recvcounts) {
                    fixed (double* pSend = send, pResult = result) {
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

            if (result != null && result.Length > rcs) {
                Debug.Assert(rcs == 0);
                result = new double[0];
            }

           

            return result;
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Gatherv"/>
        /// </summary>
        static public byte[] MPIGatherv(this byte[] send, int[] recvcounts) {
            return send.MPIGatherv(
                recvcounts,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Gatherv"/>
        /// </summary>
        static public byte[] MPIGatherv(this byte[] send, int[] recvcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);

            int outsize;
            byte[] result;
            if(rank == root) {
                outsize = recvcounts.Sum();
                result = new byte[Math.Max(outsize, 1)]; // this is only because a 0-length array maps to unsafe null
            } else {
                result = null;
                outsize = 0;
            }
            
            unsafe {
                int* displs = stackalloc int[size];
                if(rank == root)
                    for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + recvcounts[i - 1];
                    }

                fixed (int* pRcvcounts = recvcounts) {
                    int lsend = send.Length;
                    if (lsend <= 0)
                        send = new byte[0];
                    fixed (byte* pSend = send, pResult = result) {
                        Debug.Assert((rank == root) != (pResult == null));

                        csMPI.Raw.Gatherv(
                            (IntPtr)pSend,
                            lsend,
                            csMPI.Raw._DATATYPE.BYTE,
                            (IntPtr)pResult,
                            (IntPtr)pRcvcounts,
                            (IntPtr)displs,
                            csMPI.Raw._DATATYPE.BYTE,
                            root,
                            comm);
                    }
                }
            }

            if (outsize > 0) {
                Debug.Assert(outsize == result.Length);
                return result;
            } else {
                return new byte[0];
            }
        }


        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Scatterv"/>.
        /// </summary>
        static public int[] MPIScatterv(this int[] send, int[] sendcounts) {
            return send.MPIScatterv(
                sendcounts,
                root: 0,
                comm: csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// Wrapper around <see cref="IMPIdriver.Scatterv"/>.
        /// </summary>
        static public int[] MPIScatterv(this int[] send, int[] sendcounts, int root, MPI_Comm comm) {
            csMPI.Raw.Comm_Size(comm, out int size);
            csMPI.Raw.Comm_Rank(comm, out int rank);
            int[] result = new int[Math.Max(1, sendcounts[rank])];

            unsafe {
                int* displs = stackalloc int[size];
                if (rank == root) {
                    for (int i = 1; i < size; i++) {
                        displs[i] = displs[i - 1] + sendcounts[i - 1];
                    }
                    if (send.Length < displs[size - 1] + sendcounts[size - 1])
                        throw new ArgumentException("Mismatch between send counts and send buffer size.");
                }


                //if (send == null || send.Length == 0) {
                //    // Dummy to avoid null pointer exception
                //    send = new int[1];
                //}

                fixed (int* pSend = send, pSendcounts = sendcounts, pResult = result) {
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

                fixed (int* pSendcounts = sendcounts) {
                    fixed (double* pSend = send, pResult = result) {
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
