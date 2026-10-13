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
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;
using System.Threading;

namespace ilPSP {


    /// <summary>
    /// Manages buffers for temporary use;
    /// use with care -- never return the buffer objects!
    /// </summary>
    /// <remarks>
    /// Two levels:
    /// - a few thread-local slots (<see cref="NO_OF_THREADLOCAL_SLOTS"/>) per thread, which are used first;
    ///   they require no lock, and the buffers are kept (strong reference) for re-use by the same thread;
    /// - a shared pool (weak references, protected by a lock) as a fall-back, if all thread-local slots
    ///   are in use, or if a buffer is larger than <see cref="MAX_THREADLOCAL_SIZE"/>.
    ///
    /// Handles returned by <see cref="GetTempBuffer(out int, int)"/> are always positive:
    /// thread-local slots are encoded as <see cref="THREADLOCAL_HANDLE_OFFSET"/> + slot index,
    /// shared-pool buffers as their pool index (starting at 1, since some callers use 0 as a "no buffer" sentinel).
    /// A buffer must be released by <see cref="FreeTempBuffer(int)"/> **on the same thread** which obtained it.
    /// </remarks>
    public static class TempBuffer {

        // ===================================================================================
        // thread-local slots
        // ===================================================================================

        /// <summary>
        /// Number of thread-local buffers per thread.
        /// </summary>
        public const int NO_OF_THREADLOCAL_SLOTS = 3;

        /// <summary>
        /// Requests larger than this (in number of doubles; 2^20 doubles = 8 MB) always go to the shared pool,
        /// so that idle threads do not keep huge buffers alive.
        /// </summary>
        public const int MAX_THREADLOCAL_SIZE = 1 << 20;

        /// <summary>
        /// Handles larger or equal to this value refer to thread-local slots.
        /// </summary>
        const int THREADLOCAL_HANDLE_OFFSET = 1 << 30;

        // no field initializers for [ThreadStatic] fields: they would only run for the first thread.
        [ThreadStatic] static double[][] t_Buffers;
        [ThreadStatic] static bool[] t_InUse;
        [ThreadStatic] static int t_NoOfBuffersInUse;
        [ThreadStatic] static int t_MaxNoOfBuffersInUse;

        static int s_MaxNoOfBuffersInUse;

        /// <summary>
        /// Diagnostic: maximum number of buffers which were held simultaneously by one thread (over all threads, since program start).
        /// </summary>
        public static int MaxSimultaneousBuffersPerThread {
            get {
                return Volatile.Read(ref s_MaxNoOfBuffersInUse);
            }
        }

        static double[] GetThreadLocal(out int handle, int size) {
            var buffers = t_Buffers;
            var inUse = t_InUse;
            if(buffers == null) {
                buffers = new double[NO_OF_THREADLOCAL_SLOTS][];
                inUse = new bool[NO_OF_THREADLOCAL_SLOTS];
                t_Buffers = buffers;
                t_InUse = inUse;
            }

            // best fit among free slots which are large enough; otherwise, re-allocate the largest free slot
            int iFit = -1, iGrow = -1;
            for(int i = 0; i < NO_OF_THREADLOCAL_SLOTS; i++) {
                if(inUse[i])
                    continue;
                int len = buffers[i] != null ? buffers[i].Length : 0;
                if(len >= size) {
                    if(iFit < 0 || len < buffers[iFit].Length)
                        iFit = i;
                } else {
                    if(iGrow < 0 || len > (buffers[iGrow] != null ? buffers[iGrow].Length : 0))
                        iGrow = i;
                }
            }

            int iSlot;
            if(iFit >= 0) {
                iSlot = iFit;
            } else if(iGrow >= 0) {
                iSlot = iGrow;
                int oldLen = buffers[iSlot] != null ? buffers[iSlot].Length : 0;
                // some head-room, to avoid frequent re-allocation for slowly growing requests
                int newLen = Math.Max(size, Math.Min(MAX_THREADLOCAL_SIZE, oldLen + oldLen / 2));
                buffers[iSlot] = new double[newLen];
            } else {
                // all thread-local slots are in use
                handle = -1;
                return null;
            }

            inUse[iSlot] = true;
            handle = THREADLOCAL_HANDLE_OFFSET + iSlot;
            return buffers[iSlot];
        }

        static void FreeThreadLocal(int handle) {
            int iSlot = handle - THREADLOCAL_HANDLE_OFFSET;
            var inUse = t_InUse;
            if(inUse == null || iSlot < 0 || iSlot >= NO_OF_THREADLOCAL_SLOTS || !inUse[iSlot])
                // either a double free, or (more severe) the buffer is freed on a different thread than it was obtained on
                throw new ApplicationException($"TempBuffer: invalid free of thread-local handle {handle} (double free, or free on a different thread).");
            inUse[iSlot] = false;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void CountGet() {
            int n = ++t_NoOfBuffersInUse;
            if(n > t_MaxNoOfBuffersInUse) {
                t_MaxNoOfBuffersInUse = n;
                // rare: only if the per-thread maximum grows
                int glob;
                while((glob = Volatile.Read(ref s_MaxNoOfBuffersInUse)) < n) {
                    if(Interlocked.CompareExchange(ref s_MaxNoOfBuffersInUse, n, glob) == glob)
                        break;
                }
            }
        }

        // ===================================================================================
        // shared pool (fall-back)
        // ===================================================================================

        static readonly object s_PoolLock = new object();

        /// <summary>
        /// Initial pool size; the pool grows if necessary (<see cref="GetPooled"/>).
        /// </summary>
        /// <remarks>
        /// Note: <see cref="ilPSP.Environment.NumThreads"/> may not yet be initialized when this class is initialized,
        /// therefore the pool must be able to grow.
        /// </remarks>
        static int INITIAL_POOL_SIZE {
            get {
                return 8 * 8 * Math.Max(1, ilPSP.Environment.NumThreads);
            }
        }

        static WeakReference[] BigBuf = new WeakReference[INITIAL_POOL_SIZE];
        static bool[] BufLok = new bool[INITIAL_POOL_SIZE];

        static double[] GetPooled(out int _iBuf, int size) {
            lock(s_PoolLock) {
                double[] BB = null;
                int iAlloc = -1;
                int iBuf = -1;

                // try to find a buffer which already fits...
                // (index 0 is never used: some callers use 0 as "no buffer" sentinel)
                for(int iB = 1; iB < BigBuf.Length; iB++) {
                    if(BufLok[iB] == false) {
                        if(BigBuf[iB] != null && BigBuf[iB].IsAlive) {
                            double[] C = (double[])BigBuf[iB].Target;
                            if(C != null && C.Length >= size) {
                                if(BB == null) {
                                    // passenden Buffer gefunden
                                    BB = C;
                                    iBuf = iB;
                                } else {
                                    if(C.Length < BB.Length) {
                                        // noch besser passenden Buffer gefunden
                                        BB = C;
                                        iBuf = iB;
                                    }
                                }
                            } else {
                                // Buffer ist tot, oder zu klein
                                // => könnte (re-)allokiert werden, falls nötig.
                                iAlloc = iB;
                            }
                        } else {
                            // Buffer ist tot
                            // => könnte reallokiert werden, falls nötig.
                            iAlloc = iB;
                        }
                    }
                }

                if(BB != null) {
                    _iBuf = iBuf;
                    BufLok[iBuf] = true;
                    return BB;
                }

                // (re-) allocation necessary
                if(iAlloc < 0) {
                    // pool exhausted => grow it
                    int oldLen = BigBuf.Length;
                    Array.Resize(ref BigBuf, oldLen * 2);
                    Array.Resize(ref BufLok, oldLen * 2);
                    iAlloc = Math.Max(1, oldLen);
                }

                _iBuf = iAlloc;
                BB = new double[size];
                BigBuf[iAlloc] = new WeakReference(BB);
                BufLok[iAlloc] = true;
                return BB;
            }
        }

        static void FreePooled(int ibuf) {
            lock(s_PoolLock) {
                // "Sperre" aufheben
                Debug.Assert(BufLok[ibuf] == true);
                BufLok[ibuf] = false;
                Debug.Assert(BufLok[ibuf] == false);
            }
        }

        // ===================================================================================
        // public interface
        // ===================================================================================

        /// <summary>
        /// Only use this if you know what you are doing: you're about to theave tha save .NET-world.
        /// Some <c>stackalloc</c> would be ideal for this, but it may create stack overflow...
        /// </summary>
        /// <param name="_iBuf">
        /// handle, to be passed to <see cref="FreeTempBuffer(int)"/> (on the same thread); always positive
        /// </param>
        /// <param name="size">minimum length of the returned buffer</param>
        /// <returns>
        /// a buffer of length greater or equal to <paramref name="size"/>; the content is undefined
        /// </returns>
        public static double[] GetTempBuffer(out int _iBuf, int size) {
            double[] ret = null;
            if(size <= MAX_THREADLOCAL_SIZE)
                ret = GetThreadLocal(out _iBuf, size);
            else
                _iBuf = -1;

            if(ret == null)
                ret = GetPooled(out _iBuf, size);

            CountGet();
            Debug.Assert(_iBuf > 0);
            return ret;
        }

        /// <summary>
        /// Releases a buffer obtained by <see cref="GetTempBuffer(out int, int)"/>;
        /// must be called on the same thread.
        /// </summary>
        public static void FreeTempBuffer(int ibuf) {
            if(ibuf >= THREADLOCAL_HANDLE_OFFSET)
                FreeThreadLocal(ibuf);
            else
                FreePooled(ibuf);
            t_NoOfBuffersInUse--;
        }

        /// <summary>
        /// Only use this if you know what you are doing: you're about to leave the save .NET-world.
        /// </summary>
        public static MultidimensionalArray GetTempMultidimensionalarray(out int iBuf, params int[] Lengths) {
            int size = 1;
            for(int i = Lengths.Length - 1; i >= 0; i--)
                size *= Lengths[i];

            double[] buffer = GetTempBuffer(out iBuf, size);

            var R = MultidimensionalArray.CreateWrapper(buffer, Lengths);
            R.Clear();
            return R;
        }

    }
}
