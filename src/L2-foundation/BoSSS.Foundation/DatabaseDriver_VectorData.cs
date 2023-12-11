﻿using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using Newtonsoft.Json;

namespace BoSSS.Foundation.IO {
    interface IVectorDataSerializer : ISerializer {
        Guid SaveVector<T>(IList<T> vector);

        void SaveVector<T>(IList<T> vector, Guid id);

        IList<T> LoadVector<T>(Guid id, ref Partitioning part);
    }

    class VectorDataSerializer : MPIProcess, IVectorDataSerializer {
        IFileSystemDriver m_fsDriver;
        ISerializer serializer;

        public string Name => serializer.Name;

        public VectorDataSerializer(IFileSystemDriver driver, ISerializer serializer) {
            this.serializer = serializer;
            m_fsDriver = driver;
        }

        /// <summary>
        /// header file for distributed stored vectors
        /// </summary>
        [Serializable]
        public class DistributedVectorHeader {

            /// <summary>
            /// GUID to identify the vector in the storage system.
            /// </summary>
            public Guid m_Guid;

            /// <summary>
            /// partition of the vector;
            /// the i-th entry is the first index of the vector that is stored in the i-th part.
            /// the length of this array is equal to the number of parts plus 1, the last entry
            /// is the total length of the vector.
            /// </summary>
            public long[] Partitioning;

            /// <summary>
            /// A hack, which indicates that we are using the most modern version.
            /// </summary>
            public bool UseGzip;
        }

        /// <summary>
        /// Saves a vector to the database, allocating a GUID automatically.
        /// </summary>
        /// <param name="vector"></param>
        /// <returns>
        /// The GUID that was allocated to identify the vector within the storage system.
        /// </returns>
        public Guid SaveVector<T>(IList<T> vector) {
            // allocate GUID
            // =============
            Guid id;
            id = Guid.NewGuid();
            id = id.MPIBroadcast(0);

            SaveVector(vector, id);

            return id;
        }

        bool DebugSerialization = false;


        /// <summary>
        /// Determines the average size after serialization, in bytes, of items stores in a vector;
        /// </summary>
        private int GetItemSize<T>(IList<T> vector) {
            int N = Math.Min(10, vector.Count);
            var testVector = vector.GetSubVector(0, N);

            long size;
            using (MemoryStream ms = new MemoryStream(1024*1024)) {
                serializer.Serialize(ms, testVector, typeof(T));
                size = ms.ToArray().Length;
            }

            return (int)Math.Round((double)size / (double)N);
        }

        /// <summary>
        /// It was noticed that the BSON-de-serializer used in this project crashes for objects
        /// which are, after serialization, larger than about 100 megabyte;
        /// therefore, such vectors will be split up further before saving.
        /// </summary>
        private IList<T>[] RepartitionVector<T>(IList<T> vector) {

            const int SplitUpThreshold = 64 * 1024 * 1024;

            int N = vector.Count;
            int itmSz = GetItemSize(vector);

            int totSize = N * itmSz;
            if (totSize < SplitUpThreshold) {
                return new IList<T>[] { vector };
            } else {
                int NoOfParts = (int)Math.Ceiling((double)totSize / (double)SplitUpThreshold);
                var ret = new IList<T>[NoOfParts];

                for (int iPart = 0; iPart < NoOfParts; iPart++) {
                    int i0 = checked((int)(((long)N * (long)iPart) / NoOfParts));
                    int iE = checked((int)(((long)N * (long)(iPart + 1)) / NoOfParts));

                    ret[iPart] = vector.GetSubVector(i0, iE - i0);
                }

                return ret;
            }
        }


        /// <summary>
        /// saves a vector to the database, under a specified Guid
        /// </summary>
        /// <param name="vector">
        /// the part of the vector which is stored on the local process.
        /// </param>
        /// <param name="id">
        /// the Guid under which the vector should be stored; must be the same
        /// on all MPI processes.
        /// </param>
        public void SaveVector<T>(IList<T> vector, Guid id) {
            using (new FuncTrace()) {
                if (m_fsDriver == null)
                    throw new NotSupportedException("Can't save data when not connected to database.");

#if DEBUG
                // never trust the user
                Guid id0 = id.MPIBroadcast(0);
                if (id0 != id) {
                    throw new ArgumentException(
                        "Guid differs at least on one MPI process.",
                        nameof(id));
                }
#endif
                Exception e = null;
                try {
                    // split up the vector to avoid files which are to big for the JSON/BSON format:
                    IList<T>[] vectorRepart = RepartitionVector<T>(vector);
                    if (vectorRepart.Length < 1)
                        throw new ApplicationException("internal error in saving vector");

                    // determine partitioning
                    int[] NoOfPartsPerProc = vectorRepart.Length.MPIAllGather();
                    int[] _Partitioning = vectorRepart.Select(part => part.Count).ToArray().MPIAllGatherv();
                    int MyPartOffset = 0;
                    for (int i = 0; i < MyRank; i++)
                        MyPartOffset += NoOfPartsPerProc[i];

                    // loop over all parts...
                    for (int iPart = 0; iPart < vectorRepart.Length; iPart++) {

                        // save parts
                        using (Stream s = m_fsDriver.GetDistVectorDataStream(true, id, MyPartOffset + iPart)) {
                            using (var s2 = new GZipStream(s, CompressionMode.Compress)) {
                                // Use a tuple since the Json format expects one object per
                                // file (with some tricks, we could avoid this, but it's
                                // not really worth the effort)
                                Tuple<DistributedVectorHeader, IList<T>> tuple;


                                if (MyRank == 0 && iPart == 0) {
                                    // save header (only proc 0, for the 0-th part)
                                    DistributedVectorHeader header = new DistributedVectorHeader();
                                    header.m_Guid = id;
                                    header.UseGzip = !DebugSerialization;
                                    long[] part = new long[_Partitioning.Length + 1];
                                    for (int i = 0; i < _Partitioning.Length; i++) {
                                        part[i + 1] = part[i] + _Partitioning[i];
                                    }
                                    header.Partitioning = part;
                                    tuple = new Tuple<DistributedVectorHeader, IList<T>>(header, vectorRepart[iPart]);
                                } else {
                                    // save all other parts without header
                                    tuple = new Tuple<DistributedVectorHeader, IList<T>>(null, vectorRepart[iPart]);
                                }

                                // serialize it!
                                serializer.Serialize(s2, tuple, typeof(Tuple<DistributedVectorHeader, IList<T>>));

                                s2.Close();
                                s.Close();
                            }
                        }
                    }
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                } catch (Exception ee) {
                    Console.Error.WriteLine(ee.GetType().Name + " on rank " + this.MyRank + " saving vector " + id + ": " + ee.Message);
                    Console.Error.WriteLine(ee.StackTrace);
                    e = ee;
                }

                //e.ExceptionBcast();
            }
        }

        /// <summary>
        /// Loads a vector from the database
        /// </summary>
        /// <param name="id"></param>
        /// <param name="part">
        /// Optional partition of the vector among MPI processors: if null, a
        /// partition is defined by the loader logic.
        /// </param>
        public IList<T> LoadVector<T>(Guid id, ref Partitioning part) {
            using (new FuncTrace()) {

                // load header
                // -----------
                DistributedVectorHeader header = null;

                // load on p0
                Exception e = null;
                try {
                    if (MyRank == 0) {
                        // load header
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, 0))
                        using (var reader = new GZipStream(s, CompressionMode.Decompress)) {
                            header = ((Tuple<DistributedVectorHeader, IList<T>>)serializer.Deserialize(reader, typeof(Tuple<DistributedVectorHeader, IList<T>>))).Item1;
                            reader.Close();
                            s.Close();
                        }

                    }
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                // broadcast
                header = header.MPIBroadcast(0);

                // define partition, if necessary
                if (part == null) {
                    long L = header.Partitioning.Last();
                    long i0 = (L * this.MyRank) / this.Size;
                    long ie = (L * (this.MyRank + 1)) / this.Size;
                    part = new Partitioning(checked((int)(ie - i0)));
                }


                // check size
                if (part.TotalLength != header.Partitioning.Last())
                    throw new ArgumentException("IO error: length of vector to load differs from total length of list.");


                // find first part
                // ===============

                long myI0;
                int myLen;
                myI0 = part.i0;
                myLen = part.LocalLength;

                int p = 0;
                while (header.Partitioning[p] <= myI0) {
                    p++;
                }
                p--; // now p is the index of the first vector part

                // load parts
                // ==========
                List<T> ret = new List<T>(part.LocalLength);
                try {
                    int d = 0;
                    do {
                        IList<T> vecP;
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, p))
                        using (var reader = new GZipStream(s, CompressionMode.Decompress)) {
                            vecP = ((Tuple<DistributedVectorHeader, IList<T>>)serializer.Deserialize(reader, typeof(Tuple<DistributedVectorHeader, IList<T>>))).Item2;
                            reader.Close();
                            s.Close();
                        }

                        long srdInd = ((myI0 + d) - header.Partitioning[p]);
                        long CopyLen = Math.Min(vecP.Count - srdInd, myLen - d);

                        for (int i = 0; i < CopyLen; i++) {
                            ret.Add(vecP[checked((int)(srdInd + i))]);
                            d++;
                        }

                        p++; // next part;

                    } while (d < myLen);
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                // return
                // ======
                return ret;
            }
        }

        /*
        /// <summary>
        /// Loads a vector from the database
        /// </summary>
        /// <param name="id"></param>
        /// <param name="part">
        /// Optional partition of the vector among MPI processors: if null, a
        /// partition is defined by the loader logic.
        /// </param>
        public IList<T> LoadVectorSpeziale<T>(Guid id, ref Partitioning part, out Partitioning verySpecial)
        {
            using (new FuncTrace())
            {

                // load header
                // -----------
                DistributedVectorHeader header = null;

                // load on p0
                Exception e = null;
                try
                {
                    if (MyRank == 0)
                    {
                        // load header
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, 0))
                        using (var reader = new GZipStream(s, CompressionMode.Decompress))
                        {
                            header = ((Tuple<DistributedVectorHeader, IList<T>>)serializer.Deserialize(reader, typeof(Tuple<DistributedVectorHeader, IList<T>>))).Item1;
                            reader.Close();
                            s.Close();
                        }

                    }
                }
                catch (Exception ee)
                {
                    e = ee;
                }
                e.ExceptionBcast();

                // broadcast
                header = header.MPIBroadcast(0);

                // define partition, if necessary
                if (part == null)
                {
                    int L = (int)header.Partitioning.Last();
                    int i0 = (L * this.MyRank) / this.Size;
                    int ie = (L * (this.MyRank + 1)) / this.Size;
                    part = new Partitioning(ie - i0);
                }

                verySpecial = new Partitioning((int)(header.Partitioning[this.MyRank + 1] - header.Partitioning[this.MyRank]));

                // check size
                if (part.TotalLength != header.Partitioning.Last())
                    throw new ArgumentException("IO error: length of vector to load differs from total length of list.");


                // find first part
                // ===============

                long myI0;
                int myLen;
                myI0 = 0;
                myLen = (int)header.Partitioning.Last();

                int p = 0;
                while (header.Partitioning[p] <= myI0)
                {
                    p++;
                }
                p--; // now p is the index of the first vector part

                // load parts
                // ==========
                List<T> ret = new List<T>(part.LocalLength);
                try
                {
                    int d = 0;
                    do
                    {
                        IList<T> vecP;
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, p))
                        using (var reader = new GZipStream(s, CompressionMode.Decompress))
                        {
                            vecP = ((Tuple<DistributedVectorHeader, IList<T>>)serializer.Deserialize(reader, typeof(Tuple<DistributedVectorHeader, IList<T>>))).Item2;
                            reader.Close();
                            s.Close();
                        }

                        int srdInd = (int)((myI0 + d) - header.Partitioning[p]);
                        int CopyLen = Math.Min(vecP.Count - srdInd, myLen - d);

                        for (int i = 0; i < CopyLen; i++)
                        {
                            ret.Add(vecP[srdInd + i]);
                            d++;
                        }

                        p++; // next part;

                    } while (d < myLen);
                }
                catch (Exception ee)
                {
                    e = ee;
                }
                e.ExceptionBcast();

                // return
                // ======
                return ret;
            }
        }
        */

        public object Deserialize(Stream stream, Type objectType) {
            return serializer.Deserialize(stream, objectType);
        }

        public void Serialize(Stream stream, object obj, Type objectType) {
            serializer.Serialize(stream, obj, objectType);
        }
    }
}
