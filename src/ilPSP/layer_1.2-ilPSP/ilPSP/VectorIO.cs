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
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using MPI.Wrappers;

namespace ilPSP.Utils {

    /// <summary>
    /// Input/Output utility functions for vectors (all kinds of double-enumerables);
    /// </summary>
    public static class VectorIO {

        /// <summary>
        /// saves a vector (distributed over all MPI processes in the WORLD communicator) into one text file
        /// </summary>
        public static void SaveToTextFile(this IEnumerable<double> list, string filename) {
            SaveToTextFile(list, filename, csMPI.Raw._COMM.WORLD);
        }

        /// <summary>
        /// saves a vector (distributed over all MPI processes in the WORLD communicator) into one text file
        /// </summary>
        public static void SaveToTextFile(this IEnumerable<double> list, string filename, MPI_Comm comm) {
            SaveToTextFile<double>(list, filename, comm, x => x.ToString(NumberFormatInfo.InvariantInfo));
        }

        /// <summary>
        /// saves a vector (distributed over various MPI processes) into one text file
        /// </summary>
        /// <param name="list"></param>
        /// <param name="filename"></param>
        /// <param name="ToString">
        /// Optional customization of the way how the entries of <paramref name="list"/> should be transformed into strings.
        /// </param>
        public static void SaveToTextFile<T>(this IEnumerable<T> list, string filename, Func<T, string> ToString = null) {
            SaveToTextFile<T>(list, filename, csMPI.Raw._COMM.WORLD, ToString);
        }

        /// <summary>
        /// saves a vector (distributed over various MPI processes) into one text file
        /// </summary>
        /// <param name="list"></param>
        /// <param name="fileName"></param>
        /// <param name="comm"></param>
        /// <param name="ToString">
        /// Optional customization of the way how the entries of <paramref name="list"/> should be transformed into strings.
        /// </param>
        public static void SaveToTextFile<T>(this IEnumerable<T> list, string fileName , MPI_Comm comm, Func<T, string> ToString = null) {
            int Rank;
            csMPI.Raw.Comm_Rank(comm, out Rank);

            // open file
            StreamWriter stw;
            if (Rank == 0) {
                stw = new StreamWriter(fileName);
            } else {
                stw = null;
            }

            try {
                SaveToStream<T>(list, stw, comm, ToString);
            } finally {
                if (Rank == 0) {
                    stw.Flush();
                    stw.Close();
                }
            }

        }

        /// <summary>
        /// saves a vector (distributed over various MPI processes) into one text file
        /// </summary>
        /// <param name="list"></param>
        /// <param name="stw">text writer on output stream</param>
        /// <param name="comm"></param>
        /// <param name="ToString">
        /// Optional customization of the way how the entries of <paramref name="list"/> should be transformed into strings.
        /// </param>
        public static void SaveToStream<T>(this IEnumerable<T> list, TextWriter stw, MPI_Comm comm, Func<T, string> ToString = null) {
            int Rank, Size;
            csMPI.Raw.Comm_Rank(comm, out Rank);
            csMPI.Raw.Comm_Size(comm, out Size);

            if (ToString == null)
                ToString = delegate(T obj) {
                    return obj.ToString();
                };

            var sms = new SerialisationMessenger(comm);

            if (Rank == 0) {
                sms.CommitCommPaths();

                SortedDictionary<int, T[]> packets = new SortedDictionary<int, T[]>();
                packets.Add(0, list.ToArray());


                T[] rcvdata;
                int rcvRank;
                while (sms.GetNext(out rcvRank, out rcvdata)) {
                    packets.Add(rcvRank, rcvdata);
                }



                // write

                for (int r = 0; r < Size; r++) {
                    IList<T> lr = packets[r];

                    for (int n = 0; n < lr.Count; n++) {
                        stw.Write(ToString(lr[n]));
                        stw.WriteLine();
                    }
                }

            } else {
                sms.SetCommPath(0);
                sms.CommitCommPaths();

                T[] send = list as T[];
                if (send == null)
                    send = list.ToArray();

                sms.Transmitt(0, send);


                T[] dummy;
                int dummy_;
                if (sms.GetNext<T[]>(out dummy_, out dummy))
                    throw new ApplicationException("error in app");
            }
        }

        ///// <summary>
        ///// saves a vector (distributed over various MPI processes) into a binary file
        ///// </summary>
        ///// <param name="list"></param>
        ///// <param name="filename"></param>
        //public static void SaveToFile<T>(this IEnumerable<T> list, string filename) {
        //    SaveToFile(list, filename, csMPI.Raw._COMM.WORLD);
        //}

        ///// <summary>
        ///// saves a vector (distributed over various MPI processes) into a binary file
        ///// </summary>
        ///// <param name="list"></param>
        ///// <param name="filename"></param>
        ///// <param name="comm"></param>
        //public static void SaveToFile<T>(this IEnumerable<T> list, string filename, MPI_Comm comm) {
        //    int Rank, Size;
        //    csMPI.Raw.Comm_Rank(comm, out Rank);
        //    csMPI.Raw.Comm_Size(comm, out Size);


        //    var sms = new SerialisationMessenger(comm);

        //    if (Rank == 0) {
        //        sms.CommitCommPaths();

        //        SortedDictionary<int, T[]> packets = new SortedDictionary<int, T[]>();
        //        packets.Add(0, list.ToArray());


        //        T[] rcvdata;
        //        int rcvRank;
        //        while (sms.GetNext(out rcvRank, out rcvdata)) {
        //            packets.Add(rcvRank, rcvdata);
        //        }

        //        int TotLen = packets.Values.Sum(p => p.Length);
        //        T[] allData = new T[TotLen];


        //        // write
        //        int offset = 0;
        //        for (int r = 0; r < Size; r++) {
        //            var lr = packets[r];
        //            Array.Copy(lr, 0, allData, offset, lr.Length);
        //            offset += lr.Length;
        //        }

        //        using (var file = new FileStream(filename, FileMode.Create, FileAccess.Write)) {
        //            Serializer.Serialize(file, allData);
        //            file.Close();
        //        }
        //    } else {
        //        sms.SetCommPath(0);
        //        sms.CommitCommPaths();

        //        T[] send = list as T[];
        //        if (send == null)
        //            send = list.ToArray();

        //        sms.Transmitt(0, send);


        //        T[] dummy;
        //        int dummy_;
        //        if (sms.GetNext<T[]>(out dummy_, out dummy))
        //            throw new ApplicationException("error in app");
        //    }
        //}

        ///// <summary>
        ///// Loads a vector from a binary file
        ///// </summary>
        ///// <param name="filename"></param>
        ///// <param name="comm"></param>
        ///// <param name="p">
        ///// optional partition, its MPI communicator must match <paramref name="comm"/>
        ///// </param>
        ///// <returns></returns>
        //public static T[] LoadFromFile<T>(string filename, MPI_Comm comm, Partitioning p = null) {
        //    if (p != null && (comm != p.MPI_Comm))
        //        throw new ArgumentException();
        //    int Rank, Size;
        //    csMPI.Raw.Comm_Rank(comm, out Rank);
        //    csMPI.Raw.Comm_Size(comm, out Size);



        //    T[] vec = null;
        //    int GlobalLen = 0;
        //    if (Rank == 0) {
        //        using (var file = new FileStream(filename, FileMode.Open, FileAccess.Read)) {
        //            vec = Serializer.Deserialize<T[]>(file);
        //            file.Close();
        //        }

        //        GlobalLen = vec.Length;
        //    }

        //    unsafe {
        //        csMPI.Raw.Bcast((IntPtr)(&GlobalLen), 1, csMPI.Raw._DATATYPE.INT, 0, comm);
        //    }

        //    if (p == null) {
        //        int i0 = (GlobalLen * Rank) / Size;
        //        int iE = (GlobalLen * (Rank + 1)) / Size;
        //        p = new Partitioning(iE - i0, 1, comm);
        //    } else {
        //        if (p.TotalLength != GlobalLen)
        //            throw new ArgumentException("Mismatch between partition length and number of elements in file.");
        //    }

        //    Dictionary<int, T[]> SendDataPackets = new Dictionary<int, T[]>();
        //    if (Rank != 0) {
        //        for (int prc = 1; prc < Size; prc++) {
        //            T[] vec_prc = vec.GetSubVector((int)p.GetI0Offest(prc), p.GetLocalLength(prc));
        //            SendDataPackets.Add(prc, vec_prc);
        //        }
        //    }
        //    var RcvdDataPackets = SerialisationMessenger.ExchangeData(SendDataPackets, comm);

        //    if (Rank == 0) {
        //        return vec.GetSubVector((int)p.GetI0Offest(0), p.GetLocalLength(0));
        //    } else {
        //        Debug.Assert(RcvdDataPackets.Count == 1);
        //        Debug.Assert(RcvdDataPackets.ContainsKey(0));
        //        return RcvdDataPackets[0];
        //    }
        //}


        /// <summary>
        /// Loads a vector from a text file
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="Parse">Parsing function, converts a string into <typeparamref name="T"/>.</param>
        /// <param name="comm"></param>
        /// <param name="p">
        /// optional partition, its MPI communicator must match <paramref name="comm"/>
        /// </param>
        /// <returns></returns>
        public static T[] LoadFromTextFile<T>(string filename, Func<string, T> Parse, MPI_Comm comm, Partitioning p = null) {
            if (p != null && (comm != p.MPI_Comm))
                throw new ArgumentException();
            int Rank, Size;
            csMPI.Raw.Comm_Rank(comm, out Rank);
            csMPI.Raw.Comm_Size(comm, out Size);



            List<T> vec = null;
            int GlobalLen = 0;
            if (Rank == 0) {

                StreamReader rd = null;
                vec = new List<T>();
                try {
                    rd = new StreamReader(filename);

                    string line = rd.ReadLine();
                    while (line != null) {
                        T entry = Parse(line);
                        vec.Add(entry);
                        line = rd.ReadLine();
                    }

                } finally {
                    if (rd != null) {
                        rd.Close();
                    }
                }
                GlobalLen = vec.Count;
            }

            unsafe {
                csMPI.Raw.Bcast((IntPtr)(&GlobalLen), 1, csMPI.Raw._DATATYPE.INT, 0, comm);
            }

            if (p == null) {
                int i0 = (GlobalLen * Rank) / Size;
                int iE = (GlobalLen * (Rank + 1)) / Size;
                p = new Partitioning(iE - i0, comm);
            } else {
                if (p.TotalLength != GlobalLen)
                    throw new ArgumentException("Mismatch between partition length and number of elements in file.");
            }

            Dictionary<int, T[]> SendDataPackets = new Dictionary<int, T[]>();
            if (Rank != 0) {
                for (int prc = 1; prc < Size; prc++) {
                    T[] vec_prc = vec.GetSubVector((int)p.GetI0Offest(prc), p.GetLocalLength(prc));
                    SendDataPackets.Add(prc, vec_prc);
                }
            }
            var RcvdDataPackets = SerialisationMessenger.ExchangeData(SendDataPackets, comm);

            if (Rank == 0) {
                return vec.GetSubVector((int)p.GetI0Offest(0), p.GetLocalLength(0));
            } else {
                Debug.Assert(RcvdDataPackets.Count == 1);
                Debug.Assert(RcvdDataPackets.ContainsKey(0));
                return RcvdDataPackets[0];
            }
        }

        /// <summary>
        /// Loads a vector from a text file.
        /// </summary>
        public static double[] LoadFromTextFile(string filename, MPI_Comm comm, Partitioning p = null) {
            return LoadFromTextFile<double>(filename, str => double.Parse(str, NumberFormatInfo.InvariantInfo), comm);
        }

        /// <summary>
        /// Loads a vector from a text file.
        /// </summary>
        public static void LoadFromTextFile(this IList<double> vec, string filename) {
            var part = new Partitioning(vec.Count);
            var load = LoadFromTextFile<double>(filename, str => double.Parse(str, NumberFormatInfo.InvariantInfo), part.MPI_Comm, part);

            for (int i = 0; i < load.Length; i++) {
                vec[i] = load[i];
            }
        }

        ///// <summary>
        ///// Loads a vector from a binary file.
        ///// </summary>
        //public static void LoadFromFile(this IList<double> vec, string filename, MPI.Wrappers.MPI_Comm c) {
        //    var part = new Partitioning(vec.Count, 1, c);
        //    var load = LoadFromFile<double>(filename, part.MPI_Comm, part);

        //    for (int i = 0; i < load.Length; i++) {
        //        vec[i] = load[i];
        //    }
        //}

        ///// <summary>
        ///// Loads a vector from a binary file, acting on the world communicator.
        ///// </summary>
        //public static void LoadFromFile(this IList<double> vec, string filename) {
        //    LoadFromFile(vec, filename, csMPI.Raw._COMM.WORLD);
        //}
    }
}
