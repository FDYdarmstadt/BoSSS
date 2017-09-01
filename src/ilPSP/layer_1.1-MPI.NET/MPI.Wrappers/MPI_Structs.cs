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
using System.Runtime.InteropServices;

namespace MPI.Wrappers {

    /// <summary>
    /// MPI_Status struct;
    /// </summary>
    /// <remarks>
    /// In both currently considered MPI implementations
    /// (current = 16nov2010, MS-MPI, MPICH2, OpenMPI 1.4.1)
    /// the MPI_Struct contains 5 integers, but the sequence is different;
    /// The public properties use the <see cref="IMPIdriver.Flavour"/> to
    /// return the correct value.
    /// </remarks>
    [StructLayout(LayoutKind.Sequential)]
    public struct MPI_Status {
        internal int i1;
        internal int i2;
        internal int i3;
        internal int i4;
        internal int i5;
        internal int i6;

        /// <summary>
        /// 
        /// </summary>
        public int count {
            get {
                switch (csMPI.Raw.Flavour) {
                    case MPI_Flavour.MPICH:
                        return i1;
                    case MPI_Flavour.OpenMPI:
                        return i4;
                    default:
                        throw new NotImplementedException("unknown MPI flavor.");
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int cancelled {
            get {
                switch (csMPI.Raw.Flavour) {
                    case MPI_Flavour.MPICH:
                        return i2;
                    case MPI_Flavour.OpenMPI:
                        return i5;
                    default:
                        throw new NotImplementedException("unknown MPI flavor.");
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int MPI_SOURCE {
            get {
                switch (csMPI.Raw.Flavour) {
                    case MPI_Flavour.MPICH:
                        return i3;
                    case MPI_Flavour.OpenMPI:
                        return i1;
                    default:
                        throw new NotImplementedException("unknown MPI flavor.");
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int MPI_TAG {
            get {
                switch (csMPI.Raw.Flavour) {
                    case MPI_Flavour.MPICH:
                        return i4;
                    case MPI_Flavour.OpenMPI:
                        return i2;
                    default:
                        throw new NotImplementedException("unknown MPI flavor.");
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int MPI_ERROR {
            get {
                switch (csMPI.Raw.Flavour) {
                    case MPI_Flavour.MPICH:
                        return i5;
                    case MPI_Flavour.OpenMPI:
                        return i3;
                    default:
                        throw new NotImplementedException("unknown MPI flavor.");
                }
            }
        }

    }

    /// <summary>
    /// MPI operation type
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct MPI_Op {

        /// <summary>
        /// initializes this communicator form an int;
        /// </summary>
        public MPI_Op(int m) {
            this.m1 = m;
        }

        internal int m1;

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            if (obj.GetType() != typeof(MPI_Op))
                return false;
            MPI_Op other = (MPI_Op)obj;
            return (other.m1 == this.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() {
            return m1;
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="l"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static bool operator !=(MPI_Op l, MPI_Op r) {
            return (l.m1 != r.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="l"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static bool operator ==(MPI_Op l, MPI_Op r) {
            return (l.m1 == r.m1);
        }
    }


    /// <summary>
    /// MPI communicator type
    /// </summary>
    /// <remarks>
    /// Note that within the .NET - wrappers we always use FORTRAN MPI
    /// communicators, which are, by the standard, defined as the FORTRAN
    /// integer and therefore, hopefully, 4 bytes long on all platforms;<br/>
    /// The problem with the C - MPI communicator is the following: in MPICH2
    /// (and derivatives, including Microsoft MPI), the MPI communicator is
    /// defined as int, and therefore always 4 bytes, no matter if we run in
    /// 64 or 32 bit mode. On the other hand, in OpenMPI, the MPI communicator
    /// is a pointer, and its size depends on whether we run in 32 or 64 bit mode;
    /// </remarks>
    [StructLayout(LayoutKind.Sequential)]
    public struct MPI_Comm {

        /// <summary>
        /// The encapsaulted communicator
        /// </summary>
        public int m1;

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            if (obj.GetType() != typeof(MPI_Comm))
                return false;
            MPI_Comm other = (MPI_Comm)obj;
            return (other.m1 == this.m1);
        }

        /// <summary>
        /// hash value;
        /// </summary>
        public override int GetHashCode() {
            return m1;
        }

        /// <summary>
        /// bitwise in-equivalence;
        /// </summary>
        public static bool operator !=(MPI_Comm l, MPI_Comm r) {
            return (l.m1 != r.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        public static bool operator ==(MPI_Comm l, MPI_Comm r) {
            return (l.m1 == r.m1);
        }
    }

    /// <summary>
    /// MPI Datatype - type
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct MPI_Datatype {

        /// <summary>
        /// initializes this structure from an integer
        /// </summary>
        public MPI_Datatype(int i) {
            this.m1 = i;
        }

        /// <summary>
        /// The actual data-type
        /// </summary>
        public int m1;

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            if (obj.GetType() != typeof(MPI_Datatype))
                return false;
            MPI_Datatype other = (MPI_Datatype)obj;
            return (other.m1 == this.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        public override int GetHashCode() {
            return m1;
        }

        /// <summary>
        /// bitwise in-equivalence;
        /// </summary>
        public static bool operator !=(MPI_Datatype l, MPI_Datatype r) {
            return (l.m1 != r.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        public static bool operator ==(MPI_Datatype l, MPI_Datatype r) {
            return (l.m1 == r.m1);
        }
    }

    /// <summary>
    /// MPI Datatype - type
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct MPI_Request {

        /// <summary>
        /// The actual request
        /// </summary>
        public int m1;

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            if (obj.GetType() != typeof(MPI_Datatype))
                return false;
            MPI_Datatype other = (MPI_Datatype)obj;
            return (other.m1 == this.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() {
            return m1;
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="l"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static bool operator !=(MPI_Request l, MPI_Request r) {
            return (l.m1 != r.m1);
        }

        /// <summary>
        /// bitwise equivalence;
        /// </summary>
        /// <param name="l"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static bool operator ==(MPI_Request l, MPI_Request r) {
            return (l.m1 == r.m1);
        }
    }
}
