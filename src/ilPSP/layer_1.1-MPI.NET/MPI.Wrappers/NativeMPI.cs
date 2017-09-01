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


namespace MPI.Wrappers {

    /// <summary>
    /// Identifies the MPI implementation
    /// </summary>
    public enum MPI_Flavour {
        /// <summary>
        /// MPICH2 or Microsoft-MPI
        /// </summary>
        MPICH,

        /// <summary>
        /// OpenMPI (Unix-Only)
        /// </summary>
        OpenMPI
    }

    /*

    /// <summary>
    /// Direct, low-level interface to the system MPI library.
    /// </summary>
    /// <remarks>
    /// This low-level interface provides direct access to the unmanaged
    /// MPI library provided by the system. It is by nature unsafe, and
    /// should only be used by programmers experienced both in the use 
    /// of MPI from lower-level languages (e.g., C, Fortran) and with an
    /// understanding of the interaction between managed and unmanaged
    /// code, especially those issues that pertain to memory 
    /// pinning/unpinning.<br/>
    /// The purpose of this class is to build a standardized API on top
    /// of the different MPI implementations (OpenMPI and MPICH);
    /// This is a bit difficult, because this implementations are not binary equivalent;
    /// E.g. MPICH defines the MPI_Comm datatype as int, which is always 4 Byte, inependent wehter the application runs in 32- or 64 - bit mode;
    /// OpenMPI, on the other hand, defines MPI_Comm as a pointer, so it depends on the Bitness of the application/platform;
    /// We solved such issues by addressing the FORTRAN interface, where types like MPI_Comm, MPI_Op, ... are defined to be integers
    /// by the standard;
    /// </remarks>
    public class csMPI.Raw  {

        /// <summary>
        /// converts the MPI communicators in BoSSS (which are FORTRAN-communicators) into a C-Communicator
        /// </summary>
        /// <param name="EightByteComm">
        /// on exit, the value of the C - MPI communicator handle, if <see cref="GetSizeof_C_MPI_comm"/>==8,
        /// otherwise 0;
        /// </param>
        /// <param name="FourByteComm">
        /// on exit, the value of the C - MPI communicator handle, if <see cref="GetSizeof_C_MPI_comm"/>==4,
        /// otherwise 0;
        /// </param>
        /// <param name="input">
        /// (FORTRAN) MPI communicator handle to convert
        /// </param>
        /// <returns>
        /// Equal to <see cref="GetSizeof_C_MPI_comm"/>, 
        /// either 4 or 8
        /// </returns>
        /// <remarks>
        /// The length of the C - MPI communicator handle depends on the used MPI implementation 
        /// (see remarks at <see cref="T_MPI_COMM"/>) 
        /// and is given by to the return value of <see cref="GetSizeof_C_MPI_comm"/>
        /// </remarks>
        static public int MPI_Comm_f2c(MPI_Comm input, out uint FourByteComm, out ulong EightByteComm) {
            unsafe {
                byte[] C_Comm = _COMM.Comm_f2c(input);
                fixed (byte* pComm = &(C_Comm[0])) {
                    uint _FourByteComm = 0;
                    ulong _EightByteComm = 0;
                    byte* pDest;
                    switch (C_Comm.Length) {
                        case 4: pDest = (byte*)(&_FourByteComm); break;
                        case 8: pDest = (byte*)(&_EightByteComm); break;
                        default: throw new NotImplementedException("unknown size of C MPI communicator: " + C_Comm.Length + "bytes; ");
                    }
                    for (int i = 0; i < C_Comm.Length; i++)
                        pDest[i] = pComm[i];

                    FourByteComm = _FourByteComm;
                    EightByteComm = _EightByteComm;
                }
                return C_Comm.Length;
            }
        }

        //// workaround for .NET bug:
        //// https://connect.microsoft.com/VisualStudio/feedback/details/635365/runtimehelpers-initializearray-fails-on-64b-framework
        //static PlatformID[] Helper() {
        //    PlatformID[] p = new PlatformID[3];
        //    p[0] = PlatformID.Win32NT;
        //    p[1] = PlatformID.Unix;
        //    p[2] = PlatformID.Unix;
        //    return p;
        //}

        static MPI_Flavour m_Flavour = MPI_Flavour.MPICH;

        /// <summary>
        /// identifies the currently loaded MPI implementation
        /// </summary>
        static public MPI_Flavour Flavour {
            get {
                return m_Flavour;
            }
        }


        static IMPI_CommConstants m_MPI_Comm = new Platform_Native_Comm();

        /// <summary>
        /// default communicators
        /// </summary>
        public static IMPI_CommConstants _COMM {
            get {
                return m_MPI_Comm;
            }
        }

        static IMiscConstants m_MiscConstants = new Platform_Native_MiscConstants();

        /// <summary>
        /// misc MPI - constants
        /// </summary>
        public static MPI.Wrappers.IMiscConstants MiscConstants {
            get {
                return m_MiscConstants;
            }
        }

        /// <summary>
        /// predefined MPI operations
        /// </summary>
        public static IMPI_OpConstants _OP {
            get {
                return m_MPI_OP;
            }
        }

        static IMPI_OpConstants m_MPI_OP = new Platform_Native_MPI_Op();

        static IMPI_DatatypeConstants m_MPI_Datatype = new Platform_Native_MPI_Datatype();

        /// <summary>
        /// predefined MPI operations
        /// </summary>
        public static IMPI_DatatypeConstants _DATATYPE {
            get {
                return m_MPI_Datatype;
            }
        }


        /// <summary>
        /// gets an error string for an MPI error code
        /// </summary>
        /// <param name="errorCode">
        /// MPI error code (see <see cref="MPIException.ErrorCode"/>
        /// </param>
        /// <returns>
        /// the description of the error no <see cref="errorCode"/>
        /// from the MPI library
        /// </returns>
        static public string Error_string(int errorCode) {
            //unsafe {
            //    char* str = (char*)Marshal.AllocHGlobal(65000);
            //    int resultlen;
            //    instance.MPI_Error_string(errorCode, str, out resultlen);
            //    string r = Marshal.PtrToStringAnsi((IntPtr)str);
            //    Marshal.FreeHGlobal((IntPtr)str);
            //    return r;
            //    //return "MPI_ERROR_String not loaded.";
            //}
            return string.Format("Error_string not implemented yet; code is {0}", errorCode);
        }

        ///// <summary>
        ///// if true, <see cref="OpenMPI_MPI_Init"/> is taken instead 
        ///// of <see cref="MPI_Init"/>;
        ///// For some strange reason that i cannot reproduce, MPI_Init does not work when
        ///// loaded via dlsym(...);
        ///// </summary>
        //bool dirrtyMpiInitHack = false;

        ///// <summary>
        ///// see <see cref="dirrtyMpiInitHack"/>;
        ///// </summary>
        //[DllImport("mpi_f77", EntryPoint = "MPI_Init")]
        //static extern private unsafe int OpenMPI_MPI_Init(int* argc, byte*** argv);


        unsafe delegate int _MPI_Init(int* argc, byte*** argv);

        [DllImport("Platform_Native")]
        static extern int BoSSS_MPI_INIT(ref int ierr);


        /// <summary>
        /// see MPI reference;
        /// </summary>
        /// <param name="args">command line args passed to the application</param>
        /// <remarks>
        /// Do not call this directly, use <see cref="Enviroment.Bootstrap"/> instead.
        /// </remarks>
        static public void Init(string[] args) {
            //IntPtr[] argsAnsi = new IntPtr[Math.Max(args.Length, 1)];

            //for (int i = 0; i < args.Length; i++) {
            //    argsAnsi[i] = Marshal.StringToCoTaskMemAnsi(args[i]);
            //}

            //unsafe {
            //    fixed (void* argv = &(argsAnsi[0])) {
            //        void* _argv = argv;
            //        int argc = args.Length;
            //        if (instance.dirrtyMpiInitHack)
            //            MPIException.CheckReturnCode(OpenMPI_MPI_Init(&argc, (byte***)(&_argv)));
            //        else
            //            MPIException.CheckReturnCode(instance.MPI_Init(&argc, (byte***)(&_argv)));
            //    }
            //}

            //foreach (IntPtr p in argsAnsi) {
            //    Marshal.FreeCoTaskMem(p);
            //}
#if MPI_TRACE
            Console.WriteLine("MPI_Init");
#endif
            int ierr = 0;
            BoSSS_MPI_INIT(ref ierr);
            MPIException.CheckReturnCode(ierr);
        }


        [DllImport("Platform_Native")]
        static extern int BoSSS_MPI_Finalize();

        /// <summary>
        /// MPI finalize
        /// </summary>
        static public void mpiFinalize() {
#if MPI_TRACE
            Console.WriteLine("MPI_Finalize");
#endif
            MPIException.CheckReturnCode(BoSSS_MPI_Finalize());
        }



        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_COMM_RANK(ref MPI_Comm comm, out int rank, out int ierr);



        /// <summary>
        /// 
        /// </summary>
        static public void Comm_Rank(MPI_Comm comm, out int rank) {
            int ierr;
#if MPI_TRACE
            //Console.WriteLine("MPI_Comm_Rank");
#endif
            BoSSS_MPI_COMM_RANK(ref comm, out rank, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_COMM_SIZE(ref MPI_Comm comm, out int size, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Comm_Size(MPI_Comm comm, out int size) {
            int ierr;
#if MPI_TRACE
            //Console.WriteLine("MPI_Comm_Size");
#endif
            BoSSS_MPI_COMM_SIZE(ref comm, out size, out ierr);
            MPIException.CheckReturnCode(ierr);
        }


        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_INITIALIZED(ref int flag, out int ierr);


        /// <summary>
        /// see MPI reference;
        /// </summary>
        /// <returns></returns>
        static public bool Initialized() {
            int flag = 0;
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Initialized");
#endif
            BoSSS_MPI_INITIALIZED(ref flag, out ierr);
            MPIException.CheckReturnCode(ierr);
            return (flag != 0);
        }


        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_WAITANY(ref int count, [In, Out] MPI_Request[] array_of_requests, out int index, out MPI_Status status, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="count"></param>
        /// <param name="array_of_requests"></param>
        /// <param name="index"></param>
        /// <param name="status"></param>
        static public void Waitany(int count, MPI_Request[] array_of_requests, out int index, out MPI_Status status) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Waitany");
#endif
            BoSSS_MPI_WAITANY(ref count, array_of_requests, out index, out status, out ierr);
            if (index != MiscConstants.UNDEFINED)
                index--; // convert fortran index into C-index
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_WAITALL(ref int count, [In, Out] MPI_Request[] array_of_requests, [Out] MPI_Status[] array_of_statuses, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="count"></param>
        /// <param name="array_of_requests"></param>
        /// <param name="array_of_statuses"></param>
        static public void Waitall(int count, MPI_Request[] array_of_requests, MPI_Status[] array_of_statuses) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Waitall");
#endif
            BoSSS_MPI_WAITALL(ref count, array_of_requests, array_of_statuses, out ierr);
            MPIException.CheckReturnCode(ierr);
        }


        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_ALLGATHER(IntPtr sendbuf, ref int sendcount, ref MPI_Datatype sendtype,
                                     IntPtr recvbuf, ref int recvcount, ref MPI_Datatype recvtype,
                                     ref MPI_Comm comm, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Allgather(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype,
                                     IntPtr recvbuf, int recvcount, MPI_Datatype recvtype,
                                     MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Allgather");
#endif
            BoSSS_MPI_ALLGATHER(sendbuf, ref sendcount, ref sendtype, recvbuf, ref recvcount, ref recvtype, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }



        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_ALLGATHERV(IntPtr sendbuf, ref int sendcount, ref MPI_Datatype sendtype,
                                      IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, ref MPI_Datatype recvtype,
                                      ref MPI_Comm comm, out int ierr);

        /// <summary>
        /// 
        /// </summary>
        static public void Allgatherv(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype,
                                      IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, MPI_Datatype recvtype,
                                      MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Allgatherv");
#endif
            BoSSS_MPI_ALLGATHERV(sendbuf, ref sendcount, ref sendtype, recvbuf, recvcounts, displs, ref recvtype, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }




        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_IRECV(IntPtr buf, ref int count, ref MPI_Datatype datatype,
                                 ref int source, ref int tag,
                                 ref MPI_Comm comm,
                                 out MPI_Request request,
                                 out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Irecv(IntPtr buf, int count, MPI_Datatype datatype,
                                 int source, int tag,
                                 MPI_Comm comm,
                                 out MPI_Request request) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Irecv");
#endif
            BoSSS_MPI_IRECV(buf, ref count, ref datatype, ref source, ref tag, ref comm, out request, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_BARRIER(ref MPI_Comm comm, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Barrier(MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            //Console.WriteLine("MPI_Barrier");
#endif
            BoSSS_MPI_BARRIER(ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_ISSEND(IntPtr buf, ref int count, ref MPI_Datatype datatype,
                                  ref int dest, ref int tag,
                                  ref MPI_Comm comm,
                                  out MPI_Request request, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Issend(IntPtr buf, int count, MPI_Datatype datatype,
                                  int dest, int tag,
                                  MPI_Comm comm,
                                  out MPI_Request request) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Issend");
#endif
            BoSSS_MPI_ISSEND(buf, ref count, ref datatype, ref dest, ref tag, ref comm, out request, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_BCAST(IntPtr buf, ref int count, ref MPI_Datatype datatype, ref int root, ref MPI_Comm comm, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Bcast(IntPtr buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Bcast");
#endif
            BoSSS_MPI_BCAST(buf, ref count, ref datatype, ref root, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_SEND(IntPtr buf, ref int count, ref MPI_Datatype datatype, ref int dest, ref int tag, ref MPI_Comm comm, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Send(IntPtr buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Send");
#endif
            BoSSS_MPI_SEND(buf, ref count, ref datatype, ref dest, ref tag, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_RECV(IntPtr buf, ref int count, ref MPI_Datatype datatype, ref int source, ref int tag, ref MPI_Comm comm, out MPI_Status status, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Recv(IntPtr buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, out MPI_Status status) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Recv");
#endif
            BoSSS_MPI_RECV(buf, ref count, ref datatype, ref source, ref tag, ref comm, out status, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        ///// <summary>
        ///// copys <paramref name="size"/> bytes form <paramref name="src"/> to <paramref name="dest"/>;
        ///// </summary>
        ///// <param name="dest"></param>
        ///// <param name="src"></param>
        ///// <param name="size"></param>
        //[DllImport("Platform_Native", EntryPoint = "UnmanagedMemcopy")]
        //static public extern void Memcpy(IntPtr dest, IntPtr src, int size);

        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_REDUCE(IntPtr sndbuf, IntPtr rcvbuf, ref int count, ref MPI_Datatype datatype, ref MPI_Op op, ref int root, ref MPI_Comm comm, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static public void Reduce(IntPtr sndbuf, IntPtr rcvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Reduce");
#endif
            BoSSS_MPI_REDUCE(sndbuf, rcvbuf, ref count, ref datatype, ref op, ref root, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }


        [DllImport("Platform_Native")]
        static extern void BoSSS_MPI_ALLREDUCE(IntPtr sendbuf, IntPtr recvbuf, ref int count, ref MPI_Datatype MPI_datatype, ref MPI_Op MPI_Op, ref MPI_Comm MPI_Comm, out int ierr);


        /// <summary>
        /// 
        /// </summary>
        static unsafe public void Allreduce(IntPtr sndbuf, IntPtr rcvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
            int ierr;
#if MPI_TRACE
            Console.WriteLine("MPI_Allreduce");
#endif
            BoSSS_MPI_ALLREDUCE(sndbuf, rcvbuf, ref count, ref datatype, ref op, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        [DllImport("Platform_Native")]
        static extern unsafe int BoSSS_MPI_Status_C2f(MPI_Status* cSt, int[] fSt);
        [DllImport("Platform_Native")]
        static extern unsafe int BoSSS_MPI_Status_f2c(int[] fSt, MPI_Status* cSt);
        [DllImport("Platform_Native")]
        static extern unsafe int BoSSS_MPI_Request_C2f(MPI_Request cRq, int* fRq);
        [DllImport("Platform_Native")]
        static extern unsafe int BoSSS_MPI_Request_f2c(int* fRq, MPI_Request* cRq);

        [DllImport("Platform_Native")]
        static extern int BoSSS_Get_MPI_Status_Size();
        
    }

    */
}
