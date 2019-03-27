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
using System.Runtime.InteropServices;

namespace MPI.Wrappers {

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
    /// This is a bit difficult, because this implementations are not binary
    /// equivalent; E.g. MPICH defines the MPI_Comm data-type as int, which is
    /// always 4 Byte, independent of whether the application runs in 32- or
    /// 64-bit mode; OpenMPI, on the other hand, defines MPI_Comm as a pointer,
    /// so it depends on the 'bitness' of the application/platform;
    /// We solved such issues by addressing the FORTRAN interface, where types
    /// like MPI_Comm, MPI_Op, ... are defined to be integers by the standard;
    /// </remarks>
    class FortranMPIdriver : Utils.DynLibLoader, IMPIdriver {

        /// <summary>
        /// converts the MPI communicators in ilPSP (which are
        /// FORTRAN-communicators) into a C-Communicator
        /// </summary>
        /// <param name="EightByteComm">
        /// on exit, the value of the C - MPI communicator handle, if
        /// <see cref="IMPI_CommConstants.GetSizeof_C_MPI_comm"/> == 8,
        /// otherwise 0
        /// </param>
        /// <param name="FourByteComm">
        /// on exit, the value of the C - MPI communicator handle, if
        /// <see cref="IMPI_CommConstants.GetSizeof_C_MPI_comm"/> == 4,
        /// otherwise 0;
        /// </param>
        /// <param name="input">
        /// (FORTRAN) MPI communicator handle to convert
        /// </param>
        /// <returns>
        /// Equal to <see cref="IMPI_CommConstants.GetSizeof_C_MPI_comm"/>, 
        /// either 4 or 8
        /// </returns>
        /// <remarks>
        /// The length of the C - MPI communicator handle depends on the used
        /// MPI implementation (see remarks at <see cref="MPI_Comm"/>) and is
        /// given by to the return value of
        /// <see cref="IMPI_CommConstants.GetSizeof_C_MPI_comm"/>
        /// </remarks>
        public int MPI_Comm_f2c(MPI_Comm input, out uint FourByteComm, out ulong EightByteComm) {
            unsafe {
                byte[] C_Comm = _COMM.Comm_f2c(input);
                fixed (byte* pComm = &(C_Comm[0])) {
                    uint _FourByteComm = 0;
                    ulong _EightByteComm = 0;
                    byte* pDest;
                    switch (C_Comm.Length) {
                        case 4:
                            pDest = (byte*)(&_FourByteComm);
                            break;
                        case 8:
                            pDest = (byte*)(&_EightByteComm);
                            break;
                        default:
                            throw new NotImplementedException(
                                "unknown size of C MPI communicator: " + C_Comm.Length + "bytes; ");
                    }
                    for (int i = 0; i < C_Comm.Length; i++)
                        pDest[i] = pComm[i];

                    FourByteComm = _FourByteComm;
                    EightByteComm = _EightByteComm;
                }
                return C_Comm.Length;
            }
        }

        /// <summary>
        /// workaround for .NET bug:
        /// https://connect.microsoft.com/VisualStudio/feedback/details/635365/runtimehelpers-initializearray-fails-on-64b-framework
        /// </summary>
        /// <returns></returns>
        static PlatformID[] Helper() {
            PlatformID[] p = new PlatformID[6];
            p[0] = PlatformID.Win32NT;
            p[1] = PlatformID.Unix;
			p[2] = PlatformID.Unix;
			p[3] = PlatformID.Unix;
            p[4] = PlatformID.Unix;
            p[5] = PlatformID.MacOSX;
            return p;
        }

        /// <summary>
        /// Strange name mangling convention on MacOS; seems to suit OpenMPI from homebrew.
        /// </summary>
        internal static string MacOsMangling (string Name) {
            if(Name.ToUpperInvariant () == Name) {
                // one of the Fortran guys
                return Name.ToLowerInvariant () + "_";
            } else {
                return Name;
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        internal FortranMPIdriver()
            : base(
				new string[] { "msmpi.dll", "libmpi_f77.so", "libfmpich.so", "libmpi_mpifh.so", "openmpi/libmpi_mpifh.so", "/usr/local/opt/open-mpi/lib/libmpi_mpifh.dylib" },
                new string[6][][], 
				new GetNameMangling[] { Utils.DynLibLoader.Identity, Utils.DynLibLoader.Identity, Utils.DynLibLoader.Identity, Utils.DynLibLoader.Identity, Utils.DynLibLoader.Identity, MacOsMangling  },
                //new PlatformID[] {PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix },
                Helper(),
                new int[] { -1, -1, -1, -1, -1, -1 }) {

            // depending on the MPI flavor, we define the Datatype
            // ---------------------------------------------------

            if (base.CurrentLibraryName.ToLowerInvariant().Contains("msmpi")
                || base.CurrentLibraryName.ToLowerInvariant().Contains("mpich")) {
                m_MPI_Comm = new MPICH_MPI_Comm();
                m_MiscConstants = new MPICH_MiscConstants();
                m_Flavour = MPI_Flavour.MPICH;
                m_MPI_OP = new MPICH_MPI_op();
                m_MPI_Datatype = new MPICH_MPI_Datatype();
            } else if (base.CurrentLibraryName.ToLowerInvariant().Contains("libmpi")) {
                OPENMPI_Converter conv = new OPENMPI_Converter();
                m_MPI_Comm = new OpenMPI_MPI_Comm(conv);
				if (base.CurrentLibraryName.ToLowerInvariant().Contains("mpi_f77"))
					dirrtyMpiInitHack = 1;
				else if (base.CurrentLibraryName.ToLowerInvariant().Contains("mpi_mpifh"))
					dirrtyMpiInitHack = 2; 
                //throw new NotImplementedException();
                m_MiscConstants = new OPENMPI_MiscConstants(conv);
                m_Flavour = MPI_Flavour.OpenMPI;
                m_MPI_OP = new OpenMPI_MPI_Op(conv);
                m_MPI_Datatype = new OpenMPI_MPI_Datatype(conv);
            } else {
                throw new NotImplementedException("Unknown MPI;");
            }
        }

        MPI_Flavour m_Flavour;

        /// <summary>
        /// identifies the currently loaded MPI implementation
        /// </summary>
        public MPI_Flavour Flavour {
            get {
                return m_Flavour;
            }
        }

        IMPI_CommConstants m_MPI_Comm;

        /// <summary>
        /// default communicators
        /// </summary>
        public IMPI_CommConstants _COMM {
            get {
                return m_MPI_Comm;
            }
        }

        IMiscConstants m_MiscConstants;

        /// <summary>
        /// misc MPI - constants
        /// </summary>
        public IMiscConstants MiscConstants {
            get {
                return m_MiscConstants;
            }
        }

        /// <summary>
        /// predefined MPI operations
        /// </summary>
        public IMPI_OpConstants _OP {
            get {
                return m_MPI_OP;
            }
        }

        IMPI_OpConstants m_MPI_OP;

        IMPI_DatatypeConstants m_MPI_Datatype;

        /// <summary>
        /// predefined MPI operations
        /// </summary>
        public IMPI_DatatypeConstants _DATATYPE {
            get {
                return m_MPI_Datatype;
            }
        }

        ///// <summary>
        ///// the only instance of <see cref="csMPI.Raw"/>
        ///// </summary>
        //static csMPI.Raw instance = new csMPI.Raw();

#pragma warning disable 649
        unsafe delegate int _MPI_Error_string(int errorcode, char* str, out int resultlen);
        _MPI_Error_string MPI_Error_string;
#pragma warning restore 649

        /// <summary>
        /// gets an error string for an MPI error code
        /// </summary>
        /// <param name="errorCode">
        /// MPI error code (see <see cref="MPIException.ErrorCode"/>
        /// </param>
        /// <returns>
        /// the description of the error no
        /// <see cref="MPIException.ErrorCode"/> from the MPI library
        /// </returns>
        public string Error_string(int errorCode) {
            unsafe {
                char* str = (char*)Marshal.AllocHGlobal(65000);
                int resultlen;
                MPI_Error_string(errorCode, str, out resultlen);
                string r = Marshal.PtrToStringAnsi((IntPtr)str);
                Marshal.FreeHGlobal((IntPtr)str);
                return r;
                //return "MPI_ERROR_String not loaded.";
            }
        }

        #region INIT

        /// <summary>
		/// if unequal 0, <see cref="OpenMPI_MPI_Init1"/> or <see cref="OpenMPI_MPI_Init2"/> is taken instead 
        /// of <see cref="MPI_Init"/>;
        /// For some strange reason that i cannot reproduce, MPI_Init does not
        /// work when loaded via dlsym(...);
        /// </summary>
        int dirrtyMpiInitHack = 0;

        /// <summary>
        /// see <see cref="dirrtyMpiInitHack"/>;
        /// </summary>
		[DllImport("mpi_f77", EntryPoint = "MPI_Init")]
		static extern private unsafe int OpenMPI_MPI_Init1(int* argc, byte*** argv);

		/// <summary>
		/// see <see cref="dirrtyMpiInitHack"/>;
		/// </summary>
		[DllImport("mpi_mpifh", EntryPoint = "MPI_Init")]
		static extern private unsafe int OpenMPI_MPI_Init2(int* argc, byte*** argv);



#pragma warning disable 649
        unsafe delegate int _MPI_Init(int* argc, byte*** argv);
        _MPI_Init MPI_Init;
#pragma warning restore 649

        /// <summary>
        /// see MPI reference;
        /// </summary>
        /// <param name="args">
        /// command line args passed to the application
        /// </param>
        /// <remarks>
        /// Do not call this directly, use ilPSP.Enviroment.Bootstrap instead.
        /// </remarks>
        public void Init(string[] args) {
            IntPtr[] argsAnsi = new IntPtr[Math.Max(args.Length, 1)];

            for (int i = 0; i < args.Length; i++) {
                argsAnsi[i] = Marshal.StringToCoTaskMemAnsi(args[i]);
            }

            unsafe {
                fixed (void* argv = &(argsAnsi[0])) {
                    void* _argv = argv;
                    int argc = args.Length;
					if (dirrtyMpiInitHack == 1)
						MPIException.CheckReturnCode(OpenMPI_MPI_Init1(&argc, (byte***)(&_argv)));
					else if (dirrtyMpiInitHack == 2)
					    MPIException.CheckReturnCode(OpenMPI_MPI_Init2(&argc, (byte***)(&_argv)));
					else
                        MPIException.CheckReturnCode(MPI_Init(&argc, (byte***)(&_argv)));
                }
            }

            foreach (IntPtr p in argsAnsi) {
                Marshal.FreeCoTaskMem(p);
            }
        }

#pragma warning disable 649
        unsafe delegate int _MPI_Finalize();
        _MPI_Finalize MPI_Finalize;
#pragma warning restore 649

        /// <summary>
        /// MPI finalize
        /// </summary>
        public void mpiFinalize() {
            MPIException.CheckReturnCode(MPI_Finalize());
        }

        #endregion


#pragma warning disable 649
        unsafe delegate void _MPI_COMM_RANK(ref MPI_Comm comm, out int rank, out int ierr);
        _MPI_COMM_RANK MPI_COMM_RANK;
#pragma warning restore 649


        #region Group destructors
        /// <summary>
        /// 
        /// </summary>
        public void Comm_Rank(MPI_Comm comm, out int rank) {
            int ierr;
            MPI_COMM_RANK(ref comm, out rank, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        unsafe delegate void _MPI_COMM_SIZE(ref MPI_Comm comm, out int size, out int ierr);
        _MPI_COMM_SIZE MPI_COMM_SIZE;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Comm_Size(MPI_Comm comm, out int size) {
            int ierr;
            MPI_COMM_SIZE(ref comm, out size, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

        #endregion

#pragma warning disable 649
        delegate void _MPI_INITIALIZED(ref int flag, out int ierr);
        _MPI_INITIALIZED MPI_INITIALIZED;
#pragma warning restore 649

        /// <summary>
        /// see MPI reference;
        /// </summary>
        /// <returns></returns>
        public bool Initialized() {
            int flag = 0;
            int ierr;
            MPI_INITIALIZED(ref flag, out ierr);
            MPIException.CheckReturnCode(ierr);
            return (flag != 0);
        }


#pragma warning disable 649
        delegate void _MPI_WAITANY(ref int count, [In, Out] MPI_Request[] array_of_requests, out int index, out MPI_Status status, out int ierr);
        _MPI_WAITANY MPI_WAITANY;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        /// <param name="count"></param>
        /// <param name="array_of_requests"></param>
        /// <param name="index"></param>
        /// <param name="status"></param>
        public void Waitany(int count, MPI_Request[] array_of_requests, out int index, out MPI_Status status) {
            int ierr;

            // note: since the status is passed by reference/pointer, 
            // some larger buffer does not do any harm 
            // (if our internal MPI_Status structure is larger than the one actually defined by the MPI-implementation)

            MPI_WAITANY(ref count, array_of_requests, out index, out status, out ierr);
            if (index != MiscConstants.UNDEFINED)
                index--; // convert fortran index into C-index
            MPIException.CheckReturnCode(ierr);
            FixMPIStatus(ref status);
        }

#pragma warning disable 649
        delegate void _MPI_WAITALL(ref int count, [In, Out] MPI_Request[] array_of_requests, [In, Out] MPI_Status[] array_of_statii, out int ierr);
        _MPI_WAITALL MPI_WAITALL;
#pragma warning restore 649


        /// <summary>
        /// this method deals with different sizes of MPI_Status in BoSSS and the actual MPI implementation.
        /// </summary>
        /// <param name="array_of_statii"></param>
        void FixMPI_Status(MPI_Status[] array_of_statii) {
            unsafe {
                fixed (MPI_Status* pStatii = array_of_statii) {
                    int* pIStatii = (int*)pStatii;

                    int NatSz = this.MPI_Status_Size;
                    int BosSz = sizeof(MPI_Status) / sizeof(int);

                    if (BosSz == NatSz)
                        return;

                    int DifSz = BosSz - NatSz;
                    Debug.Assert(DifSz >= 0);

                    int NoOfStatii = array_of_statii.Length;

                    int* pTailNat = pIStatii + NoOfStatii * NatSz - 1;
                    int* pTailBos = pIStatii + NoOfStatii * BosSz - 1;

                    for (int i = 0; i < NoOfStatii; i++) {
                        for (int j = 0; j < DifSz; j++) {
                            *pTailBos = 0;
                            pTailBos--;
                        }

                        for (int j = 0; j < NatSz; j++) {
                            *pTailBos = *pTailNat;
                            pTailNat--;
                            pTailBos--;
                        }
                    }
                }
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="count"></param>
        /// <param name="array_of_requests"></param>
        /// <param name="array_of_statii"></param>
        public void Waitall(int count, MPI_Request[] array_of_requests, MPI_Status[] array_of_statii) {
            int ierr;
            MPI_WAITALL(ref count, array_of_requests, array_of_statii, out ierr);
            FixMPI_Status(array_of_statii);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_GATHERV(IntPtr sendbuf, ref int sendcount, ref MPI_Datatype sendtype,
                                   IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, ref MPI_Datatype recvtype,
                                   ref int root, ref MPI_Comm comm, out int ierr);
        _MPI_GATHERV MPI_GATHERV;
#pragma warning restore 649
        
        public void Gatherv(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype,
                            IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, MPI_Datatype recvtype,
                            int root, MPI_Comm comm) {
            MPI_GATHERV(sendbuf, ref sendcount, ref sendtype, recvbuf, recvcounts, displs, ref recvtype, ref root, ref comm, out int ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_SCATTERV(IntPtr sendbuf, IntPtr sendcounts, IntPtr displs, ref MPI_Datatype sendtype,
                                    IntPtr recvbuf, ref int recvcounts, ref MPI_Datatype recvtype,
                                    ref int root, ref MPI_Comm comm, out int ierr);
        _MPI_SCATTERV MPI_SCATTERV;
#pragma warning restore 649

        public void Scatterv(IntPtr sendbuf, IntPtr sendcounts, IntPtr displs, MPI_Datatype sendtype,
                             IntPtr recvbuf, int recvcount, MPI_Datatype recvtype,
                             int root, MPI_Comm comm) {
            MPI_SCATTERV(sendbuf, sendcounts, displs, ref sendtype, recvbuf, ref recvcount, ref recvtype, ref root, ref comm, out int ierr);
            MPIException.CheckReturnCode(ierr);
        }


#pragma warning disable 649
        delegate void _MPI_ALLGATHER(IntPtr sendbuf, ref int sendcount, ref MPI_Datatype sendtype,
                                     IntPtr recvbuf, ref int recvcount, ref MPI_Datatype recvtype,
                                     ref MPI_Comm comm, out int ierr);
        _MPI_ALLGATHER MPI_ALLGATHER;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Allgather(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype,
                                     IntPtr recvbuf, int recvcount, MPI_Datatype recvtype,
                                     MPI_Comm comm) {
            int ierr;
            MPI_ALLGATHER(sendbuf, ref sendcount, ref sendtype, recvbuf, ref recvcount, ref recvtype, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_GATHER(IntPtr sendbuf, ref int sendcount, ref MPI_Datatype sendtype,
                                     IntPtr recvbuf, ref int recvcount, ref MPI_Datatype recvtype, 
                                     int root, ref MPI_Comm comm, out int ierr);
        _MPI_GATHER MPI_GATHER;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Gather(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype, IntPtr recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
            int ierr;
            MPI_GATHER(sendbuf, ref sendcount, ref sendtype, recvbuf, ref recvcount, ref recvtype, root,ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }



#pragma warning disable 649
        delegate void _MPI_ALLGATHERV(IntPtr sendbuf, ref int sendcount, ref MPI_Datatype sendtype,
                                      IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, ref MPI_Datatype recvtype,
                                      ref MPI_Comm comm, out int ierr);
        _MPI_ALLGATHERV MPI_ALLGATHERV;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Allgatherv(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype,
                                      IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, MPI_Datatype recvtype,
                                      MPI_Comm comm) {
            int ierr;
            MPI_ALLGATHERV(sendbuf, ref sendcount, ref sendtype, recvbuf, recvcounts, displs, ref recvtype, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }




#pragma warning disable 649
        delegate void _MPI_IRECV(IntPtr buf, ref int count, ref MPI_Datatype datatype,
                                 ref int source, ref int tag,
                                 ref MPI_Comm comm,
                                 out MPI_Request request,
                                 out int ierr);
        _MPI_IRECV MPI_IRECV;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Irecv(IntPtr buf, int count, MPI_Datatype datatype,
                                 int source, int tag,
                                 MPI_Comm comm,
                                 out MPI_Request request) {
            int ierr;
            MPI_IRECV(buf, ref count, ref datatype, ref source, ref tag, ref comm, out request, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_BARRIER(ref MPI_Comm comm, out int ierr);
        _MPI_BARRIER MPI_BARRIER;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Barrier(MPI_Comm comm) {
            int ierr;
            MPI_BARRIER(ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_ISSEND(IntPtr buf, ref int count, ref MPI_Datatype datatype,
                                  ref int dest, ref int tag,
                                  ref MPI_Comm comm,
                                  out MPI_Request request, out int ierr);
        _MPI_ISSEND MPI_ISSEND;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Issend(IntPtr buf, int count, MPI_Datatype datatype,
                                  int dest, int tag,
                                  MPI_Comm comm,
                                  out MPI_Request request) {
            int ierr;
            MPI_ISSEND(buf, ref count, ref datatype, ref dest, ref tag, ref comm, out request, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_BCAST(IntPtr buf, ref int count, ref MPI_Datatype datatype, ref int root, ref MPI_Comm comm, out int ierr);
        _MPI_BCAST MPI_BCAST;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Bcast(IntPtr buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
            int ierr;
            MPI_BCAST(buf, ref count, ref datatype, ref root, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_SEND(IntPtr buf, ref int count, ref MPI_Datatype datatype, ref int dest, ref int tag, ref MPI_Comm comm, out int ierr);
        _MPI_SEND MPI_SEND;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Send(IntPtr buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
            int ierr;
            MPI_SEND(buf, ref count, ref datatype, ref dest, ref tag, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_RECV(IntPtr buf, ref int count, ref MPI_Datatype datatype, ref int source, ref int tag, ref MPI_Comm comm, out MPI_Status status, out int ierr);
        _MPI_RECV MPI_RECV;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Recv(IntPtr buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, out MPI_Status status) {
            int ierr;
            MPI_RECV(buf, ref count, ref datatype, ref source, ref tag, ref comm, out status, out ierr);
            FixMPIStatus(ref status);
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

#pragma warning disable 649
        delegate void _MPI_REDUCE(IntPtr sndbuf, IntPtr rcvbuf, ref int count, ref MPI_Datatype datatype, ref MPI_Op op, ref int root, ref MPI_Comm comm, out int ierr);
        _MPI_REDUCE MPI_REDUCE;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Reduce(IntPtr sndbuf, IntPtr rcvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
            int ierr;
            MPI_REDUCE(sndbuf, rcvbuf, ref count, ref datatype, ref op, ref root, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _MPI_ALLREDUCE(IntPtr sendbuf, IntPtr recvbuf, ref int count, ref MPI_Datatype MPI_datatype, ref MPI_Op MPI_Op, ref MPI_Comm MPI_Comm, out int ierr);
        _MPI_ALLREDUCE MPI_ALLREDUCE;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        unsafe public void Allreduce(IntPtr sndbuf, IntPtr rcvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
            int ierr;
            MPI_ALLREDUCE(sndbuf, rcvbuf, ref count, ref datatype, ref op, ref comm, out ierr);
            MPIException.CheckReturnCode(ierr);
        }

#pragma warning disable 649
        delegate void _COMM_GET_PARENT(out MPI_Comm parent, out int ierr);
        _COMM_GET_PARENT MPI_COMM_GET_PARENT;
#pragma warning restore 649

        /// <summary>
        /// 
        /// </summary>
        public void Comm_get_parent(out MPI_Comm parent) { 
            int ierr;
            MPI_COMM_GET_PARENT(out parent, out ierr);
            MPIException.CheckReturnCode(ierr);
        }


        int m_MPI_Status_Size = -1;

        /// <summary>
        /// the size of the native MPI_Status -- structure, in numbers of integers;
        /// Currently we know that this may be either 5 or 6 integers, depending on the MPI-implementation/version.
        /// </summary>
        public int MPI_Status_Size {
            get {
                if (m_MPI_Status_Size < 0) {
                    unsafe {
                        if (sizeof(MPI_Status) % sizeof(int) != 0)
                            throw new NotSupportedException("size of the MPI_Status must be a multiple of sizeof(int).");


                        int rcvBuf, sndBuf;
                        int count, rank, tag;
                        MPI_Comm comm;
                        int ierr;

                        comm = this._COMM.WORLD;
                        MPI_COMM_RANK(ref comm, out rank, out ierr);
                        sndBuf = 55;
                        tag = 444;

                        MPI_Datatype type = this._DATATYPE.INT;
                        MPI_Request[] Req = new MPI_Request[1];
                        MPI_Status[] Statussies = new MPI_Status[55];
                        count = 1;
                        MPI_IRECV((IntPtr)(&rcvBuf), ref count, ref type, ref rank, ref tag, ref comm, out Req[0], out ierr);
                        MPI_SEND((IntPtr)(&sndBuf), ref count, ref type, ref rank, ref tag, ref comm, out ierr);

                        fixed (MPI_Status* pStatussies = Statussies) {
                            int* pIStatussies = (int*)pStatussies;
                            int Buffersize = Statussies.Length * sizeof(MPI_Status) / sizeof(int);
                            for (int i = 0; i < Buffersize; i++) {
                                pIStatussies[i] = 666666;
                            }


                            MPI_WAITALL(ref count, Req, Statussies, out ierr);
                            m_MPI_Status_Size = 0;
                            for (int i = 0; i < Buffersize; i++) {
                                if (pIStatussies[i] != 666666) {
                                    m_MPI_Status_Size++;
                                } else {
                                    break;
                                }
                            }

                            if (m_MPI_Status_Size > sizeof(MPI_Status) / sizeof(int))
                                throw new NotSupportedException("size of native MPI status is not supported.");
                        }
                    }
                    //Console.WriteLine("MPI_Status - size: " + m_MPI_Status_Size);
                }
                return m_MPI_Status_Size;
            }
        }

        private void FixMPIStatus(ref MPI_Status st) {
            unsafe {
                Debug.Assert(sizeof(MPI_Status) <= 6 * sizeof(int), "implement higher schas");
                int NativeSize = this.MPI_Status_Size;
                if (NativeSize <= 5)
                    st.i6 = 0;
                else
                    return;

                if (NativeSize <= 4)
                    st.i5 = 0;
                else
                    return;

                if (NativeSize <= 3)
                    st.i4 = 0;
                else
                    return;

                if (NativeSize <= 2)
                    st.i3 = 0;
                else
                    return;

                if (NativeSize <= 1)
                    st.i2 = 0;
                else
                    return;
            }
        }
    }
}

