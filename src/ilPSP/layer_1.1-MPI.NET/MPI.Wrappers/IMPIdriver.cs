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

namespace MPI.Wrappers {

    /// <summary>
    /// Interface to the subset of MPI that is used in BoSSS.
    /// </summary>
    /// <remarks>
    /// Documentation partly copied from the official MPI 2.2 standard, see
    /// www.mpi-forum.org
    /// </remarks>
    public interface IMPIdriver {

        /// <summary>
        /// Default communicators
        /// </summary>
        IMPI_CommConstants _COMM {
            get;
        }

        /// <summary>
        /// Predefined data-types
        /// </summary>
        IMPI_DatatypeConstants _DATATYPE {
            get;
        }

        /// <summary>
        /// Predefined operations
        /// </summary>
        IMPI_OpConstants _OP {
            get;
        }

        /// <summary>
        /// Combining data on one processor
        /// </summary>
        void Gatherv(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype,
                     IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, MPI_Datatype recvtype,
                     int root, MPI_Comm comm);
        
        /// <summary>
        /// Sending data from one processor to many
        /// </summary>
        void Scatterv(IntPtr sendbuf, IntPtr sendcounts, IntPtr displs, MPI_Datatype sendtype,
                      IntPtr recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm);

        /// <summary>
        /// MPI_ALLGATHER can be thought of as MPI_GATHER, but where all
        /// processes receive the result, instead of just the root. The block
        /// of data sent from the j-th process is received by every process and
        /// placed in the j-th block of the buffer <paramref name="recvbuf"/>.
        /// </summary>
        void Allgather(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype, IntPtr recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);

        /// <summary>
        /// MPI_ALLGATHERV can be thought of as MPI_GATHERV, but where all
        /// processes receive the result, instead of just the root. The block of
        /// data sent from the j-th process is received by every process and
        /// placed in the j-th block of the buffer <paramref name="recvbuf"/>.
        /// These blocks need not all be the same size.
        /// </summary>
        void Allgatherv(IntPtr sendbuf, int sendcount, MPI_Datatype sendtype, IntPtr recvbuf, IntPtr recvcounts, IntPtr displs, MPI_Datatype recvtype, MPI_Comm comm);

        /// <summary>
        /// MPI_ALLREDUCE behaves the same as <see cref="Reduce"/> except that
        /// the result appears in the receive buffer of all the group members
        /// </summary>
        void Allreduce(IntPtr sndbuf, IntPtr rcvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

        /// <summary>
        /// If comm is an intracommunicator, MPI_BARRIER blocks the caller
        /// until all group members have called it. The call returns at any
        /// process only after all group members have entered the call.
        /// </summary>
        /// <param name="comm"></param>
        void Barrier(MPI_Comm comm);

        /// <summary>
        /// If comm is an intracommunicator, MPI_BCAST broadcasts a message
        /// from the process with rank root to all processes of the group,
        /// itself included. It is called by all members of the group using the
        /// same arguments for comm and root. On return, the content of root's
        /// buffer is copied to all other processes.
        /// </summary>
        void Bcast(IntPtr buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

        /// <summary>
        /// This function gives the rank of the process in the particular
        /// communicator's group. It is useful, as noted above, in conjunction
        /// with <see cref="Comm_Size"/>.
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="rank"></param>
        /// <remarks>
        /// This operation is purely local and involves no communication
        /// </remarks>
        void Comm_Rank(MPI_Comm comm, out int rank);

        /// <summary>
        /// This function indicates the number of processes involved in a
        /// communicator. For <see cref="IMPI_CommConstants.WORLD"/>, it
        /// indicates the total number of processes available (for this version
        /// of MPI, there is no standard way to change the number of processes
        /// once initialization has taken place)
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="size"></param>
        /// <remarks>
        /// This operation is purely local and involves no communication
        /// </remarks>
        void Comm_Size(MPI_Comm comm, out int size);

        /// <summary>
        /// Returns the error string associated with an error code or class.
        /// </summary>
        /// <param name="errorCode"></param>
        /// <returns></returns>
        string Error_string(int errorCode);

        /// <summary>
        /// The supported MPI implementation
        /// </summary>
        MPI_Flavour Flavour {
            get;
        }

        /// <summary>
        /// All MPI programs must contain exactly one call to an MPI
        /// initialization routine. Subsequent calls to any initialization
        /// routines are erroneous. The only MPI functions that may be invoked
        /// before the MPI initialization are MPI_GET_VERSION,
        /// <see cref="Initialized"/>, and MPI_FINALIZED.
        /// </summary>
        /// <param name="args"></param>
        void Init(string[] args);

        /// <summary>
        /// This routine may be used to determine whether <see cref="Init"/>
        /// has been called. MPI_INITIALIZED returns true if the calling
        /// process has called MPI_INIT. Whether <see cref="mpiFinalize"/> has
        /// been called does not affect the behavior of MPI_INITIALIZED. It is
        /// one of the few routines that may be called before
        /// <see cref="Init"/> is called.
        /// </summary>
        /// <returns></returns>
        bool Initialized();

        /// <summary>
        /// Start a nonblocking receive. These calls allocate a communication
        /// request object and associate it with the request handle (the
        /// argument <paramref name="request"/>). The request can be used later
        /// to query the status of the communication or wait for its
        /// completion.
        /// </summary>
        /// <param name="buf"></param>
        /// <param name="count"></param>
        /// <param name="datatype"></param>
        /// <param name="source"></param>
        /// <param name="tag"></param>
        /// <param name="comm"></param>
        /// <param name="request"></param>
        /// <remarks>
        /// A nonblocking receive call indicates that the system may start
        /// writing data into the receive buffer. The receiver should not
        /// access any part of the receive buffer after a nonblocking receive
        /// operation is called, until the receive completes.
        /// </remarks>
        void Irecv(IntPtr buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, out MPI_Request request);

        /// <summary>
        /// Start a synchronous mode, nonblocking send.
        /// </summary>
        /// <param name="buf"></param>
        /// <param name="count"></param>
        /// <param name="datatype"></param>
        /// <param name="dest"></param>
        /// <param name="tag"></param>
        /// <param name="comm"></param>
        /// <param name="request"></param>
        /// <remarks>
        /// A nonblocking send call indicates that the system may start copying
        /// data out of the send buffer.
        /// </remarks>
        void Issend(IntPtr buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, out MPI_Request request);

        /// <summary>
        /// Misc MPI constants
        /// </summary>
        IMiscConstants MiscConstants {
            get;
        }

        /// <summary>
        /// Converts the MPI communicators in ilPSP (which are
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
        int MPI_Comm_f2c(MPI_Comm input, out uint FourByteComm, out ulong EightByteComm);

        /// <summary>
        /// This routine cleans up all MPI state. Each process must call
        /// MPI_FINALIZE before it exits. Unless there has been a call to
        /// MPI_ABORT, each process must ensure that all pending nonblocking
        /// communications are (locally) complete before calling MPI_FINALIZE.
        /// Further, at the instant at which the last process calls
        /// MPI_FINALIZE, all pending sends must be matched by a receive, and
        /// all pending receives must be matched by a send.
        /// </summary>
        void mpiFinalize();

        /// <summary>
        /// Starts a blocking receive.
        /// </summary>
        /// <param name="buf"></param>
        /// <param name="count"></param>
        /// <param name="datatype"></param>
        /// <param name="source"></param>
        /// <param name="tag"></param>
        /// <param name="comm"></param>
        /// <param name="status"></param>
        void Recv(IntPtr buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, out MPI_Status status);

        /// <summary>
        /// If <paramref name="comm"/> is an intracommunicator, MPI_REDUCE
        /// combines the elements provided in the input buffer of each process
        /// in the group, using the operation <paramref name="op"/>, and
        /// returns the combined value in the output buffer of the process with
        /// rank root. The input buffer is defined by the arguments
        /// <paramref name="sndbuf"/>, <paramref name="count"/> and
        /// <paramref name="datatype"/>; the output buffer is defined by the
        /// arguments <paramref name="rcvbuf"/>, <paramref name="count"/> and
        /// <paramref name="datatype"/>; both have the same number of elements,
        /// with the same type.
        /// The routine is called by all group members using the same arguments
        /// for <paramref name="count"/>, <paramref name="datatype"/>,
        /// <paramref name="op"/>, <paramref name="root"/> and
        /// <paramref name="comm"/>. Thus, all processes provide input buffers
        /// and output buffers of the same length, with elements of the same
        /// type. Each process can provide one element, or a sequence of
        /// elements, in which case the combine operation is executed
        /// element-wise on each entry of the sequence.
        /// </summary>
        /// <param name="sndbuf"></param>
        /// <param name="rcvbuf"></param>
        /// <param name="count"></param>
        /// <param name="datatype"></param>
        /// <param name="op"></param>
        /// <param name="root"></param>
        /// <param name="comm"></param>
        void Reduce(IntPtr sndbuf, IntPtr rcvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

        /// <summary>
        /// Starts a blocking send.
        /// </summary>
        /// <param name="buf"></param>
        /// <param name="count"></param>
        /// <param name="datatype"></param>
        /// <param name="dest"></param>
        /// <param name="tag"></param>
        /// <param name="comm"></param>
        void Send(IntPtr buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

        /// <summary>
        /// Blocks until all communication operations associated with active
        /// handles in the list complete, and return the status of all these
        /// operations (this includes the case where no handle in the list is
        /// active). Both arrays have the same number of valid entries. The
        /// i-th entry in <paramref name="array_of_statuses"/> is set to the
        /// return status of the i-th operation. Requests that were created by
        /// nonblocking communication operations are deallocated and the
        /// corresponding handles in the array are set to MPI_REQUEST_NULL. The
        /// list may contain null or inactive handles. The call sets to empty
        /// the status of each such entry.
        /// </summary>
        /// <param name="count"></param>
        /// <param name="array_of_requests"></param>
        /// <param name="array_of_statuses"></param>
        void Waitall(int count, MPI_Request[] array_of_requests, MPI_Status[] array_of_statuses);

        /// <summary>
        /// Blocks until one of the operations associated with the active
        /// requests in the array has completed. If more then one operation is
        /// enabled and can terminate, one is arbitrarily chosen. Returns in
        /// index the index of that request in the array and returns in status
        /// the status of the completing communication. (The array is indexed
        /// from zero in C, and from one in Fortran.) If the request was
        /// allocated by a nonblocking communication operation, then it is
        /// deallocated and the request handle is set to MPI_REQUEST_NULL.
        /// </summary>
        /// <param name="count"></param>
        /// <param name="array_of_requests"></param>
        /// <param name="index"></param>
        /// <param name="status"></param>
        void Waitany(int count, MPI_Request[] array_of_requests, out int index, out MPI_Status status);

        /// <summary>
        /// Return the parent communicator for this process.
        /// </summary>
        void Comm_get_parent(out MPI_Comm parent);

        /// <summary>
        /// The size of an <see cref="MPI_Status"/>.
        /// </summary>
        int MPI_Status_Size {
            get;
        }
    }
}

