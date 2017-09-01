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
using log4net;

namespace MPI.Wrappers {

    /// <summary>
    /// MPI Exceptions
    /// </summary>
    public class MPIException : ApplicationException {

        /*
        /// <summary>
        /// return codes from the MPI header files
        /// </summary>
        public enum ReturnCodes {

            // MPI's error classes 
            /// <summary> Successful return code </summary>
            MPI_SUCCESS=0,     
            // Communication argument parameters 
            /// <summary> Invalid buffer pointer </summary>
            MPI_ERR_BUFFER=1,  
            /// <summary> Invalid count argument </summary>
            MPI_ERR_COUNT=2,  
            /// <summary> Invalid datatype argument </summary>
            MPI_ERR_TYPE=3,   
            /// <summary> Invalid tag argument </summary>
            MPI_ERR_TAG=4,  
            /// <summary> Invalid communicator </summary>
            MPI_ERR_COMM=5,  
            /// <summary> Invalid rank </summary>
            MPI_ERR_RANK=6,  
            /// <summary> Invalid root  </summary>
            MPI_ERR_ROOT=7,   
            /// <summary> Message truncated on receive </summary>
            MPI_ERR_TRUNCATE=14,  
            // MPI Objects (other than COMM) 
            /// <summary> Invalid group </summary>
            MPI_ERR_GROUP=8,  
            /// <summary> Invalid operation  </summary>
            MPI_ERR_OP=9,  
            /// <summary> Invalid mpi_request handle </summary>
            MPI_ERR_REQUEST=19,  
            / Special topology argument parameters 
            /// <summary> Invalid topology </summary>
            MPI_ERR_TOPOLOGY=10, 
            /// <summary> Invalid dimension argument </summary>
            MPI_ERR_DIMS=11,  
            // All other arguments.  This is a class with many kinds 
            /// <summary> Invalid argument </summary>
            MPI_ERR_ARG=12, 
            // Other errors that are not simply an invalid argument 
            /// <summary> Other error; use Error_string </summary>
            MPI_ERR_OTHER=15, 
            /// <summary> Unknown error  </summary>
            MPI_ERR_UNKNOWN=13, 
            /// <summary> Internal error code </summary>
            MPI_ERR_INTERN=16, 
            // Multiple completion has two special error classes 
            /// <summary> Look in status for error value </summary>
            MPI_ERR_IN_STATUS=17, 
            /// <summary> Pending request </summary>
            MPI_ERR_PENDING=18,    
            // New MPI-2 Error classes 
            /// <summary>  </summary>
            MPI_ERR_FILE=27, 
            /// <summary>  </summary>
            MPI_ERR_ACCESS=20, 
            /// <summary>  </summary>
            MPI_ERR_AMODE=21,   
            /// <summary>  </summary>
            MPI_ERR_BAD_FILE=22,  
            /// <summary>  </summary>
            MPI_ERR_FILE_EXISTS=25, 
            /// <summary>  </summary>
            MPI_ERR_FILE_IN_USE=26,  
            /// <summary>  </summary>
            MPI_ERR_NO_SPACE=36,  
            /// <summary>  </summary>
            MPI_ERR_NO_SUCH_FILE=37,  
            /// <summary>  </summary>
            MPI_ERR_IO=32,  
            /// <summary>  </summary>
            MPI_ERR_READ_ONLY=40,  
            /// <summary>  </summary>
            MPI_ERR_CONVERSION=23,   
            /// <summary>  </summary>
            MPI_ERR_DUP_DATAREP=24,   
            /// <summary>  </summary>
            MPI_ERR_UNSUPPORTED_DATAREP=43, 

            // MPI_ERR_INFO is NOT defined in the MPI-2 standard.  I believe that this is an oversight 
            /// <summary>  </summary>
            MPI_ERR_INFO=28,  
            /// <summary>  </summary>
            MPI_ERR_INFO_KEY=29,  
            /// <summary>  </summary>
            MPI_ERR_INFO_VALUE=30, 
            /// <summary>  </summary>
            MPI_ERR_INFO_NOKEY=31,  
            /// <summary>  </summary>
            MPI_ERR_NAME=33,    
            /// <summary> Alloc_mem could not allocate memory </summary>
            MPI_ERR_NO_MEM=34,  
            /// <summary>  </summary>
            MPI_ERR_NOT_SAME=35,  
            /// <summary>  </summary>
            MPI_ERR_PORT=38,  
            /// <summary>  </summary>
            MPI_ERR_QUOTA=39,  
            /// <summary>  </summary>
            MPI_ERR_SERVICE=41,  
            /// <summary>  </summary>
            MPI_ERR_SPAWN=42,  
            /// <summary>  </summary>
            MPI_ERR_UNSUPPORTED_OPERATION=44, 
            /// <summary>  </summary>
            MPI_ERR_WIN=45,  
            /// <summary>  </summary>
            MPI_ERR_BASE=46,  
            /// <summary>  </summary>
            MPI_ERR_LOCKTYPE=47, 

            /// <summary> Erroneous attribute key </summary>
            MPI_ERR_KEYVAL=48,  
            /// <summary>  </summary>
            MPI_ERR_RMA_CONFLICT=49, 
            /// <summary>  </summary>
            MPI_ERR_RMA_SYNC=50, 
            /// <summary>  </summary>
            MPI_ERR_SIZE=51,     
            /// <summary>  </summary>
            MPI_ERR_DISP=52,     
            /// <summary>  </summary>
            MPI_ERR_ASSERT=53,   

            /// <summary> Last valid error code for a predefined error class  </summary>
            MPI_ERR_LASTCODE=0x3fffffff,  
            /// <summary> It is also helpful to know the last valid class </summary>
            MPICH_ERR_LAST_CLASS=53   

        }
        */


        /// <summary>
        /// da logger
        /// </summary>
        static ILog logger = LogManager.GetLogger(typeof(MPIException));

        /// <summary>
        /// 
        /// </summary>
        /// <param name="errorCode"></param>
        public MPIException(int errorCode) {
            m_ErrorCode = errorCode;
            m_ErrStr = csMPI.Raw.Error_string(errorCode);
            logger.Error(m_ErrStr);
        }

        /// <summary>
        /// throws an <see>MPIException</see> if <paramref name="returnCode"/>
        /// differs from MPI_SUCCESS;
        /// </summary>
        /// <param name="returnCode">return value of the C-style MPI function</param>
        internal static void CheckReturnCode(int returnCode) {
            int ret = returnCode;
            if (ret != 0)
                throw new MPIException(ret);
        }



        /// <summary>
        /// MPI Error code
        /// </summary>
        private int m_ErrorCode;


        /// <summary>
        /// MPI return value Code
        /// </summary>
        public int ErrorCode {
            get {
                return m_ErrorCode;
            }
        }


        string m_ErrStr;

        /// <summary>
        /// MPI error string
        /// </summary>
        public string ErrStr {
            get {
                return m_ErrStr;
            }
        }



        /// <summary>
        /// the MPI error code as a string
        /// </summary>
        public override string Message {
            get {
                return "MPI error code: " + m_ErrorCode + "; '" + m_ErrStr + "';";
            }
        }



    }
}
