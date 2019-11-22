using System;
using System.Collections.Generic;
using ilPSP.LinSolvers;
using System.Runtime.InteropServices;
using System.Xml;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP;
using ilPSP.Connectors;

namespace BoSSS.Application.ExternalBinding {
    /// <summary>
    /// Initialization stuff
    /// </summary>
    public class Initializer : IForeignLanguageProxy {

        /// <summary>
        /// Constructor for export.
        /// </summary>
        [CodeGenExport]
        public Initializer() {

        }


        static bool mustFinalizeMPI;

        /// <summary>
        /// Load/lookup of native libraries
        /// </summary>
        [CodeGenExport]
        public void BoSSSInitialize() {
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out mustFinalizeMPI);
            //Console.WriteLine("elo from " +  Enviroment.MPIEnv.MPI_Rank + " of " + Enviroment.MPIEnv.MPI_Size);
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        [CodeGenExport]
        public void BoSSSFinalize() {
            if (mustFinalizeMPI)
                MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        IntPtr m_ForeignPtr;

        /// <summary>
        /// %
        /// </summary>
        public void _SetForeignPointer(IntPtr ptr) {
            if (ptr == IntPtr.Zero) {
                m_ForeignPtr = IntPtr.Zero;
            } else {

                if (m_ForeignPtr != IntPtr.Zero) {
                    throw new ApplicationException("already registered");
                }
                m_ForeignPtr = ptr;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IntPtr _GetForeignPointer() {
            return m_ForeignPtr;
        }
    }

}

