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


        /// <summary>
        /// the main purpose of this property is to guarantee that the assembly for all types are linked.
        /// </summary>
        static public Type[] ExplicitHooks {
            get {
                return new[] {
                    typeof(Foundation.Grid.Classic.GridData),
                    typeof(OpenFOAMGrid),
                    typeof(FixedOperators),
                    typeof(Initializer)
                };
            }

        }


        static bool mustFinalizeMPI;

        /// <summary>
        /// Load/lookup of native libraries
        /// </summary>
        [CodeGenExport]
        public void BoSSSInitialize() {

            // Save original stdout and stderr
            TextWriter originalOut = Console.Out;
            TextWriter originalError = Console.Error;

            try {
                // Mute stdout and stderr
                using (var muteWriter = new StringWriter()) {
                    Console.SetOut(muteWriter);
                    Console.SetError(muteWriter);

                    // Call BoSSS init but we don't want to see output 
                    mustFinalizeMPI |= BoSSS.Solution.Application.InitMPI();
                }
            } catch (Exception ex) {
                Console.WriteLine(ex);
                throw;
            } finally {
                // Restore original stdout and stderr
                Console.SetOut(originalOut);
                Console.SetError(originalError);
            }
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

