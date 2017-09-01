using System;

namespace bcl {

    /// <summary>
    /// Common interface for all subroutines of the BoSSS command line tools.
    /// </summary>
    interface IProgram {

        /// <summary>
        /// Starts the execution of the subroutine.
        /// </summary>
        void Execute();

        /// <summary>
        /// Reads a list of command line arguments specific to any class
        /// implementing this interface and stores the result so it can be used
        /// in <see cref="Execute"/>.
        /// </summary>
        /// <param name="args"></param>
        void DecodeArgs(string[] args);

        /// <summary>
        /// Provides information on how to use the subroutine for the user.
        /// </summary>
        void PrintUsage();
    }
}
