using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace bcl.nativelib_inst {
    
    class Program : IProgram {
        #region IProgram Members

        public void Execute() {
            CopyNative(StartDir);
        }

        DirectoryInfo StartDir = bcl.myEnv.BOSSS_ROOT;

        public void DecodeArgs(string[] args) {
            if (args.Length != 1)
                throw new UserInputException("wrong number of arguments");

            if (args.Length > 0) {
                DirectoryInfo di = new DirectoryInfo(bcl.BrowsePath(args[0]));
                if (!di.Exists)
                    throw new UserInputException("'" + args[0] + "' either does not exist or is no directory.");
                StartDir = di;
            }
        }

        public void PrintUsage() {
            Console.WriteLine("Usage: bcl nativelib-inst $directory");
            Console.WriteLine();
            Console.WriteLine("Deploys the recent version of native the libraries");
            Console.WriteLine("into the specified directory");
            Console.WriteLine("Arguments:");
            Console.WriteLine(" directory      directory to start the recursive installation");
            Console.WriteLine("                of native libraries; if not specified, $BOSSS_ROOT is ");
            Console.WriteLine("                used");
            Console.WriteLine();
        }

        #endregion

        /// <summary>
        /// Copies the native binaries to 
        /// <paramref name="targetDirectory"/>.
        /// </summary>
        private void CopyNative(DirectoryInfo targetDirectory) {

            bool MustContainNative = true;

            if (MustContainNative) {
                Console.WriteLine("INFO: installing native binaries in :'" + targetDirectory.FullName + "'");
                CopyDirectoryRec(bcl.myEnv.BOSSS_BIN_NATIVE, targetDirectory, "amd64");
            }
        }

        /// <summary>
        /// (tries to) do a recursive copy of a directory
        /// </summary>
        /// <param name="srcDir"></param>
        /// <param name="dstDir"></param>
        /// <param name="filter">search pattern/filter</param>
        public static void CopyDirectoryRec(DirectoryInfo srcDir, DirectoryInfo dstDir, string filter) {
            FileInfo[] srcFiles = srcDir.GetFiles();


            foreach (FileInfo srcFile in srcFiles) {
                TryCopy(srcFile.FullName, Path.Combine(dstDir.FullName, srcFile.Name));
            }

            DirectoryInfo[] subDirs;
            if (filter == null)
                subDirs = srcDir.GetDirectories();
            else
                subDirs = srcDir.GetDirectories(filter);
            foreach (DirectoryInfo srcSubDir in subDirs) {
                DirectoryInfo dstSubDir = new DirectoryInfo(Path.Combine(dstDir.FullName, srcSubDir.Name));
                if (!dstSubDir.Exists)
                    dstSubDir.Create();
                CopyDirectoryRec(srcSubDir, dstSubDir, null);
            }
        }

        /// <summary>
        /// Utility function which tries to copy a file from
        /// <paramref name="sourceFileName"/> to
        /// <paramref name="destFileName"/> overwriting existing files if
        /// required. Issues a warning (but proceeds as normal) if the copy
        /// process fails.
        /// </summary>
        /// <param name="sourceFileName">
        /// The path to the file to be copied
        /// </param>
        /// <param name="destFileName">The path to the destination</param>
        private static void TryCopy(string sourceFileName, string destFileName) {
            try {
                File.Copy(sourceFileName, destFileName, true);
            } catch (Exception e) {
                Console.WriteLine("WARNING: Unable to copy to: '"
                    + destFileName + "': " + e.GetType().Name + " says:'" + e.Message + "'");
            }
        }
    }
     
}
