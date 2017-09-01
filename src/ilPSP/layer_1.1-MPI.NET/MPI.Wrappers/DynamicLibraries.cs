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

namespace MPI.Wrappers.Utils {

    /// <summary>
    /// this class provides platform-independent handling of unmanaged dynamic libraries.
    /// </summary>
    public static class DynamicLibraries {

        /// <summary>
        /// a handle to unmanaged dynamic libraries
        /// </summary>
        public struct DynLibHandle {
            /// <summary> memory address </summary>
            public IntPtr val;
        }

        /// <summary>
        /// Unix command to open a shared object (lib*.so), aka. DLL, called <paramref name="filename"/>;
        /// </summary>
        //[DllImport("libdl.so.2", EntryPoint="dlopen", CharSet = CharSet.Ansi)]
        [DllImport("dl", CharSet = CharSet.Ansi)]
        static extern IntPtr dlopen(string filename, int flag);

        //[DllImport("fakedl", EntryPoint="my_dlerror")]
        [DllImport("dl")]
        unsafe static extern byte* dlerror();

        /// <summary>
        /// wrapper around <see cref="dlerror"/>
        /// </summary>
        /// <returns>a .NET string, instead of the the unmanaged pointer from <see cref="dlerror"/></returns>
        static string _dlerror() {
            unsafe {
                //Console.Write("fetching error string ...");
                IntPtr errStr = (IntPtr)dlerror();
                //Console.WriteLine("got error string pointer: " + errStr.ToString());
                //return "nonsense";
                if (errStr == IntPtr.Zero)
                    return null;
                else
                    return Marshal.PtrToStringAnsi((IntPtr)errStr);
            }
        }

        /// <summary>
        /// Unix command to load a symbol <paramref name="symbol"/> from a shared object,
        /// (lib*.so), aka. DLL;
        /// </summary>
        [DllImport("dl", CharSet = CharSet.Ansi)]
        unsafe static extern IntPtr dlsym(DynLibHandle handle, string symbol);

        /// <summary>
        /// Windows command to load a symbol <paramref name="lpProcName"/> from a dynamic library;
        /// </summary>
        [DllImportAttribute("kernel32.dll")]
        static extern IntPtr GetProcAddress(DynLibHandle hModule, [MarshalAsAttribute(UnmanagedType.LPStr)] string lpProcName);

        /// <summary>
        /// Unix command to close a shared object.
        /// </summary>
        [DllImport("dl")]
        static extern int dlclose(DynLibHandle handle);

        /// <summary>
        /// Windows command to open a dynamic library (*.dll), aka. DLL, called <paramref name="lpFileName"/>;
        /// </summary>
        [DllImportAttribute("kernel32.dll")]
        public static extern IntPtr LoadLibrary([InAttribute()] [MarshalAsAttribute(UnmanagedType.LPStr)] string lpFileName);

        /// <summary>
        /// Controls whether the system will handle the specified types of
        /// serious errors, or whether this process will handle them.
        /// </summary>
        /// <param name="newMode">
        /// The new error mode to be set
        /// </param>
        /// <returns>
        /// The error mode before applying the changes
        /// </returns>
        [DllImport("kernel32.dll")]
        public static extern ErrorModes SetErrorMode(ErrorModes newMode);

        /// <summary>
        /// Valid error modes for <see cref="SetErrorMode(ErrorModes)"/>
        /// </summary>
        [Flags]
        public enum ErrorModes : uint {

            /// <summary>
            /// Use the system default, which is to display all error dialog boxes.
            /// </summary>
            SYSTEM_DEFAULT = 0x0,

            /// <summary>
            /// The system does not display the critical-error-handler message
            /// box. Instead, the system sends the error to the calling process.
            /// </summary>
            SEM_FAILCRITICALERRORS = 0x0001,

            /// <summary>
            /// The system automatically fixes memory alignment faults and
            /// makes them invisible to the application.
            /// </summary>
            SEM_NOALIGNMENTFAULTEXCEPT = 0x0004,

            /// <summary>
            /// The system does not display the Windows Error Reporting dialog.
            /// </summary>
            SEM_NOGPFAULTERRORBOX = 0x0002,

            /// <summary>
            /// The OpenFile function does not display a message box when it
            /// fails to find a file. Instead, the error is returned to the caller.
            /// </summary>
            SEM_NOOPENFILEERRORBOX = 0x8000
        }

        /// <summary>
        /// last windows error, formatted as string
        /// </summary>
        /// <returns></returns>
        static string GetLastWin32Error() {
            //int errcode = Marshal.GetLastWin32Error();
            int hresult = Marshal.GetHRForLastWin32Error();
            Exception e = Marshal.GetExceptionForHR(
                hresult,
                new IntPtr(-1));

            return e.GetType().FullName + ": " + e.Message;
        }


        /// <summary>
        /// Unix command to close a shared object.
        /// </summary>
        [DllImportAttribute("kernel32.dll")]
        static extern bool FreeLibrary([InAttribute()] DynLibHandle hModule);


        /// <summary>
        /// loads a shared library (on Windows: dynamic link library, *.dll; on Unix: shared object, lib*.so)
        /// </summary>
        /// <param name="LibName"></param>
        /// <returns>library handle, for the use in functions <see cref="LoadSymbol"/> and <see cref="UnloadDynLib"/></returns>
        /// <remarks>
        /// On Windows, this function redirects to the 'LoadLibrary'-function in kernel32.dll;<br/>
        /// On Unix, this function redirects to the 'dlopen'-function in libdl.so;
        /// </remarks>
        /// <param name="errstr">
        /// on success, null;
        /// if call failed (return value is null), an error information provided by the operating system
        /// </param>
        /// <param name="loadGlobal">
        /// only effective on Linux/UNIX systems; if true, a library is loaded with the 'RTDL_GLOBAL' flag:
        /// thereby, it can be used automatically by the operating system to resolve symbols in other libraries.
        /// </param>
        public static DynLibHandle LoadDynLib(string LibName, out string errstr, bool loadGlobal) {
            PlatformID plattid = System.Environment.OSVersion.Platform;
            errstr = null;
            DynLibHandle ret;
            switch (plattid) {
                case PlatformID.Win32NT:
                    // Try to load but suppress ugly dialog box error
                    ErrorModes originalMode = SetErrorMode(ErrorModes.SEM_FAILCRITICALERRORS);
                    ret.val = LoadLibrary(LibName);
                    SetErrorMode(originalMode);

                    if (ret.val == IntPtr.Zero)
                        errstr = GetLastWin32Error();
                    break;
                case PlatformID.Unix:
                    ret.val = dlopen(LibName, loadGlobal ? (2 | 256) : 2); // 2 == RTLD_NOW, 256 == RTDL_GLOBAL
                    if (ret.val == IntPtr.Zero) {
                        errstr = _dlerror();
                    }
                    break;
                default:
                    throw new NotImplementedException("Dynamic Library Loading for " + plattid + " is not supported.");
            }
            return ret;
        }

        /// <summary>
        /// gets the address of a symbol (function or global variable) from an unmanaged library
        /// </summary>
        /// <param name="SymbName"></param>
        /// <summary>
        /// releases a dynamic library
        /// </summary>
        /// <param name="libHandle">
        /// handle, which was acquired by <see cref="LoadDynLib"/>
        /// </param>
        /// <returns>
        /// a function/symbol pointer;<br/>
        /// using the function <see cref="System.Runtime.InteropServices.Marshal.GetDelegateForFunctionPointer"/>
        /// it can be converted to a .NET delegate.
        /// </returns>
        /// <remarks>
        /// On Windows, this function redirects to the 'GetProcAddress'-function in kernel32.dll;<br/>
        /// On Unix, this function redirects to the 'dlsym'-function in libdl.so;
        /// </remarks>
        /// <param name="errstr">
        /// on success, null;
        /// if call failed (return value is null), an error information provided by the operating system
        /// </param>
        public static IntPtr LoadSymbol(DynLibHandle libHandle, string SymbName, out string errstr) {
            PlatformID plattid = System.Environment.OSVersion.Platform;
            IntPtr ret;
            errstr = null;

            Debug.Assert(Marshal.SizeOf(typeof(DynLibHandle)) == IntPtr.Size, "runtime created struct of wrong size.");

            switch (plattid) {
                case PlatformID.Win32NT:
                    ret = GetProcAddress(libHandle, SymbName);
                    if (ret == IntPtr.Zero)
                        errstr = GetLastWin32Error();
                    return ret;
                case PlatformID.Unix:
                    //Console.WriteLine("Trying to load >" + SymbName + "< form lib >" + libHandle.val + "< ...");
                    ret = dlsym(libHandle, SymbName);
                    if (ret == IntPtr.Zero) {
                        errstr = _dlerror();
                    }
                    return ret;
                default:
                    throw new NotImplementedException("Symbol Loading for " + plattid + " is not supported.");
            }
        }

        /// <summary>
        /// releases a dynamic library
        /// </summary>
        /// <param name="libHandle">
        /// handle, which was acquired by <see cref="LoadDynLib"/>
        /// </param>
        /// <remarks>
        /// On Windows, this function calls the 'FreeLibrary'-function in kernel32.dll;<br/>
        /// On Unix, this function calls the 'dlclose'-function in libdl.so;
        /// </remarks>
        public static void UnloadDynLib(DynLibHandle libHandle) {
            PlatformID plattid = System.Environment.OSVersion.Platform;
            switch (plattid) {
                case PlatformID.Win32NT:
                    FreeLibrary(libHandle);
                    break;
                case PlatformID.Unix:
                    dlclose(libHandle);
                    break;
                default:
                    throw new NotImplementedException("Dynamic Library Release for " + plattid + " is not supported.");
            }
        }
    }
}

