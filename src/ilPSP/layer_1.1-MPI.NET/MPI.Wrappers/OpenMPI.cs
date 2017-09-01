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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace MPI.Wrappers {


    //    /// <summary>
    ///// this class provides platform-independet handling of unmanaged dynamic libraries.
    ///// </summary>
    //public static class DynamicLibraries {

    //    /// <summary>
    //    /// a handle to unmanaged dynamic libraries
    //    /// </summary>
    //    public struct DynLibHandle {
    //        /// <summary> memory address </summary>
    //        public IntPtr val;
    //    }

    //    /// <summary>
    //    /// Unix command to open a shared object (lib*.so), aka. DLL, called <paramref name="filename"/>;
    //    /// </summary>
    //    [DllImport("dl", CharSet = CharSet.Ansi)]
    //    static extern DynLibHandle dlopen(string filename, int flag);

    //    [DllImport("dl")]
    //    unsafe static extern byte* dlerror();
        
    //    static string _dlerror() {
    //        unsafe {
    //            IntPtr errStr = (IntPtr)dlerror();
    //            if (errStr == IntPtr.Zero)
    //                return null;
    //            else
    //                return Marshal.PtrToStringAnsi((IntPtr)errStr);
    //        }
    //    }

    //    /// <summary>
    //    /// Unix command to load a symbol <paramref name="symbol"/> from a shared object,
    //    /// (lib*.so), aka. DLL;
    //    /// </summary>
    //    [DllImport("dl", CharSet = CharSet.Ansi)]
    //    static extern IntPtr dlsym(DynLibHandle handle, string symbol);

    //    /// <summary>
    //    /// Windows command to load a symbol <paramref name="lpProcName"/> from a dynamic library;
    //    /// </summary>
    //    [DllImportAttribute("kernel32.dll")]
    //    static extern IntPtr GetProcAddress(DynLibHandle hModule,
    //                                              [MarshalAsAttribute(UnmanagedType.LPStr)] string lpProcName);
  
    //    /// <summary>
    //    /// Unix command to close a shared object.
    //    /// </summary>
    //    [DllImport("dl")]
    //    static extern int dlclose(DynLibHandle handle);
        
    //    /// <summary>
    //    /// Windows command to open a dynamic library (*.dll), aka. DLL, called <paramref name="lpFileName"/>;
    //    /// </summary>
    //    [DllImportAttribute("kernel32.dll")]
    //    static extern DynLibHandle LoadLibrary([InAttribute()] [MarshalAsAttribute(UnmanagedType.LPStr)] string lpFileName);

    //            /// <summary>
    //    /// Unix command to close a shared object.
    //    /// </summary>
    //    [DllImportAttribute("kernel32.dll")]
    //    static extern bool FreeLibrary([InAttribute()] DynLibHandle hModule);


    //    /// <summary>
    //    /// loads a shared library (on Windows: dynamic link library, *.dll; on Unix: shared object, lib*.so)
    //    /// </summary>
    //    /// <param name="LibName"></param>
    //    /// <returns>library handle, for the use in functions <see cref="LoadSymbol"/> and <see cref="UnloadDynLib"/></returns>
    //    /// <remarks>
    //    /// On Windows, this function redicts to the 'LoadLibrary'-function in kernel32.dll;<br/>
    //    /// On Unix, this function redicts to the 'dlopen'-function in libdl.so;
    //    /// </remarks>
    //    public static DynLibHandle LoadDynLib(string LibName) {
    //        PlatformID plattid = System.Environment.OSVersion.Platform;
    //        switch (plattid) {
    //            case PlatformID.Win32NT:
    //                return LoadLibrary(LibName);
    //            case PlatformID.Unix:
    //                return dlopen(LibName, 1);
    //            default:
    //                throw new NotImplementedException("Dynamic Library Loading for " + plattid + " is not supported.");
    //        }
    //    }

    //    /// <summary>
    //    /// gets the address of a symbol (function or global variable) from an unmanaged library
    //    /// </summary>
    //    /// <param name="SymbName"></param>
    //    /// <summary>
    //    /// releases a dynamic library
    //    /// </summary>
    //    /// <param name="libHandle">
    //    /// handle, which was acquired by <see cref="LoadDynLib"/>
    //    /// </param>
    //    /// <returns>
    //    /// a function/symbol pointer;<br/>
    //    /// using the function <see cref="System.Runtime.InteropServices.Marshal.GetDelegateForFunctionPointer"/>
    //    /// it can be converted to a .NET delegate.
    //    /// </returns>
    //    /// <remarks>
    //    /// On Windows, this function redicts to the 'GetProcAddress'-function in kernel32.dll;<br/>
    //    /// On Unix, this function redicts to the 'dlsym'-function in libdl.so;
    //    /// </remarks>
    //    public static IntPtr LoadSymbol(DynLibHandle libHandle, string SymbName) {
    //        PlatformID plattid = System.Environment.OSVersion.Platform;
    //        switch (plattid) {
    //            case PlatformID.Win32NT:
    //                return GetProcAddress(libHandle, SymbName);
    //            case PlatformID.Unix:
    //                return dlsym(libHandle,SymbName);
    //            default:
    //                throw new NotImplementedException("Symbol Loading for " + plattid + " is not supported.");
    //        }
    //    }

    //    /// <summary>
    //    /// releases a dynamic library
    //    /// </summary>
    //    /// <param name="libHandle">
    //    /// handle, which was acquired by <see cref="LoadDynLib"/>
    //    /// </param>
    //    /// <remarks>
    //    /// On Windows, this function redicts to the 'FreeLibrary'-function in kernel32.dll;<br/>
    //    /// On Unix, this function redicts to the 'dlclose'-function in libdl.so;
    //    /// </remarks>
    //    public static void UnloadDynLib( DynLibHandle libHandle) {
    //        PlatformID plattid = System.Environment.OSVersion.Platform;
    //        switch (plattid) {
    //            case PlatformID.Win32NT:
    //                FreeLibrary(libHandle); break;
    //            case PlatformID.Unix:
    //                dlclose(libHandle); break;
    //            default:
    //                throw new NotImplementedException("Dynamic Library Release for " + plattid + " is not supported.");
    //        }
    //    }
    //}



    //class OpenMPI {
    //    public OpenMPI() {
            

    //    }







    //}
}
