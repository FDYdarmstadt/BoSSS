using System;
using System.Reflection;
using ilPSP;
using System.Runtime.InteropServices;

namespace BoSSS.Platform {


    /// <summary>
    /// loads 
    /// </summary>
    public abstract class DynLibLoader : IDisposable {

        public delegate string GetNameMangling(string s);


        DynamicLibraries.DynLibHandle m_LibHandle;


        static bool IsDelegate(Type t) {
            if (t.Equals(typeof(Delegate)))
                return true;
            if (t.Equals(typeof(MulticastDelegate)))
                return true;
            if (t.Equals(typeof(object)))
                return false;

            return IsDelegate(t.BaseType);
        }
                
        /// <summary>
        /// 
        /// </summary>
        /// <param name="LibNames">
        /// a list of library names
        /// </param>
        /// <param name="OsFilter">
        /// Operating system filter
        /// </param>
        /// <param name="NameMangling">
        /// 
        /// name mangling, suggestions are: <see cref="LeadingUnderscore_SmallLetters"/>, <see cref="SmallLetters_TrailingUnderscore"/>
        /// or <see cref="CAPITAL_LETTERS"/>
        /// </param>
        protected DynLibLoader(string[] LibNames, PlatformID[] OsFilter, GetNameMangling[] NameMangling) {
            if (LibNames.Length != OsFilter.Length || OsFilter.Length != NameMangling.Length)
                throw new ApplicationException("all arrays must have the same length.");

            PlatformID CurrentSys = System.Environment.OSVersion.Platform;

            for (int i = 0; i < LibNames.Length; i++) {
                if (CurrentSys != OsFilter[i])
                    continue;

                m_LibHandle = DynamicLibraries.LoadDynLib(LibNames[i]);
                if(m_LibHandle.val == IntPtr.Zero)
                    continue;

                Type myType = this.GetType();
                FieldInfo[] fields = myType.GetFields(BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public);

                // loop over all delegates ....
                foreach (FieldInfo fld in fields) {
                    if (IsDelegate(fld.FieldType)) {
                        // get function name in DLL
                        string UnmanagedName = NameMangling[i](fld.Name);
                        
                        // get function pointer
                        IntPtr FuncPtr = DynamicLibraries.LoadSymbol(m_LibHandle, UnmanagedName);
                        if (FuncPtr == IntPtr.Zero) {
                            //throw new ApplicationException("Library '" + LibNames[i] + "' not working - missing function '" + UnmanagedName + "';");
                            DynamicLibraries.UnloadDynLib(m_LibHandle);
                        }

                        // create delegate
                        fld.SetValue(this,Marshal.GetDelegateForFunctionPointer(FuncPtr,fld.FieldType));
                    }
                }

                // successfully loaded all library functions
                return;
            }

            // error
            throw new ApplicationException("unable to find/load dynamic library - none supported on actual system.");
        }


        public static string CAPITAL_LETTERS(string Name) {
            return Name.ToUpperInvariant();
        }

        public static string SmallLetters_TrailingUnderscore(string Name) {
            return Name.ToLowerInvariant() + "_";
        }

        public static string LeadingUnderscore_SmallLetters(string Name) {
            return "_" + Name.ToLowerInvariant();
        }


        /// <summary>
        /// dtor - calls dispose
        /// </summary>
        ~DynLibLoader() {
            this.Dispose();
        }
        
        #region IDisposable Members

        /// <summary>
        /// unloads the library
        /// </summary>
        public void Dispose() {
            if (m_LibHandle.val != IntPtr.Zero) {
                DynamicLibraries.UnloadDynLib(m_LibHandle);
                m_LibHandle.val = IntPtr.Zero;
            }
        }

        #endregion
    }
}
