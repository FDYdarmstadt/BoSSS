using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using BoSSS.Application.ExternalBinding;

namespace BoSSS.Application.ExternalBinding.CodeGen {

    /// <summary>
    /// Main class of the C++ binding code generator
    /// </summary>
    static class CodeGenMain {

        static Type[] TypesToExport = new Type[] { typeof(GridServer) };

        static Type[] SupportedPrimitiveTypes = new Type[] { typeof(int), typeof(int*), typeof(int**), typeof(double), typeof(double*), typeof(double**), typeof(int).MakeByRefType()};
        static string[] SupportedPrimitiveTypesSS = new [] { "int", "int*", "int**", "double", "double*", "double**", "int*" };

        static bool IsPrimitiveType(Type t) {
            return SupportedPrimitiveTypes.Contains(t);
        }
        static bool IsExportedRefType(Type t) {
            return TypesToExport.Contains(t);
        }
        static string FormatType4C(Type t) {
            if (IsExportedRefType(t))
                return t.Name.ToString();
            if (IsPrimitiveType(t))
                return SupportedPrimitiveTypesSS[Array.IndexOf(SupportedPrimitiveTypes, t)];

            throw new NotImplementedException("Unknown Type: " + t);
        }
        //static bool IsPointer(Type t) {
        //    t.IsPointer
        //}
               
        static string[] CppNamespace = new string[] { "BoSSScppWrapper" };
        
        static string FormatParameters(ParameterInfo[] ps) {
            using(var stw = new StringWriter()) {

                for(int i = 0; i < ps.Length; i++) {
                    var pi = ps[i];

                    
                    stw.Write(FormatType4C(pi.ParameterType));

                    stw.Write(" ");
                    stw.Write(pi.Name);

                    if (i < ps.Length - 1)
                        stw.Write(", ");
                }


                return stw.ToString();
            }
        }

        static List<CodeGenCppFile> Cppfiles = new List<CodeGenCppFile>();
        static List<CodeGenHeaderFile> Hfiles = new List<CodeGenHeaderFile>();


        static void ExportClass(Type t) {
            if (!t.IsClass)
                throw new ArgumentOutOfRangeException();
            if (t.IsAbstract)
                throw new NotSupportedException("Abstract classes are not supported yet.");
            if (t.IsGenericType)
                throw new NotSupportedException("Generics are not supported yet.");

            // create files
            // ============
            Console.WriteLine("Doing Class " + t.FullName);
            var Cfile = new CodeGenCppFile();
            var Hfile = new CodeGenHeaderFile();
            Cfile.FileName = t.Name;
            Hfile.FileName = t.Name;
            Cppfiles.Add(Cfile);
            Hfiles.Add(Hfile);

            // add includes to C++ file
            // ========================
            Hfile.IncludeDirectives.Add("#pragma once");
            Hfile.IncludeDirectives.Add("#include \"" + t.Name + CodeGenHeaderFile.HeaderFileSuffix + "\"");
            Cfile.IncludeDirectives.Add("#include \"" + t.Name + CodeGenHeaderFile.HeaderFileSuffix + "\"");

            // create namespaces
            // =================
            BracedSection Cnmnsp = null;
            BracedSection Hnmnsp = null;
            {
                if (CppNamespace == null || CppNamespace.Length <= 0) {
                    Cnmnsp = new BracedSection() { NoBraces = true };
                    Hnmnsp = new BracedSection() { NoBraces = true };
                    Cfile.MainCode.Add(Cnmnsp);
                    Hfile.MainCode.Add(Hnmnsp);

                } else {
                    foreach (string nmn in CppNamespace) {
                        string line = "namespace " + nmn;

                        var _Cnmnsp = new BracedSection();
                        var _Hnmnsp = new BracedSection();
                        _Cnmnsp.OutsideCode.Add(line);
                        _Hnmnsp.OutsideCode.Add(line);

                        Debug.Assert((Cnmnsp == null) == (Hnmnsp == null));
                        if (Cnmnsp == null) {
                            Cfile.MainCode.Add(_Cnmnsp);
                            Hfile.MainCode.Add(_Hnmnsp);
                            Cnmnsp = _Cnmnsp;
                            Hnmnsp = _Hnmnsp;
                        } else {
                            Cnmnsp.Children.Add(_Cnmnsp);
                            Hnmnsp.Children.Add(_Hnmnsp);
                            Cnmnsp = _Cnmnsp;
                            Hnmnsp = _Hnmnsp;
                        }
                    }

                }
            }
            Cfile.ToString();
            Hfile.ToString();


            // class declaration
            // =================
            BracedSection ClassDecl = new BracedSection();
            ClassDecl.OutsideCode.Add("class " + t.Name);
            ClassDecl.ClosingSemicolon = true;
            Hnmnsp.Children.Add(ClassDecl);

            // Init code
            // =========
            /*
int GridProxy::_InitMonoBindings()
	{
        //printf("init mono bindings.\n");
        _ClassHandle = MonoBoSSSglobals::LookupClass("GridServer", "BoSSS.Application.ExternalBinding");
        //printf("got class.\n");
        _ctor = MonoBoSSSglobals::LookupMethod(_ClassHandle, "BoSSS.Application.ExternalBinding.GridServer:.ctor", true);
        //printf("got ctor.\n");
        return 0;
	}             */


            BracedSection InitCode = new BracedSection();
            Cnmnsp.Children.Add(InitCode);
            InitCode.OutsideCode.Add("void " + t.Name + "::_InitMonoBindings()");
            InitCode.Children.Add(string.Format(
                "_ClassHandle = MonoBoSSSglobals::LookupClass(\"{0}\", \"{1}\");",
                t.Name,
                t.Namespace));

            BracedSection PublicMethodDecl = new BracedSection();
            ClassDecl.Children.Add(PublicMethodDecl);
            PublicMethodDecl.NoBraces = true;
            PublicMethodDecl.AddOutside("public:");


            BracedSection PrivateBindingDecl = new BracedSection();
            ClassDecl.Children.Add(PrivateBindingDecl);
            PrivateBindingDecl.NoBraces = true;
            PrivateBindingDecl.AddOutside("private:");
            PrivateBindingDecl.AddInner("void _InitMonoBindings();");
            PrivateBindingDecl.AddInner("MonoClass* _ClassHandle;");
            PrivateBindingDecl.AddInner("uint32_t _MonoGCHandle;");

            // Destructor
            // ==========
            {
                PublicMethodDecl.OutsideCode.Add("~" + t.Name + "();");

                BracedSection dtorCode = new BracedSection();
                Cnmnsp.Children.Add(dtorCode);
                dtorCode.AddOutside(t.Name + "::~" + t.Name + "()");
                if (t.GetInterface(typeof(IDisposable).ToString()) != null) {
                    dtorCode.AddInner("Dispose();");
                }

                dtorCode.AddInner("mono_gchandle_free(_MonoGCHandle);");
            }




            // method wrappers
            // ===============

            var methods = t.GetMethods(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly);
            var ctors = t.GetConstructors(BindingFlags.Public | BindingFlags.Instance);
            int CtorCounter = 0;
            foreach (var m in ctors) {
                if (m.IsAbstract) {
                    Console.WriteLine("Skipping abstract method: " + m.Name);
                    continue;
                }

                var Params = m.GetParameters();
                string sParams = FormatParameters(Params);

                // init code
                // ---------
                PrivateBindingDecl.AddInner("MonoMethod* _ctor_{0};", CtorCounter);
                InitCode.AddInner("_ctor_{0} = MonoBoSSSglobals::LookupMethod(_ClassHandle, \"{1}\", true);", 
                    CtorCounter,
                    t.FullName + ":" + ".ctor"  //"BoSSS.Application.ExternalBinding.GridServer:.ctor"
                    );
                
                // declaration in class
                // --------------------
                PublicMethodDecl.AddInner(t.Name + "(" + sParams + ");");

                // wrapper code
                // ------------

                BracedSection ctorImpl = new BracedSection();
                Cnmnsp.Children.Add(ctorImpl);
                ctorImpl.AddOutside("{0}::{0}({1})", t.Name, sParams);

                ctorImpl.AddInner("_InitMonoBindings();");

                // instantiate
                ctorImpl.AddInner("MonoObject* ThisObj = mono_object_new(MonoBoSSSglobals::domain, _ClassHandle);");
                ctorImpl.AddInner("_MonoGCHandle = mono_gchandle_new(ThisObj, true);");

                //argument 
                ctorImpl.AddInner("void* args[{0}];", Params.Length);
                for (int i = 0; i < Params.Length; i++) {
                    if (Params[i].ParameterType.IsPointer) {
                        ctorImpl.AddInner("args[{0}] = {1};", i, Params[i].Name);
                    } else if(Params[i].ParameterType.IsValueType) {
                        ctorImpl.AddInner("args[{0}] = &{1};", i, Params[i].Name);
                    } else {
                        throw new NotImplementedException("Todo: args of type " + Params[i].ParameterType);
                    }
                }

                // call mono
                ctorImpl.AddInner("MonoObject* exception;");
                ctorImpl.AddInner("MonoObject* retval;");
                ctorImpl.AddInner("retval = mono_runtime_invoke(_ctor_{0}, mono_gchandle_get_target(_MonoGCHandle), args, &exception);", CtorCounter);
                ctorImpl.AddInner("if (exception != NULL) {");
                ctorImpl.AddInner("    printf( \"got exception from C#\\n\");");
                ctorImpl.AddInner("}");

                CtorCounter++;
            }
            foreach (var m in methods) {
                if (m.IsAbstract) {
                    Console.WriteLine("Skipping abstract method: " + m.Name);
                    continue;
                }

                var Params = m.GetParameters();
                string sParams = FormatParameters(Params);


            }



        }



        static void Main(string[] args) {
            foreach(var t in TypesToExport) {
                ExportClass(t);
            }

            foreach(var Cf in Cppfiles) {
                Cf.WriteFile();
            }
            foreach(var Hf in Hfiles) {
                Hf.WriteFile();
            }
        }
    }
}
