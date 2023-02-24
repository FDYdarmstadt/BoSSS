using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using BoSSS.Application.ExternalBinding;
using ilPSP;
using ilPSP.Connectors;
using ilPSP.Utils;


namespace BoSSS.Application.ExternalBinding.CodeGen {

    /// <summary>
    /// Main class of the C++ binding code generator
    /// </summary>
    static class CodeGenMain {

        static Type[] s_TypesToExport = null;
       
        /// <summary>
        /// all which implements <see cref="BoSSS.Application.ExternalBinding.IForeignLanguageProxy"/>;
        /// </summary>
        static Type[] TypesToExport { 
            get {
                if (s_TypesToExport == null) {
                    var assis = BoSSS.Solution.Application.GetAllAssemblies(typeof(CodeGenMain));
                    List<Type> classes = new List<Type>();
                    foreach (var a in assis) {
                        if(a.FullName.StartsWith("BoSSS") || a.FullName.StartsWith("ilPSP") || a.FullName.StartsWith("MPI.Wrappers")) {
                            foreach(var t in a.GetTypes()) {
                                if(t.IsClass && t.GetInterfaces().Contains(typeof(IForeignLanguageProxy))) {
                                    classes.Add(t);
                                }
                            }
                        }
                    }

                    // additional hooks 
                    classes.SetAdd(typeof(FixedOperators));
                    s_TypesToExport = classes.ToArray();
                }
                return s_TypesToExport.CloneAs();
            }
        }

        static Type[] SupportedPrimitiveTypes = new Type[] { typeof(int), typeof(int*), typeof(int**), typeof(double), typeof(double*), typeof(double**), typeof(int).MakeByRefType(), typeof(IntPtr), typeof(void)};
        static string[] SupportedPrimitiveTypesSS = new [] { "int", "int*", "int**", "double", "double*", "double**", "int*", "void*", "void" };

        static bool IsPrimitiveType(Type t) {
            return SupportedPrimitiveTypes.Contains(t);
        }
        static bool IsExportedRefType(Type t) {
            Debug.Assert(TypesToExport.Contains(t) == ImplementsMagicInterface(t));
            return TypesToExport.Contains(t);
        }
        static string FormatType4C(Type t) {
            if (TypesToExport.Contains(t))
                return t.FullName.Replace(".", "::") + "*";

            if (IsPrimitiveType(t))
                return SupportedPrimitiveTypesSS[Array.IndexOf(SupportedPrimitiveTypes, t)];

            throw new NotImplementedException("Unknown Type: " + t);
        }
        
        static bool IsFromMagicInterface(MethodInfo mi) {
            var t = mi.DeclaringType;
            if (!ImplementsMagicInterface(t))
                return false;

            var imap = t.GetInterfaceMap(typeof(IForeignLanguageProxy));

            foreach(var mInt in imap.InterfaceMethods) {
                if(mi.Name == mInt.Name) {
                    if(mi.ReturnType == mInt.ReturnType) {
                        var p1 = mi.GetParameters();
                        var p2 = mInt.GetParameters();

                        if (p1.Length != p2.Length)
                            continue;
                        bool found = true;
                        for(int ip = 0; ip < p1.Length; ip++) {
                            if(p1[ip].ParameterType != p2[ip].ParameterType) {
                                found = false;
                                break;
                            }
                        }

                        if (found)
                            return true;
                        else
                            continue;
                    }
                }
            }

            return false;
        }

        static bool ImplementsMagicInterface(Type t) {
            var Is = t.GetInterfaces();
            return Array.IndexOf(Is, typeof(IForeignLanguageProxy)) >= 0;
        }

         
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
            Hfile.IncludeDirectives.Add("#pragma once");



            // create namespaces
            // =================
            string[] CppNamespace = t.Namespace.Split('.');
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

            // class declaration
            // =================
            BracedSection ClassDecl = new BracedSection();
            ClassDecl.OutsideCode.Add("class " + t.Name);
            ClassDecl.ClosingSemicolon = true;
            Hnmnsp.Children.Add(ClassDecl);

            // Init code
            // =========

            BracedSection InitCode = new BracedSection();
            Cnmnsp.Children.Add(InitCode);
            InitCode.OutsideCode.Add("void " + t.Name + "::_InitMonoBindings()");
            InitCode.Children.Add(string.Format(
                "_ClassHandle = BoSSS::Globals::LookupClass(BoSSS::Globals::_image__{0}, \"{1}\", \"{2}\");",
                t.Assembly.GetName().Name.Replace(".", "_"),
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

                BracedSection cond = new BracedSection();
                dtorCode.Children.Add(cond);
                cond.AddOutside("if (_MonoGCHandle != 0)");

                cond.AddInner("// See also _ReleaseGChandle()  and _FromMonoObject(MonoObject*) methods:");
                cond.AddInner("// For a temporary object creation through _FromMonoObject(...), we only need parts of the destruction, therefore parts of the destructor are blocked.");

                if (t.GetInterface(typeof(IDisposable).ToString()) != null) {
                    cond.AddInner("Dispose();");
                }

                cond.AddInner("_SetForeignPointer(NULL);");
                cond.AddInner("mono_gchandle_free(_MonoGCHandle);");
            }

            // Release Method
            // ==============
            {
                PublicMethodDecl.AddInner("void _ReleaseGChandle();");

                BracedSection _ReleaseGChandleImpl = new BracedSection();
                Cnmnsp.Children.Add(_ReleaseGChandleImpl);
                _ReleaseGChandleImpl.AddOutside("void {0}::_ReleaseGChandle()", t.Name);

                _ReleaseGChandleImpl.AddInner("mono_gchandle_free(_MonoGCHandle);");
                _ReleaseGChandleImpl.AddInner("_MonoGCHandle = 0; // blocks destructor functionality");
            }

            // Constructor from MonoObject
            // ===========================
            {
                PublicMethodDecl.AddInner("{0}(MonoObject* mo);", t.Name);
                BracedSection ctorImpl = new BracedSection();
                Cnmnsp.Children.Add(ctorImpl);
                ctorImpl.AddOutside("{0}::{0}(MonoObject* mo)", t.Name);
                ctorImpl.AddInner("_InitMonoBindings();");
                ctorImpl.AddInner("_MonoGCHandle = mono_gchandle_new(mo, true);");
            }

            // _FromMonoObject(...) method
            // ===========================
            {
                PublicMethodDecl.AddInner("static {0}* _FromMonoObject(MonoObject* mo);", t.Name);

                BracedSection _FromMonoObjectImpl = new BracedSection();
                Cnmnsp.Children.Add(_FromMonoObjectImpl);
                _FromMonoObjectImpl.AddOutside("{0}* {0}::_FromMonoObject(MonoObject* mo)", t.Name);

                _FromMonoObjectImpl.AddInner("{0}* tmp = new {0}(mo);", t.Name);
                _FromMonoObjectImpl.AddInner("void* LoggedRef = tmp->_GetForeignPointer();");
                _FromMonoObjectImpl.AddInner("if (LoggedRef != NULL) {");
                _FromMonoObjectImpl.AddInner("    tmp->_ReleaseGChandle();");
                _FromMonoObjectImpl.AddInner("    delete tmp;");
                _FromMonoObjectImpl.AddInner("    return (({0}*)LoggedRef);", t.Name);
                _FromMonoObjectImpl.AddInner("} else {");
                _FromMonoObjectImpl.AddInner("    tmp->_SetForeignPointer(tmp);");
                _FromMonoObjectImpl.AddInner("    return tmp;");
                _FromMonoObjectImpl.AddInner("}");
            }


            // _GetMonoObject(...) method
            // ===========================
            {
                PublicMethodDecl.AddInner("MonoObject* _GetMonoObject();", t.Name);

                BracedSection _GetMonoObjectImpl = new BracedSection();
                Cnmnsp.Children.Add(_GetMonoObjectImpl);
                _GetMonoObjectImpl.AddOutside("MonoObject* {0}::_GetMonoObject()", t.Name);
                _GetMonoObjectImpl.AddInner("return mono_gchandle_get_target(_MonoGCHandle);");
            }


            // constructor wrappers
            // ====================

            var ctors = t.GetConstructors(BindingFlags.Public | BindingFlags.Instance);
            int CtorCounter = 0;
            foreach (var m in ctors) {
                if (m.IsAbstract) {
                    Console.WriteLine("Skipping abstract method: " + m.Name);
                    continue;
                }

                if (!(m.GetCustomAttributes().Any(att => att.GetType() == typeof(CodeGenExportAttribute))))
                    continue;


                var Params = m.GetParameters();
                string sParams = FormatParameters(Params);

                // init code
                // ---------
                PrivateBindingDecl.AddInner("MonoMethod* _ctor_{0};", CtorCounter);
                InitCode.AddInner("_ctor_{0} = BoSSS::Globals::LookupMethod(_ClassHandle, \"{1}\", true);", 
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
                ctorImpl.AddInner("MonoObject* ThisObj = mono_object_new(BoSSS::Globals::_domain, _ClassHandle);");
                ctorImpl.AddInner("_MonoGCHandle = mono_gchandle_new(ThisObj, true);");

                //argument 
                CreateMethodWrapper("ctor_" + CtorCounter, Params, ctorImpl, t);

                // link c++ to .NET object
                ctorImpl.AddInner("_SetForeignPointer(this);");

                // next 
                CtorCounter++;
            }

            // method wrappers
            // ===============
            var methods = t.GetMethods(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly);
            foreach (var m in methods) {

                if (!(m.GetCustomAttributes().Any(att => att.GetType() == typeof(CodeGenExportAttribute))
                    || IsFromMagicInterface(m)))
                    continue;

                if (m.IsAbstract) {
                    throw new ApplicationException("Abstract method is not supported; Skipping abstract method:; (" + m.Name + ")");
                }
                if (m.IsGenericMethod || m.IsGenericMethodDefinition) {
                    throw new ApplicationException("Generic method is not supported; (" + m.Name + ")");
                }



                var Params = m.GetParameters();
                string sParams = FormatParameters(Params);

                string sRetType = FormatType4C(m.ReturnType);

                // init code
                // ---------
                PrivateBindingDecl.AddInner("MonoMethod* _{0};", m.Name);
                InitCode.AddInner("_{0} = BoSSS::Globals::LookupMethod(_ClassHandle, \"{1}\", true);",
                    m.Name,
                    t.FullName + ":" + m.Name
                    );

                // declaration in class
                // --------------------
                PublicMethodDecl.AddInner(sRetType + " " + m.Name + "(" + sParams + ");");

                // wrapper code
                // ------------

                BracedSection methImpl = new BracedSection();
                Cnmnsp.Children.Add(methImpl);
                methImpl.AddOutside("{0} {1}::{2}({3})", sRetType, t.Name, m.Name, sParams);
                CreateMethodWrapper(m.Name, Params, methImpl, t);

                // return value handling
                if (m.ReturnType == typeof(void)) {
                    // nothing to do
                    methImpl.AddInner("return;");
                } else if (IsPrimitiveType(m.ReturnType)) {
                    // try with type-cast
                    methImpl.AddInner("void* retptr = mono_object_unbox(retval);");
                    methImpl.AddInner("return *(({0}*) retptr);", sRetType);
                } else if (IsExportedRefType(m.ReturnType)) {
                    // wrapper object required
                    methImpl.AddInner("return {0}::_FromMonoObject(retval);", m.ReturnType.FullName.Replace(".", "::"));
                } else {
                    throw new NotSupportedException("Unknown return type of wrapper function: " + m.ReturnType);
                }
            }



        }

        private static void CreateMethodWrapper(string _Name, ParameterInfo[] Params, BracedSection methImpl, Type declaringType) {
            //argument 
            methImpl.AddInner("void* args[{0}];", Math.Max(1, Params.Length));
            for (int i = 0; i < Params.Length; i++) {
                if (Params[i].ParameterType.IsPointer) {
                    methImpl.AddInner("args[{0}] = {1};", i, Params[i].Name);
                } else if (Params[i].ParameterType.IsValueType) {
                    methImpl.AddInner("args[{0}] = &{1};", i, Params[i].Name);
                } else if (TypesToExport.Contains(Params[i].ParameterType)) {
                    methImpl.AddInner("args[{0}] = {1}->_GetMonoObject();", i, Params[i].Name);
                } else {
                    throw new NotImplementedException("Todo: args of type " + Params[i].ParameterType);
                }
            }

            // call mono
            methImpl.AddInner("MonoObject* exception;");
            methImpl.AddInner("MonoObject* retval;");
            methImpl.AddInner("retval = mono_runtime_invoke(_{0}, mono_gchandle_get_target(_MonoGCHandle), args, &exception);", _Name);
            methImpl.AddInner("if (exception != NULL) {");
            methImpl.AddInner("    fprintf( stderr, \"got exception from C# (" + declaringType.Name + "." + _Name + ") \\n\");");
            methImpl.AddInner("}");
        }

        static void CreatePrototypeDeclarations() {

            // sort types according to namespace
            // =================================
            Dictionary<string, List<Type>> namespace2type = new Dictionary<string, List<Type>>();
            foreach(var t in TypesToExport) {
                var nmn = t.Namespace;
                List<Type> tts;
                if(!namespace2type.TryGetValue(nmn,out tts)) {
                    tts = new List<Type>();
                    namespace2type.Add(nmn, tts);
                }

                tts.Add(t);
            }


            var Hfile = new CodeGenHeaderFile();
            Hfile.FileName = "Prototypes";
            Hfiles.Add(Hfile);
            Hfile.IncludeDirectives.Add("#pragma once");


            // create namespaces
            // =================
            foreach (var ky in namespace2type) {
                string[] CppNamespace = ky.Key.Split('.');
                BracedSection Hnmnsp = null;
                {
                    if (CppNamespace == null || CppNamespace.Length <= 0) {
                        Hnmnsp = new BracedSection() { NoBraces = true };
                        Hfile.MainCode.Add(Hnmnsp);

                    } else {
                        foreach (string nmn in CppNamespace) {
                            string line = "namespace " + nmn;

                            var _Cnmnsp = new BracedSection();
                            var _Hnmnsp = new BracedSection();
                            _Cnmnsp.OutsideCode.Add(line);
                            _Hnmnsp.OutsideCode.Add(line);

                            if (Hnmnsp == null) {
                                Hfile.MainCode.Add(_Hnmnsp);
                                Hnmnsp = _Hnmnsp;
                            } else {
                                Hnmnsp.Children.Add(_Hnmnsp);
                                Hnmnsp = _Hnmnsp;
                            }
                        }
                    }
                }

                foreach(var t in ky.Value) {
                    Hnmnsp.AddInner("class {0};", t.Name);
                }
            }

        }

        static void GenerateMasterInclude() {
            var Hfile = new CodeGenHeaderFile();
            Hfile.FileName = "BoSSScpp";
            Hfile.IncludeDirectives.Add("#pragma once");

            Hfile.IncludeDirectives.Add("#include <mono/metadata/mono-config.h>");
            Hfile.IncludeDirectives.Add("#include <mono/jit/jit.h>");
            Hfile.IncludeDirectives.Add("#include <mono/metadata/assembly.h>");
            Hfile.IncludeDirectives.Add("#include <mono/metadata/debug-helpers.h>");

            foreach (var Hf in Hfiles) {
                Hfile.IncludeDirectives.Add("#include \"" + Hf.FileName + CodeGenHeaderFile.HeaderFileSuffix + "\"");
            }

            //Hfile.IncludeDirectives.Add("#include \"MonoBoSSSglobals.h\"");
            Hfile.IncludeDirectives.Add("#include \"Globals.h\"");

            Hfiles.Add(Hfile);
        }


        static void GenerateGlobals() {
            var Hfile = new CodeGenHeaderFile();
            Hfile.FileName = "Globals";
            Hfile.IncludeDirectives.Add("#pragma once");
            Hfiles.Add(Hfile);

            var Cppfile = new CodeGenCppFile();
            Cppfile.FileName = "Globals";
            Cppfile.IncludeDirectives.Add("#include \"Globals.h\"");
            Cppfiles.Add(Cppfile);
            
            var Hnmnsp = new BracedSection();
            Hfile.MainCode.Add(Hnmnsp);
            Hnmnsp.AddOutside("namespace BoSSS");

            var Cnmnsp = new BracedSection();
            Cppfile.MainCode.Add(Cnmnsp);
            Cnmnsp.AddOutside("namespace BoSSS");

            var Class = new BracedSection();
            Class.AddOutside("class Globals");
            Hnmnsp.Children.Add(Class);
            Class.ClosingSemicolon = true;

            var PublicSection = new BracedSection() { NoBraces = true };
            Class.Children.Add(PublicSection);
            PublicSection.AddOutside("public:");
            PublicSection.AddInner("static MonoDomain* _domain;");
            Cnmnsp.AddInner("MonoDomain* Globals::_domain = NULL;");
            
            var PrivateSection = new BracedSection() { NoBraces = true };
            Class.Children.Add(PrivateSection);
            PrivateSection.AddOutside("private:");
                        
            Assembly[] ass = TypesToExport.Select(tt => tt.Assembly).ToSet().ToArray();
            
            PrivateSection.AddInner("static int _Initialized;");
            Cnmnsp.AddInner("int Globals::_Initialized = 0;");

            PublicSection.AddInner("static void Init(char* ManagedAssemblyDirectory);");

            int stringlen = 1024;
            var InitMethod = new BracedSection();
            {
                Cnmnsp.Children.Add(InitMethod);
                InitMethod.AddOutside("void Globals::Init(char* ManagedAssemblyDirectory)");
                InitMethod.AddInner("if (_Initialized != 0)");
                InitMethod.AddInner("   return;");
                InitMethod.AddInner("_Initialized = 0x1234;");
                InitMethod.AddInner("");
                InitMethod.AddInner("char path[{0}];", stringlen + 4);
                InitMethod.AddInner("mono_config_parse(NULL);");
                InitMethod.AddInner("_domain = mono_jit_init_version(\"BoSSSdomain\", \"v4.0.30319\");");
                InitMethod.AddInner("if (_domain == NULL) {");
                InitMethod.AddInner("    fprintf(stderr, \"Unable to setup mono domain.\");");
                InitMethod.AddInner("    throw \"Unable to setup mono domain.\";");
                InitMethod.AddInner("}");
                InitMethod.AddInner("mono_domain_set_config(_domain, ManagedAssemblyDirectory, \"\");");
            }

            foreach (var a in ass) {
                string nmn = a.GetName().Name;
                string FileName = Path.GetFileName(a.Location); 
                if(FileName.EndsWith(".DLL")) {
                    FileName = FileName.Substring(0, FileName.Length - 4);
                    FileName = FileName + ".dll";
                }
                if(FileName.EndsWith(".EXE")) {
                    FileName = FileName.Substring(0, FileName.Length - 4);
                    FileName = FileName + ".exe";
                }

                
                nmn = nmn.Replace(".", "_");
                PublicSection.AddInner("static MonoAssembly* _assembly__{0};", nmn);
                PublicSection.AddInner("static MonoImage* _image__{0};", nmn);
                Cnmnsp.AddInner("MonoAssembly* Globals::_assembly__{0} = NULL;", nmn);
                Cnmnsp.AddInner("MonoImage* Globals::_image__{0} = NULL;", nmn);



                InitMethod.AddInner("if(strlen(ManagedAssemblyDirectory) + {0} >= {1}) {{", FileName.Length + 1, stringlen);
                InitMethod.AddInner("    fprintf(stderr, \"Path to assembly {0} exceeds {1} character limit. \\n \"); ", FileName, stringlen);
                InitMethod.AddInner("}");
                InitMethod.AddInner("strcpy(path, ManagedAssemblyDirectory);");
                InitMethod.AddInner("strcat(path, \"{0}\");", FileName);

                InitMethod.AddInner("_assembly__{0} = mono_domain_assembly_open(_domain, path);", nmn);
                InitMethod.AddInner("if (_assembly__{0} == NULL) {{", nmn);
                InitMethod.AddInner("    fprintf(stderr, \"Unable to open assembly: {1} \\n \"); ", nmn, FileName);
                InitMethod.AddInner("    throw \"Unable to open assembly.\";");
                InitMethod.AddInner("}");
                InitMethod.AddInner("_image__{0} = mono_assembly_get_image(_assembly__{0});",nmn);
                InitMethod.AddInner("if (_image__{0} == NULL) {{", nmn);
                InitMethod.AddInner("    fprintf(stderr, \"Unable to get assembly image {0}.\\n\");", FileName);
                InitMethod.AddInner("    throw \"Unable to get assembly image.\";");
                InitMethod.AddInner("}");
                
            }


            PublicSection.AddInner("static MonoClass*  LookupClass(MonoImage* image, const char *_classname, const char* _namespace);");
            PublicSection.AddInner("static MonoMethod* LookupMethod(MonoClass* pClass, const char *name, mono_bool include_namespace);");

            var LookupClassMethod = new BracedSection();
            {
                Cnmnsp.Children.Add(LookupClassMethod);
                LookupClassMethod.AddOutside("MonoClass*  Globals::LookupClass(MonoImage* image, const char *_classname, const char* _namespace)");
                LookupClassMethod.AddInner("if(image == NULL) {");
                LookupClassMethod.AddInner("    fprintf(stderr, \"Image not loaded - must call MonoBoSSSglobals::Init first.\\n\");");
                LookupClassMethod.AddInner("}");
                LookupClassMethod.AddInner("MonoClass* mcls = mono_class_from_name(image, _namespace, _classname);");
                LookupClassMethod.AddInner("if (mcls == NULL) {");
                LookupClassMethod.AddInner("    fprintf(stderr, \"unable to find class %s.%s\\n\", _namespace, _classname);");
                LookupClassMethod.AddInner("    throw \"unable to find method description for method. \";");
                LookupClassMethod.AddInner("}");
                LookupClassMethod.AddInner("return mcls;");
            }

            var LookupMethodMethod = new BracedSection();
            {
                Cnmnsp.Children.Add(LookupMethodMethod);
                LookupMethodMethod.AddOutside("MonoMethod* Globals::LookupMethod(MonoClass* pClass, const char *name, mono_bool include_namespace)");

                LookupMethodMethod.AddInner("MonoMethodDesc* desc = mono_method_desc_new(name, include_namespace);");
                LookupMethodMethod.AddInner("if (desc == NULL) {");
                LookupMethodMethod.AddInner("    fprintf(stderr, \"unable to find method description for method %s\\n\", name);");
                LookupMethodMethod.AddInner("    throw \"unable to find method description for method. \";");
                LookupMethodMethod.AddInner("}");
                LookupMethodMethod.AddInner("MonoMethod* methodHandle = mono_method_desc_search_in_class(desc, pClass);");
                LookupMethodMethod.AddInner("if (methodHandle == NULL) {");
                LookupMethodMethod.AddInner("    fprintf(stderr, \"unable to find method handle for method.\\n\");");
                LookupMethodMethod.AddInner("    throw \"unable to find method handle for method. \";");
                LookupMethodMethod.AddInner("}");
                LookupMethodMethod.AddInner("mono_method_desc_free(desc);");
                LookupMethodMethod.AddInner("return methodHandle;");
            }

        }


        public static void Main(string[] args) {
            // check input
            // ===========

            if(args.Length != 1 || !Directory.Exists(args[0])) {
                Console.WriteLine("BoSSS - Code generator for C++ bindings.");
                Console.WriteLine("usage: ");
                Console.WriteLine("      ExteranlBinding.CodeGen.exe {OuputDir}");

                if(args.Length < 1)
                    Console.Error.WriteLine("Unable to get output directory.");
                else if(!Directory.Exists(args[0]))
                    Console.Error.WriteLine("Directory >{0}< does not seem to exist.");

                Console.Error.WriteLine("exiting");
                return;
            }
            string outputDir = args[0];

            // delete existing header and source files
            // =======================================
            {
                foreach(string ext in new[] { "*.h", "*.cpp" }) {
                    var fs = Directory.GetFiles(outputDir, ext);
                    foreach(var f in fs) {
                        File.Delete(f);
                    }

                }

            }


            // create type decls
            // =================
            CreatePrototypeDeclarations();

            // create wrapper code 
            // ===================
            foreach (var t in TypesToExport) {
                ExportClass(t);
            }

            // add includes to C++ file
            // ========================
            foreach(var Cf in Cppfiles) {
                foreach (var Hf in Hfiles) {
                    Cf.IncludeDirectives.Add("#include \"" + Hf.FileName + CodeGenHeaderFile.HeaderFileSuffix + "\"");
                }
            }

            // Master include
            // ==============
            GenerateMasterInclude();

            // Mono interface wrappers
            // =======================

            GenerateGlobals();


            // write code
            // ==========

            foreach(var Cf in Cppfiles) {
                Cf.WriteFile(outputDir);
            }
            foreach(var Hf in Hfiles) {
                Hf.WriteFile(outputDir);
            }

            // Write additional files
            // ======================

            File.WriteAllText(Path.Combine(outputDir, "compile.sh"), Resources.compile_sh);
            //File.WriteAllText(Path.Combine(outputDir, "MonoBoSSSglobals.h"), Resources.MonoBoSSSglobals_h);
            //File.WriteAllText(Path.Combine(outputDir, "MonoBoSSSglobals.cpp"), Resources.MonoBoSSSglobals_cpp);
            File.WriteAllText(Path.Combine(outputDir, "ExtBindingTest.cpp"), Resources.ExtBindingTest_cpp);

        }
    }
}
