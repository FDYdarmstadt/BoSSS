using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using BoSSS.Application.ExternalBinding;

namespace BoSSS.Application.ExternalBinding.CodeGen {

    /// <summary>
    /// Main class of the C++ binding code generator
    /// </summary>
    static class CodeGenMain {

        static Type[] TypesToExport = new Type[] { typeof(GridServer) };

        static string[] CppNamespace = new string[] { "BoSSScppBinding" };


        static void ExportClass(Type t) {
            if (!t.IsClass)
                throw new ArgumentOutOfRangeException();
            
            // create files
            // ============
            var Cfile = new CodeGenCppFile();
            var Hfile = new CodeGenHeaderFile();
            Cfile.FileName = t.Name;
            Hfile.FileName = t.Name;

            // add includes to C++ file
            // ========================
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
                            Cfile.MainCode.Add(Cnmnsp);
                            Hfile.MainCode.Add(Hnmnsp);
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
            Cnmnsp.Children.Add(ClassDecl);


        }



        static void Main(string[] args) {

        }
    }
}
