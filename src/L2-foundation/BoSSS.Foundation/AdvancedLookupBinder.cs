using ilPSP;
using Newtonsoft.Json.Serialization;
using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Text;

namespace BoSSS.Foundation.IO {
    /// <summary>
    /// Helps within Jupyter notebook, where <see cref="jsonFormatter"/>
    /// sometimes is not able to resolve all types and assemblies automatically.
    /// </summary>
    class AdvancedLookupBinder : DefaultSerializationBinder {

        internal AdvancedLookupBinder() {

            //this.defBinder = new DefaultSerializationBinder();

            var allAssis = System.AppDomain.CurrentDomain.GetAssemblies();

            HashSet<Assembly> assiList = new HashSet<Assembly>();
            foreach(var a in allAssis)
                GetAllAssemblies(a, assiList);

            foreach(var a in assiList) {
                try {
                    var tt = new Dictionary<string, Type>();
                    knownTypes.Add(a.GetName().Name, tt);

                    var Types_in_a = a.GetExportedTypes();
                    foreach(var t in Types_in_a) {
                        tt.Add(t.FullName, t);
                    }
                } catch  (Exception ) {

                }
            }
        }


        /// <summary>
        /// Recursive collection of all dependencies of some assembly.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="assiList">
        /// Output, list where all dependent assemblies are collected.
        /// </param>
        private static void GetAllAssemblies(Assembly a, HashSet<Assembly> assiList) {
            if(a.IsDynamic)
                return;
            if(assiList.Contains(a))
                return;
            assiList.Add(a);

            //string fileName = Path.GetFileName(a.Location);
            //var allMatch = assiList.Where(_a => Path.GetFileName(_a.Location).Equals(fileName)).ToArray();
            //if(allMatch.Length > 1) {
            //    throw new ApplicationException("internal error in assembly collection.");
            //}


            foreach(AssemblyName b in a.GetReferencedAssemblies()) {
                Assembly na;
                try {
                    na = Assembly.Load(b);
                } catch(FileNotFoundException) {
                    //string[] AssiFiles = ArrayTools.Cat(Directory.GetFiles(SearchPath, b.Name + ".dll"), Directory.GetFiles(SearchPath, b.Name + ".exe"));
                    //if(AssiFiles.Length != 1) {
                    //    //throw new FileNotFoundException("Unable to locate assembly '" + b.Name + "'.");
                    //    Console.WriteLine("Skipping: " + b.Name);
                    //    continue;
                    //}
                    //na = Assembly.LoadFile(AssiFiles[0]);
                    continue;
                }

                GetAllAssemblies(na, assiList);
            }
        }


        Dictionary<string, Dictionary<string, Type>> knownTypes = new Dictionary<string, Dictionary<string, Type>>();


        public override Type BindToType(string assemblyName, string typeName) {
            if(assemblyName != null && knownTypes.TryGetValue(assemblyName, out var typesInAssembly)) {
                if(typeName != null && typesInAssembly.TryGetValue(typeName, out Type t)) {
                    return t;
                }
            }


            return base.BindToType(assemblyName, typeName);
        }

        /*
        internal DefaultSerializationBinder defBinder;

        public override Type BindToType(string assemblyName, string typeName) {
            var t = defBinder.BindToType(assemblyName, typeName);
            return t;
            /*
            try {
            } catch(Exception) { }

            if(typeName.IsEmptyOrWhite())
                return null;

            //Console.WriteLine("Type lookup: " + assemblyName + "+" + typeName);
            var dd = knownTypes[assemblyName];
            var tt = dd[typeName];

            return tt;
        }
        */
    }
}
