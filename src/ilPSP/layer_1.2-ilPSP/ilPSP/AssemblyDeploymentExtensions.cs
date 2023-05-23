
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;

namespace ilPSP {

    /// <summary>
    /// Utility functions to copy apps
    /// </summary>
    static public class AssemblyDeploymentExtensions {


        static void GetAllDependentAssembliesRecursive(Assembly a, HashSet<Assembly> assiList, string SearchPath) {
            if(assiList.Contains(a))
                return;
            assiList.Add(a);
            //if(a.FullName.Contains("codeanalysis", StringComparison.InvariantCultureIgnoreCase))
            //    Console.Write("");

            string fileName = Path.GetFileName(a.Location);
            var allMatch = assiList.Where(_a => Path.GetFileName(_a.Location).Equals(fileName)).ToArray();
            if(allMatch.Length > 1) {
                throw new ApplicationException("internal error in assembly collection.");
            }

            foreach(AssemblyName b in a.GetReferencedAssemblies()) {
                Assembly na;
                try {
                    na = Assembly.Load(b);
                } catch(FileNotFoundException) {
                    string[] AssiFiles = ArrayTools.Cat(Directory.GetFiles(SearchPath, b.Name + ".dll"), Directory.GetFiles(SearchPath, b.Name + ".exe"));
                    if(AssiFiles.Length != 1) {
                        //throw new FileNotFoundException("Unable to locate assembly '" + b.Name + "'.");
                        //Console.WriteLine("Skipping: " + b.Name);
                        continue;
                    }
                    na = Assembly.LoadFile(AssiFiles[0]);

                }

                GetAllDependentAssembliesRecursive(na, assiList, SearchPath);
            }
        }

        /// <summary>
        /// all non-system assemblies upon which an assembly depends on
        /// </summary>
        public static Assembly[] GetAllDependentAssemblies(this Assembly a) {
            HashSet<Assembly> assiList = new HashSet<Assembly>();
            string SearchPath = Path.GetDirectoryName(a.Location);
            GetAllDependentAssembliesRecursive(a, assiList, SearchPath);

            return assiList.ToArray();
        }

        /// <summary>
        /// 
        /// </summary>
        static string[] GetManagedFileList(IEnumerable<Assembly> AllDependentAssemblies, string MainAssemblyDir) {
            List<string> files = new List<string>();

            foreach(var a in AllDependentAssemblies) {
                // new rule for .NET5: if the file is located in the same directory as the entry assembly, it should be deployed;
                // (in Jupyter, sometimes assemblies from some cache are used, therefore we cannot use the assembly location as a criterion)
                string DelpoyAss = Path.Combine(MainAssemblyDir, Path.GetFileName(a.Location));

                if(File.Exists(DelpoyAss)) {
                    files.Add(DelpoyAss);

                    string a_config = Path.Combine(MainAssemblyDir, DelpoyAss + ".config");
                    string a_runtimeconfig_json = Path.Combine(MainAssemblyDir, Path.GetFileNameWithoutExtension(DelpoyAss) + ".runtimeconfig.json");

                    foreach(var a_acc in new[] { a_config, a_runtimeconfig_json }) {
                        if(File.Exists(a_acc)) {
                            files.Add(a_acc);
                        }
                    }
                } else {
                    //Console.WriteLine("SKIPPING: " + DelpoyAss + " --- " + MainAssemblyDir);
                }
            }

            return files.ToArray();
        }

        /// <summary>
        /// copies an assembly and all dependencies to a certain <paramref name="Destination"/>
        /// </summary>
        public static void DeployAt(this Assembly entryAssembly, DirectoryInfo Destination) {
            var assiList = GetAllDependentAssemblies(entryAssembly);
            string MainAssemblyDir = Path.GetDirectoryName(entryAssembly.Location);

            string[] files = GetManagedFileList(assiList, MainAssemblyDir);

            foreach(var f in files) {
                File.Copy(f, Path.Combine(Destination.FullName, Path.GetFileName(f)));
            }

            // copy "runtimes" directory from .NET core/.NET5
            string runtimes_Src = Path.Combine(MainAssemblyDir, "runtimes");
            string runtimes_Dst = Path.Combine(Destination.FullName, "runtimes");
            if(Directory.Exists(runtimes_Src)) {
                CopyFilesRecursively(runtimes_Src, runtimes_Dst);
            }
        }

        static void CopyFilesRecursively(string sourcePath, string targetPath) {
            //Now Create all of the directories
            foreach(string dirPath in Directory.GetDirectories(sourcePath, "*", SearchOption.AllDirectories)) {
                Directory.CreateDirectory(dirPath.Replace(sourcePath, targetPath));
            }

            //Copy all the files & Replaces any files with the same name
            foreach(string newPath in Directory.GetFiles(sourcePath, "*.*", SearchOption.AllDirectories)) {
                File.Copy(newPath, newPath.Replace(sourcePath, targetPath), true);
            }
        }

    }
}
