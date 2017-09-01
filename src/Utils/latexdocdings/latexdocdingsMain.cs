using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Xml;

namespace latexdocdings {
    class latexdocdingsMain {

        /// <summary>
        /// Recursively walks through the dependencies of the .NET assembly
        /// <paramref name="rootAssemblyFile"/> 
        /// and collects them in
        /// <paramref name="assemblyPaths"/>.
        /// </summary>
        /// <param name="rootAssemblyFile">The root of the recursion.</param>
        /// <param name="assemblies">
        /// output;
        /// </param>
        static private void CollectAllDependenciesRecursive(string rootAssemblyFilePath, List<AssemblyAndOtherStuff> assemblies) {
            if(!File.Exists(rootAssemblyFilePath))
                throw new FileNotFoundException();

            string searchPath = Path.GetDirectoryName(rootAssemblyFilePath);
            string rootAssemblyFile = Path.GetFileNameWithoutExtension(rootAssemblyFilePath);



            if(assemblies.SingleOrDefault(a => a.assi.GetName().Name  == rootAssemblyFile) != null) {
                // Already added, end of recursion
                return;
            } else {
                string[] allFiles = Directory.GetFiles(searchPath);

                AssemblyAndOtherStuff a = new AssemblyAndOtherStuff();
                assemblies.Add(a);

                // load assembly
                // -------------
                a.FilePath = rootAssemblyFile;
                a.assi = Assembly.LoadFile(rootAssemblyFilePath);
                

                // load xmldoc, if present
                // -----------------------

                string xmldocFile = allFiles.SingleOrDefault(delegate(string ss) {
                    string ssName = Path.GetFileNameWithoutExtension(ss);
                    string ssLow = ss.ToLower();
                    bool fileMatch = ssName.Equals(rootAssemblyFile);
                    bool extMatch =  ssLow.EndsWith(".xml");
                    return (fileMatch && extMatch);
                });
                
                if(xmldocFile != null) {
                    XmlDocument xmldoc = new XmlDocument();
                    xmldoc.Load(xmldocFile);
                    a.xmldoc = xmldoc;
                    //Console.WriteLine("Found XML for: " + xmldocFile);
                } else {
                    a.xmldoc = null;
                    //Console.WriteLine("Missing XML for: " + rootAssemblyFile);
                }
                
                // Collect all directly referenced assemblies (recursion)
                // ------------------------------------------------------
                
                foreach(AssemblyName dependencyName in a.assi.GetReferencedAssemblies()) {
                    //string refAssemblyfilePath = Path.Combine(searchPath, dependencyName.FullName);
                    
                    string refAssemblyFile = allFiles.SingleOrDefault(delegate(string ss) {
                        string ssName = Path.GetFileNameWithoutExtension(ss);
                        string ssLow = ss.ToLower();
                        bool fileMatch = ssName.Equals(dependencyName.Name);
                        bool extMatch =  (ssLow.EndsWith(".dll") || ssLow.EndsWith(".exe"));
                        return (fileMatch && extMatch);

                    });


                    if(refAssemblyFile != null) {
                        //Console.WriteLine(refAssemblyFile);
                        CollectAllDependenciesRecursive(refAssemblyFile, assemblies);
                    } else {

                    }
                }
            }
        }


        /// <summary>
        /// Links assembly and XMLDOC -- files.
        /// </summary>
        class AssemblyAndOtherStuff {
            /// <summary>
            /// file path to the assembly
            /// </summary>
            public string FilePath;


            /// <summary>
            /// Optional XMLDOC, but can be null.
            /// </summary>
            public XmlDocument xmldoc;

            /// <summary>
            /// the assembly itself
            /// </summary>
            public Assembly assi;
        }

        /// <summary>
        /// All loaded assemblies.
        /// </summary>
        static List<AssemblyAndOtherStuff> assemblies = new List<AssemblyAndOtherStuff>();
        
        /// <summary>
        /// "meta-type" of some .NET member, as specified by the XMLDOC format.
        /// </summary>
        enum KindOfMemeber {

            /// <summary>
            /// namespace
            /// </summary>
            N, 
            
            /// <summary>
            /// type: class, interface, struct, enum, delegate
            /// </summary>
            T, 
            
            /// <summary>
            /// field
            /// </summary>
            F,
            
            /// <summary>
            /// property (including indexers or other indexed properties)
            /// </summary>
            P, 

            /// <summary>
            /// method (including such special methods as constructors, operators, and so forth)
            /// </summary>
            M,

            /// <summary>
            /// event
            /// </summary>
            E
        }


        /// <summary>
        /// .NET member; corresponds to XML-nodes '/doc/members/member' in the XMLDOC file.
        /// </summary>
        class Member {
            /// <summary>
            /// Not neccessarily unique name of the member (e.g. in the case of overloading), matching the name in the XMLDOC.
            /// </summary>
            public string FullName;

            /// <summary>
            /// Original name; i.e. for nestet types, a plus '+' instead of a dot '.' for dref is used;
            /// </summary>
            public string OrgName;

            /// <summary>
            /// The meta-type of the member.
            /// </summary>
            public KindOfMemeber kom;

            /// <summary>
            /// Corresponding XMLDOC entry.
            /// </summary>
            public XmlNode xmldocEntry;

            /// <summary>
            /// The type, or, for methods, etc. the declaring type.
            /// </summary>
            public Type type;

            /// <summary>
            /// Only used if the member is not a type.
            /// </summary>
            public MemberInfo mi;
        
        }
        
        /// <summary>
        /// All members found in scope.
        /// </summary>
        static List<Member> AllMembers = new List<Member>();


        /// <summary>
        /// Adds all type-members (<see cref="KindOfMemeber.T"/>) to the member list.
        /// </summary>
        static void FindAllMembers_T() {
            

            foreach(AssemblyAndOtherStuff a in assemblies) {
                //Console.WriteLine(assi.FullName);
                var assi = a.assi;

                if(a.xmldoc == null) {
                    //Console.WriteLine("missing XMLDOC for: " + a.assi.GetName().Name);
                    continue;
                }


                Type[] types = assi.GetTypes();
                foreach(var ty in types) {
                    Member T_mem = new Member();
                    T_mem.kom = KindOfMemeber.T;
                    T_mem.OrgName = ty.FullName;
                    if(ty.IsNested) {
                        T_mem.FullName = T_mem.OrgName.Replace('+', '.');
                    } else {
                        T_mem.FullName = T_mem.OrgName;
                    }
                    //a.xmldoc.SelectSingleNode("/doc/members/member[@name='T:BoSSS.Foundation.Basis']");
                    T_mem.xmldocEntry = a.xmldoc.SelectSingleNode("/doc/members/member[@name='T:" + T_mem.FullName + "']");

                    if(T_mem.xmldocEntry != null) {
                        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        // add only members for which we can find an XMLDOC entry;
                        // Otherwise, documentation is not supported.
                        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        AllMembers.Add(T_mem);
                    }


                    FindAllMembers_FPM(ty, T_mem.FullName, a, KindOfMemeber.F, 
                        ty.GetFields(BindingFlags.Instance | BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic));
                    FindAllMembers_FPM(ty, T_mem.FullName, a, KindOfMemeber.P, 
                        ty.GetProperties(BindingFlags.Instance | BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic));
                    FindAllMembers_FPM(ty, T_mem.FullName, a, KindOfMemeber.M, 
                        ty.GetMethods(BindingFlags.Instance | BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic));
                }
            }
        }

        /// <summary>
        /// Adds all class-members (<see cref="KindOfMemeber.F"/>, <see cref="KindOfMemeber.P"/>, <see cref="KindOfMemeber.M"/>) to the member list.
        /// </summary>
        static void FindAllMembers_FPM(Type ty, string type_name, AssemblyAndOtherStuff a, KindOfMemeber K, MemberInfo[] members) {
                        
            foreach(MemberInfo fi in members) {
                Member F_mem = new Member();
                F_mem.kom = K;
                F_mem.OrgName = fi.Name;
                F_mem.FullName = type_name + "." + fi.Name;

                //a.xmldoc.SelectSingleNode("/doc/members/member[@name='T:BoSSS.Foundation.Basis']");
                F_mem.xmldocEntry = a.xmldoc.SelectSingleNode("/doc/members/member[@name='" + K.ToString() + ":" + F_mem.FullName + "']");

                if(F_mem.xmldocEntry != null) {
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // add only members for which we can find an XMLDOC entry;
                    // Otherwise, documentation is not supported.
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    AllMembers.Add(F_mem);
                } else {
                    // skip;
                }

            }

        }



        /// <summary>
        /// Searches for a one of the documented entities (members) in <see cref="AllMembers"/>.
        /// </summary>
        static void SearchMember(string searchString) {
            var foundMembers = AllMembers.Where(member => member.FullName.EndsWith(searchString));

            foreach(var mem in foundMembers) {
                Console.WriteLine("Search for '" + searchString + "': \t"  + mem.FullName);
            }


        }

        private static Assembly MyResolveEventHandler(object sender, ResolveEventArgs args) {
            return assemblies.FirstOrDefault(a => a.assi.FullName.StartsWith(args.Name)).assi;
        }

        static string dir = @"D:\Users\kummer\Documents\BoSSS-got\src\public\L4-application\ipPoisson\bin\Release";
        //static string dir = @"C:\Users\florian\Documents\BoSSS-got\src\public\L4-application\ipPoisson\bin\Release";

        static void Main(string[] args) {
            // load assemblies
            // ===============

            CollectAllDependenciesRecursive(
                Path.Combine(dir, "ipPoisson.exe"),
                assemblies);

            AppDomain currentDomain = AppDomain.CurrentDomain;
            currentDomain.AssemblyResolve += new ResolveEventHandler(MyResolveEventHandler);

            // link all the info together
            // ==========================

            FindAllMembers_T();


            // do something cool
            // =================

            //var n2 = a.xmldoc.SelectSingleNode("/doc/members/member[@name='T:BoSSS.Foundation.Basis']");

            SearchMember("QuotientSource");
            SearchMember("VectorField");
            SearchMember("DGField");

        }
    }
}
