using System;
using System.IO;
using System.Xml;
using CWDev.SLNTools.Core;
//create visual studio solution from using Microsoft.Build.BuildEngine;

namespace bcl.create_proj {

    /// <summary>
    /// Creates a new BoSSS L4 project (i.e. a C# - project)
    /// from the src/public/L4-Application/BoSSSTemplateProject
    /// </summary>
    class Program : IProgram {

        enum VisualStudioVersions {
            
            V15
        }

        #region IProgram Members

        public void Execute() {
            // source and dest dir
            string srcDir = Path.Combine(
                bcl.myEnv.BOSSS_SRC.FullName,
                Path.Combine("public", "src",  "L4-application", "BoSSSTemplateProject"));

            string destDir;
            if (Path.IsPathRooted(m_CreationPath)) {
                destDir = m_CreationPath;
            } else {
                destDir = Path.Combine(Directory.GetCurrentDirectory(), m_CreationPath);
            }
            destDir = Path.Combine(destDir, m_newName);
            Console.WriteLine("creating in :" + destDir);


            DirectoryInfo _srcDir = new DirectoryInfo(srcDir);
            DirectoryInfo _destDir = new DirectoryInfo(destDir);

            // create destination dir
            if (_destDir.Exists)
                throw new ApplicationException("project directory '" + destDir + "' already exists.");
            _destDir.Create();

            // copy files
            // ==========
            File.Copy(Path.Combine(srcDir, "TemplateMain.cs"), Path.Combine(destDir, m_newName + "Main.cs"));
            File.Copy(Path.Combine(srcDir, "app.config"), Path.Combine(destDir, "app.config"));
            Directory.CreateDirectory(Path.Combine(destDir, "Properties"));
            File.Copy(Path.Combine(srcDir, Path.Combine("Properties", "AssemblyInfo.cs")), Path.Combine(destDir, Path.Combine("Properties", "AssemblyInfo.cs")));

            // project file
            // ============
            Guid newProjGuid = Guid.NewGuid();
            {
                XmlDocument projFile = new XmlDocument();
                projFile.Load(Path.Combine(srcDir, "BoSSSTemplateProject.csproj"));
                XmlNamespaceManager scheissDing = new XmlNamespaceManager(projFile.NameTable);
                scheissDing.AddNamespace("haß", "http://schemas.microsoft.com/developer/msbuild/2003");

                // guid
                XmlNode proj_guid = projFile.SelectSingleNode("/haß:Project/haß:PropertyGroup/haß:ProjectGuid", scheissDing);
                //Console.WriteLine(proj_guid.InnerText);
                proj_guid.InnerText = newProjGuid.ToString("B");
                //Console.WriteLine(proj_guid.InnerText);

                // RootNamespace
                XmlNode root_Namesp = projFile.SelectSingleNode("/haß:Project/haß:PropertyGroup/haß:RootNamespace", scheissDing);
                root_Namesp.InnerText = m_newName;

                // AssemblyName
                XmlNode AssemblyName = projFile.SelectSingleNode("/haß:Project/haß:PropertyGroup/haß:AssemblyName", scheissDing);
                AssemblyName.InnerText = m_newName;

                // MainFile
                var CsFiles = projFile.SelectNodes("/haß:Project/haß:ItemGroup/haß:Compile/@Include", scheissDing);
                foreach (var f in CsFiles) {
                    XmlAttribute att = f as XmlAttribute;
                    if (att.Value == "TemplateMain.cs")
                        att.Value = m_newName + "Main.cs";
                }

                // ProjectReferences
                var ProjRefs = projFile.SelectNodes("/haß:Project/haß:ItemGroup/haß:ProjectReference/@Include", scheissDing);
                foreach (var f in ProjRefs) {
                    XmlAttribute att = f as XmlAttribute;

                    FileInfo item = new FileInfo(Path.Combine(srcDir, att.Value));
                    //Console.WriteLine(item.FullName + ", " + item.Exists);

                    string newPath = GetRelPath(item, _destDir);
                    //Console.WriteLine(item.FullName + ", " + File.Exists(newPath));
                    att.Value = newPath;
                }
                // References
                var Refs = projFile.SelectNodes("/haß:Project/haß:ItemGroup/haß:Reference/haß:HintPath", scheissDing);
                foreach (var f in Refs) {
                    var n = f as XmlNode;

                    FileInfo item = new FileInfo(Path.Combine(srcDir, n.InnerText));
                    //Console.WriteLine(item.FullName + ", " + item.Exists);

                    string newPath = GetRelPath(item, _destDir);
                    //Console.WriteLine(item.FullName + ", " + File.Exists(newPath));
                    n.InnerText = newPath;
                }

                // save
                projFile.Save(Path.Combine(destDir, m_newName + ".csproj"));
            }

            // solution file
            // =============
            {
                string newSolFile = Path.Combine(destDir, m_newName + "_WithDeps.sln");

                string templateSolFile;
                switch (m_VSversion) {
                    case VisualStudioVersions.V15:
                        templateSolFile = "BoSSSTemplateProject_WithDeps.sln";
                        break;

                    default:
                        throw new NotImplementedException();
                }
                templateSolFile = Path.Combine(srcDir, templateSolFile);

                var solution = SolutionFile.FromFile(templateSolFile);

                Project template = null;
                foreach (var proj in solution.Projects) {
                    if (proj.ProjectName == "BoSSSTemplateProject") {
                        //proj.ProjectName = newName;
                        //proj.RelativePath = newName + ".csproj";
                        template = proj;
                    } else {
                        FileInfo prjFile = new FileInfo(proj.FullPath);
                        string newPath = GetRelPath(prjFile, _destDir);
                        proj.RelativePath = newPath;
                    }
                }

                if (template == null)
                    throw new ApplicationException("unable to find template project.");

                var pNew = new Project(solution,
                    newProjGuid.ToString("B"),
                    template.ProjectGuid,
                    m_newName,
                    m_newName + ".csproj",
                    template.ParentFolderGuid,
                    template.ProjectSections,
                    template.VersionControlLines,
                    template.ProjectConfigurationPlatformsLines);

                solution.Projects.Remove(template);
                solution.Projects.Add(pNew);

                solution.SaveAs(newSolFile);
            }

        }

        /// <summary>
        /// computes the relative path from location <paramref name="loc"/> to file <paramref name="f"/>;
        /// </summary>
        string GetRelPath(FileInfo f, DirectoryInfo loc) {

            string[] dirs_f = f.Directory.FullName.Split(Path.DirectorySeparatorChar);
            string[] dirs_loc = loc.FullName.Split(Path.DirectorySeparatorChar);

            int equals = 0;
            while (equals < dirs_loc.Length && equals < dirs_f.Length && dirs_f[equals] == dirs_loc[equals])
                equals++;

            string res = "";
            for (int i = dirs_loc.Length; i > equals; i--)
                res = Path.Combine(res, "..");

            for (int i = equals; i < dirs_f.Length; i++) {
                res = Path.Combine(res, dirs_f[i]);
            }

            res = Path.Combine(res, f.Name);
            return res;
        }

        /// <summary>
        /// path where the new project will be created
        /// </summary>
        string m_CreationPath;

        /// <summary>
        /// name for the project to create
        /// </summary>
        string m_newName;

        VisualStudioVersions m_VSversion = VisualStudioVersions.V15;

        public void DecodeArgs(string[] args) {
            if (args.Length < 1)
                throw new UserInputException("Missing option;");
            if (args.Length > 2)
                throw new UserInputException("Too many arguments.");

            m_newName = args[0];

            if (m_newName.IndexOfAny(Path.GetInvalidFileNameChars()) > 0 || m_newName.IndexOfAny(Path.GetInvalidPathChars()) > 0) {
                string invalids = "";
                foreach (char s in Path.GetInvalidFileNameChars())
                    invalids += ("'" + s + "' ");
                foreach (char s in Path.GetInvalidPathChars())
                    invalids += ("'" + s + "' ");
                throw new UserInputException("projName may not contain any of " + invalids + ";");
            }

            if (m_newName.IndexOfAny(new char[] { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' }) == 0)
                throw new UserInputException("projName may not start with a number;");

            if (args.Length > 1) {
                m_CreationPath = args[1];
            } else {
                m_CreationPath = ".";
            }
        }

        public void PrintUsage() {
            Console.WriteLine("Usage: bcl create-proj $projName {$location} ");
            Console.WriteLine();
            Console.WriteLine("The following are available:");
            Console.WriteLine(" projName       name of the new project");
            Console.WriteLine(" location       Optional location;");
            Console.WriteLine("                per Default, the project is created in the");
            Console.WriteLine("                current directory.");
            Console.WriteLine();
        }

        #endregion
    }
}
