using System;
using System.Diagnostics;
using System.IO;
using System.Xml;

namespace bcl {

    /// <summary>
    /// Switch for different ways to notate a filesystem-path
    /// </summary>
    enum PathKind {

        /// <summary>
        /// Windows-style
        /// </summary>
        win,

        /// <summary>
        /// Unix/Linux style
        /// </summary>
        unix
    }

    /// <summary>
    /// Main class of the BoSSS command line tools
    /// </summary>
    static class bcl {

        /// <summary>
        /// utility function to extract a path from a xml path node.
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        static public string BrowsePath(XmlNode p) {
            if (p.Attributes.Count != 1)
                throw new ApplicationException("wrong number of attributes");

            // extract value
            // =============
            string value = null;
            {
                foreach (XmlAttribute att in p.Attributes) {
                    if (att.Name.Equals("value"))
                        value = att.Value;
                    else
                        throw new ApplicationException("unknown attribute;");
                }
                if (value == null)
                    throw new ApplicationException("missing value attribute;");
            }

            // browse
            // ======
            return BrowsePath(value);

        }

        /// <summary>
        /// utility function to extract a path from a string:
        /// Unix or Windows path separators are detected and replaced by the separator of the current system,
        /// variables like '%BOSSS_ROOT%' are replaced by the boss root path of the current installation.
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        static public string BrowsePath(string value) {

            // check path
            // ==========
            if (value.Contains("\r") || value.Contains("\n") || value.Contains("\t"))
                throw new ApplicationException("path value contains illegal character");

            // autodetect for Unix/Windows - path
            // ==================================

            PathKind pk = PathKind.win;
            if (value.Contains("/"))
                pk = PathKind.unix;

            // convert to the current os
            // =========================
            if (Path.DirectorySeparatorChar == '/' && pk == PathKind.win) {
                // convert win->unix
                value = value.Replace('\\', '/');
            }
            if (Path.DirectorySeparatorChar == '\\' && pk == PathKind.unix) {
                // convert win->unix
                value = value.Replace('/', '\\');
            }

            // replace variables
            // =================

            value = value.Replace("%BOSSS_ROOT%", myEnv.BOSSS_ROOT.FullName);
            value = value.Replace("%BOSSS_ETC%", myEnv.BOSSS_ETC.FullName);

            value = value.Replace("$BOSSS_ROOT", myEnv.BOSSS_ROOT.FullName);
            value = value.Replace("$BOSSS_ETC", myEnv.BOSSS_ETC.FullName);

            // return
            // ======
            return value;
        }

        /// <summary>
        /// The <see cref="MyEnvironment"/> that has been determined in
        /// <see cref="Main"/>
        /// </summary>
        static public MyEnvironment myEnv;

        /// <summary>
        /// performs the subprogram-decoding;
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            myEnv = new MyEnvironment();
            if (args.Length < 1) {
                PrintUsage();
                return;
            }

            string command = args[0];
            string[] cmdrest = new string[args.Length - 1];
            if (cmdrest.Length > 0) {
                Array.Copy(args, 1, cmdrest, 0, cmdrest.Length);
            }

            IProgram p = null;
            try {
                p = ProgramFactory(command);

                if (cmdrest.Length == 1 && cmdrest[0].Equals("?")) {
                    PrintUsage(p);
                } else {
                    p.DecodeArgs(cmdrest);
                    p.Execute();
                }
            } catch (UserInputException uie) {
                Console.WriteLine("ERROR: bad user input: " + uie.Message);
                Console.WriteLine();
                if (p == null) {
                    PrintUsage();
                } else {
                    PrintUsage(p);
                }
                Environment.Exit(-1);
            } catch (EnviromentException ee) {
                Console.WriteLine("ERROR: BoSSS environment: " + ee.Message);
                Console.WriteLine();
                Environment.Exit(-2);
            }
        }

        /// <summary>
        /// Instantiates the correct <see cref="IProgram"/> depending on the
        /// entered <paramref name="command"/>.
        /// </summary>
        /// <param name="command">The comand/subroutine to be executed</param>
        /// <returns>
        /// An instance of a class implementing <see cref="IProgram"/>
        /// </returns>
        private static IProgram ProgramFactory(string command) {
            IProgram p;
            switch (command) {
                case "init-db":
                    // Create empty BoSSS database scheme and register it
                    p = new init_db.Program();
                    break;

                case "register-db":
                    // Register a BoSSS database in the DBE config
                    p = new register_db.Program();
                    break;

                case "config":
                    // Modify bcl config parameters
                    p = new config.Program();
                    break;

                case "deploy-at":
                    // Deploy a BoSSS application including all dependencies
                    p = new deploy_at.Program();
                    break;

                case "nativelib-inst":
                    p = new nativelib_inst.Program();
                    break;

                case "create-proj":
                    p = new create_proj.Program();
                    break;

                case "visualizers-inst":
                    p = new visualizers_inst.Program();
                    break;

                case "remove-db":
                    p = new remove_db.Program();
                    break;

                default:
                    throw new UserInputException("Invalid command");
            }
            return p;
        }

        /// <summary>
        /// Just prints a disclaimer and some info about how to use this tool
        /// </summary>
        private static void PrintDisclaimer() {
            Console.WriteLine("BoSSS command line tools.");
            Console.WriteLine("(c) 2010 FDY");
        }

        /// <summary>
        /// Prints a combination of the disclaimer (see
        /// <see cref="PrintDisclaimer"/>) and some usage information about
        /// available subroutines
        /// </summary>
        private static void PrintUsage() {
            Console.WriteLine("Usage: bcl $subroutine [$argument1 [$argument2 ...]]");
            Console.WriteLine();
            Console.WriteLine("The following subroutines are available:");
            Console.WriteLine(" config            Consistently edits the configuration in $root/etc/config.xml");
            Console.WriteLine(" deploy-at         Deploys binaries including dependencies at a desired location");
            Console.WriteLine(" init-db           Inits and registers ('bcl register-db') a new BoSSS database");
            Console.WriteLine(" push-all          Executes 'push all' for all git repositories");
            Console.WriteLine(" register-db       adds a database to the BOSSS_ROOT/etc/DBE.xml file.");
            Console.WriteLine(" nativelib-inst    deploys the native binaries in all subfolders of BOSSS_ROOT ");
            Console.WriteLine("                   (where required).");
            Console.WriteLine(" create-proj       creates a new BoSSS project from a template and adds it to ");
            Console.WriteLine("                   'MySolution.sln'.");
            Console.WriteLine(" visualizers-inst  Installs BoSSS visualizers for supported Visual Studio ");
            Console.WriteLine("                   versions (currently, VS2015).");
            Console.WriteLine(" remove-db         Removes a registered database from the BOSSS_ROOT/etc/DBE.xml file.");
            PrintDisclaimer();

            myEnv.PrintInfo();
        }

        /// <summary>
        /// Prints a combination of the disclaimer (see
        /// <see cref="PrintDisclaimer"/>) and some usage information for the
        /// subroutine (see <see cref="IProgram.PrintUsage"/>).
        /// </summary>
        /// <param name="p">The affected subroutine</param>
        private static void PrintUsage(IProgram p) {
            p.PrintUsage();
            PrintDisclaimer();
        }

        /// <summary>
        /// used by <see cref="MyEnvironment.InitAndVerify"/>
        /// </summary>
        /// <param name="root"></param>
        /// <param name="msg"></param>
        /// <returns></returns>
        static bool VerifyRootDir(DirectoryInfo root, out string msg) {
            // test for "etc"
            if (root.GetDirectories("etc").Length != 1) {
                msg = "missing \"etc\" - directory;";
                return false;
            }

            msg = "";
            return true;
        }
    }
}
