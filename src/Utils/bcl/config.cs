using System;
using System.IO;
using System.Reflection;
using System.Xml;

namespace bcl.config {

    /// <summary>
    /// Utility program to read/edit the bcl-configuration file in a consistent
    /// manner (without the need to edit the xml-file manually)
    /// </summary>
    public class Program : IProgram {

        /// <summary>
        /// The filename of the config file
        /// </summary>
        private const string CONFIG_FILE_NAME = "bcl-config.xml";

        /// <summary>
        /// The name of the root node
        /// </summary>
        private const string CONFIG_ROOT_NODE = "bclConfig";

        /// <summary>
        /// A list of options that can be handled by this program. To extend
        /// the scope of this tool, add another option and write handling
        /// functions (<see cref="Execute"/>)
        /// </summary>
        private readonly string[] validOptionNames =
            new string[] { "InstallationSource", "InstallationType" };

        /// <summary>
        /// The name of the configuration option to be read/edited
        /// </summary>
        private string optionName;

        /// <summary>
        /// The new value that should be assigned to <see cref="optionName"/>
        /// (will be null if the option should just be read)
        /// </summary>
        private string optionValue = null;

        /// <summary>
        /// The path to the config file which is
        /// <see cref="MyEnvironment.BOSSS_ETC"/>/<see cref="CONFIG_FILE_NAME"/>
        /// </summary>
        private string PathToConfigFile {
            get {
                return Path.Combine(bcl.myEnv.BOSSS_ETC.FullName, CONFIG_FILE_NAME);
            }
        }

        /// <summary>
        /// Retrieves the installationSource from the config file
        /// </summary>
        /// <returns>A string containing the configured setting</returns>
        public string GetInstallationSource() {
            XmlDocument config = GetConfig();
            XmlNode rootNode = GetRootNode(config);
            XmlNode installationSourceNode = rootNode.SelectSingleNode("installationSource");
            if (installationSourceNode == null) {
                return "";
            } else {
                return installationSourceNode.InnerText;
            }
        }

        /// <summary>
        /// Changes the path to the installation source while making
        /// sure that the location really exists
        /// </summary>
        public void SetInstallationSource() {
            // Check the path for existence
            DirectoryInfo source = new DirectoryInfo(optionValue);
            if (!source.Exists) {
                throw new UserInputException("The given path \"" + optionValue + "\" does not exist");
            }

            // Check for self reference (the TrimEnd ensures that c:\bla and
            // c:\bla\ are considered equal
            string pathToSource = source.FullName.TrimEnd('\\');
            if (pathToSource.Equals(bcl.myEnv.BOSSS_ROOT.FullName.TrimEnd('\\'))) {
                throw new UserInputException("An installation cannot be it's own source");
            }

            XmlDocument config = GetConfig();
            XmlNode rootNode = GetRootNode(config);
            XmlNode installationSourceNode = rootNode.SelectSingleNode("installationSource");

            // Modify the xml structure
            if (installationSourceNode == null) {
                installationSourceNode = config.CreateElement("installationSource");
                rootNode.AppendChild(installationSourceNode);
            } else {
                installationSourceNode.RemoveAll();
            }
            installationSourceNode.AppendChild(config.CreateTextNode(pathToSource));

            config.Save(PathToConfigFile);
        }

        /// <summary>
        /// Reads the "installationType" option from the bcl config
        /// </summary>
        /// <returns>The value stored in xml config</returns>
        public string GetInstallationType() {
            XmlDocument config = GetConfig();
            XmlNode rootNode = GetRootNode(config);
            XmlNode installationTypeNode = rootNode.SelectSingleNode("installationType");
            if (installationTypeNode == null) {
                return "";
            } else {
                return installationTypeNode.InnerText;
            }
        }

        /// <summary>
        /// Sets the "installationTyp" option in the bcl config. Valid values
        /// for <see cref="optionValue"/> are "developer" and "user"
        /// </summary>
        public void SetInstallationType() {
            string installationType;
            switch (optionValue) {
                case "developer":
                case "user":
                    installationType = optionValue;
                    break;

                default:
                    throw new UserInputException("Invalid installation type \"" + optionValue + "\"");
            }

            XmlDocument config = GetConfig();
            XmlNode rootNode = GetRootNode(config);
            XmlNode installationTypeNode = rootNode.SelectSingleNode("installationType");

            if (installationTypeNode == null) {
                installationTypeNode = config.CreateElement("installationType");
                rootNode.AppendChild(installationTypeNode);
            } else {
                installationTypeNode.RemoveAll();
            }
            installationTypeNode.AppendChild(config.CreateTextNode(installationType));

            config.Save(PathToConfigFile);
        }

        /// <summary>
        /// Create an empty configuration file
        /// Note: Overwrites the original file if it exists
        /// </summary>
        private XmlDocument InitConfigFile() {
            DirectoryInfo etcFolder = bcl.myEnv.BOSSS_ETC;
            if (!etcFolder.Exists) {
                etcFolder.Create();
            }

            XmlDocument bclConfig = new XmlDocument();
            bclConfig.AppendChild(bclConfig.CreateXmlDeclaration("1.0", null, null));
            XmlElement root = bclConfig.CreateElement(CONFIG_ROOT_NODE);
            bclConfig.AppendChild(root);

            try {
                bclConfig.Save(PathToConfigFile);
            } catch (XmlException e) {
                throw new Exception("Could not write config file for bcl in " + bcl.myEnv.BOSSS_ETC, e);
            }

            return bclConfig;
        }

        /// <summary>
        /// Tries to read the config file from <see cref="PathToConfigFile"/>
        /// </summary>
        /// <returns>An xml tree containing the current config</returns>
        private XmlDocument GetConfig() {
            FileInfo configFile = new FileInfo(PathToConfigFile);
            XmlDocument config;

            // If no config file can be found, create a fresh one
            if (configFile.Exists) {
                config = new XmlDocument();
                config.Load(PathToConfigFile);
            } else {
                config = InitConfigFile();
            }

            return config;
        }

        /// <summary>
        /// Safely extracts the root node from the config
        /// </summary>
        /// <param name="config">The xml-tree of the config file</param>
        /// <returns>The root node of the config file</returns>
        private XmlNode GetRootNode(XmlDocument config) {
            XmlNode node = config.SelectSingleNode(CONFIG_ROOT_NODE);
            if (node == null) {
                throw new Exception("Invalid config file format");
            }
            return node;
        }

        #region IProgram Member

        /// <summary>
        /// Dynamically determine the handler method for the option stored in
        /// <see cref="optionName"/>. There have to two be a handler method for
        /// every valid option stored in <see cref="validOptionNames"/>. The
        /// naming convetion a method handling the option named "DummyOption"
        /// would be GetDummyOption() for retrieving and SetDummyOption() for
        /// setting a value.
        /// </summary>
        /// <example>
        /// The configuration option "InstallationType" is being handled by
        /// the methods <see cref="GetInstallationType"/> and
        /// <see cref="SetInstallationType"/>
        /// </example>
        public void Execute() {
            if (optionValue == null) {
                // Retrieve mode (get|set)
                MethodInfo method = this.GetType().GetMethod("Get" + optionName);
                if (method == null) {
                    throw new Exception("Internal error: No retrieve-handler for the valid configuration option \"" + optionName + "\" could be found");
                }
                
                if (!method.ReturnType.Equals(Type.GetType("System.String"))) {
                    throw new Exception("Internal error: Handling function for retrieving the valid configuration option \"" + optionName + "\" is invalid");
                }

                try {
                    string result = (string)method.Invoke(this, null);
                    Console.WriteLine(result);
                } catch (TargetInvocationException e) {
                    throw e.InnerException;
                }
            } else {
                // Set mode
                MethodInfo method = this.GetType().GetMethod("Set" + optionName);
                if (method == null) {
                    throw new Exception("Internal error: No set handling for the valid configuration option \"" + optionName + "\" could be found");
                }

                try {
                    method.Invoke(this, null);
                } catch (TargetInvocationException e) {
                    throw e.InnerException;
                }
            }
        }

        /// <summary>
        /// Decodes the command line arguments
        /// Argument 1: The name of the option to be modified or read
        /// Argument 2 (optional): The new value assigned to this option
        /// </summary>
        /// <param name="args">The arguments to be decoded</param>
        public void DecodeArgs(string[] args) {
            if (args.Length < 1) {
                throw new UserInputException("Missing option");
            }
            if (args.Length > 2) {
                throw new UserInputException("Too many arguments");
            }

            optionName = args[0].ToLowerInvariant();
            bool optionNameIsValid = false;
            foreach (string name in validOptionNames) {
                if (name.Equals(optionName, StringComparison.InvariantCultureIgnoreCase)) {
                    optionNameIsValid = true;

                    // Overwrite option name so we don't have to mess with
                    // cases anymore (simply always use the one stored)
                    optionName = name;
                    break;
                }
            }
            if (!optionNameIsValid) {
                throw new UserInputException("Invalid option \"" + optionName + "\"");
            }

            if (args.Length == 2) {
                optionValue = args[1];
            }
        }

        /// <summary>
        /// See <see cref="IProgram.PrintUsage"/>.
        /// </summary>
        public void PrintUsage() {
            Console.WriteLine("Usage: bcl config $optionName [$value]");
            Console.WriteLine();
            Console.WriteLine("If no $value is specified, returns the value stored in $root/etc/config.xml for $optionName. If $value is non-empty, overwrites the value while imposing option-dependent sanity checks.");
            Console.WriteLine();
        }

        #endregion
    }
}
