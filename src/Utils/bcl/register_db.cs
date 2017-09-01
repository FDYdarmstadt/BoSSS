using System;
using System.IO;
using System.Xml;
using System.Diagnostics;
using bcl.Properties;

namespace bcl.register_db {

    /// <summary>
    /// Subroutine which can be used to register a new database in the config of
    /// the DBE (DataBase Explorer).
    /// </summary>
    class Program : IProgram {

        /// <summary>
        /// The directory that should be registered.
        /// </summary>
        DirectoryInfo databaseDirectory;

        /// <summary>
        /// Checks if the given directory has the correct directory structure
        /// for a s BoSSS database (see <see cref="init_db.Program.Execute"/>).
        /// </summary>
        /// <param name="databaseDirectory">The directory in question</param>
        /// <returns>True if the structure is valid, false otherwise</returns>
        private bool IsValidBoSSSDatabase(DirectoryInfo databaseDirectory) {
            int score = 0;
            foreach (DirectoryInfo directory in databaseDirectory.GetDirectories()) {
                switch (directory.Name) {
                    case "data":
                    case "timesteps":
                    case "grids":
                    case "sessions":
                        score++;
                        break;
                }
            }

            return (score == 4);
        }

        /// <summary>
        /// Checks whether the DBE config already contains a database entry for
        /// <see cref="databaseDirectory"/>.
        /// </summary>
        /// <param name="databasesNode">
        /// The node which contains the database entries
        /// </param>
        /// <returns>
        /// True if config already contains the path, false otherwise
        /// </returns>
        private bool DatabaseIsAlreadyRegistered(XmlNode databasesNode) {
            string path = databaseDirectory.FullName.TrimEnd('\\');

            foreach (XmlNode node in databasesNode.ChildNodes) {
                XmlNode pathNode = node.SelectSingleNode("path");
                if (pathNode == null) {
                    continue;
                }

                XmlNode pathValueNode = pathNode.Attributes.GetNamedItem("value");
                if (pathValueNode == null) {
                    continue;
                }

                if (path.Equals(Path.GetFullPath(pathValueNode.InnerText).TrimEnd('\\'))) {
                    return true;
                }
            }

            return false;
        }

        #region IProgram Member

        /// <summary>
        /// Looks for the given path, checks if it is valid BoSSS database and
        /// adds it to the DBE config (if it does not already exist)
        /// </summary>
        public void Execute() {
            // Check existence and validity
            if (!databaseDirectory.Exists) {
                throw new UserInputException("The given directory does not exist.");
            }

            if (!IsValidBoSSSDatabase(databaseDirectory)) {
                throw new UserInputException("The given directory is not a valid BoSSS database");
            }

            string DBEConfigPath = Path.Combine(bcl.myEnv.BOSSS_ETC.FullName, "DBE.xml");
            FileInfo dbeConfigFile = new FileInfo(DBEConfigPath);

            XmlDocument DBEConfig;
            if (dbeConfigFile.Exists) {
                DBEConfig = new XmlDocument();

                try {
                    DBEConfig.Load(DBEConfigPath);
                } catch (FileNotFoundException e) {
                    throw new UserInputException("Could not load DBE config from " + DBEConfigPath, e);
                }
            } else {
                DBEConfig = new XmlDocument();
                DBEConfig.LoadXml(Resources.dbe_default_xml);
                DBEConfig.Save(DBEConfigPath);
            }

            // Get Databases node from the config
            XmlNode databasesNode = DBEConfig.SelectSingleNode("DBEControl/Databases");
            if (databasesNode == null) {
                throw new EnviromentException("DBE config is corrupted");
            }
            if (DatabaseIsAlreadyRegistered(databasesNode)) {
                throw new UserInputException("Given path " + databaseDirectory.FullName + " is already registered.");
            }

            // Add the path
            XmlNode databaseNode = DBEConfig.CreateElement("Database");
            databasesNode.AppendChild(databaseNode);
            XmlNode pathNode = DBEConfig.CreateElement("path");
            databaseNode.AppendChild(pathNode);
            XmlAttribute valueAttribute = DBEConfig.CreateAttribute("value");
            valueAttribute.InnerText = databaseDirectory.FullName;
            pathNode.Attributes.Append(valueAttribute);
            DBEConfig.Save(DBEConfigPath);
        }

        /// <summary>
        /// Checks the arguments for the directory to be registered
        /// </summary>
        /// <param name="args">
        /// An array with exactly one entry containing the path
        /// </param>
        public void DecodeArgs(string[] args) {
            if (args.Length > 1) {
                throw new UserInputException("Too many arguments");
            } else if (args.Length == 0) {
                throw new UserInputException("Missing argument: database to be registered");
            } else {
                databaseDirectory = new DirectoryInfo(args[0]);
            }
        }

        /// <summary>
        /// See <see cref="IProgram.PrintUsage"/>.
        /// </summary>
        public void PrintUsage() {
            Console.WriteLine("Usage: bcl register-db $pathToDatabase");
            Console.WriteLine();
        }

        #endregion
    }
}
