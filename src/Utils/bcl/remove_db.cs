using System;
using System.IO;
using System.Xml;

namespace bcl.remove_db {

    /// <summary>
    /// Subroutine which can be used to remove a registered database from the config of
    /// DBE (DataBase Explorer)
    /// </summary>
    class Program : IProgram {

        //The directory of the database that should be removed from the DBE config 
        private DirectoryInfo databaseDirectory;

        /// <summary>
        /// Decodes and checks the arguments of remove_db. 
        /// </summary>
        /// <param name="args">
        /// Array with exactly one entry.
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
        /// Loads the DBE config XML Document, removes the database if it exists.
        /// </summary>
        public void Execute() {
            //// Check removed since user might want to remove db, _because_
            //// directory does not exist anymore
            //if (!databaseDirectory.Exists) {
            //    throw new UserInputException("The given directory does not exist.");
            //}

            //load DBE config
            string DBEConfigPath = Path.Combine(bcl.myEnv.BOSSS_ETC.FullName, "DBE.xml");
            FileInfo dbeConfigFile = new FileInfo(DBEConfigPath);

            if (!dbeConfigFile.Exists) {
                throw new EnviromentException("DBE config is missing");
            }

            XmlDocument DBEConfig = new XmlDocument();
            try {
                DBEConfig.Load(DBEConfigPath);
            } catch (XmlException e) {
                throw new EnviromentException("Could not load DBE config from" + DBEConfigPath, e);
            }
            
            //Get the Databases Node from config
            XmlNode databasesNode = DBEConfig.SelectSingleNode("DBEControl/Databases");
            if (databasesNode == null) {
                throw new EnviromentException("DBE config is corrupted");
            }

            //Remove Node from Database
            removeDatabases(databasesNode);

            //Save Changes to config
            DBEConfig.Save(DBEConfigPath);
                
        }
        /// <summary>
        /// Removes a database node with a given name. Throws exception when the database node cannot be found.   
        /// </summary>
        /// <param name="databasesNode"></param>
        private void removeDatabases(XmlNode databasesNode) {
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
                    databasesNode.RemoveChild(node);
                    return;
                }
            }

            throw new UserInputException("Given path " + path +" is not registered");
        }

        public void PrintUsage() {
            Console.WriteLine("Usage: bcl remove-db $pathToDatabase");
            Console.WriteLine();
        }
    }
}
