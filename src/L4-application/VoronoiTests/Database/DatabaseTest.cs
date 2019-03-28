using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.IO;
using NUnit.Framework;

namespace VoronoiTests.Database
{
    class DatabaseTest : RunnableTest
    {
        public static IDatabaseInfo Database {
            get {
                return database ?? (database = InitializeDatabase());
            }
        }

        static IDatabaseInfo database;

        public static IDatabaseInfo InitializeDatabase()
        {
            string databasePath = "..\\..\\bosss_TestDatabase_Voronoi";
            IDatabaseInfo emptyDatabase = CreateEmptyDatabase(databasePath);
            return emptyDatabase;
        }

        static IDatabaseInfo CreateEmptyDatabase(string basePath)
        {
            if (Directory.Exists(basePath))
            {
                throw new Exception("Database folder already exists");
            }

            Directory.CreateDirectory(basePath);
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.DistVectorDataDir));
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.GridsDir));
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.SessionsDir));
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.TimestepDir));
            IDatabaseInfo emptyDatabase = new DatabaseInfo(basePath);
            return emptyDatabase;
        }

        [TestFixtureTearDown]
        public static void DeleteDatabase()
        {
            Directory.Delete(database.Path, true);
        } 
    }
}
