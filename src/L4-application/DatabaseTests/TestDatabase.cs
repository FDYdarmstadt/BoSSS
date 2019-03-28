using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using NUnit.Framework;
using BoSSS.Platform;

namespace BoSSS.Application.DatabaseTests
{
    [TestFixture]
    class TestDatabase
    {
        protected IDatabaseInfo databaseWithFiles;

        protected IDatabaseInfo emptyDatabase;

        [TestFixtureSetUp]
        public static void InitOnce()
        {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
        }

        [SetUp]
        public void Init()
        {
            string templateDatabasePath = "..\\..\\bosss_db_test_template.zip";

            // Temporary filesystem driver to unzip the test database.
            var tmpFsDriver = new StandardFsDriver(templateDatabasePath);

            databaseWithFiles = GetDatabaseCopy(Path.Combine(
                Path.GetTempPath(), "BoSSS_DB_" + Path.GetRandomFileName()),
                tmpFsDriver.BasePath);

            emptyDatabase = CreateEmptyDatabase(Path.Combine(
                Path.GetTempPath(), "BoSSS_DB_" + Path.GetRandomFileName()));
        }

        [TearDown]
        public void CleanUp()
        {
            RemoveTestDatabase(databaseWithFiles);
            RemoveTestDatabase(emptyDatabase);
        }

        private IDatabaseInfo CreateEmptyDatabase(string basePath)
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

            return new DatabaseInfo(basePath);
        }

        private IDatabaseInfo GetDatabaseCopy(string newBasePath, string templateBasePath)
        {
            CreateEmptyDatabase(newBasePath);
            IDatabaseInfo templateDB = new DatabaseInfo(templateBasePath);

            foreach (Guid sessionID in templateDB.Controller.DBDriver.FsDriver.GetAllSessionGUIDs())
            {
                Directory.CreateDirectory(Path.Combine(newBasePath,
                    StandardFsDriver.SessionsDir, sessionID.ToString()));
            }

            foreach (string filePath in Directory.GetFiles(templateBasePath, "*",
                SearchOption.AllDirectories))
            {
                File.Copy(filePath, filePath.Replace(templateBasePath, newBasePath));
            }

            return new DatabaseInfo(newBasePath);
        }

        private void RemoveTestDatabase(IDatabaseInfo database)
        {
            Directory.Delete(database.Path, true);
        }
    }
}
