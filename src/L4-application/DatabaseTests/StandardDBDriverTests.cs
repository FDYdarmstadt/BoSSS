/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using NUnit.Framework;
using BoSSS.Platform;

namespace BoSSS.Application.DatabaseTests {

    /// <summary>
    /// Tests for the "standard" Windows file-system.
    /// </summary>
    [TestFixture]
    public class StandardDBDriverTests {

        
        public static void Main() {
            Console.WriteLine("Hello World.");
            InitOnce();

            var tst = new StandardDBDriverTests();
            tst.Init();
           // tst.TestCopySession();
            tst.TestRenameGrid();
            //tst.TestClearDatabase();
            tst.CleanUp();
        }
        

        private IDatabaseInfo m_DB1;

        private IDatabaseInfo m_DB2;

        [TestFixtureSetUp]
        public static void InitOnce() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
        }

        [SetUp]
        public void Init() {
            string templateDatabasePath = "..\\..\\bosss_db_test_template.zip";

            // Temporary filesystem driver to unzip the test database.
            var tmpFsDriver = new StandardFsDriver(templateDatabasePath);

            // database with files in it
            m_DB1 = CreateTestDatabase(Path.Combine(
                Path.GetTempPath(), "BoSSS_DB_" + Path.GetRandomFileName()),
                tmpFsDriver.BasePath);

            // empty database
            m_DB2 = CreateTestDatabase(Path.Combine(
                Path.GetTempPath(), "BoSSS_DB_" + Path.GetRandomFileName()));
        }

        [TearDown]
        public void CleanUp() {
            RemoveTestDatabase(m_DB1);
            RemoveTestDatabase(m_DB2);
        }

        [Test]
        public void TestCopyGrid() {
            var grid1 = m_DB1.Controller.Grids.First();
            var grid2 = m_DB1.Controller.CopyGrid(grid1, m_DB2);
            Assert.IsTrue(grid1.Equals(grid2), "Copied grid is different from its original.");

            var grid3 = m_DB2.Controller.Grids.First();
            Assert.IsTrue(grid1.Equals(grid3), "Copied grid has been saved improperly.");
        }

        [Test]
        public void TestSaveGridIfUnique()
        {
            var gridInfo = m_DB1.Controller.Grids.First();
            IGrid grid;
            if (gridInfo is GridProxy)
            {
                grid = gridInfo.As<GridProxy>().RealGrid;
            }
            else if (gridInfo is Foundation.Grid.Classic.GridCommons)
            {
                grid = (Foundation.Grid.Classic.GridCommons)gridInfo;
            }
            else
            {
                throw new NotSupportedException();
            }

            bool isNotUnique;
            m_DB1.Controller.DBDriver.SaveGridIfUnique(ref grid, out isNotUnique,  m_DB1);
            Assert.IsTrue(isNotUnique == true, "Same grid was not recognized.");

            var grid2 = m_DB1.Controller.CopyGrid(gridInfo, m_DB2);
            m_DB1.Controller.DBDriver.SaveGridIfUnique(ref grid, out isNotUnique, m_DB2);
            Assert.IsTrue(isNotUnique == true, "Copied grid was not recognized.");
        }

        [Test]
        public void TestUnsafelyDeleteGrid() {
            var grid1 = m_DB1.Controller.Grids.First();
            var grid2 = m_DB1.Controller.CopyGrid(grid1, m_DB2);

            m_DB2.Controller.DeleteGrid(grid2, false);

            Assert.IsTrue(m_DB2.Controller.Grids.Count() == 0,
                "Deleted grid still exists in target database.");
        }

        [Test]
        public void TestCopySession() {
            foreach (ISessionInfo session in m_DB1.Controller.Sessions) {
                m_DB1.Controller.CopySession(session, m_DB2);
            }
            IEnumerable<ISessionInfo> db2Sessions = m_DB2.Controller.Sessions;

            Assert.IsTrue(
                db2Sessions.Count() == m_DB1.Controller.Sessions.Count(),
                "Not all sessions copied");
            Assert.IsTrue(
                CountAllFiles(m_DB1.Path) == CountAllFiles(m_DB2.Path),
                "Not all files copied");
        }

        [Test]
        public void TestMoveSession() {
            int preMoveFileCount = CountAllFiles(m_DB1.Path);
            int preMoveSessCount = m_DB1.Controller.Sessions.Count();

            // The number of grids in the original database DB1 determines
            // the number of files left behind in that database.
            int gridCount = m_DB1.Controller.Sessions
                .SelectMany(si => si.GetGrids())
                .Distinct().Count();

            // move all sessions to the other database
            foreach (ISessionInfo session in m_DB1.Controller.Sessions) {
                m_DB1.Controller.MoveSession(session, m_DB2);
            }

            Assert.IsTrue(
                m_DB1.Controller.Sessions.Count() == 0,
                "Not all sessions moved from source database");
            Assert.IsTrue(
                m_DB2.Controller.Sessions.Count() == preMoveSessCount,
                "Not all sessions moved to destination database");
            Assert.IsTrue(
                 CountAllFiles(m_DB1.Path) == CountGridFiles(m_DB2),
                "More than just the grids remaining in source database");
            Assert.IsTrue(
                CountAllFiles(m_DB2.Path) == preMoveFileCount,
                "Not all files moved to destination database");

            // move the sessions back
            foreach (ISessionInfo session in m_DB2.Controller.Sessions) {
                m_DB2.Controller.MoveSession(session, m_DB1);
            }
        }

        [Test]
        public void TestDeleteSession() {
            // copy everything...
            foreach (ISessionInfo session in m_DB1.Controller.Sessions) {
                m_DB1.Controller.CopySession(session, m_DB2);
            }

            // ... and delete it afterwards
            foreach (ISessionInfo session in m_DB2.Controller.Sessions) {
                m_DB2.Controller.DeleteSession(session);
            }


            Assert.IsTrue(
                m_DB2.Controller.Sessions.Count() == 0,
                "Not all sessions deleted.");

            // Only the files associated with grids should remain in the database.
            Assert.IsTrue(
                CountGridFiles(m_DB1) == CountAllFiles(m_DB2.Path),
                "More than just the grids remaining in source database");
        }

        [Test]
        public void TestClearDatabase() {
            Assert.IsTrue(
                CountAllFiles(m_DB2.Path) == 0,
                "Database should be empty in the beginning");

            foreach (ISessionInfo session in m_DB1.Controller.Sessions) {
                m_DB1.Controller.CopySession(session, m_DB2);
            }

            Assert.IsTrue(
                CountAllFiles(m_DB2.Path) == 63,
                "Database should be non-empty after copying");

            m_DB2.Controller.ClearDatabase();
            Assert.IsTrue(
                CountAllFiles(m_DB2.Path) == 0,
                "Files remain in a database that should be empty");
        }

        [Test]
        public void TestSafelyDeleteGrid() {
            m_DB1.Controller.CopySession(m_DB1.Controller.Sessions.First(), m_DB2);

            IGridInfo grd = m_DB2.Controller.Sessions.First().GetGrids().First();
            Assert.IsFalse(
                m_DB2.Controller.DeleteGrid(grd), // Delete operation fails: no files will be changed
                "Deleted grid in use by at least one session.");

            m_DB2.Controller.DeleteSession(m_DB2.Controller.Sessions.First());

            // after the session is deleted, the grid should be successfully removed
            Assert.IsTrue(
                m_DB2.Controller.DeleteGrid(grd),
                "Removal of unused grid unsuccessful.");
            Assert.IsTrue(
                CountAllFiles(m_DB2.Path) == 0,
                "Files remain in a database that should be empty");
        }

        [Test]
        public void TestCleanUpDatabase() {
            ISessionInfo sess = m_DB1.Controller.CopySession(m_DB1.Controller.Sessions.First(), m_DB2);

            int fileCountBefore = CountAllFiles(m_DB2.Path);

            // create a bunch of files in the data subdirectory and clean up the database afterwards
            for (int i = 0; i < 10; i++) {
                string rndFile = Path.Combine(m_DB2.Path, "data", Guid.NewGuid().ToString() + ".1.data");

                if (!File.Exists(rndFile)) {
                    FileStream fs = File.Create(rndFile);
                    fs.Close();
                }
            }
            m_DB2.Controller.CleanDatabase();

            int fileCountAfter = CountAllFiles(m_DB2.Path);

            Assert.IsTrue(
                fileCountBefore == fileCountAfter,
                "Either not all or too many files deleted by clean-up operation");
        }

        [Test]
        [ExpectedException(
            typeof(System.IO.IOException),
            ExpectedMessage = "already exists",
            MatchType = MessageMatch.Contains,
            UserMessage = "Copying the same session twice should not succeed")]
        public void TestCopyFail() {
            foreach (ISessionInfo session in m_DB1.Controller.Sessions) {
                m_DB1.Controller.CopySession(session, m_DB2);
            }

            Random rnd = new Random();

            // this should cause an IOException
            m_DB1.Controller.CopySession(m_DB1.Controller.Sessions.OrderBy(
                (item) => rnd.Next()).First<ISessionInfo>(), m_DB2);
        }

        [Test]
        public void TestRenameSession() {
            string newName = "Something Galerkin something incompressible something something.";
            string newDescription = @"No one will ever read this description.
            Unless this test fails and someone has to debug it.";

            foreach (ISessionInfo session in m_DB1.Controller.Sessions) {
                session.Name = newName;
                session.Description = newDescription;
                Assert.AreEqual(session.Name, newName,
                    "Session name does not match desired value.");
                Assert.AreEqual(session.Description, newDescription,
                    "Session description does not match desired value.");
            }

        }

        [Test]
        public void TestRenameGrid() {
            string newName = "Some-dimensional somehow-ordered grid of order something.";

            foreach (IGridInfo grid in m_DB1.Controller.Grids) {
                grid.Name = newName;
                Assert.AreEqual(grid.Name, newName,
                    "Grid name does not match desired value.");
            }
        }

        /// <summary>
        /// Counts all files in a database.
        /// </summary>
        /// <param name="directoryPath">The path to the database root directory.</param>
        /// <returns>
        /// The total number of files in all subdirectories in the database.
        /// </returns>
        private int CountAllFiles(string directoryPath) {
            return Directory.GetFiles(directoryPath, "*.*", SearchOption.AllDirectories).Length;
        }

        /// <summary>
        /// Counts the number of files that are associated with all grids in a
        /// database. For every grid there should be one .grid and one or more
        /// .data file(s) in the database.
        /// </summary>
        /// <param name="db">The database in question.</param>
        /// <returns>The number of files that belong to the database's grids.</returns>
        private int CountGridFiles(IDatabaseInfo db) {
            int fileCount = db.Controller.Grids.Distinct()
                .SelectMany(grd => db.Controller.GetGridFiles(grd))
                .Count();
            return fileCount;
        }

        /// <summary>
        /// Creates an empty test database from scratch.
        /// </summary>
        private IDatabaseInfo CreateTestDatabase(string basePath) {
            if (Directory.Exists(basePath)) {
                throw new Exception("Database folder already exists");
            }

            Directory.CreateDirectory(basePath);
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.DistVectorDataDir));
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.GridsDir));
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.SessionsDir));
            Directory.CreateDirectory(Path.Combine(basePath, StandardFsDriver.TimestepDir));

            return new DatabaseInfo(basePath);
        }

        /// <summary>
        /// Creates a test database which is the copy of an already existing template.
        /// </summary>
        private IDatabaseInfo CreateTestDatabase(string newBasePath, string templateBasePath) {
            CreateTestDatabase(newBasePath);
            IDatabaseInfo templateDB = new DatabaseInfo(templateBasePath);

            foreach (Guid sessionID in templateDB.Controller.DBDriver.FsDriver.GetAllSessionGUIDs()) {
                Directory.CreateDirectory(Path.Combine(newBasePath,
                    StandardFsDriver.SessionsDir, sessionID.ToString()));
            }

            foreach (string filePath in Directory.GetFiles(templateBasePath, "*",
                SearchOption.AllDirectories)) {
                File.Copy(filePath, filePath.Replace(templateBasePath, newBasePath));
            }

            return new DatabaseInfo(newBasePath);
        }

        private void RemoveTestDatabase(IDatabaseInfo database) {
            Directory.Delete(database.Path, true);
        }
    }
}
