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

namespace BoSSS.Application.DatabaseTests  {

    /// <summary>
    /// Tests for the "standard" Windows file-system.
    /// </summary>
    /// 
    class DBDriverTests : DatabaseTest
    {
        [Test]
        public void TestCopyGrid() {
            var grid1 = databaseWithFiles.Controller.Grids.First();
            var grid2 = databaseWithFiles.Controller.CopyGrid(grid1, emptyDatabase);
            Assert.IsTrue(grid1.Equals(grid2), "Copied grid is different from its original.");

            var grid3 = emptyDatabase.Controller.Grids.First();
            Assert.IsTrue(grid1.Equals(grid3), "Copied grid has been saved improperly.");
        }

        [Test]
        public void TestUnsafelyDeleteGrid() {
            var grid1 = databaseWithFiles.Controller.Grids.First();
            var grid2 = databaseWithFiles.Controller.CopyGrid(grid1, emptyDatabase);

            emptyDatabase.Controller.DeleteGrid(grid2, false);

            Assert.IsTrue(emptyDatabase.Controller.Grids.Count() == 0,
                "Deleted grid still exists in target database.");
        }

        [Test]
        public void TestCopySession() {
            foreach (ISessionInfo session in databaseWithFiles.Controller.Sessions) {
                databaseWithFiles.Controller.CopySession(session, emptyDatabase);
            }
            IEnumerable<ISessionInfo> db2Sessions = emptyDatabase.Controller.Sessions;

            Assert.IsTrue(
                db2Sessions.Count() == databaseWithFiles.Controller.Sessions.Count(),
                "Not all sessions copied");
            Assert.IsTrue(
                CountAllFiles(databaseWithFiles.Path) == CountAllFiles(emptyDatabase.Path),
                "Not all files copied");
        }

        [Test]
        public void TestMoveSession() {
            int preMoveFileCount = CountAllFiles(databaseWithFiles.Path);
            int preMoveSessCount = databaseWithFiles.Controller.Sessions.Count();

            // The number of grids in the original database DB1 determines
            // the number of files left behind in that database.
            int gridCount = databaseWithFiles.Controller.Sessions
                .SelectMany(si => si.GetGrids())
                .Distinct().Count();

            // move all sessions to the other database
            foreach (ISessionInfo session in databaseWithFiles.Controller.Sessions) {
                databaseWithFiles.Controller.MoveSession(session, emptyDatabase);
            }

            Assert.IsTrue(
                databaseWithFiles.Controller.Sessions.Count() == 0,
                "Not all sessions moved from source database");
            Assert.IsTrue(
                emptyDatabase.Controller.Sessions.Count() == preMoveSessCount,
                "Not all sessions moved to destination database");
            Assert.IsTrue(
                 CountAllFiles(databaseWithFiles.Path) == CountGridFiles(emptyDatabase),
                "More than just the grids remaining in source database");
            Assert.IsTrue(
                CountAllFiles(emptyDatabase.Path) == preMoveFileCount,
                "Not all files moved to destination database");

            // move the sessions back
            foreach (ISessionInfo session in emptyDatabase.Controller.Sessions) {
                emptyDatabase.Controller.MoveSession(session, databaseWithFiles);
            }
        }

        [Test]
        public void TestDeleteSession() {
            // copy everything...
            foreach (ISessionInfo session in databaseWithFiles.Controller.Sessions) {
                databaseWithFiles.Controller.CopySession(session, emptyDatabase);
            }

            // ... and delete it afterwards
            foreach (ISessionInfo session in emptyDatabase.Controller.Sessions) {
                emptyDatabase.Controller.DeleteSession(session);
            }


            Assert.IsTrue(
                emptyDatabase.Controller.Sessions.Count() == 0,
                "Not all sessions deleted.");

            // Only the files associated with grids should remain in the database.
            Assert.IsTrue(
                CountGridFiles(databaseWithFiles) == CountAllFiles(emptyDatabase.Path),
                "More than just the grids remaining in source database");
        }

        [Test]
        public void TestClearDatabase() {
            Assert.IsTrue(
                CountAllFiles(emptyDatabase.Path) == 0,
                "Database should be empty in the beginning");

            foreach (ISessionInfo session in databaseWithFiles.Controller.Sessions) {
                databaseWithFiles.Controller.CopySession(session, emptyDatabase);
            }

            Assert.IsTrue(
                CountAllFiles(emptyDatabase.Path) == 63,
                "Database should be non-empty after copying");

            emptyDatabase.Controller.ClearDatabase();
            Assert.IsTrue(
                CountAllFiles(emptyDatabase.Path) == 0,
                "Files remain in a database that should be empty");
        }

        [Test]
        public void TestSafelyDeleteGrid() {
            databaseWithFiles.Controller.CopySession(databaseWithFiles.Controller.Sessions.First(), emptyDatabase);

            IGridInfo grd = emptyDatabase.Controller.Sessions.First().GetGrids().First();
            Assert.IsFalse(
                emptyDatabase.Controller.DeleteGrid(grd), // Delete operation fails: no files will be changed
                "Deleted grid in use by at least one session.");

            emptyDatabase.Controller.DeleteSession(emptyDatabase.Controller.Sessions.First());

            // after the session is deleted, the grid should be successfully removed
            Assert.IsTrue(
                emptyDatabase.Controller.DeleteGrid(grd),
                "Removal of unused grid unsuccessful.");
            Assert.IsTrue(
                CountAllFiles(emptyDatabase.Path) == 0,
                "Files remain in a database that should be empty");
        }

        [Test]
        public void TestCleanUpDatabase() {
            ISessionInfo sess = databaseWithFiles.Controller.CopySession(databaseWithFiles.Controller.Sessions.First(), emptyDatabase);

            int fileCountBefore = CountAllFiles(emptyDatabase.Path);

            // create a bunch of files in the data subdirectory and clean up the database afterwards
            for (int i = 0; i < 10; i++) {
                string rndFile = Path.Combine(emptyDatabase.Path, "data", Guid.NewGuid().ToString() + ".1.data");

                if (!File.Exists(rndFile)) {
                    FileStream fs = File.Create(rndFile);
                    fs.Close();
                }
            }
            emptyDatabase.Controller.CleanDatabase();

            int fileCountAfter = CountAllFiles(emptyDatabase.Path);

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
            foreach (ISessionInfo session in databaseWithFiles.Controller.Sessions) {
                databaseWithFiles.Controller.CopySession(session, emptyDatabase);
            }

            Random rnd = new Random();

            // this should cause an IOException
            databaseWithFiles.Controller.CopySession(databaseWithFiles.Controller.Sessions.OrderBy(
                (item) => rnd.Next()).First<ISessionInfo>(), emptyDatabase);
        }

        [Test]
        public void TestRenameSession() {
            string newName = "Something Galerkin something incompressible something something.";
            string newDescription = @"No one will ever read this description.
            Unless this test fails and someone has to debug it.";

            foreach (ISessionInfo session in databaseWithFiles.Controller.Sessions) {
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

            foreach (IGridInfo grid in databaseWithFiles.Controller.Grids) {
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
    }
}
