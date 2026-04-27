using BoSSS.Application.BoSSSpad;
using BoSSS.Application.LsTest;
using BoSSS.Application.TutorialTests;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using FreeXNSE;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using NUnit.Framework.Constraints;
using PublicTestRunner;
using SAIDT;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;
using XESTSF;

namespace ValidationTestRunner {

    /// <summary>
    /// Extends the public tests to some which are only available in the internal part of BoSSS
    /// </summary>
    class ValidationTests :ITestTypeProvider {

        public Type[] FullTest {
            get {
                var ret = new Type[] {
                    typeof(ValidationTestRunnerMain),
                    typeof(XESTSFMain),
                    typeof(SAIDTMain), // required to have the SAIDT binary available
                    typeof(FreeXNSE.FreeXNSE) // required to have the FreeXNSE binary available
                };
                return ret;
            }
        }

        public Type[] ReleaseOnlyTests {
            get {
                var ret = new Type[]{

                };
                return ret;
            }
        }

        /// <summary>
        /// MPI, DEBUG and RELEASE
        /// </summary>
        public (Type type, int NoOfProcs)[] MpiFullTests {
            get {
                var ret = new (Type, int)[] { };
                return ret;
            }
        }

        /// <summary>
        /// MPI, only RELEASE
        /// </summary>
        public (Type type, int NoOfProcs)[] MpiReleaseOnlyTests {
            get {
                var ret = new (Type, int)[] { };
                // e.g.
                //(typeof(XDGShock.Program), 4).AddToArray(ref ret); // 2nd entry is number of MPI cores
                //(typeof(XDGShock_MPITests.XDGShock_MPITestsMain), 4).AddToArray(ref ret);
                return ret;
            }
        }

        public DirectoryInfo GetRepositoryBaseDir() {
            DirectoryInfo repoRoot;
            try {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                repoRoot = dir.Parent.Parent.Parent.Parent.Parent.Parent.Parent;

                var publicDir = repoRoot.GetDirectories("public").SingleOrDefault();
                if(publicDir == null)
                    return null;

                var src = publicDir.GetDirectories("src").SingleOrDefault();
                var libs = publicDir.GetDirectories("libs").SingleOrDefault();
                var doc = publicDir.GetDirectories("doc").SingleOrDefault();

                if(src == null || !src.Exists)
                    return null;
                //throw new Exception();
                if(libs == null || !libs.Exists)
                    return null;
                //throw new Exception();
                if(doc == null || !doc.Exists)
                    return null;
                //throw new Exception();

            } catch(Exception) {
                return null;
                //throw new IOException("Unable to find repository root. 'runjobmanger' must be invoked from its default location within the BoSSS git repository.");
            }

            return repoRoot;
        }

        virtual public bool CopyManagedAssembliesCentrally => false;

        virtual public int RetryCount {
            get {
                if(PublicTestRunnerMain.BOSSS_RUNTESTFROMBACKUP_ENVVAR) 
                    return 2;
                else 
                    return 1;
            }
        }

        virtual public bool DeleteSuccessfulTestFiles => false;
    }


    static class ValidationTestRunnerMain {

        /// <summary>
        /// simpler version of <see cref="DeleteDatabaseAndDeploymentsWhenOld(string, string, string, TimeSpan)">
        /// </summary>
        public static void DeleteDatabaseAndDeploymentsWhenOld(string ProjectName, TimeSpan DeletionAge) {
            DeleteDatabaseAndDeploymentsWhenOld(ProjectName,
                $"{ProjectName}*",
                $"delete_{ProjectName}",
                DeletionAge);
        }

        /// <summary>
        /// Deletes a database <paramref name="Directory"/> if it older than specified by <paramref name="DeletionAge"/>
        /// <see cref="DeleteDatabaseAndDeploymentsWhenOld(string, string, string, TimeSpan)"/>
        /// </summary>
        public static void DeleteDatabaseAndDeploymentsWhenOld(string Directory) {
            var t = new TimeSpan(days: 30, hours: 0, minutes: 0, seconds: 12);
            DeleteDatabaseAndDeploymentsWhenOld(Directory, t);
        }

        /// <summary>
        /// Deletes a database, resp. moves it to a backup directory, 
        /// if the database in <paramref name="Directory"/> if it older than specified by <paramref name="DeletionAge"/>
        /// Furthermore, all matching job deploymnets will be deleted
        /// - this will cause a re-execution of all computations specified in the test.
        /// - otherwise, i.e. if the database is not deleted, existing sessions from the database may be taken if they match 
        ///    the respective compute job (<see cref="Job"/>)
        /// 
        /// Note: the database must be located beneath the <see cref="BatchProcessorClient.AllowedDatabasesPaths"/>
        /// of the <see cref="BoSSSshell.GetDefaultQueue"/>.
        /// </summary>
        public static void DeleteDatabaseAndDeploymentsWhenOld(string Directory, string DeployMents, string EnforceDeletionEnvVar, TimeSpan DeletionAge) {
            bool runfromBackup = File.Exists("BOSSS_RUNTESTFROMBACKUP.txt");
            if(runfromBackup) {
                Console.WriteLine("Run-From-Backup mode activated ('BOSSS_RUNTESTFROMBACKUP.txt' file exists); not messing with any database or deployment directories.");
                return;
            }
            bool masterEnforceDeletion = File.Exists("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt");
            if(runfromBackup) {
                Console.WriteLine("Run-From-Backup mode activated ('BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt' file exists); enforcing deletion.");
            }

            var nau = DateTime.Now;

            EnforceDeletionEnvVar = EnforceDeletionEnvVar.ToUpperInvariant();


            bool enforce = false;
            if(EnforceDeletionEnvVar != null) {
                string really = System.Environment.GetEnvironmentVariable(EnforceDeletionEnvVar) ?? "";
                string really2 = System.Environment.GetEnvironmentVariable("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER") ?? "";
                if((really.IsNonEmpty() && really.Trim() != "0")
                    || (really2.IsNonEmpty() && really2.Trim() != "0")
                    || masterEnforceDeletion) {

                    Console.WriteLine($"Enforcing deletion of old runs, since environment variables {EnforceDeletionEnvVar}/BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER are set to {really}/{really2} (not 0) or file 'BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt' exists!");
                    enforce = true;

                }
            }

            //Console.WriteLine("Note Deletion deactivated !!!!!!!!!!!!!!!!!!!!!!!!!");
            //return; /*

            //List<string> DirsToIgnore = new List<string>();
            //bool delete = false; // deletion is decided based on the age of the database folder
            (var delete, var DirsToIgnore) = BackupDatabase_IfOld(Directory, DeletionAge, nau, enforce, "todelete", true);

            if(delete) {
                Console.WriteLine("deletion of old deployments in execution queues");
                foreach(var q in BoSSSshell.ExecutionQueues) {
                    // delete old deployment
                    // =====================

                    {
                        Console.WriteLine("searching for deployments in: " + q.DeploymentBaseDirectory);

                        var deplDirs = (new DirectoryInfo(q.DeploymentBaseDirectory)).GetDirectories(DeployMents, SearchOption.TopDirectoryOnly);
                        foreach(var d in deplDirs) {
                            if(DirsToIgnore.Any(keepDir => keepDir.Equals(d.FullName)))
                                continue; // don't delete the db

                            if(delete && d.Exists) {
                                try {
                                    Console.WriteLine("Deleting deployment dir: " + d.FullName);
                                    d.Delete(true);
                                    Console.WriteLine("  done.");
                                } catch(Exception e) {
                                    Console.Error.WriteLine("   caught exception: " + e);
                                }
                            }
                        }
                    }


                    Console.Out.Flush();
                }
            } else {
                Console.WriteLine("No deletion of old deployments");
            }
            //*/
        }

        /// <summary>
        /// creates a backup of a database, at the end of a successful run; will be used for the backup-runs later on.
        /// </summary>
        public static void BackupDatabase(string db_Directory_name) {
            BackupDatabase_IfOld(db_Directory_name, TimeSpan.Zero, DateTime.Now, true, "bkup", false); // backup will be a copy
        }

        private static (bool database_is_expired, List<string> backup_directory_names) BackupDatabase_IfOld(string Directory, TimeSpan DeletionAge, DateTime nau, bool enforce, string backup_name_prefix, bool deleteOriginal) {
            List<string> DirsToIgnore = new List<string>();
            bool delete = false;

            // don't catch exception because we want to know when the backup did not work
            void CopyFilesRecursively(string sourcePath, string targetPath) {

                //Now Create all of the directories
                foreach(string dirPath in System.IO.Directory.GetDirectories(sourcePath, "*", SearchOption.AllDirectories)) {
                    System.IO.Directory.CreateDirectory(dirPath.Replace(sourcePath, targetPath));
                }

                //Copy all the files & Replaces any files with the same name
                foreach(string newPath in System.IO.Directory.GetFiles(sourcePath, "*.*", SearchOption.AllDirectories)) {
                    File.Copy(newPath, newPath.Replace(sourcePath, targetPath), true);
                }
            }

            foreach(var q in BoSSSshell.ExecutionQueues) {

                Console.WriteLine("Checking queue: " + q);

                // move old databases to backup
                // ============================

                {
                    //var dirsToDelete = new HashSet<DirectoryInfo>();
                    //bool delete = false;
                    foreach(var allowedPath in q.AllowedDatabasesPaths) {
                        var localBaseDir = new DirectoryInfo(allowedPath.LocalMountPath);

                        var dbDirs = localBaseDir.GetDirectories(Directory, SearchOption.TopDirectoryOnly);
                        foreach(var db in dbDirs) {
                            var age = nau - db.CreationTime;
                            if((age > DeletionAge || enforce) && db.Exists) {
                                delete = true;
                                string NewName = Path.Combine(Path.GetDirectoryName(db.FullName), backup_name_prefix + "-" + db.CreationTime.ToString("yyyyMMMdd_HHmmss") + "." + db.Name);
                                Console.WriteLine("Renaming existing database: " + db.FullName + " -> " + NewName);
                                
                                if(deleteOriginal)
                                    db.MoveTo(NewName); // just rename
                                else 
                                    CopyFilesRecursively(db.FullName, NewName); // copy
                                
                                DirsToIgnore.Add(NewName);
                            }
                        }
                    }
                }
            }

            return (delete, DirsToIgnore);
        }

        /// <summary>
        /// Runs some worksheet contained in the BoSSS handbook.
        /// </summary>
        static public void RunWorksheet(string NotebookPartialPath, bool allowErrors = false) {
            using(new BoSSS.Application.TutorialTests.NotebookRunner(NotebookPartialPath, "", allowErrors)) { }
        }

        static int Main(string[] args) {
            PublicTestRunner.PublicTestRunnerMain.TimeOutSec = 24 * 3600 * 10; // 10 days


            try {
                //BoSSS.Solution.Application.InitMPI(num_threads: 1);
                //WorksheetTests_Local_short.Run__MetaJobManager();
                //Assert.IsTrue(false, "remove me");

                return PublicTestRunner.PublicTestRunnerMain._Main(args, new ValidationTests());
            } catch(Exception e) {
                // note: this seemingly useless try-catch is here since our test runner server (FDYGITRUNNER)
                // seems to silently fail on all exceptions thrown after MPI init.

                int rank, size;
                if(csMPI.Raw.Initialized()) {
                    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                    csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                } else {
                    rank = 0;
                    size = 0;
                }

                Console.WriteLine("Got some exception: " + e);

                using(var stw = new StreamWriter("Exception-" + DateTime.Now.ToString("MMMdd_HHmmss") + "." + rank + "of" + size + ".txt")) {
                    stw.WriteLine("Got some exception: " + e);
                    stw.WriteLine(e.StackTrace);
                    stw.Flush();
                    stw.Close();
                    System.Environment.Exit(-667);
                }

                return -667;
            }
        }
    }
}
