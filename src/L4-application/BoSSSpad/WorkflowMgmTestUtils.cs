using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Helper methods for tests, but Can also be used by end-users 
    /// </summary>
    public static class WorkflowMgmTestUtils {

        /// <summary>
        /// Deletes a database <paramref name="Directory"/>
        /// 
        /// Note: the database must be located beneath the <see cref="BatchProcessorClient.AllowedDatabasesPaths"/>
        /// of the <see cref="BoSSSshell.GetDefaultQueue"/>.
        /// </summary>
        public static void DeleteDatabase(string Directory) {

            foreach(var q in BoSSSshell.ExecutionQueues) {
                foreach(var allowedPath in q.AllowedDatabasesPaths) {
                    var localBaseDir = new DirectoryInfo(allowedPath.LocalMountPath);
                    if(localBaseDir.Exists) {
                        var dbDirs = localBaseDir.GetDirectories(Directory, SearchOption.TopDirectoryOnly);
                        foreach(var db in dbDirs) {
                            Console.WriteLine("Deleting database: " + db.FullName);
                            db.Delete(true);


                        }
                    } else {
                        Console.WriteLine("Warning: missing directory: " + localBaseDir.FullName);
                    }

                }
            }
        }

        /// <summary>
        /// Deletes all deployments matching the search patter <paramref name="DirectoryWildCard"/>
        /// </summary>
        public static void DeleteDeployments(string DirectoryWildCard) {

            foreach(var q in BoSSSshell.ExecutionQueues) {

                var localBaseDir = new DirectoryInfo(q.DeploymentBaseDirectory);
                if(localBaseDir.Exists) {
                    var deplDirs = localBaseDir.GetDirectories(DirectoryWildCard, SearchOption.TopDirectoryOnly);
                    foreach(var d in deplDirs) {
                        if(BoSSS.Foundation.IO.DatabaseUtils.IsValidBoSSSDatabase(d.ToString()))
                            continue;

                        Console.WriteLine("Deleting deployment: " + d.FullName);
                        try {
                            // we can be forgiving on deletion of old deployments; 
                            // an old deployment will not harm or influence the worksheet execution
                            // (for an old database, it's a different story!
                            d.Delete(true);



                        } catch(Exception e) {
                            Console.Error.WriteLine($"{e.GetType()} during deletion of {d.FullName}: {e.Message}.");
                        }
                    }
                } else {
                    Console.WriteLine("Warning: missing directory: " + localBaseDir.FullName);
                }

            }
        }

        /// <summary>
        /// Returns all deployment directories (empty ones and BoSSS databases excluded) which match the wildcard <paramref name="DirectoryWildCard"/>
        /// </summary>
        public static DirectoryInfo[] GetAllDeployments(string DirectoryWildCard) {
            var ret = new HashSet<string>(); // hash set does not seem to make the right comparison for directory info

            foreach(var q in BoSSSshell.ExecutionQueues) {

                var localBaseDir = new DirectoryInfo(q.DeploymentBaseDirectory);
                if(localBaseDir.Exists) {
                    var deplDirs = localBaseDir.GetDirectories(DirectoryWildCard, SearchOption.TopDirectoryOnly);
                    foreach(var d in deplDirs) {
                        if(BoSSS.Foundation.IO.DatabaseUtils.IsValidBoSSSDatabase(d.ToString()))
                            continue;

                        if(d.GetFiles().Length <= 0 && d.GetDirectories().Length <= 0)
                            continue;


                        ret.Add(d.ToString());
                    }
                }

            }

            return ret.Select(str => new DirectoryInfo(str)).ToArray();
        }
    }
}
