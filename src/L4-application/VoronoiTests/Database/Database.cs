using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.IO;

namespace VoronoiTests.Database
{
    class Database
    {
        readonly IDatabaseInfo databaseInfo;

        public Database(string path)
        {
            string fullPath = Path.GetFullPath(path);
            databaseInfo = CreateEmptyDatabase(fullPath);
        }

        public IDatabaseInfo DatabaseInfo { get => databaseInfo; }

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

        public void DeleteDatabase()
        {
            if (databaseInfo != null)
            {
                string path = databaseInfo.Path;
                Directory.Delete(path, true);
            }
        }
    }
}
