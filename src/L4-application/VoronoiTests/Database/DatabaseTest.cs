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
        static Database database;

        public static IDatabaseInfo Database {
            get {
                if(database == null)
                {
                    database = InitializeDatabase();
                }
                return database.DatabaseInfo;
            }
        }

        static Database InitializeDatabase()
        {
            Database emptyDatabase = new Database("..\\..\\bosss_TestDatabase_Voronoi");
            return emptyDatabase;
        }

        public override void SetUp()
        {
            base.SetUp();
            database = InitializeDatabase();
        }

        public override void TearDown()
        {
            if(database != null)
                database.DeleteDatabase();
            base.TearDown();
        } 
    }
}
