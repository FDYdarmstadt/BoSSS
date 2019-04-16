using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using BoSSS.Foundation.IO;
using BoSSS.Application.BoSSSpad;

namespace VoronoiTests.Database.Session
{
    class BoSSSpadTests : SessionTest
    {
        public override void Run()
        {
            ToEstimatedGridConvergenceData();
        }

        [Test]
        public static void ToEstimatedGridConvergenceData()
        {
            ISessionInfo[] sessions = GetSessions();
            ITimestepInfo[] timesteps = sessions.Select(s => s.Timesteps.Last()).ToArray();
            string fieldName = "T";
            Plot2Ddata data = timesteps.ToEstimatedGridConvergenceData(fieldName);
        }

    }
}
