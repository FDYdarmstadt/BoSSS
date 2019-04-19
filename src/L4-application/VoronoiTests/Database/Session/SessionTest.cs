using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Application.SipPoisson;
using BoSSS.Foundation.IO;

namespace VoronoiTests.Database.Session
{
    class SessionTest : DatabaseTest
    {
        static ISessionInfo session;
        static ISessionInfo[] sessions;

        public static ISessionInfo GetSession()
        {
            return session ?? (session = CreateSession(500));
        }

        public static ISessionInfo[] GetSessions()
        {
            return sessions ?? (sessions = CreateSessions());
        }

        static ISessionInfo CreateSession(int numberOfVoronoiCells)
        {
            AppControl lShape = SipHardcodedControl.TestVoronoi_LDomain(numberOfVoronoiCells, db: Database);
            IApplication poisson = new SipPoissonMain();
            RunApplication(poisson, lShape);
            ISessionInfo session = poisson.CurrentSessionInfo;
            return session;
        }

        static ISessionInfo[] CreateSessions()
        {
            int numberOfSessions = 5;
            ISessionInfo[] sessions = new ISessionInfo[numberOfSessions];
            double[] numbersOfVoronoiCells = ilPSP.Utils.GenericBlas.Linspace(100, 1000, numberOfSessions);

            for (int i = 0; i < numberOfSessions; ++i)
            {
                int numberOfVoronoiCells = (int)numbersOfVoronoiCells[i];
                sessions[i] = CreateSession(numberOfVoronoiCells);
            }
            return sessions;
        }
    }
}
