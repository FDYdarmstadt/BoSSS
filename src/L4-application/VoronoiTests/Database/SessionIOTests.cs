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
using BoSSS.Foundation;

namespace VoronoiTests.Database
{
    class SessionIOTests : DatabaseTest
    {
        static SessionMethods sessionMethods;

        static SessionMethods SessionMethods {
            get {
                return sessionMethods ?? (sessionMethods = new SessionMethods(Database));
            }
        }

        static ISessionInfo CreateSession()
        {
            AppControl lShape = SipHardcodedControl.TestVoronoi_LDomain(500,db: Database);
            IApplication poisson = new SipPoissonMain();
            RunApplication(poisson, lShape);
            return poisson.CurrentSessionInfo;
        }

        [Test]
        public void AccessTimestepFields()
        {
            ISessionInfo session = CreateSession();
            IEnumerable<DGField> fields = session.Timesteps.Last().Fields;
            DGField field = fields.ElementAt(0);
        }
    }
}
