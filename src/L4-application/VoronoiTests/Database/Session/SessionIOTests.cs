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

namespace VoronoiTests.Database.Session
{
    class SessionIOTests : SessionTest
    {
        static SessionMethods sessionMethods;

        static SessionMethods SessionMethods {
            get {
                return sessionMethods ?? (sessionMethods = new SessionMethods(Database));
            }
        }

        public override void Run()
        {
            AccessTimestepFields();
        }

        [Test]
        public void AccessTimestepFields()
        {
            ISessionInfo session = GetSession();
            IEnumerable<DGField> fields = session.Timesteps.Last().Fields;
            DGField field = fields.ElementAt(0);
        }
    }
}
