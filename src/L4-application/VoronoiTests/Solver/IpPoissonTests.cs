using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Application.SipPoisson;
using NUnit.Framework;

namespace VoronoiTests.Solver
{
    class IpPoissonTests : BoSSSTestBench
    {
        public override void Run() {
            Test_LDomain();
        }

        [Test]
        void Test_LDomain()
        {
            int numberOfVoronoiCells = 20;
            AppControl lShape = VoronoiControl.TestVoronoi_LDomain(numberOfVoronoiCells, NoOfLlyodsIter: 20);
            IApplication poisson = new SipPoissonMain();
            RunApplication(poisson, lShape);
        }
    }
}
