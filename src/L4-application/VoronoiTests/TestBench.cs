using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;

namespace VoronoiTests
{
    [TestFixture]
    abstract class TestBench
    {
        [TestFixtureSetUp]
        public abstract void SetUp();

        public abstract void Run();

        [TestFixtureTearDown]
        public abstract void TearDown();
    }
}
