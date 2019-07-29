using NUnit.Framework;

namespace VoronoiTests
{
    [TestFixture]
    public abstract class TestBench
    {
        [TestFixtureSetUp]
        public abstract void SetUp();

        public abstract void Run();

        [TestFixtureTearDown]
        public abstract void TearDown();
    }
}
