using NUnit.Framework;

namespace VoronoiTests
{
    [TestFixture]
    public abstract class TestBench
    {
        [OneTimeSetUp]
        public abstract void SetUp();

        public abstract void Run();

        [OneTimeTearDown]
        public abstract void TearDown();
    }
}
