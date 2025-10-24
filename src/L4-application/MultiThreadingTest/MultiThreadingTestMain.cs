namespace MultiThreadingTest {

    public class MultiThreadingTestMain {
        static void Main(string[] args) {
            Console.WriteLine("Given input parameters: ");
            foreach(var arg in args) {
                Console.Write(arg);
            }

            BoSSS.Solution.Application.InitMPI(args);
            Console.WriteLine("Starting multithread tests");
            TPLtest.SimpleTestTask(1000);

            TPLtest.ParTestSeq(1000, 1000);

            TPLtest.ParTestOMP(1000, 1000);
            MPI.Wrappers.csMPI.Raw.mpiFinalize();

        }
    }

}