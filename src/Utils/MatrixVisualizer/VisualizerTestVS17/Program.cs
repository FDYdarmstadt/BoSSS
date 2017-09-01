using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Solution;
using ilPSP;
using ilPSP.LinSolvers;
using MatrixVisualizerVS15;
using Microsoft.VisualStudio.DebuggerVisualizers;

namespace VisualizerTest {

    class Program {

        static void Main(string[] args) {
            bool dummy;
            ilPSP.Environment.Bootstrap(args, Application.GetBoSSSInstallDir(), out dummy);

            MultidimensionalArray myMultidimensionalArray = MultidimensionalArray.Create(3, 6);
            myMultidimensionalArray[1, 2] = 3.0;
            myMultidimensionalArray[2, 5] = 2.0;
            TestShowVisualizer(myMultidimensionalArray);

            BlockDiagonalMatrix myBlockDiagonalMatrix = new BlockDiagonalMatrix(15, 3);
            myBlockDiagonalMatrix[1, 1] = 1.0;
            TestShowVisualizer(myBlockDiagonalMatrix);

            MsrMatrix myMsrMatrix = new MsrMatrix(4, 1);
            myMsrMatrix[1, 1] = 3.0;
            TestShowVisualizer(myMsrMatrix);

            MultidimensionalArray myNodeSet = MultidimensionalArray.Create(3, 2);
            myNodeSet[0, 0] = 1.0;
            myNodeSet[0, 1] = 2.0;
            myNodeSet[1, 0] = 3.0;
            myNodeSet[1, 1] = 4.0;
            myNodeSet[2, 0] = 5.0;
            myNodeSet[2, 1] = 6.0;
            NodeSet nodeSet = new NodeSet(null, myNodeSet);
            TestShowVisualizer(nodeSet);
        }

        public static void TestShowVisualizer(object objectToVisualize) {
            VisualizerDevelopmentHost myHost = new VisualizerDevelopmentHost(
                objectToVisualize, typeof(Visualizer), typeof(ObjectSource));
            myHost.ShowVisualizer();
        }
    }
}
