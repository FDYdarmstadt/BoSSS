using BoSSS.Foundation;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using MatrixVisualizerVS15;
using Microsoft.VisualStudio.DebuggerVisualizers;
using System.Diagnostics;

[assembly: DebuggerVisualizer(typeof(Visualizer), typeof(ObjectSource), Target = typeof(MsrMatrix), Description = "Matrix Visualizer")]
[assembly: DebuggerVisualizer(typeof(Visualizer), typeof(ObjectSource), Target = typeof(MultidimensionalArray), Description = "Matrix Visualizer")]
[assembly: DebuggerVisualizer(typeof(Visualizer), typeof(ObjectSource), Target = typeof(NodeSet), Description = "Matrix Visualizer")]
[assembly: DebuggerVisualizer(typeof(Visualizer), typeof(ObjectSource), Target = typeof(BlockDiagonalMatrix), Description = "Matrix Visualizer")]
[assembly: DebuggerVisualizer(typeof(Visualizer), typeof(ObjectSource), Target = typeof(BlockMsrMatrix), Description = "Matrix Visualizer")]

namespace MatrixVisualizerVS15 {

    /// <summary>
    /// Visualizer dialog for Visual Studio 15
    /// </summary>
    public class Visualizer : DialogDebuggerVisualizer {

        /// <summary>
        /// Uses <see cref="MatrixVisualizer.VisualizationProvider"/> to display the matrix
        /// </summary>
        /// <param name="windowService"></param>
        /// <param name="objectProvider"></param>
        protected override void Show(IDialogVisualizerService windowService, IVisualizerObjectProvider objectProvider) {
            using (var form = MatrixVisualizer.VisualizationProvider.GetForm(objectProvider.GetObject())) {
                windowService.ShowDialog(form);
            }
        }
    }
}
