using System.IO;
using MatrixVisualizer;
using Microsoft.VisualStudio.DebuggerVisualizers;

namespace MatrixVisualizerVS15 {

    /// <summary>
    /// Transforms the given object into something that is understood by the
    /// visualizer. The real implementation is in
    /// <see cref="VisualizationProvider"/>; this class just dispatches the
    /// calls for interoperability with different versions of Visual Studio
    /// </summary>
    public class ObjectSource : VisualizerObjectSource {

        /// <summary>
        /// Pipes data requests to the actual implementation in
        /// <see cref="VisualizationProvider.GetTarget"/>
        /// </summary>
        /// <param name="target"></param>
        /// <param name="outgoingData"></param>
        public override void GetData(object target, Stream outgoingData) {
            base.GetData(VisualizationProvider.GetTarget(target), outgoingData);
        }
    }
}
