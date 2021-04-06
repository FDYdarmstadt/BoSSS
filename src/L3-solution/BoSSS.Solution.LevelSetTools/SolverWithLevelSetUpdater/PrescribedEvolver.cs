using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    
    /// <summary>
    /// A <see cref="ILevelSetEvolver"/>-Driver around some prescribed function <see cref="ScalarFunctionTimeDep"/>
    /// (e.g. from the control file)
    /// for the level-set.
    /// </summary>
    /// <remarks>
    /// implemented by F Kummer, Jan. 2021.
    /// </remarks>
    public class PrescribedEvolver : ILevelSetEvolver {

        public PrescribedEvolver(ScalarFunctionTimeDep LevSetFunction) {
            m_function = LevSetFunction;
        }


        public IList<string> ParameterNames => null;

        public IList<string> VariableNames => null;

        public IDictionary<string, DGField> InternalFields => null;


        ScalarFunctionTimeDep m_function;

        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var ls = levelSet.DGLevelSet;
            ls.Clear();
            ls.ProjectField(m_function.SetTime(time + dt));
        }
    }
}
