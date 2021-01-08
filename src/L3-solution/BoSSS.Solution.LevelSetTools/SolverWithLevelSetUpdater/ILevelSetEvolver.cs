using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
   
    /// <summary>
    /// Driver interface, 
    /// used by the <see cref="LevelSetUpdater"/> to evolve a particular level-set 
    /// with a particular evolution method.
    /// </summary>
    /// <remarks>
    /// Over the course of the development of the BoSSS code, many different variants to 
    /// evolve interfaces were tired out, with very different implementation nature.
    /// For example, the Extension approaches, Fourier level-sets, splines, etc.;
    /// To keep track of these highly different approaches, and e.g. allow to switch between them without problems,
    /// this universal driver interface was created.
    /// - original idea of Lauritz Beck, dec20;
    /// </remarks>
    public interface ILevelSetEvolver {


        /// <summary>
        /// Parameters upon which the level-set depends
        /// </summary>
        IList<string> ParameterNames { get; }

        /// <summary>
        /// As is seems, not used at the moment
        /// </summary>
        IList<string> VariableNames { get; }

        /// <summary>
        /// Internal by-products, which can be used e.g. for the sake of logging and plotting;
        /// cf. <see cref="LevelSetUpdater.InternalFields"/>
        /// </summary>
        IDictionary<string, DGField> InternalFields { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="levelSet"></param>
        /// <param name="time"></param>
        /// <param name="dt"></param>
        /// <param name="incremental"></param>
        /// <param name="DomainVarFields">
        /// sequence determined by <see cref="VariableNames"/> (maybe)
        /// </param>
        /// <param name="ParameterVarFields"></param>
        void MovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }
}
