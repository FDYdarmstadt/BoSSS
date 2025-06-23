using BoSSS.Foundation.Grid;
using System.Collections.Generic;


namespace BoSSS.Foundation.Quadrature {



    /// <summary>
    /// Implementations of this interface describe which quadrature rule should
    /// be used at which quadrature item (cell or edge).
    /// </summary>
    public interface ICompositeQuadRule<out TQuadRule> : IEnumerable<IChunkRulePair<TQuadRule>>
        where TQuadRule : QuadRule {


        /// <summary>
        /// a composite quadrature rule is always associated with a grid/mesh
        /// </summary>
        IGridData GridData { get; }


        /// <summary>
        /// The number of quadrature items (cells or edges)
        /// </summary>
        int NumberOfItems {
            get;
        }


        /// <summary>
        /// The quadrature scaling/integration metric
        /// </summary>
        IIntegrationMetric IntegrationMetric {
            get;
        }
    }
}
