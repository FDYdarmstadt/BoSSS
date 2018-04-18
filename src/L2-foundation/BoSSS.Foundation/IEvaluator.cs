using System.Collections.Generic;

namespace BoSSS.Foundation {

    /// <summary>
    /// Evaluation of Spatial operators
    /// </summary>
    public interface IEvaluator {
        UnsetteledCoordinateMapping CodomainMapping { get; }
        

        CoordinateMapping DomainMapping { get; }


        SpatialOperator Owner { get; }


        void Evaluate<Tout>(double alpha, double beta, Tout output, double time = double.NaN, bool MPIexchange = true, double[] outputBndEdge = null) where Tout : IList<double>;
    }
}