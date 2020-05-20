using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP;
using System;

namespace BoSSS.Application.SipPoisson
{
    /// <summary>
    /// Interior Penalty Flux
    /// </summary>
    class ipFlux : BoSSS.Solution.NSECommon.SIPLaplace
    {

        public ipFlux(double penalty_const, MultidimensionalArray cj, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, cj, "T") //
        {
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap.bndFunction["T"];
        }

        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        protected override double g_Diri(ref CommonParamsBnd inp)
        {
            double v = m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp)
        {
            double v = m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            return v;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp)
        {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType)
            {
                case BoundaryType.Dirichlet:
                    return true;

                case BoundaryType.Neumann:
                    return false;

                default:
                    throw new NotImplementedException();
            }
        }
    }
}
