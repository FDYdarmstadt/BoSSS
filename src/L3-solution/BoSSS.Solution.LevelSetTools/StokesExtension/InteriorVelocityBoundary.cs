using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {
    
    /// <summary>
    /// Implements the 'internal Dirichlet boundary', resp. the 'singular source term'
    /// which couples the artificial Stokes problem (used for construction of the extension velocity)
    /// to the physical one.
    /// </summary>
    class InteriorVelocityBoundary : ILevelSetForm, IParameterHandling, ILevelSetEquationComponentCoefficient {
        int m_d;
        int m_levelSetIndex;
        DGField m_InterfaceVelocityComponent;
        string m_negativeSpecies;
        string m_positiveSpecies;
        int m_D;

        public InteriorVelocityBoundary(string positiveSpecies, string negativeSpecies, int levelSetIndex, int d, int D, DGField InterfaceVelocityComponent) {
            m_InterfaceVelocityComponent = InterfaceVelocityComponent;
            m_d = d;
            m_D = D;
            m_levelSetIndex = levelSetIndex;
            m_positiveSpecies = positiveSpecies;
            m_negativeSpecies = negativeSpecies;
        }

        public int LevelSetIndex => m_levelSetIndex;

        public string PositiveSpecies => m_positiveSpecies;

        public string NegativeSpecies => m_negativeSpecies;

        public TermActivationFlags LevelSetTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => VariableNames.VelocityVector(m_D);

        public IList<string> ParameterOrdering => VariableNames.AsLevelSetVariable("Interface", VariableNames.VelocityVector(m_D));

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            //Console.WriteLine("ivb parameter " + ParameterOrdering[0] + " " + inp.Parameters_IN[0]);
            //Console.WriteLine("ivb parameter " + ParameterOrdering[1] + " " + inp.Parameters_IN[1]);
            //Console.WriteLine("ivb parameter " + ParameterOrdering[2] + " " + inp.Parameters_IN[2]);
            //Console.WriteLine("ivb parameter out " + ParameterOrdering[0] + " " + inp.Parameters_OUT[0]);
            //Console.WriteLine("ivb parameter out " + ParameterOrdering[1] + " " + inp.Parameters_OUT[1]);
            //Console.WriteLine("ivb parameter out " + ParameterOrdering[2] + " " + inp.Parameters_OUT[2]);
            double Ret = 0;
            double pnlty = this.Penalty(inp.jCellIn, inp.jCellOut);
            Ret += (uA[m_d] - inp.Parameters_IN[m_d]) * (vA) * pnlty;
            Ret += (uB[m_d] - inp.Parameters_OUT[m_d]) * (vB) * pnlty;

            return Ret;
        }

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            //throw new NotImplementedException();
            Console.WriteLine("ivb parameter " + ParameterOrdering[0] + " " + Parameters[0].L2Norm());
            Console.WriteLine("ivb parameter " + ParameterOrdering[1] + " " + Parameters[1].L2Norm());
            Console.WriteLine("ivb parameter " + ParameterOrdering[2] + " " + Parameters[2].L2Norm());
            if(Parameters.Length != 1) {
                throw new ArgumentException();
            }
            if(!object.ReferenceEquals(Parameters[0], m_InterfaceVelocityComponent)) {
                Parameters[0].Clear();
                Parameters[0].AccLaidBack(1.0, m_InterfaceVelocityComponent);
            }

        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            return new DGField[] { m_InterfaceVelocityComponent };
        }

        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base = 4.5;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        protected double Penalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];
            double penaltySizeFactor_B = 1.0 / PosLengthScaleS[jCellOut];

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            double scaledPenalty = penaltySizeFactor * m_penalty * m_penalty_base;
            if(scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }

    }

}
