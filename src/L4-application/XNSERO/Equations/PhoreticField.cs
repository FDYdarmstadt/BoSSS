using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.XNSERO_Solver.Equations {

    /// <summary>
    /// A chemical potential which models movement of bacteria
    /// </summary>
    /// <remarks>
    /// - added February 2021, F.Kummer
    /// - intended for cooperation project with Group of B. Liebchen in TRR 146.
    /// </remarks>
    internal class PhoreticFieldBulk : BulkEquation {
        
        /// <summary>
        /// Hardcoded for the fluid A;
        /// (we do not intend to use multiphase for the particle solver at this point)
        /// </summary>
        public override string SpeciesName => "A";

        public override double MassScale => 0.0; // no time derivative

        public override string CodomainName => "PhoreticResidual";


        public PhoreticFieldBulk() {

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Phoretic);

            AddComponent(new BulkLaplace());
        }


        class BulkLaplace : BoSSS.Solution.NSECommon.SIPLaplace {

            public BulkLaplace() : base(4.0, BoSSS.Solution.NSECommon.VariableNames.Phoretic) {

            }

            /// <summary>
            /// Diffusion coefficient
            /// </summary>
            public override double Nu(double[] x, double[] p, int jCell) {
                return 1.0;
            }

            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return true;
            }
        }
    }

    /// <summary>
    /// Boundary condition for the phoretic field at the particle.
    /// </summary>
    internal class ImmersedBoundaryPhoreticField : SurfaceEquation {
        public override string FirstSpeciesName => "A";

        public override string SecondSpeciesName => "C";

        public override string CodomainName => "PhoreticResidual";


        public ImmersedBoundaryPhoreticField(LevelSetTracker lstrk, int iLevSet) {
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Phoretic);
            AddComponent(new XLaplace_Interface(lstrk, 1.0));
            AddParameter(Solution.NSECommon.VariableNames.AsLevelSetVariable(Solution.NSECommon.VariableNames.LevelSetCGidx(iLevSet), Solution.NSECommon.VariableNames.Phoretic));
        }


        /// <summary>
        /// Neumann boundary condition for the Laplace operator at the interface
        /// </summary>
        public class XLaplace_Interface : ILevelSetForm, ILevelSetEquationComponentCoefficient {

            

            protected LevelSetTracker m_LsTrk;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="lstrk"></param>
            /// <param name="_muA">
            /// Diffusion coefficient
            /// </param>
            /// <param name="__penatly_baseFactor">
            /// multiplicative safety factor for the penalty; should be between 1 and 10.
            /// </param>
            public XLaplace_Interface(LevelSetTracker lstrk, double _muA, double __penatly_baseFactor = 4.0) {
                this.m_LsTrk = lstrk;
                this.muA = _muA;
                this.penatly_baseFactor = __penatly_baseFactor;
            }


            protected double muA;
            protected double penatly_baseFactor;


            /// <summary>
            /// Neumann boundary value
            /// </summary>
            double g_Neum(ref CommonParams inp) {
                double activity = inp.Parameters_IN[0];
                return activity;
                //return 1.0;
            }


            /// <summary>
            /// Implements a Neumann boundary condition for fluid phase A
            /// </summary>
            /// <remarks>
            /// Note: when adding additional terms, the <see cref="LevelSetTerms"/> must be updated correctly.
            /// </remarks>
            public virtual double InnerEdgeForm(ref CommonParams inp,
                double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
                double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

                //Vector N = inp.Normal;
                double Acc = 0;
                double g_N = this.g_Neum(ref inp);
                Acc += muA * g_N * vA;

                return Acc;
            }

            /// <summary>
            /// Note: for a pure Neumann boundary condition, the penalty is not required;
            /// However, if, at some point a Dirichlet or Robin b.c. should be tested, 
            /// </summary>f
            protected double GetPenalty(ref CommonParams inp) {
                
                double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

                double penaltySizeFactor = 1.0 / NegCellLengthScale;
                Debug.Assert(!double.IsNaN(penaltySizeFactor));
                Debug.Assert(!double.IsInfinity(penaltySizeFactor));
                
                double penalty_muFactor;
                penalty_muFactor = muA;

                double eta = this.penatly_baseFactor * penaltySizeFactor * penalty_muFactor * m_penalty_deg;
                if(eta.IsNaNorInf())
                    throw new ArithmeticException("Inf/NaN in penalty computation.");
                return eta;
            }
                       


            MultidimensionalArray NegLengthScaleS;
            double m_penalty_deg;

            /// <summary>
            /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate(CoefficientSet, CoefficientSet, int[], int)"/>
            /// </summary>
            public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
                NegLengthScaleS = csA.CellLengthScales;

                double _p = DomainDGdeg.Max();
                double _D = csA.GrdDat.SpatialDimension;
                double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
                double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

                m_penalty_deg = Math.Max(penalty_deg_tri, penalty_deg_sqr);
            }

            /// <summary>
            /// Note: the immersed boundary is always level set 1 (indexing of level-sets starts at 0)
            /// </summary>
            public int LevelSetIndex {
                get { return 1; }
            }

            public IList<string> ArgumentOrdering {
                get { return new string[] { BoSSS.Solution.NSECommon.VariableNames.Phoretic }; }
            }

            /// <summary>
            /// The Solid Domain
            /// </summary>
            public SpeciesId PositiveSpecies {
                get { return m_LsTrk.GetSpeciesId("C"); }
            }

            /// <summary>
            /// The fluid Domain
            /// </summary>
            public SpeciesId NegativeSpecies {
                get { return m_LsTrk.GetSpeciesId("A"); }
            }

            /// <summary>
            /// Remark: in case of doubt, use <see cref="TermActivationFlags.AllOn"/>
            /// </summary>
            public TermActivationFlags LevelSetTerms {
                get {
                    //return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | 
                    return TermActivationFlags.V; // sufficient for Neumann b.c.
                }
            }

            public IList<string> ParameterOrdering {
                get { return new string[] { Solution.NSECommon.VariableNames.AsLevelSetVariable(Solution.NSECommon.VariableNames.LevelSetCGidx(1), Solution.NSECommon.VariableNames.Phoretic) }; }
            }


        }


    }
}
