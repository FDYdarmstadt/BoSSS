using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'phi'
    /// </summary>
    public class phi_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace
    {

        public phi_Diffusion(int D, double penalty_const, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, "mu") // note: in the equation for 'phi', we have the Laplacian of 'mu'
        {
            m_D = D;
            m_diff = __diff;
            m_boundaryCondMap = __boundaryCondMap;
            m_lambda = __lambda;
        }

        double m_diff;
        double m_lambda;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => null;

        protected override double g_Diri(ref CommonParamsBnd inp) {
            double UxN = 0;
            for (int d = 0; d < m_D; d++) {
                UxN += (inp.Parameters_IN[d + 1]) * inp.Normal[d];
            }

            double v;
            if (UxN >= 0) {
                v = 1.0;
            } else {
                v = 0.0;
            }

            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType) {
                case BoundaryType.Wall:
                case BoundaryType.Slip:
                case BoundaryType.SlipSymmetry:
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                case BoundaryType.Inlet:
                case BoundaryType.Outlet:
                case BoundaryType.Outflow:
                case BoundaryType.Pressure_Dirichlet:
                default:
                    return true;
            }
        }

        public override double Nu(double[] x, double[] p, int jCell) {

            //double n = 0.0;
            //double D = 0.0;

            //for (int d = 0; d < m_D; d++)
            //{
            //    n += p[1 + m_D + d].Pow2();
            //}

            //D = 0.01 * Math.Exp(-8.0 * n * m_diff.Pow2());
            //D = 1e-7;

            if (m_lambda == 0.0) {
                return -m_diff; // -D;
            } else if (m_lambda > 0.0 && m_lambda <= 1.0) {
                return -m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0), 0.0);//-D * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0),0.0);
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                // inflow / outflow
                if (g_D == 0) {
                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                } else {
                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        //Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency 
                        //Acc += (nuA * 0.0) * (_vA) * nd;        // consistency, homogenous Neumann
                    }
                    Acc *= this.m_alpha;

                }
            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }

    /// <summary>
    /// Transport flux for Cahn-Hilliard
    /// </summary>
    public class phi_Flux : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {
        public phi_Flux(int D, BoundaryCondMap<BoundaryType> __boundaryCondMap, string LevelSetName = "phi") {
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap?.bndFunction[LevelSetName];
            m_LevelSetName = LevelSetName; // depending on the context, solvers might prefer a different variable name
        }

        protected string m_LevelSetName; // depending on the context, solvers might prefer a different variable name
        protected int m_D;
        protected BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        protected Func<double[], double, double>[] m_bndFunc;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ArgumentOrdering => new string[] { m_LevelSetName };

        public IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);

        public TermActivationFlags BoundaryEdgeTerms => InnerEdgeTerms | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            // expand for treatment of input functions, for now hardcode to -1.0
            double UxN = 0;
            for (int d = 0; d < m_D; d++) {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double phi;
            if (UxN >= 0) {
                phi = _uIN[0];
            } else {
                phi = -1.0;//m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }

            return phi * UxN * _vIN;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double UxN = 0;
            for (int d = 0; d < m_D; d++) {
                UxN += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];
            }

            double phi;
            if (UxN >= 0) {
                phi = _uIN[0];
            } else {
                phi = _uOUT[0];
            }

            return phi * UxN * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            double phi = U[0];
            for (int d = 0; d < m_D; d++) {
                acc += phi * cpv.Parameters[d] * GradV[d];
            }

            return -acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }
    }

    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'mu'
    /// </summary>
    public class mu_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace {

        public mu_Diffusion(int D, double penalty_const, double __cahn, BoundaryCondMap<BoundaryType> __boundaryCondMap, string LevelSetName = "phi")
            : base(penalty_const, LevelSetName) // note: in the equation for 'mu', we have the Laplacian of 'phi'
        {
            m_D = D;
            m_cahn = __cahn * __cahn;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap?.bndFunction[LevelSetName];
            m_LevelSetName = LevelSetName;
        }

        protected string m_LevelSetName; // depending on the context, solvers might prefer a different variable name
        double m_cahn;
        int m_D;
        public override IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        protected override double g_Diri(ref CommonParamsBnd inp) {
            double UxN = 0;
            for (int d = 0; d < m_D; d++) {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double v;
            if (UxN >= 0) {
                v = 0.0;
            } else {
                // for now set inflow value to -1.0
                v = -1.0; // m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType) {
                case BoundaryType.Wall:
                case BoundaryType.Slip:
                case BoundaryType.SlipSymmetry:
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                case BoundaryType.Inlet:
                case BoundaryType.Outlet:
                case BoundaryType.Outflow:
                case BoundaryType.Pressure_Dirichlet:
                default:
                    return true;
            }
        }

        public override double Nu(double[] x, double[] p, int jCell) {
            return -m_cahn;
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                // inflow / outflow
                if (g_D != 0) {
                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                } else {
                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        //Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency 
                        //Acc += (nuA * 0.0) * (_vA) * nd;        // consistency Neumann
                    }

                    //double D0 = 1.0;
                    //Acc -= (nuA * _vA) * D0 * (_uA[0]-inp.Parameters_IN[inp.D])/ 0.01;        // Outflow boundary condition see S. Dong

                    Acc *= this.m_alpha;

                }
            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }
    
    /// <summary>
    /// nonlinear source term in the 'phi'-equation
    /// </summary>
    public class mu_Source : IVolumeForm, ISupportsJacobianComponent {

        public mu_Source(string LevelSetName = "phi") {
            m_LevelSetName = LevelSetName; // depending on the context, solvers might prefer a different variable name
        }

        protected string m_LevelSetName; // depending on the context, solvers might prefer a different variable name

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "mu", m_LevelSetName };

        public IList<string> ParameterOrdering => null; // new[] { "phi0" };

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double mu = U[0];
            double phi = U[1];
            //double phi0 = cpv.Parameters[0];

            double Acc = 0;

            Acc += -mu;

            //Acc += (3 * phi0 * phi0 - 1.0) * phi - 2 * Math.Pow(phi0, 3.0); // linearized around phi0 (Taylor expansion)
            Acc += phi.Pow(3.0) - phi; // when using Newton with Jacobian linearization is not needed

            return Acc * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var r = new IEquationComponent[] { new jacobi_mu_Source(m_LevelSetName) };
            return r;
        }

        private class jacobi_mu_Source : IVolumeForm {

            public jacobi_mu_Source(string LevelSetName) {
                m_LevelSetName = LevelSetName; // depending on the context, solvers might prefer a different variable name
            }

            protected string m_LevelSetName; // depending on the context, solvers might prefer a different variable name

            public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

            public IList<string> ArgumentOrdering => new[] { "mu", m_LevelSetName };

            public IList<string> ParameterOrdering => new[] { m_LevelSetName + "_lin" };

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                double mu = U[0];
                double phi = U[1];
                double phi0 = cpv.Parameters[0];

                double Acc = 0;

                Acc -= mu;
                Acc += 3 * phi0.Pow2() * phi - phi; // linearized around c0 (Taylor expansion)

                return Acc * V;
            }
        }
    }

    
    
    
    /// <summary>
    /// source term of phi in Model A
    /// </summary>
    public class phi_Source : IVolumeForm, ISupportsJacobianComponent
    {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public phi_Source(double _diff = 0.0, string LevelSetName = "phi") {
            m_diff = _diff;
            m_LevelSetName = LevelSetName; // depending on the context, solvers might prefer a different variable name
        }



        protected string m_LevelSetName; // depending on the context, solvers might prefer a different variable name
        //double m_lambda;
        //double m_epsilon;
        double m_diff;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { m_LevelSetName};
        public IList<string> ParameterOrdering => null; //new[] { "phi0" };

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double phi = U[0];
            //double phi0 = cpv.Parameters[0];
            // values seem shifted without this offset hack

            double Acc = 0;

            //Acc += (3.0 * phi0 * phi0 - 1.0) * phi - 2 * Math.Pow(phi0, 3.0); // linearized around (Taylor expansion)
            Acc += phi.Pow(3.0) - phi; // when using Newton with Jacobian linearization is not needed

            return m_diff * Acc * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new jacobi_phi_Source(m_diff, m_LevelSetName) };
        }

        private class jacobi_phi_Source : IVolumeForm {
            private double m_diff;

            public jacobi_phi_Source(double m_diff, string LevelSetName) {
                this.m_diff = m_diff;
                m_LevelSetName = m_LevelSetName;
            }

            protected string m_LevelSetName; // depending on the context, solvers might prefer a different variable name

            public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

            public IList<string> ArgumentOrdering => new[] { m_LevelSetName };

            public IList<string> ParameterOrdering => new[] { "phi0" };

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                double phi = U[0];
                double phi0 = cpv.Parameters[0];

                double Acc = 0;

                Acc += 3 * phi0.Pow2() * phi - phi; // linearized around c0 (Taylor expansion)

                return Acc * V;
            }
        }
    }

    /*
    /// <summary>
    /// Correction term to counter along the interface diffusion
    /// </summary>
    class phi_CurvatureCorrection : IVolumeForm, ISupportsJacobianComponent
    {

        public phi_CurvatureCorrection(int D, double _cahn = 0.0)
        {
            m_cahn = _cahn.Pow2();
            m_D = D;
        }

        int m_D;
        double m_cahn;

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "phi", VariableNames.Curvature };
        public IList<string> ParameterOrdering => null;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double Acc = 0.0;
            double[] grad = new double[m_D];
            for (int d = 0; d < m_D; d++)
            {
                grad[d] = GradU[0, d];
            }
            Acc += grad.L2Norm();
            Acc *= m_cahn * U[1];

            // sign minus should be correct, plus produces more sensual results sign of curvature?
            return Acc * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            var VolDiff = new VolumeFormDifferentiator(this, m_D);
            return new IEquationComponent[] { VolDiff };
        }
    }

    
    class curvature_Direct : IEquationComponent, IVolumeForm, ISupportsJacobianComponent
    {

        public curvature_Direct(int _D)
        {
            m_D = _D;
        }
        int m_D;

        public IList<string> ArgumentOrdering => new[] { VariableNames.Curvature };

        public IList<string> ParameterOrdering => new[] { "D" + VariableNames.Curvature };

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double Acc = 0.0;

            Acc += 1e-8 * (cpv.Parameters[0] - U[0]) * V;

            return Acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new[] { this };
        }

    }

    
  class curvature_Source : IEquationComponent, IVolumeForm, ISupportsJacobianComponent
  {

      public curvature_Source(int _D)
      {
          m_D = _D;
      }
      int m_D;

      public IList<string> ArgumentOrdering => new[] { VariableNames.Curvature };

      public IList<string> ParameterOrdering => null;

      public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

      public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
      {
          double Acc = 0.0;

          Acc += -U[0] * V;

          return Acc;
      }

      public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
      {
          return new[] { this };
      }
  }


  class curvature_Divergence : IVolumeForm, IEdgeForm, IEquationComponent, ISupportsJacobianComponent
  {
      public curvature_Divergence(int D, double penalty, double limit)
      {
          m_D = D;
          m_limit = limit;
          m_penalty = penalty;
      }

      int m_D;
      double m_limit;
      double m_penalty;

      public CellMask m_cells;

      public IList<string> ArgumentOrdering => new[] { "phi", VariableNames.Curvature };

      public IList<string> ParameterOrdering => null;

      public TermActivationFlags VolTerms => TermActivationFlags.GradUxGradV;

      public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;

      public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;

      /// <summary>
      /// a little switch...
      /// </summary>
      protected double m_alpha = 1.0;

      public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
      {
          double acc = 0;

          double[] grad = new double[m_D];
          for (int d = 0; d < cpv.D; d++)
          {
              grad[d] = GradU[0, d];
              acc -= GradU[0, d] * GradV[d] * this.m_alpha;
          }
          double norm = Math.Max(grad.L2Norm(), m_limit);
          acc *= 1 / norm;

          return acc;
      }

      public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
      {
          double Acc = 0.0;

          double pnlty = this.GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);

          Acc += Flux(_uIN, _uOUT, _Grad_uIN, _Grad_uOUT, inp.Normal) * (_vIN - _vOUT);
          Acc *= this.m_alpha;

          return Acc;
      }

      private enum FluxType
      {
          Central,

          LaxFriedrich,

          Upwind,

          Godunov
      }

      private double Flux(double[] uIN, double[] uOUT, double[,] grad_uIN, double[,] grad_uOUT, Vector normal)
      {
          var FType = FluxType.Central;

          double Acc = 0.0;

          double[] MeanGrad = new double[m_D];
          double norm;

          switch (FType)
          {

              case FluxType.Central:

                  for (int d = 0; d < m_D; d++)
                      MeanGrad[d] = 0.5 * (grad_uIN[0, d] + grad_uOUT[0, d]);

                  norm = Math.Max(MeanGrad.L2Norm(), m_limit);

                  for (int d = 0; d < m_D; d++)
                      Acc += MeanGrad[d] / norm * normal[d];

                  break;
              case FluxType.Upwind:

                  // How to calculate upwind direction
                  double P = 0.0;

                  if (P > 0)
                  {
                      for (int d = 0; d < m_D; d++)
                          MeanGrad[d] = grad_uIN[0, d];
                  }
                  else
                  {
                      for (int d = 0; d < m_D; d++)
                          MeanGrad[d] = grad_uOUT[0, d];
                  }

                  norm = Math.Max(MeanGrad.L2Norm(), m_limit);

                  for (int d = 0; d < m_D; d++)
                      Acc += MeanGrad[d] / norm * normal[d];

                  break;
              case FluxType.LaxFriedrich:

                  double[] GradIN = new double[m_D];
                  double[] GradOUT = new double[m_D];
                  double gamma = 0.0;

                  for (int d = 0; d < m_D; d++)
                  {
                      GradIN[d] = grad_uIN[0, d];
                      GradOUT[d] = grad_uOUT[0, d];
                      MeanGrad[d] = GradIN[d] - GradOUT[d];
                  }

                  double normIN = Math.Max(GradIN.L2Norm(), m_limit);
                  double normOUT = Math.Max(GradOUT.L2Norm(), m_limit);
                  gamma = 0.1 / Math.Min(normIN, normOUT);

                  for (int d = 0; d < m_D; d++)
                      Acc += 0.5 * (GradIN[d] / normIN + GradOUT[d] / normOUT) * normal[d];

                  Acc -= gamma * MeanGrad.L2Norm();

                  break;
              case FluxType.Godunov:
                  break;
              default:
                  throw new NotSupportedException();
          }

          return Acc;
      }

      public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
      {
          double Acc = 0.0;

          double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);

          double[] grad = new double[m_D];
          for (int d = 0; d < inp.D; d++)
          {
              grad[d] = _Grad_uA[0, d];
          }
          double norm = Math.Max(grad.L2Norm(), m_limit);

          bool IsDirichlet = false;

          if (IsDirichlet)
          {
              // inhom. Dirichlet b.c.
              // +++++++++++++++++++++

              double g_D = 0.0; //this.g_Diri(ref inp);

              for (int d = 0; d < inp.D; d++)
              {
                  double nd = inp.Normal[d];
                  Acc += (_Grad_uA[0, d] / norm) * (_vA) * nd;        // consistency
                  Acc += (_Grad_vA[d] / norm) * (_uA[0] - g_D) * nd;  // symmetry
              }
              Acc *= this.m_alpha;

              Acc -= (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty

          }
          else
          {

              double g_N = grad.InnerProd(inp.Normal) / norm;

              Acc += g_N * _vA * this.m_alpha;
          }
          return Acc;
      }

      /// <summary>
      /// Length scales used in <see cref="GetPenalty"/>
      /// </summary>
      protected MultidimensionalArray LengthScales;

      /// <summary>
      /// update of penalty length scales.
      /// </summary>
      public virtual void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {

          double _D = cs.GrdDat.SpatialDimension;
          double _p = DomainDGdeg.Max();

          double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
          double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

          m_penalty_deg = Math.Max(penalty_deg_tri, penalty_deg_sqr);

          this.LengthScales = cs.CellLengthScales;
      }

      /// <summary>
      /// penalty scaling through polynomial degree
      /// </summary>
      double m_penalty_deg;


      /// <summary>
      /// computation of penalty parameter according to:
      /// An explicit expression for the penalty parameter of the
      /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
      /// look at formula (7) in cited paper
      /// </summary>
      protected virtual double GetPenalty(int jCellIn, int jCellOut) {
          double cj_in = 1.0/LengthScales[jCellIn];
          double mu = m_penalty * m_penalty_deg* cj_in;
          if(jCellOut >= 0) {
              double cj_out = 1.0/LengthScales[jCellOut];
              mu = Math.Max(mu, m_penalty* m_penalty_deg * cj_out);
          }

          if(mu.IsNaNorInf())
              throw new ArithmeticException("Inf/NaN in penalty computation.");

          return mu;
      }

      public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
      {
          var EdgeDiff = new EdgeFormDifferentiator(this, m_D);
          var VolDiff = new VolumeFormDifferentiator(this, m_D);
          return new IEquationComponent[] { EdgeDiff, VolDiff };
      }
  }
  */
}
