using BoSSS.Foundation;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;
using System.Threading.Tasks;

namespace BoSSS.Application.CahnHilliard {
    public class __c_Source : BoSSS.Solution.LevelSetTools.PhasefieldLevelSet.phi_Source {

        public __c_Source(double _diff = 0.0) : base(_diff, "c") { }
    }

    public class __c_Diffusion : BoSSS.Solution.LevelSetTools.PhasefieldLevelSet.phi_Diffusion {

        public __c_Diffusion(int D, double penalty_const, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(D, penalty_const, __diff, __lambda, __boundaryCondMap) { }
    }

    public class __c_Flux : phi_Flux {
        public __c_Flux(int D, Func<DGField[]> VelocityGetter, BoundaryCondMap<BoundaryType> __boundaryCondMap) : base(D, VelocityGetter, __boundaryCondMap, "c") { }
    }

    /*
    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'c'
    /// </summary>
    public class c_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace, IVolumeForm, ISupportsJacobianComponent {

        public c_Diffusion(int D, double penalty_const, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, "mu") // note: in the equation for 'c', we have the Laplacian of 'phi'
        {
            m_D = D;
            m_diff = __diff;
            m_lambda = __lambda;
            m_boundaryCondMap = __boundaryCondMap;
        }

        double m_diff;
        //double min = double.MaxValue;
        //double max = 0.0;
        double m_lambda;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => new[] { "c0" }.Cat(VariableNames.VelocityVector(m_D)).Cat(VariableNames.LevelSetGradient(m_D));



        protected override double g_Diri(ref CommonParamsBnd inp) {
            double UxN = (new Vector(inp.Parameters_IN, 1, inp.D))*inp.Normal;

            double v;
            if (UxN >= 0) {
                // outflow
                v = 1.0;
            } else {
                // inflow
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
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                //case BoundaryType.
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    //return true;

                default:
                    throw new NotImplementedException();
            }
        }

        public override double Nu(double[] x, double[] p, int jCell) {
            double n = 0.0;

            for (int d = 0; d < m_D; d++) {
                // n += p[1 + m_D + d].Pow2();
                n += p[1 + m_D + d]*p[1 + m_D + d];
            }

            //double D = 0.0;
            //if (n.Sqrt() < 1.0 / (Math.Sqrt(2) * m_diff))
            //{
            //    D = 0.027 / (1 + n * m_diff.Pow2());
            //}
            //else
            //{
            // D = Math.Exp(-Math.Pow(2.0, 3.0) * n * Math.Pow(m_diff, 2.0));
            //}

            //if (D < min || D > max) Console.WriteLine("min/max: " + D);
            //max = Math.Max(D, max);
            //min = Math.Min(D, min);

            if (m_lambda == 0.0) {
                return -m_diff;
                //return -D;//-m_diff * U;//-gradUxN.L2Norm() * m_diff;
            } else if (m_lambda > 0.0 && m_lambda <= 1.0) {
                // return -m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0), 0.0);
                double ret = -m_diff * 1 - m_lambda * p[0]*p[0];
                if (ret > 0) {
                    return -m_diff * 1 - m_lambda * p[0]*p[0];
                } else {
                    return 0.0;
                }
                //return -(Math.Abs(p[0]) - Math.Min(p[0].Pow2(),1.0));//-1.0 * Math.Abs(p[0]);//-m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0), 0.0);
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
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency                        
                    }

                    Acc *= 0.0;            //switch, 0.0 seems stable, 1.0 explodes
                    Acc *= this.m_alpha;

                }
            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }

        // linear term return self
        public override IEquationComponent[] GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { this };
        }
    }
    */


    /// <summary>
    /// Transport flux for Cahn-Hilliard
    /// </summary>
    public class c_Flux : IVolumeForm, IEdgeForm, ISupportsJacobianComponent, IParameterHandling {
        public c_Flux(int D, Func<DGField[]> GetVelVector, BoundaryCondMap<BoundaryType> __boundaryCondMap) {
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap?.bndFunction["c"];
            m_GetVelVector = GetVelVector;
        }

        Func<DGField[]> m_GetVelVector;

        int m_D;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ArgumentOrdering => new string[] { "c" };

        public IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            if (Arguments.Length != 1)
                throw new ArgumentException();
            if (Parameters.Length != m_D)
                throw new ArgumentException();

            // Velocity Vector is provided externally -> no update required
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            if (Arguments.Length != 1)
                throw new ArgumentException();

            var Vel = m_GetVelVector();
            if (Vel == null || Vel.Length != m_D)
                throw new ArgumentException();

            return Vel;
        }


        public TermActivationFlags BoundaryEdgeTerms => InnerEdgeTerms | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double UxN = 0;
            for (int d = 0; d < m_D; d++) {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double c;
            if (UxN >= 0) {
                c = _uIN[0];
            } else {
                //c =m_bndFunc[inp.EdgeTag](inp.X, inp.time);
                c = -1.0;
            }

            return c * UxN * _vIN;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double UxN = 0;
            for (int d = 0; d < m_D; d++) {
                UxN += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];
            }

            double c;
            if (UxN >= 0) {
                c = _uIN[0];
            } else {
                c = _uOUT[0];
            }

            return c * UxN * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            double c = U[0];
            for (int d = 0; d < m_D; d++) {
                acc += c * cpv.Parameters[d] * GradV[d];
            }

            return -acc;
        }

        // linear term return self
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { this };
        }
    }


    /*
    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'c'
    /// </summary>
    public class phi_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace, ISupportsJacobianComponent {

        public phi_Diffusion(int D, double penalty_const, double __cahn, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, "c") // note: in the equation for 'phi', we have the Laplacian of 'c'
        {
            // m_cahn = __cahn * __cahn;
            m_cahn = __cahn;
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap?.bndFunction["c"];
        }

        double m_cahn;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);
        Func<double[], double, double>[] m_bndFunc;


        protected override double g_Diri(ref CommonParamsBnd inp) {
            double UxN = (new Vector(inp.Parameters_IN, 0, inp.D))*inp.Normal;

            double v;
            if (UxN >= 0) {
                // outflow
                v = 0.0;
            } else {
                // inflow
                v = m_bndFunc[inp.EdgeTag](inp.X, inp.time);
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
                    return false;

                //case BoundaryType.Flow:
                //    return true;

                default:
                    throw new NotImplementedException();
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
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency    
                    }

                    Acc *= 0.0;                 //switch, 0.0 seems stable, 1.0 explodes
                    Acc *= this.m_alpha;

                }
            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }


    }
    */

    /*
    /// <summary>
    /// nonlinear source term in the 'phi'-equation
    /// </summary>
    public class phi_Source : IVolumeForm, ISupportsJacobianComponent, IParameterHandling {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public phi_Source(bool __inc, double _cahn = 0.0) {
            m_inc = __inc;
            m_scale = _cahn;
        }


        //double m_lambda;
        //double m_epsilon;
        bool m_inc;
        double m_scale;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "mu", "c" };

        public IList<string> ParameterOrdering => new[] { "c0" };

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var c = Arguments[1];
            var c0 = Parameters[0];

            if (object.ReferenceEquals(c0, c))
                return;

            c0.Clear();
            c0.Acc(1.0, c);
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[1];
            return new[] { c };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double phi = U[0];
            double c = U[1];
            double c0 = cpv.Parameters[0];

            double Acc = 0;
            if (m_inc == false) {
                Acc += -phi;
            } else {
                Acc += -phi;

                // Acc += (3.0 * c0 * c0 - 1.0) * c - 2 * Math.Pow(c0, 3.0); // linearized around c0 (Taylor expansion)
                // Acc += (3.0 * c0 * c0 - 1.0) * c - 2 * c0*c0*c0; // linearized around c0 (Taylor expansion)
                Acc += c.Pow(3.0) - c; // for newton with jacobian no linearization is needed
            }


            return Acc * V;
        }

        // already linearized term return self
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            var r = new IEquationComponent[] { new jacobi_phi_Source(m_inc, m_scale) };
            return r;
        }


        /// <summary>
        /// nonlinear source term in the 'phi'-equation
        /// </summary>
        class jacobi_phi_Source : IVolumeForm, IParameterHandling {
            //public phi_Source(double __lambda, double __epsilon) {
            //    m_lambda = __lambda;
            //    m_epsilon = __epsilon;
            //}

            public jacobi_phi_Source(bool __inc, double _cahn = 0.0) {
                m_inc = __inc;
                m_scale = _cahn;
            }


            //double m_lambda;
            //double m_epsilon;
            bool m_inc;
            double m_scale;

            public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

            public IList<string> ArgumentOrdering => new[] { "mu", "c" };

            public IList<string> ParameterOrdering => new[] { "c0" };

            public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
                var c = Arguments[1];
                var c0 = Parameters[0];

                if (object.ReferenceEquals(c0, c))
                    return;

                c0.Clear();
                c0.Acc(1.0, c);
            }

            public DGField[] MyParameterAlloc(DGField[] Arguments) {
                var c = Arguments[1];
                return new[] { c };
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

                double phi = U[0];
                double c = U[1];
                double c0 = cpv.Parameters[0];

                double Acc = 0;
                if (m_inc == false) {
                    Acc += -phi;
                } else {
                    Acc += -phi;

                    //Acc += (m_lambda / m_epsilon) * (c0 * c0 - 1) * c; // linearized around c0
                    //Acc += ((2*c0 - 1) * (2*c0 - 1) - 1) * (c - 0.5); // linearized around c0, 0<c<1
                    // Acc += (c0 * c0 - 1) * c; // linearized around c0
                    // Acc += 3 * c0.Pow2() * c - c; // linearized around c0 (Taylor expansion)
                    Acc += 3 * c0*c0 * c - c; // linearized around c0 (Taylor expansion)
                                              //Acc += c.Pow(3) - c; // for newton with jacobian no linearization is needed
                }


                return Acc * V;
            }
        }

    }

    */

    /*
    /// <summary>
    /// linear "source" term of phi in c equation in Model A
    /// </summary>
    public class c_Source : IVolumeForm, ISupportsJacobianComponent, IParameterHandling {
        //public phi_Source(double __lambda, double __epsilon) {
        //    m_lambda = __lambda;
        //    m_epsilon = __epsilon;
        //}

        public c_Source(double _diff = 0.0) {
            m_diff = _diff;
        }



        //double m_lambda;
        //double m_epsilon;
        double m_diff;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "c" };
        public IList<string> ParameterOrdering => new[] { "c0" };

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            var c = Arguments[0];
            var c0 = Parameters[0];

            if (object.ReferenceEquals(c0, c))
                return;

            c0.Clear();
            c0.Acc(1.0, c);
        }

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            var c = Arguments[0];
            return new[] { c };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double c = U[0];
            double c0 = cpv.Parameters[0];
            // values seem shifted without this offset hack

            double Acc = 0;

            //Acc += (3.0 * c0 * c0 - 1.0) * c - 2 * Math.Pow(c0, 3.0); // linearized around (Taylor expansion)
            Acc += c.Pow(3.0) - c; // without linearization

            return m_diff * Acc * V;
        }

        // already linearized term return self
        IEquationComponent[] ISupportsJacobianComponent.GetJacobianComponents(int m_D) {
            return new IEquationComponent[] { new jacobi_c_Source(m_diff) };
        }


        /// <summary>
        /// linear "source" term of phi in c equation in Model A
        /// </summary>
        class jacobi_c_Source : IVolumeForm, IParameterHandling {
            //public phi_Source(double __lambda, double __epsilon) {
            //    m_lambda = __lambda;
            //    m_epsilon = __epsilon;
            //}

            public jacobi_c_Source(double _diff = 0.0) {
                m_diff = _diff;
            }



            //double m_lambda;
            //double m_epsilon;
            double m_diff;

            public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

            public IList<string> ArgumentOrdering => new[] { "c" };
            public IList<string> ParameterOrdering => new[] { "c0" };

            public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
                var c = Arguments[0];
                var c0 = Parameters[0];

                if (object.ReferenceEquals(c0, c))
                    return;

                c0.Clear();
                c0.Acc(1.0, c);
            }

            public DGField[] MyParameterAlloc(DGField[] Arguments) {
                var c = Arguments[0];
                return new[] { c };
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

                double c = U[0];
                double c0 = cpv.Parameters[0];
                // values seem shifted without this offset hack

                double Acc = 0;

                Acc += 3
                    * c0.Pow2()
                    * c
                    - c; // linearized around c0 (Taylor expansion)

                return m_diff * Acc * V;
            }
        }


    }

    */
}
