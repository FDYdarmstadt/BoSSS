using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP;

namespace BoSSS.Solution.LevelSetTools {
    public class CurvatureInNormalDirectionAlongLevelSet : ILevelSetForm, ISupportsJacobianComponent {

        public static double hmin = double.NaN;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public CurvatureInNormalDirectionAlongLevelSet(int iLevSet, int _d, int _D) {
            //m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;
            m_LevelSetIndex = iLevSet;
        }

        int m_D;
        int m_d;

        public virtual double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double curvature = cp.Parameters_OUT[0];
            Debug.Assert(cp.Parameters_OUT[0] == cp.Parameters_IN[0], "curvature must be continuous across interface");
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));

            double[] Normal = cp.Normal;

            double Jump = (-curvature) * Normal[m_d];

            double FlxNeg = -0.5 * Jump;
            double FlxPos = +0.5 * Jump;           

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (LevelSetIndex == 0)
                    return new string[] { VariableNames.Curvature };
                else {
                    return new string[] { VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(LevelSetIndex), VariableNames.Curvature) };
                }
            }
        }
        int m_LevelSetIndex;
        public int LevelSetIndex {
            get { return m_LevelSetIndex; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.None; }
        }
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }
    }

    /// <summary>
    /// surface tension force, in Laplace-Beltrami -- form, 
    /// boundary-line-integral term (must be used in conjunction with the surface-integral part <see cref="Curvature_LaplaceBeltrami_Surface"/>).
    /// </summary>
    public class CurvatureLaplaceBeltrami_BndLine : LinearDualValueFlux {

        public CurvatureLaplaceBeltrami_BndLine(int iLevSet, int d, bool averaging) {
            m_comp = d;
            m_averaging = averaging;
            m_iLevSet = iLevSet;
        }

        int m_comp;
        int m_iLevSet;
        bool m_averaging;

        public override IList<string> ParameterOrdering {
            get {
                if (m_iLevSet == 0)
                    return VariableNames.NormalVector(2);
                else
                    return VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.NormalVector(2));
            }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }

        static double[] tangente(double[] SurfN, double[] EdgeN) {
            Debug.Assert(SurfN.Length == EdgeN.Length);
            if (SurfN.Length != 2)
                throw new NotSupportedException();

            double[] tan = new double[] { -SurfN[1], SurfN[0] };

            if (GenericBlas.InnerProd(tan, EdgeN) < 0.0)
                tan.ScaleV(-1.0);

            return tan;
        }

        protected override void InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOuCell) {
            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormalIn = (new double[] { inp.Parameters_IN[0], inp.Parameters_IN[1] }).Normalize();
            double[] SurfaceNormalOu = (new double[] { inp.Parameters_OUT[0], inp.Parameters_OUT[1] }).Normalize();            


            double[] TangenteInn = tangente(SurfaceNormalIn, EdgeNormal);
            double[] TangenteOut = tangente(SurfaceNormalOu, EdgeNormal);

            if (m_averaging) {
                TangenteInn.ScaleV(0.5);
                TangenteInn.AccV(0.5, TangenteOut);
                TangenteInn.Normalize();
                Array.Copy(TangenteInn, TangenteOut, TangenteInn.Length);
            }

            FluxInCell = -TangenteInn[m_comp];
            FluxOuCell = +TangenteOut[m_comp];
        }



        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormalIn = (new double[] { inp.Parameters_IN[0], inp.Parameters_IN[1] }).Normalize();

            double[] TangenteInn = tangente(SurfaceNormalIn, EdgeNormal);

            FluxInCell = -TangenteInn[m_comp];
        }

        public override TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public override TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }
    }

    /// <summary>
    /// surface tension force, in Laplace-Beltrami -- form, surface-integral part (must be used in conjunction with the boundary-line-integral part <see cref="Curvature_LaplaceBeltrami_BndLine"/>).
    /// </summary>
    public class CurvatureLaplaceBeltrami_Surface : LinearFlux {


        public CurvatureLaplaceBeltrami_Surface(int iLevSet, int d) {
            m_comp = d;
            m_iLevSet = iLevSet;

            if (m_comp >= 2)
                throw new NotImplementedException("3D not implemented yet.");
        }

        //int D = 2;
        int m_comp;
        int m_iLevSet;

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return 0;
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return 0;
        }

        static MultidimensionalArray ProjMatrix(double NX, double NY) {
            var R = MultidimensionalArray.Create(2, 2);
            R[0, 0] = (1 - NX * NX);
            R[0, 1] = (-NX * NY);
            R[1, 0] = (-NY * NX);
            R[1, 1] = (1 - NY * NY);
            return R;
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            double[] SurfaceNormal = (new double[] { inp.Parameters[0], inp.Parameters[1] }).Normalize();
            //double[] SurfaceNormal = (new double[] { -inp.Xglobal[0], -inp.Xglobal[1] }).Normalize();
            double NX = SurfaceNormal[0];
            double NY = SurfaceNormal[1];
            Debug.Assert(!NX.IsNaNorInf(), "Surface normal x NAN or INf");
            Debug.Assert(!NY.IsNaNorInf(), "Surface normal y NAN or INf");

            switch (m_comp) {
                case 0:
                    output[0] = -(1 - NX * NX);
                    output[1] = -(-NX * NY);
                    return;

                case 1:
                    output[0] = -(-NY * NX);
                    output[1] = -(1 - NY * NY);
                    return;

                default:
                    throw new NotSupportedException();
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                if (m_iLevSet == 0)
                    return VariableNames.NormalVector(2);
                else
                    return VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.NormalVector(2));
            }
        }

        public override TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.None;
            }
        }

        public override TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        public override TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.None;
            }
        }
    }   
}
