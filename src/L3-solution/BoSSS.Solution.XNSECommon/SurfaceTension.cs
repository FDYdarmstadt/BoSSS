/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

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

namespace BoSSS.Solution.XNSECommon.Operator.SurfaceTension {


    public class CurvatureBasedSurfaceTension : ILevelSetForm, ISupportsJacobianComponent {

        public static double hmin = double.NaN;

        /*
        /// <summary>
        /// compactification of curvature for better numerical representation
        /// </summary>
        /// <param name="kappa">curvature, in the range $(-infty,+infty)$</param>
        /// <returns>compact curvature in the range $(-3,3)$</returns>
        static public double Comp(double kappa) {
            return 3.0 * kappa / Math.Sqrt(kappa * kappa + 8);
        }

        static double EPS = BLAS.MachineEps;

        /// <summary>
        /// inverse of <see cref="Comp"/>
        /// </summary>
        static public double CompInv(double y) {
            
            
            double abst = 1.0e-5;
            double locomp = (-3.0 + abst);
            double hicomp = (+3.0 - abst);

            if(y < locomp) {
                double yy = y - locomp;
                Debug.Assert(yy <= 0);
                yy *= -1;

                double delta = abst - EPS;
                double yyc = delta * (1 - Math.Exp(-yy / delta));

                y = locomp - yyc;
                if(y > locomp)
                    throw new Exception();
            }

            if(y > hicomp) {
                double yy = y - hicomp;
                Debug.Assert(yy >= 0);

                double delta = abst - EPS;
                double yyc = delta * (1 - Math.Exp(-yy / delta));

                y = hicomp + yyc;
                if(y < hicomp)
                    throw new Exception();
            }

            if(y <= -3.0 || y >= 3.0)
                throw new ArgumentException();


            return y * Math.Sqrt(-8.0 / (y * y - 9));
             

            //return y;
        }
        */

        //LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public CurvatureBasedSurfaceTension(int _d, int _D, double _sigma) {
            //m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;
            this.sigma = _sigma;
       }

        int m_D;
        int m_d;
        double sigma;

        //static bool rem = true;

        /*
        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] GradU_Neg, double[,] GradU_Pos) {
            double curvature = cp.ParamsPos[0];
            Debug.Assert(cp.ParamsPos[0] == cp.ParamsNeg[0], "curvature must be continuous across interface");
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));

            //double r = 0.8;
            //double r = ((cp.x[0]).Pow2() + (cp.x[1]).Pow2()).Sqrt();
            //curvature = -1.0/r;
            //if(rem) {
            //    Console.WriteLine("curvature hardcoded.");
            //    rem = false;
            //}


            double[] Normal = cp.n;

            double presJump = (-curvature * sigma) * Normal[m_d];
            
            FlxNeg = -0.5 * presJump;
            FlxPos = +0.5 * presJump;

            //{
            //    const double CC_A = 1.0;
            //    const double CC_B = 0.0;
            //    const double MU_A = 1;
            //    const double MU_B = 1;

            //    Func<double[],double> SX =  X => -((CC_A - CC_B) * (1 + X[0].Pow2()) + 2.0 * (MU_B - MU_A));
            //    Func<double[],double> SY =  X => -((CC_A - CC_B) * (1 + X[0].Pow2()) - 2.0 * (MU_B - MU_A));
            //    double surfForce;
            //    switch(this.m_d) {
            //        case 0: surfForce = SX(cp.x); break;
            //        case 1: surfForce = SY(cp.x); break;
            //        default: throw new ApplicationException();
            //    }


            //    FlxNeg -= +0.5 * surfForce * Normal[m_d];
            //    FlxPos += +0.5 * surfForce * Normal[m_d];
            //}
            


            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));
        }
        */
        

        public double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double curvature = cp.Parameters_OUT[0];
            Debug.Assert(cp.Parameters_OUT[0] == cp.Parameters_IN[0], "curvature must be continuous across interface");
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));

            //double r = 0.8;
            //double r = ((cp.x[0]).Pow2() + (cp.x[1]).Pow2()).Sqrt();
            //curvature = -1.0/r;
            //if(rem) {
            //    Console.WriteLine("curvature hardcoded.");
            //    rem = false;
            //}


            double[] Normal = cp.Normal;

            double presJump = (-curvature * sigma) * Normal[m_d];

            double FlxNeg = -0.5 * presJump;
            double FlxPos = +0.5 * presJump;

            //{
            //    const double CC_A = 1.0;
            //    const double CC_B = 0.0;
            //    const double MU_A = 1;
            //    const double MU_B = 1;

            //    Func<double[],double> SX =  X => -((CC_A - CC_B) * (1 + X[0].Pow2()) + 2.0 * (MU_B - MU_A));
            //    Func<double[],double> SY =  X => -((CC_A - CC_B) * (1 + X[0].Pow2()) - 2.0 * (MU_B - MU_A));
            //    double surfForce;
            //    switch(this.m_d) {
            //        case 0: surfForce = SX(cp.x); break;
            //        case 1: surfForce = SY(cp.x); break;
            //        default: throw new ApplicationException();
            //    }


            //    FlxNeg -= +0.5 * surfForce * Normal[m_d];
            //    FlxPos += +0.5 * surfForce * Normal[m_d];
            //}
            
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
                return new string[] { VariableNames.Curvature };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
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
    /// Implementation of <see cref="CurvatureBasedSurfaceTension"/> for use in <see cref="XSpatialOperatorMk2.SurfaceElementOperator_Ls0"/>
    /// </summary>
    public class CurvatureBasedSurfaceTension_SurfaceOperator : IVolumeForm, ISupportsJacobianComponent {

        public static double hmin = double.NaN;        

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public CurvatureBasedSurfaceTension_SurfaceOperator(int _d, int _D, double _sigma) {
            //m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;
            this.sigma = _sigma;
        }

        int m_D;
        int m_d;
        double sigma;        

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.Curvature }.Cat(VariableNames.NormalVector(m_D));
            }
        }

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double curvature = cpv.Parameters[0];
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));

            Vector Normal = Vector.Empty(m_D);
            for (int d = 0; d < m_D; d++) {
                Normal[d] = cpv.Parameters[1+d];
            }
            Normal.NormalizeInPlace();

            double presJump = (-curvature * sigma) * Normal[m_d];

            double Flx = -0.5 * presJump;

            return Flx * V;
        }
    }


    /// <summary>
    /// Represents the artificial surface force (usually only used in manufactured solutions).
    /// </summary>
    public class SurfaceTension_ArfForceSrc  : ILevelSetForm {

        public static double hmin = double.NaN;

        //LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        public SurfaceTension_ArfForceSrc(int _d, int _D) {
            //m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;
       }

        int m_D;
        int m_d;

        public double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            //throw new NotImplementedException();
            double curvature = cp.Parameters_OUT[0];
            Debug.Assert(cp.Parameters_OUT[0] == cp.Parameters_IN[0], "curvature must be continuous across interface");
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));

            double[] Normal = cp.Normal;

            double surfForce = cp.Parameters_IN[0];
            Debug.Assert(cp.Parameters_IN[0] == cp.Parameters_OUT[0]);

            double FlxNeg = -0.5 * surfForce * Normal[m_d];
            double FlxPos = +0.5 * surfForce * Normal[m_d];

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        /*
        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParams cp,
            double[] U_Neg, double[] U_Pos, double[,] GradU_Neg, double[,] GradU_Pos) {
            double curvature = cp.ParamsPos[0];
            Debug.Assert(cp.ParamsPos[0] == cp.ParamsNeg[0], "curvature must be continuous across interface");
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));
                        
            double[] Normal = cp.n;
            
            double surfForce = cp.ParamsNeg[0];
            Debug.Assert(cp.ParamsNeg[0] == cp.ParamsPos[0]);

            FlxNeg = -0.5 * surfForce * Normal[m_d];
            FlxPos = +0.5 * surfForce * Normal[m_d];
            
            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));
        }

        public override void PrimalVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParams cp,
            double[] U_Neg, double[] U_Pos) {
            FlxNeg = 0;
            FlxPos = 0;
        }

        public override void FluxPotential(out double G, double[] U) {
            G = 0;
        }

        public override void Nu(out double NuNeg, out double NuPos, ref CommonParams cp) {
            NuPos = 1.0;
            NuNeg = 1.0;
        }
        */
        

        public IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.SurfaceForceComponent(m_d) };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
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
    }
    

    /// <summary>
    /// surface tension force, in Laplace-Beltrami -- form, surface-integral part (must be used in conjunction with the boundary-line-integral part <see cref="SurfaceTension_LaplaceBeltrami_BndLine"/>).
    /// </summary>
    public class SurfaceTension_LaplaceBeltrami_Surface : LinearFlux {


        public SurfaceTension_LaplaceBeltrami_Surface(int d, double factor) {
            m_comp = d;
            m_factor = factor;

            if (m_comp >= 2)
                throw new NotImplementedException("3D not implemented yet.");
        }

        //int D = 2;
        int m_comp;
        double m_factor;

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return 0;
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return 0;
        }

        static MultidimensionalArray ProjMatrix(double NX, double NY) {
            var R = MultidimensionalArray.Create(2, 2);
            R[0, 0] = (1 - NX*NX);
            R[0, 1] = (-NX*NY);
            R[1, 0] = (-NY*NX);
            R[1, 1] = (1 - NY*NY);
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
                    output[0] = -(1 - NX * NX) * m_factor;
                    output[1] = -(-NX * NY) * m_factor;
                    return;

                case 1:
                    output[0] = -(-NY * NX) * m_factor;
                    output[1] = -(1 - NY * NY) * m_factor;
                    return;

                default:
                    throw new NotSupportedException();
            }

            //MultidimensionalArray R = ProjMatrix(NX, NY);

            //switch (m_comp) {
            //    case 0:
            //        output[0] = -(R[0, 0].Pow2() + R[1, 0].Pow2()) * m_factor;
            //        output[1] = -(R[0, 0] * R[0, 1] + R[1, 0] * R[1, 1]) * m_factor;
            //        return;

            //    case 1:
            //        output[0] = -(R[0, 1] * R[0, 0] + R[1, 1] * R[1, 0]) * m_factor;
            //        output[1] = -(R[0, 1].Pow2() + R[1, 1].Pow2()) * m_factor;
            //        return;

            //    default:
            //        throw new NotSupportedException();
            //}

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
                return VariableNames.NormalVector(2);
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

    /// <summary>
    /// surface tension force, in Laplace-Beltrami -- form, 
    /// boundary-line-integral term (must be used in conjunction with the surface-integral part <see cref="SurfaceTension_LaplaceBeltrami_Surface"/>).
    /// </summary>
    public class SurfaceTension_LaplaceBeltrami_BndLine : LinearDualValueFlux {

        public SurfaceTension_LaplaceBeltrami_BndLine(int d, double factor, bool averaging) {
            m_comp = d;
            m_factor = factor;
            m_averaging = averaging;
        }

        int m_comp;
        double m_factor;
        bool m_averaging;

        public override IList<string> ParameterOrdering {
            get {
                return VariableNames.NormalVector(2);
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

        //static bool rem = true;


        protected override void InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOuCell) {
            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormalIn = (new double[] { inp.Parameters_IN[0], inp.Parameters_IN[1] }).Normalize();
            double[] SurfaceNormalOu = (new double[] { inp.Parameters_OUT[0], inp.Parameters_OUT[1] }).Normalize();

            /*
            double[] SurfaceNormalIn;
            double[] SurfaceNormalOu;
            {
                double XD = -inp.Normale[1], YD = inp.Normale[0];
                double X0 = inp.X[0], Y0 = inp.X[1];
                double R = 0.8;
                double alpha_1 = (-X0 * XD - Y0 * YD + Math.Sqrt(R.Pow2() * XD.Pow2() + R.Pow2() * YD.Pow2() - X0.Pow2() * YD.Pow2() + 2 * X0 * XD * Y0 * YD - XD.Pow2() * Y0.Pow2())) / (XD.Pow2() + YD.Pow2());
                double alpha_2 = -(X0 * XD + Y0 * YD + Math.Sqrt(R.Pow2() * XD.Pow2() + R.Pow2() * YD.Pow2() - X0.Pow2() * YD.Pow2() + 2 * X0 * XD * Y0 * YD - XD.Pow2() * Y0.Pow2())) / (XD.Pow2() + YD.Pow2());
                Debug.Assert(!(double.IsNaN(alpha_1) || double.IsInfinity(alpha_1)));
                Debug.Assert(!(double.IsNaN(alpha_2) || double.IsInfinity(alpha_2)));

                double alpha;
                if(alpha_1.Abs() < alpha_2.Abs())
                    alpha = alpha_1;
                else
                    alpha = alpha_2;
                Debug.Assert(alpha.Abs() < inp.GridDat.Cells.h_minGlobal * 0.1);


                double X1 = X0 + XD * alpha;
                double Y1 = Y0 + YD * alpha;

                SurfaceNormalIn = (new double[] { -X1, -Y1 }).Normalize();
                SurfaceNormalOu = (new double[] { -X1, -Y1 }).Normalize();
            }
            if(rem) {
                Console.WriteLine("fake laplace beltrami normals");
                rem = false;
            }
            */
            
            //double[] SurfaceNormalIn = (new double[] { -inp.X[0], -inp.X[1] }).Normalize();
            //double[] SurfaceNormalOu = (new double[] { -inp.X[0], -inp.X[1] }).Normalize();


            double[] TangenteInn = tangente(SurfaceNormalIn, EdgeNormal);
            double[] TangenteOut = tangente(SurfaceNormalOu, EdgeNormal);

            if(m_averaging) {
                TangenteInn.ScaleV(0.5);
                TangenteInn.AccV(0.5, TangenteOut);
                TangenteInn.Normalize();
                Array.Copy(TangenteInn, TangenteOut, TangenteInn.Length);

                //var Buf = TangenteInn;
                //TangenteInn = TangenteOut;
                //TangenteOut = Buf;
            }

            FluxInCell = -TangenteInn[m_comp] * m_factor;
            FluxOuCell = +TangenteOut[m_comp] * m_factor;
        }

       

        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormalIn = (new double[] { inp.Parameters_IN[0], inp.Parameters_IN[1] }).Normalize();

            double[] TangenteInn = tangente(SurfaceNormalIn, EdgeNormal);

            FluxInCell = -TangenteInn[m_comp]*m_factor;
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
    /// IsotropicSurfaceTension_LaplaceBeltrami with max sigma as parameter
    /// </summary>
    public class IsotropicSurfaceTension_LaplaceBeltrami_Parameter : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {
        int m_comp;

        int m_D;

        /// <summary>
        /// surface tension coefficient
        /// </summary>

        /// <summary>
        /// static contact angle (for navier-slip B.C.)
        /// </summary>
        double m_theta;

        /// <summary>
        /// friction coefficient at contact line (for navier-slip B.C.)
        /// </summary>
        double m_beta;

        /// <summary>
        /// true for quasi-static computations with moving slip-wall and non-moving interface  
        /// </summary>
        //bool m_staticInt;

        IncompressibleBcType[] m_edgeTag2Type;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] velFunction;

        /// <summary>
        /// IsotropicSurfaceTension_LaplaceBeltrami with max sigma as parameter
        /// </summary>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="edgeTag2Type"></param>
        /// <param name="bcmap"></param>
        /// <param name="theta_e"></param>
        /// <param name="beta_L"></param>
        public IsotropicSurfaceTension_LaplaceBeltrami_Parameter(int d, int D, IncompressibleBcType[] edgeTag2Type, IncompressibleBoundaryCondMap bcmap,
            double theta_e, double beta_L) {
            m_comp = d;
            m_D = D;
            //m_sigma = sigma;
            m_theta = theta_e;
            m_beta = beta_L;
            m_edgeTag2Type = edgeTag2Type;
            velFunction = bcmap.bndFunction[VariableNames.Velocity_d(d)];
            //m_staticInt = _staticInt;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            // only parameter dependent, leave this empty
            return new IEquationComponent[] { };
        }

        public virtual IList<string> ParameterOrdering {
            get {
                string[] parameters = VariableNames.NormalVector(m_D).Cat(VariableNames.MaxSigma);
                return parameters;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV;
            }
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double m_sigma = cpv.Parameters[cpv.D];
            double acc = 0;

            Vector Nsurf = SurfaceNormal(cpv.Parameters, cpv.D);
            double[,] Psurf = SurfaceProjection(Nsurf);

            for(int d = 0; d < cpv.D; d++)
                acc += -m_sigma * Psurf[m_comp, d] * GradV[d];
            return -acc;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            Vector EdgeNormal = inp.Normal;
            Vector SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN, inp.D);
            Vector SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT, inp.D);

            Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            Vector Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);

            //double m_sigmaMax = Math.Max(m_sigma[inp.jCellIn], m_sigma[inp.jCellOut]);
            double m_sigmaIn = inp.Parameters_IN[inp.D];
            double m_sigmaOt = inp.Parameters_OUT[inp.D];
            double acc = 0.5 * (m_sigmaIn * Tangente_IN[m_comp] + m_sigmaOt * Tangente_OUT[m_comp]) * (_vA - _vB);

            return -acc;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            IncompressibleBcType edgType = m_edgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Pressure_Outlet: {

                    Vector EdgeNormal = inp.Normal;
                    Vector SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN, inp.D);
                    Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
                    double m_sigma = inp.Parameters_IN[inp.D];
                    Flx_InCell = -m_sigma * Tangente_IN[m_comp];

                    break;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NavierSlip_Linear: {

                    Vector EdgeNormal = inp.Normal;
                    Vector SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN, inp.D);
                    Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

                    int D = inp.D;


                    // isotropic surface tension terms
                    for (int d = 0; d < D; d++) {
                        double m_sigma = inp.Parameters_IN[inp.D];
                        Flx_InCell -= m_sigma * (EdgeNormal[d] * Tangente_IN[d]) * EdgeNormal[m_comp];
                    }


                    if (edgType == IncompressibleBcType.NavierSlip_Linear) {
                        double[] PSnI = new double[D]; // projection of surface/level-set normal onto domain boundary tangent
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                double nn = EdgeNormal[d1] * EdgeNormal[d2];
                                if (d1 == d2) {
                                    PSnI[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                                } else {
                                    PSnI[d1] += -nn * SurfaceNormal_IN[d2];
                                }
                            }
                        }
                        double[] PSnINormal_IN = PSnI.Normalize(); // line normal: tangential to domain boundary & normal on contact line


                        // Young's relation (static contact angle)
                        double m_sigma = inp.Parameters_IN[inp.D];
                        Flx_InCell -= m_sigma * Math.Cos(m_theta) * PSnINormal_IN[m_comp];

                        // dissipative contact line force
                        // beta*(u*nL)

                        double g_D = this.velFunction[inp.EdgeTag](inp.X, inp.time);

                        for (int d = 0; d < D; d++) {
                            Flx_InCell += m_beta * ((_uA[d] - g_D) * PSnINormal_IN[d]) * PSnINormal_IN[m_comp];
                        }
                    }
                    break;
                }
                default:
                    break;
            }

            return Flx_InCell * _vA;
        }

        protected static Vector SurfaceNormal(double[] param, int D) {

            Vector N = Vector.Empty(D);

            for(int d = 0; d < D; d++) {
                N[d] = param[d];
            }

            N.NormalizeInPlace();

            return N;
        }

        protected static double[,] SurfaceProjection(Vector Nsurf) {

            int D = Nsurf.Dim;
            double[,] P = new double[D, D];

            for(int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    if(dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static Vector Tangent(Vector Nsurf, Vector Nedge) {
            Debug.Assert(Nsurf.Dim == Nedge.Dim);

            int D = Nsurf.Dim;

            Vector tau = new Vector(D);
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if(d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            tau.NormalizeInPlace();
            return tau;
        }
    }

    public abstract class IBM_ContactLine {
        protected Vector FluidSurfaceNormal(ref CommonParamsVol cpv) {
            return SurfaceNormal(cpv.Parameters, cpv.D, 0);
        }

        protected Vector SolidSurfaceNormal(ref CommonParamsVol cpv) {
            return SurfaceNormal(cpv.Parameters, cpv.D, cpv.D);
        }

        protected static double Sigma(ref CommonParamsVol cpv) {
            return cpv.Parameters[cpv.Parameters.Length - 1];
        }

        static Vector SurfaceNormal(double[] parameters, int D, int offSet) {

            Vector N = new Vector(D);

            for (int d = 0; d < D; d++) {
                N[d] = parameters[d + offSet];
            }
            N = N.Normalize();
            return N;
        }

        protected static Vector Tangent(Vector Nsurf, Vector Nedge) {
            Debug.Assert(Nsurf.Dim == Nedge.Dim);

            int D = Nsurf.Dim;

            Vector tau = new Vector(D);
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            tau = tau.Normalize();
            return tau;
        }

        protected static Vector ContactLineNormal(Vector Nsurf, Vector Nedge) {
            Debug.Assert(Nsurf.Dim == Nedge.Dim);

            int D = Nsurf.Dim;

            Vector tau = new Vector(D);
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nedge[d1] * Nedge[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nsurf[d2];
                    } else {
                        tau[d1] += -nn * Nsurf[d2];
                    }
                }
            }

            tau = tau.Normalize();
            return tau;
        }
    }

    /// <summary>
    /// LaplaceBeltrami with max sigma as parameter, for contact line between two level sets
    /// </summary>
    public class SurfaceTension_LaplaceBeltrami_Contactline : IBM_ContactLine, IVolumeForm, ISupportsJacobianComponent  {
        int d;
        int D;
        int iLevSet;

        public SurfaceTension_LaplaceBeltrami_Contactline(int d, int D, int iLevSet) {
            this.d = d;
            this.D = D;
            this.iLevSet = iLevSet;
        }

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public virtual IList<string> ParameterOrdering {
            get {
                string[] parameters = VariableNames.NormalVector(D);
                parameters = parameters.Cat(VariableNames.AsLevelSetVariable(NSECommon.VariableNames.LevelSetCGidx(iLevSet), VariableNames.NormalVector(D)));
                parameters = parameters.Cat(VariableNames.MaxSigma);
                return parameters;
            }
        }

        public IList<string> ArgumentOrdering => new string[0];

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);

            Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

            double Flx_InCell = 0;
            double m_sigma = Sigma(ref cpv);

            // isotropic surface tension terms
            Flx_InCell -= m_sigma * Tangente_IN[d];
            return Flx_InCell * V;
        }        

        // only parameter dependent, leave empty
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { };
        }
    }

    /// <summary>
    /// LaplaceBeltrami with max sigma as parameter, for contact line **between two level sets**,
    /// i.e. contact line at the immersed boundary.
    /// </summary>
    public class SurfaceTension_GNBC_Contactline : IBM_ContactLine, IVolumeForm, ISupportsJacobianComponent {
        int comp;
        int D;
        int iLevSet;
        double sigma;
        double theta_e;

        public SurfaceTension_GNBC_Contactline(int d, int D, double theta_e, double sigma, int iLevSet) {
            this.comp = d;
            this.D = D;
            this.iLevSet = iLevSet;
            this.theta_e = theta_e;
            this.sigma = sigma;
        }

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public virtual IList<string> ParameterOrdering {
            get {
                string[] parameters = VariableNames.NormalVector(D);
                parameters = parameters.Cat(VariableNames.AsLevelSetVariable(NSECommon.VariableNames.LevelSetCGidx(iLevSet), VariableNames.NormalVector(D)));
                return parameters;
            }
        }

        public IList<string> ArgumentOrdering => new string[0];

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector EdgeNormal = SolidSurfaceNormal(ref cpv);
            Vector SurfaceNormal_IN = FluidSurfaceNormal(ref cpv);

            Vector Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            Vector ContactLineNormal_IN = ContactLineNormal(SurfaceNormal_IN, EdgeNormal); // normal to contact line, tangential to solid

            double Flx_InCell = 0;
            double CosTheta = 0.0;
            for(int d = 0; d < D; d++){
                CosTheta += Tangente_IN[d] * ContactLineNormal_IN[d];
            }

            // isotropic surface tension terms
            Flx_InCell += 0.5 * sigma * (CosTheta - Math.Cos(theta_e)) * ContactLineNormal_IN[comp];
            return Flx_InCell * V;
        }

        // only parameter dependent, leave empty
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { };
        }
    }

    public class IsotropicSurfaceTension_LaplaceBeltrami : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {

        int m_comp;

        int m_D;

        /// <summary>
        /// surface tension coefficient
        /// </summary>
        double m_sigma;

        /// <summary>
        /// static contact angle (for navier-slip B.C.)
        /// </summary>
        double m_theta;

        /// <summary>
        /// friction coefficient at contact line (for navier-slip B.C.)
        /// </summary>
        double m_beta;

        /// <summary>
        /// true for quasi-static computations with moving slip-wall and non-moving interface  
        /// </summary>
        //bool m_staticInt;

        IncompressibleBcType[] m_edgeTag2Type;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] velFunction;


        public IsotropicSurfaceTension_LaplaceBeltrami(int d, int D, double sigma, IncompressibleBcType[] edgeTag2Type, IncompressibleBoundaryCondMap bcmap, 
            double theta_e, double beta_L) {
            m_comp = d;
            m_D = D;
            m_sigma = sigma;
            m_theta = theta_e;
            m_beta = beta_L;
            m_edgeTag2Type = edgeTag2Type;
            velFunction = bcmap.bndFunction[VariableNames.Velocity_d(d)];
            //m_staticInt = _staticInt;
        }

        //public virtual void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {

        //    m_sigma = (MultidimensionalArray)cs.UserDefinedValues["sigmaMaxValue"];
        //}


        public virtual IList<string> ParameterOrdering {
            get {
                return VariableNames.NormalVector(m_D);
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV; 
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double acc = 0;

            double[] Nsurf = SurfaceNormal(cpv.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            //double[,] Psurf2 = new double[cpv.D, cpv.D];
            //for (int d1 = 0; d1 < cpv.D; d1++) {
            //    for (int d2 = 0; d2 < cpv.D; d2++) {
            //        for (int dd = 0; dd < cpv.D; dd++) {
            //            Psurf2[d1, d2] += Psurf[d1, dd] * Psurf[dd, d2];
            //        }
            //    }
            //}

            for (int d = 0; d < cpv.D; d++)
                acc += - m_sigma * Psurf[m_comp, d] * GradV[d];

            // stabilization
            //for(int d = 0; d < cpv.D; d++) {
            //    for(int dd = 0; dd < cpv.D; dd++) {
            //        acc += -0.1 * GradU[m_comp, d] * Nsurf[d] * GradV[dd] * Nsurf[dd];
            //    }
            //}

            return -acc;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            //Console.WriteLine("Tangente_IN = ({0}, {1})", Tangente_IN[0], Tangente_IN[1]);
            double[] Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);
            //Console.WriteLine("Tangente_OUT = ({0}, {1})", Tangente_OUT[0], Tangente_OUT[1]);

            //double m_sigmaMax = Math.Max(m_sigma[inp.jCellIn], m_sigma[inp.jCellOut]);
            double acc = 0.5 * m_sigma* (Tangente_IN[m_comp] + Tangente_OUT[m_comp]) * (_vA - _vB);

            //double surfTdiff = (Tangente_IN[m_comp] - Tangente_OUT[m_comp]);
            //if (inp.iEdge == 14 && m_comp == 0) {
            //    if (_vA != 0.0) {
            //        //Console.WriteLine("surface tension at edge {0}", inp.iEdge);
            //        //Console.WriteLine("jCellIn {0}", inp.jCellIn);
            //        //Console.WriteLine("component {0}", m_comp);
            //        Console.WriteLine("test function _vA");
            //        Console.WriteLine("acc = {0}", acc);
            //    } else if (_vB != 0.0) {
            //        //Console.WriteLine("surface tension at edge {0}", inp.iEdge);
            //        //Console.WriteLine("jCellIn {0}", inp.jCellIn);
            //        //Console.WriteLine("component {0}", m_comp);
            //        Console.WriteLine("test function _vB");
            //        Console.WriteLine("acc = {0}", acc);
            //    } else {
            //        Console.WriteLine("should be zero: acc = {0}", acc);
            //    }
            //}

            return -acc;
        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            IncompressibleBcType edgType = m_edgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Pressure_Outlet: {

                        double[] EdgeNormal = inp.Normal;
                        double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                        double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

                        Flx_InCell = -m_sigma * Tangente_IN[m_comp];
                        //for (int d = 0; d < inp.D; d++) {
                        //    Flx_InCell -= m_sigma * (EdgeNormal[d] * Tangente_IN[d]) * EdgeNormal[m_comp];
                        //}

                        break;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear: {

                        double[] EdgeNormal = inp.Normal;
                        double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                        double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

                        int D = inp.D;

                        double[] PSnI = new double[D]; // projection of surface/level-set normal onto domain boundary tangent
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                double nn = EdgeNormal[d1] * EdgeNormal[d2];
                                if (d1 == d2) {
                                    PSnI[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                                } else {
                                    PSnI[d1] += -nn * SurfaceNormal_IN[d2];
                                }
                            }
                        }
                        double PSnINorm = PSnI.L2Norm();
                        double[] PSnINormal_IN = PSnI.Normalize(); // line normal: tangential to domain boundary & normal on contact line


                        // isotropic surface tension terms
                        for (int d = 0; d < D; d++) {
                            Flx_InCell -= m_sigma * (EdgeNormal[d] * Tangente_IN[d]) * EdgeNormal[m_comp];
                        }


                        if (edgType == IncompressibleBcType.NavierSlip_Linear) {

                            // Young's relation (static contact angle)
                            Flx_InCell -= m_sigma * Math.Cos(m_theta) * PSnINormal_IN[m_comp];

                            // dissipative contact line force
                            // beta*(u*nL)

                            double g_D = this.velFunction[inp.EdgeTag](inp.X, inp.time);

                            for(int d = 0; d < D; d++) {
                                Flx_InCell += m_beta * ((_uA[d] - g_D) * PSnINormal_IN[d]) * PSnINormal_IN[m_comp];
                            }
                        }

                        //double[] PInS = new double[D];
                        //for (int d1 = 0; d1 < D; d1++) {
                        //    for (int d2 = 0; d2 < D; d2++) {
                        //        double nn = SurfaceNormal_IN[d1] * SurfaceNormal_IN[d2];
                        //        if (d1 == d2) {
                        //            PInS[d1] += (1 - nn) * EdgeNormal[d2];
                        //        } else {
                        //            PInS[d1] += -nn * EdgeNormal[d2];
                        //        }
                        //    }
                        //}
                        //double PInSNorm = PInS.L2Norm();


                        //Flx_InCell -= (m_sigma * PInSNorm) * (EdgeNormal[m_comp]);
                        //Flx_InCell -= (m_sigma * Math.Cos(m_theta)) * (PSnINormal_IN[m_comp]); // Young's relation (static contact angle)

                        //Flx_InCell += (m_beta * _uA[0] * PSnINormal_IN[m_comp]) * (PSnINormal_IN[m_comp]); // dissipative contact line force

                        break;
                    }
                default:
                    break;
            }

            return Flx_InCell * _vA;
        }

       
        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for (int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { }; // only parameter dependent, not present in jacobian
        }
    }

    /*

    /// <summary>
    /// isotropic part for the surface stress tensor
    /// </summary>
    public class BoussinesqScriven_SurfaceTension : SurfaceFlux {

        public BoussinesqScriven_SurfaceTension(int d, double sigma) {
            m_comp = d;
            m_sigma = sigma;
        }


        int m_comp;

        double m_sigma;



        public override IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }

        public override TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin, double[,] Grad_Uin) {
            return 0;
        }

        protected override void InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout, double[,] Grad_Uin, double[,] Grad_Uout, out double Flx_InCell, out double Flx_OutCell) {

            double[] EdgeNormal = inp.Normale;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            double[] Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);


            Flx_InCell = -0.5 * (Tangente_IN[m_comp] + Tangente_OUT[m_comp]) * m_sigma;
            Flx_OutCell = 0.5 * (Tangente_IN[m_comp] + Tangente_OUT[m_comp]) * m_sigma;


        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[,] GradU, double[] output) {

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            for (int d = 0; d < inp.D; d++) {
                output[d] = -Psurf[m_comp, d] * m_sigma;
            }

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


    }

    */


    /// <summary>
    /// flux formulation for terms on the interface
    /// </summary>
    public abstract class SurfaceFluxBase : IVolumeForm, IEdgeForm, BoSSS.Foundation.IEquationComponentCoefficient {


        protected int m_D;

        protected int m_comp;

        /// <summary>
        /// interface lengths in order to determine the penalty parameter.
        /// </summary>
        MultidimensionalArray m_InterLen;

        protected double m_penalty;

        /// <summary>
        /// penalty computation
        /// </summary>
        protected double penalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor_A = (this.m_InterLen[jCellIn] > 0) ? this.m_InterLen[jCellIn] : 0;
            double penaltySizeFactor_B = (jCellOut >= 0 && this.m_InterLen[jCellOut] > 0) ? this.m_InterLen[jCellOut] : 0;
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            //throw new NotImplementedException("this penalty might be unsuitable");
            double penalty_base = (double)((m_degU + 1) * (m_degU + m_D)) / ((double)m_D);
            return this.m_penalty * penaltySizeFactor * penalty_base;

        }


        protected IncompressibleBcType[] m_edgeTag2Type;

        /// <summary>
        /// required the normal vector
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return VariableNames.NormalVector(m_D);
            }
        }

        /// <summary>
        /// the velocity vector
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public virtual TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public virtual TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.UxV;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public virtual TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.UxV;
            }
        }



        public virtual double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            int D = cpv.D;

            double[] flx = new double[D];
            this.Flux(ref cpv, GradU, flx);

            double acc = 0;
            for (int d = 0; d < D; d++)
                acc += flx[d] * GradV[d];

            return -acc;
        }

        protected abstract void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux);


        public abstract double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB);


        public abstract double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA);



        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for (int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1.0 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0.0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }

        protected static double[,] Multiply(double[,] M1, double[,] M2) {

            int D = M1.GetLength(0);
            double[,] res = new double[D, D];

            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    for (int dd = 0; dd < D; dd++) {
                        res[d1, d2] += M1[d1, dd] * M2[dd, d2];
                    }
                }
            }

            return res;
        }

        /// <summary>
        /// currently used DG degree for velocity
        /// </summary>
        protected int m_degU;

        public virtual void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_degU = DomainDGdeg[0];
            for(int i = 1; i < DomainDGdeg.Length; i++)
                if(DomainDGdeg[i] != m_degU)
                    throw new ApplicationException("something is weird with the velocity DG degrees");
            if(DomainDGdeg.Length != m_D)
                throw new ApplicationException("Spatial dimension mismatch");
            if(m_D != cs.GrdDat.SpatialDimension)
                throw new ApplicationException("Spatial dimension mismatch");

            m_InterLen = (MultidimensionalArray)cs.UserDefinedValues["InterfaceLengths"]; 
        }
    }


    /// <summary>
    /// additional dynamic part - surface velocity divergence (according to the Boussinesq-Scriven model) - for the surface stress tensor
    /// </summary>
    public class BoussinesqScriven_SurfaceVelocityDivergence : SurfaceFluxBase {


        double m_lamI;

        /// <summary>
        /// local coefficients (lambda_I) for the surface divergence term.
        /// </summary>
        MultidimensionalArray m_lamIs;

        bool updateLambdaI;


        public BoussinesqScriven_SurfaceVelocityDivergence(int d, int D, double lamI, double penalty, IncompressibleBcType[] edgeTag2Type, bool updateCoeff = false) {
            m_comp = d;
            m_D = D;
            m_lamI = lamI;
            updateLambdaI = updateCoeff;
            m_penalty = penalty;
            m_edgeTag2Type = edgeTag2Type;
        }


        protected override void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux) {

            int D = inp.D;

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            double surfDiv = 0.0;
            for (int d1 = 0; d1 < D; d1++) {
                for (int dd = 0; dd < D; dd++) {
                    surfDiv += Psurf[d1, dd] * GradU[dd, d1];
                }
            }

            double lamI = GetlambdaI(inp.jCell);

            for (int d = 0; d < D; d++) {
                flux[d] = -lamI * surfDiv * Psurf[m_comp, d];
            }

        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB,
           double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double acc = 0.0;

            int D = inp.D;

            double[] Nsurf_IN = SurfaceNormal(inp.Parameters_IN);
            double[] Nsurf_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[,] Psurf_IN = SurfaceProjection(Nsurf_IN);
            double[,] Psurf_OUT = SurfaceProjection(Nsurf_OUT);

            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normal);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normal);

            double[] tauL_Aver = tauL_IN.CloneAs();
            tauL_Aver.ScaleV(0.5);
            tauL_Aver.AccV(0.5, tauL_OUT);
            tauL_Aver.Normalize();


            double divUsurf_IN = 0.0;
            double divUsurf_OUT = 0.0;
            for (int d1 = 0; d1 < inp.D; d1++) {
                for (int dd = 0; dd < inp.D; dd++) {
                    divUsurf_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                    divUsurf_OUT += Psurf_OUT[d1, dd] * _Grad_uB[dd, d1];
                }
            }

            double lamI_IN = GetlambdaI(inp.jCellIn);
            double lamI_OUT = GetlambdaI(inp.jCellOut);

            // consistency term
            //acc += (m_lamI - m_muI) * 0.5 * (divUsurf_IN + divUsurf_OUT) * tauL_Aver[m_comp] * (_vA - _vB);
            //acc += m_lamI * 0.5 * (divUsurf_IN * tauL_IN[m_comp] + divUsurf_OUT * tauL_OUT[m_comp]) * (_vA - _vB);
            acc += 0.5 * (lamI_IN * divUsurf_IN * tauL_IN[m_comp] + lamI_OUT * divUsurf_OUT * tauL_OUT[m_comp]) * (_vA - _vB);

            // symmtery term
            //for (int d = 0; d < D; d++) {
            //    for (int dd = 0; dd < D; dd++) {
            //        //acc += (m_lamI - m_muI) * 0.5 * (Psurf_IN[m_comp,dd] * _Grad_vA[dd] + Psurf_OUT[m_comp,dd] * _Grad_vB[dd]) * tauL_Aver[d] * (_uA[d] - _uB[d]);
            //        acc += m_lamI * 0.5 * (Psurf_IN[m_comp, dd] * _Grad_vA[dd] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * _Grad_vB[dd] * tauL_OUT[d]) * (_uA[d] - _uB[d]);
            //    }
            //}
            // penalty
            //double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            //acc -= m_lamI * (_uA[m_comp] - _uB[m_comp]) * (_vA - _vB) * pnlty;

            return -acc;

        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            IncompressibleBcType edgType = m_edgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.NavierSlip_Linear: {

                        double acc = 0.0;

                        int D = inp.D;

                        double[] Nedge_IN = inp.Parameters_IN;
                        double[] Nsurf_IN = SurfaceNormal(inp.Parameters_IN);
                        double[,] Psurf_IN = SurfaceProjection(Nsurf_IN);
                        double[] tauL_IN = Tangent(Nsurf_IN, inp.Normal);

                        double divUsurf_IN = 0.0;
                        for (int d1 = 0; d1 < inp.D; d1++) {
                            for (int dd = 0; dd < inp.D; dd++) {
                                divUsurf_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                            }
                        }

                        double lamI_IN = GetlambdaI(inp.jCellIn);

                        // consistency term
                        for (int d = 0; d < D; d++) {
                            acc += lamI_IN * Nedge_IN[d] * (divUsurf_IN * tauL_IN[d]) * (_vA * Nedge_IN[m_comp]);
                        }

                        // symmtery term
                        //for (int d1 = 0; d1 < D; d1++) {
                        //    for (int d2 = 0; d2 < D; d2++) {
                        //        for (int dd = 0; dd < D; dd++) {
                        //            acc += m_lamI * Nedge_IN[d1] * (Psurf_IN[m_comp, dd] * _Grad_vA[dd] * tauL_IN[d1]) * (_uA[d2] * Nedge_IN[d2]);
                        //        }
                        //    }
                        //}

                        // penalty
                        //double pnlty = this.penalty(inp.jCellIn, -1);
                        //for (int d = 0; d < D; d++) {
                        //    acc -= m_lamI * (_uA[d] * Nedge_IN[d]) * (_vA * Nedge_IN[m_comp]) * pnlty;
                        //}


                        Flx_InCell = -acc;

                        break;
                    }
                default:
                    break;
            }

            return Flx_InCell * _vA;
        }


        public override TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.GradUxV;
            }
        }

        public override TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.GradUxV;
            }
        }


        private double GetlambdaI(int jCell) {

            if (updateLambdaI) {
                return m_lamIs[jCell];
            } else {
                return m_lamI;
            }
        }


        public override void CoefficientUpdate(CoefficientSet cs, Int32[] DomainDGdeg, Int32 TestDGdeg) {
            base.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);

            m_lamIs = (MultidimensionalArray)cs.UserDefinedValues["lambda_interface"];
        }


    }



    /// <summary>
    /// additional dynamic part - surface rate of deformation (according to the Boussinesq-Scriven model) - for the surface stress tensor
    /// </summary>
    public class BoussinesqScriven_SurfaceDeformationRate_GradU : SurfaceFluxBase {


        /// <summary>
        /// surface shear viscosity
        /// </summary>
        double m_muI;

        bool semiImplicit;


        /// <summary>
        /// local coefficients (lambda_I) for the surface divergence term.
        /// </summary>
        MultidimensionalArray m_muIs;

        bool updateMuI;


        public BoussinesqScriven_SurfaceDeformationRate_GradU(int d, int D, double mu_I, double penalty, bool updateCoeff = false, bool semiImplicitOnly = false) {
            m_comp = d;
            m_D = D;
            m_penalty = penalty;
            m_muI = mu_I;
            updateMuI = updateCoeff;
            semiImplicit = semiImplicitOnly;
        }


        protected override void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux) {

            int D = inp.D;

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            double[,] GradUsurf = Multiply(Psurf, GradU);

            double muI = GetMuI(inp.jCell);

            for (int d = 0; d < D; d++) {

                if (semiImplicit) {
                    flux[d] += -muI * GradUsurf[m_comp, d];
                    continue;
                }

                for(int dd = 0; dd < D; dd++) {
                    flux[d] += -muI * GradUsurf[m_comp, dd] * Psurf[dd, d];
                }
            }

        }


        //public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

        //    int D = cpv.D;

        //    double[] Nsurf = SurfaceNormal(cpv.Parameters);
        //    double[,] Psurf = SurfaceProjection(Nsurf);

        //    double[,] GradUsurf = Multiply(Psurf, GradU);
        //    double[,] GradUsurfPsurf = new double[D, D];
        //    for (int d = 0; d < D; d++) {
        //        for (int dd = 0; dd < D; dd++) {
        //            GradUsurfPsurf[d,dd] += GradUsurf[d, dd] * Psurf[dd, d];
        //        }
        //    }

        //    double acc = 0;
        //    for(int d = 0; d < D; d++)
        //        acc += -m_muI * GradUsurfPsurf[m_comp, d] * GradV[d];

        //    // stabilization
        //    //for(int d = 0; d < D; d++) {
        //    //    for(int dd = 0; dd < D; dd++) {
        //    //        acc += -0.1 * GradU[m_comp, d] * Nsurf[d] * GradV[dd] * Nsurf[dd];
        //    //    }
        //    //}

        //    return -acc;
        //}


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, 
            double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double acc = 0.0;

            int D = inp.D;

            double[] Nsurf_IN = SurfaceNormal(inp.Parameters_IN);
            double[] Nsurf_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[,] Psurf_IN = SurfaceProjection(Nsurf_IN);
            double[,] Psurf_OUT = SurfaceProjection(Nsurf_OUT);

            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normal);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normal);

            double[] tauL_Aver = tauL_IN.CloneAs();
            tauL_Aver.ScaleV(0.5);
            tauL_Aver.AccV(0.5, tauL_OUT);
            tauL_Aver.Normalize();


            double muI_IN = GetMuI(inp.jCellIn);
            double muI_OUT = GetMuI(inp.jCellOut);

            for (int d = 0; d < D; d++) {

                for (int dd = 0; dd < D; dd++) {
                    // consistency term
                    //acc += 0.5 * m_muI * (Psurf_IN[m_comp, dd] * GradUsurf_IN[dd, d] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * GradUsurf_OUT[dd, d] * tauL_OUT[d]) * (_vA - _vB);
                    //acc += 0.5 * m_muI * (Psurf_IN[m_comp, dd] * _Grad_uA[dd, d] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * _Grad_uB[dd, d] * tauL_OUT[d]) * (_vA - _vB);
                    acc += 0.5 * (muI_IN * Psurf_IN[m_comp, dd] * _Grad_uA[dd, d] * tauL_IN[d] + muI_OUT * Psurf_OUT[m_comp, dd] * _Grad_uB[dd, d] * tauL_OUT[d]) * (_vA - _vB);

                    // symmetry term
                    //acc += 0.5 * m_muI * (Psurf_IN[d, m_comp] * _Grad_vA[dd] * tauL_IN[dd] + Psurf_OUT[d, m_comp] * _Grad_vB[dd] * tauL_OUT[dd]) * (_uA[d] - _uB[d]);
                    acc += 0.5 * (muI_IN * Psurf_IN[d, m_comp] * _Grad_vA[dd] * tauL_IN[dd] + muI_OUT * Psurf_OUT[d, m_comp] * _Grad_vB[dd] * tauL_OUT[dd]) * (_uA[d] - _uB[d]);
                }
                // symmetry term
                //for (int d2 = 0; d2 < D; d2++) {
                //    //acc += 0.5 * m_muI * (GradVsurf_IN[m_comp, d] * tauL_IN[d] + GradVsurf_OUT[m_comp, d] * tauL_OUT[d]) * (_uA[m_comp] - _uB[m_comp]);
                //    //acc += 0.5 * m_muI * (_Grad_vA[d] * tauL_IN[d] + _Grad_vB[d] * tauL_OUT[d]) * (Psurf_IN[m_comp, d2] * _uA[d2] - Psurf_OUT[m_comp, d2] * _uB[d2]);
                //    acc += 0.5 * m_muI * (_Grad_vA[d] * tauL_IN[d] + _Grad_vB[d] * tauL_OUT[d]) * (Psurf_IN[m_comp, d2] * _uA[d2] - Psurf_OUT[m_comp, d2] * _uB[d2]);
                //}
            }
            // penalty
            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            acc -= Math.Max(muI_IN, muI_OUT) * (_uA[m_comp] - _uB[m_comp]) * (_vA - _vB) * pnlty;

            return -acc;

        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }



        private double GetMuI(int jCell) {

            if (updateMuI) {
                return m_muIs[jCell];
            } else {
                return m_muI;
            }
        }


        public override void CoefficientUpdate(CoefficientSet cs, Int32[] DomainDGdeg, Int32 TestDGdeg) {
            base.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);

            m_muIs = (MultidimensionalArray)cs.UserDefinedValues["mu_interface"];
        }

    }


    /// <summary>
    /// additional dynamic part - surface rate of deformation, transposed term (according to the Boussinesq-Scriven model) - for the surface stress tensor
    /// </summary>
    public class BoussinesqScriven_SurfaceDeformationRate_GradUTranspose : SurfaceFluxBase {


        /// <summary>
        /// surface shear viscosity
        /// </summary>
        double m_muI;

        /// <summary>
        /// local coefficients (lambda_I) for the surface divergence term.
        /// </summary>
        MultidimensionalArray m_muIs;

        bool updateMuI;


        public BoussinesqScriven_SurfaceDeformationRate_GradUTranspose(int d, int D, double mu_I, double penalty, bool updateCoeff = false) {
            m_comp = d;
            m_D = D;
            m_muI = mu_I;
            updateMuI = updateCoeff;
            m_penalty = penalty;
        }
        

        protected override void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux) {

            int D = inp.D;

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            double[,] GradUTsurf = new double[D, D];
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    for (int dd = 0; dd < D; dd++) {
                        GradUTsurf[d1, d2] += Psurf[d1, dd] * GradU[d2, dd];
                    }
                }
            }

            double muI = GetMuI(inp.jCell);

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    flux[d] += -muI * GradUTsurf[m_comp, dd] * Psurf[dd, d];
                }
            }

        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB,
            double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double acc = 0.0;

            int D = inp.D;

            double[] Nsurf_IN = SurfaceNormal(inp.Parameters_IN);
            double[] Nsurf_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[,] Psurf_IN = SurfaceProjection(Nsurf_IN);
            double[,] Psurf_OUT = SurfaceProjection(Nsurf_OUT);


            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normal);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normal);

            double[] tauL_Aver = tauL_IN.CloneAs();
            tauL_Aver.ScaleV(0.5);
            tauL_Aver.AccV(0.5, tauL_OUT);
            tauL_Aver.Normalize();

            double muI_IN = GetMuI(inp.jCellIn);
            double muI_OUT = GetMuI(inp.jCellOut);

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    // consistency term
                    //acc += 0.5 * m_muI * (GradUTsurf_IN[m_comp, d] * tauL_IN[d] + GradUTsurf_OUT[m_comp, d] * tauL_OUT[d]) * (_vA - _vB);
                    //acc += 0.5 * m_muI * (Psurf_IN[m_comp, dd] * _Grad_uA[d, dd] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * _Grad_uB[d, dd] * tauL_OUT[d]) * (_vA - _vB);
                    acc += 0.5 * (muI_IN * Psurf_IN[m_comp, dd] * _Grad_uA[d, dd] * tauL_IN[d] + muI_OUT * Psurf_OUT[m_comp, dd] * _Grad_uB[d, dd] * tauL_OUT[d]) * (_vA - _vB);

                    // symmetry term 
                    //acc += 0.5 * m_muI * (GradVTsurf_IN[d] * tauL_IN[m_comp] + GradVTsurf_OUT[d] * tauL_OUT[m_comp]) * (_uA[d] - _uB[d]);
                    //acc += 0.5 * m_muI * (Psurf_IN[d, dd] * _Grad_vA[dd] * tauL_IN[m_comp] + Psurf_OUT[d, dd] * _Grad_vB[dd] * tauL_OUT[m_comp]) * (_uA[d] - _uB[d]);
                    acc += 0.5 * (muI_IN * Psurf_IN[d, dd] * _Grad_vA[dd] * tauL_IN[m_comp] + muI_OUT * Psurf_OUT[d, dd] * _Grad_vB[dd] * tauL_OUT[m_comp]) * (_uA[d] - _uB[d]);
                }
            }
            // penalty
            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            acc -= Math.Max(muI_IN, muI_OUT) * (_uA[m_comp] - _uB[m_comp]) * (_vA - _vB) * pnlty;

            return -acc;

        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }



        private double GetMuI(int jCell) {

            if (updateMuI) {
                return m_muIs[jCell];
            } else {
                return m_muI;
            }
        }


        public override void CoefficientUpdate(CoefficientSet cs, Int32[] DomainDGdeg, Int32 TestDGdeg) {
            base.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);

            m_muIs = (MultidimensionalArray)cs.UserDefinedValues["mu_interface"];
        }

    }



    public class SurfaceDeformationRate_LocalStabilization : IVolumeForm, BoSSS.Foundation.IEquationComponentCoefficient {

        int m_comp;
        int m_D;

        bool m_transposed;
  
        public SurfaceDeformationRate_LocalStabilization(int d, int D, bool transposed = false) {
            m_comp = d;
            m_D = D;
            m_transposed = transposed;
        }

        /// <summary>
        /// local stabilization coefficients.
        /// </summary>
        MultidimensionalArray m_muIs;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {

            m_muIs = (MultidimensionalArray)cs.UserDefinedValues["mu_interface"];
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return VariableNames.NormalVector(m_D);
            }
        }

        public virtual IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }


        public virtual TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        public virtual double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            int D = cpv.D;

            double[] flx = new double[D];
            this.Flux(ref cpv, GradU, flx);

            double acc = 0;
            for (int d = 0; d < D; d++)
                acc += flx[d] * GradV[d];

            return -acc;
        }

        protected void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux) {

            int D = inp.D;

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            double[,] GradUsurf = Multiply(Psurf, GradU);

            for (int d = 0; d < D; d++) {

                for (int dd = 0; dd < D; dd++) {
                    flux[d] += -m_muIs[inp.jCell] * GradUsurf[m_comp, dd] * Psurf[dd, d];
                    // transpose term
                    if(m_transposed)
                        flux[d] += -m_muIs[inp.jCell] * GradUsurf[dd, m_comp] * Psurf[dd, d];
                }
            }

        }

        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for (int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1.0 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0.0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }

        protected static double[,] Multiply(double[,] M1, double[,] M2) {

            int D = M1.GetLength(0);
            double[,] res = new double[D, D];

            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    for (int dd = 0; dd < D; dd++) {
                        res[d1, d2] += M1[d1, dd] * M2[dd, d2];
                    }
                }
            }

            return res;
        }

    }



    public class LevelSetStabilization : ILevelSetForm {

        int m_D;
        int m_d;

        double m_penalty;

        //LevelSetTracker m_LsTrk;

        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public LevelSetStabilization(int _d, int _D, double penalty) {
            //m_LsTrk = LsTrk;
            if(_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;

            this.m_penalty = penalty;
        }


        public double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.Normal;

            double acc = 0.0;

            for(int d = 0; d < m_D; d++) {
                for(int dd = 0; dd < m_D; dd++) {
                    acc += -(Grad_uA[m_d, d] - Grad_uB[m_d, d]) * Normal[d] * (Grad_vA[dd] - Grad_vB[dd]) * Normal[dd];
                }
            }

            return - m_penalty * acc;
        }


        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }


        public virtual IList<string> ParameterOrdering {
            get { return new string[] { }; }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        
    }



    public class DynamicSurfaceTension_LB_SurfaceVelocityDivergence : IEdgeForm {

        int m_comp;
        int m_D;

        double m_penalty;

        //IncompressibleBcType[] m_edgeTag2Type;

        public DynamicSurfaceTension_LB_SurfaceVelocityDivergence(int d, int D, double penalty) {    //, IncompressibleBcType[] edgeTag2Type) {
            m_comp = d;
            m_D = D;
            m_penalty = penalty;
            //m_edgeTag2Type = edgeTag2Type;
        }


        public IList<string> ParameterOrdering {
            get {
                return VariableNames.NormalVector(m_D);
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB,
           double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double acc = 0.0;

            int D = inp.D;

            double[] Nsurf_IN = SurfaceNormal(inp.Parameters_IN);
            double[] Nsurf_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[,] Psurf_IN = SurfaceProjection(Nsurf_IN);
            double[,] Psurf_OUT = SurfaceProjection(Nsurf_OUT);

            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normal);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normal);


            double divUsurf_IN = 0.0;
            double divUsurf_OUT = 0.0;
            for (int d1 = 0; d1 < inp.D; d1++) {
                for (int dd = 0; dd < inp.D; dd++) {
                    divUsurf_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                    divUsurf_OUT += Psurf_OUT[d1, dd] * _Grad_uB[dd, d1];
                }
            }

            // penalty
            double pnlty = m_penalty; // this.penalty(inp.jCellIn, inp.jCellOut);
            acc -= (divUsurf_IN * tauL_IN[m_comp] - divUsurf_OUT * tauL_OUT[m_comp]) * (_vA - _vB) * pnlty;

            return -acc;

        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            //IncompressibleBcType edgType = m_edgeTag2Type[inp.EdgeTag];

            return Flx_InCell * _vA;
        }


        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.GradUxV;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.None;
            }
        }


        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for (int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1.0 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0.0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for (int d1 = 0; d1 < D; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if (d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }


    }



    public class DynamicSurfaceTension_LB_EdgeDissipation : IEdgeForm {

        int m_comp;
        int m_D;

        double m_sigma;
        double m_beta;

        //IncompressibleBcType[] m_edgeTag2Type;


        public DynamicSurfaceTension_LB_EdgeDissipation(int d, int D, double sigma, double beta) {     //, IncompressibleBcType[] edgeTag2Type) {
            m_comp = d;
            m_D = D;
            m_sigma = sigma;
            m_beta = beta;
            //m_edgeTag2Type = edgeTag2Type;
        }


        public virtual IList<string> ParameterOrdering {
            get {
                return VariableNames.NormalVector(m_D);
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }


        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.None; 
            }
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            double[] Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);


            double[] PEnI_IN = new double[m_D];
            double[] PEnI_OUT = new double[m_D];
            for (int d1 = 0; d1 < m_D; d1++) {
                for (int d2 = 0; d2 < m_D; d2++) {
                    double nn = EdgeNormal[d1] * EdgeNormal[d2];
                    if (d1 == d2) {
                        PEnI_IN[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                        PEnI_OUT[d1] += (1 - nn) * SurfaceNormal_OUT[d2];
                    } else {
                        PEnI_IN[d1] += -nn * SurfaceNormal_IN[d2];
                        PEnI_OUT[d1] += -nn * SurfaceNormal_OUT[d2];
                    }
                }
            }

            //double acc = 0.5 * (Tangente_IN[m_comp] + Tangente_OUT[m_comp]) * m_sigma * (_vA - _vB);
            double acc = 0.0; // (PEnI_IN[m_comp] - PEnI_OUT[m_comp]) * (_vA - _vB);

            double Flx_InCell = m_sigma * PEnI_IN[m_comp];
            double Flx_OutCell = m_sigma * PEnI_OUT[m_comp];

            PEnI_IN.Normalize();
            PEnI_OUT.Normalize();

            for (int d = 0; d < m_D; d++) {
                Flx_InCell += m_beta * (_uA[d] * PEnI_IN[d]) * PEnI_IN[m_comp];
                Flx_InCell += m_beta * (_uB[d] * PEnI_OUT[d]) * PEnI_OUT[m_comp];
            }

            acc = (Flx_InCell - Flx_OutCell) * (_vA - _vB);

            return -acc;

        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            //IncompressibleBcType edgType = m_edgeTag2Type[inp.EdgeTag];

            return Flx_InCell * _vA;
        }


        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for(int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for(int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    if(dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if(d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }


    }

}
