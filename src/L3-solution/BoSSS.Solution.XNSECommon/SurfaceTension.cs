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


    public class CurvatureBasedSurfaceTension : ILevelSetForm {

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

        LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public CurvatureBasedSurfaceTension(int _d, int _D, LevelSetTracker LsTrk, double _sigma) {
            m_LsTrk = LsTrk;
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
        

        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
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
                return new string[] { "Curvature", 
                //    "NX", "NY" 
                };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }
    }


    /// <summary>
    /// Represents the artificial surface force (usually only used in manufactured solutions).
    /// </summary>
    public class SurfaceTension_ArfForceSrc  : ILevelSetForm {

        public static double hmin = double.NaN;

        LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        public SurfaceTension_ArfForceSrc(int _d, int _D, LevelSetTracker LsTrk) {
            m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;
       }

        int m_D;
        int m_d;

        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            //throw new NotImplementedException();
            double curvature = cp.ParamsPos[0];
            Debug.Assert(cp.ParamsPos[0] == cp.ParamsNeg[0], "curvature must be continuous across interface");
            Debug.Assert(!double.IsNaN(curvature) || !double.IsInfinity(curvature));

            double[] Normal = cp.n;

            double surfForce = cp.ParamsNeg[0];
            Debug.Assert(cp.ParamsNeg[0] == cp.ParamsPos[0]);

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
                return new string[] { (new string[] { "surfForceX", "surfForceY", "surfForceZ" })[this.m_d] };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
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

        public override IList<string> ArgumentOrdering {
            get {
                return new string[0]; 
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return new string[] { "NX", "NY" };
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
                return new string[] { "NX", "NY" };
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
            double[] EdgeNormal = inp.Normale;
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
            double[] EdgeNormal = inp.Normale;
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
    }

    /*

    /// <summary>
    /// surface tension force, in Laplace-Beltrami -- form, surface-integral part (must be used in conjunction with the boundary-line-integral part <see cref="SurfaceTension_LaplaceBeltrami_BndLine"/>).
    /// </summary>
    public class SurfaceTension_LaplaceBeltrami2_Surface : LinearFlux {


        public SurfaceTension_LaplaceBeltrami2_Surface(int d, double factor) {
            m_comp = d;
            m_factor = factor;
        }

        int m_comp;
        double m_factor;

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            return 0;
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return 0;
        }

        static double[] SurfaceNormal(ref CommonParamsVol inp) {

            double[] N = new double[inp.D];

            for (int d = 0; d < inp.D; d++) {
                N[d] = inp.Parameters[d];
            }

            return N.Normalize();
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {

            double[] Nsurf = SurfaceNormal(ref inp);

            double[,] Psurf = new double[inp.D, inp.D];

            for (int d = 0; d < inp.D; d++) {
                for (int dd = 0; dd < inp.D; dd++) {
                    if (dd == d)
                        Psurf[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        Psurf[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            for (int d = 0; d < inp.D; d++) {
                output[d] = -Psurf[m_comp, d] * m_factor;
            }

        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return new string[] { "NX", "NY" };
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
    public class SurfaceTension_LaplaceBeltrami2_BndLine : LinearDualValueFlux {

        public SurfaceTension_LaplaceBeltrami2_BndLine(int d, double factor, double Theta_e = Math.PI / 2.0, double beta_L = 0.0) {
            m_comp = d;
            m_factor = factor;
            m_theta = Theta_e;
            m_beta = beta_L;
        }

        int m_comp;
        double m_factor;

        double m_theta;
        double m_beta;

        public override IList<string> ParameterOrdering {
            get {
                return new string[] { "NX", "NY" };
            }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_comp) };
            }
        }

        static double[] SurfaceNormal(int D, double[] param) {

            double[] NS = new double[D];

            for (int d = 0; d < D; d++) {
                NS[d] = param[d];
            }

            return NS.Normalize();
        }

        static double[] Tangent(int D, double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

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


        protected override void InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOuCell) {

            double[] EdgeNormal = inp.Normale;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.D, inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.D, inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(inp.D, SurfaceNormal_IN, EdgeNormal);
            double[] Tangente_OUT = Tangent(inp.D, SurfaceNormal_OUT, EdgeNormal);

            FluxInCell = -Tangente_IN[m_comp] * m_factor;
            FluxOuCell = +Tangente_OUT[m_comp] * m_factor;
        }



        protected override void BorderEdgeFlux_(ref CommonParamsBnd inp, double[] Uin, out double FluxInCell) {

            //IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            //switch (edgType) { }

            double[] EdgeNormal = inp.Normale;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.D, inp.Parameters_IN);
            double[] Tangente_IN = Tangent(inp.D, SurfaceNormal_IN, EdgeNormal);

            //double Theta = Math.Atan2(SurfaceNormal_IN[1], SurfaceNormal_IN[0]) * (180 / Math.PI);
            //Console.WriteLine("Theta = {0}", Theta);

            double[] PSnI = new double[inp.D];
            for (int d1 = 0; d1 < inp.D; d1++) {
                for (int d2 = 0; d2 < inp.D; d2++) {
                    double nn = EdgeNormal[d1] * EdgeNormal[d2];
                    if (d1 == d2) {
                        PSnI[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                    } else {
                        PSnI[d1] += -nn * SurfaceNormal_IN[d2];
                    }
                }
            }
            double PSnINorm = PSnI.L2Norm();
            double[] PSnINormal_IN = PSnI.Normalize();

            //FluxInCell = -(EdgeNormal[m_comp] * m_factor * Tangente_IN[m_comp]) * (EdgeNormal[m_comp]);
            //FluxInCell = -(m_factor * Math.Cos(m_theta)) * (CLineNormal_IN[m_comp]); // Young's relation (static contact angle)

            double[] PInS = new double[inp.D];
            for (int d1 = 0; d1 < inp.D; d1++) {
                for (int d2 = 0; d2 < inp.D; d2++) {
                    double nn = SurfaceNormal_IN[d1] * SurfaceNormal_IN[d2];
                    if (d1 == d2) {
                        PInS[d1] += (1 - nn) * EdgeNormal[d2];
                    } else {
                        PInS[d1] += -nn * EdgeNormal[d2];
                    }
                }
            }
            double PInSNorm = PInS.L2Norm();

            FluxInCell = 0;

            FluxInCell -= (m_factor * PInSNorm) * (EdgeNormal[m_comp]);
            FluxInCell -= (m_factor * Math.Cos(m_theta)) * (PSnINormal_IN[m_comp]); // Young's relation (static contact angle)

            FluxInCell += (m_beta * Uin[0] * PSnINormal_IN[m_comp]) * (PSnINormal_IN[m_comp]); // dissipative contact line force
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


    public class IsotropicSurfaceTension_LaplaceBeltrami : IVolumeForm, IEdgeForm {

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

        IncompressibleBcType[] m_edgeTag2Type;


        public IsotropicSurfaceTension_LaplaceBeltrami(int d, int D, double sigma, IncompressibleBcType[] edgeTag2Type, double theta_e, double beta_L) {
            m_comp = d;
            m_D = D;
            m_sigma = sigma;
            m_theta = theta_e;
            m_beta = beta_L;
            m_edgeTag2Type = edgeTag2Type;
        }


        public virtual IList<string> ParameterOrdering {
            get {
                switch(m_D) {
                    case 2:
                        return new string[] { "NX", "NY" };
                    case 3:
                        return new string[] { "NX", "NY", "NZ" };
                    default:
                        return new string[] { };
                }
                //return new string[] { "NX", "NY" };
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
                //return new string[] { VariableNames.Velocity_d(0), VariableNames.Velocity_d(1) };
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

            for (int d = 0; d < cpv.D; d++)
                acc += - m_sigma * Psurf[m_comp, d] * GradV[d];

            return -acc;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double[] EdgeNormal = inp.Normale;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            double[] Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);

            double acc = 0.5 * (Tangente_IN[m_comp] + Tangente_OUT[m_comp]) * m_sigma * (_vA - _vB);

            return -acc;
        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            IncompressibleBcType edgType = m_edgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.NavierSlip_Linear: {

                        double[] EdgeNormal = inp.Normale;
                        double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                        double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

                        int D = inp.D;

                        double[] PSnI = new double[D];
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
                        double[] PSnINormal_IN = PSnI.Normalize();


                        // isotropic surface tension terms
                        for (int d = 0; d < D; d++) {
                            Flx_InCell -= m_sigma * (EdgeNormal[d] * Tangente_IN[d]) * EdgeNormal[m_comp];
                        }

                        // Young's relation (static contact angle)
                        Flx_InCell -= m_sigma * Math.Cos(m_theta) * PSnINormal_IN[m_comp];

                        // dissipative contact line force
                        for (int d = 0; d < D; d++) {
                            Flx_InCell += m_beta * (_uA[d] * PSnINormal_IN[d]) * PSnINormal_IN[m_comp];
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


        protected int m_comp;

        /// <summary>
        /// interface lengths in order to determine the penalty parameter.
        /// </summary>
        MultidimensionalArray m_InterLen;

        public void SetParameter(string speciesName, SpeciesId SpcId, MultidimensionalArray __InterLen) {
            this.m_InterLen = __InterLen;
        }

        protected double m_penalty;

        protected double penalty(int jCellIn, int jCellOut) {

            double muFactor = 1.0;
            double penaltySizeFactor_A = (this.m_InterLen[jCellIn] > 0) ? 1.0 / this.m_InterLen[jCellIn] : 0;
            double penaltySizeFactor_B = (jCellOut >= 0 && this.m_InterLen[jCellOut] > 0) ? 1.0 / this.m_InterLen[jCellOut] : 0;
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            //throw new NotImplementedException("this penalty might be unsuitable");
            return this.m_penalty * penaltySizeFactor * muFactor;

        }


        protected IncompressibleBcType[] m_edgeTag2Type;


        public virtual IList<string> ParameterOrdering {
            get {
                return new string[] { "NX", "NY" };
            }
        }

        public virtual IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(0), VariableNames.Velocity_d(1) };
            }
        }


        public virtual TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        public virtual TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.UxV;
            }
        }

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

        MultidimensionalArray m_LenScales;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            throw new NotImplementedException("this penalty might be unsuitable");
            m_LenScales = cs.CellLengthScales;
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


        public BoussinesqScriven_SurfaceDeformationRate_GradU(int d, double mu_I, double penalty) {
            m_comp = d;
            m_penalty = penalty;
            m_muI = mu_I;
        }


        protected override void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux) {

            int D = inp.D;

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            double[,] GradUsurf = Multiply(Psurf, GradU);

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    flux[d] += -m_muI * GradUsurf[m_comp, dd] * Psurf[dd, d];
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

            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normale);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normale);

            double[] tauL_Aver = tauL_IN.CloneAs();
            tauL_Aver.ScaleV(0.5);
            tauL_Aver.AccV(0.5, tauL_OUT);
            tauL_Aver.Normalize();


            for (int d = 0; d < D; d++) {

                for (int dd = 0; dd < D; dd++) {
                    // consistency term
                    //acc += 0.5 * m_muI * (Psurf_IN[m_comp, dd] * GradUsurf_IN[dd, d] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * GradUsurf_OUT[dd, d] * tauL_OUT[d]) * (_vA - _vB);
                    acc += 0.5 * m_muI * (Psurf_IN[m_comp, dd] * _Grad_uA[dd, d] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * _Grad_uB[dd, d] * tauL_OUT[d]) * (_vA - _vB);

                    // symmetry term
                    acc += 0.5 * m_muI * (Psurf_IN[d, m_comp] * _Grad_vA[dd] * tauL_IN[dd] + Psurf_OUT[d, m_comp] * _Grad_vB[dd] * tauL_OUT[dd]) * (_uA[d] - _uB[d]);
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
            acc -= m_muI * (_uA[m_comp] - _uB[m_comp]) * (_vA - _vB) * pnlty;

            return -acc;

        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
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


        public BoussinesqScriven_SurfaceDeformationRate_GradUTranspose(int d, double mu_I, double penalty) {
            m_comp = d;
            m_muI = mu_I;
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

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    flux[d] += -m_muI * GradUTsurf[m_comp, dd] * Psurf[dd, d];
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


            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normale);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normale);

            double[] tauL_Aver = tauL_IN.CloneAs();
            tauL_Aver.ScaleV(0.5);
            tauL_Aver.AccV(0.5, tauL_OUT);
            tauL_Aver.Normalize();


            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    // consistency term
                    //acc += 0.5 * m_muI * (GradUTsurf_IN[m_comp, d] * tauL_IN[d] + GradUTsurf_OUT[m_comp, d] * tauL_OUT[d]) * (_vA - _vB);
                    acc += 0.5 * m_muI * (Psurf_IN[m_comp, dd] * _Grad_uA[d, dd] * tauL_IN[d] + Psurf_OUT[m_comp, dd] * _Grad_uB[d, dd] * tauL_OUT[d]) * (_vA - _vB);

                    // symmetry term 
                    //acc += 0.5 * m_muI * (GradVTsurf_IN[d] * tauL_IN[m_comp] + GradVTsurf_OUT[d] * tauL_OUT[m_comp]) * (_uA[d] - _uB[d]);
                    acc += 0.5 * m_muI * (Psurf_IN[d, dd] * _Grad_vA[dd] * tauL_IN[m_comp] + Psurf_OUT[d, dd] * _Grad_vB[dd] * tauL_OUT[m_comp]) * (_uA[d] - _uB[d]);
                }
            }
            // penalty
            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            acc -= m_muI * (_uA[m_comp] - _uB[m_comp]) * (_vA - _vB) * pnlty;

            return -acc;

        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }


    }



    /// <summary>
    /// additional dynamic part - surface velocity divergence (according to the Boussinesq-Scriven model) - for the surface stress tensor
    /// </summary>
    public class BoussinesqScriven_SurfaceVelocityDivergence : SurfaceFluxBase {


        /// <summary>
        /// surface shear viscosity
        /// </summary>
        double m_muI;

        /// <summary>
        /// surface dilatational viscosity
        /// </summary>
        double m_lamI;


        public BoussinesqScriven_SurfaceVelocityDivergence(int d, double mu_I, double lam_I, double penalty, IncompressibleBcType[] edgeTag2Type) {
            m_comp = d;
            m_muI = mu_I;
            m_lamI = lam_I;
            m_penalty = penalty;
            m_edgeTag2Type = edgeTag2Type;
        }


        protected override void Flux(ref CommonParamsVol inp, double[,] GradU, double[] flux) {

            int D = inp.D;

            double[] Nsurf = SurfaceNormal(inp.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            double divUsurf = 0.0;
            for (int d1 = 0; d1 < D; d1++) {
                for (int dd = 0; dd < D; dd++) {
                    divUsurf += Psurf[d1, dd] * GradU[dd, d1];
                }
            }

            for (int d = 0; d < D; d++) {
                flux[d] = -(m_lamI - m_muI) * divUsurf * Psurf[m_comp, d];
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

            double[] tauL_IN = Tangent(Nsurf_IN, inp.Normale);
            double[] tauL_OUT = Tangent(Nsurf_OUT, inp.Normale);

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


            // consistency term
            //acc += (m_lamI - m_muI) * 0.5 * (divUsurf_IN + divUsurf_OUT) * tauL_Aver[m_comp] * (_vA - _vB);
            acc += (m_lamI - m_muI) * 0.5 * (divUsurf_IN * tauL_IN[m_comp] + divUsurf_OUT * tauL_OUT[m_comp]) * (_vA - _vB);

            // symmtery term
            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    //acc += (m_lamI - m_muI) * 0.5 * (Psurf_IN[m_comp,dd] * _Grad_vA[dd] + Psurf_OUT[m_comp,dd] * _Grad_vB[dd]) * tauL_Aver[d] * (_uA[d] - _uB[d]);
                    acc += (m_lamI - m_muI) * 0.5 * (Psurf_IN[m_comp,dd] * _Grad_vA[dd] * tauL_IN[d] + Psurf_OUT[m_comp,dd] * _Grad_vB[dd] * tauL_OUT[d]) * (_uA[d] - _uB[d]);
                }
            }
            // penalty
            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);
            acc -= (m_lamI - m_muI) * (_uA[m_comp] - _uB[m_comp]) * (_vA - _vB) * pnlty;

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
                        double[] tauL_IN = Tangent(Nsurf_IN, inp.Normale);

                        double divUsurf_IN = 0.0;
                        for (int d1 = 0; d1 < inp.D; d1++) {
                            for (int dd = 0; dd < inp.D; dd++) {
                                divUsurf_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                            }
                        }

                        // consistency term
                        for (int d = 0; d < D; d++) {
                            acc += (m_lamI - m_muI) * Nedge_IN[d] * (divUsurf_IN * tauL_IN[d]) * (_vA * Nedge_IN[m_comp]);
                        }

                        // symmtery term
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                for (int dd = 0; dd < D; dd++) {
                                    acc += (m_lamI - m_muI) * Nedge_IN[d1] * (Psurf_IN[m_comp, dd] * _Grad_vA[dd] * tauL_IN[d1]) * (_uA[d2] * Nedge_IN[d2]);
                                }
                            }
                        }

                        // penalty
                        double pnlty = this.penalty(inp.jCellIn, -1);
                        for (int d = 0; d < D; d++) {
                            acc -= (m_lamI - m_muI) * (_uA[d] * Nedge_IN[d]) * (_vA * Nedge_IN[m_comp]) * pnlty;
                        }


                        Flx_InCell = -acc;

                        break;
                    }
                default:
                    break;
            }

            return Flx_InCell * _vA;
        }


    }

}
