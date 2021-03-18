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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution.RheologyCommon {

    /// <summary>
    /// Convective part of constitutive equations for singlephase flow.
    /// </summary>
    public class ConstitutiveEqns_Convective : IVolumeForm, IEdgeForm, IEquationComponentCoefficient, ISupportsJacobianComponent {
        int Component;
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        protected double m_Weissenberg; // Weissenberg number
        protected double m_alpha; // upwind-paramter
        protected bool m_UseFDJacobian;

        /// <summary>
        /// Mapping from edge tags to boundary values.
        /// - 1st index: edge tag;
        /// - 2nd index: spatial direction, row
        /// - 3rd index: spatial direction, column
        /// </summary>
        protected Func<double[], double, double>[,,] StressFunction;

        ///// <summary>
        ///// Mapping from edge tags to boundary values.<br/>
        ///// 1st index: edge tag;<br/>
        ///// 2nd index: spatial direction
        ///// </summary>
        //protected Func<double[], double, double>[,] velFunction;

        public ConstitutiveEqns_Convective(int _Component, BoundaryCondMap<IncompressibleBcType> _BcMap, double Weissenberg, bool UseFDJacobian, double alpha = 1.0) {
            Component = _Component;
            this.m_BcMap = _BcMap;
            this.m_Weissenberg = Weissenberg;
            this.m_alpha = alpha;
            this.m_UseFDJacobian = UseFDJacobian;

            //StressFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2, 2];


            //var stressXXfuncS = m_BcMap.bndFunction[VariableNames.StressXX];
            //var stressXYfuncS = m_BcMap.bndFunction[VariableNames.StressXY];
            //var stressYYfuncS = m_BcMap.bndFunction[VariableNames.StressYY];

            //for (int et = 0; et < GridCommons.FIRST_PERIODIC_BC_TAG; et++) {
            //    StressFunction[et, 0, 0] = stressXXfuncS[et];
            //    StressFunction[et, 1, 0] = stressXYfuncS[et];
            //    StressFunction[et, 0, 1] = stressXYfuncS[et];
            //    StressFunction[et, 1, 1] = stressYYfuncS[et];
            //}

            //velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];
            //for (int d = 0; d < 2; d++)
            //    velFunction.SetColumn(m_BcMap.bndFunction[VariableNames.Velocity_d(d)], d);
        }


        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV | TermActivationFlags.UxV; }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public IList<string> ArgumentOrdering {
            get {
                string[] R;

                switch (Component) {
                    case 0:
                    R = new string[] { VariableNames.StressXX };
                    break;
                    case 1:
                    R = new string[] { VariableNames.StressXY };
                    break;
                    case 2:
                    R = new string[] { VariableNames.StressYY };
                    break;
                    default:
                    throw new NotImplementedException();
                }

                R = R.Cat(VariableNames.VelocityVector(2));

                return R;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (m_UseFDJacobian) {
                    return VariableNames.Velocity0Vector(2);
                } else {
                    return null;
                }
            }
        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _Tin, double[,] Grad_Tin, double Vin, double[] Grad_Vin) {
            Vector Normale = inp.Normal;
            double Tout;
            double Tin = _Tin[0];
            

            if (m_BcMap.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.Wall || m_BcMap.EdgeTag2Type[inp.EdgeTag] == IncompressibleBcType.FreeSlip) {
                Tout = Tin;
            } else {

                switch (Component) {
                    case 0:
                        Tout = Tin;// StressFunction[inp.EdgeTag, 0, 0](inp.X, inp.time) * m_Weissenberg; // stress_XX
                        break;
                    case 1:
                        Tout = Tin;// StressFunction[inp.EdgeTag, 0, 1](inp.X, inp.time); // stress_XY
                        break;
                    case 2:
                        Tout = Tin;//StressFunction[inp.EdgeTag, 1, 1](inp.X, inp.time); // stress_YY
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }

            //Flux In
            double res1 = 0;
            double flxIn = 0;
            double n_u1 = 0;

            if (m_UseFDJacobian) {
                for (int d = 0; d < 2; d++) {
                    n_u1 += Normale[d] * 0.5 * (inp.Parameters_IN[d] + inp.Parameters_IN[d]);
                }
            } else {
                Vector Velocity_IN = new Vector(_Tin, 1, inp.D);
                n_u1 = Normale * Velocity_IN;
            }

            double factor;
            if (n_u1 < 0)
                factor = this.m_alpha;
            else
                factor = 1.0 - this.m_alpha;

            res1 += factor * n_u1 * (Tout - Tin);
            flxIn = m_Weissenberg * res1 * Vin;

            return flxIn;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] Tin, double[] Tout, double[,] Grad_Tin, double[,] Grad_Tout, double Vin, double Vout, double[] Grad_Vin, double[] Grad_Vout) {

            Vector Normale = inp.Normal;

            //Flux In
            double res1 = 0;
            double flxIn = 0;
            double n_u1 = 0;

            if (m_UseFDJacobian) {
                for (int d = 0; d < 2; d++) {
                    n_u1 += Normale[d] * 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]);
                }
            } else {
                Vector Velocity_IN = new Vector(Tin, 1, inp.D);
                Vector Velocity_OT = new Vector(Tout, 1, inp.D);
                n_u1 = 0.5 * (Normale * (Velocity_IN + Velocity_OT));
            }

            double factor;
            if (n_u1 < 0)
                factor = this.m_alpha;
            else
                factor = 1.0 - this.m_alpha;

            res1 += factor * n_u1 * (Tout[0] - Tin[0]);
            flxIn = m_Weissenberg * res1 * Vin;


            //Flux out
            double res2 = 0;
            double flxOut = 0;
            double n_u2 = 0;

            Normale.Scale(-1.0);

            if (m_UseFDJacobian) {
                for (int d = 0; d < 2; d++) {
                    n_u2 += Normale[d] * 0.5 * (inp.Parameters_OUT[d] + inp.Parameters_IN[d]);
                }
            } else {
                Vector Velocity_IN = new Vector(Tin, 1, inp.D);
                Vector Velocity_OT = new Vector(Tout, 1, inp.D);
                n_u2 = 0.5 * (Normale * (Velocity_IN + Velocity_OT));
            }

            double factor2;
            if (n_u2 < 0)
                factor2 = this.m_alpha;
            else
                factor2 = 1.0 - this.m_alpha;

            res2 += factor2 * n_u2 * (Tin[0] - Tout[0]);
            flxOut = m_Weissenberg * res2 * Vout;

            return flxIn + flxOut;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] T, double[,] GradT, double V, double[] GradV) {

            if (m_UseFDJacobian) {
                double res = 0.0;
                res += cpv.Parameters[0] * GradT[0, 0] + cpv.Parameters[1] * GradT[0, 1];
                return m_Weissenberg * res * V;
            } else {
                Vector Velocity = new Vector(T, 1, cpv.D);
                Vector _GradT = GradT.GetRowPt(0);
                double res = Velocity * _GradT;
                return m_Weissenberg * res * V;
            }
        }

        /// <summary>
        /// Notifies new Weissenberg number
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Weissenbergnumber")) {
                m_Weissenberg = (double)cs.UserDefinedValues["Weissenbergnumber"];
            }
        }

        /// <summary>
        /// Nonlinear Form in trial variable (U, resp. gradient of U) - 
        /// using <see cref="EdgeFormDifferentiator"/> 
        /// and <see cref="VolumeFormDifferentiator"/>.
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if (SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");

            return new IEquationComponent[] {
                new EdgeFormDifferentiator(this, SpatialDimension),
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }
    }
}
