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

namespace BoSSS.Solution.RheologyCommon {
    /// <summary>
    /// Objective part of constitutive equations for single-phase flow.
    /// </summary>
    public class ConstitutiveEqns_Objective : IVolumeForm, IEdgeForm, IEquationComponent, IEquationComponentCoefficient, ISupportsJacobianComponent {

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        protected double m_Weissenberg;
        double m_StressPenalty;
        bool m_useFDJacobian;

        /// <summary>
        /// Initialize objective
        /// </summary>
        public ConstitutiveEqns_Objective(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap, double Weissenberg, double Penalty, bool UseFDJacobian) {
            this.Component = Component;
            this.m_BcMap = _BcMap;
            this.m_Weissenberg = Weissenberg;
            this.m_StressPenalty = Penalty;
            this.m_useFDJacobian = UseFDJacobian;

        }

        /// <summary>
        /// Choosing the required terms for volume form (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        /// </summary>
        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.UxV; }
        }

        /// <summary>
        /// Choosing the required terms for boundary edge form (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV; }
        }

        /// <summary>
        /// Choosing the required terms for inner edge form (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        /// </summary>
        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV; }
        }

        /// <summary>
        /// Ordering of the dependencies
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                string[] stresses;
                switch (Component) {
                    case 0:
                    stresses = new string[] { VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressXX, VariableNames.StressXY };
                    break;

                    case 1:
                    stresses = new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.StressXX, VariableNames.StressXY };
                    break;

                    case 2:
                    stresses = new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.StressXY, VariableNames.StressYY };
                    break;

                    default:
                    throw new NotImplementedException();
                }

                string[] Vels = VariableNames.VelocityVector(2);

                return stresses.Cat(Vels);
            }
        }

        /// <summary>
        /// Ordering of the parameters
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                if (m_useFDJacobian) {
                    switch (Component) {
                        case 0:
                            return new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY };
                        case 1:
                            return new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                        case 2:
                            return new string[] { VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                        default:
                            throw new NotImplementedException();
                    }
                } else {
                    return null;
                }
            }
        }

        /// <summary>
        /// Extraction of velocity gradients from arguments
        /// </summary>
        void GetVelocityGrad(out double Grad1, out double Grad2, out double Grad3, out double Grad4, double[,] U) {
            int offset = 4; // offset into trial var array
            switch (Component) {
                case 0:
                //VelGrads = new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY };
                Grad1 = U[offset + 0, 0];
                Grad2 = U[offset + 0, 1];
                Grad3 = U[offset + 0, 0];
                Grad4 = U[offset + 0, 1];
                break;

                case 1:
                //VelGrads = new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                Grad1 = U[offset + 0, 0];
                Grad2 = U[offset + 0, 1];
                Grad3 = U[offset + 1, 0];
                Grad4 = U[offset + 1, 1];
                break;

                case 2:
                //VelGrads = new string[] { VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                Grad1 = U[offset + 1, 0];
                Grad2 = U[offset + 1, 1];
                Grad3 = U[offset + 1, 0];
                Grad4 = U[offset + 1, 1];
                break;

                default:
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Calculating the integral of the volume part
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] T, double[,] Grad_T, double V, double[] GradV) {

            double Grad1;
            double Grad2;
            double Grad3;
            double Grad4;

            if (m_useFDJacobian) {

                Grad1 = cpv.Parameters[0];
                Grad2 = cpv.Parameters[1];
                Grad3 = cpv.Parameters[2];
                Grad4 = cpv.Parameters[3];

            } else {

                GetVelocityGrad(out Grad1, out Grad2, out Grad3, out Grad4, Grad_T);
            }

            double res = 0.0;
            res = ((Grad1 * T[0] + Grad2 * T[1]) + (Grad3 * T[2] + Grad4 * T[3]));

            return -m_Weissenberg * res * V;
        }

        /// <summary>
        /// Calculating the integral of the inner edge part
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] Tin, double[] Tout, double[,] Grad_Tin, double[,] Grad_Tout,
            double Vin, double Vout, double[] Grad_Vin, double[] Grad_Vout) {

            //double res = 0.0;
            //res += (Tin[0] - Tout[0]) + (Tin[1] - Tout[1]) + (Tin[2] - Tout[2]) + (Tin[3] - Tout[3]);
            //return (-m_Weissenberg) * res * (Vin - Vout);
            return 0.0;
        }

        /// <summary>
        /// Calculating the integral of the boundary edge part
        /// </summary>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Tin, double[,] Grad_Tin, double Vin, double[] Grad_Vin) {
            return 0.0;
        }

        /// <summary>
        /// update the coefficient such as the current Weissenberg number
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Weissenbergnumber")) {
                m_Weissenberg = (double)cs.UserDefinedValues["Weissenbergnumber"];
            }
        }

        /// <summary>
        /// Used by <see cref="GetJacobianComponents"/>:
        /// Since the edge forms of <see cref="ConstitutiveEqns_Objective"/> are linear (they are only penalties)
        /// one can use the implementation form the original object.
        /// Remark: the volume terms, on the other hand, are non-linear and must be differentiated.
        /// </summary>
        class OwnerCaller : IEdgeForm, IEquationComponentCoefficient {
            public OwnerCaller(ConstitutiveEqns_Objective __owner) {
                owner = __owner;
            }
            ConstitutiveEqns_Objective owner;

            public TermActivationFlags BoundaryEdgeTerms => owner.BoundaryEdgeTerms;

            public TermActivationFlags InnerEdgeTerms => owner.InnerEdgeTerms;

            public IList<string> ArgumentOrdering => owner.ArgumentOrdering.GetSubVector(0, 4); // for the boundary, the first arguments are sufficient.

            public IList<string> ParameterOrdering => null;

            public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
                return owner.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
            }

            public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
                return owner.InnerEdgeForm(ref inp, _uIN, _uOUT, _Grad_uIN, _Grad_uOUT, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            }

            public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
                owner.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if (SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");


            return new IEquationComponent[] {
                new VolumeFormDifferentiator(this, SpatialDimension),
                new OwnerCaller(this)
            };
        }

    }
}
