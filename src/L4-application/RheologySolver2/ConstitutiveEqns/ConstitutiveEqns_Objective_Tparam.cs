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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;

namespace BoSSS.Application.Rheology {
    public class ConstitutiveEqns_Objective_Tparam : IVolumeForm, IEquationComponent, IEquationComponentCoefficient {

        /// <summary>
        /// Volume integral of objective part of constitutive equations.
        /// </summary>
        ///

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        double m_Weissenberg; // relaxation factor lambda_1
        double m_ObjectiveParam;


        public ConstitutiveEqns_Objective_Tparam(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap, double Weissenberg, double ObjectiveParam) {
            this.Component = Component;
            this.m_BcMap = _BcMap;
            this.m_Weissenberg = Weissenberg;
            this.m_ObjectiveParam = ObjectiveParam;
        }

        // Choosing the required terms (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        public TermActivationFlags VolTerms
        {
            get { return TermActivationFlags.AllOn; }
        }

        /*
        public TermActivationFlags BoundaryEdgeTerms
        {
            get { return TermActivationFlags.UxV; }
        }

        public TermActivationFlags InnerEdgeTerms
        {
            get { return TermActivationFlags.UxV; }
        }
        */

        // Ordering the dependencies
        public IList<string> ArgumentOrdering
        {
            get
            {
                switch (Component) {
                    case 0:
                        return new string[] { VariableNames.VelocityX, };
                    case 1:
                        return new string[] { VariableNames.VelocityX, VariableNames.VelocityY };
                    case 2:
                        return new string[] { VariableNames.VelocityY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public IList<string> ParameterOrdering
        {
            get {
                switch (Component) {
                    case 0:
                        return new string[] { VariableNames.StressXXP, VariableNames.StressXYP };
                    case 1:
                        return new string[] { VariableNames.StressXYP, VariableNames.StressYYP, VariableNames.StressXXP, };
                    case 2:
                        return new string[] { VariableNames.StressXYP, VariableNames.StressYYP};
                    default:
                        throw new NotImplementedException();
                }
            }
        }


        // Calculating the fluxes
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double res = 0;

            switch (Component)
            {
                case 0:
                    res = 2 * (GradU[0, 0] * cpv.Parameters[0] + GradU[0, 1] * cpv.Parameters[1]);
                    break;
                case 1:
                    res = (GradU[0, 0] + GradU[1, 1]) * cpv.Parameters[0] + GradU[0, 1] * cpv.Parameters[1] + GradU[1, 0] * cpv.Parameters[2];
                    break;
                case 2:
                    res = 2 * (GradU[0, 0] * cpv.Parameters[0] + GradU[0, 1] * cpv.Parameters[1]);
                    break;
                default:
                    throw new NotImplementedException();
            }

            return (1 - m_ObjectiveParam) * -m_Weissenberg * res * V;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB,
            double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return 0.0;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0.0;
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Weissenbergnumber"))
                m_Weissenberg = (double)cs.UserDefinedValues["Weissenbergnumber"];
        }
    }
}
