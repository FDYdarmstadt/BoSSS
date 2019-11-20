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
    public class ConstitutiveEqns_Objective_withGrad : IVolumeForm, IEquationComponent {

        /// <summary>
        /// Volume integral of objective part of constitutive equations.
        /// </summary>
        ///

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        double m_Weissenberg; // relaxation factor lambda_1
        double m_ObjectiveParam;


        public ConstitutiveEqns_Objective_withGrad(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap, double Weissenberg, double ObjectiveParam) {
            this.Component = Component;
            this.m_BcMap = _BcMap;
            this.m_Weissenberg = Weissenberg;
            this.m_ObjectiveParam = ObjectiveParam;
        }

        // Choosing the required terms (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        public TermActivationFlags VolTerms
        {
            get { return TermActivationFlags.V | TermActivationFlags.UxV; }
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
                        return new string[] { VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressXX, VariableNames.StressXY,
                                              VariableNames.VelocityXGradientX, VariableNames.VelocityXGradientY, VariableNames.VelocityXGradientX, VariableNames.VelocityXGradientY};
                    case 1:
                        return new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.StressXX, VariableNames.StressXY,
                                              VariableNames.VelocityXGradientX, VariableNames.VelocityXGradientY, VariableNames.VelocityYGradientX, VariableNames.VelocityYGradientY};
                    case 2:
                        return new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.StressXY, VariableNames.StressYY,
                                              VariableNames.VelocityYGradientX, VariableNames.VelocityYGradientY, VariableNames.VelocityYGradientX, VariableNames.VelocityYGradientY};
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
                        return new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY };
                    case 1:
                        return new string[] { VariableNames.VelocityX_GradientX, VariableNames.VelocityX_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                    case 2:
                        return new string[] { VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY, VariableNames.VelocityY_GradientX, VariableNames.VelocityY_GradientY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }


        // Calculating the fluxes
        public double VolumeForm(ref CommonParamsVol cpv, double[] T, double[,] Grad_T, double V, double[] GradV) {

            return m_ObjectiveParam * -m_Weissenberg * ((T[4] * T[0] + T[5] * T[1]) + (T[6] * T[2] + T[7] * T[3])) * V;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB,
            double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return 0.0;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0.0;
        }
        
    }
}
