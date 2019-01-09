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
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Central difference scheme for divergence operator of the extra stress tensor.
    /// </summary>

    public class VelocityGradXX : IVolumeForm, IEquationComponent {

        int Component;           // spatial dimension of momentum equation
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        //double InverseReynolds;
        //double[] pen1;
        //double pen2;
        protected Func<double[], double, double>[,] VelFunction;
        protected Func<double[], double, double>[,] StressFunction;

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradUxV; }
        }

        public IList<string> ArgumentOrdering {
            get {
                switch (Component)
                {
                    case 0:
                        return new string[] { VariableNames.VelocityXGradientX, VariableNames.VelocityX };
                    case 1:
                        return new string[] { VariableNames.VelocityXGradientY, VariableNames.VelocityX };
                    case 2:
                        return new string[] { VariableNames.VelocityYGradientX, VariableNames.VelocityY };
                    case 3:
                        return new string[] { VariableNames.VelocityYGradientY, VariableNames.VelocityY };

                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public IList<string> ParameterOrdering { get; }



        public VelocityGradXX(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap) {
            this.Component = Component;
            this.m_BcMap = _BcMap;

            VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];

            VelFunction.SetColumn(m_BcMap.bndFunction[VariableNames.VelocityX], 0);
            VelFunction.SetColumn(m_BcMap.bndFunction[VariableNames.VelocityY], 1);

        }



        // Calculating the fluxes
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            switch (Component)
            {
                case 0:
                    return -U[0] * V;
                case 1:
                    return -U[0] * V;
                case 2:
                    return -U[0] * V;
                case 3:
                    return -U[0] * V;

                default:
                    throw new NotImplementedException();
            }
            //switch (Component)
            //{
            //    case 0:
            //        return (U[0] - GradU[1, 0]) * V;
            //    case 1:
            //        return (U[0] - GradU[1, 1]) * V;
            //    case 2:
            //        return (U[0] - GradU[1, 0]) * V;
            //    case 3:
            //        return (U[0] - GradU[1, 1]) * V;

            //    default:
            //        throw new NotImplementedException();
            //}
        }


        //public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB,
        //    double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB)
        //{
        //    return 0.0;
        //}

        //public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        //{
        //    return 0.0;
        //}
    }
}