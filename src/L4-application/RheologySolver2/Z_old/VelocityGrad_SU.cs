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
using ilPSP.LinSolvers;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid.Classic;
using System.Diagnostics;

namespace BoSSS.Application.Rheology
{

    /// <summary>
    /// Fluxes for convective part of extra stress tensor in the constitutive equations.
    /// </summary>
    /// 
    public class VelocityGrad_SU : IVolumeForm, IEquationComponent
    {

        int Component;           // equation index (0: xx, 1: xy, 2: yx, 3:yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        public VelocityGrad_SU(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap)
        {
            this.Component = Component;
            this.m_BcMap = _BcMap;


            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];
            for (int d = 0; d < 2; d++)
                velFunction.SetColumn(m_BcMap.bndFunction[VariableNames.Velocity_d(d)], d);
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV; }
        }

        public IList<string> ArgumentOrdering {
            get {
                switch (Component)
                {
                    case 0:
                        return new string[] { VariableNames.VelocityX };
                    case 1:
                        return new string[] { VariableNames.VelocityX };
                    case 2:
                        return new string[] { VariableNames.VelocityY };
                    case 3:
                        return new string[] { VariableNames.VelocityY };

                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public IList<string> ParameterOrdering { get; }

        // calculating the fluxes
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            switch (Component)
            {
                case 0:
                    return U[0] * GradV[0];
                case 1:
                    return U[0] * GradV[0];
                case 2:
                    return U[0] * GradV[0];
                case 3:
                    return U[0] * GradV[0];

                default:
                    throw new NotImplementedException();
            }
        }
        public double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] GradUin, double[,] GradUout,
            double Vin, double Vout, double[] GradVin, double[] GradVout)
        {
            double res = 0;

            double u = 0.5 * (inp.Parameters_IN[0] + inp.Parameters_OUT[0]);
            double v = 0.5 * (inp.Parameters_IN[1] + inp.Parameters_OUT[1]);

            Vector n = new Vector(inp.Normale[0], inp.Normale[1]);
            Vector velocityVector = new Vector(u, v);

            if (velocityVector * n > 0)
            {
                //Outflow
                res += Uout[0] * inp.Normale[0] * (Vin - Vout);
            }
            else
            {
                //Inflow
                res += Uin[0] * inp.Normale[0] * (Vin - Vout);
            }

            return res;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] GradUin, double Vin, double[] GradVin)
        {
            double res = 0;
            IncompressibleBcType edgType = m_BcMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgType)
            {
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet:

                    res += Uin[0] * inp.Normale[0] * Vin;
                    break;
            }
            return res;
        }
    }
}
