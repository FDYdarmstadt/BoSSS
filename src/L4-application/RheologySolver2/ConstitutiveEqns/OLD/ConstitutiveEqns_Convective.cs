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

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Fluxes for convective part of extra stress tensor in the constitutive equations.
    /// </summary>
    /// 
    public class ConstitutiveEqns_Convective : LinearFlux {

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        double m_Weissenberg; // Weissenberg number
        double m_alpha;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] StressFunction;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        public ConstitutiveEqns_Convective(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap, double Weissenberg, double alpha) {
            this.Component = Component;
            this.m_BcMap = _BcMap;
            this.m_Weissenberg = Weissenberg;
            this.m_alpha = 1; // alpha;

            StressFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 3];

            StressFunction.SetColumn(m_BcMap.bndFunction[VariableNames.StressXX], 0);
            StressFunction.SetColumn(m_BcMap.bndFunction[VariableNames.StressXY], 1);
            StressFunction.SetColumn(m_BcMap.bndFunction[VariableNames.StressYY], 2);

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];
            for (int d = 0; d < 2; d++)
                velFunction.SetColumn(m_BcMap.bndFunction[VariableNames.Velocity_d(d)], d);
        }

        public override IList<string> ArgumentOrdering
        {
            get
            {
                switch (Component) {
                    case 0:
                        return new string[] { VariableNames.StressXX };
                    case 1:
                        return new string[] { VariableNames.StressXY };
                    case 2:
                        return new string[] { VariableNames.StressYY };
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return VariableNames.Velocity0Vector(2);
            }
        }

        // calculating the fluxes
        protected override void Flux(ref CommonParamsVol inp, double[] T, double[] output) {
            int D = output.Length;
            Array.Clear(output, 0, D);

            double u = inp.Parameters[0];
            double v = inp.Parameters[1];

            output[0] = m_Weissenberg * u * T[0];
            output[1] = m_Weissenberg * v * T[0];

        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Tin, double[] Tout) {
            double res = 0;

            double u = 0.5 * (inp.Parameters_IN[0] + inp.Parameters_OUT[0]);
            double v = 0.5 * (inp.Parameters_IN[1] + inp.Parameters_OUT[1]);

            Vector n = new Vector(inp.Normale[0], inp.Normale[1]);
            Vector velocityVector = new Vector(u, v);

            if (velocityVector * n > 0)
            {

                //Outflow
                res += u * (m_alpha * Tin[0] + (1 - m_alpha) * Tout[0]) * inp.Normale[0];
                res += v * (m_alpha * Tin[0] + (1 - m_alpha) * Tout[0]) * inp.Normale[1];

            }
            else
            {

                //Inflow
                res += u * (m_alpha * Tout[0] + (1 - m_alpha) * Tin[0]) * inp.Normale[0];
                res += v * (m_alpha * Tout[0] + (1 - m_alpha) * Tin[0]) * inp.Normale[1];
            }

            return m_Weissenberg * res;

        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Tin) {
            double res = 0;
            IncompressibleBcType edgType = m_BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType)
            {
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:

                    // Atmospheric outlet/pressure outflow: hom. Neumann
                    double u = inp.Parameters_IN[0];
                    double v = inp.Parameters_IN[1];
                    res += u * Tin[0] * inp.Normale[0];
                    res += v * Tin[0] * inp.Normale[1];
                    res *= m_Weissenberg;
                    break;

                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet:

                    // Setup params
                    // ============
                    Foundation.CommonParams inp2;
                    inp2.GridDat = inp.GridDat;
                    inp2.Normale = inp.Normale;
                    inp2.iEdge = inp.iEdge;
                    inp2.Parameters_IN = inp.Parameters_IN;
                    inp2.X = inp.X;
                    inp2.time = inp.time;

                    // Specify Parameters_OUT
                    // ======================
                    inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                    // Outer values for Velocity and VelocityMean
                    Debug.Assert(inp.Normale.Length == 2);
                    for (int j = 0; j < 2; j++)
                    {

                        inp2.Parameters_OUT[j] = velFunction[inp.EdgeTag, j](inp.X, inp.time);

                    }

                    //Dirichlet value for Stresses
                    //============================

                    switch (Component)
                    {
                        case 0:
                            double StressXX = StressFunction[inp.EdgeTag, 0](inp.X, inp.time);
                            res += InnerEdgeFlux(ref inp2, Tin, new double[] { StressXX });

                            break;

                        case 1:
                            double StressXY = StressFunction[inp.EdgeTag, 1](inp.X, inp.time);
                            res += InnerEdgeFlux(ref inp2, Tin, new double[] { StressXY });
   
                            break;

                        case 2:
                            double StressYY = StressFunction[inp.EdgeTag, 2](inp.X, inp.time);
                            res += InnerEdgeFlux(ref inp2, Tin, new double[] { StressYY });

                            break;
                        default:
                            throw new NotImplementedException();
                    }
                    break;

                default:
                    throw new NotImplementedException("unsupported/unknown b.c. - missing implementation;");
            }

            return res;
        }
    }
}
