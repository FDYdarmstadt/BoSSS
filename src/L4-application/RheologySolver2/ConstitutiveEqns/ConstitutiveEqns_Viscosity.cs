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
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Application.Rheology {
    /// <summary>
    /// Volume integral of viscosity part of constitutive equations.
    /// </summary>
    /// 
    public class ConstitutiveEqns_Viscosity : IVolumeForm, IEdgeForm {

        int Component;           // equation index (0: xx, 1: xy, 2: yy)
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        double m_ViscosityNonNewton; // polymeric viscosity
        protected Func<double[], double, double>[,] VelFunction;
        double[] pen1;

        public ConstitutiveEqns_Viscosity(int Component, BoundaryCondMap<IncompressibleBcType> _BcMap, double beta, double[] Penalty1) {
            this.Component = Component;
            this.m_BcMap = _BcMap;
            this.m_ViscosityNonNewton = 1.0 - beta;
            this.pen1 = Penalty1;

            VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];

            VelFunction.SetColumn(m_BcMap.bndFunction[VariableNames.VelocityX], 0);
            VelFunction.SetColumn(m_BcMap.bndFunction[VariableNames.VelocityY], 1);
        }

        static string[] allArg = new string[] { VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY };
        public IList<string> ArgumentOrdering {
            get {
                //return new string[] { allArg[Component] }; 
                switch (Component) {
                    case 0:
                    return new string[] { VariableNames.VelocityX, VariableNames.VelocityX };
                    case 1:
                    return new string[] { VariableNames.VelocityX, VariableNames.VelocityY };
                    case 2:
                    return new string[] { VariableNames.VelocityY, VariableNames.VelocityY };
                    default:
                    throw new NotImplementedException();
                }
            }
        }


        // Choosing the required terms (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.AllOn; }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.AllOn; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.AllOn; }
        }


        // Calculating the integral
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] Grad_U, double V, double[] GradV) {
            
            //-2 * (1-beta) * 0.5 * (d u_i / d x_j + d u_j / d x_i)


            double res;
            switch (Component) {
                case 0:
                    res = U[0] * GradV[0] + U[1] * GradV[0];
                break;
                case 1:
                    res = (U[0] * GradV[1] + U[1] * GradV[0]);
                    break;
                case 2:
                    res = U[0] * GradV[1] + U[1] * GradV[1];
                break;
                default:
                throw new NotImplementedException();
            }

            return 2 * m_ViscosityNonNewton * 0.5 * res;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] GradUin, double[,] GradUout, double Vin, double Vout, double[] GradVin, double[] GradVout) {

            double res = 0;

            switch (Component)
            {
                case 0:
                    res = 0.5 * ((Uin[0] + Uout[0]) * inp.Normale[0] + (Uin[1] + Uout[1]) * inp.Normale[0]) // central difference fo grad(u) and grad(u)^T
                            + pen1[0] * ((Uin[0] - Uout[0]) * inp.Normale[0]) + pen1[1] * ((Uin[0] - Uout[0]) * inp.Normale[1]) // beta Penalty for grad(u)
                            + pen1[0] * ((Uin[1] - Uout[1]) * inp.Normale[0]) + pen1[1] * ((Uin[1] - Uout[1]) * inp.Normale[1]); // beta penalty for grad(u)^T

                    break;
                case 1:
                    res = 0.5 * ((Uin[0] + Uout[0]) * inp.Normale[1] + (Uin[1] + Uout[1]) * inp.Normale[0])
                            + pen1[0] * ((Uin[0] - Uout[0]) * inp.Normale[0]) + pen1[1] * ((Uin[0] - Uout[0]) * inp.Normale[1])
                            + pen1[0] * ((Uin[1] - Uout[1]) * inp.Normale[0]) + pen1[1] * ((Uin[1] - Uout[1]) * inp.Normale[1]);

                    break;
                case 2:
                    res = 0.5 * ((Uin[0] + Uout[0]) * inp.Normale[1] + (Uin[1] + Uout[1]) * inp.Normale[1])
                            + pen1[0] * ((Uin[0] - Uout[0]) * inp.Normale[0]) + pen1[1] * ((Uin[0] - Uout[0]) * inp.Normale[1])
                            + pen1[0] * ((Uin[1] - Uout[1]) * inp.Normale[0]) + pen1[1] * ((Uin[1] - Uout[1]) * inp.Normale[1]);

                    break;
                default:
                    throw new NotImplementedException();
            }

            return -2 * m_ViscosityNonNewton * 0.5 * res * (Vin - Vout);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_uA, double Vin, double[] _Grad_vA) {
            double res = 0;
            IncompressibleBcType edgType = m_BcMap.EdgeTag2Type[inp.EdgeTag];

            double n1 = 0;
            double n2 = 0;

            double Vel1 = 0;
            double Vel2 = 0;

            double VelocityX = VelFunction[inp.EdgeTag, 0](inp.X, inp.time);
            double VelocityY = VelFunction[inp.EdgeTag, 1](inp.X, inp.time);

            switch (Component) {
                case 0:
                n1 = inp.Normale[0];
                n2 = inp.Normale[0];
                Vel1 = VelocityX;
                Vel2 = VelocityX;
                break;
                case 1:
                n1 = inp.Normale[1];
                n2 = inp.Normale[0];
                    Vel1 = VelocityX;
                    Vel2 = VelocityY;
                break;
                case 2:
                n1 = inp.Normale[1];
                n2 = inp.Normale[1];
                Vel1 = VelocityY;
                Vel2 = VelocityY;
                break;
                default:
                throw new NotImplementedException();
            }

            switch (edgType) {
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:

                    // Atmospheric outlet/pressure outflow: hom. Neumann
                    res += Uin[0] * n1 + Uin[1] * n2;
                    break;

                case IncompressibleBcType.FreeSlip:
                    if (Component == 1)
                    {
                        res += Uin[0] * n1 + Uin[1] * n2;
                    }

                    break;

                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:

                    // Dirichlet value for Velocity
                    // ============================
                    if (Component == 1)
                    {
                        res += Vel1 * n1 + Vel2 * n2;// + Vel1 * n2 + Vel2 * n1;
                    }
                    else {
                        res += Vel1 * n1 + Vel2 * n2;
                    }

                    break;

                default:
                throw new NotImplementedException("unsupported/unknown b.c. - missing implementation;");
            }

            return -2 * m_ViscosityNonNewton * 0.5 * res * Vin;
        }
    }
}