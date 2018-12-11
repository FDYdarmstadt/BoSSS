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
using ilPSP;

namespace BoSSS.Application.Rheology
{
    class ConstitutiveEqns_CellWiseForm : IVolumeForm, IEdgeForm, IEquationComponentCoefficient
    {
        int ComponentRow;
        int ComponentCol;
        BoundaryCondMap<IncompressibleBcType> m_BcMap;
        double m_Weissenberg; // Weissenberg number
        double m_alpha; // upwind-paramter

        /// <summary>
        /// Mapping from edge tags to boundary values.
        /// - 1st index: edge tag;
        /// - 2nd index: spatial direction, row
        /// - 3rd index: spatial direction, column
        /// </summary>
        protected Func<double[], double, double>[,,] StressFunction;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        public ConstitutiveEqns_CellWiseForm(int _ComponentRow, int _ComponentCol, BoundaryCondMap<IncompressibleBcType> _BcMap, double Weissenberg, double alpha = 1.0)
        {
            ComponentRow = _ComponentRow;
            ComponentCol = _ComponentCol;
            this.m_BcMap = _BcMap;
            this.m_Weissenberg = Weissenberg;
            this.m_alpha = alpha;

            StressFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2,2];


            var stressXXfuncS = m_BcMap.bndFunction[VariableNames.StressXX];
            var stressXYfuncS = m_BcMap.bndFunction[VariableNames.StressXY];
            var stressYYfuncS = m_BcMap.bndFunction[VariableNames.StressYY];

            for (int et = 0; et < GridCommons.FIRST_PERIODIC_BC_TAG; et++)
            {
                StressFunction[et, 0, 0] = stressXXfuncS[et];
                StressFunction[et, 1, 0] = stressXYfuncS[et];
                StressFunction[et, 0, 1] = stressXYfuncS[et];
                StressFunction[et, 1, 1] = stressYYfuncS[et];
            }

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, 2];
            for (int d = 0; d < 2; d++)
                velFunction.SetColumn(m_BcMap.bndFunction[VariableNames.Velocity_d(d)], d);
        }


        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV; }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.StressXX , VariableNames.StressXY, VariableNames.StressYY };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return VariableNames.Velocity0Vector(2);
            }
        }

       

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            double[] Normale = inp.Normale;

            double[,] T__in = new double[2, 2];
            double[,] T_out = new double[2, 2];
            double[,] S__in = new double[2, 2];

            T__in[0, 0] = _uA[0]; // stress_XX
            T__in[1, 0] = _uA[1]; // stress_XY
            T__in[0, 1] = _uA[1]; // stress_XY
            T__in[1, 1] = _uA[2]; // stress_YY

            T_out[0, 0] = StressFunction[inp.EdgeTag,0,0](inp.X, inp.time); // stress_XX
            T_out[1, 0] = StressFunction[inp.EdgeTag, 1, 0](inp.X, inp.time); // stress_XY
            T_out[0, 1] = StressFunction[inp.EdgeTag, 0, 1](inp.X, inp.time); // stress_XY
            T_out[1, 1] = StressFunction[inp.EdgeTag, 1, 1](inp.X, inp.time); // stress_YY

            S__in[ComponentRow, ComponentCol] = _vA;

            double flxIn = MyCellBoundaryForm(Normale, T__in, T_out, S__in, inp.Parameters_IN, inp.Parameters_IN, false, 0);

            return flxIn;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
        {
            double[] Normale = inp.Normale.CloneAs();

            //double[,] T__in, double[,] T_out, double[,] S, double[] U0__in, double[] U0_out

            double[,] T__in = new double[2, 2];
            double[,] T_out = new double[2, 2];
            double[,] S__in = new double[2, 2];
            double[,] S_out = new double[2, 2];

            T__in[0, 0] = _uIN[0]; // stress_XX
            T__in[1, 0] = _uIN[1]; // stress_XY
            T__in[0, 1] = _uIN[1]; // stress_XY
            T__in[1, 1] = _uIN[2]; // stress_YY

            T_out[0, 0] = _uOUT[0]; // stress_XX
            T_out[1, 0] = _uOUT[1]; // stress_XY
            T_out[0, 1] = _uOUT[1]; // stress_XY
            T_out[1, 1] = _uOUT[2]; // stress_YY

            S__in[ComponentRow, ComponentCol] = _vIN;
            S_out[ComponentRow, ComponentCol] = _vOUT;

            double flxIn = MyCellBoundaryForm(Normale, T__in, T_out, S__in, inp.Parameters_IN, inp.Parameters_OUT, false, 0);
            Normale.ScaleV(-1.0);
            double flxOt = MyCellBoundaryForm(Normale, T_out, T__in, S_out, inp.Parameters_OUT, inp.Parameters_IN, false, 0);

            return flxIn + flxOt;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double[,] T = new double[2,2];
            double[,,] GradT = new double[2,2,2];
            double[,] S = new double[2,2];

            T[0, 0] = U[0]; // stress_XX
            T[1, 0] = U[1]; // stress_XY
            T[0, 1] = U[1]; // stress_XY
            T[1, 1] = U[2]; // stress_YY

            for(int d = 0; d < 2; d++)
            {
                GradT[0, 0, d] = GradU[0, d]; // stress_XX
                GradT[1, 0, d] = GradU[1, d]; // stress_XY
                GradT[0, 1, d] = GradU[1, d]; // stress_XY
                GradT[1, 1, d] = GradU[2, d]; // stress_YY
            }

            S[ComponentRow, ComponentCol] = V;
            Debug.Assert(GradV.L2NormPow2() == 0.0);

            return MyVolumeForm(T, GradT, S, cpv.Parameters);
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="T">
        /// Stress tensor (trial function)
        /// - 1st index: row index
        /// - 2nd index: column index
        /// </param>
        /// <param name="GradT">
        /// Gradient of stress tensor (trial function)
        /// - 1st index: row index
        /// - 2nd index: column index
        /// - 3rd index: derivative direction/gradient component
        /// </param>
        /// <param name="S">
        /// Test function tensor
        /// - 1st index: row index
        /// - 2nd index: column index
        /// </param>
        /// <param name="U0">
        /// Velocity
        /// </param>
        /// <returns></returns>

        double MyVolumeForm(double[,] T, double[,,] GradT, double[,] S, double[] U0)
        {
            double res = 0;

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double conv = 0;
                    for (int d = 0; d < 2; d++)
                    {
                        conv += U0[d] * GradT[i, j, d];
                    }

                    res += conv * S[i, j];
                }
            }
            
            return m_Weissenberg * res;
        }


        

        double MyCellBoundaryForm(double[] Normale, double[,] T__in, double[,] T_out, double[,] S, double[] U0__in, double[] U0_out, bool DomainBnd, byte EdgeTag)
        {
            double res = 0;

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double n_u = 0;
                    for (int d = 0; d < 2; d++)
                    {
                        n_u += Normale[d] * 0.5 * (U0__in[d] + U0_out[d]);
                    }

                    double factor;
                    if (n_u < 0)
                        factor = this.m_alpha;
                    else
                        factor = 1.0 - this.m_alpha;

                    res += factor * n_u * (T_out[i, j] - T__in[i, j]) * S[i, j];
                }
            }
            
            return m_Weissenberg * res;
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (cs.UserDefinedValues.Keys.Contains("Weissenbergnumber")) {
                m_Weissenberg = (double)cs.UserDefinedValues["Weissenbergnumber"];
                //Console.WriteLine("Weissenbergnumber = {0}", m_Weissenberg);
            }
        }
    }
}
