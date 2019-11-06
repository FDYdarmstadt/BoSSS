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

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Utility class which helps the user in creating the function matrices that are needed
    /// in an implementation of LinearDualValueFlux;
    /// It is only necessary to override <see cref="_DerivativeSource"/>,
    /// the function matrices and offsets (aka. intercept) are constructed from the user-defined
    /// functions by this class.
    /// </summary>
    abstract public class LinearDerivativeSource : IVolumeForm {

        
        /// <summary>
        /// implementation of <see cref="IVolumeForm.VolumeForm"/>
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this._DerivativeSource(cpv.Xglobal, cpv.Parameters, GradU)*V;
        }

        /*
        /// <summary>
        /// extracts the <paramref name="FunctionMatrix"/> by multiple calls to <see cref="_DerivativeSource"/>;
        /// </summary>
        /// <param name="x"></param>
        /// <param name="Parameters"></param>
        /// <param name="FunctionMatrix"></param>
        void DerivativeSource(double[] x, double[] Parameters, double[,] FunctionMatrix) {
            int L = FunctionMatrix.GetLength(0);
            int D = FunctionMatrix.GetLength(1);
            if (m_U == null) {
                m_U = new double[FunctionMatrix.GetLength(0), FunctionMatrix.GetLength(1)];
                int Lp = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;
                m_MyParams = new double[Lp];
            }

            if (m_MyParams.Length > 0) {
                Array.Copy(Parameters, m_MyParams, m_MyParams.Length);
            }


            for (int l = L - 1; l >= 0; l--) {
                for (int d = D - 1; d >= 0; d--) {
                    m_U[l, d] = 1;
                    FunctionMatrix[l, d] = _DerivativeSource(x, m_MyParams, m_U);
                    m_U[l, d] = 0;
                }
            }
        }
        */

        /// <summary>
        /// Implementors may override this method to implement the source term.
        /// </summary>
        /// <param name="x">spatial position</param>
        /// <param name="Parameters">values of parameter fields</param>
        /// <param name="GradientU">
        /// the gradient of all argument DG fields;<br/>
        /// 1st index: argument index, defined by <see cref="IEquationComponent.ArgumentOrdering"/><br/>
        /// 2nd index: spatial direction
        /// </param>
        /// <returns>
        /// source value
        /// </returns>
        abstract public double _DerivativeSource(double[] x, double[] Parameters, double[,] GradientU);
        

        /// <summary>
        /// to be defined by user;
        /// </summary>
        abstract public IList<string> ArgumentOrdering { get; }

        /// <summary>
        /// null by default
        /// </summary>
        virtual public IList<string> ParameterOrdering { 
            get { 
                return null; 
            } 
        }

        /// <summary>
        /// this kind of source depends on trial function gradient and test function value.
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxV;
            }
        }
    }
}
