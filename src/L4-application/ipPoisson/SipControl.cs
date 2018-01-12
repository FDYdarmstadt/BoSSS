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
using BoSSS.Solution.Control;
using BoSSS.Foundation;

namespace BoSSS.Application.SipPoisson {
    
    /// <summary>
    /// Control object for the ipPoisson solver.
    /// </summary>
    [Serializable]
    public class SipControl : AppControl {

        /// <summary>
        /// Type of <see cref="SipPoissonMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(SipPoissonMain);
        }
        
        ///// <summary>
        ///// Function which determines which part of the domain boundary is of Dirichlet type (true)
        ///// and which part of Neumann type (false).
        ///// </summary>
        //public Func<CommonParamsBnd,bool> IsDirichlet;

        ///// <summary>
        ///// Dirichlet boundary value
        ///// </summary>
        //public Func<CommonParamsBnd,double> g_Diri;

        ///// <summary>
        ///// Neumann boundary value
        ///// </summary>
        //public Func<CommonParamsBnd, double> g_Neum;
        
        /// <summary>
        /// Multiplyer for the penalty parameter, should be around 1.0.
        /// </summary>
        public double penalty_poisson = 1.3;

        /// <summary>
        /// string identifying the solver variant
        /// </summary>
        public string solver_name = "direct";
        
        /// <summary>
        /// run the solver more than once, e.g. for more reliable timing-results.
        /// </summary>
        public int NoOfSolverRuns = 2;

        /// <summary>
        /// True, if an exact solution -- in order to determine the error -- is provides.
        /// </summary>
        public bool ExactSolution_provided = false;
    }
}
