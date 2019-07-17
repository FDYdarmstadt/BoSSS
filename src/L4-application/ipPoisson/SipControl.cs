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
using System.Runtime.Serialization;
using System.IO;

namespace BoSSS.Application.SipPoisson {
    

    /// <summary>
    /// Control object for the ipPoisson solver.
    /// </summary>
    [DataContract]
    [Serializable]
    public class SipControl : AppControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public SipControl() : base() {
            base.LinearSolver.NoOfMultigridLevels = 1;
            base.CompMode = _CompMode.Steady;
            base.NoOfTimesteps = 1;
            base.LinearSolver.verbose = true;
        }



        /// <summary>
        /// Type of <see cref="SipPoissonMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(SipPoissonMain);
        }

        /// <summary>
        /// Re-sets all <see cref="AppControl.FieldOptions"/>
        /// </summary>
        public override void SetDGdegree(int p) {
            if(p < 1)
                throw new ArgumentOutOfRangeException("Symmetric interior penalty requires a DG degree of at least 1.");
            base.FieldOptions.Clear();
            base.AddFieldOption("T", p);
            base.AddFieldOption("Tex", p + 2, FieldOpts.SaveToDBOpt.unspecified); // exact solution: degree plus 2 is enough precision for comparison
        }

        /// <summary>
        /// Settings verification
        /// </summary>
        public override void Verify() {
            base.Verify();
        }

        /// <summary>
        /// Multiplier for the penalty parameter, should be around 1.0.
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        public double penalty_poisson = 1.3;

        ///// <summary>
        ///// string identifying the solver variant
        ///// </summary>
        //[DataMember]
        //public SolverCodes solver_name = SolverCodes.classic_pardiso;

        ///// <summary>
        ///// If any blocking is used (Schwarz, block Jacobi), a target for the block size.
        ///// Tests show that the ideal block size may be around 10000, but this may depend on computer, DG polynomial order, etc.
        ///// </summary>
        //[DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(99.0)]
        //public int TargetBlockSize = 10000;
        
        /// <summary>
        /// run the solver more than once, e.g. for more reliable timing-results.
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.InclusiveLowerBound(1.0)]
        public int NoOfSolverRuns = 1;

        /// <summary>
        /// True, if an exact solution -- in order to determine the error -- is provides.
        /// </summary>
        [DataMember]
        public bool ExactSolution_provided = false;

        /// <summary>
        /// Suppresses exception prompt, which disturbs local batch run with MiniBatchprocessor.
        /// </summary>
        [DataMember]
        public bool SuppressExceptionPrompt = false;

        /// <summary>
        /// Outputpath for analysis data. Set path to enable analysis output, e.g. calculation of condition number, residual plots for each multigridlevel, etc.
        /// </summary>
        [DataMember]
        public string WriteMeSomeAnalyse = null;

    }
}
