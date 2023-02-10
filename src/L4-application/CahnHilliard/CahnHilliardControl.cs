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
using BoSSS.Solution.NSECommon;

namespace BoSSS.Application.CahnHilliard {


    /// <summary>
    /// Control object for the ipPoisson solver.
    /// </summary>
    [DataContract]
    [Serializable]
    public class CahnHilliardControl : AppControlSolver {

        /// <summary>
        /// Ctor.
        /// </summary>
        public CahnHilliardControl() : base() {
            base.TimesteppingMode = _TimesteppingMode.Transient;
            base.NoOfTimesteps = 1;
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
        }

        /// <summary>
        /// Type of <see cref="CahnHilliardMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(CahnHilliardMain);
        }

        /// <summary>
        /// Re-sets all <see cref="AppControl.FieldOptions"/>
        /// </summary>
        public override void SetDGdegree(int p) {
            if (p < 1)
                throw new ArgumentOutOfRangeException("Cahn Hilliard solver, due to use of symmetric interior penalty requires a DG degree of at least 1.");
            base.FieldOptions.Clear();
            base.AddFieldOption("c", p);
            base.AddFieldOption(VariableNames.Curvature, 0);
            base.AddFieldOption("Velocity*", p);
            base.AddFieldOption("cex", p + 2); // exact solution: degree times 2
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
        public double penalty_poisson = 2.6;

        /// <summary>
        /// Model Type of Phasefield equation, see Halperin (1977)
        /// </summary>
        public enum ModelType {
            /// <summary>
            /// Order Parameter is nonconserved
            /// </summary>
            modelA,

            /// <summary>
            /// Order Parameter is conserved
            /// </summary>
            modelB,

            /// <summary>
            /// Mass is conserved
            /// </summary>
            modelC
        }

        /// <summary>
        /// Set the <see cref="ModelType"/>
        /// </summary>
        [DataMember]
        public ModelType ModTyp = ModelType.modelB;


        ///// <summary>
        ///// According to Biben (2003), Correction to account for arclength diffusion
        ///// </summary>
        //[DataMember]
        //public bool CurvatureCorrection = false;


        /// <summary>
        /// Type of algebraic correction that is performed
        /// </summary>
        public enum Correction {
            /// <summary>
            /// Mass of a phase is conserved
            /// </summary>
            Mass,

            /// <summary>
            /// Total Concentration is conserved
            /// </summary>
            Concentration,

            /// <summary>
            /// No algebraic correction
            /// </summary>
            None
        }
        

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public Correction CorrectionType = Correction.None;

        ///// <summary>
        ///// Curvature as additional Equation Component or by direct evaluation and Parameter Update
        ///// </summary>
        //[DataMember]
        //public bool UseDirectCurvature = false;

        /// <summary>
        /// True, if an exact solution -- in order to determine the error -- is provides.
        /// </summary>
        [DataMember]
        public bool ExactSolution_provided = false;

        /// <summary>
        /// Include or Exclude Convective Term in Cahn-Hilliard Equation, not implemented
        /// </summary>
        [DataMember]
        public bool includeConvection = true;

        /// <summary>
        /// Include or Exclude Convective Term in Cahn-Hilliard Equation, not implemented
        /// </summary>
        [DataMember]
        public bool includeDiffusion = true;

        /// <summary>
        /// True if Jacobian is calculated by finite differences
        /// </summary>
        [DataMember]
        public bool UseFDJacobian = true;

        /// <summary>
        /// tru to write a log file
        /// </summary>
        [DataMember]
        public bool StoreBenchmarkQinFile = false;


        ///// <summary>
        ///// Some parameter of Cahn Hilliard equation
        ///// </summary>
        //[DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        //public double kappa = 1;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation, to adjust surface vs bulk diffusion
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.InclusiveLowerBound(0.0)]
        [BoSSS.Solution.Control.InclusiveUpperBound(1.0)]
        public double lambda = 0;

        ///// <summary>
        ///// Some parameter of Cahn Hilliard equation
        ///// </summary>
        //[DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        //public double epsilon = 1;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation
        /// diff = (kappa*lambda)/epsilon
        /// </summary>
        [DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        public double diff = 1;

        ///// <summary>
        ///// Some parameter of Cahn Hilliard equation
        ///// Peclet´s Number
        ///// </summary>
        //[DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        //public double peclet = 1;

        ///// <summary>
        ///// Some parameter of Cahn Hilliard equation
        ///// Capillary Number
        ///// </summary>
        //[DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        //public double capillary = 1;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation
        /// Cahn´s Number
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        public double cahn = 1;

        ///// <summary>
        ///// Some parameter of Cahn Hilliard equation
        ///// Reynold´s Number
        ///// </summary>
        //[DataMember]
        //[BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        //public double reynold = 1;

    }
}
