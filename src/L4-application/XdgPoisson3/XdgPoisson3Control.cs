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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XNSECommon;
using System.Runtime.Serialization;

namespace BoSSS.Application.XdgPoisson3 {


    [DataContract]
    [Serializable]
    public class XdgPoisson3Control : AppControl {

        /// <summary>
        /// Ctor.
        /// </summary>
        public XdgPoisson3Control() {
            base.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            base.LinearSolver.NoOfMultigridLevels = 1;
            base.LinearSolver.SolverCode = LinearSolverConfig.Code.classic_mumps; //public string solverName = "direct";
        }

        /// <summary>
        /// Type of <see cref="SipPoissonMain"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(XdgPoisson3Main);
        }

        /// <summary>
        /// Settings verification
        /// </summary>
        public override void Verify() {
            base.Verify();
        }

        /// <summary>
        /// DG degree of Level-Set is hard-coded to 2.
        /// </summary>
        /// <param name="p"></param>
        public override void SetDGdegree(int p) {
            FieldOptions.Clear();

            FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 2,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            FieldOptions.Add("u", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
        }
        /// <summary>
        /// This is a workaround to set Dirichlet-BndConditions, if executed from wroksheet (delegates are not serializable). True will set hom. Dirichlet-BndCnd.
        /// </summary>
        [DataMember]
        public bool SetDefaultDiriBndCnd = false;

        /// <summary>
        /// Ignore all costly unnecessary stuff ...
        /// </summary>
        [DataMember]
        public bool PerformanceModeON = false;

        [DataMember]
        public XLaplaceBCs xLaplaceBCs = new XLaplaceBCs();

        [DataMember]
        public XLaplace_Interface.Mode ViscosityMode = XLaplace_Interface.Mode.SIP;

        [DataMember]
        public double MU_A;


        [DataMember]
        public double MU_B;

        /// <summary>
        /// true if exact solution is supported
        /// </summary>
        [DataMember]
        public bool ExcactSolSupported = false;

        ///// <summary>
        ///// which solver to use
        ///// </summary>
        //public string solverName = "direct";

        /// <summary>
        /// preconditioner option
        /// </summary>
        [DataMember]
        public MultigridOperator.Mode PrePreCond = MultigridOperator.Mode.SymPart_DiagBlockEquilib;

        /// <summary>
        /// XDG cell agglomeration threshold
        /// </summary>
        [DataMember]
        public double AgglomerationThreshold = 0.1;

        [DataMember]
        public double penalty_multiplyer = 2.0;

        /// <summary>
        /// Steady-State (false) or transient solution (true)?
        /// </summary>
        [DataMember]
        public bool timeDependent = false;

        [DataMember]
        public int pOff = 2;
    }
}
