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
using System.Linq;
using System.Collections.Generic;
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using MPI.Wrappers;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution.Tecplot;
using System.Globalization;
using System.IO;
using ilPSP.Utils;
using ilPSP.Connectors.Matlab;
using System.Text;
using ilPSP;
using BoSSS.Solution.AdvancedSolvers;
using NUnit.Framework;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using Code = BoSSS.Solution.Control.LinearSolverCode;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;

namespace BoSSS.Application.XdgPoisson3 {


    /// <summary>
    /// Poisson Solver with a discontinuity in diffusion (<see cref="XdgPoisson3Control.MU_A"/>, <see cref="XdgPoisson3Control.MU_B"/>) at 
    /// the Level Set.
    /// </summary>
    public class XdgPoisson3Main : BoSSS.Solution.Application<XdgPoisson3Control> {


        /// <summary>
        /// App entry point 
        /// </summary>
        static void Main(string[] args) {
            InitMPI();
            //BoSSS.Application.XdgPoisson3.Tests.IterativeSolverTest(Code.exp_gmres_levelpmg);
            //BoSSS.Application.XdgPoisson3.Tests.ParabolaTest(2, 0.6);
            //throw new Exception("remove me");
            BoSSS.Solution.Application<XdgPoisson3Control>._Main(args, false, delegate () {
                return new XdgPoisson3Main();
            });
        }

#pragma warning disable 649

        [LevelSetTracker("-:A +:B", 1)]
        LevelSetTracker LevelSetTracker;

        [InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        LevelSet Phi;

        [InstantiateFromControlFile("u", "u", IOListOption.ControlFileDetermined)]
        XDGField u;

        /// <summary>
        /// exact solution (if known)
        /// </summary>
        [InstantiateFromControlFile("uEx", "u", IOListOption.ControlFileDetermined)]
        XDGField uEx;

        /// <summary>
        /// error of numerical solution (if exact solution known)
        /// </summary>
        [InstantiateFromControlFile("uErr", "u", IOListOption.ControlFileDetermined)]
        XDGField uErr;

        [InstantiateFromControlFile("rhs", "u", IOListOption.ControlFileDetermined)]
        XDGField rhs;

        [InstantiateFromControlFile("residual", "u", IOListOption.ControlFileDetermined)]
        XDGField residual;

        [InstantiateFromControlFile(
            new string[] { "du_dx", "du_dy", "du_dz" },
            new string[] { "u", "u", "u" },
            true, true,
            IOListOption.ControlFileDetermined)]
        VectorField<XDGField> GradientU;
#pragma warning restore 649


        /*
        static void MyHandler(object sender, UnhandledExceptionEventArgs args) {
            Exception e = (Exception)args.ExceptionObject;
            Console.WriteLine("MyHandler caught : " + e.Message);
            Console.WriteLine("Runtime terminating: {0}", args.IsTerminating);
            System.Environment.Exit(-1234);
        }
        */

        protected override void SetInitial(double t) {
            //this will suppress exception prompts
            //Workaround to prevent distrubance while executing batchclient
            if (this.Control.SuppressExceptionPrompt) {
                AppDomain currentDomain = AppDomain.CurrentDomain;
                //currentDomain.UnhandledException += new UnhandledExceptionEventHandler(MyHandler);
            }

           
            base.SetInitial(t);
            this.LsTrk.UpdateTracker(t);
            base.SetInitial(t);
            this.LsTrk.UpdateTracker(t);

            this.MGColoring = new SinglePhaseField[base.MultigridSequence.Length];
            for (int iLevel = 0; iLevel < base.MultigridSequence.Length; iLevel++) {
                this.MGColoring[iLevel] = new SinglePhaseField(new Basis(this.GridData, 0), "MGColoring_level_" + iLevel);
                base.MultigridSequence[iLevel].ColorDGField(this.MGColoring[iLevel]);
            }


        }

        protected override void CreateFields() {
            base.CreateFields();
            base.LsTrk = LevelSetTracker;
            if (Control.CutCellQuadratureType != base.LsTrk.CutCellQuadratureType)
                throw new ApplicationException();
        }

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {

            if(L != null)
                throw new NotSupportedException();


            // create operator
            // ===============

            if(this.Control.SetDefaultDiriBndCnd || this.Control.xLaplaceBCs == null || this.Control.xLaplaceBCs.g_Diri == null || this.Control.xLaplaceBCs.IsDirichlet == null) {
                this.Control.xLaplaceBCs.g_Diri = ((CommonParamsBnd inp) => 0.0);
                this.Control.xLaplaceBCs.IsDirichlet = (inp => true);
            }

            double penalty_multiplyer = base.Control.penalty_multiplyer;

            int order = this.u.Basis.Degree * 2;

            double MU_A = this.Control.MU_A;
            double MU_B = this.Control.MU_B;

            Op = new XDifferentialOperatorMk2(1, 1, (A, B, C) => order, this.LsTrk.SpeciesNames, "u", "c1");
            Op.AgglomerationThreshold = this.Control.AgglomerationThreshold;
            var lap = new XLaplace_Bulk(penalty_multiplyer, "u", this.Control.xLaplaceBCs, 1.0, MU_A, MU_B, this.Control.ViscosityMode);
            Op.EquationComponents["c1"].Add(lap);      // Bulk form
            Op.EquationComponents["c1"].Add(new XLaplace_Interface( MU_A, MU_B, penalty_multiplyer, this.Control.ViscosityMode));   // coupling form
            Op.EquationComponents["c1"].Add(new RHSSource(this.rhs));
            Op.IsLinear = true;

            Op.FluxesAreNOTMultithreadSafe = true;
            Op.Commit();

        }
     


        XDifferentialOperatorMk2 Op;

  
        SinglePhaseField[] MGColoring;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            base.TerminationKey = true;
            dt = 1.0;

            bool succ = this.Op.Solve(this.u.Mapping, MgConfig:this.OpConfig,
                nsc: this.Control.NonLinearSolver, lsc: this.Control.LinearSolver,
                verbose: true,
                queryHandler: base.QueryHandler);

            if(!succ)
                throw new ArithmeticException("solver did not converge");

            if (this.Control.ExcactSolSupported) {
                this.uErr.Clear();
                this.uErr.Acc(+1.0, u);
                this.uErr.Acc(-1.0, uEx);

                double L2_ERR = this.uErr.L2Norm();

                base.QueryHandler.ValueQuery("L2_ERR", L2_ERR, true);

                this.LsTrk.UpdateTracker(0.0);
                int order = Math.Max(u.Basis.Degree, Math.Max(uErr.Basis.Degree, uEx.Basis.Degree)) * 2 + 1;
                //var scheme = new XQuadSchemeHelper(this.LsTrk, this.Control.HMFversion, this.LsTrk.SpeciesIdS.ToArray());
                var scheme = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;


                var A_scheme = scheme.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("A"));
                var B_scheme = scheme.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("B"));

                ICompositeQuadRule<QuadRule> A_rule = A_scheme.Compile(this.GridData, order);
                ICompositeQuadRule<QuadRule> B_rule = B_scheme.Compile(this.GridData, order);

                foreach (var _rule in new ICompositeQuadRule<QuadRule>[] { A_rule, B_rule }) {
                    foreach (var ckR in _rule) {
                        ckR.Rule.Weights.CheckForNanOrInf(true, true, true);
                        ckR.Rule.Nodes.CheckForNanOrInf(true, true, true);
                    }
                }

                double L2_ERR_HMF_A = this.u.GetSpeciesShadowField("A").LxError(this.Control.InitialValues_Evaluators["uEx#A"].Vectorize(), null, A_rule).Sqrt();
                double L2_ERR_HMF_B = this.u.GetSpeciesShadowField("B").LxError(this.Control.InitialValues_Evaluators["uEx#B"].Vectorize(), null, B_rule).Sqrt();
                double L2_ERR_HMF = (L2_ERR_HMF_A.Pow2() + L2_ERR_HMF_B.Pow2()).Sqrt();
                base.QueryHandler.ValueQuery("L2_ERR_HMF_A", L2_ERR_HMF_A, true);
                base.QueryHandler.ValueQuery("L2_ERR_HMF_B", L2_ERR_HMF_B, true);
                base.QueryHandler.ValueQuery("L2_ERR_HMF", L2_ERR_HMF, true);


                Console.WriteLine("Error norm (standard):       " + L2_ERR);
                Console.WriteLine("Error norm (HMF, Species A): " + L2_ERR_HMF_A);
                Console.WriteLine("Error norm (HMF, Species B): " + L2_ERR_HMF_B);
                Console.WriteLine("Error norm (HMF):            " + L2_ERR_HMF);
            }

#if TEST
            OperatorAnalysis();
#endif

            return dt;
        }

        /// <summary>
        /// Operator stability analysis
        /// </summary>
        override public IDictionary<string,double> OperatorAnalysis(OperatorAnalysisConfig config) {
            return this.Op.OperatorAnalysis(this.u.Mapping, config, this.OpConfig); 
        }
                

        MultigridOperator.ChangeOfBasisConfig[][] OpConfig {
            get {
                int p = this.u.Basis.Degree;
                int NoOfLevelWithDifferentConfig = 1;// (u.Basis.Degree - 1); // for p-MG
                var retOpConfig = NoOfLevelWithDifferentConfig.ForLoop(iLvl =>
                    new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] { 0 }, mode = this.Control.PrePreCond, DegreeS = new int[] { p - iLvl } }
                    }
                );
                return retOpConfig;
            }
        }

        protected void CustomItCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            MaxMlevel=Mgop.LevelIndex;
            //currentRes.SaveToTextFileDebug(String.Format("Res_{0}_proc",iterIndex));
            //currentSol.SaveToTextFileDebug(String.Format("Sol_{0}_proc",iterIndex));
            //Console.WriteLine("Callback executed {0} times",iterIndex);
        }

        private int m_maxMlevel;

        public int MaxMlevel {
            get{
                return m_maxMlevel;
            }
            set{
                if (value > m_maxMlevel)
                    m_maxMlevel = value;
            }
        }

     
        private static void MgConsistencyTestRec(MultigridOperator mgOp, double[] mgSolVec, double[] mgRhsVec) {

            double[] mgResidual = mgRhsVec.CloneAs();
            mgOp.OperatorMatrix.SpMV(1.0, mgSolVec, -1.0, mgResidual);
            double scale = 1.0 / (mgRhsVec.L2Norm() + mgSolVec.L2Norm());

            long DOFs = mgOp.Mapping.TotalLength;
            Debug.Assert(DOFs == mgOp.OperatorMatrix.NoOfRows);
            Debug.Assert(DOFs == mgOp.OperatorMatrix.NoOfCols);


            double mgResidual_l2Norm = mgResidual.L2Norm();
            Console.WriteLine("Multigrid Residual norm, level {0}: {1} ({2} DOFs)", mgOp.LevelIndex, mgResidual_l2Norm, DOFs);
            Assert.LessOrEqual(scale * mgResidual_l2Norm, 1.0e-8);

            if (mgOp.CoarserLevel != null) {
                double[] mgCoarseSolVec = new double[mgOp.CoarserLevel.Mapping.LocalLength];
                double[] mgCoarseRhsVec = new double[mgOp.CoarserLevel.Mapping.LocalLength];

                mgOp.CoarserLevel.Restrict(mgSolVec, mgCoarseSolVec);
                mgOp.CoarserLevel.Restrict(mgRhsVec, mgCoarseRhsVec);

                MgConsistencyTestRec(mgOp.CoarserLevel, mgCoarseSolVec, mgCoarseRhsVec);
            }
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int susamp) {
            var Fields = ArrayTools.Cat<DGField>(this.GradientU, this.Phi, this.u, this.rhs, this.residual, this.uEx, this.uErr);
            if (this.MGColoring != null && this.MGColoring.Length > 0)
                Fields = ArrayTools.Cat<DGField>(Fields, this.MGColoring);
            Tecplot.PlotFields(Fields, "XPoisson" + timestepNo.ToString(), physTime, susamp);
            Tecplot.PlotFields(Fields, "grid" + timestepNo.ToString(), physTime, 0);
        }

    }

     class RHSSource : IVolumeForm, IParameterHandling {

        public RHSSource(DGField rhsSourceField) {
            m_rhsSourceField = rhsSourceField;
        }


        DGField m_rhsSourceField;

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new string[0];

        public IList<string> ParameterOrdering => new[] { "RHSsource" };

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            return new[] { m_rhsSourceField };
        }

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            if(!object.ReferenceEquals(m_rhsSourceField,Parameters[0])) {
                Parameters[0].Clear();
                Parameters[0].Acc(1.0, m_rhsSourceField);
            }
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double rhsVal = cpv.Parameters[0];
            return -rhsVal * V;
        }
    }

}
