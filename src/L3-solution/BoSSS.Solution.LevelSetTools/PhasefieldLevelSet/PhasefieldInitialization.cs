using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{   

    /// <summary>
    /// Holds information to initialize the phasefield near its conjected equilibrium
    /// also functionalities for local p-refinement and tuning of phasefield coefficients
    /// </summary>
    partial class Phasefield : DgApplicationWithSolver<PhasefieldControl>
    {
        // remember if this is a reinitialization
        private static double Cahn_Reinit;

        ResidualLogger ResLogger;

        /// <summary>
        /// Initialize Cahn Hilliard Level Set
        /// </summary>
        public void InitCH()
        {
            CreateFields();

            SetCHCoefficents();

            InitFromSignedDistance(this.Control.cahn);

            //if (Cahn_Reinit != 0.0)
            //    //RelaxationStep();
            //    ReInit(Cahn_Reinit, this.Control.cahn);

            CreateEquationsAndSolvers(null);            

            InitLogFile(new Guid());
            WriteLogLine(0, 0.0);

            Cahn_Reinit = this.Control.cahn;
        }

        protected override void InitSolver()
        {
            XdgTimestepping.XdgTimestepping solver = new XdgTimestepping.XdgTimestepping(
                SOperator,
                CurrentState.Fields,
                CurrentResidual.Fields,
                Control.TimeSteppingScheme,
                MultigridOperatorConfig,
                MultigridSequence,
                Control.LinearSolver, Control.NonLinearSolver, 
                queryHandler:this.QueryHandler);

            LsTrk = solver.LsTrk; // register the dummy tracker which the solver created internally for the DG case

            base.Timestepping = solver;
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {

            if(L == null) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                // Creation of time-integrator (initial, no balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++

                InitSolver();
                Timestepping.RegisterResidualLogger(new ResidualLogger(this.MPIRank, this.DatabaseDriver, new Guid()));

            } else {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // restore BDF time-stepper after grid redistribution (dynamic load balancing)
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                Timestepping.DataRestoreAfterBalancing(L, CurrentState.Fields, CurrentResidual.Fields, base.LsTrk, base.MultigridSequence, this.Operator);
            }
        }

        /// <summary>
        /// Locally refine p-Order in Cutcells
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {

            throw new NotImplementedException();
        }

        int reinit = 0;

        private void InitFromSignedDistance(double cahn) {
            Console.WriteLine($"Initializing Phasefield from signed distance field:\n" +
                $"  thickness:  {cahn}");

            // we assume the current phasefield is close to the equilibrium tangenshyperbolicus form
            // also initial Phasefield from XNSE Solver must be given in signed distance form
            GridData GridDat = (GridData)(phi.GridDat);
            SinglePhaseField phiNew = new SinglePhaseField(phi.Basis);

            // compute and project 
            // step one calculate distance field phiDist = 0.5 * log(Max(1+phi, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn_old
            // step two project the new phasefield phiNew = tanh(phiDist/(sqrt(2) * Cahn_new))
            // here done in one step, with default quadscheme
            // ===================
            phiNew.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Phi
                    // -----------------------------
                    DGLevSet.Evaluate(j0, Len, NS, result);

                    // compute the pointwise values of the new level set
                    // -----------------------------

                    result.ApplyAll(x => Math.Tanh(x / (Math.Sqrt(2) * cahn)));
                }
            );

            phi.Clear();
            phi.Acc(1.0, phiNew);

            reinit++;
            //PlotCurrentState(0.0, reinit);
        }

        // Reinitialize Phasefield with changed interface thickness
        private void ReInit(double cahn_old, double cahn)
        {
            Console.WriteLine($"Reprojecting Phasefield:\n" +
                $"  old thickness:  {cahn_old}\n" +
                $"  new thickness:  {cahn}");
            // we assume the current phasefield is close to the equilibrium tangenshyperbolicus form
            SinglePhaseField phiNew = new SinglePhaseField(phi.Basis);
            GridData GridDat = (GridData)(phi.GridDat);

            // compute and project 
            // step one calculate distance field phiDist = 0.5 * log(Max(1+phi, eps)/Max(1-phi, eps)) * sqrt(2) * Cahn_old
            // step two project the new phasefield phiNew = tanh(phiDist/(sqrt(2) * Cahn_new))
            // here done in one step, with default quadscheme
            // ===================
            phiNew.ProjectField(
                (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) 
                { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes

                    // evaluate Phi
                    // -----------------------------
                    phi.Evaluate(j0, Len, NS, result);

                    // compute the pointwise values of the new level set
                    // -----------------------------

                    result.ApplyAll(x => Math.Tanh(0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * (cahn_old / cahn)));
                }                
            );

            phi.Clear();
            phi.Acc(1.0, phiNew);

            CorrectionLevSet.Clear();
            CorrectionLevSet.Acc(1.0, phi);
            this.CorrectionLsTrk.UpdateTracker(0.0);

            // update DG LevelSet
            DGLevSet.Clear();
            DGLevSet.Acc(1.0, phi);

            reinit++;
            PlotCurrentState(0.0, reinit);
        }

        /// <summary>
        /// Set Coefficients in Cahn Hilliard based on some metric
        /// </summary>
        /// <param name="_Cahn"></param>
        /// <param name="_Diff"></param>
        /// <param name="_Lambda"></param>
        private void SetCHCoefficents()
        {
            if (!this.Control.FixedConstants)
            {
                // interface thickness = f(pOrder, D), WIP
                double dInterface = 2.0 * Math.Pow(2.0, this.GridData.SpatialDimension) / this.phi.Basis.Degree;

                // dInterface * 1/4.164 * hmin
                double hmin;

                // set the interface width based on the largest CutCell
                hmin = this.CorrectionLsTrk.Regions.GetCutCellSubGrid().h_maxSubGrd;

                this.Control.cahn = dInterface * (1.0 / 4.164) * hmin / Math.Sqrt(2);

                Console.WriteLine("Setting Interface thickness to {0}", this.Control.cahn);

                // mobility coefficient, could be inverse of Pe = Re*Sc
                // following Yue et. al. (2010), see also Manzanero (2019), for now fixed mobility
                this.Control.diff = this.Control.ModTyp == PhasefieldControl.ModelType.modelA ? 1e-4 : 1e-4;//this.Control.cahn.Pow2();// (_Cahn / 4.0).Pow2() / this.Viscosity;//_Cahn; //1.0 / this.Peclet;

                // 0.0 = pure bulk diffusion, 1.0 = pure surface diffusion, not implemented in model A
                this.Control.lambda = 0.0;
            }
        }

    }
}
