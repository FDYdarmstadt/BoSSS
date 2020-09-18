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

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Holds information to initialize the phasefield near its conjected equilibrium
    /// also functionalities for local p-refinement and tuning of phasefield coefficients
    /// </summary>
    partial class Phasefield : Application
    {
        // remember if this is a reinitialization
        private static double Cahn_Reinit;

        /// <summary>
        /// Initialize Cahn Hilliard Level Set
        /// </summary>
        public void InitCH(PhasefieldControl _Control)
        {

            if (_Control != null)
                Control = _Control;
            else
                Control = new PhasefieldControl();


            CreateFields();
            CreateEquationsAndSolvers(null);

            if (Cahn_Reinit != 0.0)
                //RelaxationStep();
                ReInit(Cahn_Reinit, Cahn);

            Cahn_Reinit = Cahn;
        }

        
        /// <summary>
        /// Locally refine p-Order in Cutcells
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid)
        {

            throw new NotImplementedException();
        }

        int reinit = 0;

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
        private void SetCHCoefficents(out double _Cahn, out double _Diff, out double _Lambda)
        {
            // interface thickness = f(pOrder, D), WIP
            double dInterface =  2.0 * Math.Pow(2.0, this.GridData.SpatialDimension) /this.phi.Basis.Degree;
            
            // dInterface * 1/4.164 * hmin
            double hmin;
            //if (this.GridData is GridData)
            //{
            //    hmin = ((GridData)GridData).Cells.h_minGlobal;
            //}
            //else if (this.GridData is AggregationGridData)
            //{
            //    hmin = ((AggregationGridData)GridData).AncestorGrid.Cells.h_minGlobal;
            //}
            //else
            //{
            //    throw new NotImplementedException();
            //}

            // set the interface width based on the largest CutCell
            hmin = this.LsTrk.Regions.GetCutCellSubGrid().h_maxSubGrd;            

            _Cahn = dInterface * (1.0/4.164) * hmin;

            Console.WriteLine("Setting Interface thickness to {0}", _Cahn);

            // mobility coefficient, for now inverse of Pe = Re*Sc
            // following Yue et. al. (2010), WARNING not correct
            _Diff = this.ModTyp == ModelType.modelA ? 1.0 : _Cahn;// (_Cahn / 4.0).Pow2() / this.Viscosity;//_Cahn; //1.0 / this.Peclet;

            // 0.0 = pure bulk diffusion, 1.0 = pure surface diffusion, not implemented in model A
            _Lambda = 0.0;
        }
        
    }
}
