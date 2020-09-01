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

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Holds information to initialize the phasefield near its conjected equilibrium
    /// also functionalities for local p-refinement and tuning of phasefield coefficients
    /// </summary>
    partial class Phasefield : Application
    {

        /// <summary>
        /// Initialize Cahn Hilliard Level Set
        /// </summary>
        public void InitCH(double _dt)
        {
            CreateFields();
            CreateEquationsAndSolvers(null);
            RelaxationStep(_dt: _dt);
        }

        
        /// <summary>
        /// Locally refine p-Order in Cutcells
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid)
        {

            throw new NotImplementedException();
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
            _Diff = this.ModTyp == ModelType.modelA ? 1.0 : _Cahn.Pow2();// (_Cahn / 4.0).Pow2() / this.Viscosity;//_Cahn; //1.0 / this.Peclet;

            // 0.0 = pure bulk diffusion, 1.0 = pure surface diffusion, not implemented in model A
            _Lambda = 0.0;
        }
        
    }
}
