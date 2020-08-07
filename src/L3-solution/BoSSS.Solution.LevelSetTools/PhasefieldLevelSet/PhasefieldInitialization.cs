using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

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
        public void InitCH()
        {
            CreateFields();
            CreateEquationsAndSolvers(null);
            RelaxationStep();
        }

        
        /// <summary>
        /// Locally refine p-Order in Cutcells
        /// </summary>
        private void AdaptP()
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
            if (this.GridData is GridData)
            {
                hmin = ((GridData)GridData).Cells.h_minGlobal;
            }
            else if (this.GridData is AggregationGridData)
            {
                hmin = ((AggregationGridData)GridData).AncestorGrid.Cells.h_minGlobal;
            }
            else
            {
                throw new NotImplementedException();
            }
            
            _Cahn = dInterface * (1/4.164) * hmin;

            // mobility coefficient
            _Diff = 1.0;

            // 0.0 = bulk diffusion, 1.0 = surface diffusion
            _Lambda = 0.0;
        }
    }
}
