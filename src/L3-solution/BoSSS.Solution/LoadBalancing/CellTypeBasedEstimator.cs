using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.LoadBalancing {


    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    abstract public class CellTypeBasedEstimator : ICellCostEstimator {


        public CellTypeBasedEstimator() {
            CellClassifier = new NoOfSpeciesClassifier() {
                ConsiderAlsoNearCells = false
            };
        }


        public ICellClassifier CellClassifier { 
            get; 
            set; 
        }


        public abstract object Clone();
        
        
        public abstract int[][] GetEstimatedCellCosts();
        
        public virtual void Init(IApplication app) {
            m_app = app;
        }

        [NonSerialized]
        protected IApplication m_app;

        public abstract void UpdateEstimates(IApplication app);
    }
}
