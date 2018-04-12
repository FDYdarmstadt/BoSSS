using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;

namespace BoSSS.Foundation.XDG {
    class NonlinearLevelSetFormVectorizer :
        ILevelSetForm_V,
        ILevelSetForm_GradV //
    {
        /// <summary>
        /// ctor.
        /// </summary>
        public NonlinearLevelSetFormVectorizer(LevelSetTracker lsTrk, ILevelSetForm _OrgComponent) {
            this.m_LsTrk = lsTrk;
            this.ArgumentOrdering = _OrgComponent.ArgumentOrdering.ToArray();
            this.ParameterOrdering = _OrgComponent.ParameterOrdering != null ? _OrgComponent.ParameterOrdering.ToArray() : null;
            this.LevelSetIndex = _OrgComponent.LevelSetIndex;
            this.PositiveSpecies = _OrgComponent.PositiveSpecies;
            this.NegativeSpecies = _OrgComponent.NegativeSpecies;
            this.LevelSetTerms = _OrgComponent.LevelSetTerms;
            this.OrgComponent = _OrgComponent;
        }

        LevelSetTracker m_LsTrk;

        /// <summary>
        /// The original component that is beeing vectorized.
        /// </summary>
        ILevelSetForm OrgComponent;


        public int LevelSetIndex {
            get;
            private set;
        }

        public SpeciesId PositiveSpecies {
            get;
            private set;
        }

        public SpeciesId NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get;
            private set;
        }

        public IList<string> ArgumentOrdering {
            get;
            private set;
        }

        public IList<string> ParameterOrdering {
            get;
            private set;
        }

        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            throw new NotImplementedException();
        }

        public void LevelSetForm_GradV(LevSetIntParams inp, MultidimensionalArray Koeff_GradV) {
            throw new NotImplementedException();
        }

        public void LevelSetForm_V(LevSetIntParams inp, MultidimensionalArray Koeff_V) {
            throw new NotImplementedException();
        }
    }
}
