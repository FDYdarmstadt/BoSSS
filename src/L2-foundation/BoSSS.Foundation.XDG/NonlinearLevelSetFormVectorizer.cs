using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;

namespace BoSSS.Foundation.XDG {
    class NonlinearLevelSetFormVectorizer :
        INonlinLevelSetForm_V,
        INonlinLevelSetForm_GradV //
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

        
    }
}
