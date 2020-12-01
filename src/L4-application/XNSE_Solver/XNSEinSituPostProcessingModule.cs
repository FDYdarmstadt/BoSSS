using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// base-class for post-processing modules in XNSE
    /// </summary>
    public abstract class XNSEinSituPostProcessingModule : InSituPostProcessingModule {


        /// <summary>
        /// control object
        /// </summary>
        new protected XNSE_Control Control {
            get {
                return (XNSE_Control)(base.Control);
            }
        }


        /// <summary>
        /// current velocity solution
        /// </summary>
        protected XDGField[] CurrentVel {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.CurrentVel;
                } else if(this.SolverMain is XNSE newSolver) {
                    int D = this.SolverMain.GridData.SpatialDimension;

                    var ret = newSolver.CurrentState.Fields.Take(D).Select(f => (XDGField)f).ToArray();
                    for(int d = 0; d < D; d++) {
                        if(ret[d].Identification != VariableNames.Velocity_d(d))
                            throw new ApplicationException("Unable to identify velocity fields.");
                    }

                    return ret;
                } else {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// the level-set which represents the fluid interface
        /// </summary>
        protected LevelSet LevSet {
            get {
                return (LevelSet)(this.LsTrk.LevelSets[0]);
            }
        }


    }
}
