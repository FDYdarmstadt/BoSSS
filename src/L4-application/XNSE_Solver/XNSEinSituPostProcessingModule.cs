
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
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
        /// current velocity solution
        /// </summary>
        protected XDGField CurrentPressure {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.Pressure;
                } else if(this.SolverMain is XNSE newSolver) {
                    int D = this.SolverMain.GridData.SpatialDimension;

                    var ret = newSolver.CurrentState.Fields.ElementAt(D) as XDGField;
                    if(ret.Identification != VariableNames.Pressure)
                        throw new ApplicationException("Unable to identify pressure field.");
                    
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


        /// <summary>
        /// Cut-Cell quadrature order used for the flow solver
        /// </summary>
        protected int m_HMForder {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.m_HMForder;
                } else if(this.SolverMain is XNSE newSolver) {

                    return newSolver.QuadOrder();
                    
                } else {
                    throw new NotImplementedException();
                }
            }
        }


        /// <summary>
        /// Cut-Cell quadrature order used for the flow solver
        /// </summary>
        protected XDGField KineticEnergy {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.KineticEnergy;
                } else if(this.SolverMain is XNSE newSolver) {

                    throw new NotImplementedException("your turn, Lauritz");
                    
                } else {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// If used, the Fourier level set; otherwise null;
        /// </summary>
        protected FourierLevSetBase Fourier_LevSet {
            get { 
                if(base.SolverMain is XNSE_SolverMain oldSolver) {
                    return oldSolver.Fourier_LevSet;
                } else if(base.SolverMain is XNSE newSolver) {
                    throw new NotImplementedException("your turn, Lauritz");
                } else {
                    throw new NotSupportedException();
                }
            }
        }

        /// <summary>
        /// les gird
        /// </summary>
        protected BoSSS.Foundation.Grid.Classic.GridData GridData {
            get {
                return (Foundation.Grid.Classic.GridData)(this.SolverMain.GridData);
            }
        }



    }
}
