
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.XNSE_Solver.Legacy;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// base-class for post-processing modules in XNSE
    /// </summary>
    public abstract class XNSFEinSituPostProcessingModule : XNSEinSituPostProcessingModule {


        /// <summary>
        /// control object
        /// </summary>
        new protected XNSFE_Control Control {
            get {
                return (XNSFE_Control)(base.Control);
            }
        }


        /// <summary>
        /// current velocity solution
        /// </summary>
        new protected XDGField[] CurrentVel {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.CurrentVel;
                } else if(this.SolverMain is XNSFE newSolver) {
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
        new protected XDGField CurrentPressure {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.Pressure;
                } else if(this.SolverMain is XNSFE newSolver) {
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
        /// Cut-Cell quadrature order used for the flow solver
        /// </summary>
        new protected int m_HMForder {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.m_HMForder;
                } else if(this.SolverMain is XNSFE newSolver) {

                    return newSolver.QuadOrder();
                    
                } else {
                    throw new NotImplementedException();
                }
            }
        }


        /// <summary>
        /// Cut-Cell quadrature order used for the flow solver
        /// </summary>
        new protected XDGField KineticEnergy {
            get {
                if(this.SolverMain is XNSE_SolverMain oldSolver) {
                    
                    return oldSolver.KineticEnergy;
                } else if(this.SolverMain is XNSFE newSolver) {

                    var a =  newSolver.Timestepping.Parameters.First(t => t.Identification == BoSSS.Solution.NSECommon.VariableNames.KineticEnergy);
                    return (XDGField) a;
                } else {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// If used, the Fourier level set; otherwise null;
        /// </summary>
        new protected FourierLevSetBase Fourier_LevSet {
            get { 
                if(base.SolverMain is XNSE_SolverMain oldSolver) {
                    return oldSolver.Fourier_LevSet;
                } else if(base.SolverMain is XNSFE newSolver) {
                    if(newSolver.LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet is FourierLevelSet fls) {
                        return fls.Fourier_LevSet;
                    } else {
                        return null;
                    }
                } else {
                    throw new NotSupportedException();
                }
            }
        }

    }
}
