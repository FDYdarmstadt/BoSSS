using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ZLSinSituPostProcessing
{
    /// <summary>
    /// base-class for post-processing modules in ZLS
    /// </summary>
    public abstract class ZLSinSituPostProcessingModule : InSituPostProcessingModule
    {


        /// <summary>
        /// control object
        /// </summary>
        new protected ZLS_Control Control
        {
            get
            {
                return (ZLS_Control)(base.Control);
            }
        }


        /// <summary>
        /// current velocity solution
        /// </summary>
        protected XDGField[] CurrentVel
        {
            get
            {
                if (this.SolverMain is ZLS solver)
                {
                    int D = this.SolverMain.GridData.SpatialDimension;

                    var ret = solver.CurrentState.Fields.Take(D).Select(f => (XDGField)f).ToArray();
                    for (int d = 0; d < D; d++)
                    {
                        if (ret[d].Identification != BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))
                            throw new ApplicationException("Unable to identify velocity fields.");
                    }

                    return ret;
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// current Displacement solution
        /// </summary>
        protected XDGField[] CurrentDisplacement
        {
            get
            {
                if (this.SolverMain is ZLS solver)
                {
                    int D = this.SolverMain.GridData.SpatialDimension;

                    var ret = solver.CurrentState.Fields.TakeLast(D).Select(f => (XDGField)f).ToArray();
                    for (int d = 0; d < D; d++)
                    {
                        if (ret[d].Identification != VariableNames.DisplacementVector(D)[d])
                            throw new ApplicationException("Unable to identify Displacement fields.");
                    }

                    return ret;
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// current pressure solution
        /// </summary>
        protected XDGField CurrentPressure
        {
            get
            {
                if (this.SolverMain is ZLS solver)
                {
                    int D = this.SolverMain.GridData.SpatialDimension;

                    var ret = solver.CurrentState.Fields.ElementAt(D) as XDGField;
                    if (ret.Identification != BoSSS.Solution.NSECommon.VariableNames.Pressure)
                        throw new ApplicationException("Unable to identify pressure field.");

                    return ret;
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
        }


        /// <summary>
        /// the level-set which represents the interfaces
        /// </summary>
        protected LevelSet LevSet
        {
            get
            {
                return (LevelSet)(this.LsTrk.LevelSets[0]);
            }
        }

        protected LevelSet LevSet1
        {
            get
            {
                return (LevelSet)(this.LsTrk.LevelSets[1]);
            }
        }


        /// <summary>
        /// Cut-Cell quadrature order used for the flow solver
        /// </summary>
        protected int m_HMForder
        {
            get
            {

                if (this.SolverMain is XNSE<ZLS_Control> solver)
                {

                    return solver.QuadOrder();

                }
                else
                {
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// If used, the Fourier level set; otherwise null;
        /// </summary>
        protected FourierLevSetBase Fourier_LevSet
        {
            get
            {
                if (base.SolverMain is ZLS newSolver)
                {
                    if (newSolver.LsUpdater.LevelSets[BoSSS.Solution.NSECommon.VariableNames.LevelSetCG].DGLevelSet is FourierLevelSet fls)
                    {
                        return fls.Fourier_LevSet;
                    }
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    throw new NotSupportedException();
                }
            }
        }

        /// <summary>
        /// les gird
        /// </summary>
        protected BoSSS.Foundation.Grid.Classic.GridData GridData
        {
            get
            {
                return (BoSSS.Foundation.Grid.Classic.GridData)(this.SolverMain.GridData);
            }
        }



    }
}
