
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
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

namespace BoSSS.Application.XNSFE_Solver {

    /// <summary>
    /// base-class for post-processing modules in XNSE
    /// </summary>
    public abstract class XNSFEinSituPostProcessingModule<T> : XNSEinSituPostProcessingModule<T> where T:XNSFE_Control, new() {

        /// <summary>
        /// reference to solver application class
        /// </summary>
        protected new XNSE<T> SolverMainOverride {
            get {
                return (XNSE<T>)base.SolverMain;
            }
        }

        /// <summary>
        /// control object
        /// </summary>
        new protected T Control {
            get {
                return (T)(base.Control);
            }
        }

        /// <summary>
        /// current temperature solution
        /// </summary>
        protected XDGField CurrentTemperature {
            get {
                int D = this.SolverMainOverride.GridData.SpatialDimension;

                var ret = this.SolverMainOverride.CurrentState.Fields.ElementAt(D) as XDGField;
                if (ret.Identification != VariableNames.Temperature)
                    throw new ApplicationException("Unable to identify temperature field.");

                return ret;
            }
        }

        /// <summary>
        /// current massflux parameter
        /// </summary>
        protected DGField CurrentMassFlux {
            get {
                int D = this.SolverMainOverride.GridData.SpatialDimension;
                IReadOnlyDictionary<string, DGField> parameters = this.SolverMainOverride.LsUpdater.Parameters;

                DGField ret = null;
                for (int i = 0; i < 3; ++i) {
                    if (parameters.TryGetValue(VariableNames.MassFluxExtension, out DGField velocityField)) {
                        ret = velocityField;
                    } else {
                        throw new ApplicationException("Unable to identify mass flux extension field.");
                    }
                }

                return ret;
            }
        }

        /// <summary>
        /// current ls velocity
        /// </summary>
        protected DGField[] CurrentVelocityLevelSet {
            get {
                int D = this.SolverMainOverride.GridData.SpatialDimension;
                IReadOnlyDictionary<string, DGField> parameters = this.SolverMainOverride.LsUpdater.Parameters;

                DGField[] ret = new DGField[D];
                for (int d = 0; d < D; ++d) {
                    if (parameters.TryGetValue(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCG, VariableNames.Velocity_d(d)), out DGField velocityField)) {
                        ret[d] = velocityField;
                    } else {
                        throw new ApplicationException("Unable to identify level set velocity y field.");
                    }
                }

                return ret;
            }
        }


    }
}
