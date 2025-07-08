using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.LsTest {

    /// <summary>
    /// PDE-solver-control object which defines configuration options for nonlinear and linear equation solvers
    /// </summary>
    [Serializable]
    [DataContract]
    public class SolverWithLevelSetUpdaterTestControl : SolverWithLevelSetUpdaterControl {

        public SolverWithLevelSetUpdaterTestControl() {
            this.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
        }

        public override Type GetSolverType() {
            return typeof(SolverWithLevelSetUpdaterTestCenter);
        }

        private int _NoOfLevelSets = 1;
        /// <summary>
        /// No of level sets, standard 1
        /// </summary>
        [DataMember]
        virtual public int NoOfLevelSets {
            get {
                return _NoOfLevelSets;
            }
            set {
                _NoOfLevelSets = value;
                SetDGdegree(value);
            }
        }


        private int _DegreeOfLevelSets = 1;
        /// <summary>
        /// DG Degree
        /// </summary>
        [DataMember]
        virtual public int DegreeOfLevelSets {
            get {
                return _DegreeOfLevelSets;
            }
            set {
                _DegreeOfLevelSets = value;
                SetDGdegree(value);
            }
        }

        /// <summary>
        /// Setting time dependent (component-wise) advection velocities for each level set. (DG Background field, only the values at the respective interface should be relevant)
        /// </summary>
        /// <remarks>
        /// Note: using the setter is not recommended when working with the job management system,
        /// since these values specified here cannot be serialized.
        /// Instead, <see cref="AppControl.InitialValues"/> should be used.
        /// </remarks>
        public void SetAdvectionVelocity(int iLevSet, Func<double[], double, double>[] V) {
            int D = V.Length;
            for (int d = 0; d < D; d++)
                this.InitialValues_Evaluators_TimeDep["Var_" + VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(iLevSet), VariableNames.VelocityVector(D)[d])] = V[d];
        }

        public override void SetDGdegree(int p) {
            
            for (int n = 0; n < this.NoOfLevelSets; n++) {
                FieldOptions[VariableNames.LevelSetDGidx(n)] = new FieldOpts() {
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                };
                FieldOptions[VariableNames.LevelSetCGidx(n)] = new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                };

                // adjust the names, to defferentiate between dummy solver Var and Parameter for level set movement
                FieldOptions[VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(n), "Var_Velocity*")] = new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                }; 
            }
        }


    }


}

