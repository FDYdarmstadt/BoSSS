using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StokesHelical_Ak {
    public partial class HelicalMain : Application<HelicalControl> {

        // Level-Set tracker
        [LevelSetTracker("-:A +:B", 1)]
        public LevelSetTracker MyLsTrk;

        /// <summary>
        /// the DG representation of the level set.
        /// This one is used for level-set evolution in time; it is in general discontinuous.
        /// </summary>
        [InstantiateFromControlFile("PhiDG", "PhiDG", IOListOption.ControlFileDetermined)]
        public ScalarFieldHistory<SinglePhaseField> DGLevSet;

        /// <summary>
        /// The  continuous level set field which defines the XDG space; 
        /// it is obtained from the projection of the discontinuous <see cref="DGLevSet"/> onto the 
        /// continuous element space.
        /// </summary>
        [InstantiateFromControlFile("Phi", "Phi", IOListOption.ControlFileDetermined)]
        public LevelSet LevSet;

        /// <summary>
        /// Compute Operator Matrix. The Magic happens here: <see cref="DelComputeOperatorMatrixMod"/>
        /// </summary>
        protected virtual void DelComputeOperatorMatrixMod(BlockMsrMatrix OpMatrix, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double phystime, bool OnlyAffine) {
            // compute operator
            if(OnlyAffine) {
                if(OpMatrix != null) {
                    throw new ArgumentException();
                }
            }
            // update global dirichlet value for boundary edges
            if(!this.Control.steady) {
                    Globals.DirichletValue_uR = (X => this.Control.BoundaryValues["Dirichlet"].Values["ur"].Evaluate(X, phystime));
                    Globals.DirichletValue_uXi = (X => this.Control.BoundaryValues["Dirichlet"].Values["uxi"].Evaluate(X, phystime));
                    Globals.DirichletValue_uEta = (X => this.Control.BoundaryValues["Dirichlet"].Values["ueta"].Evaluate(X, phystime));
                    Globals.DirichletValue_psi = (X => this.Control.BoundaryValues["Dirichlet"].Values["Pressure"].Evaluate(X, phystime));
            }


            //_____________________________________________________________________________
            // create matrix and affine vector:
            //if(OpMatrix != null || OnlyAffine) {
            var mtxBuilder = diffOp_implicit.GetMatrixBuilder(Mapping, null, Mapping);
            mtxBuilder.time = phystime;
            if(!OnlyAffine) {
                // +++++++++++++++++++++++++++
                // Update Linarization AND RHS
                // +++++++++++++++++++++++++++
                mtxBuilder.ComputeMatrix(OpMatrix, OpAffine); 
            } else {
                // ++++++++++++++++++++++++++++++++
                // Update only RHS of linearization
                // ++++++++++++++++++++++++++++++++
                mtxBuilder.ComputeAffine(OpAffine);
            }

            if(OpMatrix != null) {
                OpMatrix.CheckForNanOrInfM();
            }
            // if(BC requires than PRP){ //What exactly needed ?!?
            if(base.Control.PressureReferencePoint) {
                Console.WriteLine("Ref point is used");
                SetPressureReferencePoint(new double[] { 0.5, 0.5 }, this.CurrentSolution.Mapping, 3, OpMatrix, OpAffine);
                //}
            }
            //}
            OpAffine.CheckForNanOrInfV();

            if(Control.rMin < 10e-6) {
                R0fix myR0fix = new R0fix(Mapping, Control.rMin);

                if(OpMatrix != null) {
                    var OpMatrixMod = myR0fix.GetManipulatedMatrix(OpMatrix);
                    OpMatrix.Clear();
                    OpMatrix.Acc(1.0, OpMatrixMod);
                }
                var OpAffineMod = myR0fix.GetRestrictedRHS(OpAffine);


                OpAffine.ClearEntries();
                OpAffine.AccV(1.0, OpAffineMod);
            }
            //Console.WriteLine("Auskommetieren von R0 fix");
        }

        /// <summary>
        /// Might be moved to control object? Not sure...
        /// </summary>
        bool useSchur = true;

        /// <summary>
        /// configuration options for <see cref="MultigridOperator"/>.
        /// </summary>
        MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                int pVel = this.ur.Basis.Degree;
                int pPrs = this.Pressure.Basis.Degree;
                int D = this.GridData.SpatialDimension;
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[1][];
                for(int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    if(useSchur) {

                        // using a Schur complement for velocity & pressure
                        var confMomConti = new MultigridOperator.ChangeOfBasisConfig() {
                            VarIndex = new[] { 0, 1, 2, 3 },
                            DegreeS = new[] { pVel, pVel, pVel, pPrs },
                            mode = MultigridOperator.Mode.SchurComplement

                        };

                        configs[iLevel] = new[] { confMomConti };
                    } else {
                        configs[iLevel] = new MultigridOperator.ChangeOfBasisConfig[4];
                        // configurations for velocity
                        for(int d = 0; d < D + 1; d++) {
                            configs[iLevel][d] = new MultigridOperator.ChangeOfBasisConfig() {
                                DegreeS = new int[] { Math.Max(1, pVel) },
                                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib, //this.Control.VelocityBlockPrecondMode,
                                //mode = MultigridOperator.Mode.Eye, 

                                VarIndex = new int[] { d }
                            };
                        }
                        // configuration for pressure
                        configs[iLevel][D + 1] = new MultigridOperator.ChangeOfBasisConfig() {
                            DegreeS = new int[] { Math.Max(0, pPrs) },
                            mode = MultigridOperator.Mode.Eye, // this.Control.PressureBlockPrecondMode,
                            VarIndex = new int[] { D + 1 }
                        };
                    }
                }

                return configs;
            }
        }



        protected override void SetInitial(double t) {

            base.SetInitial(t);    // Anfangswerte werden gesetzt

            LsTrk.UpdateTracker(t);
            CreateEquationsAndSolvers(null);
            After_SetInitialOrLoadRestart(t);

            double dt = Control.GetFixedTimestep();
            double Time = 0.0;

            // Initialize DGFields
            var fieldNames = new string[] { "ur", "uxi", "ueta", "psi" };
            var exactFields_n = fieldNames.Select(name => new SinglePhaseField(CurrentSolution.Fields[Array.IndexOf(fieldNames, name)].Basis, $"{name}Exact_n")).ToArray();
            var exactFields_nn = fieldNames.Select(name => new SinglePhaseField(CurrentSolution.Fields[Array.IndexOf(fieldNames, name)].Basis, $"{name}Exact_nn")).ToArray();
            var guessFields = fieldNames.Select(name => new SinglePhaseField(CurrentSolution.Fields[Array.IndexOf(fieldNames, name)].Basis, $"{name}Guess_n")).ToArray();

            if(this.Control.ExactResidual == false) {
                UpdateGuessFields(guessFields, dt);
            } else {
                UpdateExactFields(exactFields_n, exactFields_nn, Time, dt);
            }
        }

        private void UpdateGuessFields(DGField[] guessFields, double dt) {

            Func<double[], double> pressure = this.Control.InitialValues_Evaluators.Get("Pressure");
            Func<double[], double> velUR = this.Control.InitialValues_Evaluators.Get("ur");
            Func<double[], double> velUETA = this.Control.InitialValues_Evaluators.Get("ueta");
            Func<double[], double> velUXI = this.Control.InitialValues_Evaluators.Get("uxi");

            guessFields[0].ProjectField(velUR);
            guessFields[1].ProjectField(velUXI);
            guessFields[2].ProjectField(velUETA);
            guessFields[3].ProjectField(pressure);

            CoordinateVector _Un = new CoordinateVector(guessFields);
            m_Splitting_Timestepper.SetInitialExact(dt, 0.0, _Un, _Un);
        }

        private void UpdateExactFields(DGField[] exactFields_n, DGField[] exactFields_nn, double Time, double dt) {

            double t = Time;
            double a = Globals.a;
            double b = Globals.b;

            Func<double[], double, double> ExactSolution_psi;
            Func<double[], double, double> ExactSolution_uR;
            Func<double[], double, double> ExactSolution_uETA;
            Func<double[], double, double> ExactSolution_uXI;

            if(Control.steady == true) {
                // Steady Man_Sol of DDD Paper
                ExactSolution_psi = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]);
                ExactSolution_uR = (X, t) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]);
                ExactSolution_uETA = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]);
                ExactSolution_uXI = (X, t) => 0.2e1 * X[0] * X[0] * Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1])
                + Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]);
            } else {
                // Transient Man_Sol of DDD Paper
                ExactSolution_psi = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t);
                ExactSolution_uR = (X, t) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t);
                ExactSolution_uETA = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(t);
                ExactSolution_uXI = (X, t) => (0.2e1 * X[0] * X[0] * Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) 
                + Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(t);
            }


            exactFields_n[0].ProjectField(ExactSolution_uR.Vectorize(Time - dt));
            exactFields_n[1].ProjectField(ExactSolution_uXI.Vectorize(Time - dt));
            exactFields_n[2].ProjectField(ExactSolution_uETA.Vectorize(Time - dt));
            exactFields_n[3].ProjectField(ExactSolution_psi.Vectorize(Time - dt));

            CoordinateVector _Un = new CoordinateVector(exactFields_n);

            exactFields_nn[0].ProjectField(ExactSolution_uR.Vectorize(Time - 2 * dt));
            exactFields_nn[1].ProjectField(ExactSolution_uXI.Vectorize(Time - 2 * dt));
            exactFields_nn[2].ProjectField(ExactSolution_uETA.Vectorize(Time - 2 * dt));
            exactFields_nn[3].ProjectField(ExactSolution_psi.Vectorize(Time - 2 * dt));

            CoordinateVector _Unn = new CoordinateVector(exactFields_nn);
            m_Splitting_Timestepper.SetInitialExact(dt, Time, _Unn, _Un);
        }

        private void After_SetInitialOrLoadRestart(double time) {
            using(new FuncTrace()) { 
                int D = this.GridData.SpatialDimension;

                this.DGLevSet.Current.Clear();
                this.DGLevSet.Current.AccLaidBack(1.0, this.LevSet);
                this.LsTrk.UpdateTracker(time);
                this.DGLevSet.IncreaseHistoryLength(1);
                this.DGLevSet.Push();

            }
        }
    }


}
