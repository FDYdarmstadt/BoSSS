using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Application.BoSSSpad;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MathNet.Numerics.Statistics;
using MPI.Wrappers;
using NUnit.Framework;
using StokesHelical_Ak.MomentumEquations;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing.Text;
using System.IO;
using System.Linq;
using System.Reflection.Metadata;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Voronoi;
using System.Diagnostics.Metrics;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using StokesHelical_Ak.TestSpartial;
using StokesHelical_Ak.TestTransient;
using StokesHelical_Ak.NUnitTestsR0_fix;

namespace StokesHelical_Ak {
    public partial class HelicalMain : Application<HelicalControl> {

        /// <summary>
        /// Entry point of my program
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {

            InitMPI();

            if(ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                Console.Write("rm");
                foreach(var pltFile in dir.GetFiles("*.plt").Concat(dir.GetFiles("*.curve"))) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }



            var c = StokesHelical_Ak.DNS_Hagen_Poiseulle.HagenPoiseulle(degree: 4, noOfCellsR: 32, noOfCellsXi: 32);
            c.ImmediatePlotPeriod = 1;
            c.NoOfTimesteps = 5;
            var solver = new HelicalMain();
            solver.Init(c);
            solver.RunSolverMode();
            Process.Start("mpiexec");
            ////StokesHelical_Ak.DNS_Centrifuge.Centrifuge_Flow();
            ////Restart_Comparison_Regular_Grid_BDF3_with_R0fix
            //// StokesHelical_Ak.HardcodedControl.TSFP();
            ////StokesHelical_Ak.TestTransient.TestTransient.PseudoSteadyCentrifuge();
            ////StokesHelical_Ak.NUnitTestsR0_fix.R0_fix_Test.TestR0_Fix(2);
            //FinalizeMPI();
            //System.Environment.Exit(20072022);

            _Main(args, false, delegate () {
                if(ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                    var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                    Console.Write("rm");
                    foreach(var pltFile in dir.GetFiles("*.plt").Concat(dir.GetFiles("*.curve"))) {
                        Console.Write(" " + pltFile.Name);
                        pltFile.Delete();
                    }
                    Console.WriteLine(";");
                }
                return new HelicalMain();
            });

            FinalizeMPI();
            System.Environment.Exit(20072022);
        }


        // Instantiations
        //=========================================
        protected int dimension = 2;

        /// <summary>
        /// Implicit part of the spatial operator;
        ///   this is only the Stokes operator,
        /// </summary>
        DifferentialOperator diffOp_implicit;

        /// <summary>
        /// See also <see cref="diffOp_implicit"/>;
        /// In the case of operator splitting, the convective operator; otherwise empty or null.
        /// </summary>
        DifferentialOperator diffOp_explicit;

        SplittingTimestepper m_Splitting_Timestepper;

        IEvaluatorNonLin ExplicitEval;

        public double error_Conti;
        public double error_MomZ;
        public double error_MomETA;
        public double error_MomR;

        public double urErrorL2;
        public double uetaErrorL2;
        public double uxiErrorL2;
        public double psiErrorL2;

        public double urErrorLx;
        public double uetaErrorLx;
        public double uxiErrorLx;
        public double psiErrorLx;

        public double[] Unnn_restart_;
        public double[] Unn_restart_;
        public double[] Un_restart_;

        public double[] C_nnn_restart;
        public double[] C_nn_restart;

        SinglePhaseField urError;
        SinglePhaseField uxiError;
        SinglePhaseField uetaError;
        SinglePhaseField psiError;

        // Attributes for fields (Names), initialization of DG fields
        //==============================================================

        // Velocity

        /// <summary>
        /// Pressure
        /// </summary>
        [InstantiateFromControlFile("Pressure", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField Pressure;

        /// <summary>
        /// ur
        /// </summary>
        [InstantiateFromControlFile("ur", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ur;

        /// <summary>
        /// ur
        /// </summary>
        [InstantiateFromControlFile("uxi", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField uxi;

        /// <summary>
        /// \f$ e^{\eta} \f$
        /// </summary>
        [InstantiateFromControlFile("ueta", null, IOListOption.ControlFileDetermined)]
        public SinglePhaseField ueta;

        // for manufactured solutions
        SinglePhaseField Residual_MomZ;
        SinglePhaseField Residual_Conti;
        SinglePhaseField Residual_MomETA;
        SinglePhaseField Residual_MomR;

        CoordinateVector m_UnknownsVector = null;

        ///// <summary>
        ///// The dependent variables
        ///// </summary>
        //CoordinateVector UnknownsVector {
        //    get {
        //        if(m_UnknownsVector == null)
        //            m_UnknownsVector = new CoordinateVector(ur, uxi, ueta, Pressure);
        //        return m_UnknownsVector;
        //    }
        //}

        //CoordinateMapping UnknownsMap {
        //    get {
        //        return UnknownsVector.Mapping;
        //    }
        //}


        CoordinateVector m_ResidualVector = null;
        /// <summary>
        /// The residuals of the continuity and the momentum equations
        /// </summary>
        CoordinateVector ResidualVector {
            get {
                if(m_ResidualVector == null)
                    m_ResidualVector = new CoordinateVector(Residual_MomR, Residual_MomZ, Residual_MomETA, Residual_Conti);
                return m_ResidualVector;
            }
        }

        public CoordinateMapping ResidualMap {
            get {
                return ResidualVector.Mapping;
            }
        }
        Stopwatch stopwatch = new Stopwatch();
        /// <summary>
        /// Errors of the equations
        /// </summary>
        SinglePhaseField Error_MomR;
        SinglePhaseField Error_MomZ;
        SinglePhaseField Error_Conti;
        SinglePhaseField Error_MomETA;

        CoordinateVector m_ErrorVector = null;

        CoordinateVector ErrorVector {
            get {
                if(m_ErrorVector == null) {
                    m_ErrorVector = new CoordinateVector(Error_MomR, Error_MomZ, Error_MomETA, Error_Conti);
                }
                return m_ErrorVector;
            }
        }

        public CoordinateMapping ErrorMap {
            get {
                return ErrorVector.Mapping;
            }
        }

        CoordinateVector m_CurrentSolution = null;

        public CoordinateVector CurrentSolution {
            get {
                if(m_CurrentSolution == null) {
                    m_CurrentSolution = new CoordinateVector(ur, uxi, ueta, this.Pressure);
                }
                return m_CurrentSolution;
            }
        }

        CoordinateVector m_CurrentResidual = null;

        public CoordinateVector CurrentResidual {
            get {
                if(m_CurrentResidual == null) {
                    m_CurrentResidual = new CoordinateVector(this.Residual_MomR, this.Residual_MomZ, this.Residual_MomETA, this.Residual_Conti);
                }
                return m_CurrentResidual;
            }
        }

        //Pressure Reference Point
        //============================================================//
        /// <summary>
        /// modifies a matrix <paramref name="Mtx"/> and a right-hand-side <paramref name="rhs"/>
        /// in order to fix the pressure at some reference point
        /// </summary>
        /// <param name="map">row mapping for <paramref name="Mtx"/> as well as <paramref name="rhs"/></param>
        /// <param name="iVar">the index of the pressure variable in the mapping <paramref name="map"/>.</param>
        /// <param name="Mtx"></param>
        /// <param name="rhs"></param>
        static public void SetPressureReferencePoint<T>(double[] Point, UnsetteledCoordinateMapping map, int iVar, BlockMsrMatrix Mtx, T rhs)
            where T : IList<double> {
            using(new FuncTrace()) {
                var GridDat = map.GridDat;

                if(rhs.Count != map.LocalLength)
                    throw new ArgumentException();
                if(Mtx != null) {
                    if(!Mtx.RowPartitioning.EqualsPartition(map) || !Mtx.ColPartition.EqualsPartition(map))
                        throw new ArgumentException();
                }

                GridDat.LocatePoint(Point, out _, out long GlobalIndex, out _, out bool onthisProc);

                long iRowGl = -111; // 
                if(onthisProc) {
                    iRowGl = map.GlobalUniqueCoordinateIndex_FromGlobal(iVar, GlobalIndex, 0);
                }
                iRowGl = iRowGl.MPIMax(); // on all processes where this NOT set, the input is negative; thus we select the value from the respective processor, where it is set

                // clear row
                // ---------
                if(onthisProc) {
                    // ref. cell is on local MPI process
                    // set matrix row to identity
                    if(Mtx != null) {
                        Mtx.ClearRow(iRowGl);
                        Mtx.SetDiagonalElement(iRowGl, 1.0);
                    }

                    // clear RHS
                    int iRowLoc = (int)(iRowGl - map.i0);
                    rhs[iRowLoc] = 0;

                    Console.WriteLine("Local row: " + iRowLoc);
                }

                // clear column
                // ------------
                if(Mtx != null) {
                    for(long i = Mtx.RowPartitioning.i0; i < Mtx.RowPartitioning.iE; i++) {
                        if(i != iRowGl) {
                            Mtx[i, iRowGl] = 0;
                        }
                    }
                }
            }
        }
        //==============================================================//

        double[] Compute_C_restart(double[] u) {
            double[] C = new double[u.Length];
            ExplicitEval.Evaluate(1.0, 0.0, C);
            return C;
        }

        protected override void LoadRestart(out double Time, out TimestepNumber TimestepNo) {

            if(this.Control.GetBDFOrder() == 1) {
                base.LoadRestart(out Time, out TimestepNo);
                this.Control.restartTimeStep = TimestepNo.MajorNumber;
                return;
            }

            Time = 0;
            TimestepNo = 0;

            if(MPIRank == 0) {
                Console.WriteLine("----------------------------------------");
                Console.WriteLine("Loading Restart...");
            }

            if(this.Control != null && this.Control.RestartInfo != null) {

                if(this.Control.GetBDFOrder() == 3) {

                    TimestepNo = RestartFromDatabase(out Time);
                    this.CurrentSessionInfo.RestartedFrom = this.Control.RestartInfo.Item1;
                    this.CurrentSessionInfo.Save();
                    Unnn_restart_ = this.CurrentSolution.ToArray();

                    this.Control.RestartInfo = new Tuple<Guid, TimestepNumber>(this.CurrentSessionInfo.RestartedFrom, new TimestepNumber(TimestepNo.MajorNumber + 1));
                    TimestepNo = RestartFromDatabase(out Time);
                    this.CurrentSessionInfo.RestartedFrom = this.Control.RestartInfo.Item1;
                    this.CurrentSessionInfo.Save();
                    Unn_restart_ = this.CurrentSolution.ToArray();

                    this.Control.RestartInfo = new Tuple<Guid, TimestepNumber>(this.CurrentSessionInfo.RestartedFrom, new TimestepNumber(TimestepNo.MajorNumber + 1));
                    TimestepNo = RestartFromDatabase(out Time);
                    this.CurrentSessionInfo.RestartedFrom = this.Control.RestartInfo.Item1;
                    this.CurrentSessionInfo.Save();
                    Un_restart_ = this.CurrentSolution.ToArray();
                    this.Control.restartTimeStep = TimestepNo.MajorNumber;
                } else {
                    throw new Exception("ERRROR BDF Sceme not implemented");
                }

            } else {
                throw new ApplicationException("restart was not specified");
                //throw new NotImplementedException("unknown restart type.");
            }

            PostRestart(Time, TimestepNo);

            if(this.LsTrk != null) {
                if(this.LsTrk.Regions.Time != Time)
                    this.LsTrk.UpdateTracker(Time);
            }

            System.GC.Collect();
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

            if(MPIRank == 0) {
                Console.WriteLine("Finished loading restart.");
                Console.WriteLine("----------------------------------------");
            }

        }






        //Penalty Factor
        //============================================================//
        /// <summary>
        /// Berechne penalty factor
        /// </summary>
        /// <param name="penaltyFactor">penaltySafety * p *p</param>
        /// <param name="jCellIn">Zell index linke zelle</param>
        /// <param name="jCellOut"></param>
        /// <param name="cj">Oberflaeche pro Volumen, benoetigt fuer nicht kartesische zellen</param>
        /// <returns></returns>
        public double PenaltyFactor(double penaltyFactor, int jCellIn, int jCellOut, MultidimensionalArray cj) {
            double cj_in = cj[jCellIn];
            double eta = penaltyFactor * cj_in;
            if(jCellOut >= 0) {
                double cj_out = cj[jCellOut];
                eta = Math.Max(eta, penaltyFactor * cj_out);
            }
            return eta;
        }
        //============================================================//

        /// <summary>
        /// creates the field
        /// <see cref="c"/>;
        /// </summary>
        protected override void CreateFields() {
            // Start the clock
            stopwatch.Start();
            //
            Debug.Assert(MyLsTrk != null);
            base.LsTrk = MyLsTrk;
            Residual_Conti = new SinglePhaseField(Pressure.Basis, "Residual_Conti");
            m_IOFields.Add(Residual_Conti);
            Error_Conti = new SinglePhaseField(Pressure.Basis, "Error_Conti");
            m_IOFields.Add(Error_Conti);
            psiExact = new SinglePhaseField(Pressure.Basis, "presExact");
            m_IOFields.Add(psiExact);
            psiError = new SinglePhaseField(Pressure.Basis, "presError");
            m_IOFields.Add(psiError);

            Residual_MomR = new SinglePhaseField(ur.Basis, "Residual_MomR");
            Residual_MomZ = new SinglePhaseField(uxi.Basis, "Residual_MomZ");
            Residual_MomETA = new SinglePhaseField(ueta.Basis, "Residual_MomETA");
            m_IOFields.Add(Residual_MomZ);
            m_IOFields.Add(Residual_MomETA);
            m_IOFields.Add(Residual_MomR);

            Error_MomZ = new SinglePhaseField(uxi.Basis, "Error_MomZ");
            Error_MomETA = new SinglePhaseField(ueta.Basis, "Error_MomETA");
            Error_MomR = new SinglePhaseField(ur.Basis, "Error_MomR");
            m_IOFields.Add(Error_MomZ);
            m_IOFields.Add(Error_MomETA);
            m_IOFields.Add(Error_MomR);


            // Exact solution
            urExact = new SinglePhaseField(ur.Basis, "urExact");
            uetaExact = new SinglePhaseField(ueta.Basis, "uetaExact");
            uxiExact = new SinglePhaseField(uxi.Basis, "uxiExact");

            m_IOFields.Add(urExact);
            m_IOFields.Add(uetaExact);
            m_IOFields.Add(uxiExact);

            // Errors of velocities and pressure
            urError = new SinglePhaseField(ur.Basis, "urError");
            uetaError = new SinglePhaseField(ueta.Basis, "uetaError");
            uxiError = new SinglePhaseField(uxi.Basis, "uxiError");

            m_IOFields.Add(urError);
            m_IOFields.Add(uetaError);
            m_IOFields.Add(uxiError);

        }

        /// <summary>
        /// DG Fields of exact solutions
        /// </summary>
        DGField urExact;
        DGField uetaExact;
        DGField uxiExact;
        DGField psiExact;

        /// <summary>
        /// Setzt Operator Matrix aus Konti und Impulsgleichungen zusammen
        /// </summary>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            Console.WriteLine("Mit den RB auf Muss der PRP auf true sein!");//
            this.Control.PressureReferencePoint = true;
            if(Control.rMin < 10e-6) {
                this.Control.R0fixOn = true;
                if(Control.ExactResidual == true) {
                    if(Control.steady == true) {
                        // Das ist hier noch nicht so sicher!
                        Globals.activeMult = Globals.Multiplier.Bsq;
                    } else {
                        Globals.activeMult = Globals.Multiplier.one;
                        Console.WriteLine($"CreateEquationsAndSolvers need Debuggen because of Globals.Multiplier {Globals.activeMult} and rMin = {Control.rMin}");
                    }
                } else {
                    Globals.activeMult = Globals.Multiplier.Bsq;
                }
            } else {
                Globals.activeMult = Globals.Multiplier.one;
            }
            if(L != null)
                throw new NotSupportedException("No support for dynamic load balancing / dynamic mesh adaptation.");

            if(this.diffOp_implicit != null) {
                // already initialized
                return;
            }

            Globals.DirichletValue_uR = (X => this.Control.BoundaryValues["Dirichlet"].Values["ur"].Evaluate(X, 0));
            Globals.DirichletValue_uXi = (X => this.Control.BoundaryValues["Dirichlet"].Values["uxi"].Evaluate(X, 0));
            Globals.DirichletValue_uEta = (X => this.Control.BoundaryValues["Dirichlet"].Values["ueta"].Evaluate(X, 0));
            Globals.DirichletValue_psi = (X => this.Control.BoundaryValues["Dirichlet"].Values["Pressure"].Evaluate(X, 0));
            Globals.BoundaryType = X => BoundaryTypeE.Dirichlet;

            Globals.pressureStabilConti = true;
            Globals.pressureStabilEtaMom = false;

            diffOp_implicit = new DifferentialOperator(
            new string[] { "ur", "uxi", "ueta", "Pressure" },
            new string[0],   //no parameters
            new string[] { "rmom", "zmom", "etamom", "konti" },
            QuadOrderFunc.NonLinear(3));
            var tmpOp = new DependentTemporalOperator(diffOp_implicit);

            tmpOp.EquationComponents["rmom"].Add(new TransientTerm(1, 0));
            tmpOp.EquationComponents["zmom"].Add(new TransientTerm(1, 1));
            tmpOp.EquationComponents["etamom"].Add(new TransientTerm(1, 2));

            diffOp_implicit.TemporalOperator = tmpOp;
            //#########################################
            // TEMPORAL OPERATOR WITH BSQ NEEDED for transient and rMin=0
            //#########################################
            diffOp_explicit = new DifferentialOperator(
                diffOp_implicit.DomainVar, 
                new string[] { "ur0", "uxi0", "ueta0" }, 
                diffOp_implicit.CodomainVar, QuadOrderFunc.NonLinear(3));
            int order = Control.FieldOptions["ur"].Degree;

            double penalty = this.Control.penaltySafety * order * order;
            int noOfCells = this.Control.Resolution_Xi;
            Conti myContiNew = new Conti(noOfCells);
            diffOp_implicit.EquationComponents["konti"].Add(myContiNew);

            diffOp_implicit.EquationComponents["rmom"].Add(new Gradient_pressure_rMomentum(0));
            diffOp_implicit.EquationComponents["zmom"].Add(new Gradient_pressure_xiMomentum(1));

            diffOp_implicit.EquationComponents["etamom"].Add(new etaMomentum(noOfCells, this.Control.TermSwitch, penalty, PenaltyFactor));
            diffOp_implicit.EquationComponents["rmom"].Add(new rMomentum(this.Control.TermSwitch, penalty, PenaltyFactor));
            diffOp_implicit.EquationComponents["zmom"].Add(new xiMomentum(this.Control.TermSwitch, penalty, PenaltyFactor));
            // +++++++++++++++++++++++++
            // Navier-Stokes (nonlinear)
            // +++++++++++++++++++++++++
            if(Control.NavierStokes) {
                diffOp_explicit.EquationComponents["rmom"].Add(new convectiveRmom());
                diffOp_explicit.EquationComponents["etamom"].Add(new convectiveETAmom());
                diffOp_explicit.EquationComponents["zmom"].Add(new convectiveXImom());
            } else {
                Console.WriteLine("Solvinng only Stokes/convective terms deactivated");
            }
            // +++++++++++++++++++++++++
            // Exact Solution
            // +++++++++++++++++++++++++
            if(Control.ExactResidual) {
                diffOp_implicit.EquationComponents["rmom"].Add(new StokesHelical_Ak.ForcingTerms.ManSol.ForcingTermR());
                diffOp_implicit.EquationComponents["zmom"].Add(new StokesHelical_Ak.ForcingTerms.ManSol.ForcingTermXi());
                diffOp_implicit.EquationComponents["etamom"].Add(new StokesHelical_Ak.ForcingTerms.ManSol.ForcingTermEta());
                diffOp_implicit.EquationComponents["konti"].Add(new StokesHelical_Ak.ForcingTerms.ManSol.ForcingTermConti());
            }
            // Set Amplitude of HagenPoeuseulle, or rotating disc. Therefore, not in if statement
            Globals.MaxAmp = this.Control.maxAmpli;
            if(Control.HagenPoisseulle) {
                diffOp_implicit.EquationComponents["zmom"].Add(new StokesHelical_Ak.ForcingTerms.Hagen_Poiseulle.ForcingTermXi_Hagen());
                diffOp_implicit.EquationComponents["etamom"].Add(new StokesHelical_Ak.ForcingTerms.Hagen_Poiseulle.ForcingTermEta_Hagen());
            }

            diffOp_implicit.Commit();
            if(diffOp_explicit != null) {
                diffOp_explicit.Commit();
            }
            List<DGField> DepVars = new List<DGField> { };
            List<DGField> ResList = new List<DGField> { };
            DepVars.Add(ur);
            DepVars.Add(uxi);
            DepVars.Add(ueta);
            DepVars.Add(Pressure);
            ResList.Add(Residual_MomR);
            ResList.Add(Residual_MomZ);
            ResList.Add(Residual_MomETA);
            ResList.Add(Residual_Conti);

            m_Splitting_Timestepper = new SplittingTimestepper(this.DelComputeOperatorMatrixMod, 
                diffOp_explicit, diffOp_implicit.TemporalOperator, 
                this.CurrentSolution.Mapping, 
                this.Control.GetFixedTimestep(), this.Control.GetBDFOrder(), 
                this.Control.LinearSolver, this.MultigridOperatorConfig, base.MultigridSequence, 
                myR0fix, this.Control.rMin, this.Control.PressureReferencePoint);
        }

        protected override void PlotCurrentState(double physTime, BoSSS.Foundation.IO.TimestepNumber timestepNo, int superSampling = 0) {
            using(new FuncTrace()) {
                Tecplot plt1 = new Tecplot(GridData, true, false, (uint)superSampling);
                List<DGField> FieldsToPlot;

                if(Control.ExactResidual) {
                    FieldsToPlot = m_IOFields.ToList();
                    FieldsToPlot = new List<DGField> { ur, uxi, ueta, Residual_MomR, Residual_MomZ, Residual_MomETA };
                } else {
                    FieldsToPlot = m_IOFields.ToList();
                }
                string time_dep;
                if(this.Control.steady == true) {
                    time_dep = "Steady_";
                } else {
                    time_dep = "Transient_";
                }
                DateTime now = DateTime.Now;
                string dateString = now.ToString("dd.MM.yyyy");

                Assert.IsTrue(Control.PressureReferencePoint, "Pressure Refference Point should be true");

                if(Control.rMin < 10e-6) {
                    Assert.IsTrue(Control.R0fixOn, "R0_fix should be thrue");
                }

                if(Globals.activeMult == Globals.Multiplier.one && this.Control.rMin < 10e-6) {
                    Console.WriteLine("Friendly Reminder: Mutiplier One and rMin<10e-6");
                } else if(Globals.activeMult == Globals.Multiplier.Bsq && this.Control.rMin > 10e-6) {
                    throw new Exception("Mutiplier BSQ and rMin>10e-6");
                }

                string nameOf_Plot = "Hel__BDF_" + this.Control.GetBDFOrder() +"_"+this.Control.grid+ "_LinSol" + this.Control.LinearSolver.Shortname + "_rMin_" + this.Control.rMin.ToString() + "_Mult_" + Globals.activeMult.ToString() + "_PRP_" + this.Control.PressureReferencePoint + "_" + this.Control.Resolution_R.ToString() + "_x_" +
                              this.Control.Resolution_R.ToString() + "_DGd_" + this.Control.dg_degree + "_Amp_" + this.Control.maxAmpli.ToString() + "_dt_" + this.Control.dtMax + "_TiStNo_" + timestepNo;

                plt1.PlotFields(nameOf_Plot, physTime, m_IOFields);

                Console.WriteLine("Helical solution plotted \n");
            }
        }

        R0fix m_myR0fix;

        /// <summary>
        /// Cached version of the rather expensive <see cref="R0fix"/>
        /// </summary>
        R0fix myR0fix {
            get {
                if(m_myR0fix == null || !object.ReferenceEquals(m_myR0fix.GridDat, this.GridData))
                    m_myR0fix = new R0fix(this.CurrentSolution.Mapping, this.Control.rMin);
                return m_myR0fix;
            }
        }

        protected override double RunSolverOneStep(int TimestepInt, double phystime, double dt) {
            using(new FuncTrace()) {

                TimestepNumber TimestepNo = new TimestepNumber(TimestepInt, 0);
                int D = this.GridData.SpatialDimension;

                base.ResLogger.TimeStep = TimestepInt;

                dt = base.GetTimestep();
                Console.WriteLine("Instationary solve, timestep #{0}, dt = {1} ...", TimestepNo, dt);

                int p = ur.Basis.Degree;
                double dtCFL = this.GridData.ComputeCFLTime(new DGField[] { ur, uxi, ueta }, 1000) / (p * p);
                Console.WriteLine("CFL fraction: " + dt / dtCFL);

                if(Control.rMin < 10e-6) {
                    if(Globals.activeMult== Globals.Multiplier.one) {
                        Console.WriteLine("Friendly Reminder: Mutiplier One and rMin<10e-6");
                    }
                    Assert.IsTrue(Control.R0fixOn, "R0_fix should be true");
                    m_Splitting_Timestepper.Solve(phystime + dt, myR0fix, this.Control.restartTimeStep, TimestepInt, this.Control.RestartInfo != null, Unnn_restart_, Unn_restart_, Un_restart_, this.Control.GetBDFOrder());
                    // Ruecktransformation: note, fk, 11apr24: moved into Solve
                    myR0fix.CheckSolutionR0Compatibility(this.CurrentSolution);
                } else {
                    m_Splitting_Timestepper.Solve(phystime + dt, this.Control.restartTimeStep, TimestepInt, this.Control.RestartInfo != null, Unnn_restart_, Unn_restart_, Un_restart_, this.Control.GetBDFOrder());
                }
                if(this.Control.ExactResidual == true) {
                    ComputeErrorsAndResiduals(phystime + dt, TimestepNo);
                    // Beim Restart werden hier irgendwie immer neue Error_Felder angelegt
                }
                //Stopclock
                stopwatch.Stop();
                Console.WriteLine("Time elapsed: {0}ms", stopwatch.ElapsedMilliseconds);
                stopwatch.Restart();
                //Stopclock
                this.ResLogger.NextTimestep(false);
                CoordinateMapping helicalSol = new CoordinateMapping(ur, uxi, ueta, Pressure);
                CoordinateVector helicalSolVec = new CoordinateVector(helicalSol);
                helicalSolVec.SaveToTextFile($"After_PC_Solution_StokesSystem_with_dt={dt}.txt");
                helicalSolVec.SaveToTextFile("helicalSolution.txt");
                return dt;
            }
        }
        void ComputeErrorsAndResiduals(double time, TimestepNumber timestepNo) {
            using(new FuncTrace()) {
                var ResidualFields = ResidualMap.Fields;
                var ErrorFields = ErrorMap.Fields;

                double t = time;
                double nu = Globals.nu;
                double a = Globals.a;
                double b = Globals.b;

                Func<double[], double, double> ExactResidualConti;
                Func<double[], double, double> ExactResidualRmom;
                Func<double[], double, double> ExactResidualETAmom;
                Func<double[], double, double> ExactResidualXImom;

                Func<double[], double, double> ExactSolution_psi;
                Func<double[], double, double> ExactSolution_uR;
                Func<double[], double, double> ExactSolution_uETA;
                Func<double[], double, double> ExactSolution_uXI;



                if(Control.steady == true) {
                    // Steady
                    ExactSolution_psi = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]);
                    ExactSolution_uR = (X, t) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]);
                    ExactSolution_uETA = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]);
                    ExactSolution_uXI = (X, t) => 0.2e1 * X[0] * X[0] * Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1])
                    + Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]);
                    if(Globals.activeMult == Globals.Multiplier.one) {
                        // Globals.Multiplier.one
                        ExactResidualConti = (X, t) => 0;
                        ExactResidualRmom = (X, t) => (0.2e1 * X[0] * ((nu * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) + (-0.2e1 * Math.Pow(X[0], 0.3e1) + 0.2e1 * X[0]) * Math.Pow(Math.Cos(X[1]), 0.2e1)) * Math.Exp(-X[0] * X[0]) + (0.2e1 * Math.Pow(X[0], 0.3e1) - X[0]) * Math.Pow(Math.Cos(X[1]), 0.2e1) * Math.Exp(-0.2e1 * X[0] * X[0]) - nu * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) - X[0] * Math.Pow(Math.Cos(X[1]), 0.2e1)) * a * b * Math.Sqrt(a * a * X[0] * X[0] + b * b) + (-(-0.4e1 * a * a * nu * Math.Pow(X[0], 0.6e1) - 0.2e1 * a * a * Math.Pow(X[0], 0.5e1) + nu * (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(X[0], 0.4e1) - 0.2e1 * b * b * Math.Pow(X[0], 0.3e1) + 0.2e1 * ((a * a + 0.4e1) * b * b + a * a / 0.2e1) * nu * X[0] * X[0] + Math.Pow(b, 0.4e1) * nu - b * b * nu) * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) + 0.2e1 * X[0] * (-b * b * ((a * a + 0.2e1) * X[0] * X[0] + b * b - 0.1e1) * Math.Pow(Math.Cos(X[1]), 0.2e1) + Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * X[0] * X[0])) * Math.Exp(-X[0] * X[0]) - 0.2e1 * X[0] * (-(-0.4e1 * Math.Pow(X[0], 0.4e1) + (a * a + 0.4e1) * X[0] * X[0] + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(X[1]), 0.2e1) / 0.2e1 + Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * X[0] * X[0]) * Math.Exp(-0.2e1 * X[0] * X[0]) + (Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.4e1) + (0.2e1 * a * a * b * b + a * a) * X[0] * X[0] + Math.Pow(b, 0.4e1) - b * b) * (a * a * X[0] * X[0] + b * b) * nu * Math.Sin(X[1]) + Math.Pow(Math.Cos(X[1]), 0.2e1) * b * b * X[0] * (a * a * X[0] * X[0] + b * b - 0.1e1)) * Math.Pow(a * a * X[0] * X[0] + b * b, -0.2e1) * Math.Pow(X[0], -0.2e1);
                        ExactResidualETAmom = (X, t) => -Math.Pow(a * a * X[0] * X[0] + b * b, -0.5e1 / 0.2e1) * Math.Cos(X[1]) * ((((-0.2e1 * a * a * b * b * Math.Pow(X[0], 0.3e1) - 0.2e1 * Math.Pow(b, 0.4e1) * X[0]) * Math.Sin(X[1]) + (-0.4e1 * Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.8e1) + (Math.Pow(a, 0.6e1) + 0.4e1 * Math.Pow(a, 0.4e1) - 0.8e1 * a * a * b * b) * Math.Pow(X[0], 0.6e1) + (-0.4e1 * Math.Pow(b, 0.4e1) + (0.3e1 * Math.Pow(a, 0.4e1) + 0.8e1 * a * a) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(X[0], 0.4e1) + ((0.3e1 * a * a + 0.4e1) * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * X[0] * X[0] + Math.Pow(b, 0.6e1)) * nu) * Math.Exp(-X[0] * X[0]) + b * b * X[0] * Math.Sin(X[1]) * (a * a * X[0] * X[0] + b * b) * Math.Exp(-0.2e1 * X[0] * X[0]) + (a * a * b * b * Math.Pow(X[0], 0.3e1) + Math.Pow(b, 0.4e1) * X[0]) * Math.Sin(X[1]) - (Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.6e1) + (0.3e1 * Math.Pow(a, 0.4e1) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(X[0], 0.4e1) + (0.3e1 * a * a * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * X[0] * X[0] + Math.Pow(b, 0.6e1)) * nu) * Math.Sqrt(a * a * X[0] * X[0] + b * b) - 0.2e1 * X[0] * ((-0.4e1 * a * a * Math.Pow(X[0], 0.6e1) + (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(X[0], 0.4e1) + ((0.2e1 * a * a + 0.8e1) * b * b + a * a) * X[0] * X[0] + Math.Pow(b, 0.4e1) - b * b) * Math.Exp(-X[0] * X[0]) - Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.4e1) + (-0.2e1 * a * a * b * b - a * a) * X[0] * X[0] - Math.Pow(b, 0.4e1) + b * b) * a * b * nu) * Math.Pow(X[0], -0.2e1);
                        ExactResidualXImom = (X, t) => 0.2e1 * Math.Pow(a * a * X[0] * X[0] + b * b, -0.5e1 / 0.2e1) * Math.Cos(X[1]) * (-0.2e1 * (((a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) + X[0] * nu * (a * a * X[0] * X[0] + a * a / 0.2e1 + b * b)) * Math.Exp(-X[0] * X[0]) - (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) * Math.Exp(-0.2e1 * X[0] * X[0]) / 0.2e1 + Math.Sin(X[1]) * (-a * a * X[0] * X[0] / 0.2e1 - b * b / 0.2e1) - a * a * X[0] * nu / 0.2e1) * X[0] * X[0] * a * b * Math.Sqrt(a * a * X[0] * X[0] + b * b) + (-0.2e1 * Math.Pow(X[0], 0.3e1) * (X[0] - 0.1e1) * (X[0] + 0.1e1) * (a * a * X[0] * X[0] + a * a + b * b) * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) - 0.4e1 * Math.Pow(a, 0.4e1) * nu * Math.Pow(X[0], 0.10e2) + a * a * nu * (Math.Pow(a, 0.4e1) + 0.10e2 * a * a - 0.8e1 * b * b) * Math.Pow(X[0], 0.8e1) - Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.7e1) / 0.2e1 - (Math.Pow(a, 0.6e1) + (-0.6e1 * b * b + 0.2e1) * Math.Pow(a, 0.4e1) - 0.48e2 * a * a * b * b + 0.8e1 * Math.Pow(b, 0.4e1)) * nu * Math.Pow(X[0], 0.6e1) / 0.2e1 - 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(X[0], 0.5e1) - ((b * b - 0.1e1) * Math.Pow(a, 0.4e1) + (-0.6e1 * Math.Pow(b, 0.4e1) + 0.4e1 * b * b) * a * a - 0.28e2 * Math.Pow(b, 0.4e1)) * nu * Math.Pow(X[0], 0.4e1) / 0.2e1 - 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(X[0], 0.3e1) + ((b * b - 0.4e1) * a * a + 0.2e1 * Math.Pow(b, 0.4e1) - 0.10e2 * b * b) * b * b * nu * X[0] * X[0] / 0.2e1 - Math.Pow(b, 0.6e1) * X[0] / 0.2e1 + Math.Pow(b, 0.6e1) * nu / 0.2e1 - Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Exp(-X[0] * X[0]) - Math.Pow(X[0], 0.3e1) * Math.Sin(X[1]) * (a * a + 0.2e1 * b * b) * (a * a * X[0] * X[0] + b * b) * Math.Exp(-0.2e1 * X[0] * X[0]) + (-Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.5e1) - a * a * b * b * Math.Pow(X[0], 0.3e1)) * Math.Sin(X[1]) + Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.7e1) / 0.2e1 + Math.Pow(a, 0.6e1) * nu * Math.Pow(X[0], 0.6e1) / 0.2e1 + 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(X[0], 0.5e1) + Math.Pow(a, 0.4e1) * nu * (b - 0.1e1) * (b + 0.1e1) * Math.Pow(X[0], 0.4e1) / 0.2e1 + 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(X[0], 0.3e1) + (-Math.Pow(b, 0.4e1) * nu / 0.2e1 + 0.2e1 * b * b * nu) * a * a * X[0] * X[0] + Math.Pow(b, 0.6e1) * X[0] / 0.2e1 - Math.Pow(b, 0.6e1) * nu / 0.2e1 + Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Pow(X[0], -0.2e1);
                    } else {
                        // Globals.Multiplier.Bsq
                        ExactResidualConti = (X, t) => 0;
                        ExactResidualRmom = (X, t) => (0.2e1 * a * ((nu * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) + (-0.2e1 * Math.Pow(X[0], 0.3e1) + 0.2e1 * X[0]) * Math.Pow(Math.Cos(X[1]), 0.2e1)) * Math.Exp(-X[0] * X[0]) + (0.2e1 * Math.Pow(X[0], 0.3e1) - X[0]) * Math.Pow(Math.Cos(X[1]), 0.2e1) * Math.Exp(-0.2e1 * X[0] * X[0]) - nu * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) - X[0] * Math.Pow(Math.Cos(X[1]), 0.2e1)) * b * X[0] * Math.Sqrt(a * a * X[0] * X[0] + b * b) + (-(a * a * X[0] * X[0] + b * b) * (-0.4e1 * a * a * nu * Math.Pow(X[0], 0.6e1) - 0.2e1 * a * a * Math.Pow(X[0], 0.5e1) + nu * (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(X[0], 0.4e1) - 0.2e1 * b * b * Math.Pow(X[0], 0.3e1) + 0.2e1 * ((a * a + 0.4e1) * b * b + a * a / 0.2e1) * nu * X[0] * X[0] + Math.Pow(b, 0.4e1) * nu - b * b * nu) * Math.Sin(X[1]) + 0.2e1 * (-((a * a + 0.2e1) * X[0] * X[0] + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(X[1]), 0.2e1) + Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * X[0] * X[0]) * X[0]) * Math.Exp(-X[0] * X[0]) - 0.2e1 * (-b * b * (-0.4e1 * Math.Pow(X[0], 0.4e1) + (a * a + 0.4e1) * X[0] * X[0] + b * b - 0.1e1) * Math.Pow(Math.Cos(X[1]), 0.2e1) / 0.2e1 + Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * X[0] * X[0]) * X[0] * Math.Exp(-0.2e1 * X[0] * X[0]) + nu * (a * a * X[0] * X[0] + b * b) * (Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.4e1) + (0.2e1 * a * a * b * b + a * a) * X[0] * X[0] + Math.Pow(b, 0.4e1) - b * b) * Math.Sin(X[1]) + Math.Pow(Math.Cos(X[1]), 0.2e1) * b * b * X[0] * (a * a * X[0] * X[0] + b * b - 0.1e1)) * Math.Pow(a * a * X[0] * X[0] + b * b, -0.3e1);
                        ExactResidualETAmom = (X, t) => -Math.Pow(a * a * X[0] * X[0] + b * b, -0.7e1 / 0.2e1) * Math.Cos(X[1]) * ((((-0.2e1 * a * a * b * b * Math.Pow(X[0], 0.3e1) - 0.2e1 * Math.Pow(b, 0.4e1) * X[0]) * Math.Sin(X[1]) + (-0.4e1 * Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.8e1) + (Math.Pow(a, 0.6e1) + 0.4e1 * Math.Pow(a, 0.4e1) - 0.8e1 * a * a * b * b) * Math.Pow(X[0], 0.6e1) + (-0.4e1 * Math.Pow(b, 0.4e1) + (0.3e1 * Math.Pow(a, 0.4e1) + 0.8e1 * a * a) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(X[0], 0.4e1) + ((0.3e1 * a * a + 0.4e1) * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * X[0] * X[0] + Math.Pow(b, 0.6e1)) * nu) * Math.Exp(-X[0] * X[0]) + b * b * X[0] * Math.Sin(X[1]) * (a * a * X[0] * X[0] + b * b) * Math.Exp(-0.2e1 * X[0] * X[0]) + (a * a * b * b * Math.Pow(X[0], 0.3e1) + Math.Pow(b, 0.4e1) * X[0]) * Math.Sin(X[1]) - nu * (Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.6e1) + (0.3e1 * Math.Pow(a, 0.4e1) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(X[0], 0.4e1) + (0.3e1 * a * a * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * X[0] * X[0] + Math.Pow(b, 0.6e1))) * Math.Sqrt(a * a * X[0] * X[0] + b * b) - 0.2e1 * a * ((-0.4e1 * a * a * Math.Pow(X[0], 0.6e1) + (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(X[0], 0.4e1) + ((0.2e1 * a * a + 0.8e1) * b * b + a * a) * X[0] * X[0] + Math.Pow(b, 0.4e1) - b * b) * Math.Exp(-X[0] * X[0]) - Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.4e1) + (-0.2e1 * a * a * b * b - a * a) * X[0] * X[0] - Math.Pow(b, 0.4e1) + b * b) * nu * b * X[0]);
                        ExactResidualXImom = (X, t) => 0.2e1 * Math.Pow(a * a * X[0] * X[0] + b * b, -0.7e1 / 0.2e1) * Math.Cos(X[1]) * (-0.2e1 * a * (((a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) + X[0] * nu * (a * a * X[0] * X[0] + a * a / 0.2e1 + b * b)) * Math.Exp(-X[0] * X[0]) - (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) * Math.Exp(-0.2e1 * X[0] * X[0]) / 0.2e1 + Math.Sin(X[1]) * (-a * a * X[0] * X[0] / 0.2e1 - b * b / 0.2e1) - a * a * X[0] * nu / 0.2e1) * b * X[0] * X[0] * Math.Sqrt(a * a * X[0] * X[0] + b * b) + (-0.2e1 * Math.Pow(X[0], 0.3e1) * (X[0] - 0.1e1) * (X[0] + 0.1e1) * (a * a * X[0] * X[0] + a * a + b * b) * (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) - 0.4e1 * Math.Pow(a, 0.4e1) * nu * Math.Pow(X[0], 0.10e2) + a * a * nu * (Math.Pow(a, 0.4e1) + 0.10e2 * a * a - 0.8e1 * b * b) * Math.Pow(X[0], 0.8e1) - Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.7e1) / 0.2e1 - nu * (Math.Pow(a, 0.6e1) + (-0.6e1 * b * b + 0.2e1) * Math.Pow(a, 0.4e1) - 0.48e2 * a * a * b * b + 0.8e1 * Math.Pow(b, 0.4e1)) * Math.Pow(X[0], 0.6e1) / 0.2e1 - 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(X[0], 0.5e1) - nu * ((b * b - 0.1e1) * Math.Pow(a, 0.4e1) + (-0.6e1 * Math.Pow(b, 0.4e1) + 0.4e1 * b * b) * a * a - 0.28e2 * Math.Pow(b, 0.4e1)) * Math.Pow(X[0], 0.4e1) / 0.2e1 - 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(X[0], 0.3e1) + ((b * b - 0.4e1) * a * a + 0.2e1 * Math.Pow(b, 0.4e1) - 0.10e2 * b * b) * nu * b * b * X[0] * X[0] / 0.2e1 - Math.Pow(b, 0.6e1) * X[0] / 0.2e1 + Math.Pow(b, 0.6e1) * nu / 0.2e1 - Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Exp(-X[0] * X[0]) - Math.Pow(X[0], 0.3e1) * Math.Sin(X[1]) * (a * a + 0.2e1 * b * b) * (a * a * X[0] * X[0] + b * b) * Math.Exp(-0.2e1 * X[0] * X[0]) + (-Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.5e1) - a * a * b * b * Math.Pow(X[0], 0.3e1)) * Math.Sin(X[1]) + Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.7e1) / 0.2e1 + Math.Pow(a, 0.6e1) * nu * Math.Pow(X[0], 0.6e1) / 0.2e1 + 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(X[0], 0.5e1) + Math.Pow(a, 0.4e1) * nu * (b - 0.1e1) * (b + 0.1e1) * Math.Pow(X[0], 0.4e1) / 0.2e1 + 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(X[0], 0.3e1) + (0.2e1 * b * b * nu - Math.Pow(b, 0.4e1) * nu / 0.2e1) * a * a * X[0] * X[0] + Math.Pow(b, 0.6e1) * X[0] / 0.2e1 - Math.Pow(b, 0.6e1) * nu / 0.2e1 + Math.Pow(b, 0.4e1) * nu / 0.2e1);
                    }
                } else {
                    // Transient
                    ExactSolution_psi = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t);
                    ExactSolution_uR = (X, t) => (1 - Math.Exp(-X[0] * X[0])) * Math.Sin(X[1]) * Math.Cos(t);
                    ExactSolution_uETA = (X, t) => (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1]) * Math.Cos(t);
                    ExactSolution_uXI = (X, t) => (0.2e1 * X[0] * X[0] * Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1]) + Math.Pow(a * a * X[0] * X[0] + b * b, -0.1e1 / 0.2e1) * (0.1e1 - Math.Exp(-X[0] * X[0])) * Math.Cos(X[1])) * Math.Cos(t);

                    if(Globals.activeMult == Globals.Multiplier.one) {
                        // Globals.Multiplier.one
                        ExactResidualConti = (X, t) => 0;
                        ExactResidualRmom = (X, t) => (0.2e1 * X[0] * Math.Cos(t) * (((-0.2e1 * Math.Pow(X[0], 0.3e1) + 0.2e1 * X[0]) * Math.Pow(Math.Cos(X[1]), 0.2e1) * Math.Cos(t) + Math.Sin(X[1]) * nu * (a * a * X[0] * X[0] + b * b)) * Math.Exp(-X[0] * X[0]) + 0.2e1 * X[0] * (X[0] * X[0] - 0.1e1 / 0.2e1) * Math.Cos(t) * Math.Pow(Math.Cos(X[1]), 0.2e1) * Math.Exp(-0.2e1 * X[0] * X[0]) - X[0] * Math.Cos(t) * Math.Pow(Math.Cos(X[1]), 0.2e1) - Math.Sin(X[1]) * nu * (a * a * X[0] * X[0] + b * b)) * b * a * Math.Sqrt(a * a * X[0] * X[0] + b * b) + (0.2e1 * (-((a * a + 0.2e1) * X[0] * X[0] + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(X[1]), 0.2e1) + Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * X[0] * X[0]) * X[0] * Math.Pow(Math.Cos(t), 0.2e1) - (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) * (-0.4e1 * a * a * nu * Math.Pow(X[0], 0.6e1) - 0.2e1 * a * a * Math.Pow(X[0], 0.5e1) + nu * (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(X[0], 0.4e1) - 0.2e1 * b * b * Math.Pow(X[0], 0.3e1) + 0.2e1 * nu * ((a * a + 0.4e1) * b * b + a * a / 0.2e1) * X[0] * X[0] + Math.Pow(b, 0.4e1) * nu - b * b * nu) * Math.Cos(t) + Math.Sin(X[1]) * X[0] * X[0] * Math.Sin(t) * Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1)) * Math.Exp(-X[0] * X[0]) - 0.2e1 * X[0] * Math.Pow(Math.Cos(t), 0.2e1) * (-(-0.4e1 * Math.Pow(X[0], 0.4e1) + (a * a + 0.4e1) * X[0] * X[0] + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(X[1]), 0.2e1) / 0.2e1 + Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * X[0] * X[0]) * Math.Exp(-0.2e1 * X[0] * X[0]) + b * b * X[0] * Math.Pow(Math.Cos(X[1]), 0.2e1) * (a * a * X[0] * X[0] + b * b - 0.1e1) * Math.Pow(Math.Cos(t), 0.2e1) + (a * a * X[0] * X[0] + b * b) * Math.Sin(X[1]) * nu * (Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.4e1) + (0.2e1 * a * a * b * b + a * a) * X[0] * X[0] + Math.Pow(b, 0.4e1) - b * b) * Math.Cos(t) - Math.Sin(X[1]) * X[0] * X[0] * Math.Sin(t) * Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1)) * Math.Pow(a * a * X[0] * X[0] + b * b, -0.2e1) * Math.Pow(X[0], -0.2e1);
                        ExactResidualETAmom = (X, t) => -(((-0.2e1 * Math.Sin(X[1]) * b * b * X[0] * (a * a * X[0] * X[0] + b * b) * Math.Pow(Math.Cos(t), 0.2e1) + (-0.4e1 * Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.8e1) + (Math.Pow(a, 0.6e1) + 0.4e1 * Math.Pow(a, 0.4e1) - 0.8e1 * a * a * b * b) * Math.Pow(X[0], 0.6e1) + (-0.4e1 * Math.Pow(b, 0.4e1) + (0.3e1 * Math.Pow(a, 0.4e1) + 0.8e1 * a * a) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(X[0], 0.4e1) + ((0.3e1 * a * a + 0.4e1) * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * X[0] * X[0] + Math.Pow(b, 0.6e1)) * nu * Math.Cos(t) - X[0] * X[0] * Math.Sin(t) * Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1)) * Math.Exp(-X[0] * X[0]) + Math.Sin(X[1]) * b * b * X[0] * (a * a * X[0] * X[0] + b * b) * Math.Pow(Math.Cos(t), 0.2e1) * Math.Exp(-0.2e1 * X[0] * X[0]) + Math.Sin(X[1]) * b * b * X[0] * (a * a * X[0] * X[0] + b * b) * Math.Pow(Math.Cos(t), 0.2e1) - (Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.6e1) + (0.3e1 * Math.Pow(a, 0.4e1) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(X[0], 0.4e1) + (0.3e1 * a * a * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * X[0] * X[0] + Math.Pow(b, 0.6e1)) * nu * Math.Cos(t) + X[0] * X[0] * Math.Sin(t) * Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1)) * Math.Sqrt(a * a * X[0] * X[0] + b * b) - 0.2e1 * X[0] * nu * ((-0.4e1 * a * a * Math.Pow(X[0], 0.6e1) + (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(X[0], 0.4e1) + ((0.2e1 * a * a + 0.8e1) * b * b + a * a) * X[0] * X[0] + Math.Pow(b, 0.4e1) - b * b) * Math.Exp(-X[0] * X[0]) - Math.Pow(a, 0.4e1) * Math.Pow(X[0], 0.4e1) + (-0.2e1 * a * a * b * b - a * a) * X[0] * X[0] - Math.Pow(b, 0.4e1) + b * b) * Math.Cos(t) * b * a) * Math.Cos(X[1]) * Math.Pow(a * a * X[0] * X[0] + b * b, -0.5e1 / 0.2e1) * Math.Pow(X[0], -0.2e1);
                        ExactResidualXImom = (X, t) => 0.2e1 * (X[0] * X[0] * ((-0.2e1 * Math.Sin(X[1]) * (a * a * X[0] * X[0] + b * b) * Math.Cos(t) - 0.2e1 * X[0] * nu * (a * a * X[0] * X[0] + a * a / 0.2e1 + b * b)) * Math.Exp(-X[0] * X[0]) + Math.Sin(X[1]) * (a * a * X[0] * X[0] + b * b) * Math.Cos(t) * Math.Exp(-0.2e1 * X[0] * X[0]) + Math.Sin(X[1]) * (a * a * X[0] * X[0] + b * b) * Math.Cos(t) + a * a * X[0] * nu) * Math.Cos(t) * b * a * Math.Sqrt(a * a * X[0] * X[0] + b * b) + (-0.2e1 * Math.Sin(X[1]) * Math.Pow(X[0], 0.3e1) * (X[0] - 0.1e1) * (X[0] + 0.1e1) * (a * a * X[0] * X[0] + b * b) * (a * a * X[0] * X[0] + a * a + b * b) * Math.Pow(Math.Cos(t), 0.2e1) + (-0.4e1 * Math.Pow(a, 0.4e1) * nu * Math.Pow(X[0], 0.10e2) + a * a * nu * (Math.Pow(a, 0.4e1) + 0.10e2 * a * a - 0.8e1 * b * b) * Math.Pow(X[0], 0.8e1) - Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.7e1) / 0.2e1 - (Math.Pow(a, 0.6e1) + (-0.6e1 * b * b + 0.2e1) * Math.Pow(a, 0.4e1) - 0.48e2 * a * a * b * b + 0.8e1 * Math.Pow(b, 0.4e1)) * nu * Math.Pow(X[0], 0.6e1) / 0.2e1 - 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(X[0], 0.5e1) - nu * ((b * b - 0.1e1) * Math.Pow(a, 0.4e1) + (-0.6e1 * Math.Pow(b, 0.4e1) + 0.4e1 * b * b) * a * a - 0.28e2 * Math.Pow(b, 0.4e1)) * Math.Pow(X[0], 0.4e1) / 0.2e1 - 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(X[0], 0.3e1) + ((b * b - 0.4e1) * a * a + 0.2e1 * Math.Pow(b, 0.4e1) - 0.10e2 * b * b) * nu * b * b * X[0] * X[0] / 0.2e1 - Math.Pow(b, 0.6e1) * X[0] / 0.2e1 + Math.Pow(b, 0.6e1) * nu / 0.2e1 - Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Cos(t) - Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) * Math.Sin(t) * X[0] * X[0] * (X[0] * X[0] - 0.1e1 / 0.2e1)) * Math.Exp(-X[0] * X[0]) - Math.Sin(X[1]) * Math.Pow(X[0], 0.3e1) * Math.Pow(Math.Cos(t), 0.2e1) * (a * a + 0.2e1 * b * b) * (a * a * X[0] * X[0] + b * b) * Math.Exp(-0.2e1 * X[0] * X[0]) - a * a * Math.Sin(X[1]) * Math.Pow(X[0], 0.3e1) * (a * a * X[0] * X[0] + b * b) * Math.Pow(Math.Cos(t), 0.2e1) + (Math.Pow(a, 0.6e1) * Math.Pow(X[0], 0.7e1) / 0.2e1 + Math.Pow(a, 0.6e1) * nu * Math.Pow(X[0], 0.6e1) / 0.2e1 + 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(X[0], 0.5e1) + Math.Pow(a, 0.4e1) * nu * (b - 0.1e1) * (b + 0.1e1) * Math.Pow(X[0], 0.4e1) / 0.2e1 + 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(X[0], 0.3e1) + (-Math.Pow(b, 0.4e1) * nu / 0.2e1 + 0.2e1 * b * b * nu) * a * a * X[0] * X[0] + Math.Pow(b, 0.6e1) * X[0] / 0.2e1 - Math.Pow(b, 0.6e1) * nu / 0.2e1 + Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Cos(t) - X[0] * X[0] * Math.Sin(t) * Math.Pow(a * a * X[0] * X[0] + b * b, 0.2e1) / 0.2e1) * Math.Cos(X[1]) * Math.Pow(a * a * X[0] * X[0] + b * b, -0.5e1 / 0.2e1) * Math.Pow(X[0], -0.2e1);
                    } else {
                        // Globals.Multiplier.Bsq
                        ExactResidualConti = (X, t) => 0;
                        ExactResidualRmom = (X, t) => (2 * b * nu * X[0] * (a * Math.Cos(t) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) - (b * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Sin(X[1]) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) / (X[0] * Math.Pow((a * a * X[0] * X[0] + b * b), (1 / 2))))) * Math.Pow((a * a * X[0] * X[0] + b * b), -(3 / 2)) - nu * Math.Cos(t) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) - (X[0] * X[0] * X[0] * Math.Pow((a * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) - (b * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((X[0] * Math.Pow((a * a * X[0] * X[0] + b * b), (1 / 2))), -1)), 2)) * Math.Pow((a * a * X[0] * X[0] + b * b), -2) + (X[0] * X[0] * Math.Sin(t) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) - (nu * Math.Cos(t) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) + (2 * X[0] * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Sin(X[1])) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) - (2 * X[0] * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Sin(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) / (a * a * X[0] * X[0] + b * b) - (X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) + (4 * nu * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Sin(X[1]) * (X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1);
                        ExactResidualETAmom = (X, t) => (X[0] * X[0] * Math.Cos(X[1]) * Math.Sin(t) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) - nu * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) - (2 * X[0] * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) + (4 * nu * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1]) * (X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) + (a * a * X[0] * X[0] * X[0] * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -2) + (X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) - (a * a * nu * X[0] * X[0] * Math.Cos(t) * Math.Cos(X[1]) * (a * a * X[0] * X[0] + 2 * b * b) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -3) + (2 * a * b * nu * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1]) * (b * b * Math.Exp(X[0] * X[0]) - b * b * b * b * Math.Exp(X[0] * X[0]) - b * b + b * b * b * b + a * a * X[0] * X[0] + 4 * a * a * X[0] * X[0] * X[0] * X[0] - 4 * a * a * X[0] * X[0] * X[0] * X[0] * X[0] * X[0] + a * a * a * a * X[0] * X[0] * X[0] * X[0] + 8 * b * b * X[0] * X[0] - 4 * b * b * X[0] * X[0] * X[0] * X[0] - a * a * X[0] * X[0] * Math.Exp(X[0] * X[0]) - a * a * a * a * X[0] * X[0] * X[0] * X[0] * Math.Exp(X[0] * X[0]) + 2 * a * a * b * b * X[0] * X[0] - 2 * a * a * b * b * X[0] * X[0] * Math.Exp(X[0] * X[0]))) * Math.Pow((a * a * X[0] * X[0] + b * b), -(7 / 2));
                        ExactResidualXImom = (X, t) => (X[0] * X[0] * Math.Sin(t) * ((Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(1 / 2)) - (2 * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(X[1])) * Math.Pow((a * a * X[0] * X[0] + b * b), -(1 / 2)))) * Math.Pow((a * a * X[0] * X[0] + b * b), -1) - (2 * b * nu * X[0] * X[0] * ((a * a * a * X[0] * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(3 / 2)) - (b * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((X[0] * X[0]), -1) + (2 * a * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1])) * Math.Pow((a * a * X[0] * X[0] + b * b), -(1 / 2)))) * Math.Pow((a * a * X[0] * X[0] + b * b), -(3 / 2)) - (X[0] * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(1 / 2)) + (nu * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1]) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(1 / 2)) - (nu * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1]) * (12 * b * b * b * b + 2 * a * a * b * b - a * a * a * a * X[0] * X[0] + 2 * a * a * a * a * X[0] * X[0] * X[0] * X[0] - 20 * a * a * a * a * X[0] * X[0] * X[0] * X[0] * X[0] * X[0] + 8 * a * a * a * a * X[0] * X[0] * X[0] * X[0] * X[0] * X[0] * X[0] * X[0] - 28 * b * b * b * b * X[0] * X[0] + 8 * b * b * b * b * X[0] * X[0] * X[0] * X[0] - 2 * a * a * b * b * Math.Exp(X[0] * X[0]) + a * a * a * a * X[0] * X[0] * Math.Exp(X[0] * X[0]) + 8 * a * a * b * b * X[0] * X[0] - 48 * a * a * b * b * X[0] * X[0] * X[0] * X[0] + 16 * a * a * b * b * X[0] * X[0] * X[0] * X[0] * X[0] * X[0])) * Math.Pow((a * a * X[0] * X[0] + b * b), -(7 / 2)) - (X[0] * Math.Exp(-2 * X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(3 / 2)) + (b * b * nu * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(X[1]) * (2 * a * a * X[0] * X[0] + b * b) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(7 / 2)) + (2 * a * b * X[0] * X[0] * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) * (Math.Exp(-X[0] * X[0]) - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -2) + (X[0] * X[0] * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) * (a * a * Math.Exp(X[0] * X[0]) - a * a - 6 * b * b - 4 * a * a * X[0] * X[0] + 4 * a * a * X[0] * X[0] * X[0] * X[0] + 4 * b * b * X[0] * X[0])) * Math.Pow((a * a * X[0] * X[0] + b * b), -(5 / 2)) - (b * b * X[0] * Math.Exp(-X[0] * X[0]) * Math.Cos(t) * Math.Cos(t) * Math.Cos(X[1]) * Math.Sin(X[1]) * (Math.Exp(-X[0] * X[0]) - 1) * (Math.Exp(X[0] * X[0]) + 2 * X[0] * X[0] - 1)) * Math.Pow((a * a * X[0] * X[0] + b * b), -(5 / 2));
                    }
                }

                Func<double[], double, double>[] ExactResiduals;
                ExactResiduals = new[] { ExactResidualRmom, ExactResidualXImom, ExactResidualETAmom, ExactResidualConti };

                Debug.Assert(ExactResiduals.Length == ResidualFields.Count);
                Debug.Assert(ExactResiduals.Length == ErrorFields.Count);
                double[] L2errorS = new double[ExactResiduals.Length];

                for(int iVar = 0; iVar < ExactResiduals.Length; iVar++) {
                    if(ExactResiduals[iVar] != null) {

                        double L2resinorm = ResidualFields[iVar].L2Norm();
                        Console.WriteLine("L2 norm '{0}': {1}", ResidualFields[iVar].Identification, L2resinorm);
                        L2errorS[iVar] = L2resinorm;
                    }
                }

                this.error_MomR = L2errorS[0];
                this.error_MomZ = L2errorS[1];
                this.error_MomETA = L2errorS[2];
                this.error_Conti = L2errorS[3];

                // exact solution
                // ==============
                if(this.Control.steady == true) {
                    if(!Control.ExactResidual) {
                        urExact.ProjectField(Globals.DirichletValue_uR);
                        uetaExact.ProjectField(Globals.DirichletValue_uEta);
                        uxiExact.ProjectField(Globals.DirichletValue_uXi);
                        psiExact.ProjectField(Globals.DirichletValue_psi);
                    }
                } else {
                    urExact.ProjectField(ExactSolution_uR.Vectorize(time));
                    uetaExact.ProjectField(ExactSolution_uETA.Vectorize(time));
                    uxiExact.ProjectField(ExactSolution_uXI.Vectorize(time));
                    psiExact.ProjectField(ExactSolution_psi.Vectorize(time));
                }
                urError.Clear();
                uetaError.Clear();
                uxiError.Clear();
                psiError.Clear();


                if(Control.steady == true) {
                    urError.ProjectField(Globals.DirichletValue_uR.Vectorize());
                    uetaError.ProjectField(Globals.DirichletValue_uEta.Vectorize());
                    uxiError.ProjectField(Globals.DirichletValue_uXi.Vectorize());
                    psiError.ProjectField(Globals.DirichletValue_psi.Vectorize());

                } else {
                    urError.ProjectField(ExactSolution_uR.Vectorize(time));
                    uetaError.ProjectField(ExactSolution_uETA.Vectorize(time));
                    uxiError.ProjectField(ExactSolution_uXI.Vectorize(time));
                    psiError.ProjectField(ExactSolution_psi.Vectorize(time));

                }

                urError.Acc(-1, ur);
                uetaError.Acc(-1, ueta);
                uxiError.Acc(-1, uxi);
                psiError.Acc(-1, Pressure);

                double psiErrorMean = psiError.GetMeanValueTotal(null);
                psiError.AccConstant(-psiErrorMean);

                //// Hier werden die beim Restart im jedem Zeitschritt neu angelegt. WIESO!?!?!?!
                //m_IOFields.Add(urError);   // In plots
                //m_IOFields.Add(uetaError);
                //m_IOFields.Add(uxiError);
                //m_IOFields.Add(psiError);



                urErrorL2 = urError.L2Norm();
                uetaErrorL2 = uetaError.L2Norm();
                uxiErrorL2 = uxiError.L2Norm();
                psiErrorL2 = psiError.L2Norm();


                base.QueryHandler.ValueQuery("urErrorL2", urErrorL2, true);
                base.QueryHandler.ValueQuery("uetaErrorL2", uetaErrorL2, true);
                base.QueryHandler.ValueQuery("uxiErrorL2", uxiErrorL2, true);
                base.QueryHandler.ValueQuery("psiErrorL2", psiErrorL2, true);
                base.QueryHandler.ValueQuery("condNum", cond_Full_012Vars, true);


                Console.WriteLine("L2 Error of ur: " + urErrorL2);
                Console.WriteLine("L2 Error of ueta: " + uetaErrorL2);
                Console.WriteLine("L2 Error of uxi: " + uxiErrorL2);
                Console.WriteLine("L2 Error of pressure: " + psiErrorL2);


                ///
                /// Compute LxErrors:
                ///

                ScalarFunction scalarFuncUR;
                ScalarFunction scalarFuncUETA;
                ScalarFunction scalarFuncUXI;
                ScalarFunction scalarFuncPressure;



                if(Control.steady == true) {
                    scalarFuncUR = Globals.DirichletValue_uR.Vectorize();
                    scalarFuncUETA = Globals.DirichletValue_uEta.Vectorize();
                    scalarFuncUXI = Globals.DirichletValue_uXi.Vectorize();
                    scalarFuncPressure = Globals.DirichletValue_psi.Vectorize();
                    urErrorLx = ur.LxError(scalarFuncUR, (X, a, b) => Globals.f_function_(X[0]) * (b - a).Pow2(), (new CellQuadratureScheme()).SaveCompile(ur.GridDat, 2 * ur.Basis.Degree)).Sqrt();
                    uetaErrorLx = ueta.LxError(scalarFuncUETA, (X, a, b) => Globals.f_function_(X[0]) * (b - a).Pow2(), (new CellQuadratureScheme()).SaveCompile(ueta.GridDat, 2 * ueta.Basis.Degree)).Sqrt();
                    uxiErrorLx = uxi.LxError(scalarFuncUXI, (X, a, b) => Globals.f_function_(X[0]) * (b - a).Pow2(), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * uxi.Basis.Degree)).Sqrt();
                    double LxIntA = Pressure.LxError(scalarFuncPressure, (X, a, b) => Globals.f_function_(X[0]) * (b - a), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * ur.Basis.Degree));
                    double LxIntB = Pressure.LxError(scalarFuncPressure, (X, a, b) => Globals.f_function_(X[0]), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * ur.Basis.Degree));
                    double LxConst = Math.Abs(LxIntA) / LxIntB;
                    psiErrorLx = Pressure.LxError(scalarFuncPressure, (X, a, b) => Globals.f_function_(X[0]) * ((b - a) - LxConst).Pow2(), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * ur.Basis.Degree)).Sqrt();

                } else {
                    scalarFuncUR = ExactSolution_uR.Vectorize(time);
                    scalarFuncUETA = ExactSolution_uETA.Vectorize(time);
                    scalarFuncUXI = ExactSolution_uXI.Vectorize(time);
                    scalarFuncPressure = ExactSolution_psi.Vectorize(time);
                    urErrorLx = ur.LxError(scalarFuncUR, (X, a, b) => Globals.f_function_(X[0]) * (b - a).Pow2(), (new CellQuadratureScheme()).SaveCompile(ur.GridDat, 2 * ur.Basis.Degree)).Sqrt();
                    uetaErrorLx = ueta.LxError(scalarFuncUETA, (X, a, b) => Globals.f_function_(X[0]) * (b - a).Pow2(), (new CellQuadratureScheme()).SaveCompile(ueta.GridDat, 2 * ueta.Basis.Degree)).Sqrt();
                    uxiErrorLx = uxi.LxError(scalarFuncUXI, (X, a, b) => Globals.f_function_(X[0]) * (b - a).Pow2(), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * uxi.Basis.Degree)).Sqrt();
                    double LxIntA = Pressure.LxError(scalarFuncPressure, (X, a, b) => Globals.f_function_(X[0]) * (b - a), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * ur.Basis.Degree));
                    double LxIntB = Pressure.LxError(scalarFuncPressure, (X, a, b) => Globals.f_function_(X[0]), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * ur.Basis.Degree));
                    double LxConst = Math.Abs(LxIntA) / LxIntB;
                    psiErrorLx = Pressure.LxError(scalarFuncPressure, (X, a, b) => Globals.f_function_(X[0]) * ((b - a) - LxConst).Pow2(), (new CellQuadratureScheme()).SaveCompile(uxi.GridDat, 2 * ur.Basis.Degree)).Sqrt();
                }


                base.QueryHandler.ValueQuery("urErrorLx", urErrorLx, true);
                base.QueryHandler.ValueQuery("uetaErrorLx", uetaErrorLx, true);
                base.QueryHandler.ValueQuery("uxiErrorLx", uxiErrorLx, true);
                base.QueryHandler.ValueQuery("presErrorLx", psiErrorLx, true);
                base.QueryHandler.ValueQuery("condNum", cond_Full_012Vars, true);


                Console.WriteLine("LxError ur: " + urErrorLx);
                Console.WriteLine("LxError ueta: " + (uetaErrorLx));
                Console.WriteLine("LxError uxi: " + uxiErrorLx);
                Console.WriteLine("LxError pressure: " + psiErrorLx);

                ///Writing results to text file
                ///    
                string dt_ = base.GetTimestep().ToString();

                Guid sessionID = this.CurrentSessionInfo.ID;
                TextWriter file = base.DatabaseDriver.FsDriver.GetNewLog("L2error_" + dt_, sessionID);
                file.WriteLine("=========================================================================");
                file.WriteLine("DG degree: " + ur.Basis.Degree + " | " + "Number of Cells r: " + Control.Resolution_R + " | " + "Number of Cells xi: " + Control.Resolution_Xi
                                    + " | " + "Compute the case: " + Control.ControlFileText);
                file.WriteLine("L2 Error ur: " + Math.Abs(urErrorL2).ToString());
                file.WriteLine("L2 Error ueta: " + Math.Abs(uetaErrorL2).ToString());
                file.WriteLine("L2 Error uxi: " + Math.Abs(uxiErrorL2).ToString());
                file.WriteLine("L2 Error pressure: " + Math.Abs(psiErrorL2).ToString());
                file.Flush();


                file = base.DatabaseDriver.FsDriver.GetNewLog("L2errorTimesteps_" + dt_, sessionID);
                file.WriteLine("=========================================================================");
                file.WriteLine("DG degree: " + ur.Basis.Degree + " | " + "Number of Cells r: " + Control.Resolution_R + " | " + "Number of Cells xi: " + Control.Resolution_Xi
                                    + " | " + "Compute the case: " + Control.ControlFileText);
                file.WriteLine("dt= " + Control.dtFixed);
                file.WriteLine("L2 Error ur: " + Math.Abs(urErrorL2).ToString());
                file.WriteLine("L2 Error ueta: " + Math.Abs(uetaErrorL2).ToString());
                file.WriteLine("L2 Error uxi: " + Math.Abs(uxiErrorL2).ToString());
                file.WriteLine("L2 Error pressure: " + Math.Abs(psiErrorL2).ToString());
                file.Flush();


                file = base.DatabaseDriver.FsDriver.GetNewLog("L2errorTSconvergence_" + dt_, sessionID);
                file.WriteLine("=========================================================================");
                file.WriteLine("DG degree: " + ur.Basis.Degree + " | " + "Number of Cells r: " + Control.Resolution_R + " | " + "Number of Cells xi:" + Control.Resolution_Xi
                                + " | " + "Compute the case: " + Control.ControlFileText);
                file.WriteLine("dt= " + Control.dtFixed);
                file.WriteLine("L2 Error ur: " + Math.Abs(urErrorL2).ToString());
                file.WriteLine("L2 Error ueta: " + Math.Abs(uetaErrorL2).ToString());
                file.WriteLine("L2 Error uxi: " + Math.Abs(uxiErrorL2).ToString());
                file.WriteLine("L2 Error pressure: " + Math.Abs(psiErrorL2).ToString());
                file.Flush();


                // write formatted txt file for pgf plotting in latex
                if(timestepNo[0] == Control.NoOfTimesteps) {
                    Console.WriteLine("timestep number {0} written", timestepNo[0]);
                    file = base.DatabaseDriver.FsDriver.GetNewLog("L2errorTSconvergenceFormat" + dt_, sessionID);
                    file.WriteLine(Control.dtFixed + "\t" + Math.Abs(psiErrorL2).ToString());
                    file.Flush();
                }


                file = base.DatabaseDriver.FsDriver.GetNewLog("LxError_" + dt_, sessionID);
                file.WriteLine("=========================================================================");
                file.WriteLine("DG degree: " + ur.Basis.Degree + " | " + "Number of Cells r: " + Control.Resolution_R + " | " + "Number of Cells xi:" + Control.Resolution_Xi
                                + " | " + "Compute the case: " + Control.ControlFileText);
                file.WriteLine("Lx Error ur: " + Math.Abs(urErrorLx).ToString());
                file.WriteLine("Lx Error ueta: " + Math.Abs(uetaErrorLx).ToString());
                file.WriteLine("Lx Error uxi: " + Math.Abs(uxiErrorLx).ToString());
                file.WriteLine("Lx Error pressure: " + Math.Abs(psiErrorLx).ToString());
                file.Flush();


                string fileNameUR = "LxErrorURdeg_" + ur.Basis.Degree + "_" + dt_;
                file = base.DatabaseDriver.FsDriver.GetNewLog(fileNameUR, sessionID);
                file.WriteLine(Math.Abs(urErrorLx).ToString());
                file.Flush();

                string fileNameUXI = "LxErrorUXIdeg_" + ur.Basis.Degree + "_" + dt_;
                file = base.DatabaseDriver.FsDriver.GetNewLog(fileNameUXI, sessionID);
                file.WriteLine(Math.Abs(uxiErrorLx).ToString());
                file.Flush();

                string fileNameUETA = "LxErrorUETAdeg_" + ur.Basis.Degree + "_" + dt_;
                file = base.DatabaseDriver.FsDriver.GetNewLog(fileNameUETA, sessionID);
                file.WriteLine(Math.Abs(uetaErrorLx).ToString());
                file.Flush();

                string fileNamePressure = "LxErrorPSIdeg_" + Pressure.Basis.Degree + "_" + dt_;
                file = base.DatabaseDriver.FsDriver.GetNewLog(fileNamePressure, sessionID);
                file.WriteLine(Math.Abs(psiErrorLx).ToString());
                file.Flush();
            }
        }

        /// <summary>
        /// used by <see cref="ConditionNumberScalingTest"/> and <see cref="ConditionNumberScalingTest.Perform"/> to perform fully automatic condition number test
        /// </summary>
        public override IDictionary<string, double> OperatorAnalysis(OperatorAnalysisConfig config) {
            var VarGroups = new int[][] { new[] { 0 }, new[] { 1 }, new[] { 2 }, new[] { 0, 1, 2 }, new[] { 0, 1, 2, 3 } };

            bool plotStencilCondNumViz = true;

            string time_dep;
            if(this.Control.steady == true) {
                time_dep = "Steady_";
            } else {
                time_dep = "Transient_";
            }
            string nameOf_Stencil = "Stencil" + time_dep + "_R0fix=_" + Control.R0fixOn.ToString() + "_r_min_" + this.Control.rMin.ToString() + "_" + this.Control.Resolution_R.ToString() + "_x_" +
                                  this.Control.Resolution_R.ToString() + "_x_" + this.ueta.Basis.Degree.ToString() + "_dt_" + base.GetTimestep().ToString();

            if(Control.rMin < 10e-6) {
                return m_Splitting_Timestepper.OperatorAnalysisAk(VarGroups,
                    calclulateStencils: config.CalculateStencilConditionNumbers, calculateGlobals: config.CalculateGlobalConditionNumbers, plotStencilCondNumViz: plotStencilCondNumViz, nameOfStencil: nameOf_Stencil);
            } else {
                return m_Splitting_Timestepper.OperatorAnalysisAk(VarGroups,
                    calclulateStencils: config.CalculateStencilConditionNumbers, calculateGlobals: config.CalculateGlobalConditionNumbers, plotStencilCondNumViz: plotStencilCondNumViz, nameOfStencil: nameOf_Stencil);
            }
        }
    }

    class TransformToKartesian {

        public TransformToKartesian(int degree, int zRes, int phiRes) {
            Cylinder = Grid3D.Ogrid(0.5, 0.9, phiRes, phiRes / 2, GenericBlas.Linspace(0, 2 * Math.PI, zRes));
            gdat = new GridData(Cylinder);

            Basis b = new Basis(gdat, degree);
            m_VelocityVector = new VectorField<SinglePhaseField>(
                new SinglePhaseField(b, VariableNames.VelocityX), new SinglePhaseField(b, VariableNames.VelocityY), new SinglePhaseField(b, VariableNames.VelocityZ));
            m_Pressure = new SinglePhaseField(b, VariableNames.Pressure);

        }

        GridCommons Cylinder;
        public GridData gdat {
            get;
            private set;
        }

        VectorField<SinglePhaseField> m_VelocityVector;

        public VectorField<SinglePhaseField> VelocityVector {
            get {
                return m_VelocityVector;
            }
        }

        SinglePhaseField m_Pressure;

        public SinglePhaseField Pressure {
            get {
                return m_Pressure;
            }
        }

        public SinglePhaseField VelocityX {
            get {
                return m_VelocityVector[0];
            }
        }
        public SinglePhaseField VelocityY {
            get {
                return m_VelocityVector[1];
            }
        }
        public SinglePhaseField VelocityZ {
            get {
                return m_VelocityVector[2];
            }
        }

        public void SetSolution(SinglePhaseField uR_Helix, SinglePhaseField uXi_Helix, SinglePhaseField uEta_Helix, SinglePhaseField Pressure_Helix) {
            GridData gdat_Helix = (GridData)(uR_Helix.Basis.GridDat);


            VelocityVector.Clear();
            Pressure.Clear();

            var p = new HelixData(uR_Helix, uXi_Helix, uEta_Helix, Pressure_Helix);

            VelocityX.ProjectField(p.UX);
            VelocityY.ProjectField(p.UY);
            VelocityZ.ProjectField(p.UZ);

        }

        class HelixData {

            public HelixData(SinglePhaseField __uR_Helix, SinglePhaseField __uXi_Helix, SinglePhaseField __uEta_Helix, SinglePhaseField __Pressure_Helix) {
                uR_Helix = __uR_Helix;
                uXi_Helix = __uXi_Helix;
                uEta_Helix = __uEta_Helix;
                Pressure_Helix = __Pressure_Helix;
                eval = new FieldEvaluation(gdat_Helix);
            }

            SinglePhaseField uR_Helix;
            SinglePhaseField uXi_Helix;
            SinglePhaseField uEta_Helix;
            SinglePhaseField Pressure_Helix;

            GridData gdat_Helix {
                get {
                    GridData g = (GridData)(uR_Helix.Basis.GridDat);
                    Debug.Assert(object.ReferenceEquals(g, uXi_Helix.Basis.GridDat));
                    Debug.Assert(object.ReferenceEquals(g, uEta_Helix.Basis.GridDat));
                    Debug.Assert(object.ReferenceEquals(g, Pressure_Helix.Basis.GridDat));
                    return g;
                }

            }

            FieldEvaluation eval;


            public const double G = 1;
            public const double a = 1;
            public const double b = 1;
            public const double K = 1;
            public const double Z = 1;



            static double R(double x, double y) {
                return Math.Sqrt(x * x + y * y);
            }

            public static double XI(double x, double y, double z) {
                double phi = Math.Atan2(y, x);
                if(phi < 0)
                    phi = 2 * Math.PI - Math.Abs(phi);
                return (a * z + b * phi);
            }


            static public double XI_(double[] X) {
                double x = X[0];
                double y = X[1];
                double z = X[2];

                double xi = XI(x, y, z);

                return xi;
            }


            public void UR_UPHI(MultidimensionalArray r_xi, MultidimensionalArray ur_uphi) {
                int K = r_xi.GetLength(0);

                MultidimensionalArray ur_uxi_ueta = MultidimensionalArray.Create(K, 3);
                int Unalloc = eval.Evaluate(1.0, new DGField[] { uR_Helix, uXi_Helix, uEta_Helix }, r_xi, 0.0, ur_uxi_ueta);

                for(int k = 0; k < K; k++) { // loop over all nodes
                    double ur = ur_uxi_ueta[k, 0];
                    double ueta = ur_uxi_ueta[k, 1];
                    double uxi = ur_uxi_ueta[k, 2];
                    double r = r_xi[k, 0];

                    ur_uphi[k, 0] = ur;
                    ur_uphi[k, 1] = ((a * r * ueta + b * uxi) / Math.Sqrt(a * a * r * r + b * b));

                }
            }

            public void UZ_r_xi(MultidimensionalArray r_xi, MultidimensionalArray uz) {
                int K = r_xi.GetLength(0);

                MultidimensionalArray ueta_uxi = MultidimensionalArray.Create(K, 2);
                int Unalloc = eval.Evaluate(1.0, new DGField[] { uEta_Helix, uXi_Helix }, r_xi, 0.0, ueta_uxi);

                for(int k = 0; k < K; k++) { // loop over all nodes

                    double ueta = ueta_uxi[k, 0];
                    double uxi = ueta_uxi[k, 1];
                    double r = r_xi[k, 0];

                    uz[k] = ((a * r * uxi - b * ueta) / Math.Sqrt(a * a * r * r + b * b));
                }
            }

            public void UX(MultidimensionalArray input, MultidimensionalArray output) {
                int K = input.GetLength(0); // number of nodes
                var r_xi = Get_r_xi(input);
                MultidimensionalArray ur_uphi = MultidimensionalArray.Create(K, 2);
                UR_UPHI(r_xi, ur_uphi);

                for(int k = 0; k < K; k++) {
                    double x = input[k, 0];
                    double y = input[k, 1];
                    double ur = ur_uphi[k, 0];
                    double uphi = ur_uphi[k, 1];
                    double phi = Math.Atan2(y, x);
                    output[k] = Math.Cos(phi) * ur - uphi * Math.Sin(phi);
                }
            }

            public void UY(MultidimensionalArray input, MultidimensionalArray output) {
                int K = input.GetLength(0); // number of nodes
                var r_xi = Get_r_xi(input);
                MultidimensionalArray ur_uphi = MultidimensionalArray.Create(K, 2);
                UR_UPHI(r_xi, ur_uphi);

                for(int k = 0; k < K; k++) {
                    double x = input[k, 0];
                    double y = input[k, 1];
                    double ur = ur_uphi[k, 0];
                    double uphi = ur_uphi[k, 1];
                    double phi = Math.Atan2(y, x);
                    output[k] = Math.Sin(phi) * ur + uphi * Math.Cos(phi);
                }
            }


            public void UZ(MultidimensionalArray input, MultidimensionalArray output) {
                Debug.Assert(output.GetLength(0) == input.GetLength(0));
                var r_xi = Get_r_xi(input);
                UZ_r_xi(r_xi, output);
            }

            private static MultidimensionalArray Get_r_xi(MultidimensionalArray input) {
                int K = input.GetLength(0); // number of nodes
                Debug.Assert(input.GetLength(1) == 3);

                MultidimensionalArray r_xi = MultidimensionalArray.Create(K, 2);
                for(int k = 0; k < K; k++) {
                    double x = input[k, 0];
                    double y = input[k, 1];
                    double z = input[k, 2];
                    double r = R(x, y);
                    double xi = XI(x, y, z);
                    r_xi[k, 0] = r;
                    r_xi[k, 1] = xi;
                }

                return r_xi;
            }
        }


    }
}
