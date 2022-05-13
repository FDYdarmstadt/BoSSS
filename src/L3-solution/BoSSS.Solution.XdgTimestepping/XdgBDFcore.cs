using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Timestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.XdgTimestepping {
    
    /*
    class XdgTScore {

        /// <summary>
        /// How the interface motion should be integrated
        /// </summary>
        public LevelSetHandling Config_LevelSetHandling;

        public MassMatrixShapeandDependence Config_MassMatrixShapeandDependence;

        /// <summary>
        /// Quadrature order on cut cells.
        /// </summary>
        public int Config_CutCellQuadratureOrder {
            get;
            protected set;
        }

        /// <summary>
        /// Callback routine to update the operator matrix.
        /// </summary>
        public DelComputeOperatorMatrix ComputeOperatorMatrix {
            get;
            protected set;
        }

        /// <summary>
        /// stack of solution vectors 
        /// </summary>
        public CoordinateVector[] m_Stack_u;

        /// <summary>
        /// DG/XDG coefficient mapping for the test- and trial-space;
        /// After Solver call, the updated solution; on entry, the initial guess
        /// </summary>
        public CoordinateMapping CurrentStateMapping {
            get {
                return m_Stack_u[0].Mapping;
            }
        }

        public LevelSetTracker m_LsTrk;

        /// <summary>
        /// Agglomerator for the currently set level-set position . 
        /// </summary>
        protected MultiphaseCellAgglomerator m_CurrentAgglomeration;

        /// <summary>
        /// Agglomerated length scales, input for <see cref="ComputeOperatorMatrix"/>.
        /// </summary>
        internal protected Dictionary<SpeciesId, MultidimensionalArray> GetAgglomeratedLengthScales() {
            if(m_CurrentAgglomeration != null) {
                //
                // agglomerated length scales are available from 
                //
                return m_CurrentAgglomeration.CellLengthScales;
            } else {
                //
                // 'Notlösung' -- no actual agglomeration available - use length scales form a temporary agglomerator.
                //
                var agg = this.m_LsTrk.GetAgglomerator(this.Config_SpeciesToCompute, this.Config_CutCellQuadratureOrder, this.Config_AgglomerationThreshold);
                return agg.CellLengthScales;
            }
        }

        /// <summary>
        /// Les spatial operateur 
        /// </summary>
        public virtual ISpatialOperator AbstractOperator {
            get;
            protected set;
        }

        /// <summary>
        /// Species to compute, must be a subset of <see cref="LevelSetTracker.SpeciesIdS"/>
        /// </summary>
        public SpeciesId[] Config_SpeciesToCompute {
            get;
            protected set;
        }

        /// <summary>
        /// As usual the threshold for cell agglomeration.
        /// </summary>
        public double Config_AgglomerationThreshold {
            get;
            set;
        }


    }

    */

    /*
    
    /// <summary>
    /// Numerical core of the XDG BDF method, not for direct user interaction
    /// </summary>
    class XdgBDFcore : XdgTScore {

        public BDFSchemeCoeffs Tsc;

        int m_IterationCounter = 0;

        int m_CoupledIterations = 0;

        double m_CurrentDt;


        BlockMsrMatrix m_PrecondMassMatrix;


        
        /// <summary>
        /// stack of mass matrices (matrices _without_ agglomeration)
        /// </summary>
        BlockMsrMatrix[] m_Stack_MassMatrix;


        /// <summary>
        /// Callback-routine  to update the linear resp. linearized system, 
        /// see <see cref="OperatorEvalOrLin"/> resp. <see cref="NonlinearSolver.m_AssembleMatrix"/>.
        /// </summary>
        /// <param name="argCurSt">Input, current state of solution.</param>
        /// <param name="System">Output.</param>
        /// <param name="Affine">Output.</param>
        /// <param name="PrecondMassMatrix">
        /// Mass matrix including agglomeration, without any scaling,
        /// required for block-precond.
        /// </param>
        /// <param name="Linearization">
        /// - true: assemble matrix and affine vector
        /// - false: evaluate operator (<paramref name="System"/> will be null)
        /// </param>
        /// <param name="abstractOperator">
        ///  the original operator that somehow produced the matrix; yes, this API is convoluted piece-of-shit
        /// </param>
        public void AssembleMatrixCallback(out BlockMsrMatrix System, out double[] Affine, out BlockMsrMatrix PrecondMassMatrix, DGField[] argCurSt, bool Linearization, out ISpatialOperator abstractOperator) {
            using (new FuncTrace()) {

                // copy data from 'argCurSt' to 'CurrentStateMapping', if necessary 
                // -----------------------------------------------------------
                DGField[] locCurSt = CurrentStateMapping.Fields.ToArray();
                if (locCurSt.Length != argCurSt.Length) {
                    throw new ApplicationException();
                }
                int NF = locCurSt.Length;
                for (int iF = 0; iF < NF; iF++) {
                    if (object.ReferenceEquals(locCurSt[iF], argCurSt[iF])) {
                        // nothing to do
                    } else {
                        locCurSt[iF].Clear();
                        locCurSt[iF].Acc(1.0, argCurSt[iF]);
                    }

                }


                // update of level-set
                // ----------------------

                bool updateAgglom = false;
                if ((this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once && m_IterationCounter == 0)
                    || (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && m_IterationCounter == 0)
                    || (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative && CoupledIteration)) {

                    m_CoupledIterations++;
                    if(this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative)
                        Console.WriteLine("Coupled Iteration {0}:", m_CoupledIterations);

                    MoveLevelSetAndRelatedStuff(locCurSt, m_CurrentPhystime, m_CurrentDt, IterUnderrelax);

                    // note that we need to update the agglomeration
                    updateAgglom = true;
                }

                if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {
                    if (m_IterationCounter == 0) {
                        if(m_CurrentAgglomeration != null)
                            throw new ApplicationException();
                        updateAgglom = true;
                    } else {
                        if (m_CurrentAgglomeration == null)
                            throw new ApplicationException();
                    }
                    // ensure, that, when splitting is used we update the agglomerator in the very first iteration.
                }


                // update agglomeration
                // --------------------

                if (updateAgglom || m_CurrentAgglomeration == null) {

                    if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                        || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled) {
                        // Agglomeration update in the case of splitting - agglomeration does **NOT** depend on previous time-steps

                        Debug.Assert(m_IterationCounter == 0);
                        m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder, this.Config_AgglomerationThreshold,
                            AgglomerateNewborn: false, AgglomerateDecased: false, ExceptionOnFailedAgglomeration: true);

                    } else {
                        // Agglomeration update in the case of a moving interface - consider previous time-steps
                        double[] oldAggTrsh = new double[Tsc.S];
                        ArrayTools.SetAll(oldAggTrsh, this.Config_AgglomerationThreshold);
                        Debug.Assert(m_Stack_MassMatrix.Where(mm => mm != null).Count() == Tsc.S);


                        m_CurrentAgglomeration = m_LsTrk.GetAgglomerator(base.Config_SpeciesToCompute, base.Config_CutCellQuadratureOrder,
                            __AgglomerationTreshold: base.Config_AgglomerationThreshold,
                            AgglomerateNewborn: (oldAggTrsh != null), AgglomerateDecased: (oldAggTrsh != null),
                            ExceptionOnFailedAgglomeration: true,
                            oldTs__AgglomerationTreshold: oldAggTrsh);
                    }


                    // update Multigrid-XDG basis
                    Debug.Assert(object.ReferenceEquals(base.MultigridBasis[0][0].DGBasis.GridDat, m_CurrentAgglomeration.Tracker.GridDat));
                    base.MultigridBasis.UpdateXdgAggregationBasis(m_CurrentAgglomeration);
                }

                // mass matrix update
                // ------------------

                if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity) {
                    // may occur e.g. if one runs the FSI solver as a pure single-phase solver,
                    // i.e. if the Level-Set is outside the domain.

                    PrecondMassMatrix = null;
                } else {
                    // checks:
                    Debug.Assert(this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                        || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent
                        || this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent,
                        "Something is not implemented here.");
                    if (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity
                        && Config_LevelSetHandling != LevelSetHandling.None)
                        throw new NotSupportedException("Illegal configuration;");

                    Debug.Assert((m_PrecondMassMatrix == null) == (m_Stack_MassMatrix[0] == null));
                    if ((this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsNonIdentity && m_PrecondMassMatrix == null) // compute mass matrix (only once in application lifetime)
                        || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeDependent && m_IterationCounter == 0) // compute mass matrix once per timestep
                        || (this.Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsTimeAndSolutionDependent && CoupledIteration) // re-compute mass matrix in every iteration
                        ) {
                        Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));

                        //MassMatrixFactory MassFact = new MassMatrixFactory(CurrentStateMapping.BasisS.ElementAtMax(b => b.Degree), m_CurrentAgglomeration);
                        MassMatrixFactory MassFact = m_LsTrk.GetXDGSpaceMetrics(Config_SpeciesToCompute, Config_CutCellQuadratureOrder, 1).MassMatrixFactory;

                        // matrix for preconditioner (Agglom required)
                        m_PrecondMassMatrix = MassFact.GetMassMatrix(CurrentStateMapping, false);
                        m_CurrentAgglomeration.ManipulateMatrixAndRHS(m_PrecondMassMatrix, default(double[]), CurrentStateMapping, CurrentStateMapping);

                        // mass matrix for time derivative
                        m_Stack_MassMatrix[0] = new BlockMsrMatrix(CurrentStateMapping);
                        //MassFact.AccMassMatrix(m_Stack_MassMatrix[0], CurrentStateMapping, _alpha: Config_MassScale);
                        base.ComputeMassMatrixImpl(m_Stack_MassMatrix[0], m_LsTrk.RegionsHistory.Current.Time);
                    }

                    PrecondMassMatrix = m_PrecondMassMatrix;

                    CoupledIteration = true;
                }


                // operator matrix update
                // ----------------------

                // we perform the extrapolation to have valid parameters if
                // - the operator matrix depends on these values
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.m_CurrentAgglomeration.Extrapolate(CurrentStateMapping);

                // clear operator matrix (clearing and re-alloc are pretty equal, i.e. 'BlockMsrMatrix.Clear()' just releases all internal memory)
                BlockMsrMatrix OpMatrix;
                if(Linearization)
                    OpMatrix = new BlockMsrMatrix(CurrentStateMapping);
                else
                    OpMatrix = null;

                // clear affine part
                double[] OpAffine = new double[CurrentStateMapping.LocalLength];
                
                // assemble matrix & affine part
                Debug.Assert(OpMatrix == null || OpMatrix.InfNorm() == 0);
                Debug.Assert(OpAffine.L2Norm() == 0);
                Debug.Assert(object.ReferenceEquals(this.m_CurrentAgglomeration.Tracker, this.m_LsTrk));
                this.ComputeOperatorMatrix(OpMatrix, OpAffine, CurrentStateMapping, locCurSt, base.GetAgglomeratedLengthScales(), m_CurrentPhystime + m_CurrentDt, 1);


                // assemble system
                // ---------------

                double dt = m_CurrentDt;

                double[] CurrentAffine = OpAffine;
                BlockMsrMatrix CurrentOpMatrix = OpMatrix;
                BlockMsrMatrix CurrentMassMatrix = m_Stack_MassMatrix[0];
                               

                // right-hand-side, resp. affine vector
                double[] RHS = new double[CurrentAffine.Length];
                if (this.Config_LevelSetHandling == LevelSetHandling.Coupled_Once
                    || this.Config_LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // MovingMesh:
                    // (1/dt)*(M1*u1 - M0*u0) + theta1*(Op1*u1 + b1) + theta0*(Op0*u0 + b0);
                    // RHS = (1/dt)*M0*u0 - theta1*b1 - theta0*Op0*u0 - theta0*b0
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    for (int s = 1; s <= Tsc.S; s++) { // loop over BDF stages
                        if (m_Stack_MassMatrix[s] != null) {
                            m_Stack_MassMatrix[s].SpMV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s], 1.0, RHS); //   (1/dt)*M0*u0
                        } else {
                            // mass matrix is identity
                            Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                            RHS.AccV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s]);
                        }
                    }

                    RHS.AccV(-Tsc.theta1, CurrentAffine); //                                     -theta1*b1
                    if (Tsc.theta0 != 0.0) {
                        double[] evalBuffer = new double[RHS.Length];
                        this.ComputeOperatorMatrix(null, evalBuffer, m_Stack_u[1].Mapping, m_Stack_u[1].Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime, 0);
                        RHS.AccV(-Tsc.theta0, evalBuffer); // RHS -= -theta0*Op0(u0)
                    }

                } else if (this.Config_LevelSetHandling == LevelSetHandling.LieSplitting
                    || this.Config_LevelSetHandling == LevelSetHandling.StrangSplitting
                    || Config_LevelSetHandling == LevelSetHandling.FSILieSplittingFullyCoupled
                    || this.Config_LevelSetHandling == LevelSetHandling.None) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Splitting:
                    // (1/dt)M1*(u1 - u0) + theta1*(Op1*u1 + b1) + theta0*(Op1*u0 + b1);
                    // Note: only the new operator is used!!!!!!!
                    // RHS = (1/dt)*M1*u0 - theta1*b1 - theta0*Op1*u0 - theta0*b1
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    for (int s = 1; s <= Tsc.S; s++) { // loop over BDF stages
                        if (CurrentAffine != null) {
                            if (CurrentMassMatrix != null) {
                                CurrentMassMatrix.SpMV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s], 1.0, RHS); //   (1/dt)*M0*u0 
                            } else {
                                Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                                RHS.AccV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s]);
                            }
                        } else {
                            Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                            RHS.AccV(Tsc.beta[s - 1] / dt, this.m_Stack_u[s]);
                        }
                    }

                    RHS.AccV(-Tsc.theta1, CurrentAffine); //                                     -theta1*b1
                    if (Tsc.theta0 != 0.0) {
                        // For XDG & splitting, we have to evaluate the operator with 
                        // the _actual_ level set position, but 
                        // everything else (field state, etc.) from at the _old_ timestep.
                                                
                        double[] evalBuffer = new double[RHS.Length];
                        this.ComputeOperatorMatrix(null, evalBuffer, m_Stack_u[1].Mapping, m_Stack_u[1].Fields.ToArray(), base.GetAgglomeratedLengthScales(), m_CurrentPhystime, 1);
                        RHS.AccV(-Tsc.theta0, evalBuffer); // RHS -= -theta0*Op(u0) at current level-set
                    }

                } else {
                    throw new NotImplementedException();
                }
                Affine = RHS;
                Affine.ScaleV(-1.0);

                // left-hand-side
                if(Linearization) {
                    System = CurrentOpMatrix.CloneAs();
                    if(Tsc.theta1 != 1.0)
                        System.Scale(Tsc.theta1);
                } else {
                    System = null;
                }
                if (CurrentMassMatrix != null) {
                    if(Linearization) {
                        System.Acc(1.0 / dt, CurrentMassMatrix);
                    } else {
                        CurrentMassMatrix.SpMV(1.0 / dt, new CoordinateVector(CurrentStateMapping), 1.0, Affine);
                    }
                } else {
                    Debug.Assert(Config_MassMatrixShapeandDependence == MassMatrixShapeandDependence.IsIdentity);
                    if(Linearization) {
                        System.AccEyeSp(1.0 / dt);
                    } else {
                        Affine.AccV(1.0 / dt, new CoordinateVector(CurrentStateMapping));
                    }
                }

                // perform agglomeration
                // ---------------------
                Debug.Assert(object.ReferenceEquals(m_CurrentAgglomeration.Tracker, m_LsTrk));
                m_CurrentAgglomeration.ManipulateMatrixAndRHS(System, Affine, CurrentStateMapping, CurrentStateMapping);
                
                if(Linearization) {
                    m_LsTrk.CheckVectorZeroInEmptyCutCells(Affine, CurrentStateMapping, this.Config_SpeciesToCompute, m_CurrentAgglomeration, this.Config_CutCellQuadratureOrder);
                    m_LsTrk.CheckMatrixZeroInEmptyCutCells(System, CurrentStateMapping, this.Config_SpeciesToCompute, m_CurrentAgglomeration, this.Config_CutCellQuadratureOrder);
                } else {
                    m_LsTrk.CheckVectorZeroInEmptyCutCells(Affine, CurrentStateMapping, this.Config_SpeciesToCompute, m_CurrentAgglomeration, this.Config_CutCellQuadratureOrder);
                }

                // increase iteration counter         
                // --------------------------
                //OpMatrix.SaveToTextFileSparse("OpMatrix49e.txt");
                abstractOperator = AbstractOperator;
                m_IterationCounter++;
            }
        }

    }

    */
}
