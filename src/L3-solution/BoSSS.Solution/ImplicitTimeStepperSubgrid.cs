
 /*
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System;
using System.Collections.Generic;

using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using System.Linq;

namespace BoSSS.Solution.Timestepping {

#pragma warning disable 1572
#pragma warning disable 1573
#pragma warning disable 1591

    /// <summary>
    /// Base class for linear implicit time stepper (e.g. implicit euler). If
    /// the equation to be solved is denoted by dx/dt + F(x) = 0, this class
    /// provides the decomposition of the (linear)
    /// <see cref="BoSSS.Foundation.SpatialOperator"/> F(x) =
    /// <em>M</em>*x + <see cref="m_AffineOffset1"/> which is required for
    /// the execution of the timestep. In this class, x is given by
    /// <see cref="DGCoordinates"/>
    /// </summary>
    public abstract class ImplicitTimeStepperSubgrid : ITimeStepper, IDisposable {

        /// <summary>
        /// The omnipresent context
        /// </summary>
        protected Context m_Context;

        /// <summary>
        /// The sparse solver used to solve the equation system in
        /// <see cref="PerformTimeStep"/>
        /// </summary>
        /// <remarks>
        /// The Matrix of the system is stored within the sparse solver;
        /// </remarks>
        protected ISparseSolver m_Solver;

    
        /// <summary>
        /// The affine offset b from F(x) = Mx + b
        /// </summary>
        /// <remarks>
        /// If boundary conditions are time-dependent, this vector may change over 
        /// time; Re-calculation can be implemented e.g. in <see cref="BeforeTimeStep"/>.
        /// At a given time <em>t</em>, for given initial values
        /// it is 
        /// the convention is that this vector represents the inhomogeouos b.c.
        /// at time <em>t</em>+ <em>dt</em>, i.e. at the next timestep.
        /// </remarks>
        protected double[] m_AffineOffset1;
        
        /// <summary>
        /// One period of the diagonal vector
        /// </summary>
        protected double[] m_diagVecOneSec;
        /// <summary>
        /// Current time
        /// </summary>
        private double m_Time = 0.0;
        /// <summary>
        /// Coordinate Mapping defined on the Subgrid
        /// </summary>
        public SubgridCoordinateMapping m_SubgridMapping;
        /// <summary>
        /// DG Coordinates of the problem on the restricted narrow band domain
        /// </summary>
        protected double[] m_SubgridDGCoordinates;
        /// <summary>
        /// The subgrid (i.e. the narrow band)
        /// </summary>
        public SubGrid m_Subgrid;
        /// <summary>
        /// Operator matrix defined on the subgrid
        /// </summary>
        protected MsrMatrix m_CompressedMatrix;
        /// <summary>
        /// Affine part of the operator on the subgrid domain
        /// </summary>
        protected double[] m_CompressedAffine;
        /// <summary>
        /// Indicator for time dependent components of the problem
        /// </summary>
        public bool[] temporalOp;

       /// <summary>
       /// Constructor for implicit time stepping scheme using a predefined discrization of the differential operator
       /// </summary>
       /// <param name="ctx"></param>
       /// <param name="solver">Solver</param>
       /// <param name="temporalOp">true for each component which should be time dependent</param>
       /// <param name="spatialOpMtx">Matrix associated with the differential operator</param>
       /// <param name="spatialOpAffine">Affine part of the differential operator</param>
       /// <param name="subgrid">Subgrid which the problem should be computed on</param>
       /// <param name="fields">Mapping on the subgrid</param>
       /// <param name="Initialdt">Time step size</param>
        public ImplicitTimeStepperSubgrid(Context ctx, ISparseSolver solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine,
            SubGrid subgrid, SubgridCoordinateMapping fields, double Initialdt) {
            Setup1(ctx, solver, temporalOp, spatialOpMtx, spatialOpAffine,subgrid, fields, Initialdt);
        }

        /// <summary>
        /// Constructor for an empty scheme
        /// </summary>
        protected ImplicitTimeStepperSubgrid() { }

       /// <summary>
       /// Factory for implicit time stepping schemes that operate on a subgrid
       /// </summary>
       /// <typeparam name="TimeStepperType"></typeparam>
       /// <param name="ctx"></param>
        /// <param name="solver">Solver</param>
        /// <param name="temporalOp">true for each component which should be time dependent</param>
        /// <param name="spatialOpMtx">Matrix associated with the differential operator</param>
        /// <param name="spatialOpAffine">Affine part of the differential operator</param>
        /// <param name="subgrid">Subgrid which the problem should be computed on</param>
        /// <param name="fields">Mapping on the subgrid</param>
        /// <param name="Initialdt">Time step size</param>
       /// <returns></returns>
        public static TimeStepperType Factory<TimeStepperType>(Context ctx, ISparseSolver solver, bool[] temporalOp, MsrMatrix spatialOpMtx,
            IList<double> spatialOpAffine, SubGrid subgrid, SubgridCoordinateMapping fields, double Initialdt)
            where TimeStepperType : ImplicitTimeStepperSubgrid, new() {

            TimeStepperType timestepper = new TimeStepperType();
            timestepper.Setup1(ctx, solver, temporalOp, spatialOpMtx, spatialOpAffine, subgrid, fields, Initialdt);
            return timestepper;
        }
        
        /// <summary>
        /// Sparse solver which is used to perform the implicit timestep. This
        /// property can e.g. be used to change some properties of the solver
        /// </summary>
        public ISparseSolver Solver { get { return m_Solver; } }

        /// <summary>
        /// Returns an array of bools of size <paramref name="n"/> where all
        /// entries are set to "true"
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        protected static bool[] AllTrue(int n) {
            bool[] r = new bool[n];
            for (int i = 0; i < n; i++) r[i] = true;
            return r;
        }

        /// <summary>
        /// Implement this method in any subclass in order to carry out the
        /// timestep (that is, by updating <see cref="DGCoordinates"/>).
        /// </summary>
        /// <param name="dt"></param>
        abstract protected void PerformTimeStep(double dt);

        /// <summary>
        /// Sets matrix for the timestepping with time step size <paramref name="dt"/> in the solver. 
        /// </summary>
        /// <param name="OperatorMatrix"></param>
        /// <param name="dt"></param>
        abstract protected void DefineMatrix(MsrMatrix OperatorMatrix, double dt);
        public MsrMatrix SubgridMatrix {
            get { return m_CompressedMatrix; }
        }
        public double[] SubgridAffine {
            get { return m_CompressedAffine; }
        }
        /// <summary>
        /// method for setting up the timestepper, i.e. the necessary
        /// </summary>
        protected void Setup1(Context ctx, ISparseSolver solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, 
            SubGrid subgrid,SubgridCoordinateMapping fields, double InitialDeltat) {
            // check operator and arguments
            if (spatialOpMtx.NoOfRows != spatialOpMtx.NoOfCols)
                throw new ArgumentException("matrix must be quadratic.", "spatialOpMtx");
            if (spatialOpMtx.NoOfRows != fields.GlobalCount)
                throw new ArgumentException("matrix size must be equal to the GlobalCount of fields mapping", "fields,spatialOpMtx");
            if (spatialOpMtx.RowPartition.LocalLength != fields.NUpdate)
                throw new ArgumentException("number of locally stored matrix rows nust be equal to NUpdate of fields mapping.", "fields,spatialOpMtx");
            if (spatialOpAffine.Count < fields.NUpdate)
                throw new ArgumentException("length affine offset vector must be equal or larger than NUpdate of the mapping", "spatialOpAffine");

            m_Context = ctx;
            m_Solver = solver;
            m_Subgrid = subgrid;
            m_SubgridMapping = fields;
            m_AffineOffset1 = spatialOpAffine.ToArray();
            this.temporalOp = temporalOp;
            BLAS.dscal(m_AffineOffset1.Length, -1.0, m_AffineOffset1,1);
            {
                // check temporal operator
                // -----------------------
                if (m_SubgridMapping.Fields.Count != temporalOp.Length)
                    throw new ArgumentException(
                        "lenght of temporalOp must be equal to number of domain/codomain variables of the spatial differential operator",
                        "temporalOp");
                m_SubgridMapping.Compress();
                m_SubgridDGCoordinates = m_SubgridMapping.subgridCoordinates;
                
                m_SubgridMapping.SetupSubmatrix(m_AffineOffset1, spatialOpMtx, out m_CompressedAffine, out m_CompressedMatrix);
     
                bool timedep = false;
                bool fullyTimeDep = true;
                foreach (bool b in temporalOp) {
                    timedep = timedep || b;
                    fullyTimeDep = fullyTimeDep & b;
                }
                if (!timedep) {
                    throw new ArgumentException("At least one equation must be time-dependent; one entry in temporalOp must be true;", "temporalOp");
                }

                DefineMatrix(m_CompressedMatrix, InitialDeltat);
              
            }
        }
        
        #region ITimeStepper Member

        /// <summary>
        /// <see cref="ITimeStepper.Time"/>
        /// </summary>
        public double Time {
            get { return m_Time; }
        }

        /// <summary>
        /// <see cref="ITimeStepper.ResetTime"/>
        /// </summary>
        /// <param name="NewTime">The new time</param>
        public void ResetTime(double NewTime) {
            m_Time = NewTime;
        }

        /// <summary>
        /// Increments the time and performs a timestep.
        /// </summary>
        /// <param name="dt">The length of the timestep</param>
        /// <remarks>
        /// Leaves the actual implementation of the timestep to the subclass
        /// (<see cref="PerformTimeStep"/>)
        /// </remarks>
        public  void Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {
                //m_CompressedMatrix.SaveToTextFileSparse("c:\\tmp\\matrix-boese.txt");


                if (m_Context == null) {
                    throw new ApplicationException("Cannot perform a timestep before timestepper has been set-up");
                }

                BeforeTimeStep(dt);
                PerformTimeStep(dt);
                m_Time += dt;
                AfterTimeStep(dt);

            }
        }

        /// <summary>
        /// invoked by <see cref="Perform"/> before calling <see cref="PerformTimeStep"/>
        /// </summary>
        /// <param name="dt">size of timestep</param>
        virtual protected void BeforeTimeStep(double dt) { }

        /// <summary>
        /// invoked by <see cref="Perform"/> after calling <see cref="PerformTimeStep"/>
        /// </summary>
        /// <param name="dt">size of timestep</param>
        virtual protected void AfterTimeStep(double dt) { }

        /// <summary>
        /// sets the <see cref="SubgridCoordinateMapping"/> to <paramref name="m"/>
        /// The new mapping is compressed again to the narrow band only.
        /// </summary>
        /// <param name="m"></param>
        protected void SetMapping(SubgridCoordinateMapping m) {
            m_SubgridMapping = m;
            m_SubgridMapping.Compress();
            m_SubgridDGCoordinates = m_SubgridMapping.subgridCoordinates;
        }
        /// <summary>
        /// DG Coordinates on the subgrid only.
        /// </summary>
        protected double[] SubgridDGCoordinates {
            get {
                return m_SubgridDGCoordinates;
            }
            set { m_SubgridDGCoordinates=value;
            }
        }
        /// <summary>
        /// Subgrid mapping of the given mapping which the time stepper is operating on.
        /// </summary>
        public SubgridCoordinateMapping SubgridMapping {
            get {
                return m_SubgridMapping; }
        }
        /// <summary>
        /// Coordinate Mapping on the full domain
        /// </summary>
        public CoordinateMapping Mapping {
            get {
                m_SubgridMapping.Decompress();
                return (CoordinateMapping)m_SubgridMapping;
            }
        }

        /// <summary>
        /// The CoordinatesVector associated with this Mapping 
        /// </summary>
        public CoordinateVector DGCoordinates {
            get {
                m_SubgridMapping.Decompress();
                return new CoordinateVector(m_SubgridMapping); }
        }
        
        #endregion

        #region IDisposable Member

        /// <summary>
        /// Frees the umanaged memory that is maybe used by the sparse solver
        /// </summary>
        public void Dispose() {
            IDisposable solver = m_Solver as IDisposable;
            if (solver != null) {
                solver.Dispose();
            }
        }

        #endregion
    }

#pragma warning restore 1572
#pragma warning restore 1573
#pragma warning restore 1591
}