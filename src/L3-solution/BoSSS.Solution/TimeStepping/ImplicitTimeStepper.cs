/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation;
using ilPSP;
using ilPSP.LinSolvers;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Base class for linear implicit time stepper (e.g. implicit euler). If
    /// the equation to be solved is denoted by dx/dt + F(x) = 0, this class
    /// provides the decomposition of the (linear)
    /// <see cref="BoSSS.Foundation.SpatialOperator"/> F(x) =
    /// <em>M</em>*x + <see cref="m_AffineOffset1"/> which is required for
    /// the execution of the timestep. In this class, x is given by
    /// <see cref="CurrentState"/>
    /// </summary>
    public abstract class ImplicitTimeStepper : ITimeStepper, IDisposable {

        /// <summary>
        /// The sparse solver used to solve the equation system in
        /// <see cref="PerformTimeStep"/>
        /// </summary>
        /// <remarks>
        /// The Matrix of the system is stored within the sparse solver;
        /// </remarks>
        protected ISparseSolverExt m_Solver;

        ///// <summary>
        ///// 
        ///// </summary>
        //protected ISparseMatrix eM {
        //    get { return m_Solver.GetMatrix(); }
        //}
        
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

        private double m_Time = 0.0;
        private CoordinateMapping m_Mapping;
        private CoordinateVector m_DGCoordinates;

        /// <summary>
        /// Constructs an implicit timestepper from a given matrix (i.e. the
        /// spatial discretization is done elsewhere);
        /// </summary>
        public ImplicitTimeStepper(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields) {
            Setup1(solver, temporalOp, spatialOpMtx, spatialOpAffine, fields);
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="temporalOp">
        /// Indicates, for each each equation whether it is
        /// <list type="bullet">
        ///   <item>(false) a auxiliary condition  i.e. a variable where no time derivative occurs in the equation, or</item>
        ///   <item>(true) a differential equation</item>
        /// </list>
        /// At least one equation must be time-dependent;
        /// Otherwise, the <see cref="BoSSS.Solution.Solvers.LinearSolver"/> should be used;
        /// </param>
        /// <param name="spatialOp"></param>
        /// <param name="fields"></param>
        /// <param name="solver">
        /// </param>
        public ImplicitTimeStepper(ISparseSolverExt solver, bool[] temporalOp, SpatialOperator spatialOp, CoordinateMapping fields) {
            
            // check operator and arguments
            if (!spatialOp.IsCommited)
                throw new ArgumentException("operator must be committed first.", "spatialOp");
            if (spatialOp.ContainsNonlinear)
                throw new ArgumentException("spatial differential operator cannot contain nonlinear components for implicit euler.", "spatialOp");
            if (!spatialOp.ContainsLinear())
                throw new ArgumentException("spatial differential operator seems to contain no components.", "spatialOp");
            if (spatialOp.DomainVar.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("spatial differential operator must have the same number of domain and codomain variables.", "spatialOp");

            if (fields.Fields.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("the number of fields in the coordinate mapping must be equal to the number of domain/codomain variables of the spatial differential operator", "fields");

            m_Solver = solver;
            m_Mapping = fields;
            //m_AffineOffset1 = new double[fields.LocalLength];
            MsrMatrix eM;
            ImplicitTimeStepper.ComputeMatrix(spatialOp, fields, false, out eM, out m_AffineOffset1);

            //Setup2(eM, temporalOp);
            Setup1(solver, temporalOp, eM, m_AffineOffset1, fields);
        }

        /// <summary>
        /// 
        /// </summary>
        protected ImplicitTimeStepper() {}

        /// <summary>
        /// 
        /// </summary>
        public static TimeStepperType Factory<TimeStepperType>(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields)
            where TimeStepperType : ImplicitTimeStepper, new() {

            TimeStepperType timestepper = new TimeStepperType();
            timestepper.Setup1(solver, temporalOp, spatialOpMtx, spatialOpAffine, fields);
            
            //timestepper.Setup2(
            return timestepper;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="fields"></param>
        /// <param name="matrix"></param>
        /// <param name="affineOffset"></param>
        /// <param name="OnlyAffine">
        /// if true, only the <paramref name="affineOffset"/>
        /// is computed and <paramref name="matrix"/> is null on exit.
        /// </param>
        public static void ComputeMatrix(SpatialOperator spatialOp, CoordinateMapping fields, bool OnlyAffine, out MsrMatrix matrix, out double[] affineOffset) {
            // Check operator and arguments
            if (!spatialOp.IsCommited)
                throw new ArgumentException("operator must be committed first.", "spatialOp");
            if (spatialOp.ContainsNonlinear)
                throw new ArgumentException("spatial differential operator cannot contain nonlinear components for implicit euler.", "spatialOp");
            if (!spatialOp.ContainsLinear())
                throw new ArgumentException("spatial differential operator seems to contain no components.", "spatialOp");
            if (spatialOp.DomainVar.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("spatial differential operator must have the same number of domain and codomain variables.", "spatialOp");
            if (fields.Fields.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("the number of fields in the coordinate mapping must be equal to the number of domain/codomain variables of the spatial differential operator", "fields");

            // Assemble matrix and affine offset
			IPartitioning matrixPartition = fields;


            if (!OnlyAffine)
                matrix = new MsrMatrix(matrixPartition);
            else
                matrix = null;
            affineOffset = new double[fields.LocalLength];
            spatialOp.ComputeMatrixEx( 
                                      fields, null, fields,
                                      matrix, affineOffset,
                                      OnlyAffine);
        }

        /// <summary>
        /// Sparse solver which is used to perform the implicit timestep. This
        /// property can e.g. be used to change some properties of the solver
        /// </summary>
        public ISparseSolver Solver { get { return m_Solver; } }

        /// <summary>
        /// Returns an array of bool's of size <paramref name="n"/> where all
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
        /// timestep (that is, by updating <see cref="CurrentState"/>).
        /// </summary>
        /// <param name="dt"></param>
        abstract protected void PerformTimeStep(double dt);

        /// <summary>
        /// common for all constructors
        /// </summary>
        protected void Setup1(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields) {
            // check operator and arguments
            if (spatialOpMtx.NoOfRows != spatialOpMtx.NoOfCols)
                throw new ArgumentException("matrix must be quadratic.", "spatialOpMtx");
            if (spatialOpMtx.NoOfRows != fields.GlobalCount)
                throw new ArgumentException("matrix size must be equal to the GlobalCount of fields mapping", "fields,spatialOpMtx");
            if (spatialOpMtx.RowPartitioning.LocalLength != fields.LocalLength)
                throw new ArgumentException("number of locally stored matrix rows nust be equal to NUpdate of fields mapping.", "fields,spatialOpMtx");
            if (spatialOpAffine.Count < fields.LocalLength)
                throw new ArgumentException("length affine offset vector must be equal or larger than NUpdate of the mapping", "spatialOpAffine");

            m_Solver = solver;
            m_Mapping = fields;
            m_AffineOffset1 = spatialOpAffine.ToArray();

            //Setup2(spatialOpMtx, temporalOp);
            //}

            ///// <summary>
            ///// Common to all constructors; Passes the Matrix <paramref name="eM"/>
            ///// to the solver and initializes <see cref="m_diagVecOneSec"/>;
            ///// </summary>
            ///// <param name="eM"><see cref="eM"/></param>
            ///// <param name="temporalOp"><see cref="ImplicitTimeStepper"/></param>
            //private void _Setup2(MsrMatrix eM, bool[] temporalOp) {

            {
                // check temporal operator
                // -----------------------
                if (m_Mapping.Fields.Count != temporalOp.Length)
                    throw new ArgumentException(
                        "length of temporalOp must be equal to number of domain/codomain variables of the spatial differential operator",
                        "temporalOp");
                m_DGCoordinates = new CoordinateVector(m_Mapping);

                bool timedep = false;
                bool fullyTimeDep = true;
                foreach (bool b in temporalOp) {
                    timedep = timedep || b;
                    fullyTimeDep = fullyTimeDep & b;
                }
                if (!timedep) {
                    throw new ArgumentException("at least one equation must be time-dependent; one entry in temporalOp must be true;", "temporalOp");
                }

                // Construct diagonal matrix vector
                // --------------------------------
                m_diagVecOneSec = new double[m_Mapping.MaxTotalNoOfCoordinatesPerCell];
                IList<DGField> _fields = m_Mapping.Fields;
                for (int g = 0; g < temporalOp.Length; g++) {
                    if (temporalOp[g]) {
                        DGField gamma = _fields[g];
                        for (int n = 0; n < gamma.Basis.MaximalLength; n++) {
                            m_diagVecOneSec[m_Mapping.LocalUniqueCoordinateIndex(g, 0, n)] = 1.0;
                        }
                    }
                }

                // initialize linear solver
                m_Solver.DefineMatrix(spatialOpMtx);
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
        public void ResetTime(double NewTime, int timestepNumber) {
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
        public double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {


                BeforeTimeStep(dt);
                PerformTimeStep(dt);
                m_Time += dt;
                AfterTimeStep(dt);

            }
            return dt;
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
        /// sets <see cref="Mapping"/> to <paramref name="m"/>
        /// and creates a new <see cref="CurrentState"/> - mapping;
        /// </summary>
        /// <param name="m"></param>
        protected void SetMapping(CoordinateMapping m) {
            m_Mapping = m;
            m_DGCoordinates = new CoordinateVector(m);
        }

        /// <summary>
        /// The mapping of the DG fields that are affected by this timestepper
        /// </summary>
        /// <remarks>
        /// the initial value is read from those fields and 
        /// the solution is written back to it;
        /// </remarks>
        public CoordinateMapping Mapping {
            get { return m_Mapping; }
        }

        /// <summary>
        /// The DG coordinates of the <see cref="BoSSS.Foundation.DGField"/>'s in
        /// <see cref="Mapping"/> (see
        /// <see cref="BoSSS.Foundation.CoordinateMapping.Fields"/>)
        /// </summary>
        public CoordinateVector CurrentState {
            get { return m_DGCoordinates; }
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
}