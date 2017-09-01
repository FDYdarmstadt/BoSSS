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
using ilPSP.Utils;

namespace BoSSS.Solution.Solvers {


    /// <summary>
    /// Event handler used in <see cref="LinearSolver"/>;
    /// See <see cref="LinearSolver.RHSEvaluated"/>;
    /// </summary>
    /// <param name="rhs"></param>
    public delegate void OnRhsEvaluatedHandler(double[] rhs);

    /// <summary>
    /// a solver for fully-linear systems of non-evolution - Equations
    /// </summary>
    public class LinearSolver {

        /// <summary>
        /// ctor.
        /// </summary>
        public LinearSolver(ISparseSolver solver, MsrMatrix spatialOpMatrix, IList<double> spatialOpAffine, params DGField[] Unknowns)
            : this(solver, spatialOpMatrix, spatialOpAffine, new CoordinateMapping(Unknowns)) {
        }


        /// <summary>
        /// ctor.
        /// </summary>
        public LinearSolver(ISparseSolver solver, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping UnknownsMap) {
            // check operator and arguments
            // ----------------------------

            if (spatialOpMtx.NoOfRows != spatialOpMtx.NoOfCols)
                throw new ArgumentException("matrix must be quadratic.", "spatialOpMtx");
            if (spatialOpMtx.NoOfRows != UnknownsMap.GlobalCount)
                throw new ArgumentException("matrix size must be equal to the GlobalCount of fields mapping", "fields,spatialOpMtx");
            if (spatialOpMtx.RowPartitioning.LocalLength != UnknownsMap.LocalLength)
                throw new ArgumentException("number of locally stored matrix rows nust be equal to NUpdate of fields mapping.", "fields,spatialOpMtx");
            if (spatialOpAffine.Count < UnknownsMap.LocalLength)
                throw new ArgumentException("length affine offset vector must be equal or larger than NUpdate of the mapping", "spatialOpAffine");

            m_Solver = solver;
            m_Mapping = UnknownsMap;


            // finish constructor
            // ------------------
            m_AffineOffset = spatialOpAffine.ToArray();
            ConstructorCommon(spatialOpMtx);
        }

        /// <summary>
        /// empty constructor;
        /// </summary>
        protected LinearSolver() {

        }


        /// <summary>
        /// ctor.
        /// </summary>
        public LinearSolver(ISparseSolver solver, MsrMatrix spatialOpMatrix, IList<double> spatialOpAffine, CoordinateMapping UnknownsMap, SpatialOperator rhsOperator, CoordinateMapping rhsDomainFields)
            : this(solver, spatialOpMatrix, spatialOpAffine, UnknownsMap) {
            ConstructorCommon2(UnknownsMap, rhsOperator, rhsDomainFields);
        }

        /// <summary>
        /// ctor.
        /// </summary>
        public LinearSolver(ISparseSolver solver, SpatialOperator spatialOp, CoordinateMapping UnknownsMap, SpatialOperator rhsOperator, CoordinateMapping rhsDomainFields)
            : this(solver, spatialOp, UnknownsMap) {

            ConstructorCommon2(UnknownsMap, rhsOperator, rhsDomainFields);
        }


        /// <summary>
        /// common to more than one constructor
        /// </summary>
        protected void ConstructorCommon2(CoordinateMapping UnknownsMap, SpatialOperator rhsOperator, CoordinateMapping rhsDomainFields) {
            // verify input
            // ------------

            if (m_Mapping.Fields.Count != rhsOperator.CodomainVar.Count)
                throw new ArgumentException("spatial differential operator and RHS obperator must have the same number of domain variables.", "rhsOperator");

            // construct evaluator
            // -------------------

            m_rhsEvaluator = rhsOperator.GetEvaluator(rhsDomainFields, UnknownsMap);
        }

        /// <summary>
        /// ctor.
        /// </summary>
        public LinearSolver(ISparseSolver solver, SpatialOperator spatialOp, params DGField[] Unknowns)
            : this(solver, spatialOp, new CoordinateMapping(Unknowns)) {
        }

        /// <summary>
        /// ctor.
        /// </summary>
        public LinearSolver(ISparseSolver solver, SpatialOperator spatialOp, CoordinateMapping UnknownsMap) {

            // verify input
            // ------------
            if (!spatialOp.IsCommited)
                throw new ArgumentException("operator must be committed first.", "spatialOp");
            if (spatialOp.ContainsNonlinear)
                throw new ArgumentException("spatial differential operator cannot contain nonlinear components for linear solver.", "spatialOp");
            if (!spatialOp.ContainsLinear())
                throw new ArgumentException("spatial differential operator seems to contain no components.", "spatialOp");
            if (spatialOp.DomainVar.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("spatial differential operator must have the same number of domain and codomain variables.", "spatialOp");

            if (UnknownsMap.Fields.Count != spatialOp.CodomainVar.Count)
                throw new ArgumentException("the number of fields in the coordinate mapping must be equal to the number of domain/codomain variables of the spatial differential operator", "fields");

            m_Mapping = UnknownsMap;
            m_Solver = solver;

            // matrix assembly
            // ---------------
            MsrMatrix eM;
            {
                eM = new MsrMatrix(m_Mapping);
                m_AffineOffset = new double[m_Mapping.LocalLength];

                spatialOp.ComputeMatrixEx(m_Mapping, null, m_Mapping, eM, m_AffineOffset);
            }

            ConstructorCommon(eM);

        }

        /// <summary>
        /// Defines the matrix of the sparse solver <see cref="m_Solver"/>
        /// </summary>
        /// <param name="_eM"></param>
        protected virtual void ConstructorCommon(MsrMatrix _eM) {
            //eM = _eM;
            m_Solver.DefineMatrix(_eM);
            m_DGCoordinates = new CoordinateVector(m_Mapping);
        }

        /// <summary>
        /// <see cref="RhsEvaluator"/>
        /// </summary>
        SpatialOperator.Evaluator m_rhsEvaluator;

        /// <summary>
        /// optional right-hand side of the equation;
        /// This property is set (unequal to null) if the 
        /// <see cref="LinearSolver(ISparseSolver,SpatialOperator,CoordinateMapping, SpatialOperator,CoordinateMapping)"/>-constructor
        /// is used for the construction of this object.
        /// </summary>
        public SpatialOperator.Evaluator RhsEvaluator {
            get {
                return m_rhsEvaluator;
            }
        }

        ///// <summary>
        ///// Tranceiver for the fields within <see cref="Mapping"/>
        ///// </summary>
        //Transceiver m_TRX;

        ///// <summary>
        ///// constructs a linear solver from a allready pre-assembled
        ///// matrix;
        ///// The system that is solved is defined as:
        ///// <paramref name="EqSystem"/>*Unknown + <paramref name="AffineOffset"/> = rhs, 
        ///// where the mapping of "Unknown" to Dg Fields and coordinater is done by 
        ///// <paramref name="Mapping"/>; "rhs" is provided during the invocation of 
        ///// <see cref="Solve"/>;
        ///// </summary>
        ///// <param name="context">as usual</param>
        ///// <param name="EqSystem"></param>
        ///// <param name="AffineOffset"></param>
        ///// <param name="Mapping"></param>
        ///// <param name="solver">
        ///// Linear eqation system solver which should be used for implicit timestepping; 
        ///// </param>
        //public LinearSolver(Context context, MsrMatrix EqSystem, double[] AffineOffset, CoordinateMapping Mapping, ISparseSolver solver) {

        //    int th = context.IOMaster.tracer.EnterFunction("BoSSS.Solution.Solvers.LinearSolver(Context,MsrMatrix,double[],CoordinateMapping,ISparseSolver)");


        //    if (Mapping.LocalLength != EqSystem.RowPartiton.LocalLength
        //        || Mapping.GlobalCount != EqSystem.RowPartiton.TotalLength)
        //        throw new ArgumentException("mismatch between matrix and mapping", "EqSystem,Mapping");
        //    if (AffineOffset.Length != EqSystem.RowPartiton.LocalLength)
        //        throw new ArgumentException("AffineOffset.Length must be equal to EqSystem.RowPartiton.LocalLength");

        //    m_Context = context;
        //    ConstructorCommon(EqSystem, AffineOffset, Mapping, solver);

        //    context.IOMaster.tracer.LeaveFunction(th);
        //}



        ///// <summary>
        ///// some things that all constructors share in common
        ///// </summary>
        ///// <param name="EqSystem"></param>
        ///// <param name="AffineOffset"></param>
        ///// <param name="Mapping"></param>
        ///// <param name="solver"></param>
        //private void ConstructorCommon(MsrMatrix EqSystem, double[] AffineOffset, CoordinateMapping Mapping, ISparseSolver solver) {

        //    m_Mapping = Mapping;
        //    eM = EqSystem;
        //    m_AffineOffset = AffineOffset;
        //    m_Solver = solver;

        //    // feed linear solver
        //    // ------------------

        //    if (m_MatrixBuild != null)
        //        m_MatrixBuild(EqSystem, AffineOffset, Mapping);

        //    m_Solver.DefineMatrix(EqSystem);
        //    GC.Collect();


        //    // create transceiver
        //    // ------------------
        //    m_TRX = new Transceiver(m_Context.CommMaster, m_Mapping.Fields);


        //}


        ///// <summary>
        ///// constructs a linear solver from a collection of equations;
        ///// </summary>
        ///// <param name="context">as usual</param>
        ///// <param name="equations">The system of equations that should be solved;
        ///// All equations must be fully linear non-evolution equations;
        ///// </param>
        ///// <param name="solver">
        ///// Linear eqation system solver which should be used for implicit timestepping; 
        ///// </param>
        //public LinearSolver(Context context, ICollection<Equation> equations, ISparseSolver solver) :
        //    this(context, equations, solver, null) { }





        ///// <summary>
        ///// constructs a linear solver from a collection of equations;
        ///// </summary>
        ///// <param name="context">as usual</param>
        ///// <param name="equations">The system of equations that should be solved;
        ///// All equations must be fully linear non-evolution equations;
        ///// </param>
        ///// <param name="solver">
        ///// Linear eqation system solver which should be used for implicit timestepping; 
        ///// </param>
        ///// <param name="eve">
        ///// a method that is called after matrix assembly, to modify matrix and affine offset
        ///// of the operator manually
        ///// </param>
        //public LinearSolver(Context context, ICollection<Equation> equations, ISparseSolver solver, OnMatrixBuild eve) {

        //    // check input args and create mapping
        //    // -----------------------------------

        //    m_MatrixBuild = eve;
        //    m_Solver = solver;
        //    m_Context = context;
        //    int th = context.IOMaster.tracer.EnterFunction("BoSSS.Solution.Solvers.LinearSolver(Context,ICollection<Equation>,ISparseSolver)");

        //    Field[] fields = new Field[equations.Count];
        //    Equation[] eqsOrd = new Equation[equations.Count];
        //    int i = 0;
        //    foreach (Equation eq in equations) {
        //        if (Array.IndexOf<Field>(fields, eq.MyField) > 0)
        //            throw new ArgumentException("Only one equation per field is possible.");

        //        if (!eq.IsLinear())
        //            throw new ArgumentException("at least one equation contains nonlinear components");

        //        fields[i] = eq.MyField;
        //        eqsOrd[i] = eq;
        //        i++;
        //    }

        //    m_Mapping = new CoordinateMapping(context, fields);


        //    // create matrix
        //    // -------------


        //    //MsrMatrix eM = new MsrMatrix(m_Mapping, m_Mapping);
        //    {
        //        Partition part = new Partition(m_Mapping.LocalLength);
        //        eM = new MsrMatrix(part, (int)m_Mapping.GlobalCount);
        //        m_AffineOffset = new double[m_Mapping.LocalLength];

        //        {
        //            LECQuadratureEdge mxtbuilder2 = new LECQuadratureEdge(m_Context, eqsOrd, eM, m_AffineOffset, m_Mapping, m_Mapping);
        //            mxtbuilder2.Execute();
        //        }

        //        {
        //            LECVolumeQuadrature mtxBuilder = new LECVolumeQuadrature(m_Context, eqsOrd, eM, m_AffineOffset, m_Mapping, m_Mapping);
        //            mtxBuilder.Execute();
        //        }
        //    }


        //    ConstructorCommon(eM, m_AffineOffset, m_Mapping, solver);

        //    context.IOMaster.tracer.LeaveFunction(th);
        //}

        ///// <summary>
        ///// memory for the unknowns
        ///// </summary>
        //ScatteredVector m_CurrSubVec;

        ///// <summary>
        ///// matrix of the created operator
        ///// </summary>
        //public MsrMatrix eM;

        /// <summary>
        /// 
        /// </summary>
        public void Check(double[] rhs) {

            throw new NotImplementedException("only for test purpose");

            //// communicate 
            //Transceiver m_transciever = new Transceiver(m_Context.CommMaster, m_Mapping.Fields);
            //m_transciever.TransceiveStartImReturn();
            //m_transciever.TransceiveFinish();


            //double aonrm = BLAS.dnrm2(m_AffineOffset.Length, m_AffineOffset, 1);
            //Console.WriteLine("affine offset norm is: " + aonrm);



            //// evaluate the linear matrix
            //double[] res = new double[m_Mapping.LocalLength];
            //eM.Gemv<CoordinateMapping, double[]>(1.0, m_Mapping, 0, res);
            //BLAS.daxpy(res.Length, 1.0, m_AffineOffset, 1, res, 1);



            //// evaluate the reference equations
            //double[] k = new double[m_Mapping.LocalLength];

            //m_RefQuadEdge.m_Output = k;
            //m_RefQuadEdge.Execute();

            //m_RefQuadVol.m_Output = k;
            //m_RefQuadVol.Execute();


            //// compare
            //BLAS.daxpy(res.Length, -1.0, k, 1, res, 1); // after this, res should be approximately zero
            //double nrm = BLAS.dnrm2(res.Length, res, 1);
            //Console.WriteLine("check result is: " + nrm);

            //if (!File.Exists("matrixi.dat")) {
            //    FileStream fs = new FileStream("matrixi.dat", FileMode.Create);
            //    BinaryFormatter bf = new BinaryFormatter();

            //    bf.Serialize(fs, eM.m_Entries);
            //    bf.Serialize(fs, m_AffineOffset);
            //    bf.Serialize(fs, rhs);
            //    Console.WriteLine("em fl wrtn.");
            //    fs.Flush();
            //    fs.Close();

            //} else {
            //    FileStream fs = new FileStream("matrixi.dat", FileMode.Open);
            //    BinaryFormatter bf = new BinaryFormatter();


            //    MsrMatrix.MSREntry[][] matrixi = (MsrMatrix.MSREntry[][])bf.Deserialize(fs);
            //    for (int i = 0; i < matrixi.Length; i++) {
            //        MsrMatrix.MSREntry[] row = matrixi[i];
            //        if (row.Length != eM.m_Entries[i].Length)
            //            Console.WriteLine("mismatsch");


            //        for (int j = 0; j < Math.Min(row.Length, eM.m_Entries[i].Length); j++) {
            //            if (row[j].m_ColIndex != eM.m_Entries[i][j].m_ColIndex)
            //                Console.WriteLine("mixmatch");

            //            if (row[j].val != eM.m_Entries[i][j].val)
            //                Console.WriteLine("anal mixmatch");
            //        }
            //    }

            //    double[] affoff = (double[])bf.Deserialize(fs);
            //    if (affoff.Length != m_AffineOffset.Length)
            //        throw new Exception("1");
            //    for (int i = 0; i < affoff.Length; i++)
            //        if (affoff[i] != m_AffineOffset[i])
            //            throw new Exception("2");

            //    double[] r = (double[])bf.Deserialize(fs);
            //    if (r.Length != rhs.Length)
            //        throw new Exception("1");
            //    for (int i = 0; i < rhs.Length; i++)
            //        if (r[i] != rhs[i])
            //            throw new Exception("3");






            //    fs.Flush();
            //    fs.Close();
            //    Console.WriteLine("diff tested.");

            //}
        }

        ///// <summary>
        ///// debug code
        ///// </summary>
        //public void CheckParallel() {

        //    Comm.DatabaseDriver cm = m_Context.CommMaster;

        //    // communicate 
        //    Transceiver m_transciever = new Transceiver(m_Context.CommMaster, m_Mapping.Fields);
        //    m_transciever.TransceiveStartImReturn();
        //    m_transciever.TransceiveFinish();

        //    //
        //    double[] State = new double[m_Mapping.Count];
        //    m_Mapping.CopyTo(State,0);


        //    // evaluate the reference equations
        //    double[] k = new double[m_Mapping.LocalLength];

        //    m_RefQuadEdge.m_Output = k;
        //    m_RefQuadEdge.Execute();

        //    //m_RefQuadVol.m_Output = k;
        //    //m_RefQuadVol.Execute();

        //    // save compare value
        //    if (cm.Size == 1) {
        //        // save

        //        FileStream fsO = new FileStream("ref.dat",FileMode.Create);
        //        BinaryFormatter formatterO = new BinaryFormatter();

        //        formatterO.Serialize(fsO, State);
        //        formatterO.Serialize(fsO, k);
        //        formatterO.Serialize(fsO, m_Context.GridDat.CurrentGlobalIdPermutation.Values);

        //        fsO.Close();
        //        fsO.Dispose();

        //        return;
        //    }


        //    FileStream fs = new FileStream("ref.dat", FileMode.Open, FileAccess.Read, FileShare.Read);
        //    BinaryFormatter formatter = new BinaryFormatter();

        //    double[] __State = (double[])formatter.Deserialize(fs);
        //    double[] __k = (double[])formatter.Deserialize(fs);
        //    long[] __perm = (long[])formatter.Deserialize(fs);

        //    fs.Close();
        //    fs.Dispose();


        //    int i0 = (int)m_Context.GridDat.CurrentGlobalIdPermutation.i0Offset;
        //    int len = m_Context.GridDat.CurrentGlobalIdPermutation.LocalLength;

        //    if (cm.MyRank != 0)
        //        return;

        //    long[] _perm = new long[len];
        //    Array.Copy(__perm, i0, _perm, 0, len);

        //    i0 *= m_Mapping.TotalNoOfCoordinatesPerCell;
        //    len *= m_Mapping.TotalNoOfCoordinatesPerCell;

        //    double[] _State = new double[len];
        //    Array.Copy(__State, i0, _State, 0, len);

        //    double[] _k = new double[len];
        //    Array.Copy(__k, i0, _k, 0, len);



        //    // compare
        //    long[] perm = m_Context.GridDat.CurrentGlobalIdPermutation.Values;
        //    if (perm.Length != _perm.Length)
        //        throw new Exception("1");
        //    for (int i = 0; i < perm.Length; i++) {
        //        if (perm[i] != _perm[i])
        //            throw new Exception("a" + i);
        //    }


        //    if (_State.Length != m_Mapping.LocalLength)
        //        throw new Exception("2");
        //    if (_k.Length != k.Length)
        //        throw new Exception("3");

        //    //System.Diagnostics.Debugger.Break();


        //    BLAS.daxpy(_State.Length, -1.0, State, 1, _State, 1);
        //    double nrm_State = BLAS.dnrm2(_State.Length, _State, 1);
        //    Console.WriteLine("nrm_State = " + nrm_State);


        //    BLAS.daxpy(_k.Length, -1, k, 1, _k, 1);
        //    double nrm_k = BLAS.dnrm2(_k.Length, _k, 1);
        //    Console.WriteLine("nrm_k = " + nrm_k);

        //}



        ///// <summary>
        ///// 
        ///// </summary>
        ///// <param name="rhs"></param>
        ///// <returns>the L2 norm of the residual vector.</returns>
        //public double Residual(double[] rhs) {

        //    //throw new NotImplementedException();

        //    //// communicate 
        //    //Transceiver m_transciever = new Transceiver(m_Context.CommMaster, m_Mapping.Fields);
        //    //m_transciever.TransceiveStartImReturn();
        //    //m_transciever.TransceiveFinish();




        //    // evaluate the linear matrix
        //    double[] res = new double[m_Mapping.LocalLength];
        //    eM.Gemv<CoordinateVector, double[]>(1.0, this.m_DGCoordinates, 0, res);
        //    BLAS.daxpy(res.Length, 1.0, m_AffineOffset, 1, res, 1);

        //    // substract right hand side
        //    if (rhs != null) {
        //        if (rhs.Length < res.Length)
        //            throw new ArgumentException("length of right-hand-side vector is too short.");

        //        BLAS.daxpy(res.Length, -1.0, rhs, 1, res, 1);
        //    }


        //    double nrm2 = BLAS.dnrm2(res.Length, res, 1);

        //    double nrmInf = 0;
        //    for (int i = 0; i < res.Length; i++)
        //        if (Math.Abs(res[i]) > nrmInf) 
        //            nrmInf = Math.Abs(res[i]);


        //    Console.WriteLine("Residual L2 = " + nrm2 + " Linf = " + nrmInf );




        //    return nrm2;

        //}




        ///// <summary>
        ///// equal to <see cref="m_EqsSys"/>;
        ///// </summary>
        //MsrMatrix eM;


        /// <summary>
        /// see <see cref="Solver"/>
        /// </summary>
        protected ISparseSolver m_Solver;

        /// <summary>
        /// Sparse solver used to solve the linear system
        /// This proerty can be used to change some properies of the equation system solver
        /// </summary>
        public ISparseSolver Solver {
            get {
                return m_Solver;
            }
        }



        /// <summary>
        /// affine offset to the linear operator
        /// </summary>
        protected double[] m_AffineOffset;

        /// <summary>
        /// Row mapping - vector of unknowns;
        /// does alos apply to <see cref="m_AffineOffset"/>
        /// </summary>
        protected CoordinateMapping m_Mapping;


        /// <summary>
        /// the DG fields which are accected by this linear solver
        /// </summary>
        public CoordinateMapping Mapping {
            get {
                return m_Mapping;
            }
        }

        /// <summary>
        /// <see cref="DGCoordinates"/>
        /// </summary>
        protected CoordinateVector m_DGCoordinates;

        /// <summary>
        /// The DG coordinates of the <see cref="BoSSS.Foundation.DGField"/>'s in <see cref="Mapping"/>
        /// (see <see cref="BoSSS.Foundation.CoordinateMapping.Fields"/>);
        /// </summary>
        public CoordinateVector DGCoordinates {
            get {
                return m_DGCoordinates;
            }
        }

        /// <summary>
        /// solves the equation system;
        /// there is no solution output, 
        /// the solution is written to the fields in <see cref="Mapping"/>;
        /// </summary>
        /// <returns>the basic statistics of the sparse solver</returns>
        /// <remarks>
        /// If an right-hand-side operator has been specified, i.e.
        /// <see cref="RhsEvaluator"/> is not null, it is evaluated prior
        /// to the sparse solver call;
        /// </remarks>
        virtual public SolverResult Solve() {
            return Solve(null);
        }

        /// <summary>
        /// Tells a handler which right-hand-side is used as an input for the solver,
        /// during the invocation of <see cref="Solve()"/>;
        /// </summary>
        /// <remarks>
        /// This is useful if e.g. the
        /// <see cref="LinearSolver(ISparseSolver,SpatialOperator,CoordinateMapping,SpatialOperator,CoordinateMapping)"/>- 
        /// or the
        /// <see cref="LinearSolver(ISparseSolver,MsrMatrix,IList{double},CoordinateMapping,SpatialOperator, CoordinateMapping)"/>-constructor
        /// as been used to construct this object;
        /// In these cases, the <see cref="RhsEvaluator"/>-property is not null,
        /// and this event can be used to monitor the result of <see cref="RhsEvaluator"/>
        /// during the invocation of <see cref="Solve()"/>;
        /// </remarks>
        public event OnRhsEvaluatedHandler RHSEvaluated;

        /// <summary>
        /// assembles the right-hand-side (DG coordinates of <paramref name="rhsIn"/> minus the affine offset <see cref="m_AffineOffset"/>
        /// of the operator, which usually carries the boundary conditions)
        /// </summary>
        /// <param name="rhsIn">input</param>
        /// <param name="TotalRhs">output</param>
        void GetTotalRHS(IList<double> rhsIn, out double[] TotalRhs) {
            int n = m_Mapping.LocalLength;

            // check
            if (rhsIn != null)
                if (rhsIn.Count < m_Mapping.LocalLength)
                    throw new ArgumentOutOfRangeException("right hand side vector to short.");

            // clone or allocate memory 
            if (rhsIn == null) {
                TotalRhs = new double[m_AffineOffset.Length];
            } else {
                TotalRhs = rhsIn as double[];
                if (TotalRhs == null)
                    TotalRhs = rhsIn.ToArray();
            }
            rhsIn = null;
            //Console.WriteLine("Affine offset norm = " + BLAS.dnrm2(m_AffineOffset.Length,m_AffineOffset,1));

            // call evaluator
            if (m_rhsEvaluator != null)
                m_rhsEvaluator.Evaluate<double[]>(1.0, 1.0, TotalRhs);

            // notify
            if (RHSEvaluated != null)
                RHSEvaluated(TotalRhs);


            // substract affine part
            BLAS.daxpy(n, -1.0, m_AffineOffset, 1, TotalRhs, 1);
        }


        /// <summary>
        /// solves the equation system;
        /// there is no solution output, 
        /// the solution is written to the fields in <see cref="Mapping"/>;
        /// </summary>
        /// <param name="rhs">optional right-hand-side, can be null;</param>
        /// <returns>the basic statistics of the sparse solver</returns>
        /// <remarks>
        /// If an right-hand-side operator has been specified, i.e.
        /// <see cref="RhsEvaluator"/> is not null, it is evaluated prior
        /// to the sparse solver call;
        /// The result of this right-hand-side operator, added to <paramref name="rhs"/>,
        /// can be monitored by consuming the <see cref="RHSEvaluated"/>-event.
        /// </remarks>
        virtual public SolverResult Solve(IList<double> rhs) {

            // initial stuff
            // =============

            using (var tr = new ilPSP.Tracing.FuncTrace()) {


                // prepare right-hand-side
                // =======================

                // clone or allocate memory 
                double[] _rhs;
                GetTotalRHS(rhs, out _rhs);

                //// testcode
                //{
                //    Field f = new SinglePhaseField(m_Context, new Basis(m_Context, 2));

                //    int J = f.Coordinates.NoOfRows;
                //    int N = f.Coordinates.NoOfCols;
                //    int cnt = 0;
                //    for( int j = 0; j < J; j++)
                //        for (int nn = 0; nn < N; nn++) {
                //            f.Coordinates[j, nn] = _rhs[cnt];
                //            cnt++;
                //        }

                //    double InfNrm = f.LinfNorm();
                //    double TwoNorm = f.L2Norm();

                //    Console.WriteLine("rhs L2: " + TwoNorm + " inf: " + InfNrm);


                //    //Array.Copy(_rhs, f.Coordinates.C, _rhs.Length);
                //    BoSSS.Solution.Tecplot.Tecplot.PlotFields(new Field[] { f }, null, m_Context, "rhs", "---", 0, 0);
                //}

                // call solver
                // ===========

                tr.Info("entering linear solver:");
                SolverResult ret = m_Solver.Solve<CoordinateVector, double[]>(m_DGCoordinates, (double[])_rhs.Clone());
                tr.Info("no. of iterations: " + ret.NoOfIterations);
                tr.Info("converged? : " + ret.Converged);
                tr.Info("Pure solver runtime: " + ret.RunTime.TotalSeconds + " sec.");

                // finalize
                // ========
                return ret;
            }

        }
    }
}