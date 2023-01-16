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
//#define TEST

using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Kraypis;
using ilPSP.Tracing;
using System.IO;
using ilPSP.LinSolvers.PARDISO;

namespace BoSSS.Foundation.ConstrainedDGprojection {


    /// <summary>
    /// <see cref="ConstrainedDGFieldMk3.Factory"/>
    /// </summary>
    public enum ProjectionStrategy {

        /// <summary>
        /// <see cref="ConstrainedDGField_Global"/>
        /// </summary>
        globalOnly = 2,

        //globalWithPatchwisePrecond = 1,

        /// <summary>
        /// <see cref="ConstrainedDgField_Patchwise"/>
        /// </summary>
        patchwiseOnly = 0

        //patchwiseWithGlobalPrecond = 3

    }



    /// <summary>
    /// Projection of a DG field onto a continuous subspace of the DG space;
    /// Applications which require continuous approximations are:
    /// - Level Set algorithms; if a level-set field is not continuous, the zero-set is not connected, which can result in many pathological cases
    /// - Artificial Viscosity (AV): also AV fields must be continuous to provide decent results
    /// </summary>
    /// <remarks>
    /// - Explained in PhD thesis of M. Smuda, see https://tuprints.ulb.tu-darmstadt.de/17376/
    /// </remarks>
    abstract public class ConstrainedDGFieldMk3 : IDisposable {

        /// <summary>
        /// Factory to produce different variants.
        /// </summary>
        /// <param name="b"><see cref="ConstrainedDGFieldMk3.Basis"/></param>
        /// <param name="__domainLimit"><see cref="ConstrainedDGFieldMk3.domainLimit"/></param>
        /// <param name="s"></param>
        /// <returns>
        /// 
        /// </returns>
        static public ConstrainedDGFieldMk3 Factory(Basis b, CellMask __domainLimit, ProjectionStrategy s) {
            switch(s) {
                case ProjectionStrategy.globalOnly: return new ConstrainedDGField_Global(b, __domainLimit);
                case ProjectionStrategy.patchwiseOnly: return new ConstrainedDgField_Patchwise(b, __domainLimit);
                default: throw new NotImplementedException();
            }
        }



        /// <summary>
        /// ctor
        /// </summary>
        public ConstrainedDGFieldMk3(Basis b, CellMask __domainLimit) {
            m_Basis = b;
            m_grd = (GridData)b.GridDat;
            m_Mapping = new UnsetteledCoordinateMapping(b);
            m_Coordinates = new double[m_Mapping.LocalLength];
            this.internalProjection = new SinglePhaseField(m_Basis, "internalProjection");
            this.initialProjection0 = new SinglePhaseField(m_Basis, "initialProjection0");

            if(__domainLimit == null)
                __domainLimit = CellMask.GetFullMask(b.GridDat);
            if(__domainLimit.NoOfItemsLocally.MPISum() <= 0) {
                throw new ArgumentOutOfRangeException("Domain mask cannot be empty.");
            }
            this.domainLimit = __domainLimit;

        }

        /// <summary>
        /// release of internal solvers
        /// </summary>
        public abstract void Dispose();

        /// <summary>
        /// domain on which the projection is performed;
        /// null denotes the entire domain;
        /// in the case of level-set continuity projection, typically the narrow band.
        /// </summary>
        public CellMask domainLimit {
            get;
            private set;
        }

        readonly Basis m_Basis;

        /// <summary>
        /// 
        /// </summary>
        public Basis Basis {
            get {
                return m_Basis;
            }
        }

        readonly GridData m_grd;

        readonly UnsetteledCoordinateMapping m_Mapping;

        /// <summary>
        /// DG Coordinates of the current approximation;
        /// after execution of the approximation algorithm, hopefully continuous;
        /// </summary>
        protected double[] m_Coordinates;

        /// <summary>
        /// DG representation of the current solution
        /// </summary>
        protected readonly SinglePhaseField internalProjection;

        /// <summary>
        /// DG representation of the current solution;
        /// <see cref="CheckLocalProjection(out SinglePhaseField, CellMask, bool)"/> computes jump norms from this field
        /// </summary>
        public SinglePhaseField InternalDGfield {
            get {
                return internalProjection;
            }
        }

        /// <summary>
        /// DG coordinates of original discontinuous representation
        /// </summary>
        protected double[] m_Coordinates0;

        /// <summary>
        /// DG representation of the initial (discontinuous) solution
        /// </summary>
        protected readonly SinglePhaseField initialProjection0;


        /// <summary>
        /// hard-coded switch to turn some console output on/off (BAD PRACTICE)
        /// </summary>
        protected bool diagnosticOutput = false;

        /// <summary>
        /// hard-coded switch to turn debugging output for MATLAB on/off
        /// </summary>
        protected bool diagOutputMatlab = false;


        /// <summary>
        /// Projects some DG field <paramref name="orgDGField"/> onto the internal, continuous representation.
        /// </summary>
        /// <param name="orgDGField">
        /// input; unchanged on exit
        /// </param>
        abstract public void ProjectDGField(ConventionalDGField orgDGField);

        /// <summary>
        /// sets the internal DG coordinates from <paramref name="orgDGField"/>
        /// </summary>
        protected void SetDGCoordinatesOnce(DGField orgDGField, CellMask mask) {
            // get DG-coordinates (change of basis for projection on a higher polynomial degree)
            foreach(int j in mask.ItemEnum) {
                int N = Math.Min(orgDGField.Basis.GetLength(j), this.m_Basis.GetLength(j));
                for(int n = 0; n < N; n++) {
                    m_Coordinates[m_Mapping.LocalUniqueCoordinateIndex(0, j, n)] = orgDGField.Coordinates[j, n];
                }
            }

            m_Coordinates0 = m_Coordinates.CloneAs();
            m_Coordinates.ClearEntries();
            initialProjection0.CoordinateVector.Acc(1.0, m_Coordinates0);

            UpdateInternalProjection(csMPI.Raw._COMM.WORLD);
            //ProjectionSnapshots.Add(initialVal);
        }

        /// <summary>
        /// Core routines, i.e. projection on a certain patch <see cref="Patch"/>
        /// This is implemented as a class, 
        /// so that the factorization of the linear solver can be cached. (In the case of a direct solver)
        /// Then the first call to <see cref="PerformProjection"/> might be expensive,
        /// but any subsequent call is comparatively cheap.
        /// </summary>
        protected class ConstrainedProjectionInternal : IDisposable {

            /// <summary>
            /// 
            /// </summary>
            /// <param name="runOnlyLocal">
            /// - true: work on MPI SELF communicator
            /// - false: work on MPI WORLD communicator
            /// </param>
            /// <param name="__patch">
            /// <see cref="Patch"/>
            /// </param>
            /// <param name="__domain">
            /// <see cref="DomainLimit"/>
            /// </param>
            /// <param name="__owner"></param>
            public ConstrainedProjectionInternal(ConstrainedDGFieldMk3 __owner, CellMask __domain, CellMask __patch, bool runOnlyLocal) {
                if(__domain != null) {
                    if(__domain.MaskType != MaskType.Logical)
                        throw new ArgumentException();
                }


                if(runOnlyLocal) {
                    comm = MPI.Wrappers.csMPI.Raw._COMM.SELF;
                } else {
                    comm = MPI.Wrappers.csMPI.Raw._COMM.WORLD;
                }
                IsLocal = runOnlyLocal;
                this.owner = __owner;
                this.DomainLimit = __domain;
                this.Patch = __patch;
                this.Initialize();

            }

            readonly ConstrainedDGFieldMk3 owner;

            UnsetteledCoordinateMapping m_Mapping => owner.m_Mapping;
            GridData m_grd => owner.m_grd;
            Basis m_Basis => owner.m_Basis;


            /// <summary>
            /// Cell mask on which the projection will be performed;
            /// only cells within this mask will be updated.
            /// </summary>
            public CellMask Patch {
                get;
                private set;
            }

            /// <summary>
            /// Optional domain limit, can be null;
            /// on the boundary of this mask, no constraints will be enforced.
            /// If not specified, only the boundaries of the grid are considered.
            /// </summary>
            public CellMask DomainLimit {
                get;
                private set;
            }

            public MPI_Comm comm {
                get;
                private set;
            }

            public int MPIsize {
                get {
                    csMPI.Raw.Comm_Size(this.comm, out int sz);
                    return sz;
                }
            }

            public bool IsLocal {
                get;
                private set;
            }

            public void Dispose() {
                if(solver != null)
                    solver.Dispose();
            }

            NodeSet[] m_ConstrainNodes = null;


            /// <summary>
            /// Nodes at which the continuity constrains are enforced.
            /// - index: correlates with edge reference element index, <see cref="IGeometricalEdgeData.EdgeRefElements"/>, <see cref="IGeometricalEdgeData.GetRefElementIndex"/>
            /// </summary>
            NodeSet[] ConstrainNodes {
                get {
                    if(m_ConstrainNodes == null) {
                        var EdgRefElems = m_grd.iGeomEdges.EdgeRefElements;
                        int deg = m_Basis.Degree;
                        int Np;
                        switch(m_grd.SpatialDimension) {
                            case 1: Np = 1; break;
                            case 2: Np = deg + 1; break;
                            case 3:
                                // remark: (deg + 1) should be sufficient, but it seems that for numerical accuracy, we need (deg + 2)
                            Np = (deg + 2).ForLoop(i => i + 1).Sum(); // DG basis vector space dimension in 2D == number of constraints on face required.
                            //Np = m_Basis.Length;
                            break;
                            default: throw new NotImplementedException("unknown spatial dimension");
                        }

                        m_ConstrainNodes = new NodeSet[EdgRefElems.Length];
                        for(int iEref = 0; iEref < m_ConstrainNodes.Length; iEref++) {
                            int p = 1;
                            do {
                                //EdgRefElems[iEref].GetInterpolationNodes

                                m_ConstrainNodes[iEref] = EdgRefElems[iEref].GetQuadratureRule(p).Nodes;
                                EdgRefElems[iEref].GetNodeSet(p + 1, out m_ConstrainNodes[iEref], out _, out _);
                                
                                p++;
                            } while(m_ConstrainNodes[iEref].NoOfNodes < Np);
                        }
                    }
                    return m_ConstrainNodes;
                }
            }

            NodeSet GetNodeSet(int iEdg) {
                // we assume that logical and geometrical edges are equal, otherwise we throw some exception.
                if(m_grd.iLogicalEdges.EdgeToParts != null) {
                    if(m_grd.iLogicalEdges.EdgeToParts[iEdg] != null && m_grd.iLogicalEdges.EdgeToParts[iEdg][0] != iEdg) {
                        throw new NotImplementedException("missing implementation for aggregate meshes");
                    }
                }
                int iEdgeRef = m_grd.iGeomEdges.GetRefElementIndex(iEdg);
                return this.ConstrainNodes[iEdgeRef];
            }

            BitArray DomainLimitBitMask, PatchBitMask;


            void Initialize() {
                using(new FuncTrace()) {

                    EdgeMask fullEM;
                    if(IsLocal) {
                        fullEM = Patch.GetAllLocalEdgesMask();
                    } else {
                        SubGrid maskSG = new SubGrid(Patch);
                        fullEM = maskSG.AllEdgesMask;
                    }
                    EdgeMask constraintsMask = fullEM;

                    if(this.IsLocal) {
                        DomainLimitBitMask = this.DomainLimit?.GetBitMask();
                        PatchBitMask = this.Patch.GetBitMask();
                    } else {
                        DomainLimitBitMask = this.DomainLimit?.GetBitMaskWithExternal();
                        PatchBitMask = this.Patch.GetBitMaskWithExternal();
                    }


                    // assemble matrix
                    // =============== 

                    List<long> BlockI0 = new List<long>();
                    List<int> BlockLen = new List<int>();
                    long i0 = 0;
                    int nodeCount = 0;
                    foreach(int iEdg in constraintsMask.ItemEnum) {
                        if(ConsiderEdge(iEdg, out _, out _)) {

                            int NoOfNodes = GetNodeSet(iEdg).NoOfNodes;

                            BlockI0.Add(i0);
                            BlockLen.Add(NoOfNodes);
                            i0 += NoOfNodes;
                            nodeCount += NoOfNodes;
                        }
                    }

                    IBlockPartitioning rowBlockPart = new BlockPartitioning(nodeCount, BlockI0, BlockLen, this.comm, true);
                    IBlockPartitioning colBlockPart;
                    if(IsLocal) {
                        colBlockPart = m_Mapping.GetLocalBlockPartitioning();
                    } else {
                        colBlockPart = m_Mapping;
                    }
                    A = new BlockMsrMatrix(rowBlockPart, colBlockPart);
                    assembleConstrainsMatrix(A, constraintsMask);



                    // test with matlab
                    if(owner.diagOutputMatlab) {
                        MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
                        Console.WriteLine("Calling MATLAB/Octave...");
                        using(BatchmodeConnector bmc = new BatchmodeConnector()) {
                            bmc.PutSparseMatrix(A, "A");
                            bmc.Cmd("rank_A = rank(full(A))");
                            bmc.Cmd("rank_AT = rank(full(A'))");
                            bmc.GetMatrix(output, "[rank_A, rank_AT]");

                            bmc.Execute(false);
                        }

                        Console.WriteLine("A: No of Rows = {0}; rank = {1}", A.NoOfRows, output[0, 0]);
                        Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);
                    }

                    // initialize solver
                    // =================

                    AT = A.Transpose();
                    BlockMsrMatrix AAT = BlockMsrMatrix.Multiply(A, AT);
                    AAT.AssumeSymmetric = true;

                    InitSolver(AAT);
                }
            }


            /// <summary>
            /// Transpose of <see cref="A"/>
            /// </summary>
            BlockMsrMatrix AT;

            /// <summary>
            /// Constraints matrix
            /// </summary>
            BlockMsrMatrix A;

            /// <summary>
            /// internal matrix used by the linear solver
            /// </summary>
            IMutableMatrixEx solverAAT;

            public double PerformProjection() {
                using(var tr = new FuncTrace()) {
                    // solve system
                    // ============

                    var DgCoordinates0 = owner.m_Coordinates0;
                    var DgCoordinates = owner.m_Coordinates;

                    double[] RHS = new double[A.RowPartitioning.LocalLength];

                    double[] remainder = DgCoordinates0.CloneAs();
                    remainder.AccV(-1.0, DgCoordinates);

                    A.SpMV(1.0, remainder, 0.0, RHS);

                    double[] v = new double[RHS.Length]; // Lagrange multiplyer
                    double[] x = new double[DgCoordinates.Length];

                    {
                        CheckRHS(solverAAT, RHS);
                        solver.Solve(v, RHS); // solve for the Lagrange multiplyer
                        CheckSolutionResidual(tr, RHS, v);
                    }
                    // x = RHS - ATv
                    AT.SpMV(-1.0, v, 0.0, x);
                    x.AccV(1.0, remainder);

                    // x must be orthogonal to (x-remainder)
                    double[] R = remainder.CloneAs();
                    R.AccV(-1.0, x);
                    double tst = x.InnerProd(R).MPISum(comm);
                    tr.Info("Inner Product Test (should be 0.0): " + tst);


                    double l2_change = 0.0;
                    foreach(int j in Patch.ItemEnum) {
                        int Np = m_Basis.GetLength(j);
                        int ii0 = m_Mapping.LocalUniqueCoordinateIndex(0, j, 0);
                        for(int n = 0; n < Np; n++) {
                            double x_i = x[ii0 + n];
                            l2_change += x_i * x_i;
                            DgCoordinates[ii0 + n] += x_i;
                        }
                    }
                    l2_change = l2_change.MPISum(comm).Sqrt();


                    SinglePhaseField errorField;
                    (double jumpNorm, double L2err) = owner.CheckLocalProjection(this.comm, out errorField, this.Patch, true);
                    //if(owner.diagnosticOutput)
                    //    Console.WriteLine("L2 jump norm on mask: {0}", jumpNorm);
                    //if(jumpNorm > 1e-11) {
                    //    if(owner.diagnosticOutput) {
                    //        Console.WriteLine("======================");
                    //        Console.WriteLine("project mask: No of cells {0}", this.Patch.NoOfItemsLocally);
                    //    }

                    //}

                    return l2_change;
                }
            }


            private static void CheckRHS(IMutableMatrixEx solverAAT, double[] RHS) {

                // basic check on the validity of the RHS
                long __i0 = solverAAT.RowPartitioning.i0;
                int __L = solverAAT.RowPartitioning.LocalLength;

                for(int l = 0; l < __L; l++) {
                    long i = __i0 + l;
                    int nnz_i = solverAAT.GetNoOfNonZerosPerRow(i);

                    if(nnz_i <= 0) {
                        if(RHS[i] != 0)
                            throw new ArithmeticException("illegal system; nonzero RHS in absolute zero row");
                    }
                }
            }

            /// <summary>
            /// Verifies that the solution of the linear system is sufficiently good.
            /// </summary>
            void CheckSolutionResidual(FuncTrace tr, double[] RHS, double[] v) {

                // check solver residual

                double[] Resi = RHS.CloneAs();
                solverAAT.SpMV(-1.0, v, 1.0, Resi);
                double resi_l2, rhs_l2, v_l2;
                rhs_l2 = RHS.MPI_L2Norm(solverAAT.MPI_Comm);
                resi_l2 = Resi.MPI_L2Norm(solverAAT.MPI_Comm);
                v_l2 = Resi.MPI_L2Norm(solverAAT.MPI_Comm);
                double AAT_inf = solverAAT.InfNorm();

                double Denom = Math.Max(AAT_inf, Math.Max(rhs_l2, Math.Max(v_l2, Math.Sqrt(BLAS.MachineEps))));
                double RelResidualNorm = resi_l2 / Denom;

                tr.Info("solver residual " + RelResidualNorm);
                if(RelResidualNorm > 1.0e-5) {

                    //solverAAT.SaveToTextFileSparse($"Mtx{blabla}.txt");
                    //AAT.SaveToTextFileSparse($"MtxOrg{blabla}.txt");
                    //v.SaveToTextFile($"V{blabla}.txt");
                    //RHS.SaveToTextFile($"RHS{blabla}.txt");
                    //blabla++;

                    string solverName = this.solver.GetType().FullName;

                    string ErrMsg;
                    using(var stw = new StringWriter()) {
                        stw.WriteLine("High residual from solver (using {0}).", solverName);
                        stw.WriteLine("    L2 Norm of RHS:         " + rhs_l2);
                        stw.WriteLine("    L2 Norm of Solution:    " + v_l2);
                        stw.WriteLine("    L2 Norm of Residual:    " + resi_l2);
                        stw.WriteLine("    Relative Residual norm: " + RelResidualNorm);
                        stw.WriteLine("    Matrix Inf norm:        " + AAT_inf);

                        stw.WriteLine("Dumping text versions of Matrix, Solution and RHS.");
                        ErrMsg = stw.ToString();
                    }

                    throw new ArithmeticException(ErrMsg);
                }
            }

            ISparseSolver solver {
                get;
                set;
            }

            private void InitSolver(BlockMsrMatrix matrix) {
                if(solver != null)
                    throw new NotSupportedException("can only be called once");
                this.solverAAT = matrix;

                //PARDISOSolver _solver = new PARDISOSolver();
                //this.solver = _solver;
                //_solver.SymmIndefPivot = true;
                //_solver.CacheFactorization = true;
                //_solver.Parallelism = IsLocal ? Parallelism.SEQ : Parallelism.OMP;

                long DOF = matrix.RowPartitioning.TotalLength;
                var _solver = IsLocal ? SolverUtils.PatchSolverFactory() : SolverUtils.GlobalSolverFactory(DOF);

                _solver.DefineMatrix(matrix);
                this.solver = _solver;
                //if(_solver. != this.comm)
                //    throw new ApplicationException();
            }


            NodeSet GetVerticeRefNodeSet(int vertInd, int jCell) {

                double[] vertCoord = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(vertInd, -1).To1DArray();

                int D = vertCoord.Length;
                MultidimensionalArray vertCoord_glb = MultidimensionalArray.Create(1, D);
                vertCoord_glb.SetRow<double[]>(0, vertCoord);

                MultidimensionalArray vertCoord_loc = MultidimensionalArray.Create(1, 1, D);
                m_grd.TransformGlobal2Local(vertCoord_glb, vertCoord_loc, jCell, 1, 0);
                double[] vertCellCoord1 = vertCoord_loc.ExtractSubArrayShallow(0, 0, -1).To1DArray();

                return new NodeSet(m_grd.Cells.RefElements[0], vertCellCoord1, false);

            }

            /*
             * Remark: non-zero boundary conditions will never work, i.e. produce a constraint system which has no solution;
             * E.g., consider a Square (0,1)^2 and b.c. at three edges: u = 0 @ x = 0, u = x @ y = 1, u = -x @ y = -1; there is no solution for the Ansatz u = a + x*b + y*c;
             * (Fk, 17aug21)
             */


            void assembleConstrainsMatrix(BlockMsrMatrix A, EdgeMask constrainsEdges) {
                MPICollectiveWatchDog.Watch(this.comm);

                long ColumnOffset;
                if(IsLocal)
                    ColumnOffset = -m_Mapping.i0; // translate back to the partitioning living on MPI_SELF 
                else
                    ColumnOffset = 0;


                int Jup = this.m_grd.iLogicalCells.NoOfLocalUpdatedCells;

                int count = 0;
                long _nodeCount = A.RowPartitioning.i0; // offset into the row
                foreach(int edg in constrainsEdges.ItemEnum) {

                    //var edgeInfo = m_grd.Edges.Info[edg];
                    //if(edgeInfo.HasFlag(EdgeInfo.Interprocess) && !ownedInterProcEdges.Contains(edg))
                    //    continue;


                    if(ConsiderEdge(edg, out int cell1, out int cell2)) {



                        // Omit any cells outside of the patch:
                        // setting Np1 or Np2 to 0 will enforce a homogeneous boundary
                        int Np1, Np2;
                        Np1 = PatchBitMask[cell1] ? this.m_Basis.GetLength(cell1) : 0; //                 omit any cells outside of the patch:

                        if(IsLocal && cell2 >= Jup) {
                            // local patch, external cell: must enforce homogeneous boundary
                            Np2 = 0;
                        } else {
                            Np2 = PatchBitMask[cell2] ? this.m_Basis.GetLength(cell2) : 0;
                        }

                        if(Np1 <= 0 && Np2 <= 0)
                            throw new ApplicationException("error in algorithm");

                        // set continuity constraints
                        NodeSet qNodes = GetNodeSet(edg);
                        var results = m_Basis.EdgeEval(qNodes, edg, 1);

                        for(int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                            // Cell1   

                            for(int p = 0; p < Np1; p++) {
                                A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p) + ColumnOffset] = results.Item1[0, qN, p];
                            }
                            // Cell2
                            for(int p = 0; p < Np2; p++) {
                                A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p) + ColumnOffset] = -results.Item2[0, qN, p];
                            }

                            count++;
                        }
                        _nodeCount += qNodes.NoOfNodes;
                    }
                }


                if(count != A.RowPartitioning.LocalLength)
                    throw new ArgumentException("Error in algorithm");
            }





            bool ConsiderEdge(int edg, out int cell1, out int cell2) {
                cell1 = m_grd.Edges.CellIndices[edg, 0];
                cell2 = m_grd.Edges.CellIndices[edg, 1];


                if(cell2 < 0)
                    return false; // don't apply any constraint on a boundary edge

                if(DomainLimitBitMask != null) {
                    if(!DomainLimitBitMask[cell1])  // don't apply any constraint at the boundary of the domain limit
                        return false;
                    if(cell2 < DomainLimitBitMask.Length && !DomainLimitBitMask[cell2])
                        return false;
                }

                return true;
            }
        }



        /// <summary>
        /// 
        /// </summary>
        protected void UpdateInternalProjection(MPI_Comm comm) {
            MPICollectiveWatchDog.Watch();
            //Console.WriteLine("update internal projection field");
            internalProjection.Clear();
            int stride = m_Mapping.MaxTotalNoOfCoordinatesPerCell;
            internalProjection._Acc(1.0, m_Coordinates, 0, stride, true);
            if(comm == csMPI.Raw._COMM.WORLD) {
                internalProjection.MPIExchange();
            } else if(comm == csMPI.Raw._COMM.SELF) {
                // noop
            } else {
                throw new NotImplementedException("only supported for WORLD and SELF communicator");
            }
        }

        /// <summary>
        /// Computes the norm of jumps on the interior edges of the <paramref name="mask"/>
        /// L2-error norm of the projected field against the initial discontinuous field 
        /// </summary>
        public (double jumpNorm, double L2err) CheckLocalProjection(out SinglePhaseField errField, CellMask mask = null, bool onInterProc = false) {
            return CheckLocalProjection(csMPI.Raw._COMM.WORLD, out errField, mask, onInterProc);
        }

        /// <summary>
        /// Computes the norm of jumps on the interior edges of the <paramref name="mask"/>
        /// </summary>
        public double CheckLocalProjection(CellMask mask = null, bool onInterProc = false) {
            SinglePhaseField errField;
            (double jumpNorm, double L2err) = CheckLocalProjection(csMPI.Raw._COMM.WORLD, out errField, mask, onInterProc);
            return jumpNorm;
        }

        /// <summary>
        /// A hack only used for testing
        /// </summary>
        public void ScheissDrauf(SinglePhaseField f) {

            this.SetDGCoordinatesOnce(f, CellMask.GetFullMask(InternalDGfield.GridDat));
            InternalDGfield.Clear();
            InternalDGfield.AccLaidBack(1.0, f);

            this.m_Coordinates0.ClearEntries();
            this.m_Coordinates0.AccV(1.0, InternalDGfield.CoordinateVector);

            this.m_Coordinates.ClearEntries();
            this.m_Coordinates.AccV(1.0, InternalDGfield.CoordinateVector);
        }

        /// <summary>
        /// Computes the norm of jumps on the interior edges of the <paramref name="mask"/>
        /// </summary>
        protected (double jumpNorm, double L2err) CheckLocalProjection(MPI_Comm comm, out SinglePhaseField errField, CellMask mask = null, bool onInterProc = false) {
            MPICollectiveWatchDog.Watch(comm);

            bool IsLocal = IsMPIself(comm);
            EdgeMask innerEM;
            if(IsLocal) {
                innerEM = mask.GetAllInnerEdgesMask();
            } else {
                if(mask == null) {
                    mask = CellMask.GetFullMask(m_grd);
                }
                SubGrid maskSG = new SubGrid(mask);
                innerEM = maskSG.InnerEdgesMask;
            }

            UpdateInternalProjection(comm);

            // check jump norm on interior edges
            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, m_grd,
                (new EdgeQuadratureScheme(true, innerEM)).Compile(m_grd, m_Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for(int j = 0; j < Length; j++) {
                        int iEdge = j + i0;

                        var edgeInfo = m_grd.Edges.Info[iEdge];
                        if(!onInterProc && edgeInfo.HasFlag(EdgeInfo.Interprocess))
                            continue;

                        int jCell_IN = m_grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = m_grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if(jCell_OT >= 0) {

                            int iTrafo_IN = m_grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = m_grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(m_grd, iTrafo_IN, false);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(m_grd, iTrafo_OT, false);

                            internalProjection.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            internalProjection.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);

                            //if (uDiff.L2Norm() > 1e-9)
                            //    Console.WriteLine("uDiff at edge {0} between cell {1} and cell {2}: {3}", iEdge, jCell_IN, jCell_OT, uDiff.L2Norm());
                        } else {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum(comm);

            // check projection against initial projection
            errField = initialProjection0.CloneAs();
            errField.AccLaidBack(-1.0, internalProjection, mask);

            double L2err = errField.L2Norm(mask);

            // check the single cells within the projection (Debugging)
            //foreach(int cell in mask.ItemEnum) {
            //    BitArray oneBA = new BitArray(m_Basis.GridDat.Grid.NumberOfCells);
            //    oneBA[cell] = true;
            //    CellMask oneCM = new CellMask(m_Basis.GridDat, oneBA);
            //    double oneErr = errField.L2Norm(oneCM);
            //    if(oneErr > 1e-4)
            //        Console.WriteLine("cell {0}: error norm = {1}", cell, oneErr);
            //}


            return (Unorm.Sqrt(), L2err);

        }

        private static bool IsMPIself(MPI_Comm comm) {
            bool IsLocal;
            if(comm == csMPI.Raw._COMM.WORLD) {
                IsLocal = false;
            } else if(comm == csMPI.Raw._COMM.SELF) {
                IsLocal = true;
            } else {
                throw new NotSupportedException("Only supported for WORLD ans SELF communicator.");
            }

            return IsLocal;
        }

        private class RuntimeTracker : IDisposable {
            public RuntimeTracker(string recordname, bool activate) {
                m_activate = activate;
                m_RuntimeSTW = new Stopwatch();
                m_record = RuntimeRecord;
                m_name = recordname;
                if(m_activate) m_RuntimeSTW.Start();
            }

            private Stopwatch m_RuntimeSTW;
            private bool m_activate;
            private Dictionary<string, double> m_record;
            private string m_name;

            private double GetTime() {
                double thisreturn = 0;
                m_RuntimeSTW.Stop();
                thisreturn = m_RuntimeSTW.Elapsed.TotalSeconds;
                return thisreturn;
            }

            private bool AddRecordIfNew() {
                double time = GetTime();
                bool IsNew = m_record.ContainsKey(m_name) == false;
                if(IsNew)
                    m_record.Add(m_name, time);
                else {
                    m_record[m_name] += time;
                }
                return IsNew;
            }

            public void Dispose() {
                if(m_activate) AddRecordIfNew();
                m_RuntimeSTW = null;
            }
        }

        static private Dictionary<string, double> RuntimeRecord = new Dictionary<string, double>();

        public void PrintRuntimes() {
            foreach(var kv in RuntimeRecord) {
                Console.WriteLine(kv.Key + " : " + kv.Value);
            }
        }


        /// <summary>
        /// Accumulate the internal continuous field to a DG Field
        /// </summary>
        /// <param name="alpha">Scaling factor</param>
        /// <param name="DGField">output</param>
        /// <param name="mask"></param>
        public void AccToDGField(double alpha, ConventionalDGField DGField, CellMask mask = null) {
            if(!DGField.Basis.Equals(this.m_Basis))
                throw new ArgumentException("Basis does not match.");
            if(mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }
            if(mask.MaskType != MaskType.Logical)
                throw new ArgumentException("Expecting logical cell mask, but mask type is " + mask.MaskType);

            foreach(int j in mask.ItemEnum) {
                // collect coordinates for cell j
                int N = DGField.Basis.GetLength(j);
                for(int n = 0; n < N; n++) {
                    double coord_jn = m_Coordinates[m_Mapping.LocalUniqueCoordinateIndex(0, j, n)];
                    DGField.Coordinates[j, n] += coord_jn * alpha;
                }
            }
        }
    }
}
