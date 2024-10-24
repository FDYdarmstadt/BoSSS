using BoSSS.Application.BoSSSpad;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace StokesHelical_Ak {

    public class R0fix {

        public IGridData GridDat => m_map.GridDat;

        UnsetteledCoordinateMapping m_map;

        /// <summary>
        /// The prolongation operator, performs an injection from the 
        /// R0-restricted DG space to the standard DG space.
        /// </summary>
        BlockMsrMatrix Trafo;

        /// <summary>
        /// The restriction operator, performs a projection from the 
        /// standard DG space to the R0-restricted DG space.
        /// </summary>
        BlockMsrMatrix TrafoTranspose;


        CellMask R0BndyCells;
        int[] UnusedDofs_local;
        int[] UnusedDofs_global;

        public R0fix(UnsetteledCoordinateMapping map, double rMin) {
            //public R0fix(UnsetteledCoordinateMapping map, R0condType[] condType) {
            //if(map.BasisS.Count != condType.Length)
            //    throw new ArgumentException("mismatch between number of DG basis objects and number of R0-condition specifiers.");
            if(rMin < 10e-6) {
                m_map = map;
                (Trafo, UnusedDofs_local, UnusedDofs_global) = CreateTrafoMatrix(map);
                TrafoTranspose = Trafo.Transpose();
                R0BndyCells = GetR0BndyCells(map.GridDat);
            }
        }


        internal void InternalChecks() {

            Trafo.SaveToTextFileSparse("Prolo.txt");
            TrafoTranspose.SaveToTextFileSparse("Restr.txt");

        }


        /* Fk, 11apr24: deactivated, because seemingly un-used

        public R0fix(UnsetteledCoordinateMapping map, double rMin, bool old) {
            //public R0fix(UnsetteledCoordinateMapping map, R0condType[] condType) {
             //if(map.BasisS.Count != condType.Length)
             //    throw new ArgumentException("mismatch between number of DG basis objects and number of R0-condition specifiers.");
             if(rMin < 10e-6) {
                 m_map = map;
                 if(old == true) {
                     (Trafo, UnusedDofs_local, UnusedDofs_global) = CreateTrafoMatrix2(map);
                     TrafoTranspose = Trafo.Transpose();
                     R0BndyCells = GetR0BndyCells(map.GridDat);
                 } else {
                     (Trafo, UnusedDofs_local, UnusedDofs_global) = CreateTrafoMatrix(map);
                     TrafoTranspose = Trafo.Transpose();
                     R0BndyCells = GetR0BndyCells(map.GridDat);
                 }
             }
         }

        */
        public static MultidimensionalArray GetImplicitCondition_Value0(Basis B) {
            if(B.GridDat.SpatialDimension != 2)
                throw new NotSupportedException();
            if(B.GridDat.iGeomCells.GetCellType(0) != CellType.Square_Linear)
                throw new NotSupportedException();


            int p = B.Degree;

            // ============================================
            NodeSet TestNodes = new NodeSet(Square.Instance, p + 1, 2, false);
            double[] yNodes = GenericBlas.Linspace(-1, 1, p + 1);
            for(int k = 0; k < yNodes.Length; k++) {
                TestNodes[k, 0] = -1.0;
                TestNodes[k, 1] = yNodes[k];
            }
            TestNodes.LockForever();
            // return the y-Values at the test nodes
            // =========================================
            // Value0 Method
            var B_Imp = B.Evaluate(TestNodes);
            return B_Imp.ExtractSubArrayShallow(-1, -1);
        }


        public static MultidimensionalArray GetImplicitCondit_Gradient0_Stage1(Basis B) {
            if(B.GridDat.SpatialDimension != 2)
                throw new NotSupportedException();
            if(B.GridDat.iGeomCells.GetCellType(0) != CellType.Square_Linear)
                throw new NotSupportedException();
            int p = B.Degree;
            // ============================================
            if(p <= 1) {
                NodeSet TestNodes = new NodeSet(Square.Instance, 1, 2, false);
                TestNodes[0, 0] = -1.0; //x Coordinate
                TestNodes[0, 1] = 0.0;  //y Coordinate
                TestNodes.LockForever();
                var gradB = B.EvaluateGradient(TestNodes);
                return gradB.ExtractSubArrayShallow(-1, -1, 1).CloneAs().GetSolutionSpace();
            } else {
                NodeSet TestNodes = new NodeSet(Square.Instance, p, 2, false);
                double[] yNodes = GenericBlas.Linspace(-1, 1, p);
                for(int k = 0; k < yNodes.Length; k++) {
                    TestNodes[k, 0] = -1.0;
                    TestNodes[k, 1] = yNodes[k];
                }
                TestNodes.LockForever();
                var gradB = B.EvaluateGradient(TestNodes);
                return gradB.ExtractSubArrayShallow(-1, -1, 1).CloneAs().GetSolutionSpace();
            }
        }

        public static MultidimensionalArray GetImplicitCondit_Gradient0_Stage2(Basis B1) {
            if(B1.GridDat.SpatialDimension != 2)
                throw new NotSupportedException();
            if(B1.GridDat.iGeomCells.GetCellType(0) != CellType.Square_Linear)
                throw new NotSupportedException();
            // ============================================
            NodeSet CenterNode = new NodeSet(Square.Instance, 1, 2, false);
            CenterNode[0, 0] = -1.0; //x Coordinate
            CenterNode[0, 1] = 0.0;  //y Coordinate
            CenterNode.LockForever();
            // return the y-derivative at the test nodes
            // =========================================
            // Gradient0 Method Stage 2:
            var Same_Level_Cell0 = B1.Evaluate(CenterNode);

            var R1_local = GetImplicitCondit_Gradient0_Stage1(B1);

            var restricted_Cell0 = MultidimensionalArray.Create(1, R1_local.GetLength(1));
            restricted_Cell0.Multiply(1.0, Same_Level_Cell0, R1_local, 0.0, "kj", "ki", "ij");
            return restricted_Cell0.ExtractSubArrayShallow(-1, -1).CloneAs().GetSolutionSpace();
        }


        public static MultidimensionalArray GetTransformationMatrixValue0(Basis B) {
            var M = GetImplicitCondition_Value0(B); // implicit definition of the reduced DG basis
            var S = M.GetSolutionSpace(); //          convert to explicit representation

            int ExpectedRank = B.Degree + 1;
            if(S.NoOfCols != M.NoOfCols - ExpectedRank)
                throw new ArithmeticException("Gauss elimination failed - check tolerance in row echelon form.");

            return M;
        }


        // This Method gives you the Boundary Celly at r_min=0!
        public static CellMask GetR0BndyCells(IGridData g) {
            if(g.SpatialDimension != 2)
                throw new NotSupportedException();
            if(!object.ReferenceEquals(g.iGeomEdges.EdgeRefElements.Single(), Line.Instance))
                throw new NotSupportedException();


            var GrdDatTmp = ((BoSSS.Foundation.Grid.Classic.GridData)g);
            int D = g.SpatialDimension;
            int J = g.iLogicalCells.NoOfLocalUpdatedCells;  // Number of Cells?

            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.Create(2, D);

            var R0BndyCellsBitmask = new BitArray(J);
            int NoOfEdges = GrdDatTmp.Edges.Count; //Number of Inner Cells
            for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                if(GrdDatTmp.Edges.IsEdgeBoundaryEdge(iEdge)) {
                    int jCell = GrdDatTmp.Edges.CellIndices[iEdge, 0];  // Cell Index
                    int iFace = GrdDatTmp.Edges.FaceIndices[iEdge, 0];  // Face 

                    var Cell_j = GrdDatTmp.Cells.GetCell(jCell);      // CellInfosExtrahieren, Nodes, Laenge, Breite
                    var KRef = GrdDatTmp.Cells.GetRefElement(jCell);  // Referenz Zelle Groete

                    GrdDatTmp.TransformLocal2Global(KRef.GetFaceVertices(iFace), GlobalVerticesOut, jCell); // Transofrmation to Global

                    double x1 = GlobalVerticesOut[0, 0];
                    double x2 = GlobalVerticesOut[1, 0];


                    if(Math.Abs(x1) + Math.Abs(x2) < 1.0e-8) {
                        // this is an inner edge
                        R0BndyCellsBitmask[jCell] = true;
                    }
                }
            }

            var R = new CellMask(g, R0BndyCellsBitmask);
            return R;
        }

        public static (BlockMsrMatrix, int[], int[]) CreateTrafoMatrix(UnsetteledCoordinateMapping map) {
            if(map.CellDepLength)
                throw new NotSupportedException();
            var gdat = map.GridDat;
            if(map.BasisS.Count != 4)
                throw new ArgumentException("Expecting 4 DG fields: [ur, ux, ueta, pressure], but got " + map.BasisS.Count);


            // find cell at the inner, r = x = 0 array
            // =======================================
            var R0BndyCells = GetR0BndyCells(gdat);

            // get transformation matrix for each basis in the mapping
            // =======================================================

            var BasisS = map.BasisS.ToArray();
            int GAMMA = BasisS.Length;

            var TrafoMtxS = new MultidimensionalArray[GAMMA];
            var TrafoMtxSstage2 = new MultidimensionalArray[GAMMA];

            var NpS = new int[GAMMA];
            var UnusedDofs_local = new List<int>();
            var UnusedDofs_global = new List<int>();
            int NOff = 0;
            for(int iB = 0; iB < GAMMA; iB++) {  // loop over dependent variables in the PDE...
                                                 // Distinguish here between the different coefficients u_eta_;u_r, U_xi, p!
                switch(iB) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // variables ur, uxi: should be absolutely zero
                    // this can be enforced locally, i.e., in each cell
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    case 0:
                    case 1: {
                        var getTransformationMatrixValue0 = GetTransformationMatrixValue0(BasisS[iB]);
                        TrafoMtxS[iB] = getTransformationMatrixValue0.GetSolutionSpace();
                        break;
                    }
                    // ++++++++++++++++++++++++++++++++++++++++++++
                    // variables ueta, pressure: zero gradient
                    // in addition to the zero-gradient condition in each cell, 
                    // this requires a global condition which enforces that all cells have the same average value at the "r=0"-edge
                    // ++++++++++++++++++++++++++++++++++++++++++++

                    default: {
                        TrafoMtxS[iB] = GetImplicitCondit_Gradient0_Stage1(BasisS[iB]);
                        TrafoMtxSstage2[iB] = GetImplicitCondit_Gradient0_Stage2(BasisS[iB]);

                        break;
                    }
                }
                // }
                //else {Gradient 0_Stage1
                //}

                NpS[iB] = BasisS[iB].Length;
                Debug.Assert(NpS[iB] == TrafoMtxS[iB].NoOfRows);
                Debug.Assert(NpS[iB] >= TrafoMtxS[iB].NoOfCols);
                for(int n = TrafoMtxS[iB].NoOfCols; n < TrafoMtxS[iB].NoOfRows; n++) {
                    if(iB < 2) {
                        UnusedDofs_local.Add(n + NOff);                // Unused DOFs due to local Reduction
                        UnusedDofs_global.Add(n + NOff);               // Unused DOFs due to local Reduction
                    } else {
                        UnusedDofs_local.Add(n + NOff);                // Unused DOFs due to local Reduction
                        UnusedDofs_global.Add(n + NOff);                // Unused DOFs due to local Reduction
                        if(n == TrafoMtxS[iB].NoOfRows - 1) {
                            UnusedDofs_global.Add(NOff + TrafoMtxS[iB].NoOfCols - 1);    // Unused DOFs due to GLOBAL Reduction
                        }
                    }
                }
                NOff += NpS[iB];
            }
            UnusedDofs_global.Sort();


            // create transformation matrix
            // ============================
            var Trafo_local = new BlockMsrMatrix(map);
            var Trafo_global = new BlockMsrMatrix(map);

            Trafo_local.AccEyeSp(1.0);
            Trafo_global.AccEyeSp(1.0);
            long j0Cell = gdat.CellPartitioning.i0;

            //// Hier muss ich die Connection einfügen!
            foreach(int jCell in R0BndyCells.ItemEnum) { // loop over cells at x=0 -- boundary...
                long jCellGlob = j0Cell + jCell;

                long i0 = map.GetBlockI0(jCellGlob);
                int L = map.GetBlockLen(jCellGlob);
                Debug.Assert(L == NpS.Sum());

                Trafo_local.ClearBlock(i0, i0, L, L);

                long i00 = i0;
                for(int iB = 0; iB < GAMMA; iB++) { // loop over dependent variables in the PDE...
                    Debug.Assert(i00 == map.GlobalUniqueCoordinateIndex(iB, jCell, 0));
                    Trafo_local.AccBlock(i00, i00, 1.0, TrafoMtxS[iB]);
                    i00 += NpS[iB];
                }
            }

            long j0;
            if(map.MpiRank == 0) {
                int jCellFix = R0BndyCells.ItemEnum.First();
                j0 = map.GetBlockI0(jCellFix + j0Cell);
            } else {
                j0 = -1;
            }
            j0 = j0.MPIMax();
            //j0 = j0.MPIBroadcast(0);



            foreach(int jCell in R0BndyCells.ItemEnum.Skip((map.MpiRank == 0 ? 1 : 0))) { // loop over cells at x=0 -- boundary...
                long jCellGlob = j0Cell + jCell;
                long i0 = map.GetBlockI0(jCellGlob);

                int L = map.GetBlockLen(jCellGlob);
                Debug.Assert(L == NpS.Sum());

                long i00 = i0;
                long j00 = j0;

                for(int iB = 0; iB < GAMMA; iB++) { // loop over dependent variables in the PDE...
                    switch(iB) {
                        case 0:
                        case 1: {
                            // already fixed to zero, do nothing here
                            break;
                        }
                        default: {
                            Trafo_global.ClearBlock(i00, i00, NpS[iB], NpS[iB]);
                            Trafo_global.AccBlock(i00, i00, 1.0, TrafoMtxSstage2[iB]);
                            Trafo_global.AccBlock(i00, j00, 1.0, MultidimensionalArray.CreateEye(NpS[iB]));

                            // fix with respect to 0 cell
                            break;
                        }
                    }
                    i00 += NpS[iB];
                    j00 += NpS[iB];

                }
            }

            // return
            // ======
            // BIG TRAFO
            BlockMsrMatrix Trafo;
            Trafo = BlockMsrMatrix.Multiply(Trafo_local, Trafo_global);
            return (Trafo, UnusedDofs_local.ToArray(), UnusedDofs_global.ToArray());
        }

        /* Fk, 11apr24: deactivated, because seemingly un-used

        public static (BlockMsrMatrix, int[], int[]) CreateTrafoMatrix2(UnsetteledCoordinateMapping map) {
            if(map.CellDepLength)
                throw new NotSupportedException();
            var gdat = map.GridDat;

            // find cell at the inner, r = x = 0 array
            // =======================================
            var R0BndyCells = GetR0BndyCells(gdat);

            // get transformation matrix for each basis in the mapping
            // =======================================================

            var BasisS = map.BasisS.ToArray();
            int GAMMA = BasisS.Length;

            var TrafoMtxS = new MultidimensionalArray[GAMMA];
            var TrafoMtxSstage2 = new MultidimensionalArray[GAMMA];


            var NpS = new int[GAMMA];
            var UnusedDofs_local = new List<int>();
            var UnusedDofs_global = new List<int>();
            int NOff = 0;
            for(int iB = 0; iB < GAMMA; iB++) {  // loop over dependent variables in the PDE...
                                                 // Distinguish here between the different coefficients u_eta_;u_r, U_xi, p!
                switch(iB) {
                    case 0:
                    case 1: {
                        var getTransformationMatrixValue0 = GetTransformationMatrixValue0(BasisS[iB]);
                        TrafoMtxS[iB] = getTransformationMatrixValue0.GetSolutionSpace();
                        break;
                    }
                    default: {
                        TrafoMtxS[iB] = GetImplicitCondit_Gradient0_Stage1(BasisS[iB]);
                        TrafoMtxSstage2[iB] = GetImplicitCondit_Gradient0_Stage2(BasisS[iB]);

                        break;
                    }
                }
                // }
                //else {Gradient 0_Stage1
                //}

                NpS[iB] = BasisS[iB].Length;
                Debug.Assert(NpS[iB] == TrafoMtxS[iB].NoOfRows);
                Debug.Assert(NpS[iB] >= TrafoMtxS[iB].NoOfCols);
                for(int n = TrafoMtxS[iB].NoOfCols; n < TrafoMtxS[iB].NoOfRows; n++) {
                    if(iB < 2) {
                        UnusedDofs_local.Add(n + NOff);                // Unused DOFs due to local Reduction
                        UnusedDofs_global.Add(n + NOff);               // Unused DOFs due to local Reduction
                    } else {
                        UnusedDofs_local.Add(n + NOff);                // Unused DOFs due to local Reduction
                        UnusedDofs_global.Add(n + NOff);                // Unused DOFs due to local Reduction
                        if(n == TrafoMtxS[iB].NoOfRows - 1) {
                            UnusedDofs_global.Add(NOff + TrafoMtxS[iB].NoOfCols - 1);    // Unused DOFs due to GLOBAL Reduction
                        }
                    }
                }
                NOff += NpS[iB];
            }
            UnusedDofs_global.Sort();
            // create transformation matrix
            // ============================
            var Trafo_local = new BlockMsrMatrix(map);
            var Trafo_global = new BlockMsrMatrix(map);

            Trafo_local.AccEyeSp(1.0);
            Trafo_global.AccEyeSp(1.0);
            long j0Cell = gdat.CellPartitioning.i0;

            //// Hier muss ich die Connection einfügen!
            foreach(int jCell in R0BndyCells.ItemEnum) { // loop over cells at x=0 -- boundary...
                long jCellGlob = j0Cell + jCell;

                long i0 = map.GetBlockI0(jCellGlob);
                int L = map.GetBlockLen(jCellGlob);
                Debug.Assert(L == NpS.Sum());

                Trafo_local.ClearBlock(i0, i0, L, L);

                long i00 = i0;
                for(int iB = 0; iB < GAMMA; iB++) { // loop over dependent variables in the PDE...
                    Debug.Assert(i00 == map.GlobalUniqueCoordinateIndex(iB, jCell, 0));
                    Trafo_local.AccBlock(i00, i00, 1.0, TrafoMtxS[iB]);
                    i00 += NpS[iB];
                }
            }

            foreach(int jCell in R0BndyCells.ItemEnum.Skip(1)) { // loop over cells at x=0 -- boundary...
                long jCellGlob = j0Cell + jCell;
                long i0 = map.GetBlockI0(jCellGlob);

                int jCellFix = R0BndyCells.ItemEnum.First();
                long j0 = map.GetBlockI0(jCellFix + j0Cell);

                int L = map.GetBlockLen(jCellGlob);
                Debug.Assert(L == NpS.Sum());

                long i00 = i0;
                long j00 = j0;

                for(int iB = 0; iB < GAMMA; iB++) { // loop over dependent variables in the PDE...
                    switch(iB) {
                        case 0:
                        case 1: {
                            // already fixed to zero, do nothing here
                            break;
                        }
                        default: {
                            Trafo_global.ClearBlock(i00, i00, NpS[iB], NpS[iB]);
                            Trafo_global.AccBlock(i00, i00, 1.0, TrafoMtxSstage2[iB]);
                            Trafo_global.AccBlock(i00, j00, 1.0, MultidimensionalArray.CreateEye(NpS[iB]));

                            // fix with respect to 0 cell
                            break;
                        }
                    }
                    i00 += NpS[iB];
                    j00 += NpS[iB];

                }
            }

            // return
            // ======
            // BIG TRAFO
            var Trafo = BlockMsrMatrix.Multiply(Trafo_local, Trafo_global);
            return (Trafo, UnusedDofs_local.ToArray(), UnusedDofs_global.ToArray());
        }
        */

        /// <summary>
        /// Multiply from both sides with the Trafo_local Matrix. Change of Basis
        /// </summary>
        public BlockMsrMatrix GetManipulatedMatrix(BlockMsrMatrix M) {
            // Was heißt das?
            if(!M.RowPartitioning.EqualsPartition(m_map))
                throw new ArgumentException("Row partition mismatch");
            if(!M.RowPartitioning.EqualsPartition(m_map))
                throw new ArgumentException("Row partition mismatch");

            // Mutiplication!
            var R = BlockMsrMatrix.Multiply(TrafoTranspose, BlockMsrMatrix.Multiply(M, Trafo));


            var map = m_map;
            var gdat = map.GridDat;
            long j0Cell = gdat.CellPartitioning.i0;
            int GAMMA = map.BasisS.Count;

            foreach(int jCell in R0BndyCells.ItemEnum) {
                long jCellGlob = j0Cell + jCell;
                long i0 = map.GetBlockI0(jCellGlob);

                if(i0 == 0) {
                    foreach(int n in UnusedDofs_local) {
                        R[i0 + n, i0 + n] = 1.0;
                    }
                } else {
                    foreach(int n in UnusedDofs_global) {
                        R[i0 + n, i0 + n] = 1.0;
                    }
                }
            }

            return R;
        }


        /// <summary>
        /// Applies the restriction operation to the RHS
        /// </summary>
        public double[] GetRestrictedRHS<T>(T rhs) where
            T : IList<double> //
        {
            if(rhs.Count != m_map.LocalLength)
                throw new ArgumentException();

            double[] R = new double[rhs.Count];

            TrafoTranspose.SpMV(1.0, rhs, 0.0, R);

            return R;
        }

        /// <summary>
        /// Go back from the restricted solution to the original DG space.
        /// </summary>
        public double[] GetExtendedSolution<T>(T sol) where
            T : IList<double> //
        {
            if(sol.Count != m_map.LocalLength)
                throw new ArgumentException();

            double[] R = new double[sol.Count];

            Trafo.SpMV(1.0, sol, 0.0, R);

            return R;
        }

        /// <summary>
        /// <see cref="CheckSolutionR0Compatibility(SinglePhaseField, SinglePhaseField, SinglePhaseField, SinglePhaseField, bool)"/>
        /// </summary>
        public void CheckSolutionR0Compatibility<T>(UnsetteledCoordinateMapping helicalSol, T vec, bool throwException = true)
            where T : IList<double> //
        {
            if(helicalSol.BasisS.Count != 4)
                throw new ArgumentException("Expecting 4 DG fields: [ur, ux, ueta, pressure], but got " + helicalSol.BasisS.Count);

            var Cvec = new CoordinateVector(helicalSol.BasisS.Select(b => new SinglePhaseField(b)));
            Cvec.SetV(vec);

            CheckSolutionR0Compatibility(Cvec, throwException);
        }


        /// <summary>
        /// <see cref="CheckSolutionR0Compatibility(SinglePhaseField, SinglePhaseField, SinglePhaseField, SinglePhaseField, bool)"/>
        /// </summary>
        public void CheckSolutionR0Compatibility(CoordinateVector helicalSol, bool throwException = true) {
            if(helicalSol.Fields.Count != 4)
                throw new ArgumentException("Expecting 4 DG fields: [ur, ux, ueta, pressure], but got " + helicalSol.Fields.Count);
            CheckSolutionR0Compatibility(helicalSol.Mapping, throwException);
        }

        /// <summary>
        /// <see cref="CheckSolutionR0Compatibility(SinglePhaseField, SinglePhaseField, SinglePhaseField, SinglePhaseField, bool)"/>
        /// </summary>
        public void CheckSolutionR0Compatibility(CoordinateMapping helicalSol, bool throwException = true) {
            if(helicalSol.Fields.Count != 4)
                throw new ArgumentException("Expecting 4 DG fields: [ur, ux, ueta, pressure], but got " + helicalSol.Fields.Count);
            CheckSolutionR0Compatibility((SinglePhaseField)helicalSol.Fields[0], (SinglePhaseField)helicalSol.Fields[1], (SinglePhaseField)helicalSol.Fields[2], (SinglePhaseField)helicalSol.Fields[3], throwException);
        }

        /// <summary>
        /// Asserts that a solution complies with all compatibility conditions at $`r = 0`$
        /// </summary>
        /// <param name="ur">expected to be 0.0 at $`r = 0`$</param>
        /// <param name="uxi">expected to be 0.0 at $`r = 0`$</param>
        /// <param name="ueta">expected to have a gradient of 0.0 at $`r = 0`$ in $`\xi$-direction</param>
        /// <param name="Pressure">expected to have a gradient of 0.0 at $`r = 0`$ in $`\xi$-direction</param>
        /// <param name="throwException"></param>
        public void CheckSolutionR0Compatibility(SinglePhaseField ur, SinglePhaseField uxi, SinglePhaseField ueta, SinglePhaseField Pressure, bool throwException = true) {
            UnsetteledCoordinateMapping helicalSol = new UnsetteledCoordinateMapping(ur.Basis, uxi.Basis, ueta.Basis, Pressure.Basis);
            IGridData gridData = helicalSol.GridDat;
            if(!ReferenceEquals(gridData, this.GridDat))
                throw new ArgumentException("Grid mismatch/fields assigned to different grid.");

            // Initialize a 3x1 Multidimensional Array
            MultidimensionalArray nd = MultidimensionalArray.Create(3, 1);
            nd[0, 0] = -1;  // Set Point at Bottom Left Corner
            nd[1, 0] = 0;   // Set Point at Mid Left 
            nd[2, 0] = 1;   // Set Point at Top Left Corner

            // Create NodeSet with non-shallow copy
            NodeSet nodes = new NodeSet(gridData.iGeomEdges.EdgeRefElements[0], nd, true);

            // Define EdgeMask for cells at boundaries
            EdgeMask em = new EdgeMask(gridData, X => Math.Abs(X[0]) < 1e-8);

            // Determine the number of cells on the local machine
            int N = em.NoOfItemsLocally;

            // Prepare result arrays for calculations
            MultidimensionalArray result_in = MultidimensionalArray.Create(N, nodes.NoOfNodes);
            MultidimensionalArray result_out = MultidimensionalArray.Create(N, nodes.NoOfNodes);

            // List of variable names for processing
            string[] variables = { "ur", "uxi", "ueta", "Pressure" };
            SinglePhaseField[] fields = new SinglePhaseField[] { ur, uxi, ueta, Pressure };

            // Evaluate each variable and assert conditions
            for(int varIndex = 0; varIndex < helicalSol.NoOfVariables; varIndex++) {
                int offset = 0;
                foreach(Chunk c in em) {
                    int edge = c.i0;
                    int Len = c.Len;
                    fields[varIndex].EvaluateEdge(edge, Len, nodes, result_in, result_out, ResultIndexOffset: offset, ResultPreScale:0.0);
                    offset += Len;
                }
                // Assert conditions based on variable index
                if(varIndex < 2) {
                    // Check for zero value in results
                    Console.WriteLine($"R0Fix: {variables[varIndex]}, zero expected; AbsMax = {Math.Max(result_in.Min().MPIMin().Abs(), result_in.Max().MPIMax().Abs()):0.###e-00} (Min|Max Value = {result_in.Min().MPIMin():0.###e-00} | {result_in.Max().MPIMax():0.###e-00})");

                    if(throwException) {
                        Assert.Less(result_in.Max().MPIMax().Abs(), 1E-14, "R0fix_exeption: Max value for {0} should be less than {2} but is {1}", variables[varIndex], result_in.Max(), 1E-14);
                        Assert.Less(result_in.Min().MPIMin().Abs(), 1E-14, "R0fix_exeption: Min value for {0} should be less than {2} but is {1}", variables[varIndex], result_in.Min(), 1E-14);
                    }
                } else if(varIndex >= 2) {
                    // Check for zero gradient in results
                    double l_inf = Math.Max(result_in.Max().Abs(), result_in.Min().Abs()).MPIMax();
                    Console.WriteLine($"R0Fix: {variables[varIndex]}, zero variation expected; Diff = {result_in.Min().MPIMin() - result_in.Max().MPIMax():0.###e-00} (Min/Max Value = {result_in.Min().MPIMin():0.###e-00} | {result_in.Max().MPIMax():0.###e-00})");

                    if((result_in.Max().MPIMax() - result_in.Min().MPIMin()) >= 1E-13) {
                        Console.WriteLine($"ERROR: Gradient0_Check1 failed: {(result_in.Max().MPIMax() - result_in.Min().MPIMin())} is not smaller than {1E-13}");
                        Console.WriteLine("Check if Check2 fails as well:");
                        if((result_in.Max().MPIMax() - result_in.Min().MPIMin()) / (l_inf) >= 1E-10) {
                            Console.WriteLine($"Gradient0_Check2 failed: {(result_in.Max().MPIMax() - result_in.Min().MPIMin()) / (l_inf)} is not smaller than {1E-10} --> Assertion");
                            
                            if(throwException)
                                throw new Exception("Gradient Zero Check Failed!");
                        } else {
                            Console.WriteLine("No fail!");
                        }
                    }
                }
            }
        }

    }
}
