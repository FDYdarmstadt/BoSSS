/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation.Grid;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using MPI.Wrappers.Utils;
using BoSSS.Platform;
using ilPSP.Tracing;


namespace BoSSS.Foundation {
    /// <summary>
    /// This class provides extended functionalities for the class <see cref="SubgridCoordinateMapping"/> regarding modified mass matrices.
    /// </summary>
    public class SubgridCoordinateMappingExt: SubgridCoordinateMapping {
        
        /// <summary>
        /// Constructs a mapping from associated with a given <see cref="SubGrid"/> on the basis of the original mapping
        /// </summary>
        /// <param name="context">the context object</param>
        /// <param name="cm">The original mapping</param>
        /// <param name="subgrid">The subgrid</param>
        public SubgridCoordinateMappingExt(Context context, CoordinateMapping cm, SubGrid subgrid): base(context,subgrid, cm.Fields) {
        }
            /// <summary>
        /// Constructs a new mapping from an ordered list of fields and assigns a subgrid.
        /// The corresponding indices are created.
        /// </summary>
        /// <param name="context"></param>
        /// <param name="_fields">the list of <see cref="Field"/>'s for this mapping.</param>
        /// <param name="subgrid">The <see cref="SubGrid"/> that the compressed version of this coordinate mapping will be defined on.</param>
    
        public SubgridCoordinateMappingExt(Context context, SubGrid subgrid, params Field[] _fields):base(context, subgrid, (IList<Field>)_fields) {
        }

        /// <summary>
        /// Constructs a new mapping from an ordered list of fields and assigns a subgrid.
        /// The corresponding indices are created.
        /// </summary>
        /// <param name="context"></param>
        /// <param name="_fields">the list of <see cref="Field"/>'s for this mapping.</param>
        /// <param name="subgrid">The <see cref="SubGrid"/> that the compressed version of this coordinate mapping will be defined on.</param>
        public SubgridCoordinateMappingExt(Context context, SubGrid subgrid, IList<Field> _fields):base( context, subgrid,_fields) {        
        }
         /// <summary>
        /// Method for extracting a submatrix associated with the narrow band as well as the partial vector of
        /// the affine part of the matrix. 
        /// </summary>
        public void SetupSubmatrix(double[] affine, MsrMatrix original, out double[] compressedAffine, out MsrMatrix subgridMatrix, MsrMatrix massMatrix, MsrMatrix subgridMassMatrix, bool[] timeDepComp, double dt) {
            using (new FuncTrace()) {
                int entries = m_subgridIndices.Length;
                int increment = Fields.Count;

                compressedAffine = new double[entries];

                for (int pos = 0; pos < m_subgridIndices.Length; pos++) {

                    compressedAffine[pos] = affine[subgridIndices[pos]];

                }
                //------------------ List of all required quantities for setting up the new matrix with correct dimensions-----------------------------------------------------------------------
                int globalNoCells = subgrid.GlobalNoOfCells;
                int localNoCells = subgrid.LocalNoOfCells;

                int totalDoF = 0;
                for (int k = 0; k < this.Fields.Count; k++) {
                    totalDoF += this.BasisS[k].MaximalLength;
                }
                int DoFperCell = 0;
                for (int k = 0; k < this.Fields.Count; k++) {
                    DoFperCell += this.BasisS[k].MaximalLength;
                }
         // The mass matrix for the time relevant components has just the respective 
        // entries. this is required for the expression u_0/dt * MassMatrix
                int totalDoFTimeDep = 0;
                for (int k = 0; k < this.Fields.Count; k++) {
                    if (timeDepComp[k]) {
                        totalDoFTimeDep += this.BasisS[k].MaximalLength;
                    }
                }
                int DoFperCellTimeDep = 0;
                for (int k = 0; k < this.Fields.Count; k++) {
                    if (timeDepComp[k]) {
                        DoFperCellTimeDep += this.BasisS[k].MaximalLength;
                    }
                }


        //All quantities for the actual operator matrix
                int globalDoF = globalNoCells * DoFperCell;
                int localDoF = localNoCells * DoFperCell;
                Partition submatrixPartition = new Partition(localDoF, 1,  MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                Partition subgridPartition = new Partition(m_subgrid.LocalNoOfCells, 1, MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                subgridMatrix = new MsrMatrix(submatrixPartition, submatrixPartition);
        
        // All quantities for the shorter mass matrix
                //int globalDoFTimeDep = globalNoCells * DoFperCellTimeDep;
                //int localDoFTimeDep = localNoCells * DoFperCellTimeDep;
                //Partition submatrixPartitionTimeDep = new Partition(localDoFTimeDep, 1, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
              
                subgridMassMatrix = new MsrMatrix(submatrixPartition, submatrixPartition);

                long[,] neighbours;
                long[] SubgridIndex2GlobalIndex;

         //In this method, both matrices are instantiated simultaneously to avoid additional loops!      
                neighbours = m_subgrid._GridData.GetGlobalNeighbourIndices();

                int[] cellsGrid = m_subgrid.SubgridIndex2LocalCellIndex;
                int[] cellsSubgrid = m_subgrid.LocalCellIndex2SubgridIndex;
                if (cellsGrid.Length == 0) {
                    throw new ArgumentNullException("No cells are contained in the subgrid."
                + "Either use a wider range for computations or return to standard quantities on the full domain");
                }
                //Number of locally updated cells contained in the subgrid
                subgridNUpdate = 0;

                for (int s0 = 0; s0 < cellsGrid.Length; s0++) {
                    if (cellsGrid[s0] < m_Context.GridDat.NoOfLocalUpdatedCells) {
                        subgridNUpdate++;
                    }
                }

                //Global cell indices in the subgrid are required for gathering all column entries that should be in the submatrix
                //Construction of this global index array by collecting the respective IDs from all processors......
                // we may call this LocalSubgrid2GlobalSubgridIndices!
                long[] cellsGlobalSubgrid = new long[subgridNUpdate];

                for (int s1 = 0; s1 < cellsGrid.Length; s1++) {

                    if (cellsGrid[s1] < m_Context.GridDat.NoOfLocalUpdatedCells) {
                        cellsGlobalSubgrid[s1] = cellsGrid[s1] + m_Context.GridDat.GridPartition.i0;

                    }
                }

                //Global indices for the original grid that are relevant.....
                //negative entries for all cells not contained in the subgrid are
                //maintained
                long[] GlobalIDAllProcessors;
                long[] Local2Global = new long[m_Context.GridDat.NoOfLocalUpdatedCells];

                for (int s2 = 0; s2 < cellsSubgrid.Length; s2++) {

                    //   if ((cellsSubgrid[s2] < m_Context.GridDat.NoOfLocalUpdatedCells) && (cellsSubgrid[s2] >= 0)) {
                    if ((cellsSubgrid[s2] < subgridNUpdate) && (cellsSubgrid[s2] >= 0)) {
                        Local2Global[s2] = cellsSubgrid[s2] + subgridPartition.i0;
                    }
                    //else if ((cellsSubgrid[s2] < m_Context.GridDat.NoOfLocalUpdatedCells) && (cellsSubgrid[s2] < 0) )
                    else if ((cellsSubgrid[s2] < subgridNUpdate) && (cellsSubgrid[s2] < 0))
                        Local2Global[s2] = int.MinValue;

                }
                //------------- Attempt for parallelization----------------------------

                //First step: Create an array of displacement to be used in the MPI_ALLGETHERV method
                Partition p = new Partition(subgridNUpdate, 1, MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                int size = p.Size;
                int[] displacement = new int[size];
                int[] localEntries = new int[size];
                for (int i = 0; i < size; i++) {
                    displacement[i] = (int)p.GetI0Offest(i);
                    localEntries[i] = p.GetLocalLength(i);
                }

                SubgridIndex2GlobalIndex = new long[p.TotalLength];

                //Second step: Parallel communication of the global indices based on the MPI_ALLGETHERV method
                //Note that the indices are circulated among all of the processes (not only the root)  and that the length of these
                //index arrays may be variable,

                unsafe {
                    fixed (
                       void* pSend = &cellsGlobalSubgrid[0],
                       pRec = &SubgridIndex2GlobalIndex[0],
                       pdispl = &displacement[0],
                       pEntries = &localEntries[0]) {
                        csMPI.Raw.Allgatherv((IntPtr)pSend, cellsGlobalSubgrid.Length, csMPI.Raw._DATATYPE.LONG_LONG_INT,
                            (IntPtr)pRec, (IntPtr)pEntries, (IntPtr)pdispl,
                            csMPI.Raw._DATATYPE.LONG_LONG_INT, csMPI.Raw._COMM.WORLD);

                    }
                }


                Partition p1 = new Partition(Local2Global.Length, 1, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                int size1 = p1.Size;
                int[] displacement1 = new int[size1];
                int[] localEntries1 = new int[size1];
                for (int i = 0; i < size1; i++) {
                    displacement1[i] = (int)p1.GetI0Offest(i);
                    localEntries1[i] = p1.GetLocalLength(i);
                }

                GlobalIDAllProcessors = new long[p1.TotalLength];
                unsafe {
                    fixed (
                       void* p1Send = &Local2Global[0],
                       p1Rec = &GlobalIDAllProcessors[0],
                       p1displ = &displacement1[0],
                       p1Entries = &localEntries1[0]) {
                        csMPI.Raw.Allgatherv((IntPtr)p1Send, Local2Global.Length, csMPI.Raw._DATATYPE.LONG_LONG_INT,
                            (IntPtr)p1Rec, (IntPtr)p1Entries, (IntPtr)p1displ,
                            csMPI.Raw._DATATYPE.LONG_LONG_INT, csMPI.Raw._COMM.WORLD);

                    }
                }

                long[,] Globalneighbours;
                int[] displacement2 = new int[size1];
                int[] localEntries2 = new int[size1];
                for (int i = 0; i < size1; i++) {
                    displacement2[i] = (int)p1.GetI0Offest(i) * neighbours.GetLength(1);
                    localEntries2[i] = p1.GetLocalLength(i) * neighbours.GetLength(1);
                }
                Globalneighbours = new long[p1.TotalLength, neighbours.GetLength(1)];
                unsafe {
                    fixed (
                       void* p2Send = &neighbours[0, 0],
                       p2Rec = &Globalneighbours[0, 0],
                       p2displ = &displacement2[0],
                       p2Entries = &localEntries2[0]) {
                        csMPI.Raw.Allgatherv((IntPtr)p2Send, m_Context.GridDat.NoOfLocalUpdatedCells * neighbours.GetLength(1), csMPI.Raw._DATATYPE.LONG_LONG_INT,
                         (IntPtr)p2Rec, (IntPtr)p2Entries, (IntPtr)p2displ,
                         csMPI.Raw._DATATYPE.LONG_LONG_INT, csMPI.Raw._COMM.WORLD);

                    }
                }
                //--------------------------------------------------------------------------------------            

                //ATTENTION: For more than one field this process has to be looped over the total number of fields!

                int location_Subgrid_Row = 0;
                int location_Grid_Row = 0;
                int location_Subgrid_Column = 0;
                int location_Grid_Column = 0;

               

                // ---------------------------------for each cell that is contained in in the Subgrid -------------------------------v
                for (int i = 0; i < cellsGrid.Length; i++) {
                    //variables for the cell blocks and matrix positions
                    //global indices!
                    int index_Grid = cellsGrid[i];
                    int index_Subgrid = cellsSubgrid[index_Grid];
                    //vorher waren diese Indizes hier lokal definiert, verschoben in Z354,Z.355

                    for (int k = 0; k < this.Fields.Count; k++) {
                        //int location_Subgrid_Row_TimeDep = 0;
                        //int location_Grid_Row_TimeDep = 0;
                        //int location_Subgrid_Column_TimeDep = 0;
                        //int location_Grid_Column_TimeDep = 0;

                        if (index_Grid >= 0) {
                            location_Subgrid_Row = DoFperCell * index_Subgrid + (int)subgridMatrix.RowPartition.i0;
                            location_Grid_Row = DoFperCell * index_Grid + (int)original.RowPartition.i0;

                            //if (timeDepComp[k]) {

                            //    location_Subgrid_Row_TimeDep = DoFperCellTimeDep * index_Subgrid + (int)subgridMassMatrix.RowPartition.i0;
                            //    location_Grid_Row_TimeDep = DoFperCellTimeDep * index_Grid + (int)massMatrix.RowPartition.i0;


                            //}

                            //für jede entsprechende Zeile copy the corresponding column entries of this cell and all neighbours
                            for (int n = 0; n < DoFperCell; n++) {
                                location_Subgrid_Column = DoFperCell * index_Subgrid + (int)subgridMatrix.RowPartition.i0;
                                location_Grid_Column = DoFperCell * index_Grid + (int)original.RowPartition.i0;
                                
                                    //location_Subgrid_Column_TimeDep = DoFperCellTimeDep * index_Subgrid + (int)subgridMassMatrix.RowPartition.i0;
                                    //location_Grid_Column_TimeDep = DoFperCellTimeDep * index_Grid + (int)massMatrix.RowPartition.i0;

                                    for (int l = 0; l < this.m_Fields.Length;l++ )
                                        for (int m = 0; m < this.BasisS[l].MaximalLength; m++) {
                                            //Extraction of the cell's own block 
                                            double fillIn;
                                            double fillInTimeDep=0.0;
                                            if (timeDepComp[l]) {
                                                fillInTimeDep = massMatrix[location_Grid_Row, location_Grid_Column];
                                                subgridMassMatrix[location_Subgrid_Row, location_Subgrid_Column] = dt * fillInTimeDep;
                                            }
                                            fillIn = original[location_Grid_Row, location_Grid_Column];
                                            subgridMatrix[location_Subgrid_Row, location_Subgrid_Column] = (fillIn + dt*fillInTimeDep);                                          
                                            location_Subgrid_Column++;
                                            location_Grid_Column++;
                                        }

                                for (int n0 = 0; n0 < neighbours.GetLength(1); n0++) {

                                    if (Globalneighbours[index_Grid + m_Context.GridDat.GridPartition.i0, n0] >= 0) {

                                        long global = GlobalIDAllProcessors[Globalneighbours[index_Grid + m_Context.GridDat.GridPartition.i0, n0]];

                                        //int global = cellsSubgrid[neighbours[index_Grid, n0]];
                                        //Check, if the neighbour cell is part of the subgrid...
                                        if (global >= 0) {

                                            int location_Subgrid_ColumnNeighbours =
                                            (int)global * DoFperCell;
                                            int location_Grid_ColumnNeighbours =
                                            (int)Globalneighbours[index_Grid + m_Context.GridDat.GridPartition.i0, n0] * DoFperCell;
                                            for (int r = 0; r < DoFperCell; r++) {
                                                subgridMatrix[location_Subgrid_Row, location_Subgrid_ColumnNeighbours] = original[location_Grid_Row, location_Grid_ColumnNeighbours];

                                                location_Grid_ColumnNeighbours++;
                                                location_Subgrid_ColumnNeighbours++;
                                            }
                                        }
                                    }
                                }
                                location_Grid_Row++;
                                location_Subgrid_Row++;
                            }
      
                        }
                    }
                }
            }
        }

       
    }
}
