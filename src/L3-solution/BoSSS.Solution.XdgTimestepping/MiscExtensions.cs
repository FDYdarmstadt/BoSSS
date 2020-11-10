using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {
    public static class MiscExtensions {

        
        /// <summary>
        /// Verifies that matrix rows which should be zero (since the corresponding cut cells are empty) are actually empty.
        /// </summary>
        /// <remarks>
        /// Note that we can have cut-cells with zero volume. This is because the initial marking of cut-cells must be 'aggressive',
        /// because missing a cut cell is a huge problem.
        /// Later, when construction the quadrature rules for these cells, their volume might actually be zero.
        /// These cut-cells correspond to a certain set of rows in <paramref name="OpMatrix"/>, which must be zero.
        /// </remarks>
        static public void CheckMatrixZeroInEmptyCutCells(this LevelSetTracker LsTrk, IMutableMatrixEx OpMatrix, UnsetteledCoordinateMapping RowMapping, IEnumerable<SpeciesId> SpcToCheck, int CutCellQuadOrder) {



            //
            // Verification: 
            // Note that we can have cut-cells with zero volume. This is because the initial marking of cut-cells must be 'aggressive',
            // because missing a cut cell is a huge problem.
            // Later, when construction the quadrature rules for these cells, their volume might actually be zero.
            //
            // Here, we want to test that such cells have no contribution to the matrix.
            //


            int NoOfEq = RowMapping.BasisS.Count;

            foreach( var SpcId in SpcToCheck) {
                var spcMask = LsTrk.Regions.GetSpeciesMask(SpcId).GetBitMask();
                var CCvol = LsTrk.GetXDGSpaceMetrics(SpcToCheck.ToArray(), CutCellQuadOrder).CutCellMetrics.CutCellVolumes[SpcId];


                int J = LsTrk.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                for(int j = 0; j < J; j++) {

                    if(spcMask[j]) {
                        // for performance reasons, check only cells in the species mask.

                        if(CCvol[j] <= 0) {
                            // cell 'j' is an empty cut-cell:
                            // check that corresponding matrix rows are zero
                            for(int iEq = 0; iEq < NoOfEq; iEq++) {
                                int N = RowMapping.GetNumberOfModes(LsTrk, iEq, j, SpcId);
                                for(int n = 0; n < N; n++) {
                                    int iRow = RowMapping.LocalUniqueCoordinateIndex(LsTrk, iEq, j, SpcId, n) + RowMapping.i0;
                                    int NoOfNonZeros = OpMatrix.GetNoOfNonZerosPerRow(iRow);
                                    if(NoOfNonZeros > 0) {
                                        if(CCvol[j] != 0.0) {
                                            // could be a cut-cell with negative quadrature weight sum (could happen with HMF)
                                            // so be a little forgiving

                                            int[] colIdx = null;
                                            double[] entries = null;
                                            OpMatrix.GetRow(iRow, ref colIdx, ref entries);
                                            double err = entries.Max();
                                            
                                            if(err > Math.Sqrt(Math.Sqrt(BLAS.MachineEps + CCvol[j].Abs()))) { 
                                                throw new ArithmeticException($"Found non-zero row in matrix corresponding to an empty cut-cell.");
                                            }


                
                                        } else {
                                            throw new ArithmeticException("Found non-zero row in matrix corresponding to an empty cut-cell.");
                                        }
                                    } else {
                                        // everything is ok.
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Verifies that matrix rows which should be zero (since the corresponding cut cells are empty) are actually empty.
        /// </summary>
        /// <remarks>
        /// Note that we can have cut-cells with zero volume. This is because the initial marking of cut-cells must be 'aggressive',
        /// because missing a cut cell is a huge problem.
        /// Later, when construction the quadrature rules for these cells, their volume might actually be zero.
        /// These cut-cells correspond to a certain set of rows in <paramref name="OpMatrix"/>, which must be zero.
        /// </remarks>
        static public void CheckVectorZeroInEmptyCutCells<T>(this LevelSetTracker LsTrk, T Vector, UnsetteledCoordinateMapping RowMapping, IEnumerable<SpeciesId> SpcToCheck, int CutCellQuadOrder) 
            where T: IList<double> //            
        {


            int NoOfEq = RowMapping.BasisS.Count;

            foreach( var SpcId in SpcToCheck) {
                var spcMask = LsTrk.Regions.GetSpeciesMask(SpcId).GetBitMask();
                var CCvol = LsTrk.GetXDGSpaceMetrics(SpcToCheck.ToArray(), CutCellQuadOrder).CutCellMetrics.CutCellVolumes[SpcId];


                int J = LsTrk.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                for(int j = 0; j < J; j++) {
                    if(spcMask[j]) {
                        // for performance reasons, check only cells in the species mask.

                        if(CCvol[j] <= 0) {
                            // cell 'j' is an empty cut-cell:
                            // check that corresponding vector entries are zero
                            for(int iEq = 0; iEq < NoOfEq; iEq++) {
                                int N = RowMapping.GetNumberOfModes(LsTrk, iEq, j, SpcId);
                                for(int n = 0; n < N; n++) {
                                    int iRow = RowMapping.LocalUniqueCoordinateIndex(LsTrk, iEq, j, SpcId, n) + RowMapping.i0;

                                    if(Vector[iRow] != 0) {
                                        throw new ArithmeticException("Found non-zero entry in vector corresponding to an empty cut-cell.");
                                    } else {

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
}
