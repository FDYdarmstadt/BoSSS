using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// checks and verifications
    /// </summary>
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
        static public void CheckMatrixZeroInEmptyCutCells(this LevelSetTracker LsTrk, IMutableMatrixEx OpMatrix, UnsetteledCoordinateMapping RowMapping, IEnumerable<SpeciesId> SpcToCheck, MultiphaseCellAgglomerator UsedAgglom, int CutCellQuadOrder) {



            //
            // Verification: 
            // Note that we can have cut-cells with zero volume. This is because the initial marking of cut-cells must be 'aggressive',
            // because missing a cut cell is a huge problem.
            // Later, when construction the quadrature rules for these cells, their volume might actually be zero.
            //
            // Here, we want to test that such cells have no contribution to the matrix.
            //


            int J = LsTrk.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int NoOfEq = RowMapping.BasisS.Count;
            bool[] isXrow = RowMapping.BasisS.Select(b => b is XDGBasis).ToArray();

            
            List<(int j, SpeciesId sid, int iEq, int n, double CelVol, double errVal)> allErrors = null;

            foreach( var SpcId in SpcToCheck) {
                var spcMask = LsTrk.Regions.GetSpeciesMask(SpcId).GetBitMask();
                var CCvol = LsTrk.GetXDGSpaceMetrics(SpcToCheck.ToArray(), CutCellQuadOrder).CutCellMetrics.CutCellVolumes[SpcId];

                var AgglomSourcesBitMask = UsedAgglom.GetAgglomerator(SpcId).AggInfo.SourceCells.GetBitMask();


                for(int j = 0; j < J; j++) {

                    if(spcMask[j]) {
                        // for performance reasons, check only cells in the species mask.

                        if(CCvol[j] <= 0 || AgglomSourcesBitMask[j]) {
                            // cell 'j' is an empty cut-cell:
                            // check that corresponding matrix rows are zero
                            for(int iEq = 0; iEq < NoOfEq; iEq++) {
                                if(isXrow[iEq]) {
                                    int N = RowMapping.GetNumberOfModes(LsTrk, iEq, j, SpcId);
                                    for(int n = 0; n < N; n++) {
                                        int iRow = RowMapping.LocalUniqueCoordinateIndex(LsTrk, iEq, j, SpcId, n) + RowMapping.i0;
                                        int NoOfNonZeros = OpMatrix.GetNoOfNonZerosPerRow(iRow);
                                        if(NoOfNonZeros > 0) {

                                            int[] colIdx = null;
                                            double[] entries = null;
                                            OpMatrix.GetRow(iRow, ref colIdx, ref entries);
                                            double err = entries.Max();


                                            if(CCvol[j] != 0.0) {
                                                // could be a cut-cell with negative quadrature weight sum (could happen with HMF)
                                                // so be a little forgiving


                                                if(err > Math.Sqrt(Math.Sqrt(BLAS.MachineEps + CCvol[j].Abs()))) {
                                                    if(allErrors == null) {
                                                        allErrors = new List<(int j, SpeciesId sid, int iEq, int n, double CelVol, double errVal)>();
                                                    }

                                                    allErrors.Add((j, SpcId, iEq, n, CCvol[j], err));

                                                }



                                            } else {
                                                if(allErrors == null) {
                                                    allErrors = new List<(int j, SpeciesId sid, int iEq, int n, double CelVol, double errVal)>();
                                                }
                                                allErrors.Add((j, SpcId, iEq, n, CCvol[j], err));
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

            if(allErrors != null && allErrors.Count > 0) {
                Console.WriteLine("REm: remove hardcoded plot.");
                Tecplot.Tecplot.PlotFields(new DGField[] { LsTrk.LevelSets[0] as DGField }, "SeriousError", 0.0, 4);

                OpMatrix.SaveToTextFileSparse("OpMatrix.txt");
                var ErrBitMask = new BitArray(J);
                foreach(var tt in allErrors) {
                    ErrBitMask[tt.j] = true;
                }
                var ErrMask = new CellMask(LsTrk.GridDat, ErrBitMask);
                ErrMask.SaveToTextFile("CellsWithErrors.csv", WriteHeader:false);

                Console.Error.WriteLine("ERROR: Nonzeros in Operator matrix for empty cut cells:");
                Console.Error.WriteLine("    cell\tSpc\tVer\tmod\tClVol\tErr");
                foreach(var tt in allErrors) {
                    Console.Error.WriteLine($"    {tt.j}\t{LsTrk.GetSpeciesName(tt.sid)}\t{tt.iEq}\t{tt.n}\t{tt.CelVol}\t{tt.errVal}");
                }
                                                
                throw new ArithmeticException($"Found non-zero row in matrix corresponding to an empty cut-cell.");
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
        static public void CheckVectorZeroInEmptyCutCells<T>(this LevelSetTracker LsTrk, T Vector, UnsetteledCoordinateMapping RowMapping, IEnumerable<SpeciesId> SpcToCheck, MultiphaseCellAgglomerator UsedAgglom, int CutCellQuadOrder) 
            where T: IList<double> //            
        {


            int NoOfEq = RowMapping.BasisS.Count;
            bool[] isX = RowMapping.BasisS.Select(b => b is XDGBasis).ToArray();

            foreach( var SpcId in SpcToCheck) {
                var spcMask = LsTrk.Regions.GetSpeciesMask(SpcId).GetBitMask();
                var CCvol = LsTrk.GetXDGSpaceMetrics(SpcToCheck.ToArray(), CutCellQuadOrder).CutCellMetrics.CutCellVolumes[SpcId];

                var AgglomSourcesBitMask = UsedAgglom.GetAgglomerator(SpcId).AggInfo.SourceCells.GetBitMask();


                int J = LsTrk.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                for(int j = 0; j < J; j++) {
                    if(spcMask[j]) {
                        // for performance reasons, check only cells in the species mask.

                        if(CCvol[j] <= 0 || AgglomSourcesBitMask[j]) {
                            // cell 'j' is an empty cut-cell:
                            // check that corresponding vector entries are zero
                            for(int iEq = 0; iEq < NoOfEq; iEq++) {
                                if(isX[iEq]) {
                                    int N = RowMapping.GetNumberOfModes(LsTrk, iEq, j, SpcId);
                                    for(int n = 0; n < N; n++) {
                                        int iRow = RowMapping.LocalUniqueCoordinateIndex(LsTrk, iEq, j, SpcId, n);

                                        if(Vector[iRow] != 0) {
                                            if(CCvol[j] != 0.0) {
                                                // could be a cut-cell with negative quadrature weight sum (could happen with HMF)
                                                // so be a little forgiving

                                                if(Vector[iRow].Abs() > Math.Sqrt(Math.Sqrt(BLAS.MachineEps + CCvol[j].Abs()))) {
                                                    throw new ArithmeticException($"Found non-zero row in matrix corresponding to an empty cut-cell.");
                                                }

                                            } else {
                                                throw new ArithmeticException("Found non-zero entry in vector corresponding to an empty cut-cell.");
                                            }
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


        public static void CheckGauss1stOrder(LevelSetTracker LsTrk, int QuadOrder) {
            string[] species = new[] { "A", "B" };
            double[] signS = new double[] { +1, -1 };
            SpeciesId[] speciesIds = species.Select(s => LsTrk.GetSpeciesId(s)).ToArray();
            var GridData = LsTrk.GridDat;
            var helper = LsTrk.GetXDGSpaceMetrics(speciesIds, QuadOrder, 1).XQuadSchemeHelper;

            var surfScheme = helper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            for(int iSpc = 0; iSpc < 2; iSpc++) {
                // ++++++++++++++++++++++++++++++++++++++++++++++
                // check Gauss theorem for 1st order polynomials
                // ++++++++++++++++++++++++++++++++++++++++++++++

                CellQuadratureScheme volScheme = helper.GetVolumeQuadScheme(speciesIds[iSpc]);
                EdgeQuadratureScheme edgeScheme = helper.GetEdgeQuadScheme(speciesIds[iSpc]);
                double sign = signS[iSpc];


                var Vol = CellQuadrature(GridData, volScheme.Compile(GridData, QuadOrder));
                var VolBnd = CellVolumeFromGauss(LsTrk, edgeScheme.Compile(GridData, QuadOrder), surfScheme.Compile(GridData, QuadOrder), sign);

                int D = GridData.SpatialDimension;
                for(int d = 0; d < D; d++) {
                    var Err = Vol.CloneAs();
                    Err.Acc(-1.0, VolBnd.ExtractSubArrayShallow(-1, d));
                    double gaussErrTot = Err.L2Norm();

                    if(gaussErrTot >= 1.0e-5)
                        Console.Error.WriteLine(string.Format("Significant Gauß integral error for species B: {0}", gaussErrTot));
                        //throw new ApplicationException(string.Format("Significant Gauß integral error for species B: {0}", gaussErrTot));
                }

            }

        }


        /// <summary>
        /// Computes cut-cell-volume as \f$ \oint_{K_{j,\frakS}} (x,0) \cdot \vec{n} dS \f$
        /// </summary>
        static private MultidimensionalArray CellVolumeFromGauss(LevelSetTracker LsTrk, ICompositeQuadRule<QuadRule> edgeRule, ICompositeQuadRule<QuadRule> surfRule, double speciesSign) {
            var GridData = LsTrk.GridDat;
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int D = GridData.SpatialDimension;
            var ret = MultidimensionalArray.Create(J, D);
            var E2C = GridData.iLogicalEdges.CellIndices;


            // contribution from edges
            BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                new int[] { D }, GridData,
                edgeRule,
                delegate(int i0, int Length, QuadRule Qr, MultidimensionalArray EvalResult) { // Evaluate
                    var Normals = GridData.iGeomEdges.NormalsCache.GetNormals_Edge(Qr.Nodes, i0, Length);
                    var Nodes = GridData.GlobalNodes.GetValue_EdgeSV(Qr.Nodes, i0, Length);

                    EvalResult.Clear();
                    EvalResult.Multiply(1.0, Nodes, Normals, 0.0, "ikd", "ikd", "ikd");
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        int iEdge = i0 + i;
                        int jCell0 = E2C[iEdge, 0];
                        int jCell1 = E2C[iEdge, 1];

                        for(int d = 0; d < D; d++) {
                            ret[jCell0, d] += ResultsOfIntegration[i, d]; // IN-cell counts positive
                            if(jCell1 >= 0)
                                ret[jCell1, d] -= ResultsOfIntegration[i, d]; // OUT-cell counts negative
                        }
                    }
                }).Execute();

            // contribution from level-set
            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { D }, GridData,
                surfRule,
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(QR.Nodes, i0, Length);
                    var Nodes = GridData.GlobalNodes.GetValue_Cell(QR.Nodes, i0, Length);

                    EvalResult.Clear();
                    EvalResult.Multiply(1.0, Nodes, Normals, 0.0, "ikd", "ikd", "ikd");
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    var A = ret.ExtractSubArrayShallow(new int[] { i0, 0}, new int[] { i0 + Length - 1, D - 1 });
                    A.Acc(speciesSign, ResultsOfIntegration);
                }).Execute();

            // return
            return ret;
        }


        /// <summary>
        /// Computes the volume of the cut cells according to <paramref name="volRule"/>
        /// </summary>
        static private MultidimensionalArray CellQuadrature(GridData GridData, ICompositeQuadRule<QuadRule> volRule) {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { 1 }, GridData,
                volRule,
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    //Integrand.Evaluate(i0, Length, 0, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    EvalResult.SetAll(1.0);
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    var A = ret.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 });
                    var B = ResultsOfIntegration.ExtractSubArrayShallow(-1, 0);
                    A.Set(B);

                }).Execute();

            return ret;
        }

    }
}
