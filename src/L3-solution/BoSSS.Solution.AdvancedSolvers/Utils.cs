using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {
    static class Utils {

        /*
        public static double rhoDinvA(MsrMatrix A, out MsrMatrix DinvA)
        {
            double rho;

            // extract diagonal matrix
            double[] Diag = A.GetDiagVector();
            int n = Diag.Length;

            // % form (D^-1) A
            //Diag = spdiags(1./Diag, 0, n, n);
            //DinvA = Diag*A;
            DinvA = new MsrMatrix(A.RowPartitioning, A.ColPartition);
            int i0 = A.RowPartitioning.i0;

            int Lr;
            int[] ColumnIdx = null;
            double[] Values = null;

            for (int i = 0; i < n; i++)
            {
                Lr = A.GetRow(i + i0, ref ColumnIdx, ref Values);
                for (int k = 0; k < Lr; k++)
                {
                    Values[k] *= 1.0 / Diag[i];
                }
                DinvA.SetRow(i + i0, ColumnIdx, Values, Lr);
            }
#if DEBUG
            for(int i = 0; i < n; i++) {
                Debug.Assert(Math.Abs(1.0 - DinvA[i + i0, i + i0]) < BLAS.MachineEps * 10);
            }
#endif

            // estimate the largest eigen value from Arnoldi iteration
            int kk = 20;
            var rand = new Random(0);
            double[] v0 = n.ForLoop(i => rand.NextDouble());
            double[][] V; double[,] H; int kact;
            arnoldi(out V, out H, out kact, DinvA, v0, kk, false);
            kk = Math.Min(H.GetLength(0), H.GetLength(1));
            H = H.GetSubMatrix(0, kk, 0, kk);

            rho = MaxAbsEigen(H);

            return rho;
        }

        /*
        /// <summary>
        /// Maximum of the absolute value of all Eigenvalues of <paramref name="H"/>.
        /// </summary>
        static public double MaxAbsEigen(double[,] H) {
            var linalg = new ManagedLinearAlgebraProvider();

            double Eigen;
            int N = H.GetLength(0);
            double[] Matrix = H.Resize(false);
            double[] EigenVect = new double[Matrix.Length];
            double[] diagM = new double[Matrix.Length];
            System.Numerics.Complex[] EigenValues = new System.Numerics.Complex[N];
            linalg.EigenDecomp(false, N, Matrix, EigenVect, EigenValues, diagM);
            Eigen = EigenValues.Select(ev => Complex.Abs(ev)).Max();
            return Eigen;
        }
        */

        /// <summary>
        /// Arnoldi iteration 
        /// </summary>
        /// <param name="V">Output: Arnoldi vectors</param>
        /// <param name="H">Output: </param>
        /// <param name="kact">Output:</param>
        /// <param name="A">Input: (n-by-n) the matrix </param>
        /// <param name="v0">Input: n-vector</param>
        /// <param name="k">Input: number of Arnoldi steps requested</param>
        /// <param name="reorth">Input: (optional) set to 1 for reorthogonalization, (default), set to any other value to switch it off</param>
        /// <remarks>
        /// (c) Ren-Cang Li, rcli@uta.edu,  06/16/07
        /// </remarks>
        public static void arnoldi(out double[][] V, out double[,] H, out int kact, MsrMatrix A, double[] v0, int k, bool reorth = false)
        {
            //%
            //%             -----  kact=k -------
            //%      V      n-by-(k+1)  Arnoldi vectors
            //%      H      (k+1)-by-k
            //%             -----  kact=j<k -------
            //%      V      n-by-j  Arnoldi vectors
            //%      H      j-by-j

            double eps = BLAS.MachineEps;


            int n = A.RowPartitioning.LocalLength;
            if (A.ColPartition.LocalLength != A.RowPartitioning.LocalLength)
            {
                throw new ArgumentException("the sizes of input matrix incorrect");
            }

            V = (k + 1).ForLoop(i => new double[n]);
            H = new double[k + 1, k];

            double nrm2 = v0.L2NormPow2().MPISum().Sqrt();
            if (nrm2 == 0.0)
            {
                throw new ArgumentException("arnoldi: input v0 is a zero vector");
            }

            double tol = n * eps;

            V[0].SetV(v0, 1 / nrm2);   //v(:,1)=v0/nrm2;
            for (int j = 0; j < k; j++)
            {
                double[] vh = new double[n];
                A.SpMVpara(1.0, V[j], 0.0, vh);    //vh = A*V(:,j);   
                double nrmvh = vh.L2NormPow2().MPISum().Sqrt();

                //%   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                //%   by MGS
                for (int i = 0; i < j; i++)
                {
                    double hij = GenericBlas.InnerProd(V[i], vh).MPISum();
                    vh.AccV(-hij, V[i]); //vh = vh - hij*V(:,i);
                    H[i, j] = hij;
                }
                if (reorth)
                {
                    for (int i = 0; i < j; i++)
                    {
                        double tmp = GenericBlas.InnerProd(V[i], vh).MPISum();
                        vh.AccV(-tmp, V[i]);  //vh = vh - tmp*V(:,i);
                        H[i, j] = H[i, j] + tmp;
                    }
                }
                //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                H[j + 1, j] = vh.L2NormPow2().MPISum().Sqrt();
                V[j + 1].SetV(vh, 1.0 / H[j + 1, j]);

                if (H[j + 1, j] <= tol * nrmvh)
                {
                    //%             -----  kact<k -------
                    //%      V      n    -by- kact            Arnoldi vectors
                    //%      H      kact -by- kact
                    // Console.WriteLine("termination at step: " + j);
                    kact = j + 1;
                    V = V.GetSubVector(0, kact);
                    H = H.GetSubMatrix(0, kact, 0, kact);
                    return;
                }
            }
            kact = k;
            Debug.Assert(V.Length == kact + 1);
            Debug.Assert(V.Length == kact + 1);

            //%             -----  kact=k -------
            //%      V       n        -by-  (kact+1)  Arnoldi vectors
            //%      H      (kact+1)  -by-   kact
        }


        /// <summary>
        /// returns the persson sensor field for a given field
        /// </summary>
        public static void getPerssonSensorXDGField(out XDGField PerssonSensorField, XDGField FieldToTest, CellMask mask = null) {
            if (FieldToTest.Basis.Degree == 0) {
                throw new NotSupportedException("Sensor not supported for DG degree 0");
            }

            LevelSetTracker LsTrk = ((XDGBasis)FieldToTest.Basis).Tracker;
            XDGField SensorField = new XDGField(new XDGBasis(LsTrk, 0), "PerssonSensorField-" + FieldToTest.Identification);

            //setting all mean values to minus infinity
            if (mask == null) 
                mask = CellMask.GetFullMask(LsTrk.GridDat);
            foreach (int jCell in mask.ItemEnum) {
                foreach (SpeciesId speciesId in LsTrk.SpeciesIdS) {
                    SensorField.GetSpeciesShadowField(speciesId).SetMeanValue(jCell, -10.0);
                }
            }

            //Choose Quad Order
            int deg = FieldToTest.Basis.Degree;
            var AvailOrders = LsTrk.GetCachedOrders().Where(order => order >= 2 * deg);
            int order2Pick = AvailOrders.Any() ? AvailOrders.Min() : 2 * deg;

            //get the p-1 Projection Up-1
            XDGField fieldPMinus1 = GetPMinus1Projection(LsTrk, FieldToTest, order2Pick);
            var fieldPMinus1LaidBack = new XDGField(FieldToTest.Basis, "fieldPMinus1LaidBack");
            fieldPMinus1LaidBack.AccLaidBack(1.0, fieldPMinus1);

            //do the substraction U - Up-1
            var diffField = FieldToTest.CloneAs();
            diffField.Identification = "diffField";
            diffField.Acc(-1.0, fieldPMinus1LaidBack);
            int N = FieldToTest.Basis.NonX_Basis.Length; // DOFs per cell per species.

            //helper functions
            double[] GetCoords(int j, SpeciesId spc, XDGField field) {
                int iSpc = LsTrk.Regions.GetSpeciesIndex(spc, j);
                if (iSpc < 0)
                    return null;
                else
                    return field.Coordinates.GetRowPart(j, iSpc * N, N);
            }

            void SetValuePersonField(int cell, SpeciesId speciesId, double nominator, double denominator) {
                if (nominator == 0) {
                    SensorField.GetSpeciesShadowField(speciesId).SetMeanValue(cell, -10.0);
                } else {

                    if (denominator == 0) {
                        throw new ArgumentException("denominator zero but nominator not zero");
                    } else {
                        double PerssonSensor = (nominator / denominator).Sqrt();
                        double val = Math.Log10(PerssonSensor);
                        if (val == 0) {
                            Console.WriteLine("WTF");
                        }
                        SensorField.GetSpeciesShadowField(speciesId).SetMeanValue(cell, val);
                    }
                }
            }

            //compute the persson field
            foreach (SpeciesId speciesId in LsTrk.SpeciesIdS) {
                var speciesMask = LsTrk.Regions.GetSpeciesMask(speciesId);
                var cutCellMask = LsTrk.Regions.GetCutCellMask().Intersect(speciesMask);
                var UnCutSpeciesMask = speciesMask.Except(cutCellMask);

                var MMF = LsTrk.GetXDGSpaceMetrics(speciesId, order2Pick).MassMatrixFactory;
                var MMblox = MMF.GetMassMatrixBlocks(FieldToTest.Basis.NonX_Basis, speciesId);

                // 1st, do all the cut cells with non-id mass matrix
                // =================================================
                double[] tmp = new double[N];
                double[] tmp2 = new double[N];
                for (int iSub = 0; iSub < MMblox.jSub2jCell.Length; iSub++) {
                    int jCell = MMblox.jSub2jCell[iSub];
                    double[] CoordsFTT = GetCoords(jCell, speciesId, FieldToTest);
                    double[] CoordsDiff = GetCoords(jCell, speciesId, diffField);
                    if (CoordsFTT == null)
                        continue; // species not present in cell; no contribution.

                    var MM_j = MMblox.MassMatrixBlocks.ExtractSubArrayShallow(new int[] { iSub, 0, 0 }, new int[] { iSub - 1, N - 1, N - 1 });

                    MM_j.GEMV(1.0, CoordsFTT, 0.0, tmp);
                    double denominator = CoordsFTT.InnerProd(tmp);

                    MM_j.GEMV(1.0, CoordsDiff, 0.0, tmp2);
                    double nominator = CoordsDiff.InnerProd(tmp2);

                    SetValuePersonField(jCell, speciesId, nominator, denominator);

                }

                //2nd we do the NonCutCells
                foreach (int iCell in UnCutSpeciesMask.ItemEnum) {
                    double[] CoordsFTT = GetCoords(iCell, speciesId, FieldToTest);
                    double[] CoordsDiff = GetCoords(iCell, speciesId, diffField);

                    // in this iCell, we have an orthonormal basis, i.e. the mass matrix is the identity
                    double denominator = CoordsFTT.L2NormPow2();
                    double nominator = CoordsDiff.L2NormPow2();

                    SetValuePersonField(iCell, speciesId, nominator, denominator);
                }
            }


            PerssonSensorField = SensorField.CloneAs();
        }


        /// <summary>
        /// Computes the Projection onto the (p-1)-Polynomial space (relative to the input field)
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <param name="xdgfieldToTest"></param>
        /// <param name="order2Pick"></param>
        /// <returns>projection</returns>
        public static XDGField GetPMinus1Projection(LevelSetTracker LsTrk, XDGField xdgfieldToTest, int order2Pick) {
            var fieldPMinus1 = new XDGField(new XDGBasis(LsTrk, xdgfieldToTest.Basis.Degree - 1), "fieldPMinus1");
            {
                //get the MassMatrix
                MassMatrixFactory massMatrixFactory = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, order2Pick).MassMatrixFactory;
                BlockMsrMatrix massMatrix = massMatrixFactory.GetMassMatrix(xdgfieldToTest.Mapping, inverse: false);
                BlockMsrMatrix MMPmin1Pmin1Inv = massMatrixFactory.GetMassMatrix(fieldPMinus1.Mapping, inverse: true);

                //get the subMassMatrices of deg p-1
                MsrMatrix MMPmin1P;//Sub of Mass Mat with (rows p-1 ,column p)
                {
                    //compute the lists of indices used by GetSubMatrix()
                    var rowIndices = new List<long>();
                    var colIndices = new List<long>();
                    {
                        int MaxModeR = xdgfieldToTest.Mapping.MaxTotalNoOfCoordinatesPerCell / LsTrk.TotalNoOfSpecies;
                        int MaxModer = fieldPMinus1.Mapping.MaxTotalNoOfCoordinatesPerCell / LsTrk.TotalNoOfSpecies;

                        int iField;
                        int jCell;
                        int nMode;
                        for (int iRow = 0; iRow < fieldPMinus1.Mapping.TotalLength; iRow++) {
                            fieldPMinus1.Mapping.LocalFieldCoordinateIndex(iRow, out iField, out jCell, out nMode);
                            double rMode = (double)nMode;
                            int fac = (int)Math.Floor(rMode / MaxModer);
                            int row = xdgfieldToTest.Mapping.LocalUniqueCoordinateIndex(iField, jCell, nMode + fac * (MaxModeR - MaxModer));
                            rowIndices.Add(row);
                        }
                        for (int i = 0; i < massMatrix.NoOfCols; i++) {
                            colIndices.Add(i);
                        }
                    }
                    MMPmin1P = massMatrix.ToMsrMatrix().GetSubMatrix(rowIndices, colIndices);
                }

                // Here we compute M_{p-1,p-1}^{-1} * M_{p-1,p} * Coords(xdgfieldtoTest)
                var tmpVec = new double[fieldPMinus1.CoordinateVector.Length];
                MMPmin1P.SpMV(1.0, xdgfieldToTest.CoordinateVector, 0.0, tmpVec);
                MMPmin1Pmin1Inv.SpMV(1.0, tmpVec, 0.0, fieldPMinus1.CoordinateVector);


            }
            return fieldPMinus1;
        }

    }
}
