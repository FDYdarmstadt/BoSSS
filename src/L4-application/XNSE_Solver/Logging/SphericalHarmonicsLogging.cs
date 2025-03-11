using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Logging {

    /// <summary>
    /// <see cref="SphericalHarmonicsLogging{T}"/>
    /// </summary>
    [Serializable]
    public class SphericalHarmonicsLogging : SphericalHarmonicsLogging<XNSE_Control> {

    }

    /// <summary>
    /// For single droplet/bubble computations:
    /// decomposition of the fluid interface into Laplace Spherical Harmonics.
    /// This is achieved by minimizing
    /// ```math 
    /// \oint_{\mathfrak{I}} \left(
    ///   | \vec{x} | - \sum_{l,m} Y_{l,m}(\varphi,\theta) \underbrace{y_{l,m}}_{\text{?}}
    ///   \right)^2 dS
    ///   \rightarrow \text{min}
    /// ```
    /// in terms of the unknowns $` y_{l,m} `$ (which are logged into a file).
    /// These are determined by solving the linear system 
    /// ```math 
    /// \forall k,n: \ \
    /// \sum_{l,m} \left( \oint_{\mathfrak{I}} Y_{l,m} Y_{k,n} \dS \right) y_{l,m}
    /// = 
    ///   \oint_{\mathfrak{I}}
    ///     | \vec{x} | Y_{k,n} dS .
    /// ```
    /// </summary>
    [Serializable]
    public class SphericalHarmonicsLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {
        protected override string LogFileName => "SphericalHarmonics";

        /// <summary>
        /// assuming rotational symmetry, excluding the computation of modes m != 0
        /// </summary>
        [DataMember]
        public bool RotSymmetric = false;

        /// <summary>
        /// maximum mode l 
        /// </summary>
        [DataMember]
        public int MaxL = 1;


        /// <summary>
        /// Spherical Harmonics mapping from a pair to a linear index: (l,m) -> idx
        /// For the rotational symmetric case: (l,0) -> idx
        /// </summary>
        /// <param name="l">0, 1, 2, ...</param>
        /// <param name="m">-l, ... , +l</param>
        /// <returns>
        /// 0, 1, 2, 3, ...
        /// </returns>
        static internal int SH_mapping(int l, int m, bool rotSym = false) {
            if(l < 0)
                throw new ArgumentOutOfRangeException();
            if(m < -l)
                throw new ArgumentOutOfRangeException($"lo: l = {l}, m = {m}");
            if(m > l)
                throw new ArgumentOutOfRangeException($"hi: l = {l}, m = {m}");
            if(m != 0 && rotSym)
                throw new ArgumentException("m unequal 0 for rotational symmetry");

            if(rotSym)
                return l;

            int cnt = 0;
            for(int _l = 0; _l <= l; _l++) {
                if(_l < l) {
                    cnt += _l + _l + 1;
                } else {
                    cnt += m + l;
                }
                              

            }
            return cnt;
        }

        /// <summary>
        /// Vector space dimension for spherical harmonics up to (including) a given <paramref name="l"/>
        /// </summary>
        static internal int SH_dim(int l, bool rotSym = false) {
            if(l < 0)
                throw new ArgumentOutOfRangeException();

            if(rotSym)
                return l + 1;

            int cnt = 0;
            for(int _l = 0; _l <= l; _l++) {
                cnt += _l + _l + 1;
            }
            return cnt;
        }

        /// <summary>
        /// inverse mapping of <see cref="SH_mapping"/>, i.e. from a linear index to an (l,m)-pair
        /// </summary>
        /// <param name="idx"></param>
        /// <returns></returns>
        static internal (int l, int m) SH_mappingInv(int idx, bool rotSym = false) {
            if(idx < 0)
                throw new ArgumentOutOfRangeException();

            if(rotSym)
                return (idx, 0);
            
            int i0 = 0;
            for(int _l = 0; _l <= int.MaxValue; _l++) {
                int iE = i0 + _l + _l + 1;
                if(idx >= i0 && idx < iE) {
                    return (_l, idx - i0 - _l);
                }
                i0 = iE;
            }
            return (int.MinValue, int.MinValue);
        }

        /// <summary>
        /// Testcode
        /// </summary>
        void TestSHmappings() {

            int cnt = 0;
            for(int l = 0; l <= MaxL; l++) {
                for(int m = -l; m <= l; m++) {
                    if(m != 0 && RotSymmetric)
                        continue;
                    int i = SH_mapping(l, m, RotSymmetric);
                    if(i != cnt)
                        throw new ApplicationException($"cnt = {cnt}, i = {i}; ({l},{m})");
                    (int _l, int _m) = SH_mappingInv(i, RotSymmetric);
                    if(_m != m || _l != l)
                        throw new ApplicationException($"{i}: ({l},{m})!=({_l},{_m})");
                    cnt++;
                }
                Console.WriteLine();
            }
            if(cnt != SH_dim(MaxL, RotSymmetric))
                throw new ApplicationException($"cnt = {cnt}, dim = {SH_dim(MaxL - 1, RotSymmetric)}");
        }



        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using(new FuncTrace()) {
                TestSHmappings();

                int qOrder = SolverMainOverride.QuadOrder();

                var schemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, qOrder).XQuadSchemeHelper;
                var scheme = schemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask4LevSet(0).Intersect(LsTrk.Regions.GetSpeciesMask("A")));

                int dim = SH_dim(MaxL, RotSymmetric);

                MultidimensionalArray MassMatrix = MultidimensionalArray.Create(dim, dim);
                double[] RHS = new double[dim];

                double[] SH_buf = new double[dim];

                CellQuadrature.GetQuadrature(new int[] { dim, dim + 1 }, LsTrk.GridDat,
                    scheme.Compile(LsTrk.GridDat, qOrder), //  agg.HMForder),
                    delegate (int j0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        var NodesGlobal = LsTrk.GridDat.GlobalNodes.GetValue_Cell(QR.Nodes, j0, Length);

                        for(int jj = 0; jj < Length; jj++) { // loop over cells
                        for(int k = 0; k < QR.NoOfNodes; k++) {
                                Vector X = NodesGlobal.GetRowPt(jj, k);
                                double r = X.Abs();
                                (double theta, double phi) = SphericalHarmonics.GetAngular(X);

                                for(int i = 0; i < dim; i++) {
                                    (int l, int m) = SH_mappingInv(i, RotSymmetric);
                                    SH_buf[i] = SphericalHarmonics.MyRealSpherical(l, m, theta, phi);
                                }

                                for(int i = 0; i < dim; i++) {
                                    for(int j = 0; j < dim; j++) {
                                        EvalResult[jj, k, i, j] = SH_buf[i] * SH_buf[j];
                                    }
                                    EvalResult[jj, k, i, dim] = r * SH_buf[i];
                                }
                            }
                        }
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++) {
                            MassMatrix.Acc(1.0, ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, 0, 0 }, new int[] { i - 1, dim - 1, dim - 1 }));
                            double[] rhs_part = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, 0, dim }, new int[] { i - 1, dim - 1, dim - 2 }).To1DArray();
                            RHS.AccV(1.0, rhs_part);
                        }

                    }
                    ).Execute();

                // summation over all MPI processes
                double[] mpiRHS = RHS.MPISum();
                MultidimensionalArray mpiMassMatrix = MultidimensionalArray.CreateWrapper(MassMatrix.Storage.MPISum(), dim, dim);



                if(SolverMainOverride.MPIRank == 0) {

                    double[] Ylm = mpiMassMatrix.LeastSquareSolve(mpiRHS);

                    string line = PhysTime + Ylm.ToConcatString("\t", "\t", "");
                    Log.WriteLine(line);
                    Log.Flush();

                    LoggedValues.Add(Ylm);
                    while(LoggedValues.Count > 100)
                        LoggedValues.RemoveAt(0);
                }
            }
        }

        internal List<double[]> LoggedValues = new List<double[]>();





        protected override void WriteHeader(TextWriter textWriter) {
            string header = SH_dim(MaxL, RotSymmetric).ForLoop(l => SH_mappingInv(l, RotSymmetric)).ToConcatString("time\t", "\t", "");
            Log.WriteLine(header);
            Log.Flush();
        }

        public override void Setup(IApplication solverMain) {
            base.Setup(solverMain);
        }


    }


    
}
