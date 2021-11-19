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
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Logging {

    /// <summary>
    /// <see cref="XNSEinSituPostProcessingModule{T}"/>
    /// </summary>
    [Serializable]
    public class SphericalHarmonicsLogging : SphericalHarmonicsLogging<XNSE_Control> { }

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


        int MaxL = 3;


        /// <summary>
        /// Spherical Harmonics mapping: from a pair to a linear index (l,m) -> idx
        /// </summary>
        /// <param name="l">0, 1, 2, ...</param>
        /// <param name="m">-l, ... , +l</param>
        /// <returns>
        /// 0, 1, 2, 3, ...
        /// </returns>
        static internal int SH_mapping(int l, int m) {
            if(l < 0)
                throw new ArgumentOutOfRangeException();
            if(m < -l)
                throw new ArgumentOutOfRangeException($"lo: l = {l}, m = {m}");
            if(m > l)
                throw new ArgumentOutOfRangeException($"hi: l = {l}, m = {m}");

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

        static internal int SH_dim(int l) {
            if(l < 0)
                throw new ArgumentOutOfRangeException();
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
        static internal (int l, int m) SH_mappingInv(int idx) {
            if(idx < 0)
                throw new ArgumentOutOfRangeException();
            
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
                    int i = SH_mapping(l, m);
                    if(i != cnt)
                        throw new ApplicationException($"cnt = {cnt}, i = {i}; ({l},{m})");
                    (int _l, int _m) = SH_mappingInv(i);
                    if(_m != m || _l != l)
                        throw new ApplicationException($"{i}: ({l},{m})!=({_l},{_m})");
                    cnt++;
                }
                Console.WriteLine();
            }
            if(cnt != SH_dim(MaxL))
                throw new ApplicationException($"cnt = {cnt}, dim = {SH_dim(MaxL - 1)}");
        }



        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using(new FuncTrace()) {
                TestSHmappings();

                int qOrder = SolverMainOverride.QuadOrder();

                var schemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, qOrder).XQuadSchemeHelper;
                var scheme = schemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask4LevSet(0).Intersect(LsTrk.Regions.GetSpeciesMask("A")));

                int dim = SH_dim(MaxL);

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
                                (double theta, double phi) = GetAngular(X);

                                for(int i = 0; i < dim; i++) {
                                    (int l, int m) = SH_mappingInv(i);
                                    SH_buf[i] = MyRealSpherical(l, m, theta, phi);
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
            string header = SH_dim(MaxL).ForLoop(l => SH_mappingInv(l)).ToConcatString("time\t", "\t", "");
            Log.WriteLine(header);
            Log.Flush();
        }

        public override void Setup(IApplication solverMain) {
            base.Setup(solverMain);
        }

        public static double factorial(int m) {
            if(m < 0)
                throw new ArgumentException();
            double r = 1.0;
            for(int i = 1; i <= m; i++)
                r *= i;
            return r;
        }

        public static double Pow(double a, int m) {
            if(m < 0)
                throw new ArgumentException();
            double r = 1.0;
            for(int i = 0; i < m; i++)
                r *= a;
            return r;
        }


        public static (double theta, double phi) GetAngular(Vector X) {
            Vector X0 = X.Normalize();

            double theta = Math.Acos(X0.z);
            double phi = Math.Atan2(X0.y, X0.x);

            if(phi < 0)
                phi += 2 * Math.PI;

            if(phi <= 0 || phi >= Math.PI * 2 || phi.IsNaNorInf())
                throw new ApplicationException($"phi computation failed: {phi} for {X0} -- {X}");
            if(theta < 0 || theta > Math.PI || theta.IsNaNorInf())
                throw new ApplicationException($"theta computation failed: {theta} for {X0} -- {X}");

            Vector X0_recovered = new Vector(
                Math.Cos(phi) * Math.Sin(theta),
                Math.Sin(phi) * Math.Sin(theta),
                Math.Cos(theta));

            double errDist = X0.Dist(X0_recovered);
            if(errDist >= 1.0e-6)
                throw new ApplicationException($"angular coordinate computation failed: {X} -> {X0} -> ({phi},{theta}) -> {X0_recovered}");

            return (theta, phi);
        }


        static SphericalHarmonicsLogging() {
            int l_max = 10;
            Legendre_m = new double[l_max];
            for(int m = 0; m < l_max; m++) {
                Legendre_m[m] = Pow(-1, m) * factorial(2 * m) / (Pow(2, m) * factorial(m));
            }

            Spherical_Nlm = (double[,]) Array.CreateInstance(typeof(double), new int[] { l_max, 2 * l_max + 2 }, new int[] { 0, -l_max });
            for(int l = 0; l < l_max; l++) {
                for(int m = -l; m <= l; m++) {
                    Spherical_Nlm[l,m] = (1 / Math.Sqrt(2 * Math.PI)) * Math.Sqrt(((2.0 * l + 1.0) / 2.0) * factorial(l - m) / factorial(l + m)); 
                }
            }
        }

        static double[] Legendre_m;
       

        public static double MyLegendre(int l, int m, double x) {
            double f, a;

            if(l < m) {
                //f = 0;
                return 0.0;
            } else if(m < 0) {
                a = MyLegendre(l, -m,  x);
                f = a / ((Pow(-1, m)) * factorial(l - m) / factorial(l + m));
                return f;
            } else if(l == m) {
                // compute (1 - x^2)^(m/2)
                if(m % 2 == 0)
                    a = Pow(1 - x * x, m / 2);
                else
                    a = Pow(Math.Sqrt(1 - x * x), m);

                f = Legendre_m[m] * a;  //Pow(-1, m) * factorial(2 * m) / (Pow(2,m) * factorial(m)) * a;
                return f;
            } else {
                f = x * (2 * l - 1) * MyLegendre(l - 1, m, x) - (l + m - 1) * MyLegendre(l - 2, m, x);
                f = f / (l - m);
                return f;
            }
        }

        static double[,] Spherical_Nlm;

        public static double MyRealSpherical(int l, int m, double theta, double phi) {
            if(theta < 0 || theta > Math.PI) {
                throw new ArgumentException();
            }
            if(phi < 0 || phi > 2*Math.PI) {
                throw new ArgumentException();
            }

            double Nlm = Spherical_Nlm[l, Math.Abs(m)]; // Math.Sqrt((2 * l + 1) / 2 * factorial(l - m) / factorial(l + m));

            double Ylm;
            if(m >= 0) {
                Ylm =  Nlm * MyLegendre(l, m, Math.Cos(theta)) * Math.Cos(m * phi);
            } else {
                Ylm =  Nlm * MyLegendre(l, -m, Math.Cos(theta)) * Math.Sin(-m * phi);
            }
            return Ylm;
        }


    }
}
