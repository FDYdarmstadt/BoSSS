using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.EnergyCommon;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Statistic.QuadRules;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.TimeStepping;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FreeXNSE {

    public class ParameterFunctionPair {
        public double ParameterValue { get; private set; }
        public ScalarFunction Function { get; private set; }
        public ScalarFunction[] Gradient { get; private set; }
        public ParameterFunctionPair(double parameter, ScalarFunction function, ScalarFunction[] gradient) {
            ParameterValue= parameter;
            Function = function;
            Gradient = gradient;
        }

        public ParameterFunctionPair(double parameter, Func<double[], double> function, Func<double[], double>[] gradient) {
            ParameterValue = parameter;
            Function = NonVectorizedScalarFunction.Vectorize(function);
            Gradient = gradient.Select(g => NonVectorizedScalarFunction.Vectorize(g)).ToArray();
        }

        public void UpdateValue(double newvalue) {
            ParameterValue = newvalue;
        }
    }

    /// <summary>
    /// Derivative of <see cref="LevelSet"/>, mainly for the aesthetics 
    /// of having a pair of evolver and level-set.
    /// </summary>
    public class ParameterizedLevelSet : LevelSet {
        private ParameterFunctionPair[] ParameterFunctionPairs;

        ScalarFunction[] m_TestFunctions;
        private ScalarFunction[] TestFunctions {
            get {
                if(m_TestFunctions == null) {
                    m_TestFunctions = ParameterFunctionPairs.Select(p => p.Function).ToArray();
                    //m_TestFunctions = new ScalarFunction[this.GridDat.SpatialDimension * ParameterFunctionPairs.Length];
                    //for(int i = 0; i < ParameterFunctionPairs.Length; i++) {
                    //    for(int j = 0; j < this.GridDat.SpatialDimension; j++) {
                    //        var Moment = new Moment(j, i);
                    //        m_TestFunctions[i * this.GridDat.SpatialDimension + j] = Moment.Evaluate;
                    //    }
                    //}
                }
                return m_TestFunctions;
            }
        }

        internal class Moment {

            int[] p;
            internal Moment(int[] p) {
                this.p = p;
            }

            internal void Evaluate(MultidimensionalArray input, MultidimensionalArray output) {
                output.Clear();
                output.AccConstant(1.0);
                MultidimensionalArray temp = input.CloneAs();
                for(int i = 0; i < p.Length; i++) {
                    temp.ExtractSubArrayShallow(-1, i).ApplyAll(x => x.Pow(p[i]));
                    output.Multiply(1.0, output, temp.ExtractSubArrayShallow(-1, i), 0.0, "ij", "ij", "ij");
                }

            }
        }


        public ParameterizedLevelSet(ParameterFunctionPair[] ParameterFunctionPairs, Basis b, string id) : base(b, id) {
            this.ParameterFunctionPairs = ParameterFunctionPairs;
            var normalizedParameterValues = this.ParameterFunctionPairs.Select(p => p.ParameterValue).ToArray();
            //normalizedParameterValues.Normalize();
            for(int i = 0; i < normalizedParameterValues.Length; i++) {
                ParameterFunctionPairs[i].UpdateValue(normalizedParameterValues[i]);
            }
            ProjectParameterLevelSet();
        }

        private ICompositeQuadRule<QuadRule> quadruleInterface;
        private LevelSetTracker lsTrk;
        private DGField[] InterfaceVelocity;
        public void MovePhaseInterface(double dt, LevelSetTracker _lsTrk, int order, DGField[] _interfaceVelocity) {
            this.lsTrk = _lsTrk;
            var SchemeHelper = lsTrk.GetXDGSpaceMetrics(lsTrk.SpeciesIdS.ToArray(), order).XQuadSchemeHelper;
            var scheme = SchemeHelper.GetLevelSetquadScheme(0, lsTrk.Regions.GetCutCellMask4LevSet(0));
            this.quadruleInterface = scheme.Compile(lsTrk.GridDat, order);

            InterfaceVelocity = _interfaceVelocity;
            var M = ComputeMassMatrix();

            double[] oldValues = ParameterFunctionPairs.Select(p => p.ParameterValue).ToArray();
            double[] temp = new double[oldValues.Length + 1];
            double[] deltaValues = new double[oldValues.Length];

            var A = ComputeOperatorMatrix();
            M.Acc(dt, A);
            A.MatVecMul(-dt, oldValues, 0.0, temp);

            deltaValues = M.LeastSquareSolve(temp);

            //var M = ComputeMomentOperator(order);

            //double[] oldValues = ParameterFunctionPairs.Select(p => p.ParameterValue).ToArray();

            //for(int i = 0; i < oldValues.Length; i++) {
            //    Console.WriteLine(oldValues[i]);
            //}

            //double[] temp = ComputeMomentVector();
            //temp.ScaleV(dt);

            //deltaValues = M.LeastSquareSolve(temp);
            ////deltaValues = M.Solve(temp);

            double[] newValues = new double[oldValues.Length];
            for(int i = 0; i < deltaValues.Length; i++) {
                newValues[i] = oldValues[i] + deltaValues[i];
            }
            newValues.Normalize();

            for(int i = 0; i < newValues.Length; i++) {
                Console.WriteLine(newValues[i]);
                ParameterFunctionPairs[i].UpdateValue(newValues[i]);
            }

            ProjectParameterLevelSet();
        }

        MultidimensionalArray ComputeMassMatrix() {
            int N = ParameterFunctionPairs.Length;
            MultidimensionalArray M = MultidimensionalArray.Create(N + 1,N);
                 
            CellQuadrature.GetQuadrature(new int[] { N, N }, lsTrk.GridDat,
                quadruleInterface,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray F = MultidimensionalArray.Create(Length, QR.NoOfNodes, N);
                    MultidimensionalArray GlobalNodes = MultidimensionalArray.Create(Length, QR.NoOfNodes, lsTrk.GridDat.SpatialDimension);
                    lsTrk.GridDat.TransformLocal2Global(QR.Nodes, i0, Length, GlobalNodes);
                    for(int i = 0; i < Length; i++) {
                        for (int j = 0; j < N; j++) {
                            ParameterFunctionPairs[j].Function(GlobalNodes.ExtractSubArrayShallow(i, -1, -1), F.ExtractSubArrayShallow(i, -1, j));
                        }
                    }

                    EvalResult.Multiply(1.0, F, F, 0.0, "ijkl", "ijk", "ijl");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int k = 0; k < Length; k++) {
                        for (int i = 0; i < N; i++) {
                            for (int j = 0; j < N; j++) {
                                M[i, j] += ResultsOfIntegration[k, i, j];
                            }
                        }
                    }
                }
            ).Execute();



            M.ExtractSubArrayShallow(N, -1).AccConstant(1.0);

            return M;
        }

        MultidimensionalArray ComputeOperatorMatrix() {
            int N = ParameterFunctionPairs.Length;
            MultidimensionalArray M = MultidimensionalArray.Create(N + 1, N);

            CellQuadrature.GetQuadrature(new int[] { N, N }, lsTrk.GridDat,
                quadruleInterface,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray F = MultidimensionalArray.Create(Length, QR.NoOfNodes, N);
                    MultidimensionalArray GradF = MultidimensionalArray.Create(Length, QR.NoOfNodes, N, lsTrk.GridDat.SpatialDimension);
                    MultidimensionalArray U = MultidimensionalArray.Create(Length, QR.NoOfNodes, lsTrk.GridDat.SpatialDimension);
                    MultidimensionalArray UxGradF = MultidimensionalArray.Create(Length, QR.NoOfNodes, N);

                    MultidimensionalArray GlobalNodes = MultidimensionalArray.Create(Length, QR.NoOfNodes, lsTrk.GridDat.SpatialDimension);
                    lsTrk.GridDat.TransformLocal2Global(QR.Nodes, i0, Length, GlobalNodes);
                    for(int i = 0; i < Length; i++) {
                        for(int j = 0; j < N; j++) {
                            ParameterFunctionPairs[j].Function(GlobalNodes.ExtractSubArrayShallow(i, -1, -1), F.ExtractSubArrayShallow(i, -1, j));
                            for(int k = 0; k < lsTrk.GridDat.SpatialDimension; k++) {
                                ParameterFunctionPairs[j].Gradient[k](GlobalNodes.ExtractSubArrayShallow(i, -1, -1), GradF.ExtractSubArrayShallow(i, -1, j, k));
                            }
                        }
                    }

                    for(int k = 0; k < lsTrk.GridDat.SpatialDimension; k++) {
                        InterfaceVelocity[k].Evaluate(i0, Length, QR.Nodes, U.ExtractSubArrayShallow(-1, -1, k));
                    }

                    UxGradF.Multiply(1.0, U, GradF, 0.0, "ijk", "ijl", "ijkl");
                    EvalResult.Multiply(1.0, UxGradF, F, 0.0, "ijkl", "ijl", "ijk");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int k = 0; k < Length; k++) {
                        for(int i = 0; i < N; i++) {
                            for(int j = 0; j < N; j++) {
                                M[i, j] += ResultsOfIntegration[k, i, j];
                            }
                        }
                    }
                }
            ).Execute();

            return M;
        }

        double[] ComputeMomentVector() {
            int N = TestFunctions.Length;
            int D = lsTrk.GridDat.SpatialDimension;
            double[] r = new double[N];

            CellQuadrature.GetQuadrature(new int[] { N }, lsTrk.GridDat,
                quadruleInterface,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray F = MultidimensionalArray.Create(Length, QR.NoOfNodes, N);
                    MultidimensionalArray U = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                    MultidimensionalArray Normal = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                    MultidimensionalArray Norm = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                    MultidimensionalArray UxN = MultidimensionalArray.Create(Length, QR.NoOfNodes);

                    MultidimensionalArray GlobalNodes = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                    lsTrk.GridDat.TransformLocal2Global(QR.Nodes, i0, Length, GlobalNodes);
                    for(int i = 0; i < Length; i++) {
                        for(int j = 0; j < N; j++) {
                            TestFunctions[j](GlobalNodes.ExtractSubArrayShallow(i, -1, -1), F.ExtractSubArrayShallow(i, -1, j));
                        }
                    }

                    EvaluateGradient(i0, Length, QR.Nodes, Normal);
                    for(int d = 0; d < D; d++) {
                        var GradPhi_d = Normal.ExtractSubArrayShallow(-1, -1, d);
                        Norm.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                    }
                    Norm.ApplyAll(x => 1.0 / Math.Sqrt(x));
                    Normal.Multiply(1.0, Normal, Norm, 0.0, "ikd", "ikd", "ik");

                    for(int k = 0; k < D; k++) {
                        InterfaceVelocity[k].Evaluate(i0, Length, QR.Nodes, U.ExtractSubArrayShallow(-1, -1, k));
                    }

                    UxN.Multiply(1.0, U, Normal, 0.0, "ij", "ijl", "ijl");
                    EvalResult.Multiply(1.0, UxN, F, 0.0, "ijk", "ij", "ijk");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int k = 0; k < Length; k++) {
                        for(int i = 0; i < N; i++) {
                            r[i] += ResultsOfIntegration[k, i];                            
                        }
                    }
                }
            ).Execute();

            //r.ScaleV(-1.0);

            return r;
        }
        MultidimensionalArray ComputeMomentOperator(int order) {
            int N = TestFunctions.Length;
            int M = ParameterFunctionPairs.Length;

            MultidimensionalArray O = MultidimensionalArray.Create(N, M);
            double[] R1, R0;
            double delta = Math.Sqrt(BLAS.MachineEps);

            R0 = ComputeMomentOperator(this.lsTrk, order);

            for(int i = 0; i < M; i++) {
                double oldValue = ParameterFunctionPairs[i].ParameterValue;

                ParameterFunctionPairs[i].UpdateValue(oldValue + delta);

                var LevelSet1 = new ParameterizedLevelSet(this.ParameterFunctionPairs, this.Basis, this.Identification);
                var tracker = new LevelSetTracker(this.lsTrk.GridDat, this.lsTrk.CutCellQuadratureType, this.lsTrk.NearRegionWidth, this.lsTrk.SpeciesNames.ToArray(), LevelSet1);
                tracker.UpdateTracker(0.0);

                R1 = ComputeMomentOperator(tracker, order);

                ParameterFunctionPairs[i].UpdateValue(oldValue);

                R1.AccV(-1.0, R0);
                R1.ScaleV(1.0 / delta);
                O.SetColumn(i, R1);
            }

            return O;
        }


        double[] ComputeMomentOperator(LevelSetTracker levelSetTracker, int order) {
            int N = TestFunctions.Length;
            int D = levelSetTracker.GridDat.SpatialDimension;
            double[] r = new double[N];

            var SchemeHelper = levelSetTracker.GetXDGSpaceMetrics(levelSetTracker.SpeciesIdS.ToArray(), order).XQuadSchemeHelper;
            var Volscheme = SchemeHelper.GetVolumeQuadScheme(levelSetTracker.SpeciesIdS.FirstOrDefault());
            var quadruleVolume = Volscheme.Compile(levelSetTracker.GridDat, order);

            CellQuadrature.GetQuadrature(new int[] { N }, levelSetTracker.GridDat,
                quadruleVolume,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    MultidimensionalArray F = MultidimensionalArray.Create(Length, QR.NoOfNodes, N);

                    MultidimensionalArray GlobalNodes = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                    levelSetTracker.GridDat.TransformLocal2Global(QR.Nodes, i0, Length, GlobalNodes);
                    for(int i = 0; i < Length; i++) {
                        for(int j = 0; j < N; j++) {
                            TestFunctions[j](GlobalNodes.ExtractSubArrayShallow(i, -1, -1), F.ExtractSubArrayShallow(i, -1, j));
                        }
                    }

                    EvalResult.ExtractSubArrayShallow(-1, -1, -1).Acc(1.0, F);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int k = 0; k < Length; k++) {
                        for(int i = 0; i < N; i++) {
                            r[i] += ResultsOfIntegration[k, i];
                        }
                    }
                }
            ).Execute();

            return r;
        }

        

        void ProjectParameterLevelSet() {
            this.Clear();
            for(int i = 0; i < ParameterFunctionPairs.Length; i++) {
                this.ProjectField(ParameterFunctionPairs[i].ParameterValue, ParameterFunctionPairs[i].Function);
            }
        }
    }

    /// <summary>
    /// Parameterized Level Set
    /// </remarks>
    public class ParameterizedLevelSetEvolver : ILevelSetEvolver {


        /// <summary>
        /// ctor
        /// </summary>
        public ParameterizedLevelSetEvolver(string levelSetName, int hMForder, int D, IGridData grd) {
            this.SpatialDimension = D;
            this.m_HMForder = hMForder;
            this.levelSetName = levelSetName;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(this.levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            this.m_grd = grd;
        }

        IGridData m_grd;
        int SpatialDimension;
        double AgglomThreshold;
        int m_HMForder;
        string levelSetName;
        string[] parameters;

        /// <summary>
        /// should only be the interface velocity vector; typically, a phase-averaged velocity.
        /// </summary>
        public IList<string> ParameterNames => parameters;

        /// <summary>
        /// nix
        /// </summary>
        public IList<string> VariableNames => null;

        // nothing to do
        public Func<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>, bool> AfterMovePhaseInterface => null;


        /// <summary>
        /// Provides access to the internally constructed extension velocity.
        /// <see cref="ILevelSetEvolver.InternalFields"/>
        /// </summary>
        public IDictionary<string, DGField> InternalFields { 
            get {                
                return null;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(var tr = new FuncTrace()) {

                int D = levelSet.Tracker.GridDat.SpatialDimension;

                SinglePhaseField[] meanVelocity = D.ForLoop(
                    d => (SinglePhaseField)ParameterVarFields[parameters[d]]
                    );

                ((ParameterizedLevelSet)levelSet.DGLevelSet).MovePhaseInterface(dt, levelSet.Tracker, m_HMForder, meanVelocity);
            }
        }
    }
}
