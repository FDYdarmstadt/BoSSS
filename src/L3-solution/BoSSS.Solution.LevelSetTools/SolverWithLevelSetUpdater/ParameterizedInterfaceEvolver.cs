using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.ParameterizedLevelSet;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// Derivative of <see cref="LevelSet"/>, mainly for the aesthetics 
    /// of having a pair of evolver and level-set.
    /// </summary>
    public class ParameterizedLevelSet : LevelSet {
        internal double xSemiAxis;
        internal double ySemiAxis;
        internal double yCenter;

        public ParameterizedLevelSet(ParameterizedLevelSetControl control, Basis b, string id) : base(b, id) {
            this.xSemiAxis = control.xSemiAxis;
            this.ySemiAxis = control.ySemiAxis;
            this.yCenter = control.yCenter;
        }

        void Ellipsis(MultidimensionalArray X, MultidimensionalArray phi) {
            int NoOfNodes = X.GetLength(0);
            for (int i = 0; i < NoOfNodes; i++) {
                double x = X[i, 0];
                double y = X[i, 1];
                phi[i] = y - yCenter + ySemiAxis * Math.Sqrt(1 - x.Pow2() / xSemiAxis.Pow2());
            }
        }

        public void Project() {
            this.ProjectField(Ellipsis);
        }
    }

    /// <summary>
    /// 
    /// </summary>
    /// <remarks>
    /// 
    /// </remarks>
    public class ParameterizedLevelSetEvolver : ILevelSetEvolver {

        /// <summary>
        /// specialized timestepper  for the evolution of the Parameterized-LS
        /// </summary>
        ParameterizedLevelSetTimeStepper Parameterized_TimeStepper;

        IList<string> parameters;

        int m_HMForder;

        int SpatialDimension;

        ParameterizedLevelSet m_ls;

        string levelSetName;
        public IList<string> ParameterNames => parameters;

        public IList<string> VariableNames => new string[] { };

        // nothing to do
        public Action<DualLevelSet, double, double, bool, IReadOnlyDictionary<string, DGField>, IReadOnlyDictionary<string, DGField>> AfterMovePhaseInterface => null;

        /// <summary>
        /// <see cref="ILevelSetEvolver.InternalFields"/>; here, empty;
        /// </summary>
        public IDictionary<string, DGField> InternalFields {
            get { return null; }
        }


        public ParameterizedLevelSetEvolver(string interfaceName, ParameterizedLevelSet ls, ParameterizedLevelSetControl control, int hMForder, int D) {
            this.levelSetName = interfaceName;
            this.m_ls = ls;
            this.m_HMForder = hMForder;
            this.SpatialDimension = D;
            parameters = BoSSS.Solution.NSECommon.VariableNames.LevelSetGradient(D);
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(interfaceName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));


            if (control == null)
                throw new ArgumentNullException("LevelSetEvolution needs and instance of FourierLevSetControl!");

            //create specialized parameterized timestepper
            Parameterized_TimeStepper = ParameterizedLevelSetFactory.Build_Timestepper(control);
        }


        public void MovePhaseInterface(DualLevelSet levelSet, double time, double dt, bool incremental, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            ParameterizedLevelSet ls = (ParameterizedLevelSet)levelSet.DGLevelSet;
            if (!object.ReferenceEquals(ls, m_ls)) {
                throw new ApplicationException("level-set mismatch");
            }

            double forceX = ComputeForceX(levelSet, ParameterVarFields);

            //Parameterized_TimeStepper.UpdateParameterizedLevelSet();
            var Param1 = Parameterized_TimeStepper.MoveLevelSet(dt, forceX, ls.xSemiAxis, ls.ySemiAxis, ls.yCenter);
            
            
            ls.xSemiAxis = Param1.xSemi1;
            ls.ySemiAxis = Param1.ySemi1;
            ls.yCenter = Param1.yCenter1;
            ((ParameterizedLevelSet)levelSet.DGLevelSet).Project();
        }

        private double ComputeForceX(DualLevelSet levelSet, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            double forceX = 0;
            {
                var LsTrk = levelSet.Tracker;
                int D = LsTrk.GridDat.SpatialDimension;
                SinglePhaseField[] meanVelocity = D.ForLoop(
                    d => (SinglePhaseField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d))]
                    );

                VectorField<SinglePhaseField> LevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(
                    d => new SinglePhaseField(levelSet.DGLevelSet.Basis)
                    ));
                LevSetGradient.Clear();
                LevSetGradient.Gradient(1.0, levelSet.DGLevelSet);

                IEnumerable<string> requiredSpecies = LsTrk.GetSpeciesSeparatedByLevSet(levelSet.LevelSetIndex);
                IEnumerable<SpeciesId> requiredSpeciesId = requiredSpecies.Select(spc => LsTrk.GetSpeciesId(spc));
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(requiredSpeciesId, m_HMForder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(levelSet.LevelSetIndex, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, m_HMForder),
                    delegate (int j0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        var VelocityValues = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                        var LevelSetGradValues = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);

                        for (int d = 0; d < D; d++) {
                            meanVelocity[d].Evaluate(j0, Length, QR.Nodes, VelocityValues.ExtractSubArrayShallow(-1, -1, d));
                            LevSetGradient[d].Evaluate(j0, Length, QR.Nodes, LevelSetGradValues.ExtractSubArrayShallow(-1, -1, d));
                        }

                        //var Normals = LsTrk.DataHistories[levelSet.LevelSetIndex][0].GetLevelSetNormals(QR.Nodes, j0, Length);

                        for (int j = 0; j < Length; j++) { // loop over cells
                            for (int k = 0; k < QR.NoOfNodes; k++) { // loop over nodes

                                double acc = 0;
                                for (int d = 0; d < D; d++) {
                                    acc += VelocityValues[j, k, d] * LevelSetGradValues[j, k, d];//Normals[j, k, d];
                                                                                                 //acc1 += LevelSetGradValues[j, k, d];
                                }
                                EvalResult[j, k, 0] = acc;
                            }
                        }

                        //EvalResult.ExtractSubArrayShallow(-1, -1, 0).Multiply(1.0,
                        //    VelocityValues, Normals, 0.0, "jk", "jkd", "jkd");


                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        //( i0, Length, ResultsOfIntegration) => {
                        for (int i = 0; i < Length; i++)
                            forceX += ResultsOfIntegration[i, 0];
                    }
                 ).Execute();
            }

            return forceX;
        }

            
    }
}
