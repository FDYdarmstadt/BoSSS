using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XheatCommon;
using ilPSP.Tracing;
using ilPSP;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using ilPSP.Utils;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;

namespace BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Post-processing specific to <see cref="StaticDropletConvergenceTest"/>
    /// Specifically, evaluate the velocities at the contactline.
    /// </summary>
    [Serializable]
    public class SlipDropletLogging : SlipDropletLogging<XNSFE_Control> { }

    /// <summary>
    /// Post-processing specific to <see cref="StaticDropletConvergenceTest"/>
    /// </summary>
    [Serializable]
    public class SlipDropletLogging<T> : XNSFEinSituPostProcessingModule<T> where T : XNSFE_Control, new() {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "SlipDroplet";

        /// <summary>
        /// hard-coded name for the Rising bubble
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// Header for the rising bubble log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "UxN_interface", "UxN_contactline", "UxN_boundary", "UxT_boundary", "theta");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int iTimestep, double PhysTime) {
            using (new FuncTrace()) {
                double[] UxN_interface = new double[2];
                double[] UxN_contactline = new double[2];
                double[] UxN_boundary = new double[2];
                double[] UxT_boundary = new double[2];
                double[] contactAngles = new double[2];

                var returnValues = this.ComputeBenchmarkQuantities_Contactline();
                UxN_interface = returnValues.Item1;
                UxN_contactline = returnValues.Item2;
                UxN_boundary = returnValues.Item3;
                UxT_boundary = returnValues.Item4;
                contactAngles = returnValues.Item5;

                for (int i = 0; i < contactAngles.Length; i++) {
                    if (this.SolverMain.MPIRank == 0) {
                        string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", iTimestep, PhysTime,
                        UxN_interface.ElementAt(i), UxN_contactline.ElementAt(i), UxN_boundary.ElementAt(i), UxT_boundary.ElementAt(i), contactAngles.ElementAt(i));
                        Log.WriteLine(line);
                    }
                }
            }
        }

        ThermalParameters m_thermParams;
        ThermalParameters ThermParams {
            get {
                if (m_thermParams == null) {
                    m_thermParams = this.Control.ThermalParameters;
                }
                return m_thermParams;
            }
        }

        internal (double[], double[], double[], double[], double[]) ComputeBenchmarkQuantities_Contactline() {

            var Phi = (LevelSet)LsTrk.LevelSets[0];
            var LevelSetGradient = new VectorField<SinglePhaseField>(SolverMainOverride.GridData.SpatialDimension, Phi.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, (SinglePhaseField)LsTrk.LevelSets[0]);
            SinglePhaseField[] Normals = LevelSetGradient.ToArray();

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadSchemeHelper;
            SpeciesId spcId = LsTrk.SpeciesIdS[0];
            SpeciesId spcIdOther = LsTrk.SpeciesIdS[1];
            var clqs = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(spcId, spcIdOther, 0);
            var gridDat = (GridData)this.SolverMainOverride.GridData;
            var QuadDom = clqs.Domain;
            var boundaryEdge = gridDat.GetBoundaryEdgeMask().GetBitMask();
            var boundaryCutEdge = QuadDom.Intersect(new EdgeMask(gridDat, boundaryEdge, MaskType.Geometrical));
            var factory = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, LsTrk.GridDat.Grid.RefElements[0]);
            clqs = new EdgeQuadratureScheme(factory, boundaryCutEdge);

            // we proceed under the assumption to have 2 contactpoints, one at x < 0, the other at x > 0
            double[] UxN_interface = new double[2];
            double[] UxN_contactline = new double[2];
            double[] UxN_boundary = new double[2];
            double[] UxT_boundary = new double[2];
            double[] contactangle = new double[2];

            EdgeQuadrature.GetQuadrature(new int[] { 6 }, LsTrk.GridDat,
                clqs.Compile(LsTrk.GridDat, 0),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                    int D = SolverMain.Grid.SpatialDimension;

                    // contact point
                    NodeSet Enode_l = QR.Nodes;
                    int trf = LsTrk.GridDat.Edges.Edge2CellTrafoIndex[i0, 0];
                    NodeSet Vnode_l = Enode_l.GetVolumeNodeSet(LsTrk.GridDat, trf, false);
                    NodeSet Vnode_g = Vnode_l.CloneAs();
                    int cell = LsTrk.GridDat.Edges.CellIndices[i0, 0];
                    LsTrk.GridDat.TransformLocal2Global(Vnode_l, Vnode_g, cell);

                    EvalResult[0, 0, 0] = Vnode_g[0, 0];                    

                    // contact line velocity
                    MultidimensionalArray U_IN = MultidimensionalArray.Create(new int[] { 2, 1, 1, 3 });
                    MultidimensionalArray U_OUT = MultidimensionalArray.Create(new int[] { 2, 1, 1, 3 });
                    for (int d = 0; d < D; d++) {
                        this.CurrentVel[d].GetSpeciesShadowField(spcId).EvaluateEdge(i0, Length, QR.Nodes, U_IN.ExtractSubArrayShallow(0, -1, -1, d), U_OUT.ExtractSubArrayShallow(0, -1, -1, d), null, null, null, null, 0, 0);
                        this.CurrentVel[d].GetSpeciesShadowField(spcIdOther).EvaluateEdge(i0, Length, QR.Nodes, U_IN.ExtractSubArrayShallow(1, -1, -1, d), U_OUT.ExtractSubArrayShallow(1, -1, -1, d), null, null, null, null, 0, 0);
                    }

                    // contact angle
                    MultidimensionalArray normal_IN = MultidimensionalArray.Create(new int[] { 1, 1, 3 });
                    MultidimensionalArray normal_OUT = MultidimensionalArray.Create(new int[] { 1, 1, 3 });
                    for (int d = 0; d < D; d++) {
                        Normals[d].EvaluateEdge(i0, Length, QR.Nodes, normal_IN.ExtractSubArrayShallow(-1, -1, d), normal_OUT.ExtractSubArrayShallow(-1, -1, d));
                    }

                    // boundary normal
                    MultidimensionalArray normal_bndy = MultidimensionalArray.Create(new int[] { 1, 1, 3 });
                    for (int d = 0; d < D; d++) {
                        normal_bndy[0,0,d] = LsTrk.GridDat.Edges.NormalsForAffine[i0,d];
                    }

                    // contact line normal
                    var nI = new Vector(normal_IN.ExtractSubArrayShallow(0, 0, -1).To1DArray()).Normalize();
                    var nB = new Vector(normal_bndy.ExtractSubArrayShallow(0, 0, -1).To1DArray());
                    var tC = nI.CrossProduct(nB).Normalize();
                    var nC = tC.CrossProduct(nI).Normalize();
                    var tB = nB.CrossProduct(tC).Normalize();

                    // Extract Velocities
                    var uA = new Vector(U_IN.ExtractSubArrayShallow(0, 0, 0, -1).To1DArray());
                    var uB = new Vector(U_IN.ExtractSubArrayShallow(1, 0, 0, -1).To1DArray());

                    double theta = Math.Acos(-(nI * nB)) * (180 / Math.PI);

                    EvalResult[0, 0, 1] = (uA - uB) * nI;
                    EvalResult[0, 0, 2] = (uA - uB) * nC;
                    EvalResult[0, 0, 3] = (uA - uB) * nB;
                    EvalResult[0, 0, 4] = (uA - uB) * tB;
                    EvalResult[0, 0, 5] = theta;
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        if (ResultsOfIntegration[i, 0] > 0) {
                            UxN_interface[1] += ResultsOfIntegration[i, 1];
                            UxN_contactline[1] += ResultsOfIntegration[i, 2];
                            UxN_boundary[1] += ResultsOfIntegration[i, 3];
                            UxT_boundary[1] += ResultsOfIntegration[i, 4];
                            contactangle[1] += ResultsOfIntegration[i, 5];
                        } else {
                            UxN_interface[0] += ResultsOfIntegration[i, 1];
                            UxN_contactline[0] += ResultsOfIntegration[i, 2];
                            UxN_boundary[0] += ResultsOfIntegration[i, 3];
                            UxT_boundary[0] += ResultsOfIntegration[i, 4];
                            contactangle[0] += ResultsOfIntegration[i, 5];
                        }
                    }
                }
            ).Execute();

            return (UxN_interface, UxN_contactline, UxN_boundary, UxT_boundary, contactangle);
        }
    }
}