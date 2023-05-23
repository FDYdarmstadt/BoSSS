using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Tracing;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Used by <see cref="HeatedWall"/>, <see cref="TwoPhaseCouetteFlow"/>, etc.
    /// </summary>
    [Serializable]
    public class MovingContactLineLogging : MovingContactLineLogging<XNSE_Control> { }

    /// <summary>
    /// Used by <see cref="HeatedWall"/>, <see cref="TwoPhaseCouetteFlow"/>, etc.
    /// </summary>
    [Serializable]
    public class MovingContactLineLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "ContactAngle";

        /// <summary>
        /// Filename for logging contact angles
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle");
            Log.WriteLine(header);
            Log.Flush();
        }

        
        [JsonIgnore]
        [NonSerialized]
        List<double[]> contactPointsRef;

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            using(new FuncTrace()) {

                // contact angles at contact points
                //=================================

                List<double[]> contactPoints = new List<double[]>();
                List<double[]> contactVelocities = new List<double[]>();
                List<double> contactAngles = new List<double>();

                var returnValues = this.ComputeContactLineQuantities();
                contactPoints = returnValues.Item1;
                contactVelocities = returnValues.Item2;
                contactAngles = returnValues.Item3;


                List<double[]> contactPointsSorted = new List<double[]>();
                List<double[]> contactVelocitiesSorted = new List<double[]>();
                List<double> contactAnglesSorted = new List<double>();

                if(contactPointsRef == null) {
                    contactPointsRef = contactPoints;
                    contactPointsSorted = contactPoints;
                    contactVelocitiesSorted = contactVelocities;
                    contactAnglesSorted = contactAngles;
                } else {
                    // sort
                    double eps = 1e-7;

                    do {
                        contactPointsSorted.Clear();
                        contactVelocitiesSorted.Clear();
                        contactAnglesSorted.Clear();

                        eps *= 10;
                        //Console.WriteLine("sorting contact line points");
                        foreach(var ptR in contactPointsRef) {
                            //Console.WriteLine("ref point: ({0}, {1})", ptR[0], ptR[1]);
                            for(int i = 0; i < contactAngles.Count(); i++) {
                                double[] pt = contactPoints.ElementAt(i);
                                //Console.WriteLine("sorting point: ({0}, {1})", pt[0], pt[1]);
                                double xDiff = Math.Abs(pt[0] - ptR[0]);
                                //Console.WriteLine("x diff = {0}", xDiff);
                                double yDiff = Math.Abs(pt[1] - ptR[1]);
                                //Console.WriteLine("y diff = {0}", yDiff);
                                if(xDiff < eps && yDiff < eps) {
                                    //Console.WriteLine("sorted");
                                    contactPointsSorted.Add(pt.CloneAs());
                                    contactVelocitiesSorted.Add(contactVelocities.ElementAt(i).CloneAs());
                                    contactAnglesSorted.Add(contactAngles.ElementAt(i));
                                    break;
                                }
                            }
                        }

                    } while(contactPointsSorted.Count != contactPointsRef.Count && eps < 1.0);

                    contactPointsRef = contactPointsSorted;
                }

                for(int p = 0; p < contactAnglesSorted.Count; p++) {
                    string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime,
                        contactPointsSorted.ElementAt(p)[0], contactPointsSorted.ElementAt(p)[1],
                        contactVelocitiesSorted.ElementAt(p)[0], contactVelocitiesSorted.ElementAt(p)[1], contactAnglesSorted.ElementAt(p));
                    Log.WriteLine(line);
                }
                Log.Flush();

                return;
            }
        }


        ConventionalDGField[] GetMeanVelocityFromXDGField(DGField[] EvoVelocity) {

            IList<string> velocityName = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(Solution.NSECommon.VariableNames.LevelSetCG, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(3));
            IReadOnlyDictionary<string, DGField> parameters = this.SolverMainOverride.LsUpdater.Parameters;

            List<ConventionalDGField> velocity = new List<ConventionalDGField>(3);
            for (int i = 0; i < 3; ++i) {

                if (parameters.TryGetValue(velocityName[i], out DGField velocityField)) {
                    velocity.Add((ConventionalDGField)velocityField);
                }
            }
            return velocity.ToArray();

        }

        
        public Tuple<List<double[]>, List<double[]>, List<double>> ComputeContactLineQuantities() {

            List<double[]> contactPoints = new List<double[]>();
            List<double[]> contactVelocities = new List<double[]>();
            List<double> contactAngles = new List<double>();

            ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(this.CurrentVel);

            var Phi = (LevelSet)LsTrk.LevelSets[0];
            var LevelSetGradient = new VectorField<SinglePhaseField>(SolverMainOverride.GridData.SpatialDimension, Phi.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, (SinglePhaseField)LsTrk.LevelSets[0]);
            SinglePhaseField[] Normals = LevelSetGradient.ToArray();

            XQuadSchemeHelper SchemeHelper = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadSchemeHelper;
            EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"), 0);

            var gridDat = (GridData)this.SolverMainOverride.GridData;
            var QuadDom = SurfaceElement_Edge.Domain;
            var boundaryEdge = gridDat.GetBoundaryEdgeMask().GetBitMask();
            var boundaryCutEdge = QuadDom.Intersect(new EdgeMask(gridDat, boundaryEdge, MaskType.Geometrical));

            var factory = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, LsTrk.GridDat.Grid.RefElements[0]);
            SurfaceElement_Edge = new EdgeQuadratureScheme(factory, boundaryCutEdge);

            EdgeQuadrature.GetQuadrature(new int[] { 5 }, LsTrk.GridDat,
                SurfaceElement_Edge.Compile(LsTrk.GridDat, 0),
                delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                    // contact point
                    NodeSet Enode_l = QR.Nodes;
                    int trf = LsTrk.GridDat.Edges.Edge2CellTrafoIndex[i0, 0];
                    NodeSet Vnode_l = Enode_l.GetVolumeNodeSet(LsTrk.GridDat, trf, false);
                    NodeSet Vnode_g = Vnode_l.CloneAs();
                    int cell = LsTrk.GridDat.Edges.CellIndices[i0, 0];
                    LsTrk.GridDat.TransformLocal2Global(Vnode_l, Vnode_g, cell);
                    //Console.WriteLine("contact point: ({0},{1})", Vnode_g[0, 0], Vnode_g[0, 1]);

                    int D = base.SolverMainOverride.Grid.SpatialDimension;
                    for (int d = 0; d < D; d++) {
                        EvalResult[0, 0, d] = Vnode_g[0, d];
                    }

                    // contact line velocity
                    MultidimensionalArray U_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    MultidimensionalArray U_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        (meanVelocity[d] as SinglePhaseField).EvaluateEdge(i0, length, QR.Nodes, U_IN.ExtractSubArrayShallow(-1, -1, d), U_OUT.ExtractSubArrayShallow(-1, -1, d));
                    }

                    for (int d = 0; d < D; d++) {
                        EvalResult[0, 0, 2 + d] = U_IN[0, 0, d];
                    }

                    // contact angle
                    MultidimensionalArray normal_IN = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    MultidimensionalArray normal_OUT = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        Normals[d].EvaluateEdge(i0, length, QR.Nodes, normal_IN.ExtractSubArrayShallow(-1, -1, d), normal_OUT.ExtractSubArrayShallow(-1, -1, d));
                    }

                    double theta_surf = Math.Atan2(normal_IN[0, 0, 1], normal_IN[0, 0, 0]);
                    double theta_edge = Math.Atan2(LsTrk.GridDat.Edges.NormalsForAffine[i0, 1], LsTrk.GridDat.Edges.NormalsForAffine[i0, 0]);
                    double theta = (theta_surf - theta_edge) * (180 / Math.PI);

                    EvalResult[0, 0, 2 * D] = (theta > 180) ? theta - 180 : theta;
                    //Console.WriteLine("contact angle = {0}", EvalResult[0, 0, 2 * D]);

                },
                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                    int D = SolverMainOverride.Grid.SpatialDimension;
                    for (int i = 0; i < length; i++) {
                        if (ResultsOfIntegration[i, 2 * D] != 0.0) {
                            contactAngles.Add(Math.Abs(ResultsOfIntegration[i, 2 * D]));
                            double[] cp = new double[D];
                            double[] cpV = new double[D];
                            for (int d = 0; d < D; d++) {
                                cp[d] = ResultsOfIntegration[i, d];
                                cpV[d] = ResultsOfIntegration[i, 2 + d];
                            }
                            contactPoints.Add(cp);
                            contactVelocities.Add(cpV);
                        }
                    }
                }
            ).Execute();

            return new Tuple<List<double[]>, List<double[]>, List<double>>(contactPoints, contactVelocities, contactAngles);

        }


    }

    /// <summary>
    /// Used by <see cref="HeatedWall"/>, <see cref="TwoPhaseCouetteFlow"/>, etc.
    /// </summary>
    [Serializable]
    public class MovingContactLineZwoLsLogging : MovingContactLineZwoLsLogging<XNSE_Control> { }

    /// <summary>
    /// Used by <see cref="HeatedWall"/>, <see cref="TwoPhaseCouetteFlow"/>, etc.
    /// Gets the ContactPoints between levelset 0 and 1
    /// </summary>
    [Serializable]
    public class MovingContactLineZwoLsLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "ContactAngle";

        /// <summary>
        /// Filename for logging contact angles
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// 
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-angle");
            Log.WriteLine(header);
            Log.Flush();
        }


        [JsonIgnore]
        [NonSerialized]
        List<double[]> contactPointsRef;

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            using (new FuncTrace()) {

                // contact angles at contact points
                //=================================

                List<double[]> contactPoints = new List<double[]>();
                List<double[]> contactVelocities = new List<double[]>();
                List<double> contactAngles = new List<double>();

                var returnValues = this.ComputeContactLineQuantities();
                contactPoints = returnValues.Item1;
                contactVelocities = returnValues.Item2;
                contactAngles = returnValues.Item3;


                List<double[]> contactPointsSorted = new List<double[]>();
                List<double[]> contactVelocitiesSorted = new List<double[]>();
                List<double> contactAnglesSorted = new List<double>();

                if (contactPointsRef == null) {
                    contactPointsRef = contactPoints;
                    contactPointsSorted = contactPoints;
                    contactVelocitiesSorted = contactVelocities;
                    contactAnglesSorted = contactAngles;
                } else {
                    // sort
                    double eps = 1e-7;

                    do {
                        contactPointsSorted.Clear();
                        contactVelocitiesSorted.Clear();
                        contactAnglesSorted.Clear();

                        eps *= 10;
                        //Console.WriteLine("sorting contact line points");
                        foreach (var ptR in contactPointsRef) {
                            //Console.WriteLine("ref point: ({0}, {1})", ptR[0], ptR[1]);
                            for (int i = 0; i < contactAngles.Count(); i++) {
                                double[] pt = contactPoints.ElementAt(i);
                                //Console.WriteLine("sorting point: ({0}, {1})", pt[0], pt[1]);
                                double xDiff = Math.Abs(pt[0] - ptR[0]);
                                //Console.WriteLine("x diff = {0}", xDiff);
                                double yDiff = Math.Abs(pt[1] - ptR[1]);
                                //Console.WriteLine("y diff = {0}", yDiff);
                                if (xDiff < eps && yDiff < eps) {
                                    //Console.WriteLine("sorted");
                                    contactPointsSorted.Add(pt.CloneAs());
                                    contactVelocitiesSorted.Add(contactVelocities.ElementAt(i).CloneAs());
                                    contactAnglesSorted.Add(contactAngles.ElementAt(i));
                                    break;
                                }
                            }
                        }

                    } while (contactPointsSorted.Count != contactPointsRef.Count && eps < 1.0);

                    contactPointsRef = contactPointsSorted;
                }

                for (int p = 0; p < contactAnglesSorted.Count; p++) {
                    string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", TimestepNo, phystime,
                        contactPointsSorted.ElementAt(p)[0], contactPointsSorted.ElementAt(p)[1],
                        contactVelocitiesSorted.ElementAt(p)[0], contactVelocitiesSorted.ElementAt(p)[1], contactAnglesSorted.ElementAt(p));
                    Log.WriteLine(line);
                }
                Log.Flush();

                return;
            }
        }


        ConventionalDGField[] GetMeanVelocityFromXDGField(DGField[] EvoVelocity) {

            IList<string> velocityName = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(Solution.NSECommon.VariableNames.LevelSetCG, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(3));
            IReadOnlyDictionary<string, DGField> parameters = this.SolverMainOverride.LsUpdater.Parameters;

            List<ConventionalDGField> velocity = new List<ConventionalDGField>(3);
            for (int i = 0; i < 3; ++i) {

                if (parameters.TryGetValue(velocityName[i], out DGField velocityField)) {
                    velocity.Add((ConventionalDGField)velocityField);
                }
            }
            return velocity.ToArray();

        }


        public Tuple<List<double[]>, List<double[]>, List<double>> ComputeContactLineQuantities() {

            List<double[]> contactPoints = new List<double[]>();
            List<double[]> contactVelocities = new List<double[]>();
            List<double> contactAngles = new List<double>();

            ConventionalDGField[] meanVelocity = GetMeanVelocityFromXDGField(this.CurrentVel);

            var Phi = (LevelSet)LsTrk.LevelSets[0];
            var Phi2 = (LevelSet)LsTrk.LevelSets[1];

            var LevelSetGradient = new VectorField<SinglePhaseField>(SolverMainOverride.GridData.SpatialDimension, Phi.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, (SinglePhaseField)LsTrk.LevelSets[0]);
            SinglePhaseField[] Normals = LevelSetGradient.ToArray();

            var LevelSetGradient2 = new VectorField<SinglePhaseField>(SolverMainOverride.GridData.SpatialDimension, Phi2.Basis, SinglePhaseField.Factory);
            LevelSetGradient2.Gradient(1.0, (SinglePhaseField)LsTrk.LevelSets[1]);
            SinglePhaseField[] Normals1 = LevelSetGradient2.ToArray();

            XQuadSchemeHelper SchemeHelper = this.LsTrk.GetXDGSpaceMetrics(this.LsTrk.SpeciesIdS.ToArray(), this.m_HMForder).XQuadSchemeHelper;

            var gridDat = (GridData)this.SolverMainOverride.GridData;
            var ContactLineVolumeScheme = SchemeHelper.GetContactLineQuadScheme(this.LsTrk.GetSpeciesId("A"), 0);

            CellQuadrature.GetQuadrature(new int[] { 5 }, LsTrk.GridDat,
                ContactLineVolumeScheme.Compile(LsTrk.GridDat, 0),
                delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                    // contact point
                    NodeSet Vnode_l = QR.Nodes;
                    NodeSet Vnode_g = Vnode_l.CloneAs();
                    LsTrk.GridDat.TransformLocal2Global(Vnode_l, Vnode_g, i0);

                    int D = base.SolverMainOverride.Grid.SpatialDimension;
                    for (int d = 0; d < D; d++) {
                        EvalResult[0, 0, d] = Vnode_g[0, d];
                    }

                    // contact line velocity
                    MultidimensionalArray U = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        (meanVelocity[d] as SinglePhaseField).Evaluate(i0, length, QR.Nodes, U.ExtractSubArrayShallow(-1, -1, d));
                    }

                    for (int d = 0; d < D; d++) {
                        EvalResult[0, 0, 2 + d] = U[0, 0, d];
                    }

                    // contact angle
                    MultidimensionalArray normal = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        Normals[d].Evaluate(i0, length, QR.Nodes, normal.ExtractSubArrayShallow(-1, -1, d));
                    }
                    MultidimensionalArray normal1 = MultidimensionalArray.Create(new int[] { 1, 1, D });
                    for (int d = 0; d < D; d++) {
                        Normals1[d].Evaluate(i0, length, QR.Nodes, normal1.ExtractSubArrayShallow(-1, -1, d));
                    }

                    double theta_surf = Math.Atan2(normal[0, 0, 1], normal[0, 0, 0]);
                    double theta_edge = Math.Atan2(normal1[0, 0, 1], normal1[0, 0, 0]);
                    double theta = (theta_surf - theta_edge) * (180 / Math.PI);

                    EvalResult[0, 0, 2 * D] = (theta > 180) ? theta - 180 : theta;
                    //Console.WriteLine("contact angle = {0}", EvalResult[0, 0, 2 * D]);

                },
                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                    int D = SolverMainOverride.Grid.SpatialDimension;
                    for (int i = 0; i < length; i++) {
                        if (ResultsOfIntegration[i, 2 * D] != 0.0) {
                            contactAngles.Add(Math.Abs(ResultsOfIntegration[i, 2 * D]));
                            double[] cp = new double[D];
                            double[] cpV = new double[D];
                            for (int d = 0; d < D; d++) {
                                cp[d] = ResultsOfIntegration[i, d];
                                cpV[d] = ResultsOfIntegration[i, 2 + d];
                            }
                            contactPoints.Add(cp);
                            contactVelocities.Add(cpV);
                        }
                    }
                }
            ).Execute();

            return new Tuple<List<double[]>, List<double[]>, List<double>>(contactPoints, contactVelocities, contactAngles);

        }


    }
}
