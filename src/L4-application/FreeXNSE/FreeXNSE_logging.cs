using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using ilPSP.Tracing;
using ilPSP;
using log4net;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;
using Newtonsoft.Json;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using ilPSP.Utils;

namespace FreeXNSE {

    /// <summary>
    /// base-class for post-processing modules in FreeXNSE
    /// </summary>
    public abstract class FreeXNSEinSituPostProcessingModule : InSituPostProcessingModule {

        /// <summary>
        /// reference to solver application class
        /// </summary>
        protected FreeXNSE SolverMainOverride {
            get {
                return (FreeXNSE)base.SolverMain;
            }
        }

        /// <summary>
        /// control object
        /// </summary>
        new protected FreeXNSE_Control Control {
            get {
                return (FreeXNSE_Control)(base.Control);
            }
        }


        /// <summary>
        /// current velocity solution
        /// </summary>
        protected DGField[] CurrentVel {
            get {
                int D = this.SolverMainOverride.GridData.SpatialDimension;

                var ret = this.SolverMainOverride.CurrentState.Fields.Take(D).Select(f => ((XDGField)f).GetSpeciesShadowField("A")).ToArray();
                for(int d = 0; d < D; d++) {
                    if(ret[d].Identification != VariableNames.Velocity_d(d) + "#A")
                        throw new ApplicationException("Unable to identify velocity fields.");
                }

                return ret;
            }
        }

        /// <summary>
        /// current velocity solution
        /// </summary>
        protected DGField CurrentPressure {
            get {
                int D = this.SolverMainOverride.GridData.SpatialDimension;

                var ret = ((XDGField)this.SolverMainOverride.CurrentState.Fields.ElementAt(D)).GetSpeciesShadowField("A");
                if(ret.Identification != VariableNames.Pressure + "#A")
                    throw new ApplicationException("Unable to identify pressure field.");

                return ret;
            }
        }


        /// <summary>
        /// the level-set which represents the fluid interface
        /// </summary>
        protected LevelSet LevSet {
            get {
                return (LevelSet)(this.LsTrk.LevelSets[0]);
            }
        }


        /// <summary>
        /// Cut-Cell quadrature order used for the flow solver
        /// </summary>
        protected int m_HMForder {
            get {
                return this.SolverMainOverride.QuadOrder();
            }
        }


        /// <summary>
        /// les gird
        /// </summary>
        protected GridData GridData {
            get {
                return (GridData)(this.SolverMainOverride.GridData);
            }
        }



    }


    /// <summary>
    /// Used by <see cref="HeatedWall"/>, <see cref="TwoPhaseCouetteFlow"/>, etc.
    /// </summary>
    [Serializable]
    public class ContactLineLogging : FreeXNSEinSituPostProcessingModule {

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
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", "#timestep", "time", "contact-pointX", "contact-pointY", "contact-VelocityX", "contact-VelocityY", "contact-Velocity", "contact-angle");
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
                List<double> contactLineVelocities = new List<double>();
                List<double> contactAngles = new List<double>();

                var returnValues = this.ComputeContactLineQuantities();
                contactPoints = returnValues.Item1;
                contactVelocities = returnValues.Item2;
                contactLineVelocities = returnValues.Item3;
                contactAngles = returnValues.Item4;


                List<double[]> contactPointsSorted = new List<double[]>();
                List<double[]> contactVelocitiesSorted = new List<double[]>();
                List<double> contactLineVelocitiesSorted = new List<double>();
                List<double> contactAnglesSorted = new List<double>();

                if(contactPointsRef == null) {
                    contactPointsRef = contactPoints;
                    contactPointsSorted = contactPoints;
                    contactVelocitiesSorted = contactVelocities;
                    contactLineVelocitiesSorted = contactLineVelocities;
                    contactAnglesSorted = contactAngles;
                } else {
                    // sort
                    double eps = 1e-7;
                    do {
                        contactPointsSorted.Clear();
                        contactVelocitiesSorted.Clear();
                        contactLineVelocitiesSorted.Clear();
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
                                    contactLineVelocitiesSorted.Add(contactLineVelocities.ElementAt(i));
                                    contactAnglesSorted.Add(contactAngles.ElementAt(i));
                                    break;
                                }
                            }
                        }

                    } while(contactPointsSorted.Count != contactPointsRef.Count && eps < 1.0);

                    contactPointsRef = contactPointsSorted;
                }
                for(int p = 0; p < contactAnglesSorted.Count; p++) {
                    string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", TimestepNo, phystime,
                        contactPointsSorted.ElementAt(p)[0], contactPointsSorted.ElementAt(p)[1],
                        contactVelocitiesSorted.ElementAt(p)[0], contactVelocitiesSorted.ElementAt(p)[1], contactLineVelocitiesSorted.ElementAt(p), contactAnglesSorted.ElementAt(p));
                    Log.WriteLine(line);
                }
                Log.Flush();

                return;
            }
        }

        public Tuple<List<double[]>, List<double[]>, List<double>, List<double>> ComputeContactLineQuantities() {

            int D = this.GridData.SpatialDimension;

            if(D != 2) {
                throw new NotImplementedException();
            }

            List<double[]> contactPoints = new List<double[]>();
            List<double[]> contactVelocities = new List<double[]>();
            List<double> contactLineVelocities = new List<double>();
            List<double> contactAngles = new List<double>();

            DGField[] Velocity = this.CurrentVel;
            var Phi = this.LevSet;
            var LevelSetGradient = new VectorField<SinglePhaseField>(SolverMainOverride.GridData.SpatialDimension, Phi.Basis, SinglePhaseField.Factory);
            LevelSetGradient.Gradient(1.0, Phi);
            SinglePhaseField[] Normals = LevelSetGradient.ToArray();

            if(this.LsTrk.NoOfLevelSets > 1)
                throw new NotSupportedException("todo -- maybe missing level-set intersection line contributions.");

            
            XQuadSchemeHelper SchemeHelper = this.LsTrk.GetXDGSpaceMetrics(this.m_HMForder).XQuadSchemeHelper;
            EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper
                .Get_SurfaceElement_EdgeQuadScheme(this.LsTrk.GetSpeciesId("A"), this.LsTrk.GetSpeciesId("B"), 0)
                .Restrict(this.LsTrk.GridDat.GetBoundaryEdgeMask());


            EdgeQuadrature.GetQuadrature(new int[] { D + D + 1 + 1 }, LsTrk.GridDat,
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
                    for(int d = 0; d < D; d++) {
                        EvalResult[0, 0, d] = Vnode_g[0, d];
                    }

                    // contact line velocity
                    MultidimensionalArray U_IN = MultidimensionalArray.Create(new int[] { 1, QR.NoOfNodes, D });
                    MultidimensionalArray U_OUT = MultidimensionalArray.Create(new int[] { 1, QR.NoOfNodes, D });
                    for(int d = 0; d < D; d++) {
                        Velocity[d].EvaluateEdge(i0, length, QR.Nodes, U_IN.ExtractSubArrayShallow(-1, -1, d), U_OUT.ExtractSubArrayShallow(-1, -1, d), null, null, null, null, 0, 0.0);
                    }

                    for(int d = 0; d < D; d++) {
                        EvalResult[0, 0, D + d] = U_IN[0, 0, d];
                    }

                    // surface normal
                    MultidimensionalArray normal_IN = MultidimensionalArray.Create(new int[] { 1, QR.NoOfNodes, D });
                    MultidimensionalArray normal_OUT = MultidimensionalArray.Create(new int[] { 1, QR.NoOfNodes, D });
                    for(int d = 0; d < D; d++) {
                        Normals[d].EvaluateEdge(i0, length, QR.Nodes, normal_IN.ExtractSubArrayShallow(-1, -1, d), normal_OUT.ExtractSubArrayShallow(-1, -1, d));
                    }
                    Vector NormalSurf = normal_IN.GetRowPt(0, 0);
                    NormalSurf.Normalize();

                    // edge normal
                    Vector NormalEdg = LsTrk.GridDat.Edges.NormalsForAffine.GetRow(i0);
                    NormalEdg.Normalize();

                    // contactline normal
                    double[,] Pedg = FreeXNSE_utils.Projection(NormalEdg);
                    double[] NormalCL = new double[D];
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            NormalCL[d1] += Pedg[d1, d2] * NormalSurf[d2];
                        }
                    }
                    NormalCL.Normalize();

                    // contact line normal velocity
                    for(int d = 0; d < D; d++) {
                        EvalResult[0, 0, D + D] += U_IN[0, 0, d] * NormalCL[d];
                    }

                    // contact angle
                    double CosTheta = -NormalSurf.InnerProd(NormalEdg);
                    double Theta = Math.Acos(CosTheta) * 180 / Math.PI;

                    EvalResult[0, 0, D + D + 1] = Theta;
                },
                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                    int D = SolverMainOverride.Grid.SpatialDimension;
                    for(int i = 0; i < length; i++) {
                        if(ResultsOfIntegration[i, D + D + 1] != 0.0) {
                            contactAngles.Add(ResultsOfIntegration[i, D + D + 1]);
                            double[] cp = new double[D];
                            double[] cpV = new double[D];
                            for(int d = 0; d < D; d++) {
                                cp[d] = ResultsOfIntegration[i, d];
                                cpV[d] = ResultsOfIntegration[i, D + d];
                            }
                            contactPoints.Add(cp);
                            contactVelocities.Add(cpV);
                            contactLineVelocities.Add(ResultsOfIntegration[i, D + D]);
                        }
                    }
                }
            ).Execute();

            //Console.WriteLine("Count of contactpoints : {0}", contactAngles.Count);

            return new Tuple<List<double[]>, List<double[]>, List<double>, List<double>>(contactPoints, contactVelocities, contactLineVelocities, contactAngles);

        }


    }
}
