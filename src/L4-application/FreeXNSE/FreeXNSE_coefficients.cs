using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using static FreeXNSE.Bulkfriction;

namespace FreeXNSE {

    public static class Coefficientnames {
        public static string Bulkviscosityfield = "Bulkviscosityfield";
        public static string Bulkfriction = "Friction_Bulk";
        public static string Contactlinefriction = "Friction_Contactline";
        public static string Contactangle = "Contactangle";
        public static string Surfacetensionfield = "Surfacetensionfield";
        public static string Surfaceshearviscosityfield = "Surfaceshearviscosityfield";
        public static string Surfacedilatationalviscosityfield = "Surfacedilatationalviscosityfield";
    }

    /// <summary>
    /// Coefficient to give the sliplength as a function of the level set
    /// </summary>
    public class Bulkfriction : Coefficient {
        protected int D;

        private double friction;
        private double r;
        private LevelSetUpdater lsUpdate;

        public Bulkfriction(int D, double friction, double r, LevelSetUpdater lsUpdate) {
            this.D = D;
            this.friction = friction;
            this.r = r;
            this.lsUpdate = lsUpdate;
        }

        public override DelCoefficientFactory Factory => CoefficientFactory;

        public override IList<string> CoefficientsNames => new[] { Coefficientnames.Bulkfriction };

        public (string, object)[] CoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            Func<int, double[], double> Beta;

            if(friction < 0 || friction == double.PositiveInfinity) {
                Console.WriteLine("Negative or Infinite Friction parameter, setting negative - meaning no-slip");
                Beta = (iEdge, X) => -1.0;
            } else {
                if(!lsUpdate.Tracker.Equals(lstrk)) throw new ApplicationException();
                _Bulkfriction fric = new _Bulkfriction(friction, r, lstrk.GridDat, lsUpdate.LevelSets[VariableNames.LevelSetCG].DGLevelSet);
                Beta = fric.Evaluate;
            }

            (string, object)[] Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], Beta);
            return Ret;
        }

        internal class _Bulkfriction {

            double beta;
            double r;
            ILevelSet phi;
            GridData grd;
            internal _Bulkfriction(double friction, double r, GridData grd, ILevelSet phi) {
                this.phi = phi.CloneAs();
                beta = friction;
                this.r = r;
                this.grd = grd;
            }

            internal double Evaluate(int jCell, double[] X) {
                int D = grd.SpatialDimension;

                MultidimensionalArray GlobalNode = MultidimensionalArray.Create(1, D);
                GlobalNode.SetRow<double[]>(0, X);

                MultidimensionalArray LocalNode = MultidimensionalArray.Create(1, 1, D);
                grd.TransformGlobal2Local(GlobalNode, LocalNode, jCell, 1, 0);

                MultidimensionalArray Result = MultidimensionalArray.Create(1, 1);
                phi.Evaluate(jCell, 1, new NodeSet(grd.Cells.GetRefElement(jCell), LocalNode.ExtractSubArrayShallow(0,0,-1).To1DArray(), true), Result);

                double scale = Math.Pow(Math.Abs(Result[0,0]), r);
                //Console.WriteLine("At ({0}|{1}): {2}", X[0], X[1], scale);
                return beta * scale;
            }
        }                   
    }


    /// <summary>
    /// 
    /// </summary>
    public class Contactlinefriction : Coefficient {

        protected int D;
        private double friction;
        private double frictionExponent;

        public override IList<string> CoefficientsNames => new[] { Coefficientnames.Contactlinefriction };

        public Contactlinefriction(int D, double friction, double frictionExponent) {
            this.D = D;
            this.friction = friction;
            this.frictionExponent = frictionExponent;
        }
        public override DelCoefficientFactory Factory => CoefficientFactory;

        public (string, object)[] CoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            Func<double[], double>[] Beta = new Func<double[], double>[2];

            if(friction < 0 || friction == double.PositiveInfinity) {
                Console.WriteLine("Negative or Infinite Contactline friction parameter, setting negative - meaning no contactline force");
                Beta[0] = X => -1.0;
                Beta[1] = X => 0.0;
            } else {
                Beta[0] = X => friction;
                Beta[1] = X => frictionExponent;
            }

            (string, object)[] Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], Beta);
            return Ret;
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class Contactangle : Coefficient {

        protected int D;
        private double theta;

        private string name;
        public override IList<string> CoefficientsNames => new[] { name };

        public Contactangle(int D, double theta, int AdvRec = 0) {
            this.D = D;
            this.theta = theta;
            name = Coefficientnames.Contactangle;
            if(AdvRec == -1) {
                name = name + "Rec";
            } else if(AdvRec == 1) {
                name = name + "Adv";
            }
        }
        public override DelCoefficientFactory Factory => CoefficientFactory;

        public (string, object)[] CoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            Func<double[], double> Beta;

            if(theta < 0 || theta > Math.PI) {
                Console.WriteLine("Invalid contactangle, setting to pi/2");
                Beta = X => Math.PI/2.0;
            } else {
                Beta = X => theta;
            }
            
            (string, object)[] Ret = new (string, object)[1];
            Ret[0] = (name, Beta);
            return Ret;
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class Surfacetensionfield : Coefficient {

        public override IList<string> CoefficientsNames => new[] { Coefficientnames.Surfacetensionfield };

        public Surfacetensionfield() {

        }
        public override DelCoefficientFactory Factory => CoefficientFactory;

        public (string, object)[] CoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            Func<double[], double> Sigma;

            Sigma = X => 1.0; // MathNet.Numerics.SpecialFunctions.Erf(X[1]);

            (string, object)[] Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], Sigma);
            return Ret;
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class Surfaceviscosityfields : Coefficient {

        public override IList<string> CoefficientsNames => new[] { Coefficientnames.Surfaceshearviscosityfield, Coefficientnames.Surfacedilatationalviscosityfield };

        public Surfaceviscosityfields() {

        }
        public override DelCoefficientFactory Factory => CoefficientFactory;

        public (string, object)[] CoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            Func<double[], double> Lambda;
            Func<double[], double> Mu;

            Lambda = X => 2.0; // Dilatationalviscosity
            Mu = X => 1.0; // Shearviscosity


            (string, object)[] Ret = new (string, object)[2];
            Ret[0] = (CoefficientsNames[0], Mu);
            Ret[1] = (CoefficientsNames[1], Lambda);
            return Ret;
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class Bulkviscosityfield : Coefficient {

        public override IList<string> CoefficientsNames => new[] { Coefficientnames.Bulkviscosityfield };

        private LevelSetUpdater lsUpdate;
        private double viscosityscaling;

        public Bulkviscosityfield(double viscosityscaling, LevelSetUpdater lsUpdate) {
            this.lsUpdate = lsUpdate;
            this.viscosityscaling = viscosityscaling;
        }
        public override DelCoefficientFactory Factory => CoefficientFactory;

        public (string, object)[] CoefficientFactory(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {
            Func<int, double[], double> Mu;

            _Bulkviscosity visc = new _Bulkviscosity(lstrk.GridDat, lsUpdate.LevelSets[VariableNames.LevelSetCG].DGLevelSet, viscosityscaling);
            if(viscosityscaling == 0.0) {
                Mu = (jCell, X) => 1.0;
            } else {
                Mu = visc.Evaluate; // Shearviscosity
            }

            (string, object)[] Ret = new (string, object)[1];
            Ret[0] = (CoefficientsNames[0], Mu);
            return Ret;
        }

        internal class _Bulkviscosity {

            ILevelSet phi;
            GridData grd;
            double viscosityscaling;
            internal _Bulkviscosity(GridData grd, ILevelSet phi, double viscosityscaling) {
                this.phi = phi.CloneAs();
                this.grd = grd;
                this.viscosityscaling = viscosityscaling;

                ContactPoints.Clear();

                int D = grd.SpatialDimension;
                LevelSetTracker trk = new LevelSetTracker(grd, CutCellQuadratureMethod.Saye, 1, new string[] { "A", "B" }, phi);
                trk.UpdateTracker(0.0);

                XQuadSchemeHelper SchemeHelper = trk.GetXDGSpaceMetrics(trk.SpeciesIdS.ToArray(), 0).XQuadSchemeHelper;
                EdgeQuadratureScheme SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(trk.GetSpeciesId("A"), 0);

                var QuadDom = SurfaceElement_Edge.Domain;
                var boundaryEdge = grd.GetBoundaryEdgeMask().GetBitMask();
                var boundaryCutEdge = QuadDom.Intersect(new EdgeMask(grd, boundaryEdge, MaskType.Geometrical));

                var factory = trk.GetXDGSpaceMetrics(trk.SpeciesIdS.ToArray(), 0).XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, trk.GridDat.Grid.RefElements[0]);
                SurfaceElement_Edge = new EdgeQuadratureScheme(factory, boundaryCutEdge);

                EdgeQuadrature.GetQuadrature(new int[] { D }, trk.GridDat,
                    SurfaceElement_Edge.Compile(trk.GridDat, 0),
                    delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {
                        // contact point
                        NodeSet Enode_l = QR.Nodes;
                        int trf = trk.GridDat.Edges.Edge2CellTrafoIndex[i0, 0];
                        NodeSet Vnode_l = Enode_l.GetVolumeNodeSet(trk.GridDat, trf, false);
                        NodeSet Vnode_g = Vnode_l.CloneAs();
                        int cell = trk.GridDat.Edges.CellIndices[i0, 0];
                        trk.GridDat.TransformLocal2Global(Vnode_l, Vnode_g, cell);
                        //Console.WriteLine("contact point: ({0},{1})", Vnode_g[0, 0], Vnode_g[0, 1]);

                        int D = trk.GridDat.SpatialDimension;
                        for(int d = 0; d < D; d++) {
                            EvalResult[0, 0, d] = Vnode_g[0, d];
                        }                        
                    },
                    delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < length; i++) {
                            double[] cp = new double[D];
                            for(int d = 0; d < D; d++) {
                                cp[d] = ResultsOfIntegration[i, d];
                            }
                            ContactPoints.Add(cp);
                        }                        
                    }
                ).Execute();
            }

            List<double[]> ContactPoints = new List<double[]>();

            internal double Evaluate(int jCell, double[] X) {               
                double dist;
                if(!ContactPoints.IsNullOrEmpty()) {
                    dist = ContactPoints.Select(x => X.L2Dist(x)).Min();
                } else {
                    dist = 1.0;
                }
                double R = Math.Pow(Math.Abs(dist), viscosityscaling);
                return R;
            }
        }
    }

    [Serializable]
    [DataContract]
    public class ConstantGravity : Forcing {

        [DataMember]
        double[] direction;

        public ConstantGravity(double[] direction) {
            this.direction = direction.Normalize();
        }

        public override bool equals(Forcing other) {

            if(other is ConstantGravity f) {
                return this.direction.SequenceEqual(f.direction);
            }
            
            return false;
        }

        public override double[] Evaluate(double[] X, double t) {
            return direction;
        }
    }
}
