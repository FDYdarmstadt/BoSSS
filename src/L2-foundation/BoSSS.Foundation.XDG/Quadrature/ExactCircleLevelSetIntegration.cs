/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF { 
    
    
    public class ExactCircleLevelSetIntegration : IQuadRuleFactory<QuadRule> {

        static bool Rem = false;

        public ExactCircleLevelSetIntegration(int iLs, GridData c, RefElement simplex) {
            if (!Rem) {
                for (int i = 0; i < RADIUS.Length; i++) {
                    Console.WriteLine("ACHTUNG: ExactCircleLevelSetIntegration; Radius = " + RADIUS[i]);
                }
                Rem = true;
            }
            if (simplex.GetType() != typeof(Square))
                throw new ArgumentOutOfRangeException();
            this.RefElement = simplex;
            this.iLevSet = iLs;
            this._Context = c;
        }
        
        public RefElement RefElement {
            get;
            private set; 
        }

        IEnumerable<IChunkRulePair<QuadRule>> IQuadRuleFactory<QuadRule>.GetQuadRuleSet(ExecutionMask mask, int order) {
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

            var ret = new List<IChunkRulePair<QuadRule>>();
            double R = RADIUS[iLevSet];
            
            foreach (var chunk in mask) {
                for (int jCell = chunk.i0; jCell < chunk.JE; jCell++) {

                    double alpha_0, alpha_1;
                    FindAngles(jCell, RADIUS[this.iLevSet], _Context, out alpha_0, out alpha_1);

                    if (alpha_0 == alpha_1) {
                        QuadRule emptyRule = new QuadRule();
                        emptyRule.Nodes = new NodeSet(this.RefElement, new double[1,2]);
                        emptyRule.Weights = MultidimensionalArray.Create(1);
                        emptyRule.OrderOfPrecision = 123;
                        ret.Add(new ChunkRulePair<QuadRule>(new Chunk() { i0  = jCell, Len = 1 }, emptyRule));
                    } else {
                        
                        var da1Drule = this.RefElement.FaceRefElement.GetBruteForceQuadRule(8, order);

                        QuadRule Bogenrule = new QuadRule();

                        var GlobNodes = MultidimensionalArray.Create(da1Drule.NoOfNodes, 2);
                        var LocNodes = MultidimensionalArray.Create(1, da1Drule.NoOfNodes, 2);
                        Bogenrule.Nodes = new NodeSet(this.RefElement, da1Drule.NoOfNodes, 2);

                        for (int k = 0; k < da1Drule.NoOfNodes; k++) {
                            double beta = da1Drule.Nodes[k, 0];
                            double alpha = alpha_0 + (beta + 1)*0.5*(alpha_1 - alpha_0);

                            GlobNodes[k, 0] = R*Math.Cos(alpha);
                            GlobNodes[k, 1] = R*Math.Sin(alpha);

                            if (Position != null) {
                                GlobNodes[k, 0] += Position[0];
                                GlobNodes[k, 1] += Position[1];
                            }
                        }
                        _Context.TransformGlobal2Local(GlobNodes, LocNodes, jCell, 1, 0);
                        Bogenrule.Nodes.Set(LocNodes.ExtractSubArrayShallow(0, -1, -1));

                        var GlobCenter = MultidimensionalArray.Create(1, 2);
                        var LocCenter = MultidimensionalArray.Create(1, 1, 2);
                        _Context.TransformGlobal2Local(GlobCenter, LocCenter, jCell, 1, 0);

                        //double RR = Math.Sqrt((LocNodes[0, 0, 0] - LocCenter[0, 0, 0]).Pow2() + (LocNodes[0, 0, 1] - LocCenter[0, 0, 1]).Pow2());


                        //double metric = (R*(alpha_1 - alpha_0)/2.0)/_Context.ChefBasis.Scaling[jCell];
                        double metric = (R * (alpha_1 - alpha_0) / 2.0) / _Context.Cells.JacobiDet[jCell];
                        //double metric = R;
                        Bogenrule.OrderOfPrecision = 123;
                        Bogenrule.Weights = da1Drule.Weights.CloneAs();
                        Bogenrule.Weights.Scale(metric);
                        Bogenrule.Nodes.LockForever();

                        ret.Add(new ChunkRulePair<QuadRule>(new Chunk() { i0  = jCell, Len = 1 }, Bogenrule));
                         //*/

                        /*
                        var da1Drule = this.RefElement.FaceRefElement.GetBruteForceQuadRule(8, order);

                        QuadRule PswsLin = new QuadRule();

                        var GlobNodes = MultidimensionalArray.Create(da1Drule.NoOfNodes, 2);
                        var LocNodes = MultidimensionalArray.Create(1, da1Drule.NoOfNodes, 2);
                        PswsLin.Nodes = new NodeSet(this.RefElement, da1Drule.NoOfNodes, 2);

                        double x0 = Math.Cos(alpha_0) * R;
                        double y0 = Math.Sin(alpha_0) * R;
                        double x1 = Math.Cos(alpha_1) * R;
                        double y1 = Math.Sin(alpha_1) * R;
                        double dist = Math.Sqrt((x1 - x0).Pow2() + (y1 - y0).Pow2());

                        for(int k = 0; k < da1Drule.NoOfNodes; k++) {
                            double beta = da1Drule.Nodes[k, 0];
                            //double alpha = alpha_0 + (beta + 1) * 0.5 * (alpha_1 - alpha_0);

                            //GlobNodes[k, 0] = R * Math.Cos(alpha);
                            //GlobNodes[k, 1] = R * Math.Sin(alpha);

                            double a = 0.5*(beta + 1);
                            GlobNodes[k, 0] = x0*(1 - a) + x1*a;
                            GlobNodes[k, 1] = y0 * (1 - a) + y1 * a;

                        }
                        _Context.TransformGlobal2Local(GlobNodes, LocNodes, jCell, 1, 0);
                        PswsLin.Nodes.Set(LocNodes.ExtractSubArrayShallow(0, -1, -1));

                        var GlobCenter = MultidimensionalArray.Create(1, 2);
                        var LocCenter = MultidimensionalArray.Create(1, 1, 2);
                        _Context.TransformGlobal2Local(GlobCenter, LocCenter, jCell, 1, 0);

                        //double RR = Math.Sqrt((LocNodes[0, 0, 0] - LocCenter[0, 0, 0]).Pow2() + (LocNodes[0, 0, 1] - LocCenter[0, 0, 1]).Pow2());


                        //double metric = (R*(alpha_1 - alpha_0)/2.0)/_Context.ChefBasis.Scaling[jCell];
                        double metric = (dist / 2.0) / _Context.Cells.JacobiDet[jCell];
                        //double metric = R;
                        PswsLin.OrderOfPrecision = 123;
                        PswsLin.Weights = da1Drule.Weights.CloneAs();
                        PswsLin.Weights.Scale(metric);
                        PswsLin.Nodes.LockForever();

                        ret.Add(new ChunkRulePair<QuadRule>(new Chunk() { i0 = jCell, Len = 1 }, PswsLin));
                         */
                    }
                }
            }

            return ret;
        }

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        public static double[] Position;

        /// <summary>
        /// radius for each level set
        /// </summary>
        public static double[] RADIUS;

        GridData _Context;
        int iLevSet;

        public static void FindAngles(int jCell, double R, GridData _Context, out double alpha0, out double alpha1) {
            BoundingBox BB = new BoundingBox(2);
            _Context.Cells.GetCellBoundingBox(jCell, BB);

            if (Position != null) {
                BB.Min[0] -= Position[0];
                BB.Max[0] -= Position[0];
                BB.Min[1] -= Position[1];
                BB.Max[1] -= Position[1];
            }
            var X = new List<double[]>();
            {
                double x = BB.Min[0];
                double y1 = Math.Sqrt(R.Pow2() - x*x);
                double y2 = -y1;
                if (!double.IsNaN(y1) && BB.Min[1] <= y1 && y1 <= BB.Max[1]) {
                    X.Add(new double[] { x, y1 });
                }
                if (!double.IsNaN(y2) && BB.Min[1] <= y2 && y2 <= BB.Max[1]) {
                    X.Add(new double[] { x, y2 });
                }
            }
            {
                double x = BB.Max[0];
                double y1 = Math.Sqrt(R.Pow2() - x*x);
                double y2 = -y1;
                if (!double.IsNaN(y1) && BB.Min[1] <= y1 && y1 <= BB.Max[1]) {
                    X.Add(new double[] { x, y1 });
                }
                if (!double.IsNaN(y2) && BB.Min[1] <= y2 && y2 <= BB.Max[1]) {
                    X.Add(new double[] { x, y2 });
                }
            }
            {
                double y = BB.Min[1];
                double x1 = Math.Sqrt(R.Pow2() - y*y);
                double x2 = -x1;
                if (!double.IsNaN(x1) && BB.Min[0] <= x1 && x1 <= BB.Max[0]) {
                    X.Add(new double[] { x1, y });
                }
                if (!double.IsNaN(x2) && BB.Min[0] <= x2 && x2 <= BB.Max[0]) {
                    X.Add(new double[] { x2, y });
                }
            }
            {
                double y = BB.Max[1];
                double x1 = Math.Sqrt(R.Pow2() - y*y);
                double x2 = -x1;
                if (!double.IsNaN(x1) && BB.Min[0] <= x1 && x1 <= BB.Max[0]) {
                    X.Add(new double[] { x1, y });
                }
                if (!double.IsNaN(x2) && BB.Min[0] <= x2 && x2 <= BB.Max[0]) {
                    X.Add(new double[] { x2, y });
                }
            }
            
            //
            if (X.Count == 0) {
                alpha0 = 0;
                alpha1 = 0;
                return;
            }
            for (int _i = 0; _i < X.Count; _i++) {
                Debug.Assert(BB.Contains(X[_i]));
            }
            if (X.Count > 2) {
                
                for (int i = 0; i < X.Count; i++) {
                    for (int j = i+1; j < X.Count; j++) {
                        var x_i = X[i];
                        var x_j = X[j];

                        if (GenericBlas.L2Dist(x_i, x_j) <= 1.0e-12) {
                            X.RemoveAt(j);

                            j -= 1;
                        }
                    }
                }

            }



            //if (GenericBlas.L2Dist(X.GetRow(0), X.GetRow(1)) <= 1.0e-12) {
            //    alpha0 = 0;
            //    alpha1 = 0;
            //    return;
            //}



            int FoundNodes = X.Count;
            double alpha_mini = +double.MaxValue;
            double alpha_maxi = -double.MaxValue;
            for (int i = 0; i < FoundNodes; i++) {

                double alpha = Math.Atan2(X[i][1], X[i][0]);
                if (alpha < 0)
                    alpha += Math.PI*2.0;

                alpha_mini = Math.Min(alpha_mini, alpha);
                alpha_maxi = Math.Max(alpha_maxi, alpha);

                var pt = new double[] { R*Math.Cos(alpha), R*Math.Sin(alpha) };
                Debug.Assert(GenericBlas.L2Dist(pt, X[i]) < 1.0e-10);
            }
            Debug.Assert(alpha_maxi >= alpha_mini);


            double Delta_alpha = alpha_maxi - alpha_mini;
            Debug.Assert(Delta_alpha >= 0.0);
            Debug.Assert(Delta_alpha <= Math.PI*2.0);

            if (Delta_alpha < 1.0e-12 || (Delta_alpha - 2.0*Math.PI).Abs() <= 1.0e-12) {
                alpha0 = 0;
                alpha1 = 0;
                return;
            }


            if ((Math.PI*2.0 - Delta_alpha) < Delta_alpha) {
                alpha0 = alpha_maxi;
                alpha1 = alpha0 + (Math.PI*2.0 - Delta_alpha);
                return;
            } else {
                alpha0 = alpha_mini;
                alpha1 = alpha0 + Delta_alpha;
                return;
            }
        }
    }
}
