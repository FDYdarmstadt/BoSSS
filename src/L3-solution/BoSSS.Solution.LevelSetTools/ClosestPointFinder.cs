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
using BoSSS.Platform;
using BoSSS.Foundation.Grid;
using ilPSP.Tracing;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Foundation;
using ilPSP.Utils;
using System.IO;
using ilPSP;
using System.Collections;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.LevelSetTools {


    /// <summary>
    /// projects point from the Near-Field onto some zero-level set. 
    /// </summary>
    public class ClosestPointFinder {

        /*
        static int ccccc = 0;

        public static void PrintPoint(params double[] P) {
            Console.Write(ccccc);
            Console.Write(" ");
            ccccc++;
            for(int i = 0; i < P.Length; i++) {
                Console.Write(P[i].ToStringDot());
                if(i < P.Length - 1)
                    Console.Write(" ");
                else
                    Console.WriteLine();
            }
        }
         */

        /// <summary>
        /// Diagnostic function/debbugging, to be removed.
        /// </summary>
        public static void PlottiPunkti(MultidimensionalArray Pts, string FileName, MultidimensionalArray wgt = null, System.Collections.BitArray OutsiderMarker = null) {
            if(Pts.Dimension != 2)
                throw new ArgumentException();
            int NoOfPoints = Pts.GetLength(0);
            int D = Pts.GetLength(1);

            using(var outStr = new System.IO.StreamWriter(FileName)) {

                for(int k = 0; k < NoOfPoints; k++) {
                    outStr.Write(k);
                    outStr.Write('\t');

                    for(int d = 0; d < D; d++) {
                        outStr.Write(Pts[k, d].ToStringDot());
                        if(d < (D - 1))
                            outStr.Write('\t');
                    }

                    if(wgt != null) {
                        outStr.Write("\t");
                        outStr.Write(wgt[k].ToStringDot());
                    }

                    if(OutsiderMarker != null) {
                        outStr.Write("\t");
                        outStr.Write((OutsiderMarker[k] ? 1.0 : 0.0).ToStringDot());
                    }

                    outStr.WriteLine();
                }

                outStr.Flush();
                outStr.Close();
            }
        }


        /// <summary>
        /// Diagnostic function/debbugging, to be removed.
        /// </summary>
        public void PlotStays(string FileName, System.Collections.BitArray Marker) {

            List<Tuple<int,double[]>> Points = new List<Tuple<int, double[]>>();


            int NoOfPoints = this.X0_global.GetLength(0);
            int D = this.X0_global.GetLength(1);

            int Counter = 0;
            for(int k = 0; k < NoOfPoints; k++) {
                if(Marker[k]) {
                    int jsub_ofX = this.CellIndexMap[k];
                    int nn_ofX = this.NodeIndexMap[k];

                    double[] stPoint = this.X0_global.GetRow(k);
                    double[] enPoint = this.X_global.ExtractSubArrayShallow(jsub_ofX, nn_ofX, -1).To1DArray();

                    int PointsPerStay = 100; // one line ('stay') will be approximated with this number of points.
                    double delta = 1.0 / (PointsPerStay - 1);

                    for(int i = 0; i < PointsPerStay; i++) {
                        double s = delta * i;
                        double[] pt = new double[D];
                        pt.AccV(s, stPoint);
                        pt.AccV(1 - s, enPoint);
                        Points.Add(new Tuple<int, double[]>(Counter, pt));
                    }
                    
                    Counter++;
                }
            }

            using(var outStr = new System.IO.StreamWriter(FileName)) {

                foreach(var tt in Points) {
                    outStr.Write(tt.Item1);
                    outStr.Write('\t');

                    for(int d = 0; d < D; d++) {
                        outStr.Write(tt.Item2[d].ToStringDot());
                        if(d < (D - 1))
                            outStr.Write('\t');
                    }

                    

                    outStr.WriteLine();
                }

                outStr.Flush();
                outStr.Close();
            }


        }


        /// <summary>
        /// Diagnostic function/debbugging, to be removed.
        /// </summary>
        public void PlotOlt2NewStays(string FileName, MultidimensionalArray oldPoints, MultidimensionalArray newPoints) {

            List<Tuple<int,double[]>> Points = new List<Tuple<int, double[]>>();


            int NoOfPoints = this.X0_global.GetLength(0);
            int D = this.X0_global.GetLength(1);

            int Counter = 0;
            for(int k = 0; k < NoOfPoints; k++) {

                int jsub_ofX = this.CellIndexMap[k];
                int nn_ofX = this.NodeIndexMap[k];

                double[] stPoint = oldPoints.GetRow(k);
                double[] enPoint = newPoints.GetRow(k);

                int PointsPerStay = 100; // one line ('stay') will be approximated with this number of points.
                double delta = 1.0 / (PointsPerStay - 1);

                for(int i = 0; i < PointsPerStay; i++) {
                    double s = delta * i;
                    double[] pt = new double[D];
                    pt.AccV(s, stPoint);
                    pt.AccV(1 - s, enPoint);
                    Points.Add(new Tuple<int, double[]>(Counter, pt));
                }

                Counter++;

            }

            using(var outStr = new System.IO.StreamWriter(FileName)) {

                foreach(var tt in Points) {
                    outStr.Write(tt.Item1);
                    outStr.Write('\t');

                    for(int d = 0; d < D; d++) {
                        outStr.Write(tt.Item2[d].ToStringDot());
                        if(d < (D - 1))
                            outStr.Write('\t');
                    }

                    outStr.WriteLine();
                }

                outStr.Flush();
                outStr.Close();
            }
        }



        /// <summary>
        /// Diagnostic function/debbugging, to be removed.
        /// </summary>
        public void PlotCellAssociation(string FileName, System.Collections.BitArray Marker) {

            List<Tuple<int,double[]>> Points = new List<Tuple<int, double[]>>();


            int NoOfPoints = this.X0_global.GetLength(0);
            int D = this.X0_global.GetLength(1);

            int[] jSub_2_jCell = this.sgrd.SubgridIndex2LocalCellIndex;
            int[] I = this.x0I2SgrdCell;
            int JSUB = sgrd.LocalNoOfCells;
                
            int Counter = 0;
            for(int jsub = 0; jsub < JSUB; jsub++) {
                int I0 = I[jsub];
                int IE = I[jsub + 1];
                int K = IE - I0;
                int jCell = jSub_2_jCell[jsub];

                double[] CellCenter = _GridData.Cells.CellCenter.GetRow(jCell);
                
                for(int k = I0; k < IE; k++) {
                    if(Marker[k]) {

                        double[] stPoint = this.X0_global.GetRow(k);



                        int PointsPerStay = 100; // one line ('stay') will be approximated with this number of points.
                        double delta = 1.0 / (PointsPerStay - 1);

                        for(int i = 0; i < PointsPerStay; i++) {
                            double s = delta * i;
                            double[] pt = new double[D];
                            pt.AccV(s, stPoint);
                            pt.AccV(1 - s, CellCenter);
                            Points.Add(new Tuple<int, double[]>(Counter, pt));
                        }

                        Counter++;
                    }
                }
            }

            using(var outStr = new System.IO.StreamWriter(FileName)) {

                foreach(var tt in Points) {
                    outStr.Write(tt.Item1);
                    outStr.Write('\t');

                    for(int d = 0; d < D; d++) {
                        outStr.Write(tt.Item2[d].ToStringDot());
                        if(d < (D - 1))
                            outStr.Write('\t');
                    }



                    outStr.WriteLine();
                }

                outStr.Flush();
                outStr.Close();
            }


        }

        /// <summary>
        /// Original points/nodes (i.e. the input data) in global coordinates <br/>
        ///  - 1st index: subgrid cell index, into subgrid <see cref="sgrd"/> <br/>
        ///  - 2nd index: node index <br/>
        ///  - 3rd index: spatial dimension
        /// </summary>
        public MultidimensionalArray X_global;


        /// <summary>
        /// Distance of original points/nodes (<see cref="X_global"/>) to zero-level-set surface;<br/>
        ///  - 1st index: subgrid cell index, correlates with 1st index of <see cref="X_global"/><br/>
        ///  - 2nd index: node index, correlates with 2nd index of <see cref="X_global"/> <br/>
        /// </summary>
        public MultidimensionalArray Distance {
            get {
                int JSUB = sgrd.LocalNoOfCells; // no of subgrid cells
                Debug.Assert(JSUB == X_global.GetLength(0));
                int K = X_global.GetLength(1); // no of nodes
                int D = X_global.GetLength(2);

                MultidimensionalArray R = MultidimensionalArray.Create(JSUB, K);

                int[] I = this.x0I2SgrdCell;
                for (int jsub = 0; jsub < JSUB; jsub++) {
                    int I0 = I[jsub];
                    int IE = I[jsub + 1];
                    int L = IE - I0;

                    for (int nn = 0; nn < L; nn++) {

                        int jsub_ofX = this.CellIndexMap[I0 + nn];
                        int nn_ofX = this.NodeIndexMap[I0 + nn];

                        double dist = 0;
                        for (int d = 0; d < D; d++) {
                            dist += (this.X0_global[nn + I0, d] - this.X_global[jsub_ofX, nn_ofX, d]).Pow2();
                        }
                        dist = Math.Sqrt(dist);

                        R[jsub_ofX, nn_ofX] = dist*Math.Sign(x_sign[jsub_ofX, nn_ofX]);
                    }

                }

                return R;
            }
        }

        /// <summary>
        /// see <see cref="EvaluateAtCp(ScalarFunctionEx)"/>;
        /// </summary>
        public MultidimensionalArray EvaluateAtCp(DGField f) {
            return EvaluateAtCp(f.Evaluate);
        }

        /// <summary>
        /// Evaluates <see cref="f"/> at the closest points and returns the result in the input node order;
        /// </summary>
        /// <param name="f"></param>
        /// <returns>
        /// Values of <paramref name="f"/> at the interface;<br/>
        ///  - 1st index <em>j</em>: subgrid cell index, correlates with 1st index of <see cref="X_global"/><br/>
        ///  - 2nd index <em>k</em>: node index, correlates with 2nd index of <see cref="X0_global_Resorted"/> <br/>
        /// content: value of <paramref name="f"/> at closest point (to <see cref="X0_global_Resorted"/>[<em>j</em>,<em>k</em>,-]
        /// </returns>
        public MultidimensionalArray EvaluateAtCp(ScalarFunctionEx f) {
            int JSUB = sgrd.LocalNoOfCells; // no of subgrid cells
            Debug.Assert(JSUB == X_global.GetLength(0));
            int K = X_global.GetLength(1); // no of nodes
            int D = X_global.GetLength(2);
            var ctx = _GridData;
            var jsub_2_jcell = sgrd.SubgridIndex2LocalCellIndex;

            MultidimensionalArray R = MultidimensionalArray.Create(JSUB, K);

            int[] I = this.x0I2SgrdCell;
            for (int jsub = 0; jsub < JSUB; jsub++) {
                int I0 = I[jsub];
                int IE = I[jsub + 1];
                int L = IE - I0;

                if (L <= 0)
                    continue;

                var X0_Nodes = this.X0_local.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] {IE - 1, D - 1});
                
                int jCell = jsub_2_jcell[jsub];
                RefElement Kref = ctx.Cells.GetRefElement(jCell);

                MultidimensionalArray RR = MultidimensionalArray.Create(1, L);
                NodeSet nodes = new NodeSet(Kref, X0_Nodes);
                f(jCell, 1, nodes, RR);
                
                for (int nn = 0; nn < L; nn++) {

                    int jsub_ofX = this.CellIndexMap[I0 + nn];
                    int nn_ofX = this.NodeIndexMap[I0 + nn];

                    R[jsub_ofX, nn_ofX] = RR[0, nn];
                }

            }

            return R;
        }


        /// <summary>
        /// Closest points/nodes (for each original point <see cref="X_global"/>) in global coordinates; note that the ordering 
        /// does not correlate with ordering in <see cref="X_global"/>. For sorting according to original ordering, see <see cref="X0_global_Resorted"/>.
        /// <br/>
        ///  - 1st index: node index (over all cells) <br/>
        ///  - 2nd index: spatial dimension 
        /// </summary>
        public MultidimensionalArray X0_global;

        /// <summary>
        /// An error measure for the closest points defined as 
        /// \f[
        ///   \| \varphi(\vec{x}_0) \| 
        ///      + 
        ///   \| \|
        ///   \frac{ \langle \vec{x} - \vec{x}_0 , ( \nabla \varphi ) ( \vec{x}_0 ) \rangle  }
        ///        { \| \vec{x} - \vec{x}_0 \| \| ( \nabla \varphi ) ( \vec{x}_0 ) \| }
        ///   \| - 1 \|
        /// \f]
        /// The ordering correlates with <see cref="X0_global"/>.
        /// </summary>
        public MultidimensionalArray X0_ErrorMeasure;


        MultidimensionalArray m_X0_global_Resorted;


        /// <summary>
        /// closest points in global coordinates, sorted according to input nodes, i.e. indexing correlates with <see cref="X_global"/><br/>
        /// 1st index: subgrid cell index, into subgrid <see cref="sgrd"/> <br/>
        /// 2nd index: node index <br/>
        /// 3rd index: spatial dimension
        /// </summary>
        public MultidimensionalArray X0_global_Resorted {
            get {
                if (m_X0_global_Resorted == null) {
                    m_X0_global_Resorted = MultidimensionalArray.Create(this.X_global.Lengths);
                    

                    int JSUB = sgrd.LocalNoOfCells; // no of subgrid cells
                    Debug.Assert(JSUB == X_global.GetLength(0));
                    int K = X_global.GetLength(1); // no of nodes
                    int D = X_global.GetLength(2);

                    

                    int[] I = this.x0I2SgrdCell;
                    for (int jsub = 0; jsub < JSUB; jsub++) {
                        int I0 = I[jsub];
                        int IE = I[jsub + 1];
                        int L = IE - I0;

                        for (int nn = 0; nn < L; nn++) {

                            int jsub_ofX = this.CellIndexMap[I0 + nn];
                            int nn_ofX = this.NodeIndexMap[I0 + nn];

                            for (int d = 0; d < D; d++) {
                                m_X0_global_Resorted[jsub_ofX, nn_ofX, d] =  this.X0_global[nn + I0, d];
                            };

                            
                        }

                    }
                }
                return m_X0_global_Resorted;
            }

        }

        /// <summary>
        /// Closest points/nodes in local coordinates <br/>
        ///  - 1st index: node index (over all cells) <br/>
        ///  - 2nd index: spatial dimension 
        /// </summary>
        public MultidimensionalArray X0_local;


        /// <summary>
        /// Which closest point (<see cref="X0_global"/>[i,-]) maps to which input point (<see cref="X_global"/>[j,k,-])?
        /// index: node index (over all cells), correlates with 1st index of <see cref="X0_global"/>, resp. <see cref="X0_local"/><br/>
        /// content: subgrid cell index;
        /// </summary>
        /// <remarks>
        /// <see cref="X0_global"/>[i,-] is the closest point to <see cref="X_global"/>[<see cref="CellIndexMap"/>[i],<see cref="NodeIndexMap"/>[i],-]
        /// </remarks>
        public int[] CellIndexMap;

        /// <summary>
        /// <see cref="CellIndexMap"/>.
        /// </summary>
        public int[] NodeIndexMap;

        /// <summary>
        /// 
        /// </summary>
        public SubGrid sgrd;

        /// <summary>
        /// Assignment of closest points to cells.<br/>
        /// index: subgrid cell index; length: local number of cells in subgrid + 1;
        /// </summary>
        public int[] x0I2SgrdCell;


        MultidimensionalArray x_sign;

        GridData _GridData;

        /// <summary>
        /// Constructor. The call performs the computation of the closes points.
        /// </summary>
        public ClosestPointFinder(LevelSetTracker tracker, int iLevSet, SubGrid sgrd, IEnumerable<NodeSet> _Nodes, int MaxIter = 50, int MinIter = 8, double ConvThreshold = 1.0e-13, bool DiagnosticOutput = false) {
            using (var tr = new FuncTrace()) {
                _GridData = tracker.GridDat;
                LevelSet LevSet = (LevelSet)tracker.LevelSets[iLevSet];  // cast is a hack
                var RefElements = _GridData.Grid.RefElements;
                if (_Nodes.Count() != RefElements.Length)
                    throw new ArgumentException("Nodes must correspond to reference elements.", "_Nodes");
                for(int iKref = 0; iKref < RefElements.Length; iKref++) {
                    if (!object.ReferenceEquals(_Nodes.ElementAt(iKref).RefElement, RefElements[iKref]))
                        throw new ArgumentException("Nodes must correspond to reference elements.", "_Nodes");
                }

                int D = _GridData.Grid.SpatialDimension;
                int[] jSub_2_jCell = sgrd.SubgridIndex2LocalCellIndex;
                int JSUB = sgrd.LocalNoOfCells;
                this.sgrd = sgrd;
                int[] jCell_to_jSub = sgrd.LocalCellIndex2SubgridIndex;

                int NNmax = _Nodes.Max(Nodes =>  Nodes.GetLength(0)); // Maximum number of nodes
                int NNglobal = 0;
                for (int jSub = 0; jSub < JSUB; jSub++) {

                    int jCell = jSub_2_jCell[jSub];
                    int iKref = _GridData.Cells.GetRefElementIndex(jCell);
                    var Nodes = _Nodes.ElementAt(iKref);
                    int NN = Nodes.GetLength(0); // number of nodes
                    NNglobal += NN;
                }

                this.X0_ErrorMeasure = MultidimensionalArray.Create(NNglobal);
                this.X0_local = MultidimensionalArray.Create(NNglobal, D);
                this.X0_global = MultidimensionalArray.Create(NNglobal, D);
                x_sign = MultidimensionalArray.Create(JSUB, NNmax);
                
                var RefElms = _GridData.Grid.RefElements;

                var X_global = MultidimensionalArray.Create(JSUB, NNmax, D); // 
                this.X_global = X_global;
                

                int[] I = new int[JSUB + 1];
                this.x0I2SgrdCell = I;

                this.NodeIndexMap = new int[NNglobal];
                this.CellIndexMap = new int[NNglobal];
                int[] ActualjSub = new int[NNglobal];
                BitArray OutsideSgrdMarker = new BitArray(NNglobal);
                BitArray OutsiderMarker = new BitArray(NNglobal);
                
                //var stw = new System.IO.StreamWriter("myIterator2.txt");

                // ==================
                // set initial value
                // ==================
                for (int jSub = 0; jSub < JSUB; jSub++) {
                    
                    int jCell = jSub_2_jCell[jSub];
                    int iKref = _GridData.Cells.GetRefElementIndex(jCell);
                    RefElement Kref = RefElements[iKref];
                    NodeSet Nodes = _Nodes.ElementAt(iKref);
                    int NN = Nodes.GetLength(0); // number of nodes

                    I[jSub + 1] = I[jSub] + NN;
                    for (int k = I[jSub]; k < I[jSub + 1]; k++) {
                        this.CellIndexMap[k] = jSub;
                        this.NodeIndexMap[k] = k - I[jSub];
                        ActualjSub[k] = jSub;
                    }

                    var X0_loc = this.X0_local.ExtractSubArrayShallow(new int[] { I[jSub], 0 }, new int[] { I[jSub + 1] - 1, D - 1 });
                    var X0_glb = this.X0_global.ExtractSubArrayShallow(new int[] { I[jSub], 0 }, new int[] { I[jSub + 1] - 1, D - 1 });
                    //var _X0_glb = MultidimensionalArray.Create(I[jSub + 1] - I[jSub], D); 

                    bool Cont = X0_glb.IsContinious;
                    var X_glb = X_global.ExtractSubArrayShallow(jSub, -1, -1);


                    X0_loc.Set(Nodes);

                    _GridData.TransformLocal2Global(new NodeSet(Kref, X0_loc), X0_glb, jCell);

                    X_global.ExtractSubArrayShallow(jSub, -1, -1).Set(X0_glb);

                    //for (int k = I[jSub]; k < I[jSub + 1]; k++) {
                    //    stw.Write("0 ");
                    //    stw.Write(X0_glb[k - I[jSub], 0] + " ");
                    //    stw.Write(X0_glb[k - I[jSub], 1] + " ");
                    //    stw.WriteLine();
                    //}
                    
                    LevSet.Evaluate(jCell, 1, Nodes, x_sign.ExtractSubArrayShallow(new int[] { jSub, 0 }, new int[] { jSub, NN - 1 }));
                }
                


                var _X0 = new double[D];     // X0 in local coordinates (currently assigned cell)
                var _X0_G = new double[D];   // X0 in global/physical coordinates
                var _X0_L = new double[D];   // X0 in local coordinates of some other cell
                var _X0_L2 = new double[D];  // X0 in local coordinates of some other cell
                //var _X_G = new double[D];    // X: original node/input in physical coordinates
                var dummyPt1 = new double[D];
                var dummyPt2 = new double[D];

                
                double prevRadiusError = double.MaxValue;
                for(int iter = 0; iter < MaxIter; iter++) {

                    if (DiagnosticOutput) {
                        PlottiPunkti(this.X0_global, "X0_iter-" + iter + ".csv", this.X0_ErrorMeasure, OutsideSgrdMarker);
                    }
                    MultidimensionalArray X0_global_before, X0_global_mitte;
                    if(iter == MaxIter - 1) {
                        X0_global_before = this.X0_global.CloneAs();
                        X0_global_mitte = MultidimensionalArray.Create(X0_global_before.Lengths);
                    } else {
                        X0_global_before = null;
                        X0_global_mitte = null;
                    }


                    double LevelSetError = 0; // accumulator for error over all points.
                    double GradientError = 0; // accumulator for error over all points.
                    int outcnt = 0;
                    bool doResort = false; // flag which is set to true, if some point moves from one cell to another
                    int outsiders = 0;
                    OutsiderMarker.SetAll(false);

                    //int mark = -1;

                    // loop over all cells in the subgrid...
                    for(int jsub = 0; jsub < JSUB; jsub++) {
                        int I0 = I[jsub];
                        int IE = I[jsub + 1];
                        int K = IE - I0;
                        int jCell = jSub_2_jCell[jsub];
                        RefElement Kref = _GridData.Cells.GetRefElement(jCell);

                        if(K <= 0)
                            continue;

                        var X0_loc = this.X0_local.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });
                        var X0_glb = this.X0_global.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });

                        MultidimensionalArray x0_ip1_Global = MultidimensionalArray.Create(K, D);
                        var LevSetValues = MultidimensionalArray.Create(1, K);
                        var LevSetValues2 = MultidimensionalArray.Create(1, K);
                        MultidimensionalArray LevSetGrad = MultidimensionalArray.Create(1, K, D);
                        MultidimensionalArray LevSetGrad2 = MultidimensionalArray.Create(1, K, D);

                        LevSetValues.Clear();
                        LevSetGrad.Clear();

                        var old_X0_glb = X0_glb.CloneAs();

                        // suche in Gradienten-Richtung
                        // ----------------------------
                        #region gradientCorrection
                        {
                            NodeSet _X0_loc = new NodeSet(Kref, X0_loc);
                            LevSet.Evaluate(jCell, 1, _X0_loc, LevSetValues, 0, 0.0);
                            LevSet.EvaluateGradient(jCell, 1, _X0_loc, LevSetGrad);

                            // test convergence for each point

                            for(int nn = 0; nn < K; nn++) {

                                double errZeroSet = Math.Abs(LevSetValues[0, nn]);

                                double errGradient = 0, normGrad = 0, normVec = 0;
                                int jsub_ofX = this.CellIndexMap[I0 + nn];
                                int nn_ofX = this.NodeIndexMap[I0 + nn];
                                for(int d = 0; d < D; d++) {
                                    double V_d = X_global[jsub_ofX, nn_ofX, d] - X0_glb[nn, d];
                                    double N_d = LevSetGrad[0, nn, d];

                                    errGradient += V_d * N_d;
                                    normGrad += N_d * N_d;
                                    normVec += V_d * V_d;
                                }
                                errGradient = errGradient / Math.Sqrt(normVec * normGrad);

                                errGradient = Math.Abs(Math.Abs(errGradient) - 1);


                                this.X0_ErrorMeasure[I0 + nn] = errZeroSet + errGradient;
                                LevelSetError += errZeroSet;
                                GradientError += errGradient;

                                if(errZeroSet < ConvThreshold) {
                                    OutsideSgrdMarker[I0 + nn] = true;
                                }

                                //if(iter >= 40 && this.X0_ErrorMeasure[I0 + nn] > 0.02) {
                                //Console.WriteLine("Strange point at " + (I0 + nn));
                                //mark = I0 + nn;
                                //PrintPoint(X_global[jsub_ofX, nn_ofX, 0], X_global[jsub_ofX, nn_ofX, 1]);
                                //PrintPoint(X0_glb[nn, 0], X0_glb[nn, 1]);
                                //}
                            }

                            // move towards zero-levset
                            for(int nn = 0; nn < K; nn++) {

                                double sc = 0;
                                for(int d = 0; d < D; d++) {
                                    sc += LevSetGrad[0, nn, d].Pow2();
                                }


                                for(int d = 0; d < D; d++) {
                                    double xd = X0_glb[nn, d] - LevSetGrad[0, nn, d] * LevSetValues[0, nn] * (1 / sc);
                                    x0_ip1_Global[nn, d] = xd;
                                }

                                //if(I0+nn == mark)
                                //    PrintPoint(x0_ip1_Global[nn, 0], x0_ip1_Global[nn, 1]);

                                //radiusError += Math.Abs(LevSetValues[0, nn]);
                            }
                        }
                        # endregion

                        _GridData.TransformGlobal2Local(x0_ip1_Global, X0_loc, jCell, null); // transform to local
                        X0_glb.Set(x0_ip1_Global);                                            // record global 
                        
                        if(X0_global_mitte != null)
                            X0_global_mitte.ExtractSubArrayShallow(new int[] { I[jsub], 0 }, new int[] { I[jsub + 1] - 1, D - 1 })
                                .Set(x0_ip1_Global);
                        

                        // Tangential-Korrektur
                        // --------------------
                        #region tangentialCorrection
                        // Idee: der Vektor V := (x-x0) soll parallel zum Gradienten an der Stelle x0 sein (x0: closest point (output, iter), x: eingabe, const;
                        // sei N der normalenvektor; Korrektur:
                        // X0_neu = X0 - (I - (N^T*N))*V
                        {
                            NodeSet _X0_loc = new NodeSet(Kref, X0_loc);
                            LevSet.Evaluate(jCell, 1, _X0_loc, LevSetValues, 0, 0.0);
                            LevSet.EvaluateGradient(jCell, 1, _X0_loc, LevSetGrad);

                            for(int nn = 0; nn < K; nn++) {
                                double ooNorm = 0;
                                for(int d = 0; d < D; d++)
                                    ooNorm += LevSetGrad[0, nn, d].Pow2();
                                ooNorm = 1.0 / Math.Sqrt(ooNorm);

                                int jsub_ofX = this.CellIndexMap[I0 + nn];
                                int nn_ofX = this.NodeIndexMap[I0 + nn];

                                double NxV = 0;
                                for(int d = 0; d < D; d++) {
                                    double V_d = X_global[jsub_ofX, nn_ofX, d] - X0_glb[nn, d];
                                    double N_d = LevSetGrad[0, nn, d] * ooNorm;

                                    NxV += N_d * V_d;
                                }

                                //if(I0 + nn == mark)
                                //    PrintPoint(x0_ip1_Global[nn, 0] + LevSetGrad[0, nn, 0] * ooNorm*0.1, x0_ip1_Global[nn, 1] + LevSetGrad[0, nn, 1] * ooNorm*0.1);

                                for(int d = 0; d < D; d++) {
                                    double V_d = X_global[jsub_ofX, nn_ofX, d] - X0_glb[nn, d];
                                    double N_d = LevSetGrad[0, nn, d] * ooNorm;

                                    double verschiebung = (V_d - NxV * N_d);
                                    //verschiebung = 0;
                                    double X0_d = X0_glb[nn, d] + verschiebung;

                                    x0_ip1_Global[nn, d] = X0_d;

                                    //radiusError += Math.Abs(verschiebung);
                                }

                                //if(I0 + nn == mark)
                                //    PrintPoint(x0_ip1_Global[nn, 0], x0_ip1_Global[nn, 1]);

                            }


                            // adaptive under-relax:
                            // prevent the point from moving to far away from the zero-set under tangential correction
                            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




                            _GridData.TransformGlobal2Local(x0_ip1_Global, X0_loc, jCell, null); // transform to local
                            NodeSet _X0_loc2 = new NodeSet(Kref, X0_loc);
                            LevSet.Evaluate(jCell, 1, _X0_loc2, LevSetValues2, 0, 0.0);
                            LevSet.EvaluateGradient(jCell, 1, _X0_loc2, LevSetGrad2);


                            for(int nn = 0; nn < K; nn++) {
                                double newLevSetValue = LevSetValues2[0, nn];
                                double oldLevSetValue = LevSetValues[0, nn];


                                double maxGrowFactor = 1.0;

                                int jsub_ofX = this.CellIndexMap[I0 + nn];
                                int nn_ofX = this.NodeIndexMap[I0 + nn];


                                if(Math.Abs(newLevSetValue) >= maxGrowFactor * Math.Abs(oldLevSetValue)) {
                                    // tangential correction requires under-relax
                                    // ++++++++++++++++++++++++++++++++++++++++++



                                    //Console.WriteLine("Point {0}: levset value growed a lot under tangential correction: {1} to {2}", (I0 + nn), (LevSetValues[0, nn].Abs()), LevSetValues2[0, nn].Abs());

                                    //double limitLevSetValue = oldLevSetValue * maxGrowFactor;
                                    //
                                    //if(Math.Sign(newLevSetValue) != Math.Sign(oldLevSetValue))
                                    //    limitLevSetValue *= -1.0;
                                    //double alpha = (limitLevSetValue - newLevSetValue) / (oldLevSetValue - newLevSetValue);


                                    /*
                                    double lhs = 0;
                                    double rhs = 0;
                                    for(int d = 0; d < D; d++) {
                                        double xSd = x0_ip1_Global[nn, d];      // new value, including tangential correction
                                        double x1d = X0_glb[nn, d];             // old value, without tangential corredtion
                                        double yd = X_global[jsub_ofX, nn_ofX, d];  // original point
                                        double td = xSd - x1d;

                                        lhs += td * td;
                                        rhs -= x1d * td - yd * td;
                                    }

                                    double refAlpha;
                                    if(rhs.Abs() >= lhs.Abs())
                                        refAlpha = 1;
                                    else
                                        refAlpha = rhs / lhs;

                                    lhs += (-oldLevSetValue + newLevSetValue).Pow2();
                                    rhs -= oldLevSetValue * (-oldLevSetValue + newLevSetValue);


                                    double alpha;
                                    if(rhs.Abs() >= lhs.Abs())
                                        alpha = 1;
                                    else
                                        alpha = rhs / lhs;

                                    if(jCell == 583 && adaptive_iter > 0 && newLevSetValue.Abs() > 0.01)
                                        Console.WriteLine("alpha = " + alpha);
                                    alpha = Math.Max(0, Math.Min(1, alpha));
                                    */

                                    double facN0 = 0, facN1 = 0;
                                    for(int d = 0; d < D; d++) {
                                        facN0 += LevSetGrad[0, nn, d].Pow2();
                                        facN1 += LevSetGrad2[0, nn, d].Pow2();
                                    }
                                    facN0 = Math.Sqrt(1 / facN0);
                                    facN1 = Math.Sqrt(1 / facN1);


                                    


                                    double zähler = 0, nenner = 0;
                                    for(int d = 0; d < D; d++) {
                                        double y = X_global[jsub_ofX, nn_ofX, d]; // original point for which we are searching the closest point
                                        double x0 = X0_glb[nn, d];             // old value, without tangential corredtion
                                        double x1 = x0_ip1_Global[nn, d];      // new value, including tangential correction
                                        double n0 = LevSetGrad[0, nn, d] * facN0; // normal at old point
                                        double n1 = LevSetGrad2[0, nn, d] * facN1; // normal at new point

                                        zähler += n0 * (x0 - x1) + (n1 - n0) * (y - x0);
                                        nenner += (n0 - n1) * (x0 - x1);
                                    }
                                    nenner *= 2;
                                    double alpha;

                                    // Ansatz ist folgender:
                                    // X0, N0  Punkt und Normale VOR Tangentialkorrektur,
                                    // X1, N1 danach.
                                    //  X(alpha) = X0 + (X1-X0)*alpha;
                                    // Annahme: die Normale verhält in X(alpha) = N0 + (N1-N0)*alpha
                                    // dann soll q(alpha) := (Y - X(alpha))*N(alpha) minimiert werden.



                                    if(nenner.Abs() <= 1.0e-4 * zähler.Abs()) {
                                        alpha = 0.5;
                                    } else {
                                        alpha = zähler / nenner;
                                    }
                                    //if(alpha > 0 && alpha < 1)
                                    //    Console.WriteLine("öha! : " + alpha);
                                    //if(jCell == 583 && iter == 49 && newLevSetValue.Abs() > 0.01)
                                    //    Console.WriteLine("alpha = {0} = {1}/{2}", alpha, zähler, nenner);
                                    //alpha = 1;
                                    alpha = Math.Max(-0.0, Math.Min(0.7, alpha));


                                    //Debug.Assert(0.0 <= alpha && alpha <= 1.0);
                                    Debug.Assert(object.ReferenceEquals(x0_ip1_Global.Storage, X0_glb.Storage) == false);
                                    for(int d = 0; d < D; d++) {
                                        x0_ip1_Global[nn, d] = X0_glb[nn, d] * (1 - alpha) + x0_ip1_Global[nn, d] * alpha;
                                        //x0_ip1_Global[nn, d] = X0_glb[nn, d] * alpha + x0_ip1_Global[nn, d] * (1 - alpha);
                                    }
                                }
                            }

                        }
                        #endregion

                        _GridData.TransformGlobal2Local(x0_ip1_Global, X0_loc, jCell, null); // transform to local
                        X0_glb.Set(x0_ip1_Global);                                            // record global 


                        #region AltSearch_1
                        /*
                        if(iter == MaxIter -1) {

                            NodeSet _X0_loc = new NodeSet(Kref, X0_loc);
                            LevSet.Evaluate(jCell, 1, _X0_loc, LevSetValues, 0, 0.0);
                            LevSet.EvaluateGradient(jCell, 1, _X0_loc, LevSetGrad);

                            for(int nn = 0; nn < K; nn++) {
                                double minError = double.MaxValue;
                                double mindist = double.MaxValue;

                                for(int mm = 0; mm < K; mm++) {

                                    if(OutsideSgrdMarker[mm + I0])
                                        continue;

                                    double errZeroSet = Math.Abs(LevSetValues[0, mm]);

                                    double errGradient = 0, normGrad = 0, normVec = 0;
                                    int jsub_ofX = this.CellIndexMap[I0 + nn];
                                    int nn_ofX = this.NodeIndexMap[I0 + nn];
                                    for(int d = 0; d < D; d++) {
                                        double V_d = X_global[jsub_ofX, nn_ofX, d] - X0_glb[mm, d];
                                        double N_d = LevSetGrad[0, mm, d];

                                        errGradient += V_d * N_d;
                                        normGrad += N_d * N_d;
                                        normVec += V_d * V_d;
                                    }
                                    errGradient = errGradient / Math.Sqrt(normVec * normGrad);

                                    errGradient = Math.Abs(Math.Abs(errGradient) - 1);

                                    double _X0_ErrorMeasure = errZeroSet + errGradient;

                                    double _quasiDist = Math.Sqrt(normVec) + errZeroSet;


                                    //if(_X0_ErrorMeasure < minError) {
                                    if(_quasiDist < mindist) {
                                        for(int d = 0; d < D; d++) {
                                            x0_ip1_Global[nn, d] = X0_glb[mm, d];
                                        }
                                        minError = _X0_ErrorMeasure;
                                        mindist = _quasiDist;


                                        if(OutsideSgrdMarker[I0 + mm] == false && OutsideSgrdMarker[I0 + nn] == true)
                                            outsiders--;



                                        OutsideSgrdMarker[I0 + nn] = OutsideSgrdMarker[I0 + mm];
                                        //outsiders[I0 + nn] = outsiders[I0 + mm];
                                    }

                                }
                            }

                        }
                        _GridData.TransformGlobal2Local(x0_ip1_Global, X0_loc, jCell); // transform to local
                        X0_glb.Set(x0_ip1_Global);                                            // record global 
                        */
                        #endregion

                        //if(iter == MaxIter - 1) {
                        //    PlotOlt2NewStays("whatjump.csv", X0_global_before, X0_global_mitte);
                        //    PlotOlt2NewStays("whatjump2.csv", X0_global_mitte, X0_global);
                        //}
                        

                        // test if any nodes have left their cell, and find new cells
                        // ----------------------------------------------------------
                        #region nodesCellTest
                        {
                            int[] Neighs = null, dummy = null;
                            //double h_min_Neighs = 0.0;   // minimum cell diameter

                            for(int nn = 0; nn < K; nn++) { // loop over all nodes in sub-cell 'jsub'...

                                for(int d = 0; d < D; d++)
                                    _X0[d] = X0_loc[nn, d];

                                if(!Kref.IsWithin(_X0)) {

                                    // the iteration procedure jumped into another cell...
                                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++

                                    outcnt++;
                                    
                                    if(iter < Math.Max(5000, MaxIter / 4)) { // after some iterations, we do not change anymore,
                                        //              since some points may 'oscillate' between two cells

                                        doResort = true;  // it is necessary to sort the X0-vertices again.

                                        if(Neighs == null) {
                                            _GridData.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out Neighs, out dummy);
                                        }

                                        for(int d = 0; d < D; d++) {
                                            _X0_G[d] = X0_glb[nn, d];
                                            //_X_G[d] = old_X0_glb[nn, d];
                                        }

                                        int new_jCell = int.MinValue, new_jSub = int.MinValue;
                                        double dist = double.MaxValue;
                                        for(int nc = 0; nc < Neighs.Length; nc++) {
                                            int trial_jCell = Neighs[nc];

                                            int trial_jSub = jCell_to_jSub[trial_jCell];


                                            if(_GridData.Cells.IsInCell(_X0_G, trial_jCell, _X0_L)) {

                                                if(trial_jSub >= 0) {
                                                    // neighbour is in subgrid -> we are ok with that
                                                    new_jCell = trial_jCell;
                                                    new_jSub = trial_jSub;
                                                    dist = 0.0;
                                                    OutsideSgrdMarker[I0 + nn] = false;
                                                    break;
                                                } else {
                                                    // neighbour is NOT in subgrid -> we can't (resp. don't want to) assign to this cell
                                                    OutsideSgrdMarker[I0 + nn] = true;
                                                }
                                            }
                                        }

                                        if(dist != 0.0) {
                                            // contained in none of the neighbour cells: choose the closest one.

                                            outsiders++;

                                            for(int nc = -1; nc < Neighs.Length; nc++) {
                                                int trial_jCell;
                                                if(nc < 0)
                                                    trial_jCell = jCell;
                                                else
                                                    trial_jCell = Neighs[nc];

                                                int trial_jSub = jCell_to_jSub[trial_jCell];
                                                if(trial_jSub < 0)
                                                    // neighbour is not in subgrid -> continue
                                                    continue;

                                                double dist_to_trialCell = _GridData.Cells.ClosestPointInCell(_X0_G, trial_jCell, _X0_L2, dummyPt1, dummyPt2);

                                                if(dist_to_trialCell < dist) {
                                                    new_jCell = trial_jCell;
                                                    new_jSub = trial_jSub;
                                                    dist = dist_to_trialCell;
                                                    Array.Copy(_X0_L2, _X0_L, D);
                                                }
                                            }
                                            OutsiderMarker[I0 + nn] = true;
                                        }

                                        ActualjSub[I0 + nn] = new_jSub;
                                        for(int d = 0; d < D; d++)
                                            X0_loc[nn, d] = _X0_L[d];

                                        Debug.Assert(Kref.IsWithin(_X0_L) == _GridData.Cells.IsInCell(_X0_G, new_jCell));
                                        Debug.Assert(Kref.IsWithin(_X0_L) ^ OutsiderMarker[I0 + nn]);
                                        Debug.Assert(_GridData.Cells.IsInCell(_X0_G, new_jCell) ^ OutsiderMarker[I0 + nn]);
                                    } else {
                                        OutsiderMarker[I0 + nn] = true;
                                        outsiders++;
                                    }

                                } else {
                                    // the iteration procedure remained within the same cell...
                                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++

                                    OutsideSgrdMarker[I0 + nn] = false;
                                    OutsiderMarker[I0 + nn] = false;
                                    Debug.Assert(ActualjSub[I0 + nn] == jsub);
                                }
                            }
                        }
                        #endregion
                    }

                    // alternative search
                    // ------------------
                    
                    #region AltSearch_2
                    if(iter == MaxIter - 1) {

                        
                        //MultidimensionalArray LevSetGrad = MultidimensionalArray.Create(NNglobal, D);
                        MultidimensionalArray LevSetValues = MultidimensionalArray.Create(NNglobal);
                        
                        for(int jsub = 0; jsub < JSUB; jsub++) { // loop over all cells in the sub-grid
                            int I0 = I[jsub];
                            int IE = I[jsub + 1];
                            int K = IE - I0; // number of points in sub-cell
                            int jCell = jSub_2_jCell[jsub];
                            RefElement Kref = _GridData.Cells.GetRefElement(jCell);

                            if(K > 0) {

                                var X0_loc = this.X0_local.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });
                                var X0_glb = this.X0_global.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });
                                MultidimensionalArray x0_ip1_Global = MultidimensionalArray.Create(K, D);

                                var _LevSetValues = LevSetValues.ExtractSubArrayShallow(new int[] { I0 }, new int[] { IE - 1 }).ResizeShallow(1, K);
                                //var _LevSetGrad = LevSetGrad.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 }).ResizeShallow(1, K, D);

                                NodeSet _X0_loc = new NodeSet(Kref, X0_loc);
                                LevSet.Evaluate(jCell, 1, _X0_loc, _LevSetValues, 0, 0.0);
                                //LevSet.EvaluateGradient(jCell, 1, _X0_loc, _LevSetGrad);
                            }
                        }

                        MultidimensionalArray newX0_global = this.X0_global.CloneAs();
                        MultidimensionalArray newX0_local = this.X0_local.CloneAs();
                        int[] newActualjSub = ActualjSub.CloneAs();
                        BitArray newOutsiderMarker = OutsiderMarker.CloneAs();

                        for(int jsub = 0; jsub < JSUB; jsub++) {
                            int I0 = I[jsub];
                            int IE = I[jsub + 1];
                            int K = IE - I0;
                            int jCell = jSub_2_jCell[jsub];


                            int[] Neighs, dummy;
                            _GridData.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out Neighs, out dummy);
                            
                            for(int nn = 0; nn < K; nn++) {  // loop over all nodes in 'jCell' ...
                                double mindist = double.MaxValue;
                                int jsub_ofX = this.CellIndexMap[I0 + nn];
                                int nn_ofX = this.NodeIndexMap[I0 + nn];
                                int idxNewPoint = -1;
                                int jsubNewPoint = -1;
                                
                                //idxNewPoint = nn + I0;
                                //jsubNewPoint = jsub;

                                
                                for(int ii = -1; ii < Neighs.Length; ii++) {
                                    int jCellNeigh;
                                    if(ii < 0)
                                        jCellNeigh = jCell;
                                    else
                                        jCellNeigh = Neighs[ii];

                                    int jSubNeigh = jCell_to_jSub[jCellNeigh];
                                    if(jSubNeigh < 0)
                                        continue;

                                    int I0neigh = I[jSubNeigh];
                                    int IEneigh = I[jSubNeigh + 1];
                                    int K_neigh = IEneigh - I0neigh;

                                    for(int mm = 0; mm < K_neigh; mm++) {

                                        if(OutsideSgrdMarker[mm + I0neigh])
                                            continue;


                                        double errZeroSet = Math.Abs(LevSetValues[I0neigh + mm]);

                                        double normVec = 0;
                                        for(int d = 0; d < D; d++) {
                                            double V_d = X_global[jsub_ofX, nn_ofX, d] - X0_global[mm + I0neigh, d];
                                            normVec += V_d * V_d;
                                        }
                                        normVec = Math.Sqrt(normVec);

                                        double _quasiDist = normVec + errZeroSet;


                                        //if(_X0_ErrorMeasure < minError) {
                                        if(_quasiDist < mindist) {

                                            mindist = _quasiDist;
                                            idxNewPoint = mm + I0neigh;
                                            jsubNewPoint = jSubNeigh;
                                            //outsiders[I0 + nn] = outsiders[I0 + mm];
                                        }
                                    }
                                }
                                

                                if(idxNewPoint < 0) {
                                    throw new ApplicationException("Schade.");
                                } else {

                                    for(int d = 0; d < D; d++) {
                                        newX0_global[nn + I0, d] = X0_global[idxNewPoint, d];
                                        newX0_local[nn + I0, d] = this.X0_local[idxNewPoint, d];
                                    }
                                    //Debug.Assert(OutsideSgrdMarker[idxNewPoint] == false);
                                    
                                    newActualjSub[nn + I0] = ActualjSub[idxNewPoint];
                                    newOutsiderMarker[nn + I0] = OutsiderMarker[idxNewPoint];
                                }
                            }




                            //_GridData.TransformGlobal2Local(x0_ip1_Global, X0_loc, jCell); // transform to local
                            //X0_glb.Set(x0_ip1_Global);

                        }

                        //for(int i = 0; i < NNglobal; i++) {
                        //    bool jSubGleich = ActualjSub[i] == newActualjSub[i];
                        //    bool outmGleich = OutsiderMarker[i] == newOutsiderMarker[i];
                        //    bool pnktGleich = ArrayTools.AreEqual(X0_global.GetRow(i), newX0_global.GetRow(i));
                        //    Debug.Assert(jSubGleich);
                        //    Debug.Assert(outmGleich);
                        //    Debug.Assert(pnktGleich);
                        //}

                        //OutsideSgrdMarker.SetAll(false);
                        ActualjSub = newActualjSub;
                        X0_global = newX0_global;
                        X0_local = newX0_local;
                        OutsiderMarker = newOutsiderMarker;
                    }
                    #endregion
                   

                    // resorting
                    // ---------

                    #region resortingNodes
                    {
                        {
#if DEBUG
                            var _pt = new double[D];
                            var _ptG = new double[D];
                            for(int i = 0; i < NNglobal; i++) {
                                for(int d = 0; d < D; d++) {
                                    _pt[d] = this.X0_local[i, d];
                                    _ptG[d] = this.X0_global[i, d];
                                }

                                int jCell = jSub_2_jCell[ActualjSub[i]];
                                var splx = RefElms[_GridData.Cells.GetRefElementIndex(jCell)];

                                Debug.Assert(splx.IsWithin(_pt) ^ OutsiderMarker[i]);
                                Debug.Assert(_GridData.Cells.IsInCell(_ptG, jCell) ^ OutsiderMarker[i]);
                            }
#endif
                        }


                        // resort nodes, if necessary
                        if(doResort) {

                            // determine order
                            int[] NewIndex = new int[ActualjSub.Length];
                            for(int i =  NewIndex.Length - 1; i >= 0; i--)
                                NewIndex[i] = i;

                            Array.Sort(NewIndex, delegate(int a, int b) {
                                int diff;
                                diff = ActualjSub[a] - ActualjSub[b];

                                if(diff != 0)
                                    return diff;

                                diff = this.CellIndexMap[a] - this.CellIndexMap[b];

                                if(diff != 0)
                                    return diff;

                                diff = this.NodeIndexMap[a] - this.NodeIndexMap[b];

                                Debug.Assert(a == b || diff != 0);

                                return diff;
                            });

                            // sort data
                            int[] new_R_NodeIndexMap = new int[NNglobal];
                            int[] new_R_CellIndexMap = new int[NNglobal];
                            int[] new_ActualjSub = new int[NNglobal];
                            var new_X0_ErrorMeasure = MultidimensionalArray.Create(NNglobal);
                            var new_R_X0_local = MultidimensionalArray.Create(NNglobal, D);
                            var new_R_X0_global = MultidimensionalArray.Create(NNglobal, D);
                            var new_OutsideSgrdMarker = new System.Collections.BitArray(NNglobal);
                            var new_OutsiderMarker = new System.Collections.BitArray(NNglobal);

                            for(int hh = 0; hh < (NNglobal); hh++) {
                                int tt = NewIndex[hh];

                                new_R_NodeIndexMap[hh] = this.NodeIndexMap[tt];
                                new_R_CellIndexMap[hh] = this.CellIndexMap[tt];
                                new_ActualjSub[hh] = ActualjSub[tt];
                                new_OutsideSgrdMarker[hh] = OutsideSgrdMarker[tt];
                                new_X0_ErrorMeasure[hh] = this.X0_ErrorMeasure[tt];
                                new_OutsiderMarker[hh] = OutsiderMarker[tt];

                                int jCell = jSub_2_jCell[new_ActualjSub[hh]];

                                for(int d = 0; d < D; d++) {
                                    new_R_X0_local[hh, d] = this.X0_local[tt, d];
                                    new_R_X0_global[hh, d] = this.X0_global[tt, d];
                                }

                                //Debug.Assert(_Context.GridDat.IsInCell(R.X0_global.GetRow(tt), jCell));
                                //Debug.Assert(_Context.GridDat.IsInCell(new_R_X0_global.GetRow(hh), jCell));
                            }
                            this.X0_ErrorMeasure = new_X0_ErrorMeasure;
                            this.NodeIndexMap = new_R_NodeIndexMap;
                            this.CellIndexMap = new_R_CellIndexMap;
                            ActualjSub = new_ActualjSub;
                            this.X0_global = new_R_X0_global;
                            this.X0_local = new_R_X0_local;
                            OutsideSgrdMarker = new_OutsideSgrdMarker;
                            OutsiderMarker = new_OutsiderMarker;

                            // new array bla bla
                            I[0] = 0;
                            int currentCount = 0;
                            int p = 0;
                            for(int jSub = 0; jSub < JSUB; jSub++) {
                                while(p < new_ActualjSub.Length) {
                                    if(new_ActualjSub[p] <= jSub) {

                                        if(new_ActualjSub[p] == jSub)
                                            currentCount++;

                                        p++;

                                    } else {
                                        break;
                                    }
                                }

                                I[jSub + 1] = I[jSub] + currentCount;
                                currentCount = 0;
                            }
                            Debug.Assert(I[JSUB] == NNglobal);

                            // test
#if DEBUG
                            for(int i = 1; i < NNglobal; i++) {
                                Debug.Assert(ActualjSub[i - 1] <= ActualjSub[i]);
                            }

                            for(int jsub = 0; jsub < JSUB; jsub++) {
                                int I0 = I[jsub];
                                int IE = I[jsub + 1];
                                int K = IE - I0;

                                for(int nn = 0; nn < K; nn++) {
                                    Debug.Assert(ActualjSub[I0 + nn] == jsub);
                                }
                            }

#endif
                            // 
                        }

                        {
#if DEBUG
                            for(int jsub = 0; jsub < JSUB; jsub++) {

                                var pt = new double[D];
                                var ptG = new double[D];

                                int I0 = I[jsub];
                                int IE = I[jsub + 1];
                                int K = IE - I0;
                                int jCell = jSub_2_jCell[jsub];
                                var splx = RefElms[_GridData.Cells.GetRefElementIndex(jCell)];

                                if(K <= 0)
                                    continue;

                                var X0_loc = this.X0_local.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });
                                var X0_glb = this.X0_global.ExtractSubArrayShallow(new int[] { I0, 0 }, new int[] { IE - 1, D - 1 });

                                for(int nn = 0; nn < K; nn++) {
                                    for(int d = 0; d < D; d++) {
                                        pt[d] = X0_loc[nn, d];
                                        ptG[d] = X0_glb[nn, d];
                                    }

                                    Debug.Assert(splx.IsWithin(pt) ^ OutsiderMarker[I0 + nn]);
                                    Debug.Assert(_GridData.Cells.IsInCell(ptG, jCell) ^ OutsiderMarker[I0 + nn]);
                                }

                            }
#endif
                        }
                    }
                    #endregion

                    // alternative search
                    // ------------------
                    #region altSearc_3
                    /*
                    if(iter == MaxIter - 2) {

                        System.Collections.BitArray new_OutsideSgrdMarker = OutsideSgrdMarker.CloneAs();


                        for(int jsub = 0; jsub < JSUB; jsub++) {
                            int I0 = I[jsub];
                            int IE = I[jsub + 1];
                            int K = IE - I0;
                            int jCell = jSub_2_jCell[jsub];
                            RefElement Kref = _GridData.Cells.GetRefElement(jCell);

                            var X0_loc = this.X0_local.ExtractSubArrayShallow(new int[] { I[jsub], 0 }, new int[] { I[jsub + 1] - 1, D - 1 });
                            var X0_glb = this.X0_global.ExtractSubArrayShallow(new int[] { I[jsub], 0 }, new int[] { I[jsub + 1] - 1, D - 1 });
                            MultidimensionalArray x0_ip1_Global = MultidimensionalArray.Create(K, D);
                            MultidimensionalArray LevSetGrad = MultidimensionalArray.Create(1, K, D);
                            var LevSetValues = MultidimensionalArray.Create(1, K);

                            NodeSet _X0_loc = new NodeSet(Kref, X0_loc);
                            LevSet.Evaluate(jCell, 1, _X0_loc, LevSetValues, 0, 0.0);
                            LevSet.EvaluateGradient(jCell, 1, _X0_loc, LevSetGrad);

                            
                            //{
                            //    bool zumindestEiner = false;
                            //    for(int mm = 0; mm < K; mm++) {

                            //        if(OutsideSgrdMarker[mm + I0])
                            //            continue;
                            //        zumindestEiner = true;
                            //    }

                            //    if(!zumindestEiner)
                            //        Console.WriteLine("    Zelle ohne Punkt: " + jCell);
                            //}

                            int[] Neighs, dummy;
                            _GridData.Cells.GetCellNeighbours(jCell, GridData.CellData.GetCellNeighbours_Mode.ViaVertices, out Neighs, out dummy);
                            





                            for(int nn = 0; nn < K; nn++) {
                                double minError = double.MaxValue;
                                double mindist = double.MaxValue;

                                for(int mm = 0; mm < K; mm++) {

                                    if(OutsideSgrdMarker[mm + I0])
                                        continue;
                                    double errZeroSet = Math.Abs(LevSetValues[0, mm]);

                                    double errGradient = 0, normGrad = 0, normVec = 0;
                                    int jsub_ofX = this.CellIndexMap[I0 + nn];
                                    int nn_ofX = this.NodeIndexMap[I0 + nn];
                                    for(int d = 0; d < D; d++) {
                                        double V_d = X_global[jsub_ofX, nn_ofX, d] - X0_glb[mm, d];
                                        double N_d = LevSetGrad[0, mm, d];

                                        errGradient += V_d * N_d;
                                        normGrad += N_d * N_d;
                                        normVec += V_d * V_d;
                                    }
                                    errGradient = errGradient / Math.Sqrt(normVec * normGrad);

                                    errGradient = Math.Abs(Math.Abs(errGradient) - 1);

                                    double _X0_ErrorMeasure = errZeroSet + errGradient;

                                    double _quasiDist = Math.Sqrt(normVec) + errZeroSet;


                                    //if(_X0_ErrorMeasure < minError) {
                                    if(_quasiDist < mindist) {
                                        for(int d = 0; d < D; d++) {
                                            x0_ip1_Global[nn, d] = X0_glb[mm, d];
                                        }
                                        minError = _X0_ErrorMeasure;
                                        mindist = _quasiDist;


                                        if(OutsideSgrdMarker[I0 + mm] == false && OutsideSgrdMarker[I0 + nn] == true)
                                            outsiders--;



                                        new_OutsideSgrdMarker[I0 + nn] = OutsideSgrdMarker[I0 + mm];
                                        //outsiders[I0 + nn] = outsiders[I0 + mm];
                                    }

                                }


                            }

                            


                            _GridData.TransformGlobal2Local(x0_ip1_Global, X0_loc, jCell); // transform to local
                            X0_glb.Set(x0_ip1_Global);

                        }

                        OutsideSgrdMarker = new_OutsideSgrdMarker;
                    } //*/
                    #endregion

                    // termination criterion
                    // =====================

                    double scali = (NNglobal);
                    if(prevRadiusError <= ConvThreshold * scali && LevelSetError <= ConvThreshold * scali && outsiders <= 0 && iter >= MinIter) {
                        //Console.WriteLine("Convergence criterion reached.");
                        break;
                    }
                    prevRadiusError = LevelSetError;
                    if (DiagnosticOutput) {
                        Console.WriteLine("iter #" + iter + ", LevelSetError = " + LevelSetError.ToStringDot() + ", GradientError = " + GradientError + ", outcnt = " + outcnt + ", not assignable to Neighbour: " + outsiders + " (of " + this.X0_local.GetLength(0) + ")");
                    }
                }
                if (DiagnosticOutput) {
                    PlottiPunkti(this.X0_global, "X0_iter-" + MaxIter + ".csv", this.X0_ErrorMeasure, OutsideSgrdMarker);
                }
                //PlotPairings("outOfSubgrid.csv", OutsideSgrdMarker);
                //PlotPairings("OutOfCell.csv", OutsiderMarker);
                //PlotCellAssociation("CellAssoc.csv", OutsiderMarker);

                int noPts = this.X0_global.GetLength(0);
                var ba = new System.Collections.BitArray(noPts);
                ba.SetAll(true);
                if (DiagnosticOutput) {
                    PlotStays("HR.csv", ba);   //GenerateHighResiMarker());

                    PlotCellAssociation("CellAssoc_HR.csv", GenerateHighResiMarker());
                }
            }
        }

        private System.Collections.BitArray GenerateHighResiMarker() {
            System.Collections.BitArray HighResi;
            HighResi = new System.Collections.BitArray(this.X0_ErrorMeasure.GetLength(0));
            for(int i = 0; i < HighResi.Length; i++) {
                HighResi[i] = this.X0_ErrorMeasure[i] > 0.0005;
            }
            return HighResi;
        }

        
        private static bool LocatePoint(GridData _Context, int[] jCell_to_jSub, int[] ActualjSub, int I0, MultidimensionalArray X0_loc, double[] _X0_G, double[] _X0_L, int[] Neighs, int nn) {
            bool found = false;
            int D = _Context.SpatialDimension;

            foreach (int jNeigh in Neighs) {
                if (_Context.Cells.IsInCell(_X0_G, jNeigh, _X0_L)) {
                    //newAss++;
                    found = true;

                    int new_jSub = jCell_to_jSub[jNeigh];
                    if (new_jSub < 0)
                        //Console.WriteLine("WARN");
                        throw new Exception("");

                    ActualjSub[I0 + nn] = new_jSub;

                    for (int d = 0; d < D; d++)
                        X0_loc[nn, d] = _X0_L[d];

                    break;
                }
            }

            return found;
        }
        // */


        public void WriteX0ToFile(string filename) {
            using (var stw = new StreamWriter(filename)) {

                int L = this.X0_global.GetLength(0);
                int D = this.X0_global.GetLength(1);

                GridData grdDat = (GridData)(sgrd.GridData);

                int JSUB = sgrd.LocalNoOfCells;
                int[] jSub_2_jCell = sgrd.SubgridIndex2LocalCellIndex;

                for (int jsub = 0; jsub < JSUB; jsub++) {
                    int I0 = this.x0I2SgrdCell[jsub];
                    int IE = this.x0I2SgrdCell[jsub + 1];
                    int K = IE - I0;
                    int jCell = jSub_2_jCell[jsub];

                    for (int l = I0; l < IE; l++) {
                        stw.Write(l);
                        stw.Write(" ");

                        // subgrid of x0
                        stw.Write(jsub);
                        stw.Write(" ");

                        // cell index of x0
                        stw.Write(jCell);
                        stw.Write(" ");
                        
                        // inside the assigned cell ? 
                        var X0 = X0_global.GetRow(l);
                        bool inside = grdDat.Cells.IsInCell(X0, jCell);
                        stw.Write((inside ? 1.0 : 0.0));
                        stw.Write(" ");
                        //stw.Write(inside);
                        //stw.Write(" ");
                        //if (!inside)
                        //    Debugger.Break();
                        
                        for (int d = 0; d < D; d++) {
                            stw.Write(X0[d].ToStringDot());
                            stw.Write(" ");
                        }
                        
                        stw.WriteLine();
                    }
                }
            }
        }
    }
}
