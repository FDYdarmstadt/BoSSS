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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Connectors.Matlab;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.Grid;

namespace BoSSS.Application.SipPoisson {

    /// <summary>
    /// predefined control-objects
    /// </summary>
    static public class SipHardcodedControl {

        /// <summary>
        /// Test on a curved grid.
        /// </summary>
        public static SipControl TestCurved() {
            var R = new SipControl();
            R.ProjectName = "ipPoison/curved";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = 2, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 15 });
            R.InitialValues_Evaluators.Add("RHS", X => 0.0);
            R.InitialValues_Evaluators.Add("Tex", X => (Math.Log(X[0].Pow2() + X[1].Pow2()) / Math.Log(4.0)) + 1.0);
            R.ExactSolution_provided = true;

            R.GridFunc = delegate() {
                var grd = Grid2D.CurvedSquareGrid(GenericBlas.Linspace(1, 2, 3), GenericBlas.Linspace(0, 1, 11), CellType.Square_9, true);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(X => 1);
                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     double x = X[0], y = X[1];
                     return Math.Sqrt(x * x + y * y);
                 });


            return R;
        }

        /// <summary>
        /// Creates Nodes, yes it really does!
        /// </summary>
        /// <param name="res"></param>
        /// <param name="stetch">
        /// Factor which determines how much the intervals in the output grow; 1.0 is no stretching.
        /// </param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        static double[] CreateNodes(int res, double stetch, double min, double max) {
            if(stetch == 1.0)
                return GenericBlas.Linspace(min, max, res + 1);
            else
                return Grid1D.ExponentialSpaceing(min, max, res + 1, stetch); // without proper preconditioning,
            // a stretched grid is much more expensive than
            // an equidistant grid !!!
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// </summary>
        public static SipControl TestCartesian1(int xRes = 32, double xStretch = 1.0, int yRes = 16, double yStretch = 1.01, int pDG = 2) {
            var RR = new SipControl();
            RR.ProjectName = "ipPoison/cartesian";
            RR.savetodb = false;

            RR.FieldOptions.Add("T", new FieldOpts() { Degree = pDG, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            RR.FieldOptions.Add("Tex", new FieldOpts() { Degree = pDG * 2 });
            RR.InitialValues_Evaluators.Add("RHS", X => 1.0);
            RR.InitialValues_Evaluators.Add("Tex", X => (0.5 * X[0].Pow2() - 10 * X[0]));
            RR.ExactSolution_provided = true;

            RR.GridFunc = delegate() {
                double[] xNodes = CreateNodes(xRes, xStretch, 0, 10);
                double[] yNodes = CreateNodes(yRes, yStretch, -1, +1);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });
                
                return grd;
            };


            RR.AddBoundaryValue(BoundaryType.Dirichlet.ToString());
            RR.AddBoundaryValue(BoundaryType.Neumann.ToString());


            RR.GridPartType = BoSSS.Foundation.Grid.GridPartType.none;


            return RR;
        }

        /// <summary>
        /// Test on a Cartesian grid, with an exact polynomial solution.
        /// </summary>
        public static SipControl TestCartesian3D(int xRes = 32, double xStretch = 1.0, int yRes = 16, double yStretch = 1.0, int zRes = 16, double zStretch = 1.0) {
            var R = new SipControl();
            R.ProjectName = "ipPoison/cartesian";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = 6, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 6 });
            R.InitialValues_Evaluators.Add("RHS", X => 1.0);
            R.InitialValues_Evaluators.Add("Tex", X => (0.5 * X[0].Pow2() - 10 * X[0]));
            R.ExactSolution_provided = true;

            R.GridFunc = delegate() {
                double[] xNodes = CreateNodes(xRes, xStretch, 0, 10);
                double[] yNodes = CreateNodes(yRes, yStretch, -1, +1);
                double[] zNodes = CreateNodes(zRes, zStretch, -1, +1);

                var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };
            
            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString());
            R.AddBoundaryValue(BoundaryType.Neumann.ToString());


            return R;
        }


     


        /// <summary>
        /// Test on a Cartesian grid, with a sinusodial solution.
        /// </summary>
        /// <param name="Res">
        /// Grid resolution
        /// </param>
        /// <param name="Dim">
        /// spatial dimension
        /// </param>
        /// <param name="deg">
        /// polynomial degree
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        public static SipControl TestCartesian2(int Res, int Dim, SolverCodes solver_name = SolverCodes.exp_softpcg_mg, int deg = 3) {
            if(Dim != 2 && Dim != 3)
                throw new ArgumentOutOfRangeException();
            
            var R = new SipControl();
            R.ProjectName = "ipPoison/cartesian";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg*2 });
            R.InitialValues_Evaluators.Add("RHS", X => -Math.Sin(X[0]));
            R.InitialValues_Evaluators.Add("Tex", X => Math.Sin(X[0]));
            R.ExactSolution_provided = true;
            R.NoOfMultigridLevels = int.MaxValue;
            R.solver_name = solver_name;
            //R.TargetBlockSize = 100;

            R.TracingNamespaces = "BoSSS,ilPSP";
            

            R.GridFunc = delegate() {
                GridCommons grd = null;
                if(Dim == 2) {
                    double[] xNodes = GenericBlas.Linspace(0, 10, Res*5 + 1);
                    double[] yNodes = GenericBlas.SinLinSpacing(-1, +1, 0.6, Res + 1);
                    
                    grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                } else if(Dim == 3) {
                    double[] xNodes = GenericBlas.Linspace(0, 10, Res*5 + 1);
                    double[] yNodes = GenericBlas.SinLinSpacing(-1, +1, 0.6, Res + 1);
                    double[] zNodes = GenericBlas.SinLinSpacing(-1, +1, 0.6, Res + 1);

                    grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                } else {
                    throw new NotSupportedException();
                }
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.EdgeTagNames.Add(2, BoundaryType.Neumann.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-6)
                        ret = 1;
                    else
                        ret = 2;
                    return ret;
                });

                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     double x = X[0], y = X[1];

                     if(Math.Abs(X[0] - (0.0)) < 1.0e-8)
                         return 0.0;

                     throw new ArgumentOutOfRangeException();
                 });

            R.AddBoundaryValue(BoundaryType.Neumann.ToString(), "T",
                 delegate (double[] X) {
                     if(Math.Abs(X[1] - 1.0) < 1.0e-8 || Math.Abs(X[1] + 1.0) < 1.0e-8) // y = -1, y = +1
                         return 0;

                     if(X.Length > 2 && (Math.Abs(X[2] - 1.0) < 1.0e-8 || Math.Abs(X[2] + 1.0) < 1.0e-8)) // z = -1, z = +1
                         return 0;

                     if(Math.Abs(X[0] - (+10.0)) < 1.0e-8)
                         return Math.Cos(10.0);

                     throw new ArgumentOutOfRangeException();
                 });



           
            return R;
        }


        /// <summary>
        /// Poisson Equation on a (-1,1)x(-1,1), Dirichlet everywhere
        /// </summary>
        public static SipControl Square(int xRes = 5, int yRes = 5, int deg = 5) {

            //Func<double[], double> exRhs = X => 2 * X[0] * X[0] + 2 * X[1] * X[1] - 4;
            //Func<double[], double> exSol = X => (1.0 - X[0] * X[0]) * (1.0 - X[1] * X[1]);

            //Func<double[], double> exSol = X => (1.0 - X[1]);
            //Func<double[], double> exRhs = X => 0.0;

            Func<double[], double> exSol = X => -Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5);
            Func<double[], double> exRhs = X => (Math.PI * Math.PI * 0.5 * Math.Cos(X[0] * Math.PI * 0.5) * Math.Cos(X[1] * Math.PI * 0.5)); // == - /\ exSol


            var R = new SipControl();
            R.ProjectName = "ipPoison/square";
            R.savetodb = false;
            //R.DbPath = "D:\\BoSSS-db";

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = 4 });
            R.InitialValues_Evaluators.Add("RHS", exRhs);
            R.InitialValues_Evaluators.Add("Tex", exSol);
            R.ExactSolution_provided = true;

            R.GridFunc = delegate() {
                double[] xNodes = GenericBlas.Linspace(-1,1,xRes);
                double[] yNodes = GenericBlas.Linspace(-1,1,yRes);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte ret = 1;
                    return ret;
                });


                return grd;
            };

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T", exSol);

            R.NoOfSolverRuns = 1;

            R.AdaptiveMeshRefinement = true;
            R.NoOfTimesteps = 5; 


            return R;
        }

        /*

        /// <summary>
        /// Test on a 2D Voronoi mesh
        /// </summary>
        /// <param name="Res">
        /// number of randomly chosen Delaunay vertices
        /// </param>
        /// <param name="deg">
        /// polynomial degree
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        public static SipControl TestVoronoiOld(int Res, SolverCodes solver_name = SolverCodes.classic_pardiso, int deg = 3) {
            
             if (System.Environment.MachineName.ToLowerInvariant().EndsWith("rennmaschin")
                //|| System.Environment.MachineName.ToLowerInvariant().Contains("jenkins")
                ) {
                // This is Florians Laptop;
                // he is to poor to afford MATLAB, so he uses OCTAVE
                BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;
                BatchmodeConnector.MatlabExecuteable = "C:\\cygwin64\\bin\\bash.exe";
            } 

            
            var R = new SipControl();
            R.ProjectName = "SipPoisson-Voronoi";
            R.SessionName = "testrun";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg*2 });
            R.InitialValues_Evaluators.Add("RHS", X => 1.0);
            R.InitialValues_Evaluators.Add("Tex", X => 0.0);
            R.ExactSolution_provided = false;
            R.NoOfMultigridLevels = int.MaxValue;
            R.solver_name = solver_name;
            //R.TargetBlockSize = 100;


            

            bool IsIn(double xi, double yi) {
                
                //for(int l = 0; l < bndys.Length; l++) {
                //    Debug.Assert(bndys[l].Normal.Length == 2);
                //    if (bndys[l].PointDistance(xi, yi) > 0.0)
                //        return false;
                //}
                if (xi > 1.0)
                    return false;
                if (yi > 1.0)
                    return false;
                if (xi < 0 && yi < 0)
                    return false;
                if (xi < -1)
                    return false;
                if (yi < -1)
                    return false;

                return true;
            }

            bool IsInV(Vector X) {
                Debug.Assert(X.Dim == 2);
                return IsIn(X.x, X.y);
            }

            int Mirror(ref double[] _x, ref double[] _y, AffineManifold[] bndys) {
                if (_x.Length != _y.Length)
                    throw new ArgumentException();
                var x = _x.ToList();
                var y = _y.ToList();
                int N = _x.Length;



                // filter all points that are outside of the domain
                for (int n = 0; n < N; n++) {
                    if (!IsIn(x[n], y[n])) {
                        x.RemoveAt(n);
                        y.RemoveAt(n);
                        N--;
                        n--;
                    }
                }
                Debug.Assert(x.Count == N);
                Debug.Assert(y.Count == N);
                for (int n = 0; n < N; n++) {
                    Debug.Assert(IsIn(x[n], y[n]));
                }

                // mirror each point 
                for(int n = 0; n < N; n++) {
                    double xn = x[n];
                    double yn = y[n];
                    for(int l = 0; l < bndys.Length; l++) {
                        var bndy_l = bndys[l];

                        double dist = bndy_l.PointDistance(xn, yn);
                        
                        if(dist < 0) {
                            double xMirr = xn - bndy_l.Normal[0] * dist*2;
                            double yMirr = yn - bndy_l.Normal[1] * dist*2;

                            Debug.Assert(bndy_l.PointDistance(xMirr, yMirr) > 0);

                            if(!IsIn(xMirr, yMirr)) {
                                x.Add(xMirr);
                                y.Add(yMirr);
                            }
                        }
                    }
                }

                // return
                _x = x.ToArray();
                _y = y.ToArray();
                return N;
            }

            
            IGrid GridFunc() {
                GridCommons grd = null;

                var Matlab = new BatchmodeConnector();

                // boundaries for L-domain
                AffineManifold[] Boundaries = new AffineManifold[6];
                Boundaries[0] = new AffineManifold(new[] { 0.0, 1.0 }, new[] { 0.0, 1.0 });
                Boundaries[1] = new AffineManifold(new[] { 1.0, 0.0 }, new[] { 1.0, 0.0 });
                Boundaries[2] = new AffineManifold(new[] { -1.0, 0.0 }, new[] { -1.0, 0.0 });
                Boundaries[3] = new AffineManifold(new[] { -1.0, 0.0 }, new[] { 0.0, 0.0 });
                Boundaries[4] = new AffineManifold(new[] { 0.0, -1.0 }, new[] { 0.0, -1.0 });
                Boundaries[5] = new AffineManifold(new[] { 0.0, -1.0 }, new[] { 0.0, 0.0 });

                
                // generate Delaunay vertices
                Random rnd = new Random(0);
                double[] xNodes = Res.ForLoop(idx => rnd.NextDouble()*2 - 1);
                double[] yNodes = Res.ForLoop(idx => rnd.NextDouble()*2 - 1);
                int ResFix = Mirror(ref xNodes, ref yNodes, Boundaries);

                
                var Nodes = MultidimensionalArray.Create(xNodes.Length, 2);
                Nodes.SetColumn(0, xNodes);
                Nodes.SetColumn(1, yNodes);
                
                Matlab.PutMatrix(Nodes, "Nodes");
               
                // compute Voronoi diagramm
                Matlab.Cmd("[V, C] = voronoin(Nodes);");

                // output (export from matlab)
                int[][] OutputVertexIndex = new int[Nodes.NoOfRows][];
                Matlab.GetStaggeredIntArray(OutputVertexIndex, "C");
                Matlab.GetMatrix(null, "V");

                 // run matlab
                Matlab.Execute(false);

                // import here
                MultidimensionalArray VertexCoordinates = (MultidimensionalArray)(Matlab.OutputObjects["V"]);

                // correct indices (1-based index to 0-based index)
                foreach(int[] cell in OutputVertexIndex) {
                    int K = cell.Length;
                    for (int k = 0; k < K; k++) {
                        cell[k]--;
                    }
                }

                // tessellation
                List<Cell> cells = new List<Cell>();
                List<int[]> aggregation = new List<int[]>();
                for(int jV = 0; jV < ResFix; jV++) { // loop over Voronoi Cells
                    Debug.Assert(IsInV(Nodes.GetRowPt(jV)));

                    int[] iVtxS = OutputVertexIndex[jV];
                    int NV = iVtxS.Length;

                    List<int> Agg2Pt = new List<int>();

                    for(int iTri = 0; iTri < NV - 2; iTri++) { // loop over triangles of voronoi cell
                        int iV0 = iVtxS[0];
                        int iV1 = iVtxS[iTri + 1];
                        int iV2 = iVtxS[iTri + 2];

                        Vector V0 = VertexCoordinates.GetRowPt(iV0);
                        Vector V1 = VertexCoordinates.GetRowPt(iV1);
                        Vector V2 = VertexCoordinates.GetRowPt(iV2);

                        double[] D1 = V1 - V0;
                        double[] D2 = V2 - V0;
                        Debug.Assert(D1.CrossProduct2D(D2).Abs() > 1.0e-8);
                        if(D1.CrossProduct2D(D2) < 0) {
                            Vector T = V2;
                            int t = iV2;
                            V2 = V1;
                            iV2 = iV1;
                            V1 = T;
                            iV1 = t;
                        }

                        //double[] Center = V0.Plus(V1).Plus(V2).Mul(1.0 / 3.0);
                        //Debug.Assert(IsIn(Center[0], Center[1]));

                        Cell Cj = new Cell();
                        Cj.GlobalID = cells.Count;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                        Cj.NodeIndices = new int[] { iV0, iV1, iV2 };
                        Cj.TransformationParams.SetRowPt(0, V0);
                        Cj.TransformationParams.SetRowPt(1, V1);
                        Cj.TransformationParams.SetRowPt(2, V2);

                        Agg2Pt.Add(cells.Count);

                        cells.Add(Cj);
                    }

                    aggregation.Add(Agg2Pt.ToArray());
                }

                // return grid
                grd = new Grid2D(Triangle.Instance);
                grd.Cells = cells.ToArray();
                grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
                grd.DefineEdgeTags(X => (byte)1);

                //grd.Plot2DGrid();

                // create aggregation grid
                var agrd = new AggregationGrid(grd, aggregation.ToArray());
                return agrd;

            };
            R.GridFunc = GridFunc;

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     //double x = X[0], y = X[1];

                     return 0.0;
                     //if(Math.Abs(X[0] - (0.0)) < 1.0e-8)
                     //    return 0.0;
                     //
                     //throw new ArgumentOutOfRangeException();
                 });

            


           
            return R;
        }

        //*/

        /// <summary>
        /// Test on a 2D Voronoi mesh
        /// </summary>
        /// <param name="Res">
        /// number of randomly chosen Delaunay vertices
        /// </param>
        /// <param name="deg">
        /// polynomial degree
        /// </param>
        /// <param name="solver_name">
        /// Name of solver to use.
        /// </param>
        public static SipControl TestVoronoi(int Res, SolverCodes solver_name = SolverCodes.classic_pardiso, int deg = 3) {
            
             if (System.Environment.MachineName.ToLowerInvariant().EndsWith("rennmaschin")
                //|| System.Environment.MachineName.ToLowerInvariant().Contains("jenkins")
                ) {
                // This is Florians Laptop;
                // he is to poor to afford MATLAB, so he uses OCTAVE
                BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;
                BatchmodeConnector.MatlabExecuteable = "C:\\cygwin64\\bin\\bash.exe";
            }



            var R = new SipControl();
            R.ProjectName = "SipPoisson-Voronoi";
            R.SessionName = "testrun";
            R.savetodb = false;

            R.FieldOptions.Add("T", new FieldOpts() { Degree = deg, SaveToDB = FieldOpts.SaveToDBOpt.TRUE });
            R.FieldOptions.Add("Tex", new FieldOpts() { Degree = deg*2 });
            R.InitialValues_Evaluators.Add("RHS", X => -1.0 + X[0]*X[0] );
            R.InitialValues_Evaluators.Add("Tex", X => 0.0);
            R.ExactSolution_provided = false;
            R.NoOfMultigridLevels = int.MaxValue;
            R.solver_name = solver_name;
            R.NoOfMultigridLevels = 0;
            //R.TargetBlockSize = 100;


            
            
            bool IsIn(double xi, double yi) {
                
                //for(int l = 0; l < bndys.Length; l++) {
                //    Debug.Assert(bndys[l].Normal.Length == 2);
                //    if (bndys[l].PointDistance(xi, yi) > 0.0)
                //        return false;
                //}
                if (xi > 1.0)
                    return false;
                if (yi > 1.0)
                    return false;
                if (xi < 0 && yi < 0)
                    return false;
                if (xi < -1)
                    return false;
                if (yi < -1)
                    return false;

                return true;
            }

           

            Vector[] DomainBndyPolygon = new[] {
                new Vector(+0,+0),
                new Vector(-1,+0),
                new Vector(-1,+1),
                new Vector(+1,+1),
                new Vector(+1,-1),
                new Vector(+0,-1)
            };
            
            /*
            
            bool IsIn(double xi, double yi) {
                double myEps = 0.0;
                if (xi > 1.0 + myEps)
                    return false;
                if (yi > 1.0 + myEps)
                    return false;
                if (xi < -1 - myEps)
                    return false;
                if (yi < -1 - myEps)
                    return false;

                return true;
            }

           

            Vector[] DomainBndyPolygon = new[] {
                new Vector(-1,+1),
                new Vector(+1,+1),
                new Vector(+1,-1),
                new Vector(-1,-1)
            };
            //*/

            bool IsInV(Vector X) {
                Debug.Assert(X.Dim == 2);
                return IsIn(X.x, X.y);
            }

            bool Idenity(Vector A, Vector B) {
                bool t2 = (A - B).AbsSquare() < 1.0e-10;
                bool t1 = (A - B).AbsSquare() < 1.0e-15;
                Debug.Assert(t1 == t2);
                return t1;
            }
            
            IGrid GridFunc() {
                                
                // generate Delaunay vertices
                Random rnd = new Random(0);
                int RR = Res ;
                var Node = MultidimensionalArray.Create(RR, 2);

                bool useMirror = false;
                double scl = useMirror ? 2.0 : 4.0;

                Node.SetColumn(0, RR.ForLoop(idx => rnd.NextDouble() * scl - 0.5*scl));
                Node.SetColumn(1, RR.ForLoop(idx => rnd.NextDouble() * scl - 0.5*scl));

                // generate mesh
                return Voronoi.VoronoiMeshGen.FromPolygonalDomain(Node, DomainBndyPolygon, useMirror, IsInV, Idenity);

            };
            R.GridFunc = GridFunc;

            R.AddBoundaryValue(BoundaryType.Dirichlet.ToString(), "T",
                 delegate (double[] X) {
                     //double x = X[0], y = X[1];

                     return 0.0;
                     //if(Math.Abs(X[0] - (0.0)) < 1.0e-8)
                     //    return 0.0;
                     //
                     //throw new ArgumentOutOfRangeException();
                 });

            


           
            return R;
        }

        /*
        static void VoronoiTestCode() { 
            Vector[] DomainBndy = new[] {
                new Vector(-1, 0), // 6
                new Vector(-1, 1), // 5
                new Vector(1, 1), // 4
                new Vector(1, -1), // 3
                new Vector(0, -1), // 2
                new Vector(0, 0), // 1
            };

            bool IsIn(double xi, double yi) {

                //for(int l = 0; l < bndys.Length; l++) {
                //    Debug.Assert(bndys[l].Normal.Length == 2);
                //    if (bndys[l].PointDistance(xi, yi) > 0.0)
                //        return false;
                //}
                if (xi > 1.0)
                    return false;
                if (yi > 1.0)
                    return false;
                if (xi < 0 && yi < 0)
                    return false;
                if (xi < -1)
                    return false;
                if (yi < -1)
                    return false;

                return true;
            }

            bool IsInV(Vector X) {
                Debug.Assert(X.Dim == 2);
                return IsIn(X.x, X.y);
            }

            Gnuplot gp = new Gnuplot();

            gp.PlotXY(DomainBndy.Select(X => X.x).ToArray(), DomainBndy.Select(X => X.y).ToArray(), "domain",
                new PlotFormat("-xk"));

            Vector[] VoronoiCell;

            int iTestCase = 3;
            switch (iTestCase) {
                case 1:
                VoronoiCell = new Vector[] {
                    new Vector(0.5, 0.5),
                    new Vector(0.5, -2),
                    new Vector(-2,-2),
                    new Vector(-2, 0.5)
                };
                break;

                case 2:
                VoronoiCell = new Vector[] {
                    new Vector(-0.7, 0.5),
                    new Vector(-0.2, 0.0),
                    new Vector(-0.8,0.0)
                };
                break;

                case 3:
                VoronoiCell = new Vector[] {
                    new Vector(-1.5, 0.5),
                    new Vector(-0.2, 0.5),
                    new Vector(-0.2, 0.0),
                    new Vector(-1.5, 0.0)
                };
                break;

                default:
                throw new ArgumentOutOfRangeException();
            }

            var org = VoronoiCell.CloneAs();
            var Test = Voronoi.PolygonClipping.WeilerAthertonClipping(DomainBndy, IsInV, VoronoiCell);
            
            gp.PlotXY(org.Select(X => X.x).ToArray(), org.Select(X => X.y).ToArray(), "org",
                new PlotFormat("-ob"));
            gp.PlotXY(Test.Select(X => X.x).ToArray(), Test.Select(X => X.y).ToArray(), "intersect",
                new PlotFormat("-or"));


            gp.SetXRange(-2.2, 2.2);
            gp.SetYRange(-2.2, 2.2);
            gp.Execute();
            Console.ReadKey();

            return;
        }
        */

    }
}
