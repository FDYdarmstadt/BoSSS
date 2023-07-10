using MPI.Wrappers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Connectors;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Application.CahnHilliard;
using static BoSSS.Application.CahnHilliard.CahnHilliardMain;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.Grid.RefElements;
using System.IO;
using BoSSS.Solution.LevelSetTools.StokesExtension;
using BoSSS.Solution.NSECommon;

// TODO: test mit norm, jump norm zum laufen bringen
// erst dotnet, dann mono, dann openfoam

namespace BoSSS.Application.ExternalBinding {
    

    public class FixedOperators : IForeignLanguageProxy {

        /// <summary>
        /// 
        /// </summary>
        [CodeGenExport]
        public FixedOperators() {

        }


        IntPtr m_ForeignPtr;

        /// <summary>
        /// %
        /// </summary>
        public void _SetForeignPointer(IntPtr ptr) {
            if (ptr == IntPtr.Zero) {
                m_ForeignPtr = IntPtr.Zero;
            } else {

                if (m_ForeignPtr != IntPtr.Zero) {
                    throw new ApplicationException("already registered");
                }
                m_ForeignPtr = ptr;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IntPtr _GetForeignPointer() {
            return m_ForeignPtr;
        }


        SinglePhaseField c;
        SinglePhaseField[] velocity;
        SinglePhaseField mu;
        
        
        OpenFoamMatrix _mtx;
        OpenFoamPatchField _ptch;
        double rho0 = 1;
        double rho1 = 2;
        ScalarFunction InitFunc()
        {
            // double rMin = 2.0e-3 / sqrt(noOfTotalCells) * 3.0 / sqrt(2);
            // double radius = 0.5e-3;
            double radius = 15.0/2.0;
            // double radius = rMin * 1.3;
            return ((_3D)((x, y, z) => Math.Tanh((-Math.Sqrt(Math.Pow(x, 2) + Math.Pow(z, 2)) + Math.Pow(radius, 1)) * Math.Sqrt(2)))).Vectorize();
        }

        /// <summary>
        /// Returns the auxiliary field appearing during the solution of the Cahn-Hilliard system
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField GetFlux() {
            int D = 3;
            OpenFoamDGField Flux = new OpenFoamDGField(_ptch.Grid, 0, D);
            int nCells = _ptch.Grid.NumberOfCells;

            for (int j = 0; j < nCells; j++) {
                double cMean = c.GetMeanValue(j);
                double C0 = (cMean + 1.0)/2.0;
                double C1 = 1.0 - C0;
                for (int d = 0; d < D; d++) {
                    double massFlux = (rho0 * C0 + rho1 * C1) / 2.0 * velocity[d].GetMeanValue(j);
                    if (d == 1 && velocity[d].GetMeanValue(j) < 1e-10)
                        Console.WriteLine("Hello from BoSSS " + velocity[d].GetMeanValue(j));
                    Flux.SetDGcoordinate(d, j, 0, massFlux);
                }
            }
            return Flux;
        }

        /// <summary>
        /// Returns the auxiliary field appearing during the solution of the Cahn-Hilliard system
        /// </summary>
        [CodeGenExport]
        public OpenFoamDGField GetMu() {
            int D = 3;
            OpenFoamDGField Mu0 = new OpenFoamDGField(_ptch.Grid, 0, D);
            int nCells = _ptch.Grid.NumberOfCells;

            for (int j = 0; j < nCells; j++) {
                Mu0.SetDGcoordinate(0, j, 0, mu.GetMeanValue(j));
            }
            return Mu0;
        }

        /// <summary>
        /// Return the deformation parameter of the ellipse passing through the points given by coordinates xs and ys. If more than 6 points are given, a least-square fit is performed.
        /// </summary>
        double GetDeformParameter(double[] xs, double[] ys){
            // credit to https://stackoverflow.com/a/48002645, on which this is largely based
            // fit ellipse passing through given points xs and ys
            int len = xs.Length;
            var D = MultidimensionalArray.Create(len, 6);
            var C = MultidimensionalArray.Create(6, 6);
            for (int i = 0; i < len; i++){
                D[i, 0] = xs[i] * xs[i];
                D[i, 1] = xs[i] * ys[i];
                D[i, 2] = ys[i] * ys[i];
                D[i, 3] = xs[i];
                D[i, 4] = ys[i];
                D[i, 5] = 1;
            }
            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){
                    C[i,j] = 0.0;
                }
            }
            var DT = D.TransposeTo();
            var S = DT * D;
            C[0,2] = 2;
            C[2,0] = 2;
            C[1,1] = -1;
            var SInv = S.InvertTo();
            var SInvC = SInv * C;
            double[,] SInvCArr = SInvC.To2DArray();
            Matrix<double> SInvCMat = Matrix<double>.Build.DenseOfArray(SInvCArr);
            Evd<double> eigen = SInvCMat.Evd(Symmetricity.Asymmetric);
            double[] eigenvals = eigen.EigenValues.Real().ToArray();
            var eigenvectors = eigen.EigenVectors.ToArray();
            var n = eigenvals.IndexOfMax();
            var eigenvec = new double[6];
            for (int i = 0; i < 6; i++){
                eigenvec[i] = eigenvectors[i, n];
            }

            // calculate length of ellipse axes based on its raw parameters
            double a = eigenvec[0];
            double b = eigenvec[1]/2.0;
            double c = eigenvec[2];
            double d = eigenvec[3]/2.0;
            double f = eigenvec[4]/2.0;
            double g = eigenvec[5];

            double num=b*b-a*c;
            double cx=(c*d-b*f)/num;
            double cy=(a*f-b*d)/num;

            double angle=0.5*Math.Atan(2*b/(a-c))*180/Math.PI;
            double up = 2.0*(a*f*f+c*d*d+g*b*b-2.0*b*d*f-a*c*g);
            double down1=(b*b-a*c)*( (c-a)*Math.Sqrt(1.0+4.0*b*b/((a-c)*(a-c)))-(c+a));
            double down2=(b*b-a*c)*( (a-c)*Math.Sqrt(1.0+4.0*b*b/((a-c)*(a-c)))-(c+a));
            double longAxisLength = Math.Sqrt(Math.Abs(up/down1));
            double shortAxisLength = Math.Sqrt(Math.Abs(up/down2));
            Console.WriteLine("Long Axis: " + longAxisLength);
            Console.WriteLine("Short Axis: " + shortAxisLength);
            double deformationParameter = Math.Abs((longAxisLength - shortAxisLength) / (longAxisLength + shortAxisLength));

            return deformationParameter;
        }

        /// <summary>
        /// 1D-bisection method for finding zeros of Field c.
        /// </summary>
        double Bisect(SinglePhaseField c, double xMin, double xMax, int axis, int iter=0, double otherCoord = 0, int iters = 15){
            double[] MakeCoord(double loc){
                double yM = 0.05e-3;
                if (axis == 0){
                    return new double[]{loc, yM, otherCoord};
                } else if (axis == 2){
                    return new double[]{otherCoord, yM, loc};
                } else {
                    throw new ApplicationException("Direction is not supported");
                }
            }
            double xInt = (xMin + xMax)/2;
            double cMin = c.ProbeAt(MakeCoord(xMin));
            double cMax = c.ProbeAt(MakeCoord(xMax));
            double cInt = c.ProbeAt(MakeCoord(xInt));

            if (iter > iters){
                return Newton(c, xInt, axis, otherCoord); // after iters iterations, we should be close enough so that Newton can do the rest
            }
            if (cMin * cInt < 0){ // zero lies between xMin and xInt
                return Bisect(c, xMin, xInt, axis, iter+1, otherCoord, iters);
            }
            else if (cMax * cInt < 0){ // zero lies between xMax and xInt
                return Bisect(c, xInt, xMax, axis, iter+1, otherCoord, iters);
            } else {
                throw new ApplicationException("Something went wrong during bisect search");
            }
        }

        /// <summary>
        /// 1D-Newton method for finding zeros of Field c.
        /// </summary>
        double Newton(SinglePhaseField c, double xRoot, int axis, double otherCoord = 0, double delta = 1e-8, double tolerance = 1e-8, double trustRegion = 1e-4, int maxIter = 1000){
            int iter = 0;
            double error = 1e10;
            int direction;
            double[] MakeCoord(double loc){
                double yM = 0.05e-3;
                if (axis == 0){
                    return new double[]{loc, yM, otherCoord};
                } else if (axis == 2){
                    return new double[]{otherCoord, yM, loc};
                } else {
                    throw new ApplicationException("Direction is not supported");
                }
            }
            while (error > tolerance){
                if (iter > maxIter){
                    throw new ApplicationException("Unable to converge in " + maxIter + " iterations");
                }
                // Console.WriteLine("Iteration " + iter);
                // Console.WriteLine(xRoot);
                double f = c.ProbeAt(MakeCoord(xRoot));
                double fPlus = c.ProbeAt(MakeCoord(xRoot+delta));
                double fMinus = c.ProbeAt(MakeCoord(xRoot-delta));
                double fDiff = (fPlus - fMinus)/(2*delta);
                double newtonStep;
                if (Math.Abs(fDiff) < 1e-25){
                    direction = Math.Sign(f*xRoot);
                    newtonStep = trustRegion * direction;
                    // Console.WriteLine("Case 1");
                } else {
                    newtonStep = -f / fDiff;
                    if (Math.Abs(newtonStep) > trustRegion) {
                        // direction = Math.Sign(f*xRoot);
                        direction = Math.Sign(newtonStep);
                        newtonStep = trustRegion * direction;
                        // Console.WriteLine("Case 2");
                    }
                }
                xRoot = xRoot + newtonStep;
                error = Math.Abs(f);
                // Console.WriteLine(f);
                // Console.WriteLine(fDiff);
                // Console.WriteLine(newtonStep);
                // Console.WriteLine();
                iter++;
            }
            return xRoot;
        }

        /// <summary>
        /// Return the deformation of the ellipse described by the field c (c = 0 isosurface)
        /// </summary>
        public double GetDeformParameter(SinglePhaseField c) {
            try { // since this is basically post-processing, it should not crash the simulation
                double xMax = 0.75e-3;
                double xMin = -0.75e-3;
                double zMax = 0.75e-3;
                double zMin = -0.75e-3;
                double R = xMax/2;

                // find intersections of levelset with some axes
                double xIntersectionPos = Bisect(c, 0.0, xMax*0.99, 0);
                double xIntersectionNeg = Bisect(c, xMin*0.99, 0.0, 0);
                double zIntersectionPos = Bisect(c, 0.0, zMax*0.99, 2);
                double zIntersectionNeg = Bisect(c, zMin*0.99, 0.0, 2);
                double xIntersectionPos2 = Bisect(c, 0.0, xMax*0.99, 0, otherCoord: zIntersectionPos / 2.0);
                double xIntersectionNeg2 = Bisect(c, xMin*0.99, 0.0, 0, otherCoord: zIntersectionPos / 2.0);
                double zIntersectionNeg2 = Bisect(c, zMin*0.99, 0.0, 2, otherCoord: xIntersectionNeg / 2.0);
                double zIntersectionPos2 = Bisect(c, 0.0, zMax*0.99, 2, otherCoord: xIntersectionNeg / 2.0);
                Console.WriteLine("intersections of levelset: ");
                Console.WriteLine("    x = [ " + xIntersectionPos + ", " + xIntersectionNeg + ", " + 0 + ", " + 0 + ", " + xIntersectionPos2 + ", " + xIntersectionNeg2 + ", " + xIntersectionNeg / 2.0 + ", " + xIntersectionNeg / 2.0 + " ]");
                Console.WriteLine("    y = [ " + 0 + ", " + 0 + ", " + zIntersectionPos + ", " + zIntersectionNeg + ", " + zIntersectionPos / 2.0 + ", "+ zIntersectionPos / 2.0 + ", " + zIntersectionNeg2 + ", " + zIntersectionPos2 + " ]");
                var xs = new double[] { xIntersectionPos, xIntersectionNeg, 0, 0, xIntersectionPos2, xIntersectionNeg2, xIntersectionNeg / 2.0 , xIntersectionNeg / 2.0 };
                var ys = new double[] { 0, 0, zIntersectionPos, zIntersectionNeg,zIntersectionPos / 2.0,  zIntersectionPos / 2.0, zIntersectionNeg2, zIntersectionPos2 };
                var deformParams = new List<double>();
                double deformParam = GetDeformParameter(xs, ys);
                for (int i = 0; i < xs.Length; i++){
                    var _xs = xs.Where((val, indx) => indx != i).ToArray();
                    var _ys = ys.Where((val, indx) => indx != i).ToArray();
                    deformParams.Add(GetDeformParameter(_xs, _ys));
                }
                Console.WriteLine("Max deformParam: " + deformParams.Max());
                Console.WriteLine("Min deformParam: " + deformParams.Min());
                Console.WriteLine("deformParam relative uncertainty: " + (deformParams.Max() - deformParams.Min()) / deformParam);
                return deformParam;
            } catch (Exception e) {
                    Console.WriteLine(e);
                    return 0;
                }
        }

        /// <summary>
        /// Return the norm of the field c
        /// </summary>
        public double Norm(OpenFoamMatrix mtx = null) {
            if (mtx != null) {
                OpenFoamDGField C = mtx.Fields[0];
                c = C.Fields[0] as SinglePhaseField;
                if (InitFunc() != null){
                    c.Clear();
                    c.ProjectField(InitFunc());
                }
            }
            return c.L2Norm();
        }
        /// <summary>
        /// Return the relative change of the norm of the field c
        /// </summary>
        public double NormRelChange() {
            double postNorm = c.L2Norm();
            SinglePhaseField c0 = new SinglePhaseField(c.Basis);
            c0.Clear();
            c0.ProjectField(InitFunc());
            double preNorm = c0.L2Norm();
            return (postNorm - preNorm)/postNorm;
        }

        /// <summary>
        /// Return the jump norm of the field c
        /// </summary>
        public double JumpNorm() {
            GridData grd = (GridData)c.GridDat;
            int D = grd.SpatialDimension;

            c.MPIExchange();

            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { D + 1 }, grd,
                (new EdgeQuadratureScheme()).Compile(grd, c.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++) {
                        int iEdge = j + i0;
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if (jCell_OT >= 0) {

                            int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN, false);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT, false);

                            c.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            c.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);
                        } else {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum();

            return Unorm.Sqrt();
        }

        /// <summary>
        /// Solves the Cahn-Hilliard equation
        /// This method only contains arguments that can be made available to OpenFOAM given the limitations of the mono-C-interface.
        /// For usage exclusively within BoSSS, the method <see cref="CahnHilliardInternal"/> is generally more convenient.
        /// </summary>
        [CodeGenExport]
        public void CahnHilliard(OpenFoamMatrix mtx, OpenFoamSurfaceField Flux, OpenFoamDGField U, OpenFoamPatchField ptch, OpenFoamPatchField ptchU, double deltaT) {

            // TODO sync from OpenFOAM
             double epsilon = 1e-5; // capillary width
             // double r = 5e-4; // radius
             double r = 1; // dimensional form
             double cahn = 5e-5;

             double shearRate = 8.9235;
             // double shearRate = 89.235;
             double u = shearRate * r;
             double kappa = 5e-11;
             double sigma = 0.063;
             // double sigma = 1.0;
             // double diff = 3 * kappa * sigma / (2 * Math.Sqrt(2) * r * u * epsilon);
            //                        // double M = 1; // mobility parameter
            double M = Math.Sqrt(epsilon); // mobility parameter
            double lam = 3 / (2 * Math.Sqrt(2)) * sigma * epsilon; // Holger's lambda
            double diff = M * lam * 10;

            bool convection = true;
            bool stat = false;

            // TODO only for 1D test
            // cahn = 0.1;
            // diff = 0.1;
            // convection = false;
            // stat = true;
            CahnHilliardParameters chParams = new CahnHilliardParameters(_dt: deltaT, _diffusion: diff, _cahn: cahn, _stationary: stat, _endT: deltaT*1.1, _convection: convection);
            CahnHilliardInternal(mtx, Flux, U, ptch, ptchU, null, chParams);
        }
        // public static bool FirstTimeStep = true;

        static void VerifyTrackerState(LevelSetTracker lsTrk) {
            var NoOfCells_a = lsTrk.Regions.GetSpeciesMask("a").NoOfItemsLocally.MPISum();
            var NoOfCells_b = lsTrk.Regions.GetSpeciesMask("a").NoOfItemsLocally.MPISum();
            var NoOfCells_cut = lsTrk.Regions.GetCutCellMask().NoOfItemsLocally.MPISum();

            if (NoOfCells_a > 0 && NoOfCells_b > 0 && NoOfCells_cut == 0)
                throw new ArithmeticException("Level-Set yields positive and negative cells, but no cut-cells can be identified; probably the level-set was just polynomial degree 0 (this does not work)");
        }

        static bool FirstSolve = true;


        /// <summary>
        /// Solves the Cahn-Hilliard equation
        /// This method also contains arguments that cannot be made available to OpenFOAM due to limitations of the mono-C-interface.
        /// </summary>
        public void CahnHilliardInternal(OpenFoamMatrix mtx, OpenFoamSurfaceField Flux, OpenFoamDGField U, OpenFoamPatchField ptch, OpenFoamPatchField ptchU, ScalarFunction func = null, CahnHilliardParameters chParams = new CahnHilliardParameters()) {
            try {
            // {
                _mtx = mtx;
                _ptch = ptch;

                /*
                double pow(double x, int e)
                { // inefficient but who cares
                    // return Math.Tanh(x);
                    double acc = 1.0;
                    for (int i = 0; i < e; i++)
                        acc *= x;
                    return acc;
                }
                */


                // grid, etc
                // =========

                GridData grd = mtx.ColMap.GridDat as GridData;
                // PlotGrid("grid.plt", grd);

                Basis b = mtx.ColMap.BasisS[0];
                Basis bMu = mtx.ColMap.BasisS[0];
                int noOfTotalCells = b.GridDat.Grid.NumberOfCells;

                OpenFoamDGField C = mtx.Fields[0];
                var u = U.Fields[0] as SinglePhaseField;
                var v = U.Fields[1] as SinglePhaseField;
                var w = U.Fields[2] as SinglePhaseField;
                velocity = new[] { u, v, w };
                // u.ProjectField(((_3D)((x, y, z) => 0.05)).Vectorize());
                // v.ProjectField(((_3D)((x, y, z) => 0.0)).Vectorize());
                // w.ProjectField(((_3D)((x, y, z) => 0.0)).Vectorize());
                // OpenFOAMGrid ofGrid = ptch.Grid;
                // OpenFoamDGField fields = new(ofGrid, b.Degree, 2);
                // OpenFoamMatrix fullMtx = new(ofGrid, fields);

                var map = new UnsetteledCoordinateMapping(b, bMu);



                // op.LinearizationHint = LinearizationHint.GetJacobiOperator;


                /*
                var RealLevSet = new LevelSet(b, "Levset");
                */

                int J = b.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                c = C.Fields[0] as SinglePhaseField;

                if (c.L2Norm() < 1e-20) {
                    if (func == null) {
                        Console.WriteLine("No initialization function given - using droplet");
                        func = InitFunc();
                    }
                    ScalarFunction UInitFunc() {
                        return ((_3D)((x, y, z) => 0.01*z)).Vectorize();
                    }
                    Console.WriteLine("Zero order parameter field encountered - initializing with given function");
                    c.Clear();
                    u.Clear();
                    v.Clear();
                    w.Clear();
                    u.ProjectField(UInitFunc());

                    // var cP0 = new SinglePhaseField(new Basis(c.GridDat, 0));
                    // cP0.ProjectField(func);
                    // c.AccLaidBack(1.0, cP0);

                    c.ProjectField(func);
                }
                mu = new SinglePhaseField(b);
                // for (int j = 0; j < J; j++)
                // {
                //     int N = b.GetLength(j);
                //     for (int n = 0; n < N; n++) {
                //         c.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                //         // mu.Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                //     }
                // }

                // foreach (var val in new double[]{0, 0.5, 1, 1.5, 1.99, 2.01 ,10})
                //     Console.WriteLine("Correct Value: " + Math.Tanh(val) + "My value: " + tanh(val));

                // // var LSGrad = new MultidimensionalArray(noOfTotalCells);
                // var nsma = new MultidimensionalArray(2);
                // nsma.Allocate(noOfTotalCells, 3);
                // var ns = new NodeSet(Cube.Instance, nsma, false);
                // var LSGrad = new MultidimensionalArray(3);
                // LSGrad.Allocate(noOfTotalCells, ns.NoOfNodes, 3);
                // RealLevSet.EvaluateGradient(0, noOfTotalCells, ns, LSGrad); <- This causes a mono exception
                // var LSGradX = new SinglePhaseField(b);
                // var LSGradY = new SinglePhaseField(b);
                // var LSGradZ = new SinglePhaseField(b);
                // for (int i = 0; i < noOfTotalCells; i++)
                // {
                //     LSGradX.SetMeanValue(i, LSGrad[i, 0, 0]);
                //     LSGradY.SetMeanValue(i, LSGrad[i, 0, 1]);
                //     LSGradZ.SetMeanValue(i, LSGrad[i, 0, 2]);
                //     // LSGradX.SetMeanValue(i, LSGrad[0,i]);
                //     // LSGradY.SetMeanValue(i, LSGrad[1,i]);
                //     // LSGradZ.SetMeanValue(i, LSGrad[2,i]);
                // }
                // RealLevSet.EvaluateGradient
                var domfields = (IReadOnlyDictionary<string, DGField>)(new Dictionary<string, DGField>() { { "c", c }, { "mu", mu } });
                var paramfields = (IReadOnlyDictionary<string, DGField>)(new Dictionary<string, DGField>(){
                        {"VelocityX", u},
                        {"VelocityY", v},
                        {"VelocityZ", w},
                        {"c0", c},
                        // {"LevelSetGradient[0]", LSGradX},
                        // {"LevelSetGradient[1]", LSGradY},
                        // {"LevelSetGradient[2]", LSGradZ}
                        {"LevelSetGradient[0]", c},
                        {"LevelSetGradient[1]", c},
                        {"LevelSetGradient[2]", c}
                        });
                Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> GetNamedInputFields = delegate (DGField[] fields) {

                    return (domfields, paramfields);
                };

                /*
                RealLevSet.Clear();
                RealLevSet.Acc(1.0, c);
                LevelSetUpdater lsu = new LevelSetUpdater(grd, XQuadFactoryHelper.MomentFittingVariants.Classic,
                                                         2, new string[] { "a", "b" },
                                                         GetNamedInputFields,
                                                         RealLevSet, "c", ContinuityProjectionOption.ConstrainedDG);
                // note on continuity projection: 
                // - For the Stokes Extension, we only require level-set-surface integrals; no cut-edge integrals are required;
                // - therefore, the level-set does not need to be strictly continuous.
                // => ContinuityProjectionOption.None should be ok.
                lsu.EnforceContinuity();
                
                var RealTracker = lsu.Tracker;
                // mu.Laplacian(-cahn, c);
                // mu.Acc(-1.0, c);
                // mu.ProjectPow(1.0, c, 3.0);
                // RealLevSet.Clear();
                // RealLevSet.Acc(1.0, c);

                foreach (var _ls in RealTracker.LevelSets) {
                    var _dgls = _ls as LevelSet;
                    if (_dgls != null) {
                        if (_dgls.L2Norm() == 0) {
                            throw new ArithmeticException("level-set is exactly zero");
                        }
                    }
                }
                

                RealTracker.UpdateTracker(0);
                VerifyTrackerState(RealTracker);
                */

                // SubGrid subgr = RealTracker.Regions.GetNearFieldSubgrid(6);
                // SubGridBoundaryModes subgrbnd = 0;
                // CellMask subgrMask = subgr.VolumeMask;
                // CellMask fullMask = CellMask.GetFullMask(grd);
                // SubGrid fullSubGrd = new SubGrid(fullMask);

                CellMask mask = null;
                //SubGrid sgrid = null;
                // CellMask mask = fullMask;
                // SubGrid sgrid = fullSubGrd;
                // CellMask mask = subgrMask;
                // SubGrid sgrid = subgr;

                //int noOfNarrowBandCells = 0;
                //if (mask != null) {
                //    foreach (bool cellInNarrowBand in mask.GetBitMask()) {
                //        if (cellInNarrowBand) {
                //            noOfNarrowBandCells++;
                //        }
                //    }
                //    if (noOfNarrowBandCells == 0) {
                //        Console.WriteLine("Solving only in a narrow band containing " + noOfNarrowBandCells + " of " + noOfTotalCells + " cells");
                //        // throw new ApplicationException("No interface found");
                //        // mask = fullMask;
                //        // sgrid = fullSubGrd;
                //        // Console.WriteLine("No narrow band cells detected, solving on the whole domain");
                //    } else {
                //        Console.WriteLine("Solving only in a narrow band containing " + noOfNarrowBandCells + " of " + noOfTotalCells + " cells");
                //    }
                //}

                //System.Collections.BitArray subGridCellMask = mask?.GetBitMask();


                // Assembly of Cahn Hilliard operator
                // ==================================
                int nParams = 7;
                SpatialOperator CahnHillOp;
                {
                    double lambda = 0.0;
                    double penalty_const = 2.6;
                    double diff = chParams.Diffusion;
                    double cahn = chParams.Cahn;


                    // var op = new SpatialOperator(2, nParams, 2, QuadOrderFunc.Linear(), "c", "mu", "VelocityX", "VelocityY", "VelocityZ", "c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_mu");
                    CahnHillOp = new SpatialOperator(2, nParams, 2, QuadOrderFunc.NonLinearWithoutParameters(3), "c", "mu", "VelocityX", "VelocityY", "VelocityZ", "c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_mu");
                    // var op = new XSpatialOperatorMk2(2, nParams, 2, QuadOrderFunc.Linear(), new List<string>{"a", "b"}, "c", "mu", "VelocityX", "VelocityY", "VelocityZ","c0", "LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "Res_c", "Res_mu");
                    // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "mu", "c0", "VelocityX", "VelocityY", "VelocityZ", "c_Res", "mu_Res");
                    // var op = new SpatialOperator(2, 4, 2, QuadOrderFunc.Linear(), "c", "mu", "c0","LevelSetGradient[0]", "LevelSetGradient[1]", "LevelSetGradient[2]", "c_Res", "mu_Res");


                    var CHCdiff = new CahnHilliardCDiff(ptch, penalty_const, diff, lambda, mask);
                    var CHCconv = new CahnHilliardCConv(() => velocity, null, mask);
                    var CHMudiff = new CahnHilliardMuDiff(ptch, penalty_const, cahn, mask);
                    var CHMusource = new CahnHilliardMuSource(mask);
                    CahnHillOp.EquationComponents["Res_c"].Add(CHCdiff);
                    CahnHillOp.EquationComponents["Res_c"].Add(CHCconv);
                    CahnHillOp.EquationComponents["Res_mu"].Add(CHMusource);
                    CahnHillOp.EquationComponents["Res_mu"].Add(CHMudiff);
                    // op.LinearizationHint = LinearizationHint.GetJacobiOperator;

                    double[] MassScales = new double[2];
                    MassScales[0] = 1.0;
                    // MassScales[1] = 1.0;
                    CahnHillOp.TemporalOperator = new ConstantTemporalOperator(CahnHillOp, MassScales);

                    CahnHillOp.LinearizationHint = LinearizationHint.GetJacobiOperator;

                    CahnHillOp.Commit();
                }



             

                // Boundary Condition map
                // ======================

                IncompressibleBoundaryCondMap BCmap;
                {

                    // TODO
                    var leftBVC = new AppControl.BoundaryValueCollection();
                    leftBVC.type = "Pressure_Outlet";
                    var rightBVC = new AppControl.BoundaryValueCollection();
                    rightBVC.type = "Pressure_Outlet";
                    var topBVC = new AppControl.BoundaryValueCollection();
                    topBVC.type = "Wall";
                    // double shearRate = 0.89235;
                    double shearRate = 0.01;
                    topBVC.Evaluators.Add("VelocityX", (X,t) => shearRate*X[2]);
                    var bottomBVC = new AppControl.BoundaryValueCollection();
                    bottomBVC.type = "Wall";
                    bottomBVC.Evaluators.Add("VelocityX", (X,t) => shearRate*X[2]);
                    var fbBVC = new AppControl.BoundaryValueCollection();
                    // fbBVC.type = IncompressibleBcType.SlipSymmetry.ToString();
                    fbBVC.type = IncompressibleBcType.FreeSlip.ToString();
                    // fbBVC.type = IncompressibleBcType.Pressure_Outlet.ToString(); // leads to y-velocity after 1st timestep
                    // fbBVC.type = IncompressibleBcType.Wall.ToString();
                    var bcmapcollection = new Dictionary<string, AppControl.BoundaryValueCollection>() {
                            { "left", leftBVC},
                            { "right", rightBVC},
                            { "bottom", bottomBVC},
                            { "top", topBVC },
                            { "frontAndBack", fbBVC},
                            // { "front", fbBVC},
                            // { "back", fbBVC},
                    };
                    // string[] bndFuncName = new string[]{"left", "right", "bottom", "top"};
                    BCmap = new IncompressibleBoundaryCondMap(grd, bcmapcollection, PhysicsMode.Incompressible);
                }

                

                // perform Stokes Extension
                // ========================
                StokesExtension stokesExt;
                SinglePhaseField[] uStokes;
                {
                    stokesExt = new StokesExtension(3, BCmap, 3, 0.0, true, useBCMap: true);
                    uStokes = velocity.Select(Vel_d => Vel_d.CloneAs()).ToArray();
                  /*  stokesExt.SolveExtension(0, RealTracker, uStokes, velocity); */
                    // stokesExt.SolveExtension(0, RealTracker, velocity, velocity);
                }


                /*
                lsu.InitializeParameters(domfields, paramfields);
                */

                // var tp = new Tecplot(grd.Grid.GridData, 3);
                // Tecplot("plot.1", 0.0, 3, c, mu, RealLevSet, u, v, w, uStokes[0], uStokes[1], uStokes[2]);
                uStokes[0].Identification = VariableNames.Velocity0X;
                uStokes[1].Identification = VariableNames.Velocity0Y;
                uStokes[2].Identification = VariableNames.Velocity0Z;
                u.Identification = VariableNames.VelocityX;
                v.Identification = VariableNames.VelocityY;
                w.Identification = VariableNames.VelocityZ;

                

                // Timestepper initialization
                // ==========================
                XdgTimestepping TimeStepper;
                {
                    SinglePhaseField Res_c = new SinglePhaseField(b);
                    SinglePhaseField Res_mu = new SinglePhaseField(b);
                    List<DGField> ParameterMap = new List<DGField>();
                    for (int i = 0; i < nParams; i++) {
                        ParameterMap.Add(new SinglePhaseField(b));
                    }
                    for (int j = 0; j < J; j++) {
                        int N = b.GetLength(j);
                        for (int n = 0; n < N; n++)
                            ParameterMap[3].Coordinates[j, n] += C.GetDGcoordinate(0, j, n);
                    }
                    for (int j = 0; j < 3; j++) {
                        ParameterMap[j] = new[] { u, v, w }[j];
                    }



                    NonLinearSolverConfig nls = new NonLinearSolverConfig();
                    var ls = new DirectSolver.Config();
                    nls.SolverCode = NonLinearSolverCode.Newton;
                    // nls.SolverCode = NonLinearSolverCode.Picard;
                    // nls.ConvergenceCriterion = 1e-5;
                    // nls.MaxSolverIterations = 100;
                    nls.verbose = true;


                    TimeStepper = new XdgTimestepping(CahnHillOp,
                                                      new SinglePhaseField[] { c, mu },
                                                      new SinglePhaseField[] { Res_c, Res_mu },
                                                      // TimeSteppingScheme.CrankNicolson,
                                                      TimeSteppingScheme.ImplicitEuler,
                                                      null,
                                                      null,
                                                      LinearSolver: ls,
                                                      NonLinearSolver: nls
                                                      // RealTracker,
                                                      // ,lsu: (() => lsu)
                                                      // _LevelSetHandling: LevelSetHandling.LieSplitting,
                                                      // _LevelSetHandling: LevelSetHandling.Coupled_Once,
                                                      // _AgglomerationThreshold: 0.0
                                                      );
                }

                // temporal integration
                // ====================

                {
                    double endTime = chParams.endT;
                    double dt = chParams.dt;
                    double time = 0.0;
                    int t = 0;

                    double cMean0 = c.IntegralOver(null);
                    SinglePhaseField cVar = new SinglePhaseField(b);
                    cVar.Clear();
                    cVar.Acc(1.0, c);
                    cVar.ProjectPow(1.0, cVar, 2.0);
                    cVar.AccConstant(-Math.Pow(cMean0, 2.0));
                    double cVarVal = cVar.IntegralOver(null);
                    Console.WriteLine("cMean: " + cMean0);
                    Console.WriteLine("cVar: " + cVarVal);
                    double deformParameter;
                    if (FirstSolve){
                        Tecplot("plot.1", 0.0, 3, c, mu, u, v, w);
                        TimeStepper.Solve(time, dt);
                        Tecplot("plot.2", 0.0, 3, c, mu, u, v, w);
                        deformParameter = GetDeformParameter(c);
                        Console.WriteLine("Deformation Parameter: " + deformParameter);
                        t++;
                        time += dt;
                        FirstSolve = false;
                    }
                    var RealLevSet = new LevelSet(b, "Levset");
                    RealLevSet.Clear();
                    RealLevSet.Acc(1.0, c);
                    LevelSetUpdater lsu;
                    lsu = new LevelSetUpdater(grd, XQuadFactoryHelper.MomentFittingVariants.Classic,
                                                             2, new string[] { "a", "b" },
                                                             GetNamedInputFields,
                                                             RealLevSet, "c", ContinuityProjectionOption.None);
                    // note on continuity projection:
                    // - For the Stokes Extension, we only require level-set-surface integrals; no cut-edge integrals are required;
                    // - therefore, the level-set does not need to be strictly continuous.
                    // => ContinuityProjectionOption.None should be ok.
                    lsu.EnforceContinuity();
                    lsu.InitializeParameters(domfields, paramfields);
                    var RealTracker = lsu.Tracker;

                    foreach (var _ls in RealTracker.LevelSets)
                    {
                        var _dgls = _ls as LevelSet;
                        if (_dgls != null)
                        {
                            if (_dgls.L2Norm() == 0)
                            {
                                throw new ArithmeticException("level-set is exactly zero");
                            }
                        }
                    }

                    RealTracker.UpdateTracker(0);
                    VerifyTrackerState(RealTracker);

                    while (time < endTime) {
                        RealLevSet.Clear();
                        RealLevSet.Acc(1.0, c);
                        RealTracker.UpdateTracker(time);
                        VerifyTrackerState(RealTracker);
                        uStokes = velocity.Select(Vel_d => Vel_d.CloneAs()).ToArray();
                        // Tecplot("plotb4." + (t + 2), time, 3, c, mu, u, v, w);
                        if (chParams.Convection){
                            stokesExt.SolveExtension(0, RealTracker, uStokes, velocity);
                        }
                        TimeStepper.Solve(time, dt);

                        double cMean = c.IntegralOver(null);
                        cVar.Clear();
                        cVar.Acc(1.0, c);
                        cVar.ProjectPow(1.0, cVar, 2.0);
                        cVar.AccConstant(-Math.Pow(cMean, 2.0));
                        cVarVal = cVar.IntegralOver(null);
                        Console.WriteLine("cMean change compared to starting configuration: " + (cMean - cMean0));
                        Console.WriteLine("relative cMean change compared to starting configuration: " + (cMean - cMean0)/cMean0);
                        Console.WriteLine("cVar: " + cVarVal);
                        deformParameter = GetDeformParameter(c);
                        Console.WriteLine("Deformation Parameter: " + deformParameter);

                        Tecplot("plot." + (t + 2), time, 3, c, mu, u, v, w);

                        time += dt;
                        t++;

                    }
                }
            } catch (Exception e) {
               Console.WriteLine(e.GetType());
               Console.WriteLine(e.Message);
               Console.WriteLine(e.StackTrace);
               Console.WriteLine(e);
               throw new AggregateException(e);
            }
            // }
        }

        /// <summary>
        /// 
        /// </summary>
        [CodeGenExport]
        // public void Laplacian(OpenFoamMatrix mtx) {
        public void Laplacian(OpenFoamMatrix mtx, OpenFoamPatchField ptch) {

            // grid, etc
            // =========

            GridData grd = mtx.ColMap.GridDat as GridData;
            // PlotGrid("grid.plt", grd);

            var b = mtx.ColMap.BasisS[0];
            var map = new UnsetteledCoordinateMapping(b);

            var L = new Laplace(1.3, ptch);
            var op = new SpatialOperator(1, 0, 1, QuadOrderFunc.Linear(), "T", "c1");

            op.EquationComponents["c1"].Add(L);
            op.Commit();

            // evaluate operator
            // =================

            var eval = op.GetMatrixBuilder(map, null, map);
            eval.ComputeMatrix(mtx, mtx.RHSbuffer);
            mtx.RHSbuffer.ScaleV(-1); // convert LHS affine vector to RHS

            Console.WriteLine("Computed Laplacian Matrix, norm is " + mtx.InfNorm());
            // Console.WriteLine("Computed Laplacian Matrix, RHS is ");
            // foreach (var elem in mtx.RHSbuffer)
            //     Console.Write(elem + " ");
            // Console.WriteLine();
        }

        static public void PlotGrid(string filename, IGridData grd) {


            string SanitizeName(string s) {
                char[] ot = s.ToCharArray();
                for(int k = 0; k < ot.Length; k++) {
                    if(char.IsWhiteSpace(ot[k])) {
                        ot[k] = '_';
                    }

                    if (ot[k] == '(')
                        ot[k] = 'L';
                    if (ot[k] == ')')
                        ot[k] = 'R';
                }
                return new string(ot);
            }



            var et2Name = grd.EdgeTagNames;
            Console.WriteLine($"Grid containing {et2Name.Count} EdgeTag names: ");
            int i = 0;
            foreach (var t in et2Name) { // loop over all different edge tag names...
                string name = t.Value;
                byte tag2color = t.Key;

                string sname = SanitizeName(name);

                if (name.Equals(sname)) {
                    Console.WriteLine($"   {i}: {name} -- tag = {tag2color}");
                } else {
                    Console.WriteLine($"   {i}: {name} -- tag = {tag2color}   (marked as '{sname}') in output file.");
                }
                i++;
            }


            var B0 = new Basis(grd, 0);
            SinglePhaseField[] bndyMarkers = new SinglePhaseField[et2Name.Count + 1];

            int[,] Edge2GeomCell = grd.iGeomEdges.CellIndices;
            int[] G2L = grd.iGeomCells.GeomCell2LogicalCell;
            byte[] EdgeTags = grd.iGeomEdges.EdgeTags;

            i = 0;
            foreach(var t in et2Name) { // loop over all different edge tag names...
                string name = t.Value;
                byte tag2color = t.Key;
                string sname = SanitizeName(name);


                var FI = new SinglePhaseField(B0, "Marker-" + sname);
                bndyMarkers[i] = FI;
                i++;

                for (int e = 0; e < EdgeTags.Length; e++) { // loop over edges...
                    byte tag_e = EdgeTags[e];

                    if(tag_e == tag2color) {
                        // mar cells next to edge e

                        foreach(int jG in Edge2GeomCell.GetRow(e)) {
                            if (jG < 0)
                                continue;

                            // convert geometrical cell index to logical cell index
                            int jL;
                            if (G2L == null)
                                jL = jG;
                            else
                                jL = G2L[jG];

                            // color respective cell
                            FI.SetMeanValue(jL, tag2color);
                        }

                    }
                }

            }

            var dummy = new SinglePhaseField(B0, "DummyData");
            bndyMarkers[bndyMarkers.Length - 1] = dummy;

            Tecplot(filename, 0.0, 0, bndyMarkers);
        }

        static public void Tecplot(string filename, double time, int supersampling, params BoSSS.Foundation.DGField[] flds) {
            if (flds == null || flds.Length <= 0) {
                Console.WriteLine("No DG fields specified - not writing any output files.");
                return;
            }

            if (supersampling > 3) {
                Console.WriteLine("Plotting with a supersampling greater than 3 is deactivated because it would very likely exceed this machines memory.");
                Console.WriteLine("Higher supersampling values are supported by external plot application.");
                supersampling = 3;
            }

            string directory = "~/";
            string FullPath;
            // if (directory == null || directory.Length <= 0) {
            //     directory = CurrentDocDir ?? "";
            //     FullPath = Path.Combine(directory, filename);
            // } else {
            FullPath = filename;
            // }

            Console.WriteLine("Writing output file {0}...", FullPath);


            BoSSS.Solution.Tecplot.Tecplot.PlotFields(flds, FullPath, time, supersampling);

        }

        /// <summary>
        /// SIP-form for the Laplacian
        /// </summary>
        class Laplace : BoSSS.Solution.NSECommon.SIPLaplace {

            OpenFoamPatchField _ptch;

            /// <summary>
            /// 
            /// </summary>
            public Laplace(double penalty_const, OpenFoamPatchField ptch)
                : base(penalty_const, "T")
            {
                _ptch = ptch;
            }

            /// <summary>
            /// always true
            /// </summary>
            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                return _ptch.IsDirichlet(inp.EdgeTag % GridCommons.FIRST_PERIODIC_BC_TAG);
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) {
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Index of a boundary edge should not be zero");
                }
                return _ptch.Values[inp.EdgeTag % GridCommons.FIRST_PERIODIC_BC_TAG - 1][0];
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                return 0;
            }
        }

        class CahnHilliardMuSource : mu_Source {

            // OpenFoamPatchField _ptch;

            public CahnHilliardMuSource(CellMask Subgrid = null)
                : base("c"){
            }
            // public bool GetIsDiri(ref CommonParamsBnd inp){
            //     return this.IsDirichlet(ref inp);
            // }
        }

        class CahnHilliardMuDiff : mu_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardMuDiff(OpenFoamPatchField ptch, double penalty_const, double __cahn, CellMask Subgrid = null)
                : base(3, penalty_const, __cahn, null, "c"){
                _ptch = ptch;
            }

            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                // return !_ptch.IsDirichlet(inp.EdgeTag);
                // return _ptch.IsDirichlet(inp.EdgeTag);
                return true;
                // return false;
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) {
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Index of a boundary edge should not be zero");
                }
                // Console.WriteLine("diriMu " + _ptch.Values[inp.EdgeTag - 1][0]);
                return 0.0;
                // return _ptch.Values[inp.EdgeTag - 1][0];
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                // throw new Exception("g_neum should not be called");
                return 0;
            }
        }

        class CahnHilliardCDiff : __c_Diffusion {

            OpenFoamPatchField _ptch;

            public CahnHilliardCDiff(OpenFoamPatchField ptch, double penalty_const, double __diff, double __lambda, CellMask Subgrid = null)
                : base(3, penalty_const, __diff, __lambda, null){
                this._ptch = ptch;
            }
            protected override bool IsDirichlet(ref CommonParamsBnd inp) {
                // return true;
                return _ptch.IsDirichlet(inp.EdgeTag % GridCommons.FIRST_PERIODIC_BC_TAG);
            }
            public bool GetIsDiri(int et){
                return _ptch.IsDirichlet(et % GridCommons.FIRST_PERIODIC_BC_TAG);
            }

            /// <summary>
            /// Dirichlet boundary value
            /// </summary>
            override protected double g_Diri(ref Foundation.CommonParamsBnd inp) { // TODO generalize
                // throw new Exception("g_diri should not be called");
                if (inp.EdgeTag == 0){
                    throw new ApplicationException("Edge Tag of a boundary edge should not be zero");
                }
                throw new Exception("g_diri should not be called");
                // Console.WriteLine("diriC " + _ptch.Values[inp.EdgeTag - 1][0]);
                // Console.WriteLine(_ptch.Values.Count);
                // Console.WriteLine(inp.EdgeTag % GridCommons.FIRST_PERIODIC_BC_TAG);
                // Console.WriteLine(inp.EdgeTag);
                return _ptch.Values[inp.EdgeTag % GridCommons.FIRST_PERIODIC_BC_TAG - 1][0];
                // return 1.0;
            }

            /// <summary>
            /// Neumann boundary value
            /// </summary>
            override protected double g_Neum(ref Foundation.CommonParamsBnd inp) {
                // throw new Exception("g_neum should not be called");
                return 0;
            }
        }
        class CahnHilliardCConv : __c_Flux {
            public CahnHilliardCConv(Func<DGField[]> VelocityGetter, OpenFoamSurfaceField Flux, CellMask Subgrid = null)
                : base(3, VelocityGetter, null){
                // this.m_Flux = Flux;
            }

            // OpenFoamSurfaceField m_Flux;

            public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN)
            {
                // expand for treatment of input functions, for now hardcode to -1.0
                double UxN = 0;
                UxN += (inp.Parameters_IN[0]) * inp.Normal[0];

                double phi = -1.0;
                // if (UxN >= 0)
                // {
                //     phi = _uIN[0];
                // }
                // else
                // {
                //     phi = -1.0;//m_bndFunc[inp.EdgeTag](inp.X, inp.time);
                // }

                return phi * UxN * _vIN;
                // // expand for treatment of input functions, for now hardcode to -1.0
                // double referenceValue = base.BoundaryEdgeForm(ref inp, _uIN, _Grad_uIN, _vIN, _Grad_vIN);
                // // Console.WriteLine("Test1");
                // if (m_Flux == null) {
                //     Console.WriteLine("Warning: flux not passed to BoSSS - using low-quality velocity field");
                //     return base.BoundaryEdgeForm(ref inp, _uIN, _Grad_uIN, _vIN, _Grad_vIN);
                // }
                // double UxN = m_Flux.GetFlux(ref inp);
                // // Console.Write("Test2");
                // // Console.Write(UxN);

                // double phi;
                // if (UxN >= 0)
                // {
                // // Console.Write("Test3");
                //     phi = _uIN[0];
                // }
                // else
                // {
                // // Console.Write("Test4");
                //     phi = -1.0;//m_bndFunc[inp.EdgeTag](inp.X, inp.time);
                // }

                // using (var fs = new FileStream("./out.txt", FileMode.Append))
                // using (var sw = new StreamWriter(fs))
                // {
                //     sw.WriteLine("flux: " + (double)(phi * UxN * (_vIN)));
                //     sw.WriteLine("velocityfield (boundary): " + (double)(referenceValue));
                // }
                // // Console.WriteLine("Difference between flux and velocityfield (boundary): " + (double)(phi * UxN * _vIN - referenceValue));
                // // Console.Write("Test5");
                // return phi * UxN * _vIN;
            }

            // public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
            // {
            //     // double referenceValue = base.InnerEdgeForm(ref inp, _uIN, _uOUT, _Grad_uIN, _Grad_uOUT, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            //     // Console.WriteLine("Test1inner");
            //     if (m_Flux == null) {
            //         // Console.WriteLine("Warning: flux not passed to BoSSS - using low-quality velocity field");
            //         return base.InnerEdgeForm(ref inp, _uIN, _uOUT, _Grad_uIN, _Grad_uOUT, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            //     }
            //     double UxN = 0.5 * (m_Flux.GetFluxIN(ref inp) + m_Flux.GetFluxOUT(ref inp));
            //     // Console.Write("Testinner2");

            //     double UxNReference = 0;
            //     for (int d = 0; d < 3; d++)
            //         UxNReference += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];

            //     double phi;
            //     if (UxN >= 0)
            //     {
            //         // Console.Write("Test3inner");
            //         phi = _uIN[0];
            //     }
            //     else
            //     {
            //         // Console.Write("Test4inner");
            //         phi = _uOUT[0];
            //     }

            //     using (var fs = new FileStream("./out.txt", FileMode.Append))
            //     using (var sw = new StreamWriter(fs))
            //     {
            //         // sw.WriteLine("flux: " + (double)(phi * UxN * (_vIN - _vOUT)));
            //         // sw.WriteLine("velocityfield (inner): " + (double)(referenceValue));
            //         sw.WriteLine("flux: " + (double)(UxN));
            //         sw.WriteLine("velocityfield (inner): " + (double)(UxNReference));
            //     }
            //     // Console.Write("Test5inner");
            //     return phi * UxN * (_vIN - _vOUT);
            // }

        }

        // class CahnHilliardCSource : c_Source {

        //     OpenFoamPatchField _ptch;

        //     public CahnHilliardCSource(OpenFoamPatchField ptch)
        //         : base(1.0){
        //         this._ptch = ptch;
        //     }
        // }
    }
}
