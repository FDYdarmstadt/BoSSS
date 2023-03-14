﻿using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG.Quadrature
{
    class SayeGaussRule_Cube :
        SayeComboIntegrand<LinearPSI<Cube>, LinearSayeSpace<Cube>>,
        ISayeGaussRule,
        ISayeGaussComboRule
    {
        LevelSetTracker.LevelSetData lsData;

        IRootFindingAlgorithm rootFinder;

        IGridData grid;

        int iKref;

        public enum QuadratureMode { Surface, PositiveVolume, NegativeVolume, Combo };

        QuadratureMode mode;

        public SayeGaussRule_Cube(
            LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm RootFinder,
            QuadratureMode Mode)
        {
            this.lsData = _lsData;
            this.rootFinder = RootFinder;
            mode = Mode;
            grid = lsData.GridDat;
            iKref = lsData.GridDat.iGeomCells.RefElements.IndexOf(this.RefElement, (A, B) => object.ReferenceEquals(A, B));
        }

        #region ISayeGaussRule

        public int order { get; set; }

        public QuadRule Evaluate(int Cell)
        {
            LinearSayeSpace<Cube> startArg = CreateStartSetup();
            return Evaluate(Cell, startArg);
        }

        private LinearSayeSpace<Cube> CreateStartSetup()
        {
            bool IsSurfaceIntegral = (mode == QuadratureMode.Surface);

            LinearPSI<Cube> psi = new LinearPSI<Cube>(Cube.Instance);
            int domainSign = 1;
            if (mode == QuadratureMode.NegativeVolume) {
                domainSign = -1;
            }
            Tuple<LinearPSI<Cube>, int> psi_i_s_i = new Tuple<LinearPSI<Cube>, int>(psi, domainSign);
            LinearSayeSpace<Cube> arg = new LinearSayeSpace<Cube>(Cube.Instance, psi_i_s_i, IsSurfaceIntegral);
            arg.Reset();
            return arg;
        }

        public RefElement RefElement {
            get {
                return Cube.Instance;
            }
        }

        #endregion

        #region ISayeGaussComboRule

        public QuadRule[] ComboEvaluate(int Cell)
        {
            LinearSayeSpace<Cube> startArg = CreateStartSetup();
            return ComboEvaluate(Cell, startArg);
        }

        #endregion

        #region SayeComboIntegrand

        protected override ISayeQuadRule GetEmptyQuadratureRule()
        {
            ISayeQuadRule emptyRule = new NodesAndWeights(RefElement);
            return emptyRule;
        }

        #endregion

        #region Build SayeIntegrand

        protected override int EvaluateSign(MultidimensionalArray gradient, int heightDirection, int s_i, bool surface, int sigma)
        {
            int m = Math.Sign(gradient[heightDirection]);
            if (surface == true || sigma * s_i == m)
            {
                return sigma * m;
            }
            else
            {
                return 0;
            }
        }

        protected override LinearPSI<Cube>[] ExtractSubPsis(LinearPSI<Cube> psi, LinearSayeSpace<Cube> arg, int heightDirection)
        {
            double[] bounds = arg.GetBoundaries(heightDirection);
            double x_U = bounds[1];
            double x_L = bounds[0];
            LinearPSI<Cube> Psi_U = psi.ReduceDim(heightDirection, x_U);
            LinearPSI<Cube> Psi_L = psi.ReduceDim(heightDirection, x_L);
            return new LinearPSI<Cube>[] { Psi_L, Psi_U };
        }

        protected override MultidimensionalArray Gradient(LinearPSI<Cube> psi, NodeSet Node, int Cell)
        {
            NodeSet nodeOnPsi = psi.ProjectOnto(Node);
            MultidimensionalArray gradient = ReferenceGradient(nodeOnPsi, Cell);
            return gradient;
        }

        protected override double EvaluateAt(LinearPSI<Cube> Psi, NodeSet Node, int cell)
        {
            Debug.Assert(Node.NoOfNodes == 1);

            NodeSet nodeOnPsi = Psi.ProjectOnto(Node);
            double value = lsData.GetLevSetValues(nodeOnPsi, cell, 1)[0, 0];
            return value;
        }


        double boundSafety = 1.1;
        /// <summary>
        /// First order approximation of  delta >= sup_x|psi(x) - psi(x_center)|  
        /// </summary>
        /// <param name="Arg"></param>
        /// <param name="psi"></param>
        /// <param name="x_center"></param>
        /// <param name="cell"></param>
        /// <returns></returns>
        protected override double EvaluateBounds(LinearSayeSpace<Cube> Arg, LinearPSI<Cube> psi, NodeSet x_center, int cell)
        {
            double[] arr = new double[RefElement.SpatialDimension];
            Arg.Diameters.CopyTo(arr, 0);
            psi.SetInactiveDimsToZero(arr);
            MultidimensionalArray diameters = MultidimensionalArray.CreateWrapper(arr, 3, 1);

            NodeSet nodeOnPsi = psi.ProjectOnto(x_center);
            MultidimensionalArray grad = ReferenceGradient(nodeOnPsi, cell);


            MultidimensionalArray jacobian = grid.Jacobian.GetValue_Cell(nodeOnPsi, cell, 1);
            jacobian = jacobian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 });

            LevelSet levelSet = lsData.LevelSet as LevelSet;
            MultidimensionalArray hessian = MultidimensionalArray.Create(1, 1, 3, 3);
            levelSet.EvaluateHessian(cell, 1, nodeOnPsi, hessian);
            hessian = hessian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 }).CloneAs();
            hessian.ApplyAll(x => Math.Abs(x));

            MultidimensionalArray delta1 = diameters.TransposeTo() * jacobian.TransposeTo() * hessian * jacobian * diameters; 


            grad.ApplyAll(x => Math.Abs(x));
            double delta = delta1[0, 0];
            for (int i = 0; i < 3; ++i) {
                delta += grad[i] * diameters[i, 0];
            }
            delta *= boundSafety;
            return delta;
        }

        protected override int FindPromisingHeightDirection(LinearPSI<Cube> Psi, NodeSet Node, int Cell)
        {
            NodeSet nodeOnPsi = Psi.ProjectOnto(Node);
            MultidimensionalArray gradient = ReferenceGradient(nodeOnPsi, Cell);

            int heightDirection = 0;
            double max = double.MinValue;
            for (int i = 0; i < RefElement.SpatialDimension; ++i)
            {
                if ( ( Math.Abs(gradient[i]) > max ) && ( !Psi.DirectionIsFixed(i) ) )
                {
                    max = Math.Abs(gradient[i]);
                    heightDirection = i;
                }
            }
            return heightDirection;
        }
               
        double GlobalSurfaceCurvature(MultidimensionalArray gradient, MultidimensionalArray Hessian) {
            double curvature = 0;
            double norm = gradient.L2Norm();
            for(int i = 0; i < 3; ++i) {
                curvature += Hessian[i, i] * norm;
                curvature -= gradient[i]  / norm * (gradient[0] * Hessian[i, 0] + gradient[1] * Hessian[i, 1] + gradient[2] * Hessian[i, 2]);
            }
            curvature /= (norm * norm);
            curvature *= 0.5;
            return curvature;
        }

        double GlobalLineCurvature(MultidimensionalArray gradient, MultidimensionalArray Hessian, int inactiveDim) {
            double norm = 0;
            for (int i = 0; i < 3; ++i) {
                if (i != inactiveDim) {
                    norm += gradient[i] * gradient[i]; 
                }
            }
            norm = Math.Sqrt(norm);
            
            double curvature = 0;
            for (int i = 0; i < 3; ++i) {
                if(i != inactiveDim) {
                    curvature += Hessian[i, i] * norm;
                    double gradH = 0;
                    for (int j = 0; j < 3; ++j) {
                        if (j != inactiveDim) {
                            gradH += gradient[j] * Hessian[i, j];
                        }
                    }
                    curvature -= gradient[i] / norm * (gradH);
                }
            }
            curvature /= (norm * norm);
            curvature *= 0.5;
            return curvature;
        }

        //double safety = 2;

        protected override bool HeightDirectionIsSuitable(
            LinearSayeSpace<Cube> arg,
            LinearPSI<Cube> psi,
            NodeSet x_center,
            int heightDirection,
            MultidimensionalArray gradient,
            int cell) {

            //throw new NotImplementedException();

            //Determine bounds
            //-----------------------------------------------------------------------------------------------------------------
            NodeSet nodeOnPsi = psi.ProjectOnto(x_center);

            MultidimensionalArray jacobian = grid.Jacobian.GetValue_Cell(nodeOnPsi, cell, 1);
            jacobian = jacobian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 });

            LevelSet levelSet = lsData.LevelSet as LevelSet;
            MultidimensionalArray hessian = MultidimensionalArray.Create(1, 1, 3, 3);
            levelSet.EvaluateHessian(cell, 1, nodeOnPsi, hessian);
            hessian = hessian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 }).CloneAs();

            hessian.ApplyAll(x => Math.Abs(x));

            //abs(Hessian) * 0,5 * diameters.^2 = delta ,( square each entry of diameters) , 
            //this bounds the second error term from taylor series
            //+ + + + 
            double[] arr = new double[] { arg.Diameters[0], arg.Diameters[1], arg.Diameters[2] };
            psi.SetInactiveDimsToZero(arr);
            MultidimensionalArray diameters = MultidimensionalArray.CreateWrapper(arr, 3, 1);
            MultidimensionalArray delta = hessian * jacobian * diameters;

            delta = delta.ExtractSubArrayShallow(-1, 0);

            //Check if suitable
            //-----------------------------------------------------------------------------------------------------------------

            //|gk| > δk
            //Gradient should not be able to change sign
            if (Math.Abs(gradient[heightDirection]) > delta[heightDirection]) {
                bool suitable = true;
                // ||Grad + maxChange|| should be smaller than 20 * 
                // Sum_j( g_j + delta_j)^2 / (g_k - delta_k)^2 < 20
                double sum = 0;

                for (int j = 0; j < delta.Length; ++j) {
                    if (!psi.DirectionIsFixed(j)) {
                        sum += Math.Pow(Math.Abs(gradient[j]) + delta[j], 2);
                    }
                }
                sum /= Math.Pow(Math.Abs(gradient[heightDirection]) - delta[heightDirection], 2);

                suitable &= sum < 20;

                return suitable;
            }
            return true;
        }

        protected override LinearSayeSpace<Cube> Subdivide(LinearSayeSpace<Cube> Arg)
        {
            LinearSayeSpace<Cube> newArg = Arg.Subdivide();
            return newArg;

        }

        protected override bool SubdivideSuitable(int numOfSubdivions)
        {
            if (numOfSubdivions > 10)
                return false;
            return true;
        }

        protected override LinearSayeSpace<Cube> DeriveNewArgument(LinearSayeSpace<Cube> OldSquare)
        {
            return OldSquare.DeriveNew();
        }

        protected MultidimensionalArray ReferenceGradient(NodeSet Node, int Cell)
        {
            MultidimensionalArray gradient = lsData.GetLevelSetGradients(Node, Cell, 1);
            gradient = gradient.ExtractSubArrayShallow(new int[] { 0, 0, -1 }).CloneAs();

            MultidimensionalArray jacobian = grid.Jacobian.GetValue_Cell(Node, Cell, 1);
            jacobian = jacobian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 });
            jacobian.TransposeInPlace();

            double[] tmp_grad = gradient.Storage;
            jacobian.MatVecMulInplace(1, tmp_grad, false);
            
            return gradient;
        }

        #endregion

        #region Evaluate Saye Integrand

        protected override QuadRule CreateZeroQuadrule()
        {
            QuadRule zeroRule = QuadRule.CreateEmpty(RefElement, 1, RefElement.SpatialDimension);
            zeroRule.Nodes.LockForever();
            return zeroRule;
        }

        protected override SayeQuadRule SetLowOrderQuadratureNodes(LinearSayeSpace<Cube> arg)
        {
            Console.WriteLine($"Low order Quadrature required in cell: {cell}, " +
                $"center: {lsData.GridDat.Cells.CellCenter[cell, 0]}, {lsData.GridDat.Cells.CellCenter[cell, 1]}, {lsData.GridDat.Cells.CellCenter[cell, 2]}");
            return SetGaussQuadratureNodes(arg, 1);
        }

        protected override SayeQuadRule SetGaussQuadratureNodes(LinearSayeSpace<Cube> arg)
        {
            return SetGaussQuadratureNodes(arg, order);
        }

        SayeQuadRule SetGaussQuadratureNodes(LinearSayeSpace<Cube> arg, int gaussOrder) {
            //Aquire needed data
            //------------------------------------------------------------------------------------------------------------
            MultidimensionalArray nodes_GaussRule_3D;
            MultidimensionalArray weights_GaussRule_3D;
            MultidimensionalArray centerArr = arg.GetCellCenter().ExtractSubArrayShallow(0, -1);

            double[] diameters = arg.Diameters;
            double jacobianDet;

            //Gaussrule 2d or Gaussrule 3d?
            //2d Gaussrule embedded in 3d space on [-1,1]^3
            //------------------------------------------------------------------------------------------------------------
            if (arg.Dimension == 2) {
                QuadRule gaussRule_2D = Square.Instance.GetQuadratureRule(gaussOrder);
                MultidimensionalArray nodes_GaussRule_2D = gaussRule_2D.Nodes;

                nodes_GaussRule_3D = MultidimensionalArray.Create(gaussRule_2D.NoOfNodes, 3);
                weights_GaussRule_3D = gaussRule_2D.Weights.CloneAs();

                //Embed 2D nodes in 3d space
                for (int i = 0; i < gaussRule_2D.NoOfNodes; ++i) {
                    int nodePositionCounter = 0;
                    for (int j = 0; j < 3; ++j) {
                        if (arg.DimActive(j)) {
                            nodes_GaussRule_3D[i, j] = nodes_GaussRule_2D[i, nodePositionCounter];
                            ++nodePositionCounter;
                        } else {
                            nodes_GaussRule_3D[i, j] = 0;
                        }
                    }
                }
                //Set rest
                jacobianDet = 1;
                for (int j = 0; j < 3; ++j) {
                    if (arg.DimActive(j)) {
                        jacobianDet *= diameters[j];
                    }
                }
            }
            //3d Gauss quadrature rule on [-1,1]^3
            //------------------------------------------------------------------------------------------------------------
            else {
                Debug.Assert(arg.Dimension == 3);
                //Extract Rule
                QuadRule gaussRule_3D = Cube.Instance.GetQuadratureRule(gaussOrder);
                nodes_GaussRule_3D = ((MultidimensionalArray)gaussRule_3D.Nodes).CloneAs();
                weights_GaussRule_3D = gaussRule_3D.Weights.CloneAs();
                //Set rest

                jacobianDet = diameters[0] * diameters[1] * diameters[2];
            }

            //AffineTransformation of nodes, scale weights
            //------------------------------------------------------------------------------------------------------------
            //Scale Nodes
            for (int i = 0; i < 3; ++i) {
                if (arg.DimActive(i)) {
                    nodes_GaussRule_3D.ColScale(i, diameters[i]);
                }
            }
            //Move Nodes
            int[] index = new int[] { 0, -1 };
            for (int i = 0; i < nodes_GaussRule_3D.Lengths[0]; ++i) {
                index[0] = i;
                nodes_GaussRule_3D.AccSubArray(1, centerArr, index);
            }

            //Scale Weights
            weights_GaussRule_3D.Scale(jacobianDet);

            //Set return data
            //------------------------------------------------------------------------------------------------------------
            SayeQuadRule transformed_GaussRule_2D = new SayeQuadRule(nodes_GaussRule_3D, weights_GaussRule_3D, RefElement);
            return transformed_GaussRule_2D;
        }

        protected override void RestrictToBound(MultidimensionalArray X, double bound, int direction)
        {
            X[0, direction] = bound;
        }

        protected override double[] FindRoots(LinearPSI<Cube> psi, MultidimensionalArray X, int heightDirection, double[] bounds, int cell)
        {
            MultidimensionalArray XonPsi = psi.ProjectOnto(X);
            //XonPsi = XonPsi.ExtractSubArrayShallow( 0, -1 );
            //double[] start = XonPsi.To1DArray();
            //double[] end = XonPsi.To1DArray();
            Vector start = XonPsi.GetRowPt(0);
            Vector end = XonPsi.GetRowPt(0);

            start[heightDirection] = -1;
            end[heightDirection] = 1;

            HMF.LineSegment line = new HMF.LineSegment(3, RefElement, start, end);
            LevelSet levelSet = lsData.LevelSet as LevelSet;
            line.ProjectBasisPolynomials(levelSet.Basis);
            double[] roots = rootFinder.GetRoots(line, levelSet, cell, this.iKref);

            return roots;
        }

        protected override SayeQuadRule BuildQuadRule(MultidimensionalArray X, double X_weight, int heightDirection, double length)
        {
            QuadRule gaussRule_1D = Line.Instance.GetQuadratureRule(order);

            double[,] nodesArr = new double[gaussRule_1D.NoOfNodes, 3];
            for (int i = 0; i < gaussRule_1D.NoOfNodes; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    nodesArr[i, j] = X[0, j];
                    if (j == heightDirection)
                    {
                        nodesArr[i, j] += length / 2 * gaussRule_1D.Nodes[i, 0];
                    }
                }
            }

            MultidimensionalArray weights = gaussRule_1D.Weights.CloneAs();
            weights.Scale(length / 2);
            weights.Scale(X_weight);

            MultidimensionalArray nodes = new MultidimensionalArray(2);
            nodes.InitializeFrom(nodesArr);
            SayeQuadRule transformed_GaussRule_1D = new SayeQuadRule(nodes, weights, RefElement);

            return transformed_GaussRule_1D;
        }

        protected override SayeQuadRule BuildSurfaceQuadRule(MultidimensionalArray X, double X_weight, int heightDirection, int cell)
        {
            double weight = X_weight;

            NodeSet node = new NodeSet(RefElement, X.To2DArray(), true);
            MultidimensionalArray gradient = ReferenceGradient(node, cell);
            weight *= gradient.L2Norm() / Math.Abs(gradient[heightDirection]);

            MultidimensionalArray jacobian = grid.Jacobian.GetValue_Cell(node, cell, 1).ExtractSubArrayShallow(0, 0, -1 , -1);
            //Scale weight
            if (IsScalingMatrix(jacobian)) {
                weight /= jacobian[heightDirection, heightDirection];
            } else {
                throw new NotImplementedException("To do");
            }

            MultidimensionalArray weightArr = new MultidimensionalArray(1);
            weightArr.Allocate(1);
            weightArr[0] = weight;
            return new SayeQuadRule(node, weightArr, RefElement);
        }

        protected bool IsScalingMatrix(MultidimensionalArray matrix) {
            double offDiag = 0.0;
            for(int i = 0; i < 3; ++i) {
                for(int j = i + 1; j > 3; ++i) {
                    offDiag += matrix[i, j];
                    offDiag += matrix[j, i];
                }
            }
            if(offDiag == 0.0)
                return true;
            else return false;
        }

        public override double[] GetBoundaries(LinearSayeSpace<Cube> arg, int heightDirection)
        {
            return arg.Boundaries[heightDirection];
        }

        public override NodeSet NodeOnRay(MultidimensionalArray X, int direction, double distance)
        {
            MultidimensionalArray nodeArr = X.CloneAs();
            nodeArr[0, direction] += distance;
            NodeSet node = new NodeSet(RefElement, nodeArr, true);
            return node;
        }

        #endregion

    }
}
