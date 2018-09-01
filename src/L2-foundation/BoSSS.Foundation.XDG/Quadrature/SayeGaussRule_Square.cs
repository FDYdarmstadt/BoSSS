
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;
using ilPSP;
using System.Collections;


namespace BoSSS.Foundation.XDG.Quadrature
{

    /// <summary>
    /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
    /// </summary>
    public class SayeGaussRule_Volume2D :
        SayeFactory_Square
    {
        public SayeGaussRule_Volume2D(LevelSetTracker.LevelSetData _lsData, IRootFindingAlgorithm RootFinder) :
            base(_lsData, RootFinder, QuadratureMode.Volume)
        {
        }
    }

    /// <summary>
    /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 2D case
    /// </summary>
    public class SayeGaussRule_LevelSet2D :
        SayeFactory_Square
    {
        public SayeGaussRule_LevelSet2D(LevelSetTracker.LevelSetData _lsData, IRootFindingAlgorithm RootFinder) :
            base(_lsData, RootFinder, QuadratureMode.Surface)
        {
        }
    }

    public class SayeFactory_Square :
        SayeIntegrand<LinearPSI<Square>, SayeSquare>,
        IQuadRuleFactory<QuadRule>
    {
        LevelSetTracker.LevelSetData lsData;

        IRootFindingAlgorithm rootFinder;

        int order;

        IGridData grid;

        int iKref;

        public enum QuadratureMode { Surface, Volume };

        QuadratureMode mode; 

        public SayeFactory_Square(
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder,
            QuadratureMode Mode)
        {
            this.lsData = _lsData;
            this.rootFinder = RootFinder;
            mode = Mode;
        }

        SayeSquare CreateStartSetup()
        {
            bool IsSurfaceIntegral = (mode == QuadratureMode.Surface);

            LinearPSI<Square> psi = new LinearPSI<Square>(Square.Instance);
            Tuple<LinearPSI<Square>, int> psi_i_s_i = new Tuple<LinearPSI<Square>, int>(psi, 1);
            SayeSquare arg = new SayeSquare(psi_i_s_i, IsSurfaceIntegral);
            arg.Reset();
            return arg;
        }

        #region IQaudRuleFactory<QuadRule>

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int Order)
        {
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

            order = Order;
            QuadRule gaussRule_1D = Line.Instance.GetQuadratureRule(Order);
            var result = new List<ChunkRulePair<QuadRule>>();
            grid = mask.GridData;
            iKref = mask.GridData.iGeomCells.RefElements.IndexOf(this.RefElement, (A, B) => object.ReferenceEquals(A, B));

            int number = 0;
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    SayeSquare arg = CreateStartSetup();
                    QuadRule sayeRule = this.Evaluate(cell, arg);
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), sayeRule);
                    result.Add(sayePair);
                    ++number;
                }
            }
            stopWatch.Stop();
            long ts = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Number Of Cells " + number);
            Console.WriteLine("RunTime " + ts + "ms");
            return result;
        }

        public RefElement RefElement {
            get {
                return Square.Instance;
            }
        }

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region Build SayeIntegrand

        protected override int EvaluateSign(MultidimensionalArray gradient, int heightDirection, int s_i, bool surface, int sigma)
        {
            int m = Math.Sign(gradient[heightDirection]);
            if(surface == true || sigma * s_i == m )
            {
                return sigma * m;
            }
            else
            {
                return 0;
            }
        }

        protected override LinearPSI<Square>[] ExtractSubPsis(LinearPSI<Square> psi, SayeSquare arg, int heightDirection)
        {
            double[] bounds = arg.GetBoundaries(heightDirection);
            double x_U = bounds[1];
            double x_L = bounds[0];
            LinearPSI<Square> Psi_U = psi.ReduceDim(heightDirection, x_U);
            LinearPSI<Square> Psi_L = psi.ReduceDim(heightDirection, x_L);
            return new LinearPSI<Square>[] { Psi_L, Psi_U };
        }

        protected override MultidimensionalArray Gradient(LinearPSI<Square> psi, NodeSet Node, int Cell)
        {
            //Move Node onto psi 
            MultidimensionalArray gradient = ScaledReferenceGradient(Node, Cell);

            return gradient;
        }

        protected override double EvaluateAt(LinearPSI<Square> Psi, NodeSet Node, int cell)
        {
            Debug.Assert(Node.NoOfNodes == 1);
            
            NodeSet nodeOnPsi = Psi.ProjectOnto(Node);
            double value = lsData.GetLevSetValues(nodeOnPsi, cell, 1)[0,0];
            return value;
        }

        /// <summary>
        /// First order approximation of  delta >= sup_x|psi(x) - psi(x_center)|  
        /// </summary>
        /// <param name="Arg"></param>
        /// <param name="psi"></param>
        /// <param name="x_center"></param>
        /// <param name="cell"></param>
        /// <returns></returns>
        protected override double EvaluateBounds(SayeSquare Arg, LinearPSI<Square> psi, NodeSet x_center, int cell)
        {
            double[] arr = new double[Arg.Dimension];
            Arg.Diameters.CopyTo(arr, 0);
            psi.SetInactiveDimsToZero(arr);
            MultidimensionalArray diameters = MultidimensionalArray.CreateWrapper(arr, new int[] { 2});

            NodeSet nodeOnPsi = psi.ProjectOnto(x_center);


            MultidimensionalArray grad = ScaledReferenceGradient(nodeOnPsi, cell);
            
            grad.ApplyAll(x => Math.Abs(x));     
            double delta = grad.InnerProduct(diameters);

            return delta;
        }

        protected override int FindPromisingHeightDirection(LinearPSI<Square> Psi, NodeSet Node, int Cell)
        {
            MultidimensionalArray gradient = ScaledReferenceGradient(Node, Cell);
            if (Math.Abs(gradient[0]) > Math.Abs(gradient[1]))
                return 0;
            return 1;
        }

        protected override bool HeightDirectionIsSuitable(SayeSquare arg, LinearPSI<Square> psi, NodeSet x_center, 
            int heightDirection, MultidimensionalArray gradient, int cell)
        {
            //Determine bounds
            //-----------------------------------------------------------------------------------------------------------------
            double[] arr = arg.Diameters;
            psi.SetInactiveDimsToZero(arr);
            MultidimensionalArray diameters = MultidimensionalArray.CreateWrapper(arr, new int[] {2, 1 });

            NodeSet nodeOnPsi = psi.ProjectOnto(x_center);
            MultidimensionalArray hessian = lsData.GetLevelSetReferenceHessian(nodeOnPsi, cell, 1);
            hessian = hessian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1});

            MultidimensionalArray jacobian = grid.InverseJacobian.GetValue_Cell(nodeOnPsi, cell, 1);
            jacobian = jacobian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 });
            hessian = jacobian * hessian; 

            hessian.ApplyAll(x => Math.Abs(x));

            //abs(Hessian) * diameters = delta 
            MultidimensionalArray delta = hessian * diameters;
            delta = delta.ExtractSubArrayShallow(new int[] { -1, 0 });
            
            //Check if suitable
            //-----------------------------------------------------------------------------------------------------------------

            //|gk| > δk
            if( Math.Abs(gradient[heightDirection]) > delta[heightDirection])
            {
                bool suitable = true;
                // Sum_j( g_j + delta_j)^2 / (g_k - delta_k)^2 < 20
                double sum = 0;

                for (int j = 0; j < delta.Length; ++j)
                {
                    sum += Math.Pow(gradient[j] + delta[j], 2);
                }
                sum /= Math.Pow(gradient[heightDirection] - delta[heightDirection], 2);

                suitable &= sum < 20;

                return suitable;
            }
            return false;
        }

        protected override SayeSquare Subdivide(SayeSquare Arg)
        {
            SayeSquare newArg = Arg.Subdivide();
            return newArg;
        }

        protected override bool SubdivideSuitable(SayeSquare arg)
        {
            return true;
        }

        protected override SayeSquare DeriveNewArgument(SayeSquare OldSquare)
        {
            return OldSquare.DeriveNew();
        }

        #endregion

        private MultidimensionalArray ScaledReferenceGradient(NodeSet Node, int Cell)
        {
            MultidimensionalArray gradient = lsData.GetLevelSetReferenceGradients(Node, Cell, 1);
            gradient = gradient.ExtractSubArrayShallow(new int[] { 0, 0, -1 });

            MultidimensionalArray jacobian = grid.InverseJacobian.GetValue_Cell(Node, Cell, 1);
            jacobian = jacobian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1 });

            double[] tmp_grad = gradient.Storage;
            jacobian.gemv(1, tmp_grad, 0, tmp_grad);

            return gradient;
        }

        #region Evaluate Saye Integrand

        protected override SayeQuadRule SetLowOrderQuadratureNodes(SayeSquare arg)
        {
            throw new NotImplementedException();
        }

        protected override SayeQuadRule SetGaussQuadratureNodes(SayeSquare arg)
        {
            //Aquire needed data
            //------------------------------------------------------------------------------------------------------------
            QuadRule gaussRule_2D = Square.Instance.GetQuadratureRule(order);
            MultidimensionalArray nodes_GaussRule_2D = ((MultidimensionalArray)gaussRule_2D.Nodes).CloneAs();
            MultidimensionalArray weights_GaussRule_2D = gaussRule_2D.Weights.CloneAs();

            double[] diameters = arg.Diameters;
            MultidimensionalArray centerArr = arg.GetCellCenter().ExtractSubArrayShallow(new int[] { 0, -1 });
            double jacobian = diameters[0] * diameters[1];

            //AffineTransformation of nodes, scale weights
            //------------------------------------------------------------------------------------------------------------
            //Scale Nodes
            for(int i = 0; i < 2; ++i){
                nodes_GaussRule_2D.ColScale(i, diameters[i]);
            }
            //Scale Weights
            weights_GaussRule_2D.Scale(jacobian);
            //Move Nodes
            int[] index = new int[] { 0, -1 };
            for (int i = 0; i < gaussRule_2D.NoOfNodes; ++i)
            {
                index[0] = i;
                nodes_GaussRule_2D.AccSubArray(1, centerArr, index);
            }

            //Set return data
            //------------------------------------------------------------------------------------------------------------
            SayeQuadRule transformed_GaussRule_2D = new SayeQuadRule(nodes_GaussRule_2D, weights_GaussRule_2D);
            return transformed_GaussRule_2D;
        }

        protected override void RestrictToBound(MultidimensionalArray X, double bound, int direction)
        {
            X[0, direction] = bound;
        }

        protected override double[] FindRoots(LinearPSI<Square> psi, MultidimensionalArray X, int heightDirection, double[] bounds, int cell)
        {
            MultidimensionalArray XonPsi = psi.ProjectOnto(X);
            XonPsi = XonPsi.ExtractSubArrayShallow(new int[] { 0, -1 });
            double[] start = XonPsi.To1DArray();
            double[] end = XonPsi.To1DArray();

            start[heightDirection] = -1;
            end[heightDirection] = 1;

            HMF.LineSegment line = new HMF.LineSegment(2, RefElement, start, end);
            LevelSet levelSet = lsData.LevelSet as LevelSet;
            line.ProjectBasisPolynomials(levelSet.Basis);
            double[] roots = rootFinder.GetRoots(line, levelSet, cell,  this.iKref);

            return roots;
        }

        protected override SayeQuadRule BuildQuadRule(MultidimensionalArray X, double X_weight, int heightDirection, double length)
        {
            QuadRule gaussRule_1D = Line.Instance.GetQuadratureRule(order);

            double[,] nodesArr = new double[gaussRule_1D.NoOfNodes, 2];
            for (int i = 0; i < gaussRule_1D.NoOfNodes; ++i)
            {
                for(int j = 0; j < 2; ++j)
                {
                    nodesArr[i, j] = X[0,j];
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
            SayeQuadRule transformed_GaussRule_1D = new SayeQuadRule(nodes, weights);
            
            return transformed_GaussRule_1D;
        }

        protected override SayeQuadRule BuildSurfaceQuadRule(MultidimensionalArray X, double X_weight, int heightDirection, int cell)
        {
            double weight = X_weight;

            NodeSet node = new NodeSet(RefElement, X.To2DArray());
            MultidimensionalArray gradient = lsData.GetLevelSetGradients(node, cell, 1);
            gradient = gradient.ExtractSubArrayShallow(new int[] { 0, 0, -1 });
            
            //Scale weight
            weight *= gradient.L2Norm() / Math.Abs(gradient[heightDirection]);

            MultidimensionalArray weightArr = new MultidimensionalArray(1);
            weightArr.Allocate(1);
            weightArr[0] = weight;
            return new SayeQuadRule(node, weightArr);
        }

        public override double[] GetBoundaries(SayeSquare arg, int heightDirection)
        {
            return arg.Boundaries[heightDirection];
        }

        public override NodeSet NodeOnRay(MultidimensionalArray X, int direction, double distance)
        {
            MultidimensionalArray nodeArr = X.CloneAs();
            nodeArr[0, direction] += distance;
            NodeSet node = new NodeSet(RefElement, nodeArr);
            return node;
        }

        #endregion
    }

    public class SayeSquare :
        SayeArgument<LinearPSI<Square>>
    {
        static Square refElement = Square.Instance;

        BitArray removedDims = new BitArray(2);
        protected bool subdivided = false;

        public List<Tuple<LinearPSI<Square>, int>> psiAndS = new List<Tuple<LinearPSI<Square>, int>>();
        public enum Dim { x, y};

        public SayeSquare(Tuple<LinearPSI<Square>, int> PsiAndS, bool _Surface)
        {
            dim = 2;
            surface = _Surface;
            StandardSetup();
            psiAndS.Add(PsiAndS);
        }

        private SayeSquare(double[][] _Boundaries)
        {
            dim = 2;
            SetBoundaries(_Boundaries);
        }

        private void StandardSetup()
        {
            boundaries = new double[2][];
            boundaries[0] = new double[] { -1, 1 };
            boundaries[1] = new double[] { -1, 1 };
            diameters = new double[] { 1, 1 };
        }

        public SayeSquare DeriveNew()
        {
            double[][] newBoundary = boundaries.Copy();
            SayeSquare arg = new SayeSquare(newBoundary);
            arg.subdivided = this.subdivided;
            arg.Surface = this.surface;
            return arg;
        }

        public SayeSquare Subdivide()
        {
            subdivided = true;
            double[][] newBoundary = boundaries.Copy();
                
            //Figure out new Boundaries
            int k;
            if (diameters[0] > diameters[1])
            {
                k = 0;
            }
            else
            {
                k = 1;
            }
            double[] maxBounds = boundaries[k];
            double newBound = (maxBounds[0] + maxBounds[1]) / 2;
            maxBounds[1] = newBound;
            SetBoundaries(boundaries);

            newBoundary[k][0] = newBound;

            SayeSquare sibling = new SayeSquare(newBoundary);
            sibling.psiAndS = new List<Tuple<LinearPSI<Square>, int>>( this.psiAndS);
            sibling.removedDims = this.removedDims.CloneAs();
            sibling.subdivided = true;
            sibling.Surface = this.surface;

            return sibling;
        }

        private void RecalculateCenter()
        {
            double[] centerArr = new double[2];
            centerArr[0] = ( boundaries[0][0] + boundaries[0][1] ) / 2.0;
            centerArr[1] = ( boundaries[1][0] + boundaries[1][1] ) / 2.0;
            this.center = new NodeSet(refElement, centerArr);
            this.center.LockForever();
        }

        double[][] boundaries;

        private void SetBoundaries(double[][] _Boundaries)
        {
            boundaries = _Boundaries;
            diameters = new double[2];
            diameters[0] = (_Boundaries[0][1] - _Boundaries[0][0]) / 2;
            diameters[1] = (_Boundaries[1][1] - _Boundaries[1][0]) / 2;
            RecalculateCenter();
        }

        public double[][] Boundaries
        {
            get => boundaries;       
        }

        public double[] GetBoundaries(int dimension)
        {
            return boundaries[dimension];
        }

        double[] diameters;

        public double[] Diameters
        {
            get => diameters;
        }

        #region ISayeArgument

        NodesAndWeightsLinkedList data = new NodesAndWeightsLinkedList(refElement.SpatialDimension, refElement);

        public void Reset()
        {
            data.Reset();
        }

        public override ISayeQuadRule NodesAndWeights 
        {
            get => data;
        }

        int heightDirection;
        public override int HeightDirection 
        {
            get 
            {
                return heightDirection;
            }
            set 
            {
                heightDirection = value;
            }
        }

        bool surface;

        public override bool Surface 
        {
            get => surface;
            set => surface = value;
        }

        public override void RemoveDimension(int k)
        {
            --dim;
            removedDims[k] = true;
            if (Dimension == 1)
            {
                heightDirection = 0;
                while (removedDims[heightDirection])
                {
                    ++heightDirection;
                }
            }
        }

        NodeSet center;

        public override NodeSet GetCellCenter()
        {
            NodeSet _center;
            if (!subdivided)
            {
                _center = refElement.Center;
            }
            else
            {
                _center = this.center;
            }
            
            return _center;
        }

        int dim;

        public override int Dimension => dim;
        
        public override IList<Tuple<LinearPSI<Square>, int>> PsiAndS 
        {
            get => psiAndS;
        }

        public override int n => psiAndS.Count;

        #endregion
    }
}

