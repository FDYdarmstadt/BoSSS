using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
            base(_lsData, RootFinder)
        {
        }

        public override SayeSquare StartSetup()
        {
            LinearPSI<Square> psi = new LinearPSI<Square>(Square.Instance);
            Tuple<LinearPSI<Square>, int> psi_i_s_i = new Tuple<LinearPSI<Square>, int>(psi, 1);
            SayeSquare arg = new SayeSquare(psi_i_s_i, false);
            arg.Reset();
            return arg;
        }
    }

    /// <summary>
    /// Gauss rules for \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ in the 2D case
    /// </summary>
    class SayeGaussRule_LevelSet2D :
        SayeFactory_Square
    {
        public SayeGaussRule_LevelSet2D(LevelSetTracker.LevelSetData _lsData, IRootFindingAlgorithm RootFinder) :
            base(_lsData, RootFinder)
        {
        }

        public override SayeSquare StartSetup()
        {
            throw new NotImplementedException();
        }
    }

    public abstract class SayeFactory_Square :
        SayeIntegrand<LinearPSI<Square>, SayeSquare>,
        IQuadRuleFactory<QuadRule>
    {
        LevelSetTracker.LevelSetData lsData;

        IRootFindingAlgorithm rootFinder;

        public abstract SayeSquare StartSetup();

        int order;

        public SayeFactory_Square(
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder)
        {
            this.lsData = _lsData;
            this.rootFinder = RootFinder;
        }

        #region IQaudRuleFactory<QuadRule>

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int Order)
        {
            order = Order;
            QuadRule gaussRule_1D = Line.Instance.GetQuadratureRule(Order);
            var result = new List<ChunkRulePair<QuadRule>>();

            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    SayeSquare arg = StartSetup();
                    QuadRule sayeRule = this.Evaluate(cell, arg);
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(chunk, sayeRule);
                    result.Add(sayePair);
                }
            }
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

        protected override int EvaluateSign()
        {
            //WIP, Not Finished!
            return 1;
        }

        protected override LinearPSI<Square>[] ExtractSubPsis(LinearPSI<Square> psi, SayeSquare arg, int heightDirection)
        {
            double[] bounds = arg.GetBoundaries(heightDirection);
            double x_U = bounds[0];
            double x_L = bounds[1];
            LinearPSI<Square> Psi_U = psi.ReduceDim(heightDirection, x_U);
            LinearPSI<Square> Psi_L = psi.ReduceDim(heightDirection, x_L);
            return new LinearPSI<Square>[] { Psi_U, Psi_L };
        }

        protected override MultidimensionalArray Gradient(LinearPSI<Square> psi, NodeSet Node, int Cell)
        {
            MultidimensionalArray gradient = lsData.GetLevelSetReferenceGradients(Node, Cell, 1);
            gradient = gradient.ExtractSubArrayShallow(new int[] { 0,0,-1});
            return gradient;
        }

        protected override double EvaluateAt(LinearPSI<Square> Psi, NodeSet Node, int cell)
        {
#if DEBUG
            if (Node.NoOfNodes != 1)
            {
                throw new NotSupportedException();
            }
#endif
            NodeSet nodeOnPsi = Psi.ProjectOnto(Node);
            double value = lsData.GetLevSetValues(Node, cell, 1)[0,0];
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
            double[] arr = Arg.Diameters;
            psi.SetInactiveDimsToZero(arr);
            MultidimensionalArray diameters = MultidimensionalArray.CreateWrapper(arr, new int[] { 2});

            NodeSet nodeOnPsi = psi.ProjectOnto(x_center);
            MultidimensionalArray grad = lsData.GetLevelSetReferenceGradients(nodeOnPsi, cell, 1);
            grad.ApplyAll(x => Math.Abs(x));

            grad = grad.ExtractSubArrayShallow(new int[] { 0, 0, -1 });

            double delta = grad.InnerProduct(diameters);

            return delta;
        }

        protected override int FindPromisingHeightDirection(LinearPSI<Square> Psi, NodeSet Node, int Cell)
        {
            MultidimensionalArray gradient = lsData.GetLevelSetReferenceGradients(Node, Cell, 1);
            if (Math.Abs(gradient[0, 0, 0]) > Math.Abs(gradient[0, 0, 1]))
                return 0;
            return 1;
        }

        protected override bool HeightDirectionIsSuitable(SayeSquare arg, LinearPSI<Square> psi, NodeSet x_center, int heightDirection,
            MultidimensionalArray gradient, int cell)
        {
            //Determine bounds
            //-----------------------------------------------------------------------------------------------------------------
            double[] arr = arg.Diameters;
            psi.SetInactiveDimsToZero(arr);
            MultidimensionalArray diameters = MultidimensionalArray.CreateWrapper(arr, new int[] {2, 1 });

            NodeSet nodeOnPsi = psi.ProjectOnto(x_center);
            MultidimensionalArray hessian = lsData.GetLevelSetReferenceHessian(nodeOnPsi, cell, 1);
            hessian = hessian.ExtractSubArrayShallow(new int[] { 0, 0, -1, -1});
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

        #region Evaluate Saye Integrand

        protected override QuadRule SetLowOrderQuadratureNodes(SayeSquare arg)
        {
            throw new NotImplementedException();
        }

        protected override QuadRule SetGaussQuadratureNodes(SayeSquare arg)
        {
            QuadRule gaussRule_2D = Square.Instance.GetQuadratureRule(order);
            double[] diameters = arg.Diameters;

            /*
            double[] center = arg.GetCellCenter().To1DArray();

            
            //AffineTransformation
            gaussRule_2D.Nodes.ColScale();
            gaussRule_2D.Nodes.RowScale();
            
            //Scale Weights
            
            gaussRule_2D.Nodes.AccVector(1,center);
            */

            return gaussRule_2D;
        }

        protected override double FindRoot(LinearPSI<Square> psi, MultidimensionalArray X, double[] bounds, int cell)
        {
            throw new NotImplementedException();
        }

        protected override QuadRule BuildQuadRule()
        {
            throw new NotImplementedException();
        }

        public override double[] GetBoundaries(SayeSquare arg, int heightDirection)
        {
            throw new NotImplementedException();
        }

        public override NodeSet NodeOnRay(MultidimensionalArray X, int direction, double distance)
        {
            throw new NotImplementedException();
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

        public SayeSquare()
        {
            dim = 2;
            StandardSetup();
        }

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
            sibling.Surface = this.surface;
            this.RecalculateCenter();
            sibling.RecalculateCenter();

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

        NodesAndWeightsLinkedList data = new NodesAndWeightsLinkedList(refElement.SpatialDimension);

        public void Reset()
        {
            data.Reset();
        }

        public override INodesAndWeights NodesAndWeights 
        {
            get => data;
        }

        public override int HeightDirection 
        {
            get;
            set;
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

