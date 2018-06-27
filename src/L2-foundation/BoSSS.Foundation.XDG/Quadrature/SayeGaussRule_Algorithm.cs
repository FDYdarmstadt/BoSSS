using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using System.Collections;
using System.Diagnostics;


namespace BoSSS.Foundation.XDG.Quadrature
{
    public interface INodesAndWeights
    {
        IEnumerable<MultidimensionalArray> IntegrationNodes 
        {
            get;
        }
        void AddRule(QuadRule rule);
        void ExpandRule(QuadRule rule);
        QuadRule GetQuadRule();
    }

    public interface IPsi
    {
       
    }

    public abstract class SayeArgument<S>
        where S : IPsi
    {
        public enum Mode
        {
            Standart,
            GaussQuadrature,
            LowOrderQuadrature,
            DomainIsEmpty
        };

        public Mode Status = Mode.Standart;

        public abstract void RemoveDimension(int k);

        public abstract bool Surface {
            get;
            set;
        }
        public abstract int Dimension 
        {
            get;
        }
        public abstract IList<Tuple<S, int>> PsiAndS 
        {
            get;
        }
        public abstract int n 
        {
            get;
        }
        public abstract NodeSet GetCellCenter();
        public abstract int HeightDirection {
            get;
            set;
        }
        public abstract INodesAndWeights NodesAndWeights {
            get;
        }
    }

    public abstract class SayeIntegrand<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        int cell;
        public SayeIntegrand()
        {
        }

        public QuadRule Evaluate(int Cell, T arg)
        {
            //Setup algorithm
            //----------------------------------------------------------------------
            cell = Cell;
            TreeNode<T> recursionTree = new TreeNode<T>();
            TreeNode<T> fullSpace = recursionTree.AddChild(arg);

            //Build Integrand
            //----------------------------------------------------------------------
            SayeRecursion(fullSpace);

            //Evaluate Integrand
            //----------------------------------------------------------------------
            //Fill nodesAndWeights
            recursionTree.UnrollFunc(IntegrandEvaluation);
            //ConvertToQuadRule
            QuadRule toRuleThemAll = fullSpace.Value.NodesAndWeights.GetQuadRule();

            return toRuleThemAll;
        }

        #region Evaluate Integrand
        
        private void IntegrandEvaluation(TreeNode<T> node)
        {
            T nodeArg = node.Value;

            //Check a bunch of stuff
            if(nodeArg == null)
            {
                return;
            }
            QuadRule newRule;
            switch (nodeArg.Status)
            { 
                case SayeArgument<S>.Mode.Standart:
                    VolumeIntegrand(node);
                    break;
                case SayeArgument<S>.Mode.GaussQuadrature:
                    newRule = SetGaussQuadratureNodes(nodeArg);
                    nodeArg.NodesAndWeights.AddRule(newRule);
                    break;
                case SayeArgument<S>.Mode.LowOrderQuadrature:
                    newRule = SetLowOrderQuadratureNodes(nodeArg);
                    nodeArg.NodesAndWeights.AddRule(newRule);
                    break;
                case SayeArgument<S>.Mode.DomainIsEmpty:
                    throw new NotSupportedException();
                default:
                    throw new NotSupportedException();
            }     
        }
        
        private void VolumeIntegrand(TreeNode<T> node)
        {
            T arg = node.Value;
            //Rekursionsanker
            if (arg.Dimension == 1)
            {
                MultidimensionalArray newNode = arg.GetCellCenter();
                VolumeIntegrand(newNode, arg, true);
            }
            else
            {
                foreach (TreeNode<T> childNode in node.Children)
                {
                    T childArg = childNode.Value;
                    foreach(MultidimensionalArray integrationNode in childArg.NodesAndWeights.IntegrationNodes)
                    {
                        VolumeIntegrand(integrationNode, arg, false);
                    }
                }
            }
        }

        //Algorithm 1
        private void VolumeIntegrand(MultidimensionalArray X, T arg, bool isNew)
        {
            int heightDirection = arg.HeightDirection;
            double[] roots = new double[arg.n + 2];
            double[] bounds = GetBoundaries(arg, heightDirection);
            roots[0] = bounds[0];
            roots[1] = bounds[1];
            for (int i = 0; i < arg.n; ++i)
            {
                S psi_i = arg.PsiAndS[i].Item1;
                roots[i] = FindRoot(psi_i, X, bounds,this.cell);
            }
            Array.Sort(roots);
            for (int j = 0; (j < roots.Length - 1) && (roots[j + 1] < double.MaxValue); ++j)
            {
                double L = roots[j + 1] - roots[j];
                NodeSet x_c = NodeOnRay(X, heightDirection, roots[j] + L / 2.0);

                bool updateIntegrand = true;

                foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
                {
                    S psi_i = psiAndS.Item1;
                    int s_i = psiAndS.Item2;
                    updateIntegrand &= s_i * EvaluateAt(psi_i, x_c) >= 0;
                }
                if (updateIntegrand)
                {
                    QuadRule newRule = BuildQuadRule();
                    if (isNew)
                    {
                        arg.NodesAndWeights.AddRule(newRule);
                        isNew = false;
                    }
                    else
                    {
                        arg.NodesAndWeights.ExpandRule(newRule);
                    }
                        
                }    
            }
        }

        //Algorithm 2
        private void SurfaceIntegrand(int arg)
        {

        }

        #endregion

        #region Evaluate Integrand: Abstract Functions 

        protected abstract QuadRule SetLowOrderQuadratureNodes(T arg);

        protected abstract QuadRule SetGaussQuadratureNodes(T arg);

        //Return 1 Root. If there isn't a root, return double.MaxValue 
        protected abstract double FindRoot(S psi, MultidimensionalArray X, double[] bounds, int cell);

        protected abstract QuadRule BuildQuadRule();

        public abstract double[] GetBoundaries(T arg, int heightDirection);

        public abstract NodeSet NodeOnRay(MultidimensionalArray X, int direction, double distance);

        #endregion

        #region BuildIntegrand

        //Algorithm 3
        private void SayeRecursion(TreeNode<T> treeNode)
        {
            T arg = treeNode.Value;
            Debug.Assert(arg.Surface == false || (arg.Dimension > 1 && arg.n == 1));
            //Check treeNode : Prune
            //-----------------------------------------------------------------------------------------------
            if (arg.Dimension == 1)
            {
                return;
            }
            NodeSet x_center = arg.GetCellCenter();

            for (int i = arg.n - 1 ; i >= 0;  --i)
            {
                Tuple<S, int> psiAndS = arg.PsiAndS[i];
                S psi_i = psiAndS.Item1;
                double delta = EvaluateBounds(arg, psi_i, x_center);
                if (delta <= EvaluateAt(psi_i, x_center))
                {
                    int s_i = psiAndS.Item2;
                    if (s_i * EvaluateAt(psi_i, x_center) >= 0)
                    {
                        //Remove ψi from the list and decrement n by one.
                        arg.PsiAndS.RemoveAt(i);
                    }
                    else
                    {
                        //The domain of integration is empty; return 0.
                        arg.Status = SayeArgument<S>.Mode.DomainIsEmpty;
                        return;
                    }
                }
            }

            if (arg.n == 0)
            {
                //Apply a tensor - product Gaussian quadrature scheme
                arg.Status = SayeArgument<S>.Mode.GaussQuadrature;
                return;
            }

            //Find new subspace
            //-----------------------------------------------------------------------------------------------
            S psi_1 = arg.PsiAndS[0].Item1;
            int k = FindPromisingHeightDirection(psi_1, x_center);
            T subspaceArg = DeriveNewArgument(arg);
            foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
            {
                S psi_i = psiAndS.Item1;
                int s_i = psiAndS.Item2;
                MultidimensionalArray g = Gradient(psi_i, x_center);
                if (HeightDirectionIsSuitable(arg, psi_i, x_center, k, g))
                {
                    S[] subPsis = ExtractSubPsis(psi_i, arg, k);
                    S psi_U = subPsis[0];
                    S psi_L = subPsis[1];
                    int s_U = EvaluateSign();
                    int s_L = EvaluateSign();
                    Tuple<S, int> newPsiAndS_U = new Tuple<S, int>(psi_U, s_U);
                    Tuple<S, int> newPsiAndS_L = new Tuple<S, int>(psi_L, s_L);
                    subspaceArg.PsiAndS.Add(newPsiAndS_U);
                    subspaceArg.PsiAndS.Add(newPsiAndS_L);
                }
                else
                {
                    //The height function direction ek is not suitable for ψi. If already subdivided too
                    //many times, revert to a low - order method(see discussion).Otherwise split U
                    if (SubdivideSuitable(arg))
                    {
                        //Subdivide
                        T siblingArg = Subdivide(arg);
                        TreeNode<T>sibling = treeNode.AddSibling(siblingArg);
                        //Recalculate
                        SayeRecursion(treeNode);
                        SayeRecursion(sibling);
                    }
                    else
                    {
                        arg.Status = SayeArgument<S>.Mode.LowOrderQuadrature;
                    }
                    return;
                }
            }
            subspaceArg.Surface = false;
            arg.HeightDirection = k;
            subspaceArg.RemoveDimension(k);
            TreeNode<T> subSpaceNode = treeNode.AddChild(subspaceArg);
            SayeRecursion(subSpaceNode);
        }

        #endregion

        #region BuildIntegrand: Wrappers for Readability

        double EvaluateAt(S Psi, NodeSet Point)
        {
            return EvaluateAt(Psi, Point, this.cell);
        }

        int FindPromisingHeightDirection(S psi, NodeSet Point)
        {
            return FindPromisingHeightDirection(psi, Point, this.cell);
        }

        bool HeightDirectionIsSuitable(T arg, S psi, NodeSet Point, int heightDirection, MultidimensionalArray gradient)
        {
            return HeightDirectionIsSuitable(arg, psi, Point, heightDirection, gradient,  this.cell);
        }

        MultidimensionalArray Gradient(S psi, NodeSet Node)
        {
            return Gradient(psi, Node, this.cell);
        }
        double EvaluateBounds(T arg, S psi, NodeSet x_center)
        {
            return EvaluateBounds(arg, psi, x_center, this.cell);
        }

        #endregion

        #region BuildIntegrand: Abstract Functions

        protected abstract MultidimensionalArray Gradient(S psi, NodeSet Node, int Cell);

        protected abstract double EvaluateAt(S Psi, NodeSet Point, int cell);

        protected abstract int EvaluateSign();

        protected abstract S[] ExtractSubPsis(S psi_i, T arg, int heightDirection);

        protected abstract bool SubdivideSuitable(T arg);

        protected abstract T Subdivide(T arg);

        protected abstract int FindPromisingHeightDirection(S psi, NodeSet Point, int cell);

        protected abstract bool HeightDirectionIsSuitable(T arg, S psi, NodeSet Point, int heightDirection, 
            MultidimensionalArray gradient, int cell);

        protected abstract double EvaluateBounds(T arg, S psi, NodeSet x_center, int cell);

        protected abstract T DeriveNewArgument(T arg);

        #endregion
        
    }
}
