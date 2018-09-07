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
    public interface ISayeQuadRule
    {
        IEnumerable<Tuple<MultidimensionalArray,double>> IntegrationNodes 
        {
            get;
        }
        void AddRule(SayeQuadRule rule, bool deriveFromExistingNode);
        void RemoveActiveNode();
        QuadRule GetQuadRule();
    }

    public interface IPsi
    {
       
    }

    interface ISortedBoundedList<T>
        : IEnumerable<T>
        where T : IComparable<T>
    {
        void SetBounds(T min, T max);
        void Add(params T[] arr);
        T this[int index] 
        {
            get;
        }
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
        public abstract ISayeQuadRule NodesAndWeights {
            get;
        }
    }

    public abstract class SayeIntegrand<S, T>
        where S : IPsi
        where T : SayeArgument<S>
    {
        protected int cell;
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

        protected void IntegrandEvaluation(TreeNode<T> node)
        {
            T nodeArg = node.Value;

            //Check a bunch of stuff
            if (nodeArg == null)
            {
                return;
            }
            SayeQuadRule newRule;
            switch (nodeArg.Status)
            {
                case SayeArgument<S>.Mode.Standart:
                    SetStandartNodes(node);
                    break;
                case SayeArgument<S>.Mode.GaussQuadrature:
                    newRule = SetGaussQuadratureNodes(nodeArg);
                    nodeArg.NodesAndWeights.AddRule(newRule, false);
                    break;
                case SayeArgument<S>.Mode.LowOrderQuadrature:
                    newRule = SetLowOrderQuadratureNodes(nodeArg);
                    nodeArg.NodesAndWeights.AddRule(newRule, false);
                    break;
                case SayeArgument<S>.Mode.DomainIsEmpty:
                    break;
                default:
                    throw new NotSupportedException();
            }
        }

        private void SetStandartNodes(TreeNode<T> node)
        {
            T arg = node.Value;
            //Rekursionsanker
            if (arg.Dimension == 1)
            {
                MultidimensionalArray newNode = ((MultidimensionalArray)arg.GetCellCenter()).CloneAs();
                VolumeIntegrand(newNode, 1.0, arg, true);
            }
            else
            {
                Action<MultidimensionalArray, double, T, bool> integrand;
                if (arg.Surface)
                {
                    integrand = SurfaceIntegrand;
                }
                else
                {
                    integrand = VolumeIntegrand;
                }
                foreach (TreeNode<T> childNode in node.Children)
                {
                    T childArg = childNode.Value;
                    foreach (Tuple<MultidimensionalArray, double>
                        integrationNode in childArg.NodesAndWeights.IntegrationNodes)
                    {
                        integrand(integrationNode.Item1, integrationNode.Item2, arg, false);
                    }
                }
            }
        }

        //Algorithm 1
        //page: A1005
        private void VolumeIntegrand(MultidimensionalArray X, double X_weight, T arg, bool isNew)
        {
            //Calculate the set of roots and union with the interval endpoints ( line 1)
            int heightDirection = arg.HeightDirection;

            //Sort R into ascending order such that ... (line 2) : Implemented with a sortedList
            ISortedBoundedList<double> roots = new SayeSortedList(arg.n * 2);
            double[] bounds = GetBoundaries(arg, heightDirection);
            roots.SetBounds(bounds[0], bounds[1]);
            RestrictToBound(X, bounds[0], arg.HeightDirection);
            for (int i = 0; i < arg.n; ++i) 
            {
                S psi_i = arg.PsiAndS[i].Item1;
                double[] rootsInPsi = FindRoots(psi_i, X, heightDirection, bounds, this.cell);
                roots.Add(rootsInPsi);
            }
            
            //For j = 1 to l - 1 do (line 4)
            bool xIsUnchanged = true;
            for (int j = 0; j < roots.Count() - 1; ++j)
            {
                //Define L and x_c(line 5)
                double L = roots[j + 1] - roots[j];
                NodeSet x_c = NodeOnRay(X, heightDirection, roots[j] - roots[0] + L / 2.0);
                
                //If s_i * psi_i >= 0 for all i (line 6)
                bool updateIntegrand = true;
                foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
                {
                    S psi_i = psiAndS.Item1;
                    int s_i = psiAndS.Item2;
                    updateIntegrand &= s_i * EvaluateAt(psi_i, x_c) >= 0;
                }
                //Update I = ...(line 7)
                if (updateIntegrand)
                {
                    SayeQuadRule newRule = BuildQuadRule(x_c, X_weight, heightDirection, L);
                    bool deriveFromExistingNode = !isNew && xIsUnchanged;
                    arg.NodesAndWeights.AddRule(newRule, deriveFromExistingNode);
                    xIsUnchanged = false;
                }    
            }
            if (xIsUnchanged && !isNew)
            {
                arg.NodesAndWeights.RemoveActiveNode();
            }
        }

        //Algorithm 2
        //page: A1005
        private void SurfaceIntegrand(MultidimensionalArray X, double X_weight, T arg, bool isNew)
        {
            Debug.Assert(arg.n == 1);
            //Calculate the roots of phi in the interval R = ... (line 1)
            S phi = arg.PsiAndS[0].Item1;
            int heightDirection = arg.HeightDirection;
            double[] bounds = GetBoundaries(arg, heightDirection);
            ISortedBoundedList<double> roots = new SayeSortedList(1);
            roots.SetBounds(bounds[0], bounds[1]);
            double[] newRoots = FindRoots(phi, X, heightDirection, bounds, this.cell);
            roots.Add(newRoots);

            //if there is a root, insert node 
            Debug.Assert(roots.Count() <= 3);
            if (roots.Count() > 2)
            {
                X[0, heightDirection] = roots[1];
                SayeQuadRule surfaceQuadNode = BuildSurfaceQuadRule(X, X_weight, heightDirection, this.cell);
                arg.NodesAndWeights.AddRule(surfaceQuadNode, true);
            }
            //else remove node
            else
            {
                arg.NodesAndWeights.RemoveActiveNode();
            }
        }

        #endregion

        #region Evaluate Integrand: Abstract Functions 

        protected abstract void RestrictToBound(MultidimensionalArray X, double bound, int direction);

        protected abstract SayeQuadRule SetLowOrderQuadratureNodes(T arg);

        protected abstract SayeQuadRule SetGaussQuadratureNodes(T arg);

        //Return 1 Root. If there isn't a root, return double.MaxValue 
        protected abstract double[] FindRoots(S psi, MultidimensionalArray X, int heightDirection, double[] bounds, int cell);

        protected abstract SayeQuadRule BuildQuadRule( MultidimensionalArray X, double X_weight, int heightDirection, double length);

        protected abstract SayeQuadRule BuildSurfaceQuadRule(MultidimensionalArray X, double X_weight, int heightDirection, int cell);

        public abstract double[] GetBoundaries(T arg, int heightDirection);

        public abstract NodeSet NodeOnRay(MultidimensionalArray X, int direction, double distance);

        #endregion

        #region BuildIntegrand

        //Algorithm 3
        //page: A1006
        protected void SayeRecursion(TreeNode<T> treeNode)
        {
            T arg = treeNode.Value;

            //Assert ... (line 1)
            Debug.Assert(arg.Surface == false || (arg.Dimension > 1 && arg.n == 1));
            //Check treeNode : Prune
            //-----------------------------------------------------------------------------------------------
            // If d = 1 ... (line 2)
            if (arg.Dimension == 1)
            {
                return;
            }
            //Define x_c (line 3)
            NodeSet x_center = arg.GetCellCenter();

            //for i = n downto 1 do (line 4)
            for (int i = arg.n - 1 ; i >= 0;  --i)
            {
                Tuple<S, int> psiAndS = arg.PsiAndS[i];
                S psi_i = psiAndS.Item1;
                //Evaluate bounds on the value of psi_i on U such that sup_(x element U)|psi_i(x) - psi_i(x_c)| <= delta (line 5)
                double delta = EvaluateBounds(arg, psi_i, x_center);
                // if |psi_i(x_c)| >= delta then (line 6)
                if (delta <= Math.Abs(EvaluateAt(psi_i, x_center)))
                {
                    int s_i = psiAndS.Item2;
                    //if s_i * psi_i >= 0 (line 7)
                    if (s_i * EvaluateAt(psi_i, x_center) >= 0)
                    {
                        //Remove ψi from the list and decrement n by one. (line 8)
                        arg.PsiAndS.RemoveAt(i);
                    }
                    else
                    {
                        //The domain of integration is empty; return 0. (line 9,10)
                        arg.Status = SayeArgument<S>.Mode.DomainIsEmpty;
                        return;
                    }
                }
            }
            //if n = 0 then (line 11)
            if (arg.n == 0)
            {
                //If this is a Surface Integral, this Subcell is empty
                if (arg.Surface)
                {
                    arg.Status = SayeArgument<S>.Mode.DomainIsEmpty;
                }
                else
                {
                    //Apply a tensor - product Gaussian quadrature scheme (line 12)
                    arg.Status = SayeArgument<S>.Mode.GaussQuadrature;
                }
                return;
            }

            //Find new subspace
            //-----------------------------------------------------------------------------------------------
            //Set k = argmax_j|d_xj Psi_1(x_c)| and initialize Psi-tilde = empty (line 14)
            S psi_1 = arg.PsiAndS[0].Item1;
            int k = FindPromisingHeightDirection(psi_1, x_center);
            T subspaceArg = DeriveNewArgument(arg);
            //for i = 1 to n do (line 15)
            foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
            {
                S psi_i = psiAndS.Item1;
                int s_i = psiAndS.Item2;
                //Evaluate g := gradient(Psi_i(x_c)) (line 16)
                MultidimensionalArray g = Gradient(psi_i, x_center);
                //Determine bounds and check bounds (line 17,18)
                if (HeightDirectionIsSuitable(arg, psi_i, x_center, k, g))
                {
                    //Define Psi_i^L and Psi_i^U (line 19)
                    S[] subPsis = ExtractSubPsis(psi_i, arg, k);
                    S psi_U = subPsis[1];
                    S psi_L = subPsis[0];
                    //Evaluate signs s_i^L and s_i^U (line 20)
                    int s_U = EvaluateSign(g, k, s_i, arg.Surface, 1);
                    int s_L = EvaluateSign(g, k, s_i, arg.Surface, -1);
                    //Add {Psi_i^L, s_i^L}  and {Psi_i^U, s_i^U} to the collection Psi-tilde (line 21)
                    Tuple<S, int> newPsiAndS_U = new Tuple<S, int>(psi_U, s_U);
                    Tuple<S, int> newPsiAndS_L = new Tuple<S, int>(psi_L, s_L);
                    subspaceArg.PsiAndS.Add(newPsiAndS_U);
                    subspaceArg.PsiAndS.Add(newPsiAndS_L);
                }
                else
                {
                    //The height function direction ek is not suitable for ψi. If already subdivided too
                    //many times, revert to a low - order method(see discussion).Otherwise split U (line 23)
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
            //Define a new integrand (line 24)
            arg.HeightDirection = k;
            //return I(...) (line 25)
            subspaceArg.Surface = false;
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

        protected abstract int EvaluateSign(MultidimensionalArray gradient, int heightDirection, int s_i, bool surface, int sigma);

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
