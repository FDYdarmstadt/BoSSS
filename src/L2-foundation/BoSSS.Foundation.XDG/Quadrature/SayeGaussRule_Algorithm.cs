using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG.Quadrature
{

    public interface IPsi
    {
        double EvaluateAt(double[] point);
    }

    public class SayeArgument<S>
        where S : IPsi
    {

        public int dimension;
        public LinkedList<Tuple<S, int>> PsiAndS;
        public bool gaussQuadrature = false;
        public bool lowOrderQuadrature = false;
        int cell;

        public SayeArgument()
        {
        }

        public SayeArgument(int Cell)
        {
            cell = Cell;
        }

        public int n {
            get {
                return PsiAndS.Count;
            }
        }
    }

    public abstract class SayeIntegrand<S>
        where S : IPsi
    {
        QuadRule rule;
        public enum Dim { x, y, z };
        public SayeIntegrand()
        {
        }

        public QuadRule Evaluate(int Cell, QuadRule Rule1d)
        {
            //Setup algorithm
            rule = Rule1d;
            TreeNode<SayeArgument> recursionTree = new TreeNode<SayeArgument>();
            SayeArgument arg = new SayeArgument(Cell);
            TreeNode<SayeArgument> fullSpace = recursionTree.AddChild(arg);

            //Execute algorithm
            SayeRecursion(fullSpace);
            QuadRule result = EvaluateSayeRecursion();
            return result;
        }

        //Algorithm 2
        private QuadRule EvaluateSayeRecursion()
        {
            throw new NotImplementedException();
        }

        //Algorithm 3
        private void SayeRecursion(TreeNode<SayeArgument<S>> treeNode)
        {
            SayeArgument<S> arg = treeNode.Value;

            //Check treeNode : Prune
            //------------------------------------------------------------------------------------------------------------------------
            if (arg.dimension == 1)
            {
                return;
            }
            SetCellCenter(arg);

            foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
            {
                S psi_i = psiAndS.Item1;
                double delta = EvaluateBounds(psi_i);
                if (delta <= psi_i.EvaluateAt(x_center))
                {
                    int s_i = psiAndS.Item2;
                    if (s_i * psi_i.EvaluateAt(x_center) >= 0)
                    {
                        //Remove ψi from the list and decrement n by one.
                    }
                    else
                    {
                        //The domain of integration is empty; return 0.
                        return;
                    }
                }
            }
            if (arg.n == 0)
            {
                //Apply a tensor - product Gaussian quadrature scheme
                arg.gaussQuadrature = true;
                return;
            }

            //Find new subspace
            //------------------------------------------------------------------------------------------------------------------------
            SayeArgument<S> subspaceArg = new SayeArgument<S>();
            Dim k = FindPromisingHeightDirection();
            foreach (Tuple<S, int> psiAndS in arg.PsiAndS)
            {
                if (HeightDirectionIsSuitable())
                {

                }
                else
                {
                    //The height function direction ek is not suitable for ψi. If already subdivided too
                    //many times, revert to a low - order method(see discussion).Otherwise split U
                    if (SubdivideSuitable(arg))
                    {
                        TreeNode<SayeArgument<S>>[] subdividedNodes = SubdivideNode(treeNode);
                        foreach (TreeNode<SayeArgument<S>> subdividedNode in subdividedNodes)
                        {
                            SayeRecursion(subdividedNode);

                        }
                    }
                    else
                    {
                        arg.lowOrderQuadrature = true;
                    }
                    return;
                }
            }
            TreeNode<SayeArgument<S>> subSpaceNode = treeNode.AddChild(subspaceArg);
            SayeRecursion(subSpaceNode);
        }

        //Abstract funcs

        protected abstract void SetCellCenter(SayeArgument<S> arg);

        protected abstract bool SubdivideSuitable(SayeArgument<S> arg);

        protected abstract TreeNode<SayeArgument<S>>[] SubdivideNode(TreeNode<SayeArgument<S>> node);

        protected abstract Dim FindPromisingHeightDirection();

        protected abstract bool HeightDirectionIsSuitable();

        protected abstract double EvaluateBounds(S psi);

    }

    public class TreeNode<T>
    {
        LinkedList<TreeNode<T>> children = new LinkedList<TreeNode<T>>();
        TreeNode<T> parent;
        T value;

        public TreeNode()
        {
            value = default(T);
            parent = null;
        }

        public TreeNode(T Value)
        {
            value = Value;
            parent = null;
        }

        public TreeNode(T Value, TreeNode<T> Parent)
        {
            value = Value;
            parent = Parent;
        }

        public TreeNode<T> AddSibling(T Sibling)
        {
            if (parent != null)
            {
                TreeNode<T> sibling = parent.AddChild(Sibling);
                return sibling;
            }
            return null;
        }

        public TreeNode<T> AddChild(T Child)
        {
            TreeNode<T> childNode = new TreeNode<T>(Child, this);
            children.AddFirst(childNode);
            return childNode;
        }

        public LinkedList<TreeNode<T>> Children {
            get {
                return children;
            }
        }

        public T Value {
            get { return value; }
        }

        public int SumFuncOnTree(Func<T, int, int> sumFunc)
        {
            int n = 0; //Holds the sum of function( children)
            foreach (TreeNode<T> child in this.Children)
            {
                n += child.SumFuncOnTree(sumFunc);
            }
            return sumFunc(this.value, n);
        }
    }

}
