using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG.Quadrature
{
    enum Dim { x, y, z }

    struct Argument
    {
        public int quadOrder;
        public Dim[] activeDims;
        public Dim heightDirection;
        public int nodes;

        public static int CalculateNumberOfNodes(Argument node, int numberOfSubNodes)
        {
            if (numberOfSubNodes == 0)
                return node.nodes;

            return node.quadOrder * numberOfSubNodes;
        }
    }

    class TreeNode<T>
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
