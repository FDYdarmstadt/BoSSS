using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using BoSSS.Foundation.Grid.RefElements;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Foundation.XDG.Quadrature
{
    public class LinearPSI<S>
        : IPsi
        where S : RefElement

    {
        BitArray fixedDims;
        double[] fixedValues;
        private static S refElement;

        public LinearPSI(S _RefElement)
        {
            refElement = _RefElement;
            fixedDims = new BitArray(refElement.SpatialDimension);
            fixedValues = new double[refElement.SpatialDimension];
        }

        public LinearPSI()
        {
            Debug.Assert(refElement != null);
            fixedDims = new BitArray(refElement.SpatialDimension);
            fixedValues = new double[refElement.SpatialDimension];
        }

        public static int SpatialDimension
        {
            get => refElement.SpatialDimension;
        }

        private LinearPSI(BitArray FixedDims, double[] FixedValues)
        {
            fixedDims = FixedDims;
            fixedValues = FixedValues;
        }

        public LinearPSI<S> ReduceDim(int dimension, double value)
        {
            BitArray FixedDims = new BitArray(refElement.SpatialDimension);
            double[] FixedValues = new double[refElement.SpatialDimension];
            for (int i = 0; i < fixedDims.Count; ++i)
            {
                if (this.fixedDims[i] || (dimension == i))
                {
                    FixedDims[i] = true;
                    if (dimension == i)
                    {
                        FixedValues[i] = value;
                    }
                    else
                    {
                        FixedValues[i] = this.fixedValues[i];
                    }
                }
            }
            return new LinearPSI<S>(FixedDims, FixedValues);
        }

        public NodeSet ProjectOnto(NodeSet nodes)
        {
            NodeSet projection = new NodeSet(refElement, nodes.NoOfNodes, 2);
            for (int j = 0; j < nodes.NoOfNodes; ++j)
            {
                for (int i = 0; i < fixedDims.Count; ++i)
                {
                    if (fixedDims[i])
                    {
                        projection[j, i] = fixedValues[i];
                    }
                    else
                    {
                        projection[j, i] = nodes[j, i];
                    }
                }
            }
            projection.LockForever();
            return projection;
        }

        public void SetInactiveDimsToZero( double[] arr)
        {
            Debug.Assert(arr.Length == SpatialDimension);

            for (int i = 0; i < SpatialDimension; ++i)
            {
                if (fixedDims[i])
                    arr[i] = 0;
            }
        }
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

        public int SumOverTree(Func<T, int, int> sumFunc)
        {
            int n = 0; //Holds the sum of function( children)
            foreach (TreeNode<T> child in this.Children)
            {
                n += child.SumOverTree(sumFunc);
            }
            if (this.value == null)
                return n;
            return sumFunc(this.value, n);
        }

        public void UnrollFunc(Action<TreeNode<T>> Func)
        {
            foreach (TreeNode<T> child in this.Children)
            {
                child.UnrollFunc(Func);
            }
            Func(this);
        }
    }

    public class NodesAndWeightsLinkedList : INodesAndWeights
    {
        static LinkedList<Tuple<MultidimensionalArray, double>> data = new LinkedList<Tuple<MultidimensionalArray, double>>();

        public void Reset()
        {
            data.Clear();
        }

        int spatialDim;

        int length;

        LinkedListNode<Tuple<MultidimensionalArray, double>> startNode;
        static LinkedListNode<Tuple<MultidimensionalArray, double>> activeNode;

        public NodesAndWeightsLinkedList(int SpatialDim)
        {
            spatialDim = SpatialDim;
        }

        public IEnumerable<MultidimensionalArray> IntegrationNodes 
        {
            get 
            {
                int counter = 0;
                activeNode = startNode;
                LinkedListNode<Tuple<MultidimensionalArray, double>> nextNode;
                while (counter < length - 1)
                {
                    nextNode = activeNode.Next;
                    yield return activeNode.Value.Item1;
                    activeNode = nextNode;
                    ++counter;
                }
                yield return activeNode.Value.Item1; ;
            }
        }

        /*
        bool Equals(BitArray A, BitArray B)
        {
            if (A.Length != B.Length)
                return false;
            for(int i = 0; i < A.Length; ++i)
            {
                if (A[i] != B[i])
                    return false;
            }
            return true;
        }
        */

        public void AddRule(QuadRule rule)
        {
            length = rule.NoOfNodes;
            //SetStartNode
            MultidimensionalArray node = rule.Nodes.ExtractSubArrayShallow(new int[] { 0, -1 });
            double weight = rule.Weights[0];
            startNode = new LinkedListNode<Tuple<MultidimensionalArray, double>>( new Tuple<MultidimensionalArray, double>(node, weight));
            data.AddLast(startNode);
            LinkedListNode<Tuple<MultidimensionalArray, double>> workingNode = startNode; 
            for (int i = 1; i < rule.NoOfNodes; ++i)
            {
                node = rule.Nodes.ExtractSubArrayShallow(new int[] { i, -1 });
                weight = rule.Weights[i];
                Tuple<MultidimensionalArray, double> newValue = new Tuple<MultidimensionalArray, double>(node, weight);
                workingNode = data.AddAfter(workingNode, newValue);
            }
            activeNode = startNode;
        }

        public void ExpandRule(QuadRule rule)
        {
            //If StartNode == null then startNode = static activeNode
            if(startNode == null)
            {
                startNode = activeNode;
            }
        }

        public QuadRule GetQuadRule()
        {
            //convert static NodeList to Quadrule
            throw new NotImplementedException();
        }
    }

    public static class Extensions
    {
        public static T[][] Copy<T>(this T[][] arr) 
        {
            T[][] copy = new T[arr.Length][];
            for(int i = 0; i < arr.Length; ++i)
            {
                copy[i] = new T[arr[i].Length];
                Array.Copy(arr[i], copy[i], arr[i].Length);
            }
            return copy;
        }
    }
}
