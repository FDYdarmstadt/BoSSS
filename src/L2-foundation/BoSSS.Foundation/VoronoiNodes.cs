using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid.Voronoi;

namespace BoSSS.Foundation.Grid.Voronoi
{
    static class VoronoiIDProvider
    {
        static long globalID = 0;

        static Queue<long> returnedIds = new Queue<long>();

        public static long GetID()
        {
            long id;
            if(returnedIds.Count > 0)
            {
                id = returnedIds.Dequeue();
            }
            else
            {
                id = globalID;
                ++globalID;
            }
            return id;
        }

        public static void ReturnID(long id)
        {
            returnedIds.Enqueue(id);
        }
    }

    public class VoronoiNodes : MultidimensionalArrayOrList<VoronoiNode>
    {
        long[] globalIds;

        public long[] GlobalIds {
            get { return globalIds; }
        }

        MultidimensionalArrayOrList<VoronoiNode> nodes;

        public MultidimensionalArray Positions {
            get { return Array; }
        }

        public List<VoronoiNode> Nodes {
            get { return List; }
        }

        //Velocity of each cell jCell
        public MultidimensionalArray Velocity;

        void InitializeVelocity(int count, int dimension)
        {
            this.Velocity = MultidimensionalArray.Create(count, dimension);
        }

        public VoronoiNodes(MultidimensionalArray positions)
            :base(positions)
        {
            int dimension = positions.GetLength(1);
            SetGlobalIds();
            InitializeVelocity(Count, dimension);
        }

        void SetGlobalIds()
        {
            globalIds = new long[Count];
            for(int i = 0; i < Count; ++i)
            {
                globalIds[i] = VoronoiIDProvider.GetID();
            }
        }

        public VoronoiNodes(MultidimensionalArray positions, long[] globalIds)
            : base(positions)
        {
            if (positions.GetLength(0) != globalIds.Length)
            {
                throw new Exception("Dimension mismatch of nodes and globalIds.");
            }
            this.globalIds = globalIds;
            int dimension = positions.GetLength(1);
            InitializeVelocity(Count, dimension);
        }

        public VoronoiNodes(List<VoronoiNode> nodes)
            : base(nodes)
        {
            int dimension = nodes[0].Dim;
            InitializeVelocity(Count, dimension);
        }

        protected override List<VoronoiNode> ToList(MultidimensionalArray positions)
        {
            List<VoronoiNode> nodeList = new List<VoronoiNode>(positions.GetLength(0));
            for (int i = 0; i < positions.GetLength(0); ++i)
            {
                VoronoiNode node = new VoronoiNode( new Vector(positions.GetRow(i)), GlobalIds[i]);
                nodeList.Add(node);
            }
            return nodeList;
        }

        protected override MultidimensionalArray ToArray(List<VoronoiNode> nodes)
        {
            InitializeVelocity(nodes.Count, nodes[0].Dim);
            MultidimensionalArray positions = MultidimensionalArray.Create(nodes.Count, nodes[0].Dim);
            globalIds = new long[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                Vector position = nodes[i].Position;
                positions.SetRowPt(i, position);
                globalIds[i] = nodes[i].GlobalID;
            }
            return positions;
        }
    }

    public interface INode
    {
        Vector Position { get; set; }
    }

    public class VoronoiNode : INode
    {
        public long GlobalID { get; }

        public Vector Position { get; set; }

        public int Dim {
            get { return Position.Dim; }
        }

        public VoronoiNode()
        {
            GlobalID = VoronoiIDProvider.GetID();
        }

        public VoronoiNode(Vector position)
        {
            GlobalID = VoronoiIDProvider.GetID();
            Position = position;
        }

        public VoronoiNode(Vector position, long globalId)
        {
            GlobalID = globalId;
            Position = position;
        }
    }

    public abstract class MultidimensionalArrayOrList<T>
    {
        enum DataState { MultiDimensionalArray, List };
        DataState state;

        List<T> list;
        MultidimensionalArray array;

        Func<MultidimensionalArray, List<T>> ArrayToListCast;
        Func<List<T>, MultidimensionalArray> ListToArrayCast;

        protected int Count {
            get { return DataStateMatches(DataState.List) ? list.Count : array.GetLength(0); }
        }

        protected List<T> List {
            get {
                if (!DataStateMatches(DataState.List))
                    UpdateList();
                return list;
            }
            set {
                state = DataState.List;
                list = value;
            }
        }

        bool DataStateMatches(DataState incoming)
        {
            return incoming == state;
        }

        void UpdateList()
        {
            list = ToList(array);
        }

        protected MultidimensionalArray Array {
            get {
                if (!DataStateMatches(DataState.MultiDimensionalArray))
                    UpdateArray();
                return array;
            }
            set {
                state = DataState.MultiDimensionalArray;
                array = value;
            }
        }

        void UpdateArray()
        {
            array = ToArray(list);
        }

        protected abstract MultidimensionalArray ToArray(List<T> list);

        protected abstract List<T> ToList(MultidimensionalArray array);

        public MultidimensionalArrayOrList(List<T> list)
        {
            state = DataState.List;
            this.list = list;
        }

        public MultidimensionalArrayOrList(MultidimensionalArray array)
        {
            state = DataState.MultiDimensionalArray;
            this.array = array;
        }
    }
}
