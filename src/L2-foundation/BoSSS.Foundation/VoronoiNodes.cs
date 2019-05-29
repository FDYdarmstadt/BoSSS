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
        static ulong globalID = 0;

        static Queue<ulong> returnedIds = new Queue<ulong>();

        public static ulong GetID()
        {
            ulong id;
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

        public static void ReturnID(ulong id)
        {
            returnedIds.Enqueue(id);
        }
    }

    public class VoronoiNodes : MultidimensionalArrayOrList<VoronoiNode>
    {
        ulong[] globalIds;

        public ulong[] GlobalIds {
            get { return globalIds; }
        }

        MultidimensionalArrayOrList<VoronoiNode> nodes;

        public MultidimensionalArray Positions {
            get { return Array; }
        }

        public IList<VoronoiNode> Nodes {
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
            globalIds = new ulong[Count];
            for(int i = 0; i < Count; ++i)
            {
                globalIds[i] = VoronoiIDProvider.GetID();
            }
        }

        public VoronoiNodes(MultidimensionalArray positions, ulong[] globalIds)
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

        public VoronoiNodes(IList<VoronoiNode> nodes)
            : base(nodes)
        {
            int dimension = nodes[0].Dim;
            InitializeVelocity(Count, dimension);
        }

        protected override IList<VoronoiNode> ToList(MultidimensionalArray positions)
        {
            List<VoronoiNode> nodeList = new List<VoronoiNode>(positions.GetLength(0));
            for (int i = 0; i < positions.GetLength(0); ++i)
            {
                VoronoiNode node = new VoronoiNode( new Vector(positions.GetRow(i)), GlobalIds[i]);
                nodeList.Add(node);
            }
            return nodeList;
        }

        protected override MultidimensionalArray ToArray(IList<VoronoiNode> nodes)
        {
            InitializeVelocity(nodes.Count, nodes[0].Dim);
            MultidimensionalArray positions = MultidimensionalArray.Create(nodes.Count, nodes[0].Dim);
            globalIds = new ulong[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                Vector position = nodes[i].Position;
                positions.SetRowPt(i, position);
                globalIds[i] = nodes[i].GlobalID;
            }
            return positions;
        }
    }

    public class VoronoiNode
    {
        public ulong GlobalID { get; }

        public Vector Position { get; set; }

        public int Dim {
            get { return Position.Dim; }
        }

        public VoronoiNode()
        {
            GlobalID = VoronoiIDProvider.GetID();
        }

        public VoronoiNode(Vector position)
            : this()
        {
            Position = position;
        }

        public VoronoiNode(Vector position, ulong globalId)
        {
            GlobalID = globalId;
            Position = position;
        }
    }

    public abstract class MultidimensionalArrayOrList<T>
    {
        enum DataState { MultiDimensionalArray, List };
        DataState state;

        IList<T> list;
        MultidimensionalArray array;

        public int Count {
            get { return DataStateMatches(DataState.List) ? list.Count : array.GetLength(0); }
        }

        protected IList<T> List {
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

        protected abstract MultidimensionalArray ToArray(IList<T> list);

        protected abstract IList<T> ToList(MultidimensionalArray array);

        public MultidimensionalArrayOrList(IList<T> list)
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
