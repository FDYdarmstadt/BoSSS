using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi
{
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
            InitializeVelocity(Count, dimension);
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
                VoronoiNode node = new VoronoiNode
                {
                    Position = new Vector(positions.GetRow(i)),
                    GlobalID = GlobalIds[i]
                };
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

    public class VoronoiNode
    {
        public long GlobalID { get; set; }

        public Vector Position { get; set; }

        public int Dim {
            get { return Position.Dim; }
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
