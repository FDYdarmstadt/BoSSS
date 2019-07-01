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

    public class VoronoiNodes : DataChameleon<MultidimensionalArray[], IList<VoronoiNode>>
    {
        ulong[] globalIds;

        public ulong[] GlobalIds {
            get { return globalIds; }
        }

        public MultidimensionalArray Positions {
            get { return Red[0]; }
        }

        public IList<VoronoiNode> Nodes {
            get { return Blue; }
        }

        //Velocity of each cell jCell
        public MultidimensionalArray Velocity {
            get { return Red[1]; }
        }

        public VoronoiNodes(MultidimensionalArray positions)
            : base(new MultidimensionalArray[]
            {
                positions,
                InitializeVelocityFrom(positions)
            })
        {
            SetGlobalIds();
        }

        public VoronoiNodes(MultidimensionalArray positions, ulong[] globalIds)
            : base(new MultidimensionalArray[]
            {
                positions,
                InitializeVelocityFrom(positions)
            })
        {
            if (positions.GetLength(0) != globalIds.Length)
            {
                throw new Exception("Dimension mismatch of nodes and globalIds.");
            }
            this.globalIds = globalIds;
        }

        public VoronoiNodes(MultidimensionalArray positions, MultidimensionalArray velocity, ulong[] globalIds)
            : base(new MultidimensionalArray[]
            {
                positions,
                velocity
            })
        {
            if (positions.GetLength(0) != globalIds.Length)
            {
                throw new Exception("Dimension mismatch: Cannot create Nodes.");
            }
            if (velocity.GetLength(0) != globalIds.Length)
            {
                throw new Exception("Dimension mismatch: Cannot create Nodes.");
            }
            this.globalIds = globalIds;
        }

        public VoronoiNodes(IList<VoronoiNode> nodes)
            : base(nodes)
        { }

        static MultidimensionalArray InitializeVelocityFrom(MultidimensionalArray positions)
        {
            return MultidimensionalArray.Create(positions.Lengths);
        } 

        void SetGlobalIds()
        {
            int numberOfNodes = Positions.GetLength(0);
            globalIds = new ulong[numberOfNodes];
            for(int i = 0; i < numberOfNodes; ++i)
            {
                globalIds[i] = VoronoiIDProvider.GetID();
            }
        }

        protected override IList<VoronoiNode> ToBlue(MultidimensionalArray[] positionsAndVelocity)
        {
            int numberOfNodes = positionsAndVelocity[0].GetLength(0);
            List<VoronoiNode> nodeList = new List<VoronoiNode>(numberOfNodes);
            for (int i = 0; i < numberOfNodes; ++i)
            {
                VoronoiNode node = new VoronoiNode(
                    new Vector(positionsAndVelocity[0].GetRow(i)),
                    new Vector(positionsAndVelocity[1].GetRow(i)),
                    GlobalIds[i]);
                nodeList.Add(node);
            }
            return nodeList;
        }

        protected override MultidimensionalArray[] ToRed(IList<VoronoiNode> nodes)
        {
            MultidimensionalArray positions = MultidimensionalArray.Create(nodes.Count, nodes[0].Dim);
            MultidimensionalArray velocities = MultidimensionalArray.Create(nodes.Count, nodes[0].Dim);
            globalIds = new ulong[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                Vector position = nodes[i].Position;
                positions.SetRowPt(i, position);
                Vector velocity = nodes[i].Velocity;
                velocities.SetRowPt(i, velocity);
                globalIds[i] = nodes[i].GlobalID;
            }
            return new[] { positions, velocities };
        }

        public int Count {
            get { return State == Color.Red ? Red[0].GetLength(0) : Blue.Count; }
        }
    }

    public class VoronoiNode
    {
        public ulong GlobalID { get; }

        public Vector Position { get; set; }

        public Vector Velocity{get; set;}

        public int Dim {
            get { return Position.Dim; }
        }

        public VoronoiNode()
        {
            GlobalID = VoronoiIDProvider.GetID();
        }

        public VoronoiNode(Vector position, Vector velocity)
            : this()
        {
            Position = position;
            Velocity = velocity;
        }

        public VoronoiNode(Vector position, Vector velocity, ulong globalId)
        {
            GlobalID = globalId;
            Position = position;
            Velocity = velocity; 
        }
    }

    public abstract class DataChameleon<TRed, TBlue>
    {
        protected enum Color { Red, Blue };

        protected Color State {
            get;
            private set;
        }

        TBlue blue;

        TRed red;

        protected TBlue Blue {
            get {
                if (!DataStateMatches(Color.Blue))
                    UpdateBlue();
                State = Color.Blue;
                return blue;
            }
            set {
                State = Color.Blue;
                blue = value;
            }
        }

        protected TRed Red {
            get {
                if (!DataStateMatches(Color.Red))
                    UpdateRed();
                State = Color.Red;
                return red;
            }
            set {
                State = Color.Red;
                red = value;
            }
        }

        public DataChameleon(TBlue blue)
        {
            State = Color.Blue;
            this.blue = blue;
        }

        public DataChameleon(TRed red)
        {
            State = Color.Red;
            this.red = red;
        }
        
        bool DataStateMatches(Color incoming)
        {
            return incoming == State;
        }

        void UpdateBlue()
        {
            blue = ToBlue(red);
        }

        void UpdateRed()
        {
            red = ToRed(blue);
        }

        protected abstract TRed ToRed(TBlue data);

        protected abstract TBlue ToBlue(TRed data);
    }
}
