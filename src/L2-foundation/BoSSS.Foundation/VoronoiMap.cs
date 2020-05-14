using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Voronoi
{
    public class Map
    {
        protected IList<int> map;

        public Map(IList<int> map)
        {
            this.map = map;
        }

        public int Length {
            get { return map.Count; }
        }

        public void SetMapping(int entry, int mapIndice)
        {
            map[mapIndice] = entry;
        }

        public int GetMapping( int mapIndice)
        {
            return map[mapIndice];
        }

        public T GetCorrespondingEntry<T>(int mapIndice, IList<T> list)
        {
            return list[map[mapIndice]];
        }

        public void SetCorrespondingEntry<T>(T entry, int mapIndice, IList<T> list)
        {
            list[map[mapIndice]] = entry;
        }

        public int Max()
        {
            int max = map[0];
            for (int i = 1; i < map.Count; ++i)
            {
                if (map[i].CompareTo(max) > 0)
                {
                    max = map[i];
                }
            }
            return max;
        }

        public int Min()
        {
            int min = map[0];
            for (int i = 1; i < map.Count; ++i)
            {
                if (map[i].CompareTo(min) < 0)
                {
                    min = map[i];
                }
            }
            return min;
        }

        public bool IsBijective()
        {
            HashSet<int> indices = new HashSet<int>();
            //Check if all links are unique
            for (int i = 0; i < map.Count; ++i)
            {
                if (!indices.Contains(map[i]))
                {
                    indices.Add(map[i]);
                }
                else
                {
                    return false;
                }
            }
            return true;
        }

        public bool Combine()
        {
            throw new NotImplementedException();
        }
    }

    public enum Connection { Created, Removed, Remained, Mirror};

    public class ConnectionMap : Map
    {
        IList<Connection> connections;

        public ConnectionMap(IList<Connection> connections, IList<int> map)
            : base(map)
        {
            this.connections = connections;
        }

        public ConnectionMap(Connection connection, IList<int> map)
            : base(map)
        {
            connections = new Connection[map.Count];
            for(int i = 0; i < map.Count; ++i)
            {
                connections[i] = connection;
            }
        }

        public static ConnectionMap CreateEmpty(Connection nodeConnection, int length)
        {
            Connection[] connections = new Connection[length];
            for (int i = 0; i < length; ++i)
            {
                connections[i] = nodeConnection;
            }
            int[] map = new int[length];
            return new ConnectionMap(connections, map);
        }

        public Connection GetConnection(int mapIndice)
        {
            return connections[mapIndice];
        }

        public void SetConnection(Connection entry, int mapIndice)
        {
            connections[mapIndice] = entry;
        }

        public void AddReverse(ConnectionMap towardsThis)
        {
            for (int i = 0; i < towardsThis.Length; ++i)
            {
                SetMapping(i, towardsThis.GetMapping(i));
                SetConnection(Connection.Remained, i);
            }
        }
    }
}
