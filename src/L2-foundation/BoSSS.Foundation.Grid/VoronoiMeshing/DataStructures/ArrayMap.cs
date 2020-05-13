using System;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{
    class ArrayMap<TValue> : IDictionary<int, TValue>
    {
        BitArray inUse;

        TValue[] values;

        public int Capacity { get; private set; }

        public ICollection<int> Keys => throw new NotImplementedException();

        public ICollection<TValue> Values => throw new NotImplementedException();

        public int Count { get; private set; }

        public bool IsReadOnly => throw new NotImplementedException();

        public TValue this[int key] {
            get {
                return values[key];
            }
            set {
                Add(key, value);
            }
        }

        public ArrayMap(int capacity)
        {
            this.Capacity = capacity;
            Setup();
        }

        void Setup()
        {
            this.Count = 0;
            values = new TValue[Capacity];
            inUse = new BitArray(Capacity);
        }

        public void Add(int key, TValue value)
        {
            if (key < Capacity && key > -1)
            {
                if (!inUse[key])
                {
                    ++Count;
                    inUse[key] = true;
                }
                values[key] = value;
            }
            else
            {
                throw new Exception("Key out of bounds.");
            }
        }

        public void Add(KeyValuePair<int, TValue> item)
        {
            Add(item.Key, item.Value);
        }

        public bool ContainsKey(int key)
        {
            if (key < Capacity && key > -1)
            {
                return inUse[key];
            }
            else
            {
                throw new Exception("Key out of bounds.");
            }
        }

        public bool Contains(KeyValuePair<int, TValue> item)
        {
            return ContainsKey(item.Key);
        }

        public bool Remove(int key)
        {
            if (key < Capacity && key > -1)
            {
                if (ContainsKey(key))
                {
                    --Count;
                    values[key] = default(TValue);
                    inUse[key] = false;
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                throw new Exception("Key out of bounds.");
            }
        }

        public bool Remove(KeyValuePair<int, TValue> item)
        {
            return Remove(item.Key);
        }

        public bool TryGetValue(int key, out TValue value)
        {
            if (ContainsKey(key))
            {
                value = values[key];
                return true;
            }
            else
            {
                value = default(TValue);
                return false;
            }
        }

        IEnumerable<KeyValuePair<int, TValue>> Enumerate()
        {
            for (int i = 0; i < Capacity; ++i)
            {
                if (ContainsKey(i))
                {
                    yield return new KeyValuePair<int, TValue>(i, values[i]);
                }
            }
        }

        public IEnumerator<KeyValuePair<int, TValue>> GetEnumerator()
        {
            return (IEnumerator<KeyValuePair<int, TValue>>)Enumerate();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public void Clear()
        {
            Setup();
        }

        public void CopyTo(KeyValuePair<int, TValue>[] array, int arrayIndex)
        {
            throw new NotImplementedException();
        }

        
    }
}
