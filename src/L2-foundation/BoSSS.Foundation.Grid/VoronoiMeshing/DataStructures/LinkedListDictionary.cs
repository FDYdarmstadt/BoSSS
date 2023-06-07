using System;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{

    /// <summary>
    /// LinkedList implementation of Dictionary. Needs minimal ammount of memory, 
    /// but has O(n) complexity for Contains Key/Value and TryGetValue.
    /// Should be used for very small (less than 10) dictionary.Count. 
    /// </summary>
    /// <typeparam name="TKey"></typeparam>
    /// <typeparam name="TValue"></typeparam>
    class LinkedListDictionary<TKey, TValue> : IDictionary<TKey, TValue>
    {
        LinkedList<KeyValuePair<TKey, TValue>> keysAndValues;

        public LinkedListDictionary()
        {
            keysAndValues = new LinkedList<KeyValuePair<TKey, TValue>>();
        }

        public ICollection<TKey> Keys {
            get{
                TKey[] keys = new TKey[keysAndValues.Count];
                int i = 0;
                foreach(var node in keysAndValues)
                {
                    keys[i] = node.Key;
                    ++i;
                }
                return keys;
            }
        }

        public ICollection<TValue> Values {
            get {
                TValue[] values = new TValue[keysAndValues.Count];
                int i = 0;
                foreach (var node in keysAndValues)
                {
                    values[i] = node.Value;
                    ++i;
                }
                return values;
            }
        }

        public int Count => keysAndValues.Count;

        public bool IsReadOnly => false;

        public TValue this[TKey key] {
            get {
                if(TryGetValue(key, out TValue value))
                {
                    return value;
                }
                else
                {
                    throw new IndexOutOfRangeException("Key not in dictionary.");
                }
            }
            set {
                Add(key, value);
            }
        }

        public bool ContainsKey(TKey key)
        {
            foreach(var node in keysAndValues)
            {
                TKey nodeKey = node.Key;
                if (nodeKey.Equals(key))
                {
                    return true;
                }
            }
            return false;
        }

        public void Add(TKey key, TValue value)
        {
            keysAndValues.AddLast(new KeyValuePair<TKey, TValue>(key, value));
        }

        public bool Remove(TKey key)
        {
            var node = keysAndValues.First;
            while( node != null)
            {
                TKey nodeKey = node.Value.Key;
                if (nodeKey.Equals(key))
                {
                    keysAndValues.Remove(node);
                    return true;
                }
                node = node.Next;
            }
            return false;
        }

        public bool TryGetValue(TKey key, out TValue value)
        {
            foreach(var node in keysAndValues)
            {
                if (node.Key.Equals(key))
                {
                    value = node.Value;
                    return true;
                }
            }
            value = default(TValue);
            return false;
        }

        public void Add(KeyValuePair<TKey, TValue> item)
        {
            keysAndValues.AddLast(item);
        }

        public void Clear()
        {
            keysAndValues.Clear();
        }

        public bool Contains(KeyValuePair<TKey, TValue> item)
        {
            return keysAndValues.Contains(item);
        }

        public void CopyTo(KeyValuePair<TKey, TValue>[] array, int arrayIndex)
        {
            keysAndValues.CopyTo(array, arrayIndex);
        }

        public bool Remove(KeyValuePair<TKey, TValue> item)
        {
            return keysAndValues.Remove(item);
        }

        public IEnumerator<KeyValuePair<TKey, TValue>> GetEnumerator()
        {
            return keysAndValues.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
