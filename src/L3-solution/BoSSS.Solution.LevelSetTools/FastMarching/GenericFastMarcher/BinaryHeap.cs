/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.FastMarcher {
    public class BinaryHeap<T> : IBoundedPriorityQueue<T> where T : IComparable<T> {
        int Nmax;
        int N;
        protected HeapItem<T>[] TheHeap;
        protected int[] Positions;
        LinkedList<int> unused;

        public BinaryHeap(int size) {
            Nmax = size;
            TheHeap = new HeapItem<T>[size];
            Positions = new int[size];
            N = 0;
            unused = new LinkedList<int>();
            for (int i = 1; i < Nmax + 1; ++i) {
                unused.AddLast(i);
            }
        }

        public int insert(T key) {
            //get ID
            int id = unused.First.Value;
            unused.RemoveFirst();

            //Write into arrays
            Positions[id - 1] = N;
            HeapItem<T> newItem = new HeapItem<T>(id, key);
            TheHeap[N] = newItem;

            //Swap until invariant is fullfilled (going upwards in binary Tree)
            int index = N + 1;
            while (index / 2 > 0 && TheHeap[index / 2 - 1].key.CompareTo(key) > 0) {
                TheHeap[index - 1] = TheHeap[index / 2 - 1];
                TheHeap[index / 2 - 1] = newItem;

                Positions[TheHeap[index - 1].ID - 1] = index - 1;
                Positions[newItem.ID - 1] = index / 2 - 1;
                index = index / 2;
            }
            //increase length
            ++N;
            return id;
        }

        public T extractMinimum() {
            if (N == 0) {
                throw new NotSupportedException("Heap is empty.");
            }
            T key = TheHeap[0].key;
            int id = TheHeap[0].ID;

            //Replace first Item with last.
            TheHeap[0] = TheHeap[N - 1];
            TheHeap[N - 1] = null;
            HeapItem<T> newItem = TheHeap[0];

            //decrease N and add freed ID to unused list. 
            --N;
            unused.AddFirst(id);

            //Update Positions only when Heap not empty
            if (N == 0) {
                return key;
            }
            Positions[newItem.ID - 1] = 0;

            //Swap until invariant is fullfilled (going downwards in binary Tree)
            int index = 1;
            int newIndex;
            while (index * 2 < N + 1) {
                if (index * 2 == N) {
                    if (TheHeap[index * 2 - 1].key.CompareTo(newItem.key) > 0) {
                        break;
                    }
                    newIndex = index * 2;
                } else {
                    if (TheHeap[index * 2 - 1].key.CompareTo(newItem.key) < 0 || TheHeap[index * 2].key.CompareTo(newItem.key) < 0) {
                        if (TheHeap[index * 2 - 1].key.CompareTo(TheHeap[index * 2].key) < 0) {
                            newIndex = index * 2;
                        } else {
                            newIndex = index * 2 + 1;
                        }
                    } else {
                        break;
                    }
                }

                TheHeap[index - 1] = TheHeap[newIndex - 1];
                TheHeap[newIndex - 1] = newItem;

                Positions[TheHeap[index - 1].ID - 1] = index - 1;
                Positions[newItem.ID - 1] = newIndex - 1;
                index = newIndex;
            }

            return key;
        }

        public T findMinimum(T key) {
            return TheHeap[0].key;
        }

        public void decreaseKey(int ID, T key) {

            int index = Positions[ID - 1] + 1;
            HeapItem<T> newItem = new HeapItem<T>(ID, key);

            if (TheHeap[index - 1].key.CompareTo(key) > 0) {
                TheHeap[index - 1] = newItem;
                while (index / 2 > 0 && TheHeap[index / 2 - 1].key.CompareTo(key) > 0) {
                    TheHeap[index - 1] = TheHeap[index / 2 - 1];
                    TheHeap[index / 2 - 1] = newItem;

                    Positions[TheHeap[index - 1].ID - 1] = index - 1;
                    Positions[newItem.ID - 1] = index / 2 - 1;
                    index = index / 2;
                }
            }
        }

        public int number() {
            return N;
        }

        public void changeKey(int ID, T key) {
            int index = Positions[ID - 1] + 1;
            if (TheHeap[index - 1].key.CompareTo(key) == 0) {
                return;
            }
            if (TheHeap[index - 1].key.CompareTo(key) > 0) {
                decreaseKey(ID, key);
            } else {
                increaseKey(ID, key);
            }
        }

        public void increaseKey(int ID, T key) {

            int index = Positions[ID - 1] + 1;
            HeapItem<T> newItem = new HeapItem<T>(ID, key);

            if (TheHeap[index - 1].key.CompareTo(key) < 0) {
                TheHeap[index - 1] = newItem;
                int newIndex;
                while (index * 2 < N + 1) {
                    if (index * 2 == N) {
                        if (TheHeap[index * 2 - 1].key.CompareTo(newItem.key) > 0) {
                            break;
                        }
                        newIndex = index * 2;
                    } else {
                        if (TheHeap[index * 2 - 1].key.CompareTo(newItem.key) < 0 || TheHeap[index * 2].key.CompareTo(newItem.key) < 0) {
                            if (TheHeap[index * 2 - 1].key.CompareTo(TheHeap[index * 2].key) < 0) {
                                newIndex = index * 2;
                            } else {
                                newIndex = index * 2 + 1;
                            }
                        } else {
                            break;
                        }
                    }

                    TheHeap[index - 1] = TheHeap[newIndex - 1];
                    TheHeap[newIndex - 1] = newItem;

                    Positions[TheHeap[index - 1].ID - 1] = index - 1;
                    Positions[newItem.ID - 1] = newIndex - 1;
                    index = newIndex;
                }
            }
        }
    }

    public class HeapItem<T> {
        public int ID;
        public T key;

        public HeapItem(int id, T Key) {
            ID = id;
            key = Key;
        }
    }
}
