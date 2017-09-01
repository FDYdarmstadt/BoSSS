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
using ilPSP;
using System.Diagnostics;

namespace BoSSS.Foundation.Caching {

    /// <summary>
    /// Options for Different Caching Strategies
    /// </summary>
    public enum CacheStrategy {
        /// <summary>
        /// Cache object with lowest number of accesses will be discadred first.
        /// </summary>
        DiscardLeastFrequentlyUsed,

        /// <summary>
        /// Cache object wich were not accessed for the longest time will be discarded first.
        /// </summary>
        DiscardLeastRecentyUsed
    }

    
    /// <summary>
    /// Cache for almost everything. All kinds of cache-logic objects
    /// may place their numerical results here and can have temy back, if they are not thrown away.
    /// </summary>
    public static class Cache {

        
        /// <summary>
        /// Reference counter
        /// </summary>
        static internal uint RefCounter = 1;
        
        /// <summary>
        /// Maximal memory, in bytes, which is used for caching.
        /// </summary>
        public static int MaxMem {
            get;
            set;
        }

        static CacheStrategy m_Strategy;

        /// <summary>
        /// Caching strategy; may affect performance, or not.
        /// </summary>
        public static CacheStrategy Strategy {
            get {
                return m_Strategy;
            }
            set {
                m_Strategy = value;
                ResetCache();
            }
        }

        static Cache() {
            MaxMem = 128*1024*1024; // a few megs
            m_Strategy = CacheStrategy.DiscardLeastRecentyUsed;
            ResetCache();
        }

        static void ResetCache() {
            Tail = new CacheBank() {
                iBank = int.MinValue
            };
            Head = new CacheBank() {
                iBank = int.MinValue
            };

            Tail.prev = Head;
            Head.next = Tail;

            emptyBanks = new List<int>();
            Banks = new List<CacheBank>();
            CurrentSize = 0;
            m_NoOfUsedBanks = 0;

        }


        /// <summary>
        /// Number of cache hits.
        /// </summary>
        static public int Hits {
            get;
            private set;
        }

        /// <summary>
        /// Number of cache misses.
        /// </summary>
        static public int Misses {
            get;
            private set;
        }


        /// <summary>
        /// Returns the item stored under reference <paramref name="Reference"/>
        /// -- or --
        /// null, if the object has been removed from the cache.
        /// </summary>
        public static object GetItem(ulong Reference) {
            int iBank = (int)(Reference & 0x00000000FFFFFFFF);
            int iref = (int)((Reference & 0xFFFFFFFF00000000) >> 32);

            CacheBank cb = Banks[iBank];
            if(cb != null && iref == cb.Reference) {
                Access(iBank);
                Hits++;
                return cb.Item;
            } else {
                Misses++;
                return null;
            }
        }


        /// <summary>
        /// Returns true if the item stored under reference <paramref name="Reference"/> is still present in the cache, false otherwise.
        /// In contrast to <see cref="GetItem"/>, this is cheaper and does not affect cache statistics.
        /// </summary>
        public static bool IsAlive(ulong Reference) {
            int iBank = (int)(Reference & 0x00000000FFFFFFFF);
            int iref = (int)((Reference & 0xFFFFFFFF00000000) >> 32);

            CacheBank cb = Banks[iBank];
            return (cb != null && iref == cb.Reference);
        }

        /// <summary>
        /// Caches some object <paramref name="Item"/>.
        /// </summary>
        /// <param name="HonestItemSizeInBytes">The size of <paramref name="Item"/> in bytes.</param>
        /// <param name="Item">Object which should be cached.</param>
        /// <returns>
        /// A reference code for later access to the object <paramref name="Item"/>.
        /// </returns>
        public static ulong CacheItem(object Item, int HonestItemSizeInBytes) {
            RefCounter++;

            Debug.Assert((!(Item is MultidimensionalArray)) || (((MultidimensionalArray)Item).Length * sizeof(double) == HonestItemSizeInBytes));

            CacheBank cb = new CacheBank() {
                Item = Item,
                MemSize = HonestItemSizeInBytes,
                UseCount = 1,
                Reference = RefCounter
            };

            while(CurrentSize + cb.MemSize > MaxMem) {
                if(!RemoveAtTail())
                    break;
            }

            int iBank = Insert(cb);
            Debug.Assert(object.ReferenceEquals(Banks[iBank], cb));
            Debug.Assert(cb.iBank == iBank);

            ulong R = ((ulong)iBank);
            R |= (((ulong)cb.Reference) << 32);
            return R;
        }

        /// <summary>
        /// Removes some cached item from the cache.
        /// </summary>
        public static void ForgetItem(ulong Reference) {
            int iBank = (int)(Reference & 0x00000000FFFFFFFF);
            int iref = (int)((Reference & 0xFFFFFFFF00000000) >> 32);

            CacheBank cb = Banks[iBank];
            if((cb != null && iref == cb.Reference)) {
                // object is still in cache, so there is something to forget.
                Remove(iBank);
            }

        }

        public static ulong ReCacheItem(MultidimensionalArray Item, int HonestItemSizeInBytes, ulong Reference) {
            return CacheItem(Item, HonestItemSizeInBytes); 
        }


        static List<int> emptyBanks;
        static List<CacheBank> Banks;
        static int CurrentSize;
        static int m_NoOfUsedBanks;

        public static int NoOfUsedBanks {
            get {
                return m_NoOfUsedBanks;
            }
        }
        
        static CacheBank Tail; // dummy 'bank', marking the tail of the cache-ranking list (next objects to discard)

        static CacheBank Head; // dummy 'bank', marking the tail of the cache-ranking list (last object to discard)
        
        static void Remove(int iBank) {
            CacheBank toRemove = Banks[iBank];
            Debug.Assert(toRemove.iBank == iBank);
            Banks[iBank] = null;
            toRemove.prev.next = toRemove.next;
            toRemove.next.prev = toRemove.prev;
            CurrentSize = CurrentSize - toRemove.MemSize;
            emptyBanks.Add(iBank);
            m_NoOfUsedBanks--;
            Debug.Assert(m_NoOfUsedBanks >= 0);
        }
        
        static bool RemoveAtTail() {
            Debug.Assert(Banks.Count > 0);
            if(object.ReferenceEquals(Tail.prev, Head)) {
                // cache is already empty
                Debug.Assert(m_NoOfUsedBanks == 0);
                Debug.Assert(CurrentSize == 0);
                Debug.Assert(object.ReferenceEquals(Head.next, Tail));
                return false;
            } else {
                Remove(Tail.prev.iBank);
                return true;
            }
        }

        static CacheBank Access(int iBank) {
            CacheBank a = Banks[iBank];
            if(a == null)
                return a;

            a.UseCount++;

            // update the bank ranking:
            switch(Strategy) {
                case CacheStrategy.DiscardLeastRecentyUsed:
                //
                  
                // remove 'a', where-ever it is...
                a.prev.next = a.next;
                a.next.prev = a.prev;

                // ... and put it on top:
                a.next = Head.next;
                a.prev = Head;
                Head.next.prev = a;
                Head.next = a;

                break;
                
                case CacheStrategy.DiscardLeastFrequentlyUsed:
                //

                // bubble up until the sorting is ok.
                BubbleUp(a);

                break;

                default:
                throw new NotImplementedException();
            }

            return a;
        }

        private static void BubbleUp(CacheBank a) {
            while(
                (!object.ReferenceEquals(a.prev, Head))  // reached top of ranking
                && (a.UseCount >= a.prev.UseCount))     // >= prefers least-recently used in case of doubt
                {
                CacheBank p = a.prev;
                CacheBank pp = p.prev;
                CacheBank n = a.next;

                pp.next = a;
                a.next = p;
                p.next = n;

                a.prev = pp;
                p.prev = a;
                n.prev = p;
            }
        }

        static int Insert(CacheBank cb) {
            if(emptyBanks.Count > 0) {
                cb.iBank = emptyBanks[emptyBanks.Count - 1];
                emptyBanks.RemoveAt(emptyBanks.Count - 1);
                Debug.Assert(Banks[cb.iBank] == null);
                Banks[cb.iBank] = cb;
            } else {
                cb.iBank = Banks.Count;
                Banks.Add(cb);
            }
            m_NoOfUsedBanks++;
            CurrentSize += cb.MemSize;// update the bank ranking:
            switch(Strategy) {
                case CacheStrategy.DiscardLeastRecentyUsed:
                // insert at head
                cb.next = Head.next;
                cb.prev = Head;
                Head.next.prev = cb;
                Head.next = cb;

                break;

                case CacheStrategy.DiscardLeastFrequentlyUsed:
                // insert at tail
                cb.next = Tail;
                cb.prev = Tail.prev;
                Tail.prev.next = cb;
                Tail.prev = cb;

                // bubble up until the sorting is ok.
                BubbleUp(cb);

                break;

                default:
                throw new NotImplementedException();
            }
            
            Debug.Assert(Banks[cb.iBank].iBank == cb.iBank);
            return cb.iBank;
        }
        
        class CacheBank {
            public uint Reference; // reference by which the cache-logic identifies that some stored item is still theirs
            public int iBank;      // bank index, constant for object lifetime

            
            public CacheBank next; // double-linked list for cache bank ranking
            public CacheBank prev; // double-linked list for cache bank ranking
            
            public object Item; // payload
            
            public int MemSize; // size of the payload

            internal int UseCount; // how often the item was accessed
        }





    }
}
