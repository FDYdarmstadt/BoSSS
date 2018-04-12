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
    //A Heap that wraps arround the BinaryHeap to allow for easy manipulation of Phi
    class MarchingHeap : IFastMarchingQueue<IMarchingNode> {
        BinaryHeap<WrapperINode> Heap; 
        
        public MarchingHeap(int Size) {
            Heap = new BinaryHeap<WrapperINode>(Size);
        }

        public void changeKey(int ID, IMarchingNode key) {
            WrapperINode Key = new WrapperINode(key);
            Heap.changeKey(ID, Key);
        }

        public int insert(IMarchingNode key) {
            WrapperINode Key = new WrapperINode(key);
            int ID = Heap.insert(Key);
            return ID;
        }

        public IMarchingNode extractMinimum() {
            return Heap.extractMinimum().ToNode();
        }

        public IMarchingNode findMinimum(IMarchingNode key) {
            WrapperINode Key = new WrapperINode(key);
            return Heap.findMinimum(Key).ToNode();
        }

        public void decreaseKey(int ID, IMarchingNode key) {
            WrapperINode Key = new WrapperINode(key);
            Heap.decreaseKey(ID, Key);
        }

        public int number() {
            return Heap.number();
        }
    }

    class WrapperINode : IComparable<WrapperINode> {
        double value;
        IMarchingNode node;

        public WrapperINode(IMarchingNode Node) {
            node = Node;
            value = Node.Value;
        }

        public IMarchingNode ToNode() {
            return node;
        }

        public int CompareTo(WrapperINode OtherNode) {
            if (value > OtherNode.value) {
                return 1;
            } else {
                if (value < OtherNode.value) {
                    return -1;
                }
                return 0;
            }
        }
    }
}
