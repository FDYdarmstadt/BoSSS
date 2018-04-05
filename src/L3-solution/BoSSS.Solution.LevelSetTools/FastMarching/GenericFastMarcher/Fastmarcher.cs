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
using System.Collections;

namespace BoSSS.Solution.LevelSetTools.FastMarcher {
 
    public class Fastmarcher {

        HashSet<IMarchingNode> accepted; //Holds all the nodes that are already accepted
        HashSet<IMarchingNode> considered; //Holds all the nodes that are considered
        IFastMarchingQueue<IMarchingNode> consideredQueue; //Always knows the node with the minimal value

        /// <summary>
        /// A generic Fastmarcher.
        /// This algorithm runs on a grap with nodes of type IMarchingNode. Each node holds a value and implements 
        /// a method that can calculate this value from neighboring nodes. The algorithm starts from a set of initial nodes with 
        /// given value. It then picks the lowest value for each node.  
        /// See https://en.wikipedia.org/wiki/Fast_marching_method.
        /// </summary>
        /// <param name="emptyHeap">
        /// This heap determines the marching order.
        /// </param>  
        public Fastmarcher(IFastMarchingQueue<IMarchingNode> emptyHeap) {
            if (emptyHeap.number() != 0) {
                throw new NotSupportedException("Heap must be empty");
            }
            consideredQueue = emptyHeap;
        }

        /// <summary>
        /// Will solve the fastmarching problem, i.e, assign a value to each node in the graph. 
        /// The method needs a graph as input that consists of IMarchingNodes.
        /// </summary>
        /// <param name="initialNodes">
        /// The initially accepted nodes of the fastmarching graph.
        /// </param> 
        public void march(IMarchingNode[] initialNodes) {

            //Initialize lists
            initializeAccepted(initialNodes);
            initializeConsidered(initialNodes);

            while (consideredQueue.number() != 0) {
                //update accepted nodes
                IMarchingNode node = consideredQueue.extractMinimum();
                node.Accept();
                accepted.Add(node);

                //update considered nodes
                if (consideredQueue.number() != 0) {
                    updateConsidered(node);
                }
            }
        }

        void updateConsidered(IMarchingNode node) {
            //add neighbours that are not in accepted to considered
            //Todo? : Speed up by only going through neighbors that have not been accepted by handling the neighborslist
            foreach(IMarchingNode neighbor in node.Neighbors) {
                if (!accepted.Contains(neighbor)) {
                    neighbor.CalculateValue();
                    if (considered.Contains(neighbor)) {
                        consideredQueue.changeKey(neighbor.QueueID, neighbor);
                    } else {
                        neighbor.QueueID = consideredQueue.insert(neighbor);
                        considered.Add(neighbor);
                    }
                }
            }     
        }

        void initializeConsidered(IMarchingNode[] initialNodes) {
            considered = new HashSet<IMarchingNode>();
            foreach (IMarchingNode node in accepted) {
                updateConsidered(node);
            }
        }

        void initializeAccepted(IMarchingNode[] initialNodes) {
            accepted = new HashSet<IMarchingNode>();
            foreach (IMarchingNode node in initialNodes) {
                accepted.Add(node);
            }
        }

    }
}

