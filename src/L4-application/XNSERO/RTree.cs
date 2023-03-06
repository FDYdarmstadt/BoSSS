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

using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSERO_Solver {
    public class RTree {

        private readonly int SpatialDimension;
        private readonly int MaxEntriesPerNode = 4;
        private readonly List<TreeNode> Nodes = new List<TreeNode>();
        private readonly double Tolerance;

        /// <summary>
        /// Constructor for an R-Tree (Rectangle-tree)
        /// </summary>
        /// <param name="SpatialDimension"></param>
        /// <param name="Tolerance">Tolerance parameter for periodic boundaries. Only effects performance, not accuracy. Should be at least in the order of the minimal grid length.</param>
        public RTree(int SpatialDimension, double Tolerance) {
            this.SpatialDimension = SpatialDimension;
            this.Tolerance = Tolerance;
            if (this.SpatialDimension != 2)
                throw new NotImplementedException("R-Tree only implemented for 2D");
        }

        /// <summary>
        /// Initializes an R-Tree based on the minimal bounding rectangle (MBR) of the particles.
        /// </summary>
        /// <param name="Particles"></param>
        /// <param name="Timestep"></param>
        public void InitializeTree(Particle[] Particles, double Timestep) {
            using (new FuncTrace()) {
                TreeNode firstNode = new(true, -1, 0, new Vector(SpatialDimension));
            Nodes.Add(firstNode);
            for(int p = 0; p < Particles.Length; p++) {
                List<Vector> virtualDomainOrigin = Particles[p].Motion.OriginInVirtualPeriodicDomain;
                TreeNode particleNode = new(false, 0, FindSmallestEmptyID(), new Vector(SpatialDimension));
                InsertNodeToList(particleNode);
                particleNode.ParticleID = p;
                particleNode.MBR = CalculateParticleMBR(Particles[p], Timestep, new Vector(SpatialDimension));
                InsertParticle(particleNode, Nodes[0]);
                    for (int i = 0; i < virtualDomainOrigin.Count; i++) {
                        Vector virtualPosition = virtualDomainOrigin[i] + Particles[p].Motion.GetPosition();
                        if (Particles[p].Motion.IsInsideOfPeriodicDomain(virtualPosition, 1.5 * (Particles[p].GetLengthScales().Max()))) 
                        {
                            TreeNode particleNodePeriodic = new TreeNode(false, 0, FindSmallestEmptyID(), virtualDomainOrigin[i]);
                            InsertNodeToList(particleNodePeriodic);
                            particleNodePeriodic.ParticleID = p;
                            particleNodePeriodic.MBR = CalculateParticleMBR(Particles[p], Timestep, virtualDomainOrigin[i]);
                            InsertParticle(particleNodePeriodic, Nodes[0]);
                        }

                    }
                }
            }
        }

        /// <summary>
        /// Updates the MBRs of the tree. Note: this does not change the node structure, i.e. each child-node keeps its parent.
        /// To ensure good performance the tree should be reinitialized every 10-50 time-steps. 
        /// </summary>
        /// <param name="Particles"></param>
        public void UpdateTree(Particle[] Particles, double Timestep) {
            using (new FuncTrace()) {
                for (int i = 0; i < Nodes.Count; i++) {
                    if (Nodes[i].ParticleID > -1)
                        Nodes[i].MBR = CalculateParticleMBR(Particles[Nodes[i].ParticleID], Timestep, Nodes[i].Origin);
                }
                for (int i = 0; i < Nodes.Count; i++) {
                    if (Nodes[i].IsLeaf)
                        UpdateMBR(Nodes[i]);
                }
            }
        }

        /// <summary>
        /// Insert node to the Node-list at the correct position.
        /// </summary>
        /// <param name="CurrentNode"></param>
        private void InsertNodeToList(TreeNode CurrentNode) {
            if (CurrentNode.NodeID < Nodes.Count)
                Nodes.Insert(CurrentNode.NodeID, CurrentNode);
            else
                Nodes.Add(CurrentNode);
        }

        /// <summary>
        /// Search for overlapping MBRs of particles. If two particle MBRs overlap we should check for collisions.
        /// </summary>
        /// <param name="Particle">The current particle</param>
        /// <param name="ParticleID">The position of the particle in the particle-array in the main solver.</param>
        /// <param name="Timestep">Timestep size.</param>
        /// <returns>A list of all particle IDs with an overlapping MBR</returns>
        public List<int> SearchForOverlap(Particle Particle, int ParticleID, double Timestep) {
            using (new FuncTrace()) {
                //es kann passieren, dass manche Überlappungen mehrmals gefunden werden, z.B. weil zwei höhere Knoten mit zwei Partikeln überlappen und dann beide Pfade untersucht werden.
                List<Vector> virtualDomainOrigin = Particle.Motion.OriginInVirtualPeriodicDomain;
                MinimalBoundingRectangle particleMBR = CalculateParticleMBR(Particle, Timestep, new Vector(SpatialDimension));
                List<int> overlappingParticles = new List<int>();
                for (int i = 0; i < Nodes[0].Children.Count; i++) {
                    if (particleMBR.Intersect(Nodes[0].Children[i].MBR).Volume > 0) {
                        if (Nodes[0].Children[i].ParticleID > -1 && Nodes[0].Children[i].ParticleID != ParticleID) {
                            overlappingParticles.Add(Nodes[0].Children[i].ParticleID);
                        } else {
                            overlappingParticles.AddRange(SearchForOverlapRecursive(particleMBR, Nodes[0].Children[i], ParticleID));
                        }
                    }
                }
                for (int j = 0; j < virtualDomainOrigin.Count(); j++) {
                    MinimalBoundingRectangle particleMBRPeriodic = CalculateParticleMBR(Particle, Timestep, virtualDomainOrigin[j]);
                    for (int i = 0; i < Nodes[0].Children.Count; i++) {
                        if (particleMBRPeriodic.Intersect(Nodes[0].Children[i].MBR).Volume > 0) {
                            if (Nodes[0].Children[i].ParticleID > -1 && Nodes[0].Children[i].ParticleID != ParticleID) {
                                overlappingParticles.Add(Nodes[0].Children[i].ParticleID);
                            } else {
                                overlappingParticles.AddRange(SearchForOverlapRecursive(particleMBRPeriodic, Nodes[0].Children[i], ParticleID));
                            }
                        }
                    }
                }
                return overlappingParticles;
            }
        }

        /// <summary>
        /// Recursive search for overlapping particle MBRs
        /// </summary>
        /// <param name="ParticleMBR"></param>
        /// <param name="CurrentNode"></param>
        /// <param name="ParticleID"></param>
        /// <returns></returns>
        private List<int> SearchForOverlapRecursive(MinimalBoundingRectangle ParticleMBR, TreeNode CurrentNode, int ParticleID) {
            List<int> overlappingParticles = new List<int>();
            for (int i = 0; i < CurrentNode.Children.Count; i++) {
                if (ParticleMBR.Intersect(CurrentNode.Children[i].MBR).Volume > 0) {
                    if (CurrentNode.Children[i].ParticleID > -1 && CurrentNode.Children[i].ParticleID != ParticleID) {
                        overlappingParticles.Add(CurrentNode.Children[i].ParticleID);
                    }
                    else
                        overlappingParticles.AddRange(SearchForOverlapRecursive(ParticleMBR, CurrentNode.Children[i], ParticleID));
                }
            }
            return overlappingParticles;
        }

        /// <summary>
        /// Insert particle node to the tree. Two criteria for insertion are followed:
        /// 1. Minimize overlap between different MBRs
        /// 2. Fill up each node to <see cref="MaxEntriesPerNode"/>
        /// </summary>
        /// <param name="ParticleNode">A node containing a single particle (and nothing else).</param>
        /// <param name="CurrentNode">The current node. If it is a leaf the <paramref name="ParticleNode"/> will be added to this node,
        /// otherwise the children are considered.</param>
        private void InsertParticle(TreeNode ParticleNode, TreeNode CurrentNode) {
            if (CurrentNode.IsLeaf) {
                ParticleNode.ParentNodeID = CurrentNode.NodeID;
                CurrentNode.Children.Add(ParticleNode);
                UpdateMBR(CurrentNode);
                if (CurrentNode.Children.Count > MaxEntriesPerNode)
                    Split(CurrentNode);
            } else {
                double minIntersect = double.MaxValue;
                int minIntersectEntry = 0;
                for (int i = 0; i < CurrentNode.Children.Count; i++) {
                    List<TreeNode> temp = new List<TreeNode>();
                    if (!CurrentNode.Children.IsNullOrEmpty()) {
                        temp.Add(CurrentNode.Children[0]);
                        for (int j = 1; j < CurrentNode.Children.Count; j++) {
                            if (!CurrentNode.Children[j].Children.IsNullOrEmpty() && CurrentNode.Children[j].Children.Count < temp[j - 1].Children.Count)
                                temp.Insert(j - 1, CurrentNode.Children[j]);
                            else
                                temp.Add(CurrentNode.Children[j]);
                        }
                        CurrentNode.Children.Clear();
                        CurrentNode.Children.AddRange(temp);
                    }
                    TreeNode child = CurrentNode.Children[i];
                    child.TempMBR = child.MBR.Sum(ParticleNode.MBR);
                    double tempIntersect = 0;
                    for(int j = 0; j < CurrentNode.Children.Count; j++) {
                        if (j != i) 
                            tempIntersect += child.TempMBR.Intersect(CurrentNode.Children[j].MBR).Volume;
                    }
                    if(minIntersect > tempIntersect) {
                        minIntersect = tempIntersect;
                        minIntersectEntry = i;
                    }
                }
                CurrentNode.Children[minIntersectEntry].MBR = CurrentNode.Children[minIntersectEntry].TempMBR;
                InsertParticle(ParticleNode, CurrentNode.Children[minIntersectEntry]);
                if (CurrentNode.Children.Count > MaxEntriesPerNode)
                    Split(CurrentNode);
            }
        }

        /// <summary>
        /// Split <paramref name="Node"/> into two new nodes, because <see cref="MaxEntriesPerNode"/> was reached.
        /// Note: if the root node is split the tree size increases, because the two new nodes will also get a new parent node, which is then the new root.
        /// </summary>
        /// <param name="Node">The node to be split.</param>
        private void Split(TreeNode Node) {
            Vector diff = Node.MBR.LeftUp - Node.MBR.RightLow;
            int longestDimension = 0;
            for(int d = 0; d < SpatialDimension; d++) {
                diff[d] = diff[d].Abs();
            }
            for(int d = 0; d < SpatialDimension; d++) {
                if (diff[d] == diff.Max()) {
                    longestDimension = d;
                    break;
                }
            }

            List<TreeNode> temp = new List<TreeNode>();
            temp.Add(Node.Children[0]);
            for (int i = 1; i < Node.Children.Count; i++) {
                if (Node.Children[i].MBR.LeftUp[longestDimension] < temp[i - 1].MBR.LeftUp[longestDimension])
                    temp.Insert(i - 1, Node.Children[i]);
                else
                    temp.Add(Node.Children[i]);
            }
            Node.Children.Clear();

            TreeNode NewNode = new TreeNode(Node.IsLeaf, Node.ParentNodeID, FindSmallestEmptyID(), new Vector(SpatialDimension));
            InsertNodeToList(NewNode);
            InsertNode(temp[0], Node);
            InsertNode(temp.Last(), NewNode);
            if (NewNode.ParentNodeID == -1) {
                // highest Node, we need to create a new parent (tree size increase)
                TreeNode NewParentNode = new TreeNode(false, -1, Node.NodeID, new Vector(SpatialDimension));
                Node.NodeID = FindSmallestEmptyID();
                Nodes.RemoveAt(0);
                InsertNodeToList(Node);
                Nodes.Insert(0, NewParentNode);
                for (int i = 0; i < Node.Children.Count; i++) {
                    Node.Children[i].ParentNodeID = Node.NodeID;
                }
                InsertNode(Node, NewParentNode);
                InsertNode(NewNode, NewParentNode);
            } else {
                for (int i = 0; i < Nodes.Count; i++) {
                    if (NewNode.ParentNodeID == Nodes[i].NodeID)
                        InsertNode(NewNode, Nodes[i]);
                }
            }

            for (int i = 1; i < temp.Count - 1; i++) {
                for (int j = 0; j < Nodes.Count; j++) {
                    if (Node.ParentNodeID == Nodes[j].NodeID)
                        InsertNodeMinIntersect(temp[i], Nodes[j]);
                }
            }

        }

        /// <summary>
        /// Insert <paramref name="NewNode"/> into <paramref name="ParentNode"/>.
        /// </summary>
        /// <param name="NewNode"></param>
        /// <param name="ParentNode"></param>
        private void InsertNode(TreeNode NewNode, TreeNode ParentNode) {
            NewNode.ParentNodeID = ParentNode.NodeID;
            ParentNode.Children.Add(NewNode);
            UpdateMBR(ParentNode);
            if (ParentNode.Children.Count > MaxEntriesPerNode)
                Split(ParentNode);
        }

        /// <summary>
        /// Insert <paramref name="NewNode"/> into one of the child nodes of <paramref name="ParentNode"/>, based on the minimizing overlap strategy.
        /// </summary>
        /// <param name="NewNode"></param>
        /// <param name="ParentNode"></param>
        private void InsertNodeMinIntersect(TreeNode NewNode, TreeNode ParentNode) {
            double minIntersect = double.MaxValue;
            int minIntersectEntry = 0;
            for (int i = 0; i < ParentNode.Children.Count; i++) {
                TreeNode currentNode = ParentNode.Children[i];
                currentNode.TempMBR = currentNode.MBR.Sum(NewNode.MBR);
                double tempIntersect = 0;
                for (int j = 0; j < ParentNode.Children.Count; j++) {
                    if (j != i)
                        tempIntersect += currentNode.TempMBR.Intersect(ParentNode.Children[j].MBR).Volume;
                }
                if (minIntersect > tempIntersect) {
                    minIntersect = tempIntersect;
                    minIntersectEntry = i;
                } else if(minIntersect == tempIntersect && currentNode.Children.Count < ParentNode.Children[minIntersectEntry].Children.Count) 
                    minIntersectEntry = i;
            }
            InsertNode(NewNode, ParentNode.Children[minIntersectEntry]);
        }

        /// <summary>
        /// Removes <paramref name="NodeToRemove"/> and updates all higher nodes.
        /// </summary>
        /// <param name="NodeToRemove"></param>
        private void RemoveNode(TreeNode NodeToRemove) {
            for (int i = Nodes.Count - 1; i >= 0; i--) {
                if (NodeToRemove.NodeID == Nodes[i].NodeID) {
                    Nodes.RemoveAt(i);
                    break;
                } 
            }
            for(int i = 0; i < Nodes.Count; i++) {
                if(NodeToRemove.ParentNodeID == Nodes[i].NodeID) {
                    for(int j = 0; j < Nodes[i].Children.Count; j++) {
                        if(Nodes[i].Children[j].NodeID == NodeToRemove.NodeID) {
                            Nodes[i].Children.RemoveAt(j);
                            if (Nodes[i].Children.Count < 1) {
                                RemoveNode(Nodes[i]);
                            } else {
                                UpdateMBR(Nodes[i]);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Update the MBR of <paramref name="CurrentNode"/> and its parents.
        /// </summary>
        /// <param name="CurrentNode"></param>
        private void UpdateMBR(TreeNode CurrentNode) {
            CurrentNode.MBR = new MinimalBoundingRectangle(CurrentNode.Children[0].MBR.LeftUp, CurrentNode.Children[0].MBR.RightLow);
            for (int i = 1; i < CurrentNode.Children.Count; i++) {
                CurrentNode.MBR = CurrentNode.MBR.Sum(CurrentNode.Children[i].MBR);
            }
            if(CurrentNode.ParentNodeID > -1)
                UpdateMBR(Nodes[CurrentNode.ParentNodeID]);
        }

        /// <summary>
        /// Finds the smallest available node ID.
        /// </summary>
        /// <returns></returns>
        private int FindSmallestEmptyID() {
            for (int i = 1; i < Nodes.Count; i++) {
                if (Nodes[i].NodeID - Nodes[i - 1].NodeID > 1) {
                    return Nodes[i - 1].NodeID + 1;
                }
            }
            return Nodes.Count;
        }

        /// <summary>
        /// MBR = minimal bounding rectangle
        /// </summary>
        private MinimalBoundingRectangle CalculateParticleMBR(Particle Particle, double Timestep, Vector Origin) {
            Vector[] vertices = new Vector[2];
            for (int i = 0; i < 2; i++)
                vertices[i] = new Vector(SpatialDimension);
            int noOfSubParticles = Particle.NoOfSubParticles;
            int noOfTimesteps = 2;
            MinimalBoundingRectangle mbr = new(new Vector(SpatialDimension), new Vector(SpatialDimension));
            for (int j = 0; j < noOfTimesteps; j++) {
                for (int i = 0; i < noOfSubParticles; i++) {
                    if (i > 0)
                        throw new Exception("No support for sub-particles");
                    Vector Position = j == 0 ? new Vector(Particle.Motion.GetPosition(0)) : new Vector(PredictParticePositionNextTimestep(Particle, Timestep));
                    Position += Origin;
                    Vector Angle = j == 0 ? new Vector(Particle.Motion.GetAngle(0)) : PredictParticleAngleNextTimestep(Particle,Timestep);
                    Vector LeftUp = new(SpatialDimension);
                    Vector RightLow = new(SpatialDimension);
                    for (int d = 0; d < SpatialDimension; d++) {
                        double sign = d == 0 ? -1 : 1;
                        Vector supportVector = new(SpatialDimension);
                        supportVector[d] = sign;
                        LeftUp[d] = Particle.GetSupportPoint(supportVector, Position, Angle, i)[d] + sign * Tolerance;
                        supportVector[d] = -sign;
                        RightLow[d] = Particle.GetSupportPoint(supportVector, Position, Angle, i)[d] - sign * Tolerance;
                    }
                    if (j == 0)
                        mbr = new MinimalBoundingRectangle(LeftUp, RightLow);
                    else
                        mbr = mbr.Sum(new MinimalBoundingRectangle(LeftUp, RightLow));
                }
            }
            return mbr;
        }

        private static Vector PredictParticePositionNextTimestep(Particle Particle, double Timestep) {
            return Particle.Motion.GetPosition(0) + (Particle.Motion.GetTranslationalVelocity(0) + 4 * Particle.Motion.GetTranslationalVelocity(1) + Particle.Motion.GetTranslationalVelocity(2)) * Timestep / 3;
        }

        private static Vector PredictParticleAngleNextTimestep(Particle Particle, double Timestep) {
            return new Vector(Particle.Motion.GetAngle(0) + (Particle.Motion.GetRotationalVelocity(0) + 4 * Particle.Motion.GetRotationalVelocity(1) + Particle.Motion.GetRotationalVelocity(2)) * Timestep / 3);
        }
    }

    class TreeNode {
        public bool IsLeaf;
        public List<TreeNode> Children;
        public int ParticleID;
        public int ParentNodeID;
        public int NodeID;
        public MinimalBoundingRectangle MBR;
        public MinimalBoundingRectangle TempMBR;
        public Vector Origin;

        public TreeNode(bool IsLeaf, int ParentNodeID, int NodeID, Vector Origin) {
            this.IsLeaf = IsLeaf;
            this.ParentNodeID = ParentNodeID;
            this.ParticleID = -1;
            this.NodeID = NodeID;
            Children = new List<TreeNode>();
            this.Origin = Origin;
        }

        public void AddChild(TreeNode Child) {
            Children.Add(Child);
        }

    }

    class MinimalBoundingRectangle {
        public Vector LeftUp = new Vector(3);
        public Vector RightLow = new Vector(3);
        public double Volume;
        public int SpatialDimension;

        public MinimalBoundingRectangle(Vector LeftUp, Vector RightLow) {
            this.LeftUp = new Vector(LeftUp);
            this.RightLow = new Vector(RightLow);
            this.SpatialDimension = LeftUp.Dim;
            Vector diff = this.LeftUp - this.RightLow;
            if (SpatialDimension == 2)
                this.Volume = diff[0].Abs() * diff[1].Abs();
            else if (SpatialDimension == 3)
                this.Volume = diff[0].Abs() * diff[1].Abs() * diff[2].Abs();
            else
                throw new Exception("Dimension mismatch");
        }

        public MinimalBoundingRectangle Sum(MinimalBoundingRectangle AdditionalMBR) {
            Vector sumLeftUp = new Vector(LeftUp);
            Vector sumRightLow = new Vector(RightLow);
            if (LeftUp[0] > AdditionalMBR.LeftUp[0])
                sumLeftUp[0] = AdditionalMBR.LeftUp[0];
            if (RightLow[0] < AdditionalMBR.RightLow[0])
                sumRightLow[0] = AdditionalMBR.RightLow[0];
            for (int i = 1; i < LeftUp.Count; i++) {
                if (LeftUp[i] < AdditionalMBR.LeftUp[i])
                    sumLeftUp[i] = AdditionalMBR.LeftUp[i];
                if (RightLow[i] > AdditionalMBR.RightLow[i])
                    sumRightLow[i] = AdditionalMBR.RightLow[i];
            }
            return new MinimalBoundingRectangle(sumLeftUp, sumRightLow);
        }

        public MinimalBoundingRectangle Intersect(MinimalBoundingRectangle MBR) {
            Vector intersectLeftUp = new Vector(LeftUp.Dim);
            Vector intersectRightLow = new Vector(RightLow.Dim);
            double dx = Math.Min(RightLow[0], MBR.RightLow[0]) - Math.Max(LeftUp[0], MBR.LeftUp[0]);
            double dy = Math.Min(LeftUp[1], MBR.LeftUp[1]) - Math.Max(RightLow[1], MBR.RightLow[1]);
            double dz = LeftUp.Dim == 3 ? Math.Min(LeftUp[1], MBR.LeftUp[1]) - Math.Max(RightLow[1], MBR.RightLow[1]) : 0;
            if(dx >= 0 && dy >= 0 && dz >= 0) {
                intersectLeftUp[0] = Math.Max(LeftUp[0], MBR.LeftUp[0]);
                intersectLeftUp[1] = Math.Min(LeftUp[1], MBR.LeftUp[1]);
                if(intersectLeftUp.Dim == 3)
                    intersectLeftUp[2] = Math.Min(LeftUp[2], MBR.LeftUp[2]);

                intersectRightLow[0] = Math.Min(RightLow[0], MBR.RightLow[0]);
                intersectRightLow[1] = Math.Max(RightLow[1], MBR.RightLow[1]);
                if (intersectRightLow.Dim == 3)
                    intersectRightLow[2] = Math.Max(RightLow[2], MBR.RightLow[2]);
            }
            return new MinimalBoundingRectangle(intersectLeftUp, intersectRightLow);
        }
    }
}
