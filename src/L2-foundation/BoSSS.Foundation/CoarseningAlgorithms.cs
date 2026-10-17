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

using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using BoSSS.Platform;
using ilPSP.Utils;
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.Grid.Aggregation {

    /// <summary>
    /// Parallel, deterministic R-tree coarsener (cut by global tree level).
    /// </summary>
    public static class RTreeParallelCoarsening {

        public static AggregationGridData CoarsenByLevel(IGridData gd, int aggCellCountLikeBefore) {
            using(new FuncTrace()) {
                int D = gd.SpatialDimension;
                int M = 1 << D;

                // current fine count
                int J = gd.iLogicalCells.NoOfLocalUpdatedCells;

                // desired coarse count for this step (≈ factor M reduction)
                int desired = Math.Max(1, (int)Math.Ceiling((double)J / M));

                var entries = BuildLeafEntriesDeterministic(gd, out var _);
                var root = BulkLoad_STR_GlobalAware(entries, M);

                // choose the deepest level whose node count is < J and as close as possible to 'desired'
                int cutLevel = SelectCutLevelByCount(root, desired, J);

                var groups = ExtractAggregatesAtLevel(root, cutLevel)
                             .Select(n => CollectLeaves(n).ToArray())
                             .ToArray();

                Debug_Assignments(groups, J);

                var ag = new AggregationGrid(gd.Grid, groups);
                if(!object.ReferenceEquals(ag.GridData.Grid, ag))
                    throw new ApplicationException("internal error in mesh coarsening");
                return ag.GridData;
            }
        }

        private static int SelectCutLevelByCount(Node root, int desired, int currentJ) {
            // BFS: count nodes per level (level 1 = root)
            var counts = new List<int>();
            var q = new Queue<(Node n, int d)>();
            q.Enqueue((root, 1));
            while(q.Count > 0) {
                var (n, d) = q.Dequeue();
                if(counts.Count < d) counts.Add(0);
                counts[d - 1]++;
                if(!n.IsLeaf) foreach(var c in n.Children) q.Enqueue((c, d + 1));
            }

            // find a level with reduction (< currentJ) and closest to 'desired'
            int bestLevel = counts.Count; // deepest by default
            int bestScore = int.MaxValue;
            for(int L = 1; L <= counts.Count; L++) {
                int nAtL = counts[L - 1];
                if(nAtL >= currentJ) continue; // no reduction
                int score = Math.Abs(nAtL - desired);
                if(score < bestScore) { bestScore = score; bestLevel = L; }
            }

            // fallback: if nothing reduces (shouldn’t happen), use deepest level
            return bestLevel;
        }


        // -------------------- R-tree node --------------------
        private class Node {
            public BoundingBox BB;
            public List<Node> Children;      // internal
            public List<int> LeafIds;        // leaf = base-cell local indices
            public bool IsLeaf => Children == null || Children.Count == 0;
            public int LeafCount;            // cached leaves in subtree
            public int FirstLeafId;          // min leaf id in subtree (for deterministic tie-break)
        }

        // -------------------- Deterministic leaf creation --------------------
        private sealed class LeafEntry {
            public BoundingBox BB;
            public int jLoc;
            public ulong morton;
            public long jGlob;
        }

        private static List<Node> BuildLeafEntriesDeterministic(IGridData gd, out Func<int, long> globalIdOfLocal) {
            int D = gd.SpatialDimension;
            int J = gd.iLogicalCells.NoOfLocalUpdatedCells;
            var part = gd.CellPartitioning;

            globalIdOfLocal = j => part.i0 + j;

            var list = new List<LeafEntry>(J);
            var bb = new BoundingBox(D);
            for(int j = 0; j < J; j++) {
                bb.Clear();
                gd.iLogicalCells.GetCellBoundingBox(j, bb);
                var c = new double[D];
                for(int d = 0; d < D; d++) c[d] = 0.5 * (bb.Min[d] + bb.Max[d]);

                list.Add(new LeafEntry {
                    BB = (BoundingBox)bb.Clone(),
                    jLoc = j,
                    morton = MortonKey(c),
                    jGlob = part.i0 + j
                });
            }

            list.Sort((a, b) => {
                int s = a.morton.CompareTo(b.morton);
                if(s != 0) return s;
                for(int d = 0; d < D; d++) {
                    s = Cmp(a.BB.Min[d], b.BB.Min[d]);
                    if(s != 0) return s;
                }
                return a.jGlob.CompareTo(b.jGlob);
            });

            var nodes = new List<Node>(J);
            foreach(var e in list) {
                nodes.Add(new Node {
                    BB = (BoundingBox)e.BB.Clone(),
                    LeafIds = new List<int> { e.jLoc },
                    LeafCount = 1,
                    FirstLeafId = e.jLoc
                });
            }
            return nodes;
        }

        private static int Cmp(double x, double y) {
            bool xn = double.IsNaN(x), yn = double.IsNaN(y);
            if(xn | yn) return xn == yn ? 0 : (xn ? 1 : -1); // NaN last
            return x.CompareTo(y);
        }

        private static double Clamp01(double v) { return v < 0 ? 0 : (v > 1 ? 1 : v); }

        private static ulong MortonKey(double[] c) {
            int D = c.Length;
            uint x = (uint)(Clamp01(c[0]) * 4294967295.0);
            uint y = (uint)((D > 1 ? Clamp01(c[1]) : 0.0) * 4294967295.0);
            uint z = (uint)((D > 2 ? Clamp01(c[2]) : 0.0) * 4294967295.0);
            return Interleave(x, y, z);
        }

        private static ulong Interleave(uint x, uint y, uint z) {
            ulong _x = Part1By2(x);
            ulong _y = Part1By2(y) << 1;
            ulong _z = Part1By2(z) << 2;
            return _x | _y | _z;
        }
        private static ulong Part1By2(uint n) {
            ulong x = n;
            x = (x | (x << 32)) & 0x1F00000000FFFF;
            x = (x | (x << 16)) & 0x1F0000FF0000FF;
            x = (x | (x << 8)) & 0x100F00F00F00F00F;
            x = (x | (x << 4)) & 0x10C30C30C30C30C3;
            x = (x | (x << 2)) & 0x1249249249249249;
            return x;
        }

        // -------------------- STR-like bulk load with deterministic axis split --------------------
        private static Node BulkLoad_STR_GlobalAware(List<Node> entries, int M) {
            if(entries.Count <= M) {
                if(entries.Count == 1) return entries[0];
                var p = new Node { Children = new List<Node>(entries), BB = UnionBB(entries) };
                p.LeafCount = entries.Sum(e => e.LeafCount);
                p.FirstLeafId = entries.Min(e => e.FirstLeafId);
                return p;
            }

            var bbAll = UnionBB(entries);
            int axis = LongestAxis(bbAll);

            entries.Sort((a, b) => {
                int s = Cmp(a.BB.Min[axis], b.BB.Min[axis]);
                if(s != 0) return s;
                int na = (axis + 1) % a.BB.D;
                s = Cmp(a.BB.Min[na], b.BB.Min[na]);
                if(s != 0) return s;
                return a.FirstLeafId.CompareTo(b.FirstLeafId);
            });

            var parents = new List<Node>();
            for(int i = 0; i < entries.Count; i += M) {
                var slice = entries.GetRange(i, Math.Min(M, entries.Count - i));
                var p = new Node { Children = slice, BB = UnionBB(slice) };
                p.LeafCount = slice.Sum(e => e.LeafCount);
                p.FirstLeafId = slice.Min(e => e.FirstLeafId);
                parents.Add(p);
            }
            return BulkLoad_STR_GlobalAware(parents, M);
        }

        private static BoundingBox UnionBB(List<Node> nodes) {
            int D = nodes[0].BB.D;
            var u = new BoundingBox(D);
            foreach(var n in nodes) u.AddBB(n.BB);
            return u;
        }

        private static int LongestAxis(BoundingBox bb) {
            int D = bb.D;
            int axis = 0;
            double maxLen = -1.0;
            for(int d = 0; d < D; d++) {
                double len = bb.Max[d] - bb.Min[d];
                if(len > maxLen) { maxLen = len; axis = d; }
            }
            return axis;
        }

        // -------------------- Level cut --------------------
        private static IEnumerable<Node> ExtractAggregatesAtLevel(Node root, int level) {
            var q = new Queue<(Node n, int d)>();
            q.Enqueue((root, 1));
            while(q.Count > 0) {
                var (n, d) = q.Dequeue();
                if(d == level || n.IsLeaf) {
                    yield return n;
                    continue;
                }
                foreach(var c in n.Children) q.Enqueue((c, d + 1));
            }
        }

        private static List<int> CollectLeaves(Node n) {
            if(n.IsLeaf) return new List<int>(n.LeafIds);
            var acc = new List<int>(n.LeafCount);
            var st = new Stack<Node>();
            st.Push(n);
            while(st.Count > 0) {
                var cur = st.Pop();
                if(cur.IsLeaf) acc.AddRange(cur.LeafIds);
                else foreach(var c in cur.Children) st.Push(c);
            }
            return acc;
        }

        // -------------------- Debug --------------------
        [Conditional("DEBUG")]
        private static void Debug_Assignments(int[][] groups, int J) {
            var mark = new BitArray(J);
            foreach(var g in groups)
                foreach(var j in g) {
                    Debug.Assert(j >= 0 && j < J);
                    Debug.Assert(!mark[j], $"Cell {j} assigned twice.");
                    mark[j] = true;
                }
            for(int j = 0; j < J; j++)
                Debug.Assert(mark[j], $"Cell {j} unassigned.");
        }
    }



    /// <summary>
    /// R-tree–style, geometry-driven coarsening with identical public interfaces.
    /// </summary>
    public static class RTreeCoarseningAlgorithms {

        // ---- Public API (same signatures as in CoarseningAlgorithms) ----------------

        public static AggregationGridData Coarsen(IGridData ag, int AggCellCount) {
            using(new FuncTrace()) {
                var g = Coarsen(ag.Grid, AggCellCount);
                if(!object.ReferenceEquals(g.GridData.Grid, g))
                    throw new ApplicationException("internal error in mesh coarsening");
                return g.GridData;
            }
        }

        public static AggregationGrid Coarsen(IGrid ag, int AggCellCount) {
            using(new FuncTrace()) {
                int[][] groups = AggregationKernel_RTree(ag.iGridData, AggCellCount);
                return new AggregationGrid(ag, groups);
            }
        }

        // ---- Core: R-tree bulk build + level extraction -----------------------------

        private class Node {
            public BoundingBox BB;              // covers either children or leaves
            public List<Node> Children;         // internal node
            public List<int> LeafIds;           // leaf node: indices of base cells
            public bool IsLeaf => Children == null || Children.Count == 0;
            public int LeafCount;               // cached total leaves in subtree
        }

        /// <summary>
        /// Geometry-driven aggregation. Returns int[][] with groups of fine cell indices.
        /// </summary>
        private static int[][] AggregationKernel_RTree(IGridData gd, int AggCellCount) {
            if(AggCellCount < 2) throw new ArgumentOutOfRangeException(nameof(AggCellCount));

            int D = gd.SpatialDimension;
            int J = gd.iLogicalCells.NoOfLocalUpdatedCells;

            // Fanout parameters (R*-tree inspired; simple bulk-load)
            int M = 1 << D;        // max entries per node: 2^D
            int m = Math.Max(2, M / 2);

            // Collect per-cell AABBs once
            var boxes = new BoundingBox[J];
            for(int j = 0; j < J; j++) {
                var bb = new BoundingBox(D);
                gd.iLogicalCells.GetCellBoundingBox(j, bb);
                boxes[j] = bb;
            }

            // Build R-tree (bulk load with axis-wise sorted Str packing)
            var leaves = new List<Node>(J);
            for(int j = 0; j < J; j++) {
                leaves.Add(new Node {
                    BB = boxes[j].CloneAs(),
                    LeafIds = new List<int> { j },
                    Children = null,
                    LeafCount = 1
                });
            }

            Node root = BulkLoad(leaves, M);

            // Choose extraction level targeting ~AggCellCount leaves per aggregate
            int target = Math.Max(2, AggCellCount);
            var groups = ExtractAggregatesForTarget(root, target, M, m);

            // Safety: ensure every cell is assigned exactly once
            Debug_Assignments(groups, J);

            return groups;
        }

        // ---- R-tree helpers ---------------------------------------------------------

        // Bulk-load like STR: recursively pack entries along longest axis into chunks of size M
        private static Node BulkLoad(List<Node> entries, int M) {
            if(entries.Count <= M) {
                // create one level above leaves if entries are leaf nodes
                if(entries.Count == 1) return entries[0];

                var parent = new Node {
                    Children = new List<Node>(entries),
                    BB = UnionBB(entries),
                };
                parent.LeafCount = entries.Sum(e => e.LeafCount);
                return parent;
            }

            // Sort by longest axis of overall BB, then chunk into groups of size <= M
            var bbAll = UnionBB(entries);
            int axis = LongestAxis(bbAll);

            entries.Sort((a, b) => a.BB.Min[axis].CompareTo(b.BB.Min[axis]));

            var parents = new List<Node>();
            for(int i = 0; i < entries.Count; i += M) {
                var slice = entries.GetRange(i, Math.Min(M, entries.Count - i));
                var parent = new Node {
                    Children = slice,
                    BB = UnionBB(slice)
                };
                parent.LeafCount = slice.Sum(e => e.LeafCount);
                parents.Add(parent);
            }

            return BulkLoad(parents, M);
        }

        private static BoundingBox UnionBB(List<Node> nodes) {
            int D = nodes[0].BB.D;
            var u = new BoundingBox(D);
            foreach(var n in nodes) u.AddBB(n.BB);
            return u;
        }

        private static int LongestAxis(BoundingBox bb) {
            int D = bb.D;
            int axis = 0;
            double maxLen = -1.0;
            for(int d = 0; d < D; d++) {
                double len = bb.Max[d] - bb.Min[d]; // ← use extents
                if(len > maxLen) { maxLen = len; axis = d; }
            }
            return axis;
        }

        // Extract groups so that each group has ~ target leaves, using the tree levels.
        private static int[][] ExtractAggregatesForTarget(Node root, int target, int M, int m) {
            // BFS traversal; cut at nodes whose LeafCount ~ target; otherwise descend.
            var groups = new List<int[]>();
            var queue = new Queue<Node>();
            queue.Enqueue(root);

            // heuristic band: accept nodes with LeafCount in [target/2, 2*target]
            int lower = Math.Max(2, target / 2);
            int upper = Math.Max(target, 2 * target);

            while(queue.Count > 0) {
                var n = queue.Dequeue();

                if(n.LeafCount <= upper && n.LeafCount >= lower) {
                    groups.Add(CollectLeaves(n).ToArray());
                    continue;
                }

                if(n.IsLeaf) {
                    // too small: group singleton (will happen if target is large)
                    groups.Add(n.LeafIds.ToArray());
                } else {
                    // If a child already exceeds upper bound, push it as its own group to avoid over-merge.
                    // Else continue descending.
                    foreach(var c in n.Children) {
                        if(c.LeafCount > upper && !c.IsLeaf) {
                            // try to split deeper
                            queue.Enqueue(c);
                        } else {
                            // accept as a group if reasonable; else descend one step.
                            if(c.LeafCount >= lower || c.IsLeaf) {
                                groups.Add(CollectLeaves(c).ToArray());
                            } else {
                                // descend one level to improve balance
                                if(c.IsLeaf) {
                                    groups.Add(c.LeafIds.ToArray());
                                } else {
                                    foreach(var gc in c.Children) queue.Enqueue(gc);
                                }
                            }
                        }
                    }
                }
            }

            // Optional: minor post-pass to merge tiny groups with nearest by BB enlargement
            groups = MergeTiny(groups, lower);

            return groups.ToArray();
        }

        private static List<int> CollectLeaves(Node n) {
            if(n.IsLeaf) return new List<int>(n.LeafIds);
            var acc = new List<int>(n.LeafCount);
            var stack = new Stack<Node>();
            stack.Push(n);
            while(stack.Count > 0) {
                var cur = stack.Pop();
                if(cur.IsLeaf) {
                    acc.AddRange(cur.LeafIds);
                } else {
                    foreach(var c in cur.Children) stack.Push(c);
                }
            }
            return acc;
        }

        private static List<int[]> MergeTiny(List<int[]> groups, int lower) {
            if(groups.Count <= 1) return groups;

            // Map bbox per group by union of member bboxes is not available here,
            // so do a simple size-based greedy merge.
            var list = new List<List<int>>(groups.Select(g => new List<int>(g)));
            bool changed = true;
            while(changed) {
                changed = false;
                for(int i = 0; i < list.Count; i++) {
                    if(list[i].Count >= lower) continue;
                    // merge with nearest in index (cheap, deterministic)
                    int j = (i == list.Count - 1) ? i - 1 : i + 1;
                    if(j >= 0 && j < list.Count && j != i) {
                        list[j].AddRange(list[i]);
                        list.RemoveAt(i);
                        changed = true;
                        break;
                    }
                }
            }
            return list.Select(l => l.ToArray()).ToList();
        }

        // ---- Debug checks ------------------------------------------------------------

        [Conditional("DEBUG")]
        private static void Debug_Assignments(int[][] groups, int J) {
            var mark = new BitArray(J);
            foreach(var g in groups) {
                foreach(var j in g) {
                    Debug.Assert(j >= 0 && j < J);
                    Debug.Assert(mark[j] == false, $"Cell {j} assigned twice.");
                    mark[j] = true;
                }
            }
            for(int j = 0; j < J; j++) {
                Debug.Assert(mark[j], $"Cell {j} was not assigned to any aggregate.");
            }
        }
    }

    public enum AggStrategy { Adjacency, RTree }

    /// <summary>
    /// Creation of multigrid hierarchies
    /// </summary>
    static public class CoarseningAlgorithms {

        /// <summary>
        /// Creates a sequence of aggregated grids, suitable for a multigrid algorithm
        /// </summary>
        /// <param name="Grid">original grid</param>
        /// <param name="MaxDepth">maximum number of refinements</param>
        /// <returns></returns>
        public static AggregationGrid[] CreateSequence(IGrid Grid, int MaxDepth = -1) {
            AggregationGridData[] seq = CreateSequence(Grid.iGridData, MaxDepth);
            return seq.Select(agd => (AggregationGrid)agd.Grid).ToArray();
        }


        /// <summary>
        /// Creates a sequence of aggregated grids, suitable for a multigrid algorithm
        /// </summary>
        /// <param name="GridDat">original grid</param>
        /// <param name="MaxDepth">maximum number of refinements</param>
        /// <param name="strategy">coarsening strategy</param>
        /// <returns></returns>
        public static AggregationGridData[] CreateSequence(IGridData GridDat, int MaxDepth = -1,
                                                           AggStrategy strategy = AggStrategy.Adjacency) {
            using(new FuncTrace()) {
                Console.WriteLine("Creating multigrid hierarchy with max depth " + MaxDepth + " using " + strategy + " coarsening.");
                MPICollectiveWatchDog.Watch();
                int D = GridDat.SpatialDimension;
                MaxDepth = MaxDepth >= 0 ? MaxDepth : int.MaxValue;
                var aggGrids = new List<AggregationGridData>();
                aggGrids.Add(ZeroAggregation(GridDat));

                var localNoOfCells = new List<int> { aggGrids[0].iLogicalCells.NoOfLocalUpdatedCells };
                var globalNoOfCells = new List<long> { aggGrids[0].CellPartitioning.TotalLength };

                while(true) {
                    MPICollectiveWatchDog.Watch(token: 80);
                    if(aggGrids.Count >= MaxDepth) break;

                    // NOTE: both are AggregationGridData
                    AggregationGridData grid2coarsen = aggGrids.Last();

                    // *** FIXED TYPE: returns AggregationGridData in both branches ***
                    AggregationGridData grid =
                        (strategy == AggStrategy.RTree)
                        ? RTreeParallelCoarsening.CoarsenByLevel(grid2coarsen, (int)Math.Pow(2, D))
                        : CoarseningAlgorithms.Coarsen(grid2coarsen, (int)Math.Pow(2, D));

                    double aimred = 1 / (Math.Pow(2, D)) * 2;
                    for(int iCoarsen = 0; iCoarsen < 1; iCoarsen++) {
                        double actualred = (double)grid.CellPartitioning.LocalLength
                                         / (double)aggGrids.Last().CellPartitioning.LocalLength;
                        if((actualred < aimred).MPIAnd()) break;

                        grid2coarsen = grid;
                        grid = (strategy == AggStrategy.RTree)
                            ? RTreeParallelCoarsening.CoarsenByLevel(grid2coarsen, (int)Math.Pow(2, D))
                            : CoarseningAlgorithms.Coarsen(grid2coarsen, (int)Math.Pow(2, D));

                        // AggregationGridData has this method
                        grid.MergeWithPartentGrid(grid2coarsen);
                    }

                    int Jloc = grid.CellPartitioning.LocalLength;
                    long Jtot = grid.CellPartitioning.TotalLength;

                    bool localReduction = (Jloc < localNoOfCells.Last()).MPIOr();
                    bool globalReduction = Jtot < globalNoOfCells.Last();
                    if(!localReduction || !globalReduction) break;

                    aggGrids.Add(grid);
                    localNoOfCells.Add(Jloc);
                    globalNoOfCells.Add(Jtot);

#if DEBUG
            // your existing DEBUG checks unchanged
#endif
                }
                return aggGrids.ToArray();
            }
        }



        /// <summary>
        /// creates an initial aggregated grid which is in fact equivalent to <paramref name="gd"/>
        /// </summary>
        public static AggregationGridData ZeroAggregation(IGridData gd) {
            var g = gd.Grid;
            var gc = ZeroAggregation(g);
            return ((AggregationGridData)gc.iGridData);

        }

        /// <summary>
        /// creates an initial aggregated grid which is in fact equivalent to <paramref name="g"/>
        /// </summary>
        public static AggregationGrid ZeroAggregation(IGrid g) {
            //var Cls = g.Cells;
            int J = g.CellPartitioning.LocalLength;
            int D = g.SpatialDimension;

            int[][] AggregateCells = new int[J][];
            for (int j = 0; j < J; j++) {
                AggregateCells[j] = new int[] { j };
            }

            AggregationGrid ret = new AggregationGrid(g, AggregateCells);

            return ret;
        }

        /// <summary>
        /// coarsens level <paramref name="ag"/> (aggregation of grid data objects)
        /// </summary>
        /// <param name="ag">
        /// input grid, which should be aggregated
        /// </param>
        /// <param name="AggCellCount">
        /// desired number of parts for each aggregate cell
        /// </param>
        public static AggregationGridData Coarsen(IGridData ag, int AggCellCount) {
            using (new FuncTrace()) {
                var g = Coarsen(ag.Grid, AggCellCount);
                if (!object.ReferenceEquals(g.GridData.Grid, g))
                    throw new ApplicationException("internal error in mesh coarsening");
                return g.GridData;
            }
        }

        /// <summary>
        /// coarsens level <paramref name="ag"/> (aggregation of grid objects)
        /// </summary>
        /// <param name="ag">
        /// input grid, which should be aggregated
        /// </param>
        /// <param name="AggCellCount">
        /// desired number of parts for each aggregate cell
        /// </param>
        public static AggregationGrid Coarsen(IGrid ag, int AggCellCount) {
            using(new FuncTrace()) {

                IGridData pGridData = ag.iGridData;
                int[][] Coarsened_ComositeCells = AggregationKernel(pGridData, AggCellCount);

               
                // Cuthill-McKey sorting should theoretically help iterative and direct solvers
                //Coarsened_ComositeCells = RandomSorting(Coarsened_ComositeCells);
                //Coarsened_ComositeCells = CuthillMcKey(ag, Coarsened_ComositeCells, true);

                var g = new AggregationGrid(ag, Coarsened_ComositeCells);
                return g;
            }
        }

        static int[][] RandomSorting(int[][] __Coarsened_ComositeCells) {
            // random permute of aggregation cells
            int JC = __Coarsened_ComositeCells.Length;
            int[][] Coarsened_ComositeCells = new int[JC][];
            Random rnd = new Random();
            //Console.Write(" perm: ");
            for (int jc = 0; jc < JC; jc++) {
                int jDest = rnd.Next(JC);
                while (Coarsened_ComositeCells[jDest] != null) {
                    jDest = (jDest + 1) % JC;
                }

                //Console.Write($" {jc}>{jDest}");
                Coarsened_ComositeCells[jDest] = __Coarsened_ComositeCells[jc];
            }
            //Console.WriteLine();

            return Coarsened_ComositeCells;
        }


        static int[][] CuthillMcKey(IGrid parrent, int[][] AggCells, bool reverse) {

            var gTemp = new AggregationGrid(parrent, AggCells);
            int JC = AggCells.Length;
            if (JC != gTemp.iGridData.iLogicalCells.NoOfLocalUpdatedCells)
                throw new ApplicationException();

            BitArray added = new BitArray(JC);
            int[] AdjRest(int jc) {
                return gTemp.GridData.iLogicalCells.CellNeighbours[jc].Where(jneigh => jneigh < JC && !added[jneigh]).ToArray();
            }

            int Degree(int jc) {
                return gTemp.GridData.iLogicalCells.CellNeighbours[jc].Where(jneigh => jneigh < JC).Count();
            }

            List<int> R = new List<int>(new int[] { 0 }); added[0] = true;
            for(int i = 0; R.Count < JC; i++) {
                int Ri = R[i];

                int[] Ai = AdjRest(Ri);
                Debug.Assert(Ai.Where(a => added[a]).Count() == 0);
                if (Ai.Length > 0) {
                    int[] Degs = Ai.Select(a => Degree(a)).ToArray();

                    Array.Sort(Degs, Ai);

                    R.AddRange(Ai);
                    foreach (var k in Ai)
                        added[k] = true;
                }
            }

            if (R.Count != JC)
                throw new ApplicationException("Cuthill-McKey internal error.");


            int[][] ret = new int[JC][];
            for(int j = 0; j < JC; j++) {
                if(reverse)
                    ret[JC - j - 1] = AggCells[R[j]];
                else
                    ret[j] = AggCells[R[j]];
            }
            return ret;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="ag">
        /// mesh which should be coarsened
        /// </param>
        /// <param name="AggCellCount">
        /// number of cells from the grid <paramref name="ag"/> that should be aggregated in each aggregate cell
        /// </param>
        /// <returns>
        /// - 1st index: local cell mesh of the aggregate mesh
        /// - 2nd index: enumeration of cell parts
        /// </returns>
        private static int[][] AggregationKernel(IGridData ag, int AggCellCount) {
            int Jloc = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int D = ag.SpatialDimension;
            if (AggCellCount < 2)
                throw new ArgumentOutOfRangeException();





            // sort cells of parent grid by size:
            // we want to aggregate the smallest cells at first.
            int[] Perm = Jloc.ForLoop(j => j).OrderBy(j => ag.iLogicalCells.GetCellVolume(j)).ToArray();

            BitArray UsedCellMarker = new BitArray(Jloc);

            List<int[]> Coarsened_ComositeCells = new List<int[]>();

            // caching Bounding-Boxes and cell sizes (quite expensive to compute)
            BoundingBox[] Bbxes = new BoundingBox[Jloc];
            double[] CellVol = new double[Jloc];
            for(int j = 0; j < Jloc; j++) {
                Bbxes[j] = new BoundingBox(D);
                ag.iLogicalCells.GetCellBoundingBox(j, Bbxes[j]);
                CellVol[j] = ag.iLogicalCells.GetCellVolume(j);
            }


            //
            List<int> aggCell = new List<int>();
            List<int> NeighCandidates = new List<int>();

            // loop over aggregated cells of parent grid...
            //int[][] Neighbourship = ag.iLogicalCells.CellNeighbours;
            
           

            for (int i = 0; i < Jloc; i++) {
                int jCell = Perm[i]; // pick next cell

                int[] Neighbourship_jCell = ag.GetCellNeighboursViaEdges(jCell, OmmitPeriodic: true).Select(ttt => ttt.jCellLoc).ToArray();
                //int[] Neighbourship_jCell = Neighbourship[jCell];

                //Debug.Assert(Neighbourship_jCell.Contains(jCell) == false);
                if (!UsedCellMarker[jCell]) { // if the cell is not already agglomerated to another cell

                    aggCell.Clear();
                    aggCell.Add(jCell);
                    UsedCellMarker[jCell] = true;

                    for (int iPass = 1; iPass < AggCellCount; iPass++) {


                        // list of all neighbor cells which were not already aggregated to some other cell:
                        NeighCandidates.Clear();
                        foreach (int j in aggCell) {
                            Debug.Assert(j < Jloc);
                            Debug.Assert(UsedCellMarker[j] == true);
                            int[] Neighbourship_j = ag.GetCellNeighboursViaEdges(j, OmmitPeriodic: true).Select(ttt => ttt.jCellLoc).ToArray();
                            foreach (int jNeigh in Neighbourship_j) {
                            //foreach (int jNeigh in Neighbourship_jCell) {
                                if (jNeigh >= Jloc)
                                    continue;
                                if (UsedCellMarker[jNeigh] == true)
                                    continue;

                                Debug.Assert(aggCell.Contains(jNeigh) == false); // for all cells which are already in 'aggCell', the marker should be true
                                if (!NeighCandidates.Contains(jNeigh))
                                    NeighCandidates.Add(jNeigh);
                            }
                        }


                        int NN = NeighCandidates.Count;

                        if (NN > 0) {
                            double[] sizes = new double[NN];
                            BoundingBox[] aggBB = new BoundingBox[NN];
                            double[] aggBBaspect = new double[NN];
                            for (int iNeig = 0; iNeig < NN; iNeig++) { // loop over all candidates...
                                int jCellNeigh = NeighCandidates[iNeig];

                                aggBB[iNeig] = new BoundingBox(D); //ag.CompositeCellBB[jCell].CloneAs();
                                foreach (int jTaken in aggCell) {
                                    aggBB[iNeig].AddBB(Bbxes[jTaken]);
                                }

                                aggBB[iNeig].AddBB(Bbxes[jCellNeigh]);
                                //sizes[iNeig] = ag.iLogicalCells.GetCellVolume(jCell) + ag.iLogicalCells.GetCellVolume(jCellNeigh);
                                sizes[iNeig] = CellVol[jCell] + CellVol[jCellNeigh];
                                aggBBaspect[iNeig] = aggBB[iNeig].AspectRatio;
                            }

                            double[] RelSizes = sizes.CloneAs(); RelSizes.ScaleV(1.0 / sizes.Max());
                            double[] RelAspects = aggBBaspect.CloneAs(); RelAspects.ScaleV(1.0 / aggBBaspect.Max());

                            double[] Quality = new double[NN];
                            Quality.AccV(0.7, RelSizes);
                            Quality.AccV(0.3, RelAspects);

                            int iChoice = Quality.IndexOfMin(q => q);
                            int jNeigChoice = NeighCandidates[iChoice];

                            Debug.Assert(aggCell.Contains(jNeigChoice) == false);
                            aggCell.Add(jNeigChoice);
                            UsedCellMarker[jNeigChoice] = true;
                        } else {
                            // No neighbour available => unable to coarsen
                            break;
                        }
                    }

                    // add agglom cell
                    {
                        int[] aggCell_Fix = aggCell.ToArray();
#if DEBUG
                        foreach (int j in aggCell_Fix) {
                            Debug.Assert(j >= 0);
                            Debug.Assert(j < Jloc);
                            Debug.Assert(UsedCellMarker[j] == true);
                        }

#endif
                        Coarsened_ComositeCells.Add(aggCell_Fix);
                    }
                } else {
                    // cell already done.
                    continue;
                }
            }
            Debug.Assert(UsedCellMarker.ToBoolArray().Where(b => !b).Count() == 0, "some cell was not processed.");
            return Coarsened_ComositeCells.ToArray();
        }


        static bool IsAnc(IGridData gdat, AggregationGridData aGdat) {
            if (object.ReferenceEquals(aGdat.ParentGrid, gdat))
                return true;

            if (aGdat.ParentGrid is AggregationGridData)
                return IsAnc(gdat, (AggregationGridData)(aGdat.ParentGrid));

            return false;
        }


        /// <summary>
        /// assigns the cell index of the aggregated cell <em>j</em> to all (fine) grid cells that 
        /// the aggregated cell <em>j</em> consists of.
        /// </summary>
        static public void ColorDGField(this AggregationGridData ag, DGField f) {
            IGridData Anc = f.GridDat;
            if (!IsAnc(Anc, ag))
                throw new ArgumentException("Field 'f' must be defined on an ancestor grid of 'ag'.");


            f.Clear();
            int Jag = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int Janc = Anc.iLogicalCells.NoOfLocalUpdatedCells;

            Debug.Assert(Anc.iGeomCells.Count == ag.iGeomCells.Count);
           
            int[] jG2jL = Anc.iGeomCells.GeomCell2LogicalCell;

            int[] Colors = new int[Jag];
            BitArray Marked = new BitArray(Janc);
            for (int j = 0; j < Jag; j++) { // loop over logical/aggregate cells

                // determine colors of neighbor cells
                int[] Neighs = ag.iLogicalCells.CellNeighbours[j];
                var NeighColors = Neighs.Where(jN => jN < Jag).Select(jN => Colors[jN]);

                // select color for cell 'j'
                int iCol = 1;
                for(iCol = 1; iCol < 2 * Jag; iCol++) {
                    if(!NeighColors.Contains(iCol))
                        break;
                }
                Colors[j] = iCol;

                // color all logical cells in ancestor grid
                // idea: convert to geometrical and back to logical
                //    in this way we can e.g. skip multiple grid levels
                foreach (int jGeom in ag.iLogicalCells.AggregateCellToParts[j]) {
                    int jLogAnc;
                    if (jG2jL != null) {
                        jLogAnc = jG2jL[jGeom];
                    } else {
                        jLogAnc = jGeom;
                    }
                    if (!Marked[jLogAnc]) {
                        f.SetMeanValue(jLogAnc, iCol);
                        Marked[jLogAnc] = true;
                    } else {
                        // nop
                    }
                }
            }
        }
        
        

            /*

        /// <summary>
        /// Mapping from composite/aggegated cells to base grid cells.<br/>
        /// 1st index: composite/aggregated cell index <em>j</em> <br/>
        /// 2nd index: list of cells which from the aggregated cell <em>j</em>. <br/>
        /// </summary>
        public int[][] AggregateCells;


        int[] m_jCellToAggCell = null;
        int[] m_jCellToAggk = null;

        /// <summary>
        /// mapping from base grid cells to composite cells<br/>
        /// index: local cell index.
        /// </summary>
        public int[] jCellToAggCell {
            get {
                InitInverseMapping();
                return m_jCellToAggCell;
            }
        }

        /// <summary>
        /// mapping from base grid cells to cell index within aggregate cells<br/>
        /// index: local cell index.
        /// </summary>
        public int[] jCellToAggk {
            get {
                InitInverseMapping();
                return m_jCellToAggk;
            }
        }

        private void InitInverseMapping() {
            if(m_jCellToAggCell == null) {
                m_jCellToAggCell = new int[this.GridDat.Cells.NoOfLocalUpdatedCells];
                m_jCellToAggk = new int[this.GridDat.Cells.NoOfLocalUpdatedCells];
#if DEBUG
                m_jCellToAggCell.SetAll(-1234);
                m_jCellToAggk.SetAll(-123450);

#endif
                int JAGG = this.AggregateCells.Length;
                for(int jagg = 0; jagg < JAGG; jagg++) {
                    int[] AggCell = AggregateCells[jagg];

                    for(int k = AggCell.Length - 1; k >= 0; k--) {
                        int jCell = AggCell[k];
                        Debug.Assert(m_jCellToAggCell[jCell] < 0);
                        m_jCellToAggCell[jCell] = jagg;
                        m_jCellToAggk[jCell] = k;
                    }
                }
#if DEBUG
                foreach(int jagg in m_jCellToAggCell) {
                    Debug.Assert(jagg >= 0 && jagg < JAGG);
                }

                for(int j = 0; j < m_jCellToAggCell.Length; j++) {
                    int jAgg = m_jCellToAggCell[j];
                    int k = m_jCellToAggk[j];
                    Debug.Assert(this.AggregateCells[jAgg][k] == j);
                }
#endif
            }
        }


        /// <summary>
        /// Bounding boxes of the composite/aggregated cells.
        /// 
        /// Index: composite/aggregated cell index, correlates with first index of <see cref="AggregateCells"/>.
        /// </summary>
        public BoundingBox[] CompositeCellBB;


        */


    }
}
