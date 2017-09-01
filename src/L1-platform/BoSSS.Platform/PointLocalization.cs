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
using System.Text;
using System.Diagnostics;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP;

namespace BoSSS.Platform.Utils.Geom {
    
    /// <summary>
    /// This class constructs a binary geometrical tree (see <see cref="GeomBinTreeBranchCode"/>)
    /// from a cloud of points.
    /// </summary>
    public class PointLocalization {

        
        /// <summary>
        /// Easy-to-use constructor.
        /// </summary>
        /// <param name="_points">
        /// a cloud of points
        /// </param>
        /// <param name="BoundingBoxExtensionFactor">
        /// extends the bounding box of the points (<see cref="PointsBB"/>) by this factor (minimum is zero)
        /// </param>
        /// <param name="Permutation">
        /// output; permutation of <paramref name="_points"/> by sorting in tree; if null, an internal buffer is allocated
        /// </param>
        public PointLocalization(MultidimensionalArray _points, double BoundingBoxExtensionFactor = 0.01, int[] Permutation = null)
            : this(_points, FindBoundingBox(_points, BoundingBoxExtensionFactor), Permutation) {
        }

        static BoundingBox FindBoundingBox(MultidimensionalArray _points, double extFactor) {
            BoundingBox PointsBB = new BoundingBox(_points);
            PointsBB.ExtendByFactor(extFactor);
            return PointsBB;
        }

        /// <summary>
        /// Expert's constructor.
        /// </summary>
        /// <param name="_points"></param>
        /// <param name="BB">predefined bounding box, must contain all points in <paramref name="_points"/></param>
        /// <param name="Permutation">
        /// output; permutation of <paramref name="_points"/> by sorting in tree; if null, an internal buffer is allocated
        /// </param>
        public PointLocalization(MultidimensionalArray _points, BoundingBox BB, int[] Permutation) {
            using (new FuncTrace()) {
                int D = _points.GetLength(1);
                int N = _points.GetLength(0);
                
                PointsBB = BB;

                // find code
                // =========
                double[] pt = new double[D];
                GeomBinTreeBranchCode[] _Codes = new GeomBinTreeBranchCode[N];
                for (int n = 0; n < N; n++) {
                    _points.GetRow(n, pt);
                    _Codes[n] = GeomBinTreeBranchCode.CreateFormPoint(PointsBB, pt);
                }

                // sort
                // ====

                int[] IndexIntoQuadNodes;
                if (Permutation == null)
                    IndexIntoQuadNodes = new int[N];  // we only sort an index list into the arrays 'QuadNodes', 'quadweights' and 'Locations'
                else {
                    if (Permutation.Length != N)
                        throw new ArgumentException("false Länge", "Permutation");
                    IndexIntoQuadNodes = Permutation;
                }
                for (int i = 0; i < N; i++) IndexIntoQuadNodes[i] = i;
                Array.Sort<int>(IndexIntoQuadNodes, delegate(int a, int b) {
                    GeomBinTreeBranchCode _a = _Codes[a];
                    GeomBinTreeBranchCode _b = _Codes[b];
                    return _a.CompareTo(_b);
                }); // radix sort would be better here, because it would have linear runtime


                // store sorted
                // ============
                this.Codes = new GeomBinTreeBranchCode[N];
                this.Points = MultidimensionalArray.Create(N, D);
                for (int n = 0; n < N; n++) {
                    int it = IndexIntoQuadNodes[n];
                    for (int d = 0; d < D; d++)
                        this.Points[n, d] = _points[it, d];
                    this.Codes[n] = _Codes[it];
                }
            }
        }

        /// <summary>
        /// initializes <see cref="__tree"/> on demand
        /// </summary>
        private TreeNode[] InitTree() {
            if (__tree == null) {
                // construct tree
                // ==============

                TreeNodeTmp _tree = new TreeNodeTmp();
                //Console.WriteLine("entering tree ");
                //DateTime startRec = DateTime.Now;
                _tree.InitRecursive(default(GeomBinTreeBranchCode), 0x80000000, this, 0, Codes.Length, 0);
                //TimeSpan dt = DateTime.Now - startRec;
                //Console.WriteLine("tree finished: " + dt.ToString());
                __tree = new TreeNode[_tree.NoOfElements()];
                int dummy = 0;
                TreeNode.CollectInArray(_tree, __tree, ref dummy);

                // test
                // ====
                //vt = new VerifyTree(this);
                //vt.verify();
                //VerifyRec(0, this.Codes.Length, default(BoundingBoxCode),0);
            }
            return __tree;
        }

        TreeNode[] __tree;

        /// <summary>
        /// only for debug/testing purpose; verifies if the tree is build correctly
        /// </summary>
        /// <param name="i0">for recursion start, set to 0</param>
        /// <param name="Len">for recursion start, set to <see cref="Codes"/>.Length </param>
        /// <param name="cd">for recursion start, set to default value of <see cref="BoundingBoxCode"/></param>
        /// <param name="idx">for recursion start, set to 0</param>
        [Conditional("DEBUG")]
        public void VerifyRec(int i0, int Len, BoundingBoxCode cd, int idx) {
            TreeNode[] tree = InitTree();
            
            TreeNode nt = tree[idx];

            for (int n = 0; n < Codes.Length; n++) {
                bool innen = cd.IsInside(this.Codes[n]);
                bool innen2 = n >= i0 && n < (i0 + Len);
                if (innen != innen2)
                    throw new ApplicationException();
            }


            if (nt.Left >= 0) {
                int i0L = i0;
                int LenL = nt.iMid - i0L;
                VerifyRec(i0L,LenL,cd.GetSubBox(false),nt.Left);
            }
            if (nt.Right >= 0) {
                int i0R = nt.iMid;
                int LenR = i0 + Len - i0R;
                VerifyRec(i0R, LenR, cd.GetSubBox(true), nt.Right);
            }
        }



        /*
        // test code
        VerifyTree vt;
        class VerifyTree {

            int Level;

            //bool[] touchedPointsPerLevel;

            public VerifyTree(PointLocalization _owner) {
                Level = 0;
                bb = (BoundingBox)_owner.PointsBB.Clone();
                owner = _owner;
                code.Code = 0;

                i0 = 0;
                L = owner.Codes.Length;

                //touchedPoints = new bool[owner.Codes.Length];
                

                left = new VerifyTree(false,this);
                right = new VerifyTree(true,this);
            }

            private VerifyTree(bool leftorRight, VerifyTree parrent) {
                owner = parrent.owner;

                this.Level = parrent.Level + 1;
                this.code = parrent.code;

                GeomBinTreeBranchCode add;
                add.Code = (0x80000000 >> parrent.Level);

                if (leftorRight) {
                    // right
                    this.code.Code += add.Code;
                    i0 = Array.BinarySearch<GeomBinTreeBranchCode>(owner.Codes, 0, owner.Codes.Length, this.code);
                    if (i0 < 0) i0 = ~i0;

                    add.Code -= 1;
                    add.Code += this.code.Code;
                    int iE = Array.BinarySearch<GeomBinTreeBranchCode>(owner.Codes, 0, owner.Codes.Length, add);
                    if (iE < 0) iE = ~iE;
                    L = iE - i0;
                } else {
                    // left 
                    i0 = Array.BinarySearch<GeomBinTreeBranchCode>(owner.Codes, 0, owner.Codes.Length, this.code);
                    if (i0 < 0) i0 = ~i0;

                    add.Code += this.code.Code;
                    int iE = Array.BinarySearch<GeomBinTreeBranchCode>(owner.Codes, 0, owner.Codes.Length, add);
                    if (iE < 0) iE = ~iE;
                    L = iE - i0;
                }

                this.bb = new BoundingBox(owner.Points.GetLength(1));
                owner.PointsBB.SubBoxFromCode(this.bb, this.code, (uint)this.Level);

                if (Level < 8) {
                    left = new VerifyTree(false, this);
                    right = new VerifyTree(true, this);
                }
            }

            public VerifyTree Get(GeomBinTreeBranchCode cd, int level) {
                if (this.code == cd && this.Level == level)
                    return this;

                VerifyTree ret = null;
                if (left != null)
                    ret = left.Get(cd, level);
                if (ret != null)
                    return ret;

                if (right != null)
                    ret = right.Get(cd, level);
                return ret;
            }


            public void verify() {
                int N = owner.Codes.Length;

                bool[] geomCheck = new bool[N];
                bool[] bitmaskCheck = new bool[N];
                bool[] indexCheck = new bool[N];

                uint mask = 0;
                for (int i = 0; i < Level; i++) {
                    mask += 0x80000000 >> i;
                }

                int failedTests = 0;

                if (this.code.Code == 0x80000000 && Level == 2)
                    Console.Write("");

                for (int n = 0; n < N; n++) {
                    double[] _pt = new double[owner.Points.GetLength(1)];
                    for (int d = 0; d < _pt.Length; d++) _pt[d] = owner.Points[n, d];

                    geomCheck[n] = this.bb.IsInsideStrict(_pt);

                    if (n >= i0 && n < (i0 + L))
                        indexCheck[n] = true;

                    uint code = owner.Codes[n].Code & mask;
                    if (code == this.code.Code)
                        bitmaskCheck[n] = true;

                    if (!(geomCheck[n] == indexCheck[n] && indexCheck[n] == bitmaskCheck[n]))
                        failedTests++;


                    //if (n >= i0 && n < (i0 + L)) {
                    //    // point should be inside, i.e. loc must be true
                    //    if (!loc)
                    //        throw new ApplicationException("error in alg");
                    //} else {
                    //    // loc should be false
                    //    if (loc)
                    //        throw new ApplicationException("error in alg 2");
                    //}                    
                }

                if(failedTests != 0)
                    Console.WriteLine("Fialed Tests: " + failedTests);

                if (left != null) left.verify();
                if (right != null) right.verify();

                if (left != null && right != null) {
                    if (left.L + right.L != this.L)
                        throw new ApplicationException();

                    if (left.i0 != this.i0)
                        throw new ApplicationException();

                    if (left.i0 + left.L != right.i0)
                        throw new ApplicationException();
                }
            }

            VerifyTree left;

            VerifyTree right;

            PointLocalization owner;

            GeomBinTreeBranchCode code;

            BoundingBox bb;

            int i0;

            public int I0 {
                get { return i0; }
            }

            int L;

            public int Len {
                get { return L; }
            }
        }
        // end test code */

        /// <summary>
        /// points, sorted according to their codes;
        /// <br/>
        /// 1st index: point index; <br/>
        /// 2nd index: spatial coordinate;
        /// </summary>
        public MultidimensionalArray Points;

        /// <summary>
        /// for each point in <see cref="Points"/>, its branch code.
        /// </summary>
        public GeomBinTreeBranchCode[] Codes;

        /// <summary>
        /// bounding box of all points, possibly extended a bit
        /// </summary>
        public BoundingBox PointsBB;
        
        /// <summary>
        /// Finds all points (in <see cref="Points"/>) within
        /// the radius <paramref name="eps"/> around point <paramref name="pt"/>;
        /// </summary>
        /// <param name="IndicesIntoPoints">
        /// Output; on exit, the indices of all found points (1st index of <see cref="Points"/>);
        /// </param>
        /// <param name="eps"></param>
        /// <param name="pt"></param>
        public void FindNearPoints(ICollection<int> IndicesIntoPoints, double eps, params double[] pt) {
            if (helper == null)
                helper = new BoundingBox(this.PointsBB.D);
            InitTree();

            IndicesIntoPoints.Clear();
//            GeomBinTreeBranchCode cd;
//            cd.Code = 0;
            BoundingBoxCode cd; cd.Branch.Code = 0; cd.SignificantBits = 0;
            uint MaxDepth = 0;
            FindNearPointsRecursice(IndicesIntoPoints, eps * eps, pt, 0, Codes.Length, this.PointsBB.D, cd, 0x80000000, ref MaxDepth, false, 0);
            //Console.WriteLine(MaxDepth);

            
        }

        BoundingBox helper;

        void FindNearPointsRecursice(ICollection<int> IndicesIntoPoints, double epsPow2, double[] pt, int i0, int Len, int _D, BoundingBoxCode BoxCode, uint shifti, ref uint MaxDepth, bool forceTst, int treeNodeIdx) {
            TreeNode[] tree = __tree;
            uint Depth = BoxCode.SignificantBits;
            {
                /*/ test and debug:
                VerifyTree l = vt.Get(BoxCode, (int)Depth);
                if (l != null) {
                    if (l.I0 != i0)
                        throw new ApplicationException();
                    if (l.Len != Len)
                        throw new ApplicationException();
                }

                PointsBB.SubBoxFromCode(helper, BoxCode, Depth);
                for (int n = 0; n < Codes.Length; n++) {
                    double[] _pt = new double[Points.GetLength(1)];
                    for (int d = 0; d < pt.Length; d++) _pt[d] = Points[n, d];

                    bool loc = helper.IsInsideStrict(_pt);

                    if (n >= i0 && n < (i0 + Len)) {
                        if (!loc)
                            throw new ApplicationException("error in alg");
                    } else {
                        if (loc)
                            throw new ApplicationException("error in alg 2");
                    }
                }
                // end */


                if (Len <= 0)
                    // no points in branch => terminate recursion
                    return;
                                
                PointsBB.SubBoxFromCode(helper, BoxCode);
                //helper.ExtendByFactor(0.01);
                double di = helper.Distance(pt);
                if (di * di > epsPow2) {
                    //if (helper.IsInside(pt))
                    //    throw new ApplicationException("ohne Worte");


                    //for (int n = i0; n < (i0 + Len); n++) {
                    //    double distPow2 = 0;
                    //    for (int d = _D - 1; d >= 0; d--) {
                    //        double delta = pt[d] - this.Points[n, d];
                    //        distPow2 += delta * delta;
                    //    }

                    //    if (distPow2 <= epsPow2)
                    //        throw new ApplicationException("ohne Worte");
                    //}
                    
                    
                    return; // no interesting points in box => terminate recursion
                }
            }

            MaxDepth = Math.Max(MaxDepth, Depth);
            
            if (Len <= 2 || Depth == 32 || forceTst) {
                // only a few points remaining, or maximum tree depth reached - test all of them

                for (int n = i0; n < (i0 + Len); n++) {
                    double distPow2 = 0;
                    for (int d = _D - 1; d >= 0; d--) {
                        double delta = pt[d] - this.Points[n, d];
                        distPow2 += delta * delta;
                    }

                    if (distPow2 <= epsPow2)
                        IndicesIntoPoints.Add(n);
                }

            } else {
                // do recursion;

                BoundingBoxCode BoxB = BoxCode; BoxB.Branch.Code += shifti;
                shifti = shifti >> 1;
                Depth++;
                BoxB.SignificantBits++;
                BoxCode.SignificantBits++;
                                
                int idxMid = tree[treeNodeIdx].iMid;


                int LenLeft = idxMid - i0;
                int LenRigt = Len + i0 - idxMid;

                bool forceLeft = (LenLeft == Len);
                bool forceRigt = (LenRigt == Len);

                FindNearPointsRecursice(IndicesIntoPoints, epsPow2, pt, i0,     LenLeft, _D, BoxCode, shifti, ref MaxDepth, forceLeft, tree[treeNodeIdx].Left );
                FindNearPointsRecursice(IndicesIntoPoints, epsPow2, pt, idxMid, LenRigt, _D, BoxB,    shifti, ref MaxDepth, forceRigt, tree[treeNodeIdx].Right);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="branchCode">the branch code; only the bits 32 to 32-<paramref name="BranchDepth"/> are taken into account</param>
        /// <param name="BranchDepth">num</param>
        /// <param name="i0"></param>
        /// <param name="Len"></param>
        public void GetPointsInBranch(GeomBinTreeBranchCode branchCode, int BranchDepth, out int i0, out int Len) {
            if (BranchDepth < 0 || BranchDepth > 32)
                throw new ArgumentException("must be between 0 (including) and 32 (including).", "BranchDepth");

            {
                // ancient implementation
                //i0 = 0;
                //Len = Codes.Length;
                //GetPointsInBranchRecursive(branchCode, BranchDepth, ref i0, ref Len, 0x80000000, 0);
            }

            {
                // faster ????

                uint mask = 0;
                uint m = 0x80000000;
                for (int i = 0; i < BranchDepth; i++) {
                    mask += m >> i;
                }
                
                GeomBinTreeBranchCode loCode;
                loCode.Code = branchCode.Code & mask;

                GeomBinTreeBranchCode hiCode;
                hiCode.Code = branchCode.Code | (~mask);

                //int _i0, _Len;
                if (hiCode < this.Codes[0] || loCode > this.Codes[this.Codes.Length - 1]) {
                    i0 = 0;
                    Len = 0;
                } else {

                    int loi = Array.BinarySearch(this.Codes, loCode); // use the .NET binary search ...
                    if (loi < 0)
                        loi = ~loi;

                    int hii = Array.BinarySearch(this.Codes, hiCode); // use the .NET binary search ...
                    if (hii < 0) {
                        hii = ~hii;
                        hii--;
                    }

                    i0 = loi;
                    Len = hii - loi + 1;
                }


            }

            Test_GetPointsInBranch(i0, Len, branchCode, BranchDepth);


        }
        
        /// <summary>
        /// old implementation of <see cref="GetPointsInBranch"/>
        /// </summary>
        void GetPointsInBranchRecursive(GeomBinTreeBranchCode branchCode, int BranchDepth, ref int i0, ref int Len, uint mask, int idxIntoTree) {
            if (BranchDepth <= 0)
                return;
            BranchDepth--;

            TreeNode[] tree = InitTree();

            int iMid = tree[idxIntoTree].iMid;
            //if (iMid < 0) {
            //    iMid = Array.BinarySearch<GeomBinGeomTreeCode>(this.Codes, 0, this.Codes.Length, ???);
            //    if (iMid < 0) iMid = ~iMid;
            //}

            if ((branchCode.Code & mask) == 0) {
                // left
                Len = iMid - i0;
                if( Len > 0)
                    GetPointsInBranchRecursive(branchCode, BranchDepth, ref i0, ref Len, mask >> 1, tree[idxIntoTree].Left);
            } else {
                // right
                Len = i0 + Len - iMid;
                i0 = iMid;
                if (Len > 0)
                    GetPointsInBranchRecursive(branchCode, BranchDepth, ref i0, ref Len, mask >> 1, tree[idxIntoTree].Right);
            }
        }


        [Conditional("DEBUG")]
        void Test_GetPointsInBranch(int i0, int Len, GeomBinTreeBranchCode branchCode, int BranchDepth) {
            // logischer test 
            {
                /*
                uint mask = 0;
                uint m = 0x80000000;
                for (int i = 0; i < BranchDepth; i++) {
                    mask += m >> i;
                }


                GeomBinTreeBranchCode loCode;
                loCode.Code = branchCode.Code & mask;

                GeomBinTreeBranchCode hiCode;
                hiCode.Code = branchCode.Code | (~mask);

                
                int L = Codes.Length;
                for (int l = 0; l < L; l++) {
                    var Code = Codes[l].Code;

                    bool inside = (l >= i0 && l < (i0 + Len));
                    bool gt = (Code >= loCode.Code) && (Code <= hiCode.Code);

                    if (inside != gt)
                        Console.WriteLine("schas");
                }
                

                /*
                int InsideErr = 0, OutsideErr = 0;

                int L = Codes.Length;
                for (int l = 0; l < L; l++) {
                    var Code = Codes[l].Code;

                    if ((Code & mask) == (branchCode.Code & mask)) {
                        // inside
                        if (l < i0 || l >= (i0 + Len))
                            InsideErr++;
                    } else {
                        if (l >= i0 && l < (i0 + Len))
                            OutsideErr++;
                    }
                }

                if (InsideErr != 0 || OutsideErr != 0)
                    throw new ApplicationException("internal error - should not happen.");
                 */
            }

            // geometrischer test
            {
                BoundingBox bb = new BoundingBox(this.Points.GetLength(1));
                BoundingBoxCode bbcode;
                bbcode.SignificantBits = (uint)BranchDepth;
                bbcode.Branch = branchCode;
                this.PointsBB.SubBoxFromCode(bb, bbcode);
                int D = this.Points.GetLength(1);

                double[] pt = new double[D];
                for (int i = i0; i < (i0 + Len); i++) {
                    for (int d = 0; d < D; d++) {
                        pt[d] = this.Points[i, d];
                    }

                    if (!bb.Contains(pt))
                        throw new ApplicationException("test failed.");
                }
            }
        }


        /// <summary>
        /// 
        /// </summary>
        struct TreeNode {
            
            /// <summary>
            /// index into <see cref="Points"/> (1st index) and <see cref="Codes"/> which separates the
            /// 'left' and 'right' tree branch
            /// </summary>
            public int iMid;
            
            /// <summary>
            /// index of 1st child in array (bit of corresponding level is 0)
            /// </summary>
            public int Left;
            
            /// <summary>
            /// index of 2nd child in array (bit of corresponding level is 1)
            /// </summary>
            public int Right;

            /// <summary>
            /// converts a tree, represented as <see cref="TreeNodeTmp"/> (heap objects)
            /// into a compact array of value types
            /// </summary>
            /// <param name="treeTmp"></param>
            /// <param name="targ"></param>
            /// <param name="cnt"></param>
            public static void CollectInArray(TreeNodeTmp treeTmp, TreeNode[] targ, ref int cnt) {
                int ithis = cnt;
                targ[ithis].iMid = treeTmp.iMid;
                if (treeTmp.Left != null) {
                    cnt++;
                    targ[ithis].Left = cnt;
                    CollectInArray(treeTmp.Left, targ, ref cnt);
                } else {
                    targ[ithis].Left = int.MinValue;
                }
                if (treeTmp.Right != null) {
                    cnt++;
                    targ[ithis].Right = cnt;
                    CollectInArray(treeTmp.Right, targ, ref cnt);
                } else {
                    targ[ithis].Right = int.MinValue;
                }
            }

        }
        
        class TreeNodeTmp {
            public int iMid;
            public TreeNodeTmp Left;
            public TreeNodeTmp Right;

            /// <summary>
            /// computes the total number of elements in the tree
            /// </summary>
            public int NoOfElements() {
                int ret = 1;
                if (Left != null) ret += Left.NoOfElements();
                if (Right != null) ret += Right.NoOfElements();
                return ret;
            }

            [Conditional("DEBUG")]
            public void Test(int i0, int iMid, PointLocalization owner, BoundingBoxCode b) {
                for (int i = 0; i < owner.Codes.Length; i++) {
                    bool innen = b.IsInside(owner.Codes[i]);
                    bool innen2 = i >= i0 && i < iMid;
                    if (innen != innen2)
                        throw new ApplicationException();
                }
            }

            /*
            static int MyBinarySearch(GeomBinTreeBranchCode[] array, int i0, int Length, GeomBinTreeBranchCode searchValue) {
                if(Length <= 0)
                    return -1;

                if (array[i0] >= searchValue)
                    return i0;
                if (array[i0 + Length - 1] < searchValue)
                    return i0 + Length;

                if (Length >= 2) {
                    int ifnd = int.MinValue;
                    MyBinarySearchRec(array, i0, Length, i0, Length, ref searchValue, out ifnd);
                    return ifnd;
                } else {
                    if (array[i0] >= searchValue)
                        return i0;
                    else
                        return i0 + 1;
                }
            }

            static void MyBinarySearchRec(GeomBinTreeBranchCode[] array, int i0, int Length, int _i0, int _Len, ref GeomBinTreeBranchCode searchValue, out int ifnd) {
                if( Length < 2)
                    throw new ApplicationException();
                int iMid = i0 + Length / 2;

                if (array[iMid] >= searchValue && array[iMid - 1] < searchValue) {
                    ifnd = iMid;
                } else {
                    if (array[iMid] >= searchValue) {
                        MyBinarySearchRec(array, i0, Length, iMid, i0 + _Len - iMid, ref searchValue, out ifnd);
                    } else {
                        MyBinarySearchRec(array, i0, Length, i0, iMid - i0, ref searchValue, out ifnd);
                    }

                }

                //if( iMid == i0 && array[iMid] >= searchValue && ) {
                //    return i0;
                //}
                //if( 
            }*/

            public void InitRecursive(GeomBinTreeBranchCode BoxCode, uint shifti, PointLocalization owner, int i0, int Len, uint RecDepth) {
                GeomBinTreeBranchCode BoxB = BoxCode; BoxB.Code |= shifti;


                int __iMid;
                for (__iMid = i0; __iMid < (i0 + Len); __iMid++) {
                    if (owner.Codes[__iMid].Code >= BoxB.Code)
                        break;
                    // if there is a performance problem, this linear search should be exchanged by a binary one
                }
                this.iMid = __iMid;

                //{
                //    GeomBinTreeBranchCode search = BoxB;
                //    search.Code -= 1;
                //    int iMid2 = MyBinarySearch(owner.Codes, 0, owner.Codes.Length, search);
                //    if (iMid2 != iMid)
                //        throw new ApplicationException();
                //}


                
                //if (iMid < 0) iMid = ~iMid;

                int LenLeft = iMid - i0;
                int LenRigt = Len + i0 - iMid;

                {
                    /*
                    // test left
                    {
                        BoundingBoxCode cdLeft;
                        cdLeft.Branch = BoxCode;
                        cdLeft.SignificantBits = RecDepth + 1;
                        Test(i0, iMid, owner, cdLeft);
                    }

                    // test right
                    {
                        BoundingBoxCode cdRight; 
                        cdRight.Branch = BoxB;
                        cdRight.SignificantBits = RecDepth + 1;
                        Test(iMid, i0 + Len, owner, cdRight);
                    }
                    */
                }

                if (shifti == 1)
                    return; // reached max. supported tree depth
                shifti = shifti >> 1;


                if (LenLeft > 0) {
                    Left = new TreeNodeTmp();
                    Left.InitRecursive(BoxCode, shifti, owner, i0, LenLeft, RecDepth + 1);
                }
                if (LenRigt > 0) {
                    Right = new TreeNodeTmp();
                    Right.InitRecursive(BoxB, shifti, owner, iMid, LenRigt, RecDepth + 1);
                }
            }
        }
    }
}
