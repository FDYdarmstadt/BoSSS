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
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Platform.Utils.Geom {

    /// <summary>
    /// primitive bounding box
    /// </summary>
    public class BoundingBox : ICloneable {

        /// <summary>
        /// creates an empty bounding box
        /// </summary>
        /// <param name="D">dimension of the bounding box</param>
        public BoundingBox(int D) {
            if (D < 1)
                throw new ApplicationException("zero or lower dimensional bounding boxes are not known.");
            if (D > 3)
                throw new ArgumentException("more than 3D is not supported");

            Min = new double[D];
            Max = new double[D];
            Clear();
        }

        /// <summary>
        /// the aspect ratio of the bounding box.
        /// </summary>
        public double AspectRatio {
            get {
                double MinDist = double.MaxValue;
                double Maxdist = -double.MaxValue;

                for (int d = this.D - 1; d >= 0; d--) {
                    double Dist = this.Max[d] - this.Min[d];
                    Debug.Assert(Dist >= 0);
                    Maxdist = Math.Max(Maxdist, Dist);
                    MinDist = Math.Min(MinDist, Dist);
                }

                return Maxdist / MinDist;
            }
        }

        /// <summary>
        /// Distance between <see cref="Min"/> and <see cref="Max"/>.
        /// </summary>
        public double Diameter {
            get {
                return GenericBlas.L2Dist(this.Min, this.Max);
            }
        }

        /// <summary>
        /// creates a box that contains a given cloud of points;
        /// </summary>
        /// <param name="_points">
        /// cloud of points
        ///  - 1st index: point index;
        ///  - 2nd index: spatial coordinate;
        /// </param>
        public BoundingBox(params double[][] _points)
            : this(_points[0].GetLength(0)) //
        {
            foreach (var pete in _points)
                this.AddPoint(pete);
        }


        /// <summary>
        /// Creates a box that contains a given cloud of points;
        /// </summary>
        /// <param name="_points">
        /// cloud of points;
        ///  - 1st index: point index;
        ///  - 2nd index: spatial coordinate;
        /// </param>
        public BoundingBox(double[,] _points) : this(_points.GetLength(1)) {
            this.AddPoints(_points);
        }

        /// <summary>
        /// Creates a box that contains a given cloud of points;
        /// </summary>
        /// <param name="_points">
        /// cloud of points;
        /// - 1st index: point index;
        /// - 2nd index: spatial coordinate;
        /// </param>
        public BoundingBox(MultidimensionalArray _points) : this(_points.GetLength(1)) {
            this.AddPoints(_points);
        }

        /// <summary>
        /// extends the bounding box by some factor into each direction
        /// </summary>
        /// <param name="extFactor">
        /// must be larger or equal to 0; a factor of 0 keeps the bounding box unchanged.
        /// </param>
        public void ExtendByFactor(double extFactor) {

            if (extFactor < 0)
                throw new ArgumentException("bounding box extension factor must be larger or equal to 0.0.");
            for (int d = 0; d < D; d++) {
                double diam_d = this.Max[d] - this.Min[d];

                diam_d *= extFactor;
                this.Min[d] -= diam_d;
                this.Max[d] += diam_d;
            }
        }

        /// <summary>
        /// the lower bound of the bounding box
        /// </summary>
        public double[] Min;

        /// <summary>
        /// the upper bound of the bounding box
        /// </summary>
        public double[] Max;


        /// <summary>
        /// Minimum witdh across all spatial dimensions.
        /// </summary>
        public double h_min {
            get {
                double _min = double.MaxValue;
                int D = this.D;
                Debug.Assert(Min.Length == Max.Length);
                for (int d = 0; d < D; d++) {
                    _min = Math.Min(Max[d] - Min[d], _min);
                }
                return _min;
            }
        }

        /// <summary>
        /// Maximum witdh across all spatial dimensions.
        /// </summary>
        public double h_max {
            get {
                double _max = 0;
                int D = this.D;
                Debug.Assert(Min.Length == Max.Length);
                for (int d = 0; d < D; d++) {
                    _max = Math.Max(Max[d] - Min[d], _max);
                }
                return _max;
            }
        }


        /// <summary>
        /// Empties this box.
        /// </summary>
        public void Clear() {
            for (int d = D - 1; d >= 0; d--) {
                this.Min[d] = double.MaxValue;
                this.Max[d] = -double.MaxValue;
            }
        }

        /// <summary>
        /// Adds a list of points to this bounding box.
        /// </summary>
        public void AddPoints(double[,] _points) {
            int N = _points.GetLength(0);
            int _D = _points.GetLength(1);
            if (this.D != _D)
                throw new ArgumentException("wrong spatial dimension of point array.");

            for (int d = 0; d < D; d++) {
                double mind = this.Min[d];
                double maxd = this.Max[d];

                for (int n = 0; n < N; n++) {
                    mind = Math.Min(mind, _points[n, d]);
                    maxd = Math.Max(maxd, _points[n, d]);
                }

                Min[d] = mind;
                Max[d] = maxd;
            }
        }

        /// <summary>
        /// adds a list of points to this bounding box
        /// </summary>
        /// <param name="_points"></param>
        public void AddPoints(MultidimensionalArray _points) {
            int N = _points.GetLength(0);
            int _D = _points.GetLength(1);
            if (this.D != _D)
                throw new ArgumentException("wrong spatial dimension of point array.");

            for (int d = 0; d < D; d++) {
                double mind = this.Min[d];
                double maxd = this.Max[d];

                for (int n = 0; n < N; n++) {
                    mind = Math.Min(mind, _points[n, d]);
                    maxd = Math.Max(maxd, _points[n, d]);
                }

                Min[d] = mind;
                Max[d] = maxd;
            }
        }

        /// <summary>
        /// adds a point to the bounding box, i.e. enlarges it when the point is outside and does not change it
        /// when the point is inside.
        /// </summary>
        public void AddPoint(params double[] pt) {
            int D = pt.Length;
            if (D != Min.Length)
                throw new ArgumentException("wrong spatial dimension of point");
            for (int d = D - 1; d >= 0; d--) {
                Min[d] = Math.Min(Min[d], pt[d]);
                Max[d] = Math.Max(Max[d], pt[d]);
            }
        }

        /// <summary>
        /// adds a point to the bounding box, i.e. enlarges it when the point is outside and does not change it
        /// when the point is inside (3D - version);
        /// </summary>
        /// <param name="pt"></param>
        public void AddPoint(Vector3D pt) {
            if (D != 3)
                throw new ArgumentException("wrong spatial dimension of point");
            Min[0] = Math.Min(Min[0], pt.x);
            Max[0] = Math.Max(Max[0], pt.x);
            Min[1] = Math.Min(Min[1], pt.y);
            Max[1] = Math.Max(Max[1], pt.y);
            Min[2] = Math.Min(Min[2], pt.z);
            Max[2] = Math.Max(Max[2], pt.z);
        }

        /// <summary>
        /// adds a point to the bounding box, i.e. enlarges it when the point is outside and does not change it
        /// when the point is inside (2D - version);
        /// </summary>
        /// <param name="pt"></param>
        public void AddPoint(Vector2D pt) {
            if (D != 2)
                throw new ArgumentException("wrong spatial dimension of point");
            Min[0] = Math.Min(Min[0], pt.x);
            Max[0] = Math.Max(Max[0], pt.x);
            Min[1] = Math.Min(Min[1], pt.y);
            Max[1] = Math.Max(Max[1], pt.y);
        }

        /// <summary>
        /// increases the size of this bounding box to contain the box <paramref name="other"/>
        /// </summary>
        public void AddBB(BoundingBox other) {
            if (other.D != this.D)
                throw new ArgumentException("mismatch in spatial dimension.");

            for (int d = D - 1; d >= 0; d--) {
                this.Min[d] = Math.Min(this.Min[d], other.Min[d]);
                this.Max[d] = Math.Max(this.Max[d], other.Max[d]);
            }
        }

        /// <summary>
        /// true if the volume of the bounding box is zero (in any measure of dimension higher or equal to 1)
        /// </summary>
        public bool IsEmpty {
            get {
                for (int d = Max.Length - 1; d >= 0; d--)
                    if (Max[d] < Min[d])
                        return true;
                return false;
            }
        }

        /// <summary>
        /// spatial dimension (length of <see cref="Min"/>, <see cref="Max"/>);
        /// </summary>
        public int D {
            get {
                if (Min.Length != Max.Length)
                    throw new ApplicationException("illegal modification doe to bounding box.");
                return Min.Length;
            }
        }

        /// <summary>
        /// tests whether a point is inside this bounding box or not
        /// </summary>
        public bool Contains(params double[] pt) {
            int D = Min.Length;
            for (int d = 0; d < D; d++) {
                if (pt[d] < Min[d] || pt[d] > Max[d])
                    return false;
            }
            return true;
        }

        /// <summary>
        /// tests whether a point is inside this bounding box or not
        /// </summary>
        public bool ContainsStrict(params double[] pt) {
            int D = Min.Length;
            for (int d = 0; d < D; d++) {
                if (pt[d] < Min[d] || pt[d] >= Max[d])
                    return false;
            }
            return true;
        }

        /// <summary>
        /// true, of the bounding box <paramref name="other"/> is fully contained in this box
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Contains(BoundingBox other) {
            int _D = D;

            for (int d = 0; d < _D; d++) {
                if (other.Min[d] < this.Min[d])
                    return false;
                if (other.Max[d] > this.Max[d])
                    return false;
            }

            return true;
        }


        #region ICloneable Members

        /// <summary>
        /// creates a non-shallow copy of this bounding box
        /// </summary>
        public object Clone() {
            BoundingBox ret = new BoundingBox(this.D);
            ret.Max = (double[])this.Max.Clone();
            ret.Min = (double[])this.Min.Clone();
            return ret;
        }

        #endregion


        /// <summary>
        /// computes the leave bounding-box of a geometric binary tree
        /// </summary>
        /// <param name="subBB">
        /// on exit, the bounding box associated with the branch code <paramref name="code"/>;
        /// The output is not returned (via return value) to avoid heap allocation.
        /// </param>
        /// <param name="code">
        /// branch code 
        /// </param>
        public void LeaveBoxFromCode(BoundingBox subBB, GeomBinTreeBranchCode code) {
            BoundingBoxCode cd;
            cd.Branch = code;
            cd.SignificantBits = 32;
            SubBoxFromCode(subBB, cd);
        }


        /// <summary>
        /// computes the leave bounding-box of a geometric binary tree
        /// </summary>
        /// <param name="subBB">
        /// on exit, the bounding box associated with the branch code <paramref name="bbCode"/>;
        /// The output is not returned (via return value) to avoid heap allocation.
        /// </param>
        /// <param name="bbCode">
        /// bounding box code
        /// </param>
        public void SubBoxFromCode(BoundingBox subBB, BoundingBoxCode bbCode) {
            if (subBB.D != this.D)
                throw new ArgumentException("must have the same dimension as this bounding box.", "subBB");
            int _D = D;
            uint TreeDepth = bbCode.SignificantBits;
            if (TreeDepth > 32)
                throw new ArgumentException("tree Depth is limited to 32.");

            Array.Copy(this.Max, subBB.Max, _D);
            Array.Copy(this.Min, subBB.Min, _D);

            uint cd = bbCode.Branch.Code;
            uint mxx = 0x80000000;
            //for (int i = 0; i < 30; i++)
            //    mxx *= 2;

            int cur_d = 0;
            for (int i = (int)TreeDepth - 1; i >= 0; i--) {
                double mean = (subBB.Min[cur_d] + subBB.Max[cur_d]) * 0.5;

                // split bb;
                //uint kaes = cd / mxx;
                //switch (kaes) {
                //    case 0: subBB.Max[cur_d] = mean; break;
                //    case 1: subBB.Min[cur_d] = mean; break;
                //    default: throw new ApplicationException("error in algorithm");
                //}
                if ((cd & mxx) != 0) {
                    subBB.Min[cur_d] = mean;
                } else {
                    subBB.Max[cur_d] = mean;
                }

                // increment
                mxx = mxx >> 1;
                //mxx = Math.Max(1, mxx / 2);
                //cd = cd % mxx;
                cur_d++;
                if (cur_d >= _D)
                    cur_d = 0;
            }
        }

        /// <summary>
        /// True, if this bounding box overlaps with another one.
        /// </summary>
        public bool Overlap(BoundingBox other) {
            if (this.D != other.D)
                throw new ArgumentException("Cannot compare boxes of different spatial dimension.");

            int _d = this.D;
            for (int d = 0; d < _d; d++) {
                if (other.Min[d] > this.Max[d])
                    return false;
                if (other.Max[d] < this.Min[d])
                    return false;
            }

            return true;
        }


        /// <summary>
        /// finds the minimum distance of point <paramref name="pt"/> to one point within
        /// the box;
        /// </summary>
        /// <param name="pt"></param>
        /// <returns></returns>
        public double Distance(params double[] pt) {

            int _D = pt.Length;
            if (_D != D)
                throw new ArgumentException("wrong spatial dimension of point");
            if (_D <= 0 || _D > 3)
                throw new NotSupportedException("Only supported for spatial dimensions 1, 2 and 3.");

            // 'special' cases
            // ===============

            double dist = 0;
            int MatchCnt = 0;

            for (int d = 0; d < _D; d++) {
                if (pt[d] >= this.Min[d] && pt[d] <= this.Max[d])
                    MatchCnt++;
                //dist += 0;
                else {
                    dist += Math.Min(
                        Math.Abs(pt[d] - this.Min[d]),
                        Math.Abs(pt[d] - this.Max[d])).Pow2();
                }
            }
            dist = Math.Sqrt(dist);

            if (MatchCnt == _D)
                return 0.0; // point 'pt' is inside the box
            if (MatchCnt == (_D - 1))
                return dist; // there is a line through point 'pt', which is perpendicular to one face of the box

            if (MatchCnt == 1 && _D == 3) {
                return dist;
            }

            // 'regular' cases: compare with all edge points
            // =============================================
            Vector3D vec; // we also embedd the 1D and the 2D case in a 3D Vector (Vector3D is a stack object, while an array would be a heap object)
            vec.x = 0;
            vec.y = 0;
            vec.z = 0;

            Vector3D _pt;
            _pt.x = 0;
            _pt.y = 0;
            _pt.z = 0;
            for (int d = 0; d < _D; d++)
                _pt[d] = pt[d];

            int NoOfVertice = 1;
            for (int d = _D - 1; d >= 0; d--)
                NoOfVertice *= 2; // number of vertices of the box

            dist = double.MaxValue;
            for (int c = 0; c < NoOfVertice; c++) {

                // construct vertex ...
                int div = NoOfVertice / 2;
                int cc = c;
                for (int d = 0; d < _D; d++) {
                    int kase = cc / div;
                    switch (kase) {
                        case 0:
                            vec[d] = this.Min[d];
                            break;
                        case 1:
                            vec[d] = this.Max[d];
                            break;
                        default:
                            throw new ApplicationException("error in alg");
                    }

                    cc = cc % div;
                    div /= 2;
                }

                // compare distance
                dist = Math.Min(dist, Vector3D.Dist(_pt, vec));
            }
            return dist;
        }
    }

    /// <summary>
    /// This structure represents the so-called _branch code_ of a binary geometric tree:
    /// 
    /// We construct a binary geometric tree by cutting a bounding box symmetrically along the x-Axis,
    /// then along the y-Axis, (then the z-Axis, if present), and again by x-Axis, y-Axis, ... <br/>
    /// The first cut corresponds with the most-significant bit (2^31, in integer representation) of <see cref="Code"/>, the second 
    /// cut corresponds with the 2^30-bit, ...
    /// </summary>
    public struct GeomBinTreeBranchCode : IComparable<GeomBinTreeBranchCode> {

        /// <summary>
        /// variable to store the branch code
        /// </summary>
        public uint Code;

        /// <summary>
        /// finds the branch (denoted by the branch code of the binary geometric tree) of some point 
        /// </summary>
        /// <param name="bb"></param>
        /// <param name="pt"></param>
        /// <returns></returns>
        public static GeomBinTreeBranchCode CreateFormPoint(BoundingBox bb, params double[] pt) {
            int D = pt.Length;
            if (D != bb.D)
                throw new ArgumentException("mismatch between dimension of bounding box and point.");
            if (!bb.Contains(pt))
                throw new ArgumentException("point cannot be outside of given bounding box.");

            GeomBinTreeBranchCode cd;
            cd.Code = 0;

            BoundingBox red = (BoundingBox)bb.Clone();
            double sep = (red.Min[0] + red.Max[0]) * 0.5;

            //// we compute the most-significant really arithmetically (=2^31), to 
            //// be get the same algorithms on little and big endian machines, however they are weird.
            uint mxx = 0x80000000;
            //for (int i = 0; i < 30; i++)
            //    mxx *= 2;

            int cur_d = 0;
            for (int level = 31; level >= 0; level--) {

                // decide in which side of the bounding box our point is located
                sep = (red.Min[cur_d] + red.Max[cur_d]) * 0.5;
                bool side = false;
                if (pt[cur_d] >= sep)
                    side = true;

                // reduce bounding box
                if (side) {
                    red.Min[cur_d] = sep;
                } else {
                    red.Max[cur_d] = sep;
                }

                // note branch
                if (side) {
                    cd.Code |= mxx;
                }
                mxx /= 2;

                // rotate spatial dimension
                cur_d++;
                if (cur_d >= D)
                    cur_d = 0;
            }


            BoundingBox test = new BoundingBox(D);
            bb.LeaveBoxFromCode(test, cd);
            if (!test.ContainsStrict(pt))
                throw new ApplicationException();

            return cd;
        }

        #region IComparable<GeomBinGeomTreeCode> Members

        /// <summary>
        /// compares the <see cref="Code"/>
        /// </summary>
        public int CompareTo(GeomBinTreeBranchCode other) {
            long diff = (long)this.Code - (long)other.Code;  // cast to long is the safest way
            if (diff < 0)
                return -1;
            if (diff > 0)
                return 1;
            return 0;
        }


        #endregion

        /// <summary>
        /// comparison operator
        /// </summary>
        public static bool operator <(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b) {
            return (a.Code < b.Code);
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public static bool operator >(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b) {
            return (a.Code > b.Code);
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public static bool operator <=(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b) {
            return (a.Code <= b.Code);
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public static bool operator >=(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b) {
            return (a.Code >= b.Code);
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public static bool operator ==(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b) {
            return (a.Code == b.Code);
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public static bool operator !=(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b) {
            return (a.Code == b.Code);
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public override bool Equals(object obj) {
            try {
                return (((GeomBinTreeBranchCode)obj).Code == this.Code);
            } catch (InvalidCastException) {
                return false;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="commonBits">
        /// on exit, the number of common branches in codes <paramref name="a"/> and <paramref name="b"/>
        /// </param>
        /// <returns></returns>
        public static GeomBinTreeBranchCode Combine(GeomBinTreeBranchCode a, GeomBinTreeBranchCode b, out int commonBits) {
            GeomBinTreeBranchCode ret;
            ret.Code = 0;

            uint mask = 0x80000000;
            for (commonBits = 0; commonBits < 32; commonBits++) {
                uint tst_a = mask & a.Code;
                uint tst_b = mask & b.Code;

                mask = mask >> 1;

                if (tst_a == tst_b)
                    ret.Code |= tst_a;
                else
                    return ret;
            }

            return ret;
        }

        /// <summary>
        /// comparison operator
        /// </summary>
        public override int GetHashCode() {
            return (int)this.Code;
        }
    }

    /// <summary>
    /// This represents the data that identifies a bounding box in a binary geometric tree;
    /// </summary>
    public struct BoundingBoxCode {

        /// <summary>
        /// Branch code of the bounding box; only the last <see cref="SignificantBits"/> bits are 
        /// of interest (here, the <em>last bit</em> is the most significant one, i.e. the bit that corresponds with 2^31).
        /// </summary>
        public GeomBinTreeBranchCode Branch;

        /// <summary>
        /// Number of subdivisions of the Root Bounding Box;
        /// </summary>
        /// <remarks>
        /// number of significant bits of the branch code <see cref="Branch"/>; <br/>
        /// If this is 0, this structure represents the root box of the binary geometric tree.<br/>
        /// If 1, the most significant bit of <see cref="Branch"/> will be taken into account,
        /// if 2, the two most significant bits of <see cref="Branch"/> will be ...
        /// </remarks>
        public uint SignificantBits;

        /// <summary>
        /// Returns the bitmask that marks the significant bits of the branch code <see cref="Branch"/>;
        /// </summary>
        public uint GetBitMask() {
            uint M0 = 0x80000000;

            uint Ret = 0;
            for (int i = 0; i < SignificantBits; i++) {
                Ret |= M0;
                M0 = M0 >> 1;
            }

            return Ret;
        }

        ///// <summary>
        ///// 
        ///// </summary>
        ///// <param name="container"></param>
        ///// <param name="bbout"></param>
        //public void GetBox(BoundingBox container, BoundingBox bbout) {
        //    container.SubBoxFromCode
        //}

        /// <summary>
        /// true, if the branch notated by <paramref name="code"/> is contained in this box
        /// </summary>
        public bool IsInside(GeomBinTreeBranchCode code) {
            uint mask = GetBitMask();
            return ((this.Branch.Code & mask) == (code.Code & mask));
        }

        /// <summary>
        /// true, if the box notated by <paramref name="code"/> is fully contained in this box
        /// </summary>
        public bool IsInside(BoundingBoxCode code) {
            if (code.SignificantBits < this.SignificantBits)
                return false;
            return IsInside(code.Branch);
        }

        /// <summary>
        /// gets the code of either the 'left' or the 'right' sub-box of this box
        /// </summary>
        /// <param name="leftOrRight">
        /// if false, the 'left' branch, i.e. the box with lower coordinate values;<br/>
        /// if false, the 'right' branch, i.e. the box with higher coordinate values;<br/>
        /// </param>
        /// <returns></returns>
        public BoundingBoxCode GetSubBox(bool leftOrRight) {
            if (this.SignificantBits >= 32)
                throw new ApplicationException("maximum tree depth reached.");

            BoundingBoxCode Ret = this;
            Ret.SignificantBits++;
            if (leftOrRight) {
                // right branch
                uint upper = 0x80000000 >> ((int)this.SignificantBits);
                Ret.Branch.Code |= upper;
            }
            return Ret;
        }
    }

    /// <summary>
    /// 2D Rotation about the origin
    /// </summary>
    public class Rotation2D {

        /// <summary>
        /// Initializes the Rotation Matrix
        /// </summary>
        /// <param name="phi"></param>
        public Rotation2D(double phi) {
            cos = Math.Cos(phi);
            sin = Math.Sin(phi);
        }

        private double cos;
        private double sin;

        /// <summary>
        /// Rotate the Coordinate
        /// </summary>
        /// <param name="X">input coordinate</param>
        /// <returns>rotated coordinate </returns>
        public double[] Transform(double[] X) {
            if (X.Length != 2) { throw new ArithmeticException("Only 2D!"); };
            double x = X[0] * cos - X[1] * sin;
            double y = X[0] * sin + X[1] * cos;
            return new double[2] { x, y };
        }

        /// <summary>
        /// Rotate the Coordinate about the angle phi
        /// </summary>
        /// <param name="X">input coordinate</param>
        /// <param name="phi">rotation angle</param>
        /// <returns>rotated coordinate </returns>
        public double[] Transform(double[] X, double phi) {
            if (X.Length != 2) { throw new ArithmeticException("Only 2D!"); };
            double x = X[0] * Math.Cos(phi) - X[1] * Math.Sin(phi);
            double y = X[0] * Math.Sin(phi) + X[1] * Math.Cos(phi);
            return new double[2] { x, y };
        }
    }

    public static class SimpleGeoTools {
        /// <summary>
        /// Calculates the cross product of two vectors in 2D
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double CrossProduct2D(double[] a, double[] b) {
            return a[0] * b[1] - a[1] * b[0];
        }

        /// <summary>
        /// Calculates the shortest distance from a point to a line
        /// </summary>
        /// <param name="X">The point of interest</param>
        /// <param name="anyPointOnLine">Any point on the line</param>
        /// <param name="directionVector">Direction vector of the line</param>
        /// <returns>Distance to the line</returns>
        public static double DistanceFromPointToLine(double[] X, double[] anyPointOnLine, double[] directionVector) {
            double[] X_minus_pointOnLine = new double[] { X[0] - anyPointOnLine[0], X[1] - anyPointOnLine[1] };
            double distance = CrossProduct2D(directionVector, X_minus_pointOnLine) / Math.Sqrt(Math.Pow(directionVector[0], 2) + Math.Pow(directionVector[1], 2));

            return -distance;
        }

        /// <summary>
        /// Smoothes, e.g., an initial condition over the range h/p using
        /// a tanh profile
        /// </summary>
        /// <param name="distance">The distance to e.g. a shock</param>
        /// <param name="cellSize">The cell size h</param>
        /// <param name="dgDegree">The DG polynomial order</param>
        /// <param name="smoothingWidth">Width of the smoothing</param>
        /// <returns></returns>
        public static double SmoothJump (double distance, double cellSize, int dgDegree, double smoothingWidth) {
            // Smoothing should be in the range of h/p
            double maxDistance = smoothingWidth * cellSize / Math.Max(dgDegree, 1);

            return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
        }
    }
}