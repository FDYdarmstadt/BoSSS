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
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Foundation {

    /// <summary>
    /// Extensions for handling non-orthonormal, DG approximation of (absolute) degree <see cref="Degree"/>.
    /// Use with caution!!!
    /// </summary>
    public partial class Basis {

        /// <summary>
        /// creates a new basis, with a nonstandard polynomial set
        /// </summary>
        /// <param name="_grd">the <see cref="GridData"/> that stores grid information</param>
        /// <param name="degree">highest polynomial degree of basis polynomials</param>
        /// <param name="polys">List of polynomials</param>
        public Basis(IGridData _grd, int degree, List<PolynomialList> polys) {
            List<PolynomialList> P = polys;
            var grd = ((GridData)_grd).CopyWithNewBasis(); // Make sure the GridData is using this custom polynomial set
            grd.ChefBasis.SetCustomPolynomials(P.Max(p => p.MaxAbsoluteDegree), P.ToArray());
            this.m_GridDat = grd;

            this.Polynomials = P;
            this.IsOrthonormal = false;

            this.MinimalLength = Polynomials.Min(pl => pl.Count);
            this.MaximalLength = Polynomials.Max(pl => pl.Count);
            this.Degree = this.Polynomials.Max(pl => pl.MaxAbsoluteDegree);

        }

    }   

    
    /// <summary>
    /// Provides routines to transform coefficients between different Polynomial bases
    /// </summary>
    public abstract class BasisTransformator {

        protected Basis m_origin;

        /// <summary>
        /// Origin Basis of this BasisTransformer
        /// </summary>
        public Basis Origin {
            get {
                return m_origin;
            }
        }

        protected Basis m_dest;

        /// <summary>
        /// Destination Basis of this BasisTransformer
        /// </summary>
        public Basis Destination { 
            get {
                return m_dest;
            } 
        }

        protected double m_offset; // relative value, the edges of the refelements are elongated by

        /// <summary>
        /// the DG basis functions for the reference elements
        /// <list type="bullet">
        ///     <item>1st index: reference element index (see <see cref="BoSSS.Foundation.Grid.IGeometricalCellsData.RefElements"/>)</item>
        ///     <item>2nd index: Polynomial index</item>
        /// </list>
        /// </summary>
        public List<PolynomialList> Polynomials {
            get;
            protected set;
        }

        public BasisTransformator(Basis origin, double offset) {
            m_origin = origin;
            m_offset = offset;

            GeneratePolynomials();
            m_dest = new Basis(m_origin.GridDat, m_origin.Degree, this.Polynomials);

            ConstructTransform();
            ConstructInverseTransform();
#if DEBUG
            TestTransformQuality();
#endif
        }

        /// <summary>
        /// Origin -> Destination
        /// - index: reference element index
        /// </summary>
        public MultidimensionalArray[] Origin2Dest {
            private set;
            get;
        }

        /// <summary>
        /// Destination -> Origin
        /// - index: reference element index
        /// </summary>
        public MultidimensionalArray[] Dest2Origin {
            private set;
            get;
        }

        /// <summary>
        /// From origin basis to destination basis $` T = \langle B_i,B_j \rangle^{-1} \langle B_j, P_k \rangle = M^{-1} \cdot S $`
        /// </summary>
        protected virtual void ConstructTransform() {
            var Krefs = ((GridCommons)m_origin.GridDat.Grid).RefElements;
            int N = Krefs.Length;

            MultidimensionalArray[] S = new MultidimensionalArray[N];
            MultidimensionalArray[] M = new MultidimensionalArray[N];
            Origin2Dest = new MultidimensionalArray[N];

            for (int iKref = 0; iKref < N; iKref++) {
                int I = m_dest.Polynomials[iKref].Count;
                int J = m_origin.Polynomials[iKref].Count;
                var Kref = Krefs[iKref];
                var rule = Kref.GetQuadratureRule(Math.Min(m_dest.Degree + m_origin.Degree, Kref.HighestKnownOrder));
                S[iKref] = MultidimensionalArray.Create(I, J);
                M[iKref] = MultidimensionalArray.Create(I, I);
                Origin2Dest[iKref] = MultidimensionalArray.Create(I, J);

                var BsvalOrg = m_origin.Evaluate(rule.Nodes);
                var BsvalDst = m_dest.Evaluate(rule.Nodes);
                var BsvalDstXwgt_T = MultidimensionalArray.Create(I, rule.NoOfNodes);

                BsvalDstXwgt_T.Multiply(1.0, BsvalDst, rule.Weights, 0.0, "nk", "kn", "k");
                M[iKref].GEMM(1.0, BsvalDstXwgt_T, BsvalDst, 0.0);
                S[iKref].GEMM(1.0, BsvalDstXwgt_T, BsvalOrg, 0.0);


             

                // Calculate Transformation
                if (!m_dest.IsOrthonormal) {
                    //M[iKref].InvertSymmetrical();
                    //Origin2Dest[iKref].DGEMM(1.0, M[iKref], S[iKref], 0.0);
                    M[iKref].SolveSymmetricEx(Origin2Dest[iKref], S[iKref]);
                } else {
                    Origin2Dest[iKref].Acc(1.0, S[iKref]); // in this case M is unity
                }

                /*
                 * Test code:

                double[] coordsOrg = new double[m_origin.Length];
                double[] coordsDst = new double[m_dest.Length];

                coordsOrg.FillRandom();
                Origin2Dest[iKref].GEMV(1.0, coordsOrg, 0.0, coordsDst);
                double[] valuesOrg = new double[rule.Nodes.NoOfNodes];
                double[] valuesDst = new double[rule.Nodes.NoOfNodes];

                BsvalOrg.GEMV(1.0, coordsOrg, 0.0, valuesOrg);
                BsvalDst.GEMV(1.0, coordsDst, 0.0, valuesDst);

                double err = valuesOrg.L2Dist(valuesDst);
                Console.WriteLine("Error is " + err);
                */
            }
        }


        /// <summary>
        /// From destination basis to origin basis: $` T = \langle P_i,P_j \rangle^{-1} \langle P_j, B_k \rangle $`
        /// </summary>
        protected virtual void ConstructInverseTransform() {
            var Krefs = ((GridCommons)m_origin.GridDat.Grid).RefElements;
            int N = Krefs.Length;

            MultidimensionalArray[] S = new MultidimensionalArray[N];
            MultidimensionalArray[] M = new MultidimensionalArray[N];
            Dest2Origin = new MultidimensionalArray[N];

            for (int iKref = 0; iKref < N; iKref++) {
                int I = m_dest.Polynomials[iKref].Count;
                int J = m_origin.Polynomials[iKref].Count;
                var Kref = Krefs[iKref];
                var rule = Kref.GetQuadratureRule(Math.Min(m_dest.Degree + m_origin.Degree, Kref.HighestKnownOrder));
                S[iKref] = MultidimensionalArray.Create(J, I);
                M[iKref] = MultidimensionalArray.Create(J, J);
                Dest2Origin[iKref] = MultidimensionalArray.Create(J, I);

                // not well suited for 3D and higher Polynomialdegrees!!!!, condition number of mass matrix becomes very high
                var B = m_dest.Evaluate(rule.Nodes);
                var P = m_origin.Evaluate(rule.Nodes);

                MultidimensionalArray PVal;
                MultidimensionalArray BVal;
                for (int k = 0; k < rule.Nodes.NoOfNodes; k++) {
                    double weight = rule.Weights[k];
                    PVal = P.ExtractSubArrayShallow(k, -1);
                    BVal = B.ExtractSubArrayShallow(k, -1);
                    // Quadratur M
                    {   
                        if(!m_origin.IsOrthonormal)
                            M[iKref].Multiply(weight, PVal, PVal, 1.0, "ij", "i", "j");
                    }
                    // Quadratur S
                    {
                        S[iKref].Multiply(weight, PVal, BVal, 1.0, "ij", "i", "j");
                    }
                }

                // Calculate Transformation
                if (!m_origin.IsOrthonormal) {
                    //M[iKref].InvertSymmetrical();
                    //Dest2Origin[iKref].DGEMM(1.0, M[iKref], S[iKref], 0.0);
                   
                    M[iKref].SolveSymmetricEx(Dest2Origin[iKref], S[iKref]);
                } else {
                    Dest2Origin[iKref].Acc(1.0, S[iKref]); // in this case M is unity
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        protected virtual void TestTransformQuality() {
            var Krefs = ((GridCommons)m_origin.GridDat.Grid).RefElements;
            int N = Krefs.Length;

            for (int iKref = 0; iKref < N; iKref++) {
                int I = m_dest.Polynomials[iKref].Count;
                int J = m_origin.Polynomials[iKref].Count;
                MultidimensionalArray TinvT = MultidimensionalArray.Create(J, J);
                TinvT.GEMM(1.0, Dest2Origin[iKref], Origin2Dest[iKref], 0.0);
                TinvT.AccEye(-1.0);
                /*
                Console.WriteLine("============================================================");
                Console.WriteLine("Testing Quality of BasisTransformation, degree:{0}", Math.Max(m_dest.Degree, m_origin.Degree));
                Console.WriteLine("     |T^-1*T - I|_2:         {0}", TinvT.L2Norm());
                Console.WriteLine("     |T|_2:                  {0}", Origin2Dest[iKref].L2Norm());
                Console.WriteLine("     |T^-1|_2:               {0}", Dest2Origin[iKref].L2Norm());
                Console.WriteLine("     k_2(T)=|T^-1|_2*|T|_2:  {0}", Dest2Origin[iKref].L2Norm()*Origin2Dest[iKref].L2Norm());
                Console.WriteLine("============================================================");
                */
            }
        }

        protected abstract void GeneratePolynomials();
    }


    /// <summary>
    /// Calculates the Transformation from and to Bernstein Basis Coefficients
    /// </summary>
    public class BernsteinTransformator : BasisTransformator{        

        /// <summary>
        /// Indices of Coefficients lying on faces of the refelement
        ///     1st index: Refelement
        ///     2nd index: Face
        /// </summary>
        public List<int>[][] FaceCoefficients {
            get;
            private set;
        }

        public BernsteinTransformator(Basis origin, double offset = 0.0) : base(origin, offset) {           
        }
        
        private static double nchoosek(int n, int i) {
            if (i > n || i <0 || n<0) throw new ArgumentException();
            if (n == i || i == 0) return 1;

            if (i > n - i) {
                return nchoosek(n, n - i);
            }
            return nchoosek(n - 1, i) + nchoosek(n - 1, i - 1);
        }

        // returns the i-th bernstein polynomial of order n in the d-th variable
        private Polynomial Bernstein(int index, int order, int d, int D, double a = -1.0, double b = 1.0) {
            if (index > order) throw new ArgumentOutOfRangeException();
            Polynomial Poly; // <Exponent, Coefficient>
            Polynomial P1 = new Polynomial(); //(b-u)^(n-i)
            Polynomial P2 = new Polynomial(); //(u-a)^(i)

            double factor = 1.0 / Math.Pow((b-a), order) * nchoosek(order, index);

            int[] exp = new int[D];
            for (int k = 0; k <= order - index; k++) {
                exp.Clear();
                exp[d] = k;
                P1.AddCoeff(nchoosek(order - index, k) * Math.Pow(-1.0, k) * Math.Pow(b, order - index - k), exp);
            }
            for (int j = 0; j <= index; j++) {
                exp.Clear();
                exp[d] = j;
                P2.AddCoeff(nchoosek(index, j) * Math.Pow(-a, index - j), exp);
            }

            Poly = P1 * P2;
            Poly *= factor;

            return Poly;
        }

        protected override void GeneratePolynomials() {
            int D = m_origin.GridDat.SpatialDimension;
            if (D < 1 && D > 3) throw new NotSupportedException("Wrong spatial dimension, only works with 1, 2 or 3");
            int p = m_origin.Degree;
            var Krefs = ((GridCommons)m_origin.GridDat.Grid).RefElements;
            FaceCoefficients = new List<int>[Krefs.Length][];

            if (this.Polynomials.IsNullOrEmpty()) this.Polynomials = new List<PolynomialList>();

            for (int iKref = 0; iKref < Krefs.Length; iKref++) {
                var Kref = Krefs[iKref];
                FaceCoefficients[iKref] = new List<int>[Kref.NoOfFaces];
                for (int r = 0; r < Kref.NoOfFaces; r++) {
                    FaceCoefficients[iKref][r] = new List<int>();
                }
                if (Kref is BoSSS.Foundation.Grid.RefElements.Cube || Kref is BoSSS.Foundation.Grid.RefElements.Square) {

                    List<Polynomial> polys = new List<Polynomial>();
                    double[] a = new double[D]; // Minima per dimension
                    double[] b = new double[D]; // Maximum per dimension

                    NodeSet Vertices = Kref.Vertices;
                    NodeSet OffsetVertices = new NodeSet(Vertices.RefElement, Vertices.NoOfNodes, D, false);
                    for (int d = 0; d < D; d++) {
                        for (int n = 0; n < Vertices.NoOfNodes; n++) {
                            OffsetVertices[n, d] = Vertices[n, d];
                            OffsetVertices[n, d] += m_offset * (Vertices[n, d] - Kref.Center[0, d]);
                        }
                    }
                    OffsetVertices.LockForever();
                    for (int d = 0; d < D; d++) {
                        a[d] = OffsetVertices.ExtractSubArrayShallow(-1, d).To1DArray().Min();
                        b[d] = OffsetVertices.ExtractSubArrayShallow(-1, d).To1DArray().Max();
                    }

                    int[] i = new int[D + 1];
                    i.Clear();
                    i[D] = -1;
                    int k = 0;
                    while (i[D] == -1) {
                        Polynomial P = new Polynomial();
                        P.AddCoeff(1.0, new int[D]);

                        for (int d = 0; d < D; d++) {
                            P *= Bernstein(i[d], p, d, D, a[d], b[d]);
                        }
                        polys.Add(P);

                        i[0]++;
                        while (i[k] == p + 1) {
                            i[k] = 0;
                            i[++k]++;
                            if (i[k] != p + 1)
                                k = 0;
                        }
                    }

                    var PolyList = new PolynomialList(polys);
                    Polynomials.Add(PolyList);

                    // collect vertices on scaled refelement face
                    NodeSet FaceVertices = Kref.FaceCenters;
                    NodeSet OffsetFaceVertices = new NodeSet(FaceVertices.RefElement, FaceVertices.NoOfNodes, D, false);
                    for (int d = 0; d < D; d++) {
                        for (int n = 0; n < FaceVertices.NoOfNodes; n++) {
                            OffsetFaceVertices[n, d] = FaceVertices[n, d];
                            OffsetFaceVertices[n, d] += m_offset * (FaceVertices[n, d] - Kref.Center[0, d]);
                        }
                    }
                    OffsetFaceVertices.LockForever();
                    MultidimensionalArray R = MultidimensionalArray.Create(Kref.NoOfFaces, PolyList.Count);
                    PolyList.Evaluate(OffsetFaceVertices, R);
                    for (int j = 0; j < PolyList.Count; j++) {
                        for (int r = 0; r < Kref.NoOfFaces; r++) {
                            if (Math.Abs(R[r, j]) > 1e-12) FaceCoefficients[iKref][r].Add(j);
                        }
                    }

                    for (int r = 0; r < Kref.NoOfFaces; r++) {
                        if (FaceCoefficients[iKref][r].Count != (int)Math.Pow(p + 1, D - 1)) 
                            Console.WriteLine("Warning wrong number of coinciding coefficients on face {0}", r);
                    }

                } else if (Kref is BoSSS.Foundation.Grid.RefElements.Tetra || Kref is BoSSS.Foundation.Grid.RefElements.Triangle || Kref is BoSSS.Foundation.Grid.RefElements.Line) {
                    // these elements are all simplices, we use barycentric coordinates to construct the bernstein polynomials
                    List<Polynomial> polys = new List<Polynomial>();

                    NodeSet Vertices = Kref.Vertices;
                    NodeSet OffsetVertices = new NodeSet(Vertices.RefElement, Vertices.NoOfNodes, D, false);
                    for (int d = 0; d < D; d++) {
                        for (int n = 0; n < Vertices.NoOfNodes; n++) {
                            OffsetVertices[n, d] = Vertices[n, d];
                            OffsetVertices[n, d] += m_offset * (Vertices[n, d] - Kref.Center[0, d]);
                        }
                    }
                    OffsetVertices.LockForever();

                    // compute barycentric coordinates and express as polynomials
                    Polynomial[] Cartesian = new Polynomial[D + 1];
                    for (int d = 0; d <= D; d++) {
                        int[] exp = new int[D];
                        exp.Clear();
                        if (d < D) exp[d] = 1;
                        var P = new Polynomial();
                        P.AddCoeff(1.0, exp);
                        Cartesian[d] = P;
                    }
                    Polynomial[] Barycentric = new Polynomial[D + 1];
                    MultidimensionalArray V = MultidimensionalArray.Create(D + 1, D + 1);
                    V.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { D - 1, D }).Acc(1.0, OffsetVertices.TransposeTo());
                    V.ExtractSubArrayShallow(D, -1).SetAll(1.0);
                    V.InvertInPlace();
                    for (int d = 0; d <= D; d++) {
                        var P = new Polynomial();
                        P.AddCoeff(0.0, new int[D]); // we need to initialize the polynomial to be of Spatial Dimension D, this is controlled via the Length of the exponents
                        for (int e = 0; e <= D; e++) {
                            P += V[d, e] * Cartesian[e];
                        }
                        Barycentric[d] = P;
                    }

                    int[] i = new int[D + 1];
                    i.Clear();
                    int k = 0;
                    int cancel = 0;
                    int count = 0;
                    while (cancel == 0) {
                        i[D] = p - i.Take(D).Sum();
                        double factor = 1.0;
                        Polynomial P = new Polynomial();
                        P.AddCoeff(1.0, new int[D]);
                        for (int d = 1; d < p+1; d++) { factor *= d; } // p!
                        for (int d = 0; d <= D; d++) {                            
                            for (int j = 1; j <= i[d]; j++) {
                                factor *= 1.0 / j;
                                P *= Barycentric[d];
                            }
                        }
                        P *= factor;
                        polys.Add(P);

                        // as we use barycentric coordinates this simpler routine can be used to obtain coefficients on faces
                        int FaceSum = Enumerable.Range(0, D+1).Sum();
                        for (int d = 0; d <= D; d++) {
                            if (i[d] == 0) {
                                var coincidingFaces = Kref.VertexIndicesToFaces[d];
                                int oppositeFace = FaceSum - coincidingFaces.Sum();
                                FaceCoefficients[iKref][oppositeFace].Add(count); // add the index of the coefficient positioned on that face 
                            }
                        }

                        //for (int d = 0; d <= D; d++) {
                        //    Console.Write("i{0}:{1}, ", d, i[d]);
                        //}
                        //Console.WriteLine();

                        count++;
                        i[0]++;
                        // The sum over all indices should equal the requested degree, therefore we need to upper limit of the loop
                        while (i[k] == (p + 1) - i.Take(D).Skip(k + 1).Sum()) {
                            if (k == D - 1 && i[k] == p + 1) {
                                cancel = 1;
                                break;
                            }

                            i[k] = 0;
                            i[++k]++;
                            if (i[k] != (p + 1) - i.Take(D).Skip(k + 1).Sum()) {                                
                                k = 0;
                            }                            
                        }
                    }

                    var PolyList = new PolynomialList(polys);
                    Polynomials.Add(PolyList);

                    for (int r = 0; r < Kref.NoOfFaces; r++) {
                        if (D == 2) {
                            if (FaceCoefficients[iKref][r].Count != p + 1) Console.WriteLine("Warning wrong number of coinciding coefficients on face {0}", r);
                        } else if (D == 3) {
                            if (FaceCoefficients[iKref][r].Count != (p + 1) * (p + 2) / 2) Console.WriteLine("Warning wrong number of coinciding coefficients on face {0}", r);
                        } else if (D == 1) {
                            if (FaceCoefficients[iKref][r].Count != 1) Console.WriteLine("Warning wrong number of coinciding coefficients on face {0}", r);
                        }
                    }
                    
                } else {
                    throw new NotSupportedException();
                }
            }

            
        }
    }

}
