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
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// Base class that implements/defines features which all reference
    /// elements (Triangles, Lines, Tetras which build up the mesh) have in
    /// common.
    /// </summary>
    /// <remarks>
    /// Geometrically, the reference element is defined to be the convex hull of
    /// <see cref="Vertices"/>.
    /// </remarks>
    public abstract partial class RefElement : IEquatable<RefElement> {

        /// <summary>
        /// See <see cref="Volume"/>.
        /// </summary>
        private double m_Volume = -1;

        /// <summary>
        /// the <see cref="SpatialDimension"/>-dimensional volume of this
        /// reference element
        /// </summary>
        public double Volume {
            get {
                var qr = this.GetQuadratureRule(0);
                if (m_Volume < 0) {
                    m_Volume = qr.Weights.Sum();
                }
                return m_Volume;
            }
        }

        /*

        /// <summary>
        /// finds the polynomial in <see cref="OrthonormalPolynomials"/> with
        /// matching Guid <paramref name="g"/>
        /// </summary>
        public Polynomial GetPolynomialByGuid(Guid g) {
            foreach (Polynomial p in this.OrthonormalPolynomials) {
                if (p.MyGuid.Equals(g))
                    return p;
            }

            throw new ArgumentException("unknown Polynomial GUID");
        }
         */

        /// <summary>
        /// maximum polynomial degree that is available for a
        /// <see cref="Basis"/>-object
        /// </summary>
        public int HighestSupportedPolynomialDegree {
            get {
                int deg = 0;
                foreach (Polynomial p in this.OrthonormalPolynomials) {
                    deg = Math.Max(deg, p.AbsoluteDegree);
                }
                return deg;
            }
        }

        /// <summary>
        /// must be initialized by constructor of derived class;
        /// </summary>
        protected int m_NoOfFaces = -123445;

        /// <summary>
        /// the number of edges, or neighbors;
        /// </summary>
        public int NoOfFaces {
            get {
                return m_NoOfFaces;
            }
        }

        NodeSet m_Center;

        /// <summary>
        /// Center (of gravity) of this reference element.
        /// </summary>
        public NodeSet Center {
            get {
                if(m_Center == null) {
                    m_Center = new NodeSet(this, new double[Math.Max(1, this.SpatialDimension)]);
                }
                return m_Center;
            }
        }

        /// <summary>
        /// Number of vertices this simplex consists of
        /// </summary>
        public int NoOfVertices {
            get {
                return Vertices.GetLength(0);
            }
        }

        /// <summary>
        /// must be initialized by constructor of derived class;
        /// </summary>
        protected RefElement m_FaceRefElement;


        /// <summary>
        /// a hack to ensure certain reference-equalities.
        /// </summary>
        internal void HackFaceRefElement(RefElement FaceRefelement) {
            if(this.m_FaceRefElement.GetType() != FaceRefelement.GetType())
                throw new ArgumentException();
            this.m_FaceRefElement = FaceRefelement;
        }


        /// <summary>
        /// The reference element that forms the faces of this element.
        /// </summary>
        public RefElement FaceRefElement {
            get {
                return m_FaceRefElement;
            }
        }


        /// <summary>
        /// transforms some vertices from local coordinates on one face to the
        /// local coordinate system of this simplex
        /// </summary>
        /// <param name="FaceIndex">specifies the edge</param>
        /// <param name="InpVerticesOnFace">
        /// Input; Vertices in the local coordinate system of the face;
        /// 1st index: vertex index; 2nd index: spatial coordinate index, range
        /// is 0 (including) to max(<see cref="SpatialDimension"/>-1, 1)
        /// (excluding)
        /// </param>
        /// <param name="OutVerticesInSimplex">
        /// On exit, the <paramref name="InpVerticesOnFace"/> transformed to
        /// the local coordinate system of this simplex;
        /// 1st index: vertex index;
        /// 2nd index: spatial coordinate vector, in the range of 0 (including)
        /// to <see cref="SpatialDimension"/> (excluding);
        /// </param>
        public abstract void TransformFaceCoordinates(int FaceIndex, MultidimensionalArray InpVerticesOnFace, MultidimensionalArray OutVerticesInSimplex);

        /// <summary>
        /// Transforms some vectors from the local coordinate system on a face
        /// to the local coordinate system of the simplex
        /// </summary>
        /// <param name="FaceIndex">Specifies the face</param>
        /// <param name="gradientsOnFace">
        /// Input; Vectors in the local coordinate system of the face;
        /// <list type="bullet">
        ///     <item>
        ///         1st index: Vector index
        ///     </item>
        ///     <item>
        ///         2nd index: Spatial coordinate index; Range is 0 (including)
        ///         to max(<see cref="SpatialDimension"/>-1, 1) (excluding)
        ///     </item>
        /// </list>
        /// </param>
        /// <param name="gradientsInSimplex">
        /// Output; Vectors in the local coordinate system of the element;
        /// <list type="bullet">
        ///     <item>
        ///         1st index: Vector index    
        ///     </item>
        ///     <item>
        ///         2nd index: Spatial coordinate index; Range is 0 (including)
        ///         to <see cref="SpatialDimension"/> (excluding)
        ///     </item>
        /// </list>
        /// </param>
        public void TransformFaceVectors(int FaceIndex, MultidimensionalArray gradientsOnFace, MultidimensionalArray gradientsInSimplex) {
            AffineTrafo trafo = GetFaceTrafo(FaceIndex).CloneAs();
            trafo.Affine.ClearEntries();
            trafo.Transform(gradientsOnFace, gradientsInSimplex);
        }

        /// <summary>
        /// see <see cref="GetEmbeddedFaceTrafo"/>
        /// </summary>
        AffineTrafo[] m_EmbeddedFaceTransformation;

        /// <summary>
        /// returns an affine-linear transformation which transforms from the
        /// local coordinate system of the <paramref name="FaceIndex"/>-th face
        /// to the local coordinate system of the simplex.
        /// </summary>
        /// <remarks>
        /// Let D be the spatial dimension of the simplex, then the spatial
        /// dimension of the face is D-1; but, for the transformation, the face
        /// coordinate system is embedded into
        /// \f$ \mathbb{R}^{D}\f$  and we additionally map
        /// the D-th standard basis vector to
        /// \f$ \vec{c}_e + \vec{n}_e\f$  (center of edge
        /// \f$ e\f$  plus normal), i.e.
        /// \f$ 
        /// \mathbb{R}^{D} \ni (0,\ldots,0,1) \mapsto \vec{c}_e + \vec{n}_e
        /// \f$ .
        /// </remarks>
        public AffineTrafo GetEmbeddedFaceTrafo(int FaceIndex) {
            if (m_EmbeddedFaceTransformation == null) {
                m_EmbeddedFaceTransformation = new AffineTrafo[NoOfFaces];

                int[,] EdgIdx = FaceToVertexIndices;
                var vtx_Cod = this.Vertices;
                var vtx_Dom = this.FaceRefElement.Vertices;

                int D = this.SpatialDimension;
                int D_dom = this.FaceRefElement.SpatialDimension;
                Debug.Assert(D_dom == D - 1);

                var centers = this.FaceCenters;
                var Normals = this.FaceNormals;

                double[,] Cod = new double[D + 1, D];
                double[,] Dom = new double[D + 1, D];

                for(int i = 0; i < NoOfFaces; i++) { // loop over all edges.

                    // map Edge Simplex to volume simplex
                    for(int l = 0; l < D; l++) {
                        for(int d = 0; d < D_dom; d++)
                            Dom[l, d] = vtx_Dom[l, d];
                        Dom[l, D - 1] = 0;

                        for(int d = 0; d < D; d++)
                            Cod[l, d] = vtx_Cod[EdgIdx[i, l], d];
                    }

                    // map (0,0,1) -> Edge_center + Edge_normal
                    for(int d = 0; d < D_dom; d++)
                        Dom[D, d] = 0;
                    Dom[D, D - 1] = 1;
                    for(int d = 0; d < D; d++)
                        Cod[D, d] = centers[i, d] + Normals[i, d];


                    // compute trafo:
                    m_EmbeddedFaceTransformation[i] = AffineTrafo.FromPoints(Dom, Cod);
                }

            }
            return m_EmbeddedFaceTransformation[FaceIndex];
        }

        /// <summary>
        /// see <see cref="GetInverseEmbeddedFaceTrafo"/>
        /// </summary>
        private AffineTrafo[] m_InverseEmbeddedFaceTransformation;

        /// <summary>
        /// the transformation from the coordinate system of the simplex to
        /// the coordinate system of face <paramref name="FaceIndex"/>.
        /// </summary>
        /// <remarks>
        /// Let D be the spatial dimension of the simplex, then the spatial
        /// dimension of the face is D-1; but, for the transformation, the face
        /// coordinate system is embedded into
        /// \f$ \mathbb{R}^{D}\f$  and we additionally map
        /// \f$ \vec{c}_e + \vec{n}_e\f$  (center of edge
        /// \f$ e\f$  plus normal) to the D-th standard
        /// basis vector, i.e.
        /// \f$ 
        /// \vec{c}_e + \vec{n}_e \mapsto (0,\ldots,0,1) \in \mathbb{R}^{D}
        /// \f$ .
        /// </remarks>
        public AffineTrafo GetInverseEmbeddedFaceTrafo(int FaceIndex) {
            if (m_InverseEmbeddedFaceTransformation == null) {
                m_InverseEmbeddedFaceTransformation = new AffineTrafo[this.NoOfFaces];

                for (int l = 0; l < this.NoOfFaces; l++) {
                    var Et = this.GetEmbeddedFaceTrafo(l);
                    m_InverseEmbeddedFaceTransformation[l] = Et.Invert();
                }
            }
            return m_InverseEmbeddedFaceTransformation[FaceIndex];
        }

        /// <summary>
        /// <see cref="GetInverseFaceTrafo"/>
        /// </summary>
        private AffineTrafo[] m_InverseFaceTransformation;

        /// <summary>
        /// the transformation from the coordinate system of the simplex to
        /// the coordinate system of face <paramref name="FaceIndex"/>.
        /// </summary>
        public AffineTrafo GetInverseFaceTrafo(int FaceIndex) {
            if (m_InverseFaceTransformation == null) {
                m_InverseFaceTransformation = new AffineTrafo[this.NoOfFaces];

                int D = this.SpatialDimension;

                for (int l = 0; l < this.NoOfFaces; l++) {
                    var Et = this.GetInverseEmbeddedFaceTrafo(l);

                    m_InverseFaceTransformation[l] = new AffineTrafo(D, D - 1);
                    var tr = m_InverseFaceTransformation[l];

                    for (int i = 0; i < D - 1; i++) {
                        for (int j = 0; j < D; j++) {
                            tr.Matrix[i, j] = Et.Matrix[i, j];
                        }
                        tr.Affine[i] = Et.Affine[i];
                    }
                }
            }
            return m_InverseFaceTransformation[FaceIndex];
        }

        /// <summary>
        /// <see cref="GetFaceTrafo"/>
        /// </summary>
        private AffineTrafo[] m_FaceTransformation;

        /// <summary>
        /// returns an affine-linear transformation which transforms from the
        /// local coordinate system of the <paramref name="FaceIndex"/>-th face
        /// the the local coordinate system of the simplex.
        /// </summary>
        virtual public AffineTrafo GetFaceTrafo(int FaceIndex) {
            if (m_FaceTransformation == null) {
                m_FaceTransformation = new AffineTrafo[this.NoOfFaces];

                int D = this.SpatialDimension;
                int DMinusOne = Math.Max(D - 1, 1);

                for (int l = 0; l < this.NoOfFaces; l++) {
                    var Et = this.GetEmbeddedFaceTrafo(l);

                    m_FaceTransformation[l] = new AffineTrafo(DMinusOne, D);
                    var tr = m_FaceTransformation[l];

                    for (int i = 0; i < D; i++) {
                        for (int j = 0; j < DMinusOne; j++) {
                            tr.Matrix[i, j] = Et.Matrix[i, j];
                        }
                        tr.Affine[i] = Et.Affine[i];
                    }
                }
            }
            return m_FaceTransformation[FaceIndex];
        }

        /// <summary>
        /// See <see cref="FaceTrafoGramianSqrt"/>
        /// </summary>
        double[] m_FaceTrafoGramianSqrt;

        /// <summary>
        /// Square-root of the Gramian determinant of the face-to-volume
        /// transformation, see <see cref="GetFaceTrafo"/>. <br/>
        /// index: face index
        /// </summary>
        public virtual double[] FaceTrafoGramianSqrt {
            get {
                if (m_FaceTrafoGramianSqrt == null) {
                    m_FaceTrafoGramianSqrt = new double[this.NoOfFaces];
                    int D = this.SpatialDimension;

                    MultidimensionalArray Transpose = MultidimensionalArray.Create(D - 1, D);
                    MultidimensionalArray Product = MultidimensionalArray.Create(D - 1, D - 1);

                    for (int l = 0; l < this.NoOfFaces; l++) {
                        var Trafo = GetFaceTrafo(l);
                        Debug.Assert(Trafo.Matrix.NoOfRows == D);
                        Debug.Assert(Trafo.Matrix.NoOfCols == D - 1);

                        Trafo.Matrix.TransposeTo(Transpose);
                        Product.Clear();
                        Product.GEMM(1.0, Transpose, Trafo.Matrix, 0.0);
                        m_FaceTrafoGramianSqrt[l] = Math.Sqrt(Product.Determinant());
                    }
                }

                return m_FaceTrafoGramianSqrt;
            }
        }

        /// <summary>
        /// See <see cref="VertexIndicesToFaces"/>
        /// </summary>
        int[][] m_VertexIndicesToFaces;

        /// <summary>
        /// For some vertex (index) of this reference element, the faces (by
        /// their index) which share this vertex;<br/>
        /// 1st index: vertex index;<br/>
        /// 2nd index: collection, no special order<br/>
        /// </summary>
        public int[][] VertexIndicesToFaces {
            get {
                if (m_VertexIndicesToFaces == null) {
                    int NF = this.NoOfFaces;

                    m_VertexIndicesToFaces = new int[NoOfVertices][];

                    var V2F = this.FaceToVertexIndices;
                    for (int iFace = 0; iFace < this.NoOfFaces; iFace++) {
                        for (int k = 0; k < V2F.GetLength(1); k++) {
                            int I = V2F[iFace, k];

                            if (m_VertexIndicesToFaces[iFace] == null)
                                m_VertexIndicesToFaces[iFace] = new int[0];

                            iFace.AddToArray(ref m_VertexIndicesToFaces[I]);
                        }
                    }

                }
                return m_VertexIndicesToFaces;
            }
        }

        /// <summary>
        /// See <see cref="VertexIndicesToCoFaces"/>
        /// </summary>
        int[][] m_VertexIndicesToCoFaces;

        /// <summary>
        /// For some vertex (index) of this reference element, the co-faces
        /// (by their index) which share this vertex;<br/>
        /// 1st index: vertex index;<br/>
        /// 2nd index: collection, no special order<br/>
        /// </summary>
        public int[][] VertexIndicesToCoFaces {
            get {
                if (m_VertexIndicesToCoFaces == null) {
                    int NF = this.NoOfFaces;

                    m_VertexIndicesToCoFaces = new int[NoOfVertices][];

                    var V2F = this.CoFaceVerticeIndices;
                    for (int iCoFace = 0; iCoFace < this.NoOfFaces; iCoFace++) {
                        for (int k = 0; k < V2F.GetLength(1); k++) {
                            int I = V2F[iCoFace, k];

                            if (m_VertexIndicesToCoFaces[iCoFace] == null)
                                m_VertexIndicesToCoFaces[iCoFace] = new int[0];

                            iCoFace.AddToArray(ref m_VertexIndicesToCoFaces[I]);
                        }
                    }

                }
                return m_VertexIndicesToCoFaces;
            }
        }

        /// <summary>
        /// <see cref="FaceToVertexIndices"/>
        /// </summary>
        int[,] m_FaceToVertexIndices;

        /// <summary>
        /// For all faces, the vertex indices of the reference element that make up the face;<br/>
        /// <list type="bullet">
        ///     <item>
        ///         1st index: face index, in the range of 0 to <see cref="NoOfFaces"/>.
        ///     </item>
        ///     <item>
        ///         2nd index: face vertex index, in the range of 0 to <see cref="FaceRefElement"/>.
        ///     </item>
        /// </list>
        /// </summary>
        public int[,] FaceToVertexIndices {
            get {
                if (m_FaceToVertexIndices == null)
                    InitFaceVerticesIndex();
                return (int[,])m_FaceToVertexIndices.Clone();
            }
        }

        private void InitFaceVerticesIndex() {
            m_FaceToVertexIndices = new int[this.NoOfFaces, this.FaceRefElement.NoOfVertices];

            for (int FaceIndex = 0; FaceIndex < this.NoOfFaces; FaceIndex++) {

                MultidimensionalArray EdgeVertices = FaceRefElement.Vertices.CloneAs();
                MultidimensionalArray VolumeVertices = MultidimensionalArray.Create(EdgeVertices.GetLength(0), this.SpatialDimension);
                TransformFaceCoordinates(FaceIndex, EdgeVertices, VolumeVertices);

                int D = this.SpatialDimension;
                int f = 0;
                for (int i = 0; i < EdgeVertices.GetLength(0); i++) {

                    for(int d = 0; d < D; d++) {
                        
                        //Console.WriteLine("    Face {0}: {1} \t{2}", FaceIndex, d, VolumeVertices[i, d]);
                    }

                    for (int j = 0; j < this.Vertices.GetLength(0); j++) {

                        double dist = 0;

                        for (int d = 0; d < D; d++) {
                            double delta = VolumeVertices[i, d] - Vertices[j, d];
                            dist += delta * delta;

                            //Console.WriteLine("    {0}: {1} \t{2}", d, VolumeVertices[i, d],Vertices[j, d]);
                        }

                        //Console.Write("{3}: {0}-{1} dist = {2}", i, j, dist, FaceIndex);

                        if(dist < 1.0e-6) {
                            m_FaceToVertexIndices[FaceIndex, f] = j;
                            f++;
                            //Console.WriteLine(" BINGO! f = " + f);
                        } else {
                            //Console.WriteLine(" leider nein!"); 
                        }


                    }
                }

                if(f != this.FaceRefElement.NoOfVertices) {
                    //Console.WriteLine("f = " + f);
                    throw new ApplicationException("internal error: unable to find indices of all edge vertices;");
                }

            }

        }

        /// <summary>
        /// Geometric center of all faces; <br/>
        /// 1st index: face index, in the range of 0 to <see cref="NoOfFaces"/>; <br/>
        /// 2nd index: spatial dimension, in the range of 0 (including) to
        /// <see cref="SpatialDimension"/> (excluding)
        /// </summary>
        public NodeSet FaceCenters {
            get {
                if (m_FaceCenters == null)
                    InitFaceCenters();
                return m_FaceCenters;
            }
        }

        /// <summary>
        /// Center of each face.
        /// </summary>
        public NodeSet GetFaceCenter(int iFace) {
            if(iFace < 0 || iFace >= this.NoOfFaces)
                throw new ArgumentException("face index out of range.");

            if(m_FaceCenterPF == null) {
                m_FaceCenterPF = new NodeSet[this.NoOfFaces];
                for(int iF = 0; iF < m_FaceCenterPF.Length; iF++) {
                    m_FaceCenterPF[iF] = new NodeSet(this, this.FaceCenters.GetRow(iF));
                }
            }
            return m_FaceCenterPF[iFace];
        }

        private NodeSet[] m_FaceCenterPF;

        
        /// <summary>
        /// See <see cref="FaceCenters"/>
        /// </summary>
        private NodeSet m_FaceCenters;

        /// <summary>
        /// Normal vectors on all faces;<br/>
        /// 1st index: face index, in the range of 0 to <see cref="NoOfFaces"/>;<br/>
        /// 2nd index: spatial dimension, in the range of 0 (including) to
        /// <see cref="SpatialDimension"/> (excluding)
        /// </summary>
        public MultidimensionalArray FaceNormals {
            get {
                if (m_FaceNormals == null)
                    InitFaceNormals();
                return m_FaceNormals;
            }
        }

        /// <summary>
        /// See <see cref="FaceNormals"/>
        /// </summary>
        private MultidimensionalArray m_FaceNormals;

        private void InitFaceCenters() {
            int[,] ev = FaceToVertexIndices;
            int D = SpatialDimension;
            int N = ev.GetLength(1);

            m_FaceCenters = new NodeSet(this, ev.GetLength(0), D);

            // loop over all edges ...
            for (int m = 0; m < m_FaceCenters.GetLength(0); m++) {

                // loop over the vertices of all edges ...
                for (int n = 0; n < N; n++) {

                    for (int d = 0; d < D; d++)
                        m_FaceCenters[m, d] += Vertices[ev[m, n], d];

                }

                for (int d = 0; d < D; d++)
                    m_FaceCenters[m, d] *= 1.0 / ((double)N);
            }

            m_FaceCenters.LockForever();
        }

        private void InitFaceNormals() {
            var ec = FaceCenters;
            int D = SpatialDimension;

            // RefElements are zero-centered -> vectors to edge centers are
            // already normal to the faces
            m_FaceNormals = MultidimensionalArray.Create(ec.Lengths);
            m_FaceNormals.Set(ec);

            // loop over all edges ...
            for (int m = 0; m < m_FaceCenters.GetLength(0); m++) {

                double len = 0;
                for (int d = 0; d < D; d++) {
                    len += m_FaceNormals[m, d] * m_FaceNormals[m, d];
                }
                len = 1.0 / Math.Sqrt(len);

                for (int d = 0; d < D; d++) {
                    m_FaceNormals[m, d] *= len;
                }
            }
        }

        /// <summary>
        /// The dimension of the reference element in the sense of measure-theory.
        /// </summary>
        public virtual int SpatialDimension {
            get {
                return Vertices.GetLength(1);
            }
        }

        /// <summary>
        /// The vertices of the element in local coordinates;<br/>
        /// Geometrically, the element is defined to be the convex hull of
        /// these points;
        /// </summary>
        /// <remarks>
        /// Indices are defined as follows:
        /// <list type="bullet">
        ///   <item>
        ///     1st index: Vertex index
        ///   </item>
        ///   <item>
        ///     2nd index: spatial dimension, 0 for 1D and 0,1 for 2D and 0,1,2
        ///     for 3D
        ///   </item>
        /// </list>
        /// The average (center of gravity) of these points must be 0. The
        /// number of vertices must be at least <i>D</i>+1, where <i>D</i>
        /// denotes the spatial dimension. Furthermore, if the vertices are
        /// denoted as <i>p</i>[0], <i>p</i>[1], ... <i>p</i>[<i>D</i>],
        /// the vectors
        /// (<i>p</i>[1] - <i>p</i>[0], ... , <i>p</i>[<i>D</i>] - <i>p</i>[0])
        /// must form a basis of the <i>D</i>-dimensional space.
        /// </remarks>
        public NodeSet Vertices {
            get {
                m_Vertices.LockForever();
                return m_Vertices;
            }
        }

        /// <summary>
        /// See <see cref="Vertices"/>
        /// </summary>
        protected NodeSet m_Vertices;


        NodeSet[] m_FaceVertices;

        /// <summary>
        /// Returns the vertices of the <paramref name="iFace"/>--th face.
        /// </summary>
        public NodeSet GetFaceVertices(int iFace) {
            if(m_FaceVertices == null) {
                m_FaceVertices = new NodeSet[this.NoOfFaces];
                int D = this.SpatialDimension;
                int K = this.FaceRefElement.NoOfVertices;
                for(int iF = 0; iF < m_FaceVertices.Length; iF++) {
                    m_FaceVertices[iF] = new NodeSet(this, K, D);
                    this.TransformFaceCoordinates(iF, this.FaceRefElement.Vertices, m_FaceVertices[iF]);
                    this.m_FaceVertices[iF].LockForever();
                }
            }
            return m_FaceVertices[iFace];
        }



        /// <summary>
        /// tests whether a point <paramref name="pt"/>, in local coordinates,
        /// lies within this simplex, i.e. is element of the convex hull of the
        /// <see cref="Vertices"/>
        /// </summary>
        /// <param name="pt"></param>
        /// <returns></returns>
        public bool IsWithin(double[] pt) {
            return IsWithin(pt, 0.0);
        }

        /// <summary>
        /// Tests whether a point <paramref name="pt"/>, in local coordinates,
        /// lies within this simplex (with a tolerance given by
        /// <paramref name="tolerance"/>), i.e. is element of the convex hull
        /// of the <see cref="Vertices"/>
        /// </summary>
        /// <param name="pt"></param>
        /// <param name="tolerance">
        /// A (small) positive number giving the tolerance (again, in the local
        /// coordinate system) of the check.
        /// </param>
        /// <returns></returns>
        public abstract bool IsWithin(double[] pt, double tolerance);

        /// <summary>
        /// Computes, for point <paramref name="pt_in"/>, the closest point on
        /// the boundary of the reference element (<paramref name="pt_out"/>).
        /// </summary>
        /// <param name="pt_in">
        /// input; some point in local cell coordinates
        /// </param>
        /// <param name="pt_out">
        /// output; closest point on reference element boundary, in local
        /// coordinates.
        /// </param>
        /// <returns>the distance</returns>
        public double ClosestPoint(double[] pt_in, double[] pt_out) {
            if (pt_in.Length != this.SpatialDimension)
                throw new ArgumentException("wrong spatial dimension.");
            if (pt_out.Length != this.SpatialDimension)
                throw new ArgumentException("wrong spatial dimension.");

            int D = this.SpatialDimension;
            double distance = double.MaxValue;

            // case 1: check all vertices
            // ==========================
            {
                var Vtx = this.Vertices;
                for (int iVtx = Vtx.GetLength(0) - 1; iVtx >= 0; iVtx--) {
                    double dist = 0;
                    for (int d = 0; d < D; d++)
                        dist += (pt_in[d] - Vtx[iVtx, d]).Pow2();
                    dist = dist.Sqrt();

                    if (dist < distance) {
                        for (int d = 0; d < D; d++) {
                            pt_out[d] = Vtx[iVtx, d];
                        }

                        distance = dist;
                    }
                }

            }

            // case 2: some point on the faces may be closer
            // =============================================
            {
                //MultidimensionalArray _pt_in = MultidimensionalArray.CreateWrapper(pt_in,1,D);
                double[] W = new double[D];
                double[] R = new double[D - 1];
                MultidimensionalArray _W = MultidimensionalArray.CreateWrapper(W, 1, D);
                MultidimensionalArray _R = MultidimensionalArray.CreateWrapper(R, 1, D - 1);
                var edgSplx = this.FaceRefElement;


                int _NoOfFaces = this.NoOfFaces;
                for (int iFace = 0; iFace < _NoOfFaces; iFace++) {

                    // project the oint onto the face:
                    double alfa = 0;
                    for (int d = 0; d < D; d++) {
                        alfa += (pt_in[d] - this.FaceCenters[iFace, d]) * this.FaceNormals[iFace, d];
                    }

                    if (Math.Abs(alfa) < distance) {
                        for (int d = 0; d < D; d++) {
                            W[d] = pt_in[d] - alfa * this.FaceNormals[iFace, d];
                        }

                        this.GetInverseFaceTrafo(iFace).Transform(_W, _R);

                        if (edgSplx.IsWithin(R, 1.0e-8)) {
                            distance = Math.Abs(alfa);
                            Array.Copy(W, pt_out, D);
                        }
                    }
                }
            }

            return distance;
        }

        /// <summary>
        /// provides a quadrature rule of desired or higher order;
        /// </summary>
        /// <param name="DesiredOrder">
        /// The desired polynomial degree that should integrated exactly
        /// by the given quadrature rule;
        /// highest valid value is provided
        /// by <see cref="HighestKnownOrder"/>;
        /// </param>
        /// <returns>
        /// the real order of the returned quadrature rule
        /// <see cref="BoSSS.Foundation.Quadrature.QuadRule.OrderOfPrecision"/>
        /// may be greater or equal to <paramref name="DesiredOrder"/>;
        /// </returns>
        public virtual QuadRule GetQuadratureRule(int DesiredOrder) {
            if (DesiredOrder > HighestKnownOrder) {
                throw new ArgumentOutOfRangeException("no quadrature rule for desired order " + DesiredOrder + " available for simplex " + this.GetType().Name + ".", "DesiredOrder");
            }

#if DEBUG
            foreach(var qr in this.m_QuadRules)
                Debug.Assert(qr.Nodes.IsLocked);
#endif

            Quadrature.QuadRule best = null;
            int OrderOfPrecision = -1;
            foreach (Quadrature.QuadRule qr in m_QuadRules) {
                if (qr.OrderOfPrecision >= DesiredOrder) {
                    if(OrderOfPrecision < 0) {
                        OrderOfPrecision = qr.OrderOfPrecision;
                        best = qr;
                    } else {
                        if(qr.NoOfNodes < best.Weights.Length) {
                            best = qr;
                            OrderOfPrecision = qr.OrderOfPrecision;
                        }
                    }
                }
            }

            if (OrderOfPrecision < 0)
                throw new ArgumentOutOfRangeException("no quadrature rule for desired order available.", "DesiredOrder");

            return best;
        }

        /// <summary>
        /// creates a rule for integrating over the boundary of this reference element
        /// </summary>
        public virtual CellBoundaryQuadRule GetBoundaryQuadRule(int DesiredOrder) {
            
            int NF = this.NoOfFaces;
            var FaceRule = this.FaceRefElement.GetQuadratureRule(DesiredOrder);
            int D = this.SpatialDimension;

            var R = CellBoundaryQuadRule.CreateEmpty(this, FaceRule.NoOfNodes*NF, D, NF);

            int ncnt = 0;
            for (int iface = 0; iface < NoOfFaces; iface++) {
                var Nodes = R.Nodes.ExtractSubArrayShallow(new int[] {ncnt, 0 }, new int[] { ncnt + FaceRule.NoOfNodes - 1, D - 1});
                this.TransformFaceCoordinates(iface, FaceRule.Nodes, Nodes);

                var Weights = R.Weights.ExtractSubArrayShallow(new int[] { ncnt }, new int[] { ncnt + FaceRule.NoOfNodes - 1 });
                Weights.Acc(this.FaceTrafoGramianSqrt[iface], FaceRule.Weights);

                R.NumbersOfNodesPerFace[iface] = FaceRule.NoOfNodes;

                ncnt += FaceRule.NoOfNodes;
            }
            R.Nodes.LockForever();

            return R;
        }

        /// <summary>
        /// Provides a brute-force quadrature rule (suitable for integrating 
        /// </summary>
        /// <param name="NoOfSubDiv">
        /// depth of the subdivision-tree which is used to construct the 
        /// brute-force - rule (see <see cref="GetSubdivisionTree"/>);
        /// </param>
        /// <param name="BaseRuleOrder">
        /// quad rule order of the base rule which is "multiplied"
        /// </param>
        /// <returns>
        /// </returns>
        public virtual QuadRule GetBruteForceQuadRule(int NoOfSubDiv, int BaseRuleOrder) {

            // generate subdivision of this simplex
            // ------------------------------------
            RefElement.SubdivisionTreeNode subdiv = this.GetSubdivisionTree(NoOfSubDiv);
            RefElement.SubdivisionTreeNode[] leaves = subdiv.GetLeaves();

            // get "base rule", which will be "multiplied"
            // -------------------------------------------
            Quadrature.QuadRule BaseRule = GetQuadratureRule(BaseRuleOrder);
            int D = BaseRule.Nodes.GetLength(1);
            double[] vtx = new double[D];

            // return values memalloc
            // ----------------------
            Quadrature.QuadRule ret = QuadRule.CreateEmpty(
                this, BaseRule.Weights.Length * leaves.Length, D);
            ret.OrderOfPrecision = BaseRule.OrderOfPrecision;

            // generate brute-force rule
            // -------------------------

            // loop over all leaves of the subdivision ...
            int cnt = 0;
            for (int i = 0; i < leaves.Length; i++) {

                double det = leaves[i].Trafo2Root.Matrix.Determinant();

                // loop over all nodes of the base quad. rule
                int N = BaseRule.NoOfNodes;
                leaves[i].Trafo2Root.Transform(
                    BaseRule.Nodes,
                    ret.Nodes.ExtractSubArrayShallow(
                        new int[] { cnt, 0 },
                        new int[] { cnt + N - 1, D - 1 }));
                for (int n = 0; n < N; n++) {
                    ret.Weights[cnt] = BaseRule.Weights[n] * det;
                    cnt++;
                }
            }
            ret.Nodes.LockForever();

            // return
            // ------
            
            return ret;
        }

        /// <summary>
        /// order of the most precise quadrature rule which is available for
        /// this reference element
        /// </summary>
        public virtual int HighestKnownOrder {
            get {
                int highest = 0;
                foreach (var qr in m_QuadRules) {
                    highest = Math.Max(qr.OrderOfPrecision, highest);
                }
                return highest;
            }
        }

        /// <summary>
        /// the maximum distance between two vertices of the simplex
        /// (<see cref="Vertices"/>);
        /// </summary>
        /// <returns></returns>
        public double GetMaxDiameter() {
            double h = 0;

            int NoVtx = Vertices.GetLength(0);
            int D = Vertices.GetLength(1);

            for (int i = 0; i < NoVtx; i++) {
                for (int j = i + 1; j < NoVtx; j++) {
                    double distPow2 = 0;
                    double dxd;
                    for (int d = 0; d < D; d++) {
                        dxd = Vertices[i, d] - Vertices[j, d];
                        distPow2 += dxd * dxd;
                    }

                    double dist = Math.Sqrt(distPow2);
                    h = Math.Max(h, dist);
                }
            }

            return h;
        }

        /// <summary>
        /// collection of all available quadrature rules for this simplex;
        /// </summary>
        public ICollection<QuadRule> m_QuadRules = new List<QuadRule>();

        /// <summary>
        /// Orthonormal polynomials with respect to the reference element, sorted by degree.
        /// </summary>
        protected Polynomial[] OrthonormalPolynomials {
            set;
            get;
        }

        /// <summary>
        /// Returns a complete basis of orthonormal polynomials up to degree <see cref="MaxDeg"/>.
        /// </summary>
        public PolynomialList GetOrthonormalPolynomials(int MaxDeg) {
            if(MaxDeg < 0 || MaxDeg > this.HighestSupportedPolynomialDegree) {
                throw new NotSupportedException("polynomial degree " + MaxDeg + " is not supported for " + this.GetType().Name
                        + " (highest supported degree is " + this.HighestSupportedPolynomialDegree + ").");
            }

#if DEBUG
            for(int i = 1; i < this.OrthonormalPolynomials.Length; i++) {
                if(OrthonormalPolynomials[i - 1].AbsoluteDegree > OrthonormalPolynomials[i].AbsoluteDegree) {
                    throw new ApplicationException("Internal error, polynomials for " + this.GetType().Name +" not sorted as they should.");
                }
            }
#endif
            
            var R = new PolynomialList(this.OrthonormalPolynomials.TakeWhile(p => p.AbsoluteDegree <= MaxDeg));
            Debug.Assert(R.Count == GetNoOfOrthonormalPolynomialsUptoDegree(MaxDeg));

            return R;
        }

        /// <summary>
        /// Number of polynomials with an absolute degree equal to <paramref name="p"/>.
        /// </summary>
        public int GetNoOfOrthonormalPolynomialsForDegree(int p) {
            if(p < 0)
                throw new ArgumentOutOfRangeException();

            int n;
            switch(this.SpatialDimension) {
                case 0:
                return 0;

                case 1:
                n = 1;
                break;

                case 2:
                n = p + 1;
                break;

                case 3:
                n = (p * p + 3 * p + 2) / 2;
                break;

                default:
                throw new NotImplementedException("Unknown spatial dimension");
            }

            Debug.Assert(this.OrthonormalPolynomials.Where(poly => poly.AbsoluteDegree == p).Count() == n);

            return n;
        }


        /// <summary>
        /// Number of polynomials with an absolute degree smaller or equal to <paramref name="p"/>.
        /// </summary>
        public int GetNoOfOrthonormalPolynomialsUptoDegree(int p) {
            if(p < 0)
                throw new ArgumentOutOfRangeException();

            int n = 0;
            for(int i = 0; i <= p; i++)
                n += GetNoOfOrthonormalPolynomialsForDegree(i);

            return n;
        }



        /// <summary>
        /// Provides transformations for a subdivision of this simplex.
        /// </summary>
        /// <returns>
        /// An list of transformations, which map the vertices of this simplex
        /// (<see cref="Vertices"/>) to the vertices of the sub-simplexes.
        /// </returns>
        /// <remarks>
        /// By a subdivision (of an simplex), we mean a complete, disjoint
        /// partitioning into a number of <i>n</i> simplexes of similar shape.
        /// <br/>
        /// Note for implementers: The return value of this method should be
        /// equal every time;
        /// </remarks>
        public abstract AffineTrafo[] GetSubdivision();

        /// <summary>
        /// returns the vertices for the first subdivision level;
        /// </summary>
        /// <remarks>
        /// array indices: the 1st index corresponds with the return value of
        /// <see cref="GetSubdivision"/>;<br/>
        /// 2nd and 3rd index are defined like for <see cref="Vertices"/>
        /// </remarks>
        public NodeSet[] SubdivisonVertices {
            get {
                if (m_SubdivisonVertices != null)
                    return m_SubdivisonVertices;

                AffineTrafo[] SubDiv = GetSubdivision();

                m_SubdivisonVertices = new NodeSet[SubDiv.Length];

                int D = this.SpatialDimension;
                double[] vtx = new double[D];
                int[] idx = new int[2];

                for (int i = 0; i < SubDiv.Length; i++) {

                    int J = Vertices.GetLength(0); // No. of Vertice

                    m_SubdivisonVertices[i] = new NodeSet(this, Vertices.GetLength(0), Vertices.GetLength(1));

                    for (int j = 0; j < J; j++) {
                        idx[1] = j;
                        m_Vertices.ExtractVector(vtx, 0, 0, D, idx);
                        double[] Tvtx = SubDiv[i].Transform(vtx);
                        IMatrixExtensions.SetRow<double[]>(m_SubdivisonVertices[i], j, Tvtx);
                    }

                    m_SubdivisonVertices[i].LockForever();
                }

                return m_SubdivisonVertices;
            }
        }

        /// <summary>
        /// cache for <see cref="SubdivisonVertices"/>
        /// </summary>
        NodeSet[] m_SubdivisonVertices;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Levels"></param>
        /// <returns></returns>
        public SubdivisionTreeNode GetSubdivisionTree(int Levels) {
            // check arg.
            if (Levels < 0) {
                throw new ArgumentException(
                    "number of subdivision levels must be greater or equal to 0.", "Levels");
            }

            // gen. level 0
            SubdivisionTreeNode ret = new SubdivisionTreeNode(this);
            // generate subdiv
            ret.SubdivideRecursive(Levels);

            // return
            return ret;
        }

        /// <summary>
        /// One node of a subdivision tree, as provided by <see cref="GetSubdivisionTree(int)"/>
        /// </summary>
        public class SubdivisionTreeNode {

            /// <summary>
            /// Creates the root object of the tree
            /// </summary>
            /// <param name="s"></param>
            internal SubdivisionTreeNode(RefElement s) {
                m_Parrent = null;
                Children = null;
                m_RefElement = s;

                Vertices = s.Vertices;
                InitializeTrafo();
            }

            void BuildExpensiveVertexLists() {
                // vertex list
                double len = m_RefElement.GetMaxDiameter() * 1.0e-5;
                List<double[]> _gvl = new List<double[]>();
                CollectGlobalVerticesRecursive(_gvl, len);

                double[,] gvl = new double[_gvl.Count, this.Vertices.GetLength(1)];
                for (int i = 0; i < gvl.GetLength(0); i++)
                    for (int d = 0; d < gvl.GetLength(1); d++)
                        gvl[i, d] = _gvl[i][d];

                SetGlobalVerticeRecursive(gvl);
            }

            /// <summary>
            /// recursive ctor.
            /// </summary>
            private SubdivisionTreeNode(SubdivisionTreeNode parrent, AffineTrafo trafo) {
                m_Parrent = parrent;
                m_RefElement = parrent.m_RefElement;

                double[,] rvtx = m_RefElement.Vertices.To2DArray();

                int J = rvtx.GetLength(0); // No. of vertices
                int D = rvtx.GetLength(1); // spatial dim.

                this.Vertices = new NodeSet(m_RefElement, J, D);

                for (int j = 0; j < J; j++) {
                    double[] vtx = ArrayTools.GetRow(rvtx, j);
                    double[] t1 = trafo.Transform(vtx);
                    double[] t = m_Parrent.Trafo2Root.Transform(t1);

                    this.Vertices.SetRow(j, t);
                }
                this.Vertices.LockForever();

                InitializeTrafo();
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="__GlobalVertice"></param>
            /// <param name="EqualDiamPow2">
            /// two vertices are considered equal, if their squared Euclidean distance is 
            /// smaller or equal to this number;
            /// </param>
            void CollectGlobalVerticesRecursive(List<double[]> __GlobalVertice, double EqualDiamPow2) {
                // define global vertices indices
                // =============================

                int NoVtx = this.Vertices.GetLength(0);
                int D = this.Vertices.GetLength(1);

                // for all vertices v of this subdivision part...
                m_GlobalVerticeInd = new int[NoVtx];
                for (int n = 0; n < NoVtx; n++) {
                    int ind = -1;

                    // for all vertices vtx_m which are currently in the global vertices list...
                    for (int m = 0; m < __GlobalVertice.Count; m++) {

                        // .. compute Euclidean distance:
                        double[] vtx_m = __GlobalVertice[m];

                        double distPow2 = 0;
                        for (int d = 0; d < D; d++) {
                            double dxd = Vertices[n, d] - vtx_m[d];
                            distPow2 += dxd * dxd;
                        }

                        if (distPow2 <= EqualDiamPow2) {
                            // vertex is already in globals list
                            ind = m;
                            break;
                        }
                    }

                    // vertex not found -> add to list
                    if (ind < 0) {
                        // found new vertex
                        double[] vtx = new double[D];
                        for (int d = 0; d < D; d++)
                            vtx[d] = Vertices[n, d];

                        ind = __GlobalVertice.Count;
                        __GlobalVertice.Add(vtx);
                    }

                    m_GlobalVerticeInd[n] = ind;
                }

                // do recursion
                // ============
                if (Children != null) {
                    foreach (SubdivisionTreeNode c in Children) {
                        c.CollectGlobalVerticesRecursive(__GlobalVertice, EqualDiamPow2);
                    }
                }
            }

            void SetGlobalVerticeRecursive(double[,] __GlobalVertice) {
                m_GlobalVertice = __GlobalVertice;
                if (Children != null) {
                    foreach (SubdivisionTreeNode c in Children) {
                        c.SetGlobalVerticeRecursive(__GlobalVertice);
                    }
                }
            }

            /// <summary>
            /// implementation of <see cref="GetLeaves"/>;
            /// </summary>
            /// <param name="outp"></param>
            private void GetLeavesRecursive(List<SubdivisionTreeNode> outp) {
                if (Children == null)
                    outp.Add(this);
                else {
                    foreach (var c in Children)
                        c.GetLeavesRecursive(outp);
                }
            }

            /// <summary>
            /// Returns all leaves of the subdivision tree,
            /// i.e. all nodes which have no child's.
            /// </summary>
            /// <returns></returns>
            public SubdivisionTreeNode[] GetLeaves() {
                if (m_Parrent != null)
                    throw new ApplicationException("must be called from the parrent object");
                List<SubdivisionTreeNode> ret = new List<SubdivisionTreeNode>();

                GetLeavesRecursive(ret);

                Debug.Assert(ret.ListEquals(ret[0].GetLevel()));

                return ret.ToArray();
            }

            /// <summary>
            /// The refinement level of this node in the tree.
            /// </summary>
            public int RefinementLevel {
                get {
                    if (m_Parrent == null)
                        return 0;
                    else
                        return m_Parrent.RefinementLevel + 1;
                }
            }

            /// <summary>
            /// Returns all tree nodes on the actual refinement level (see <see cref="RefinementLevel"/>).
            /// </summary>
            public SubdivisionTreeNode[] GetLevel() {
                var root = this.TreeRoot;
                var R = new List<SubdivisionTreeNode>();

                root.GetLevelRecursive(R, this.RefinementLevel);

                return R.ToArray();
            }



            NodeSet[] m_FaceVertices;

            /// <summary>
            /// Returns the vertices of the <paramref name="iFace"/>--th face.
            /// </summary>
            public NodeSet GetFaceVertices(int iFace) {
                if(m_FaceVertices == null) {
                    var Kref = m_RefElement;
                    m_FaceVertices = new NodeSet[Kref.NoOfFaces];
                    int NoOfVtxPerFace = Kref.FaceRefElement.NoOfVertices;
                    int D = Kref.SpatialDimension;
                    Debug.Assert(Kref.FaceToVertexIndices.GetLength(0) == Kref.NoOfFaces);
                    Debug.Assert(Kref.FaceToVertexIndices.GetLength(1) == NoOfVtxPerFace);

                    for(int _iFace = 0; _iFace < m_FaceVertices.Length; _iFace++) {
                        m_FaceVertices[_iFace] = new NodeSet(Kref, NoOfVtxPerFace, D);

                        for(int k = 0; k < NoOfVtxPerFace; k++) {
                            int k1 = Kref.FaceToVertexIndices[_iFace, k];
                            for(int d = 0; d < D; d++) {
                                m_FaceVertices[_iFace][k, d] = Vertices[k1, d];
                            }
                        }
                        m_FaceVertices[_iFace].LockForever();
                    }
                }

                return m_FaceVertices[iFace];
            }



            /// <summary>
            /// Cache for <see cref="GetNeighbor(int)"/>.
            /// </summary>
            int[] m_NeighborIndex;

            /// <summary>
            /// Cache for <see cref="GetNeighbor(int)"/>.
            /// </summary>
            int[] m_NeighborFace;


            /// <summary>
            /// Neighbor ship information on the respective refinement level.
            /// </summary>
            /// <param name="iFace">
            /// A face index of the reference element.
            /// </param>
            /// <returns>
            /// An pair
            /// - item 1: neighbor index, i.e. index into the return value of <see cref="GetLevel"/>.
            /// - item 2: face index of the neighbor 
            /// Negative return values indicate that there is no neighbor at <paramref name="iFace"/>.
            /// </returns>
            public Tuple<int,int> GetNeighbor(int iFace) {
                // Check Args
                // ==========
                var Kref = this.m_RefElement;
                int NoOfFaces = Kref.NoOfFaces;
                if (iFace < 0 || iFace >= Kref.NoOfFaces)
                    throw new ArgumentException();

                // check caches
                // ============
                Debug.Assert((m_NeighborFace == null) == (m_NeighborIndex == null));
                if(m_NeighborIndex == null) {
                    m_NeighborIndex = new int[NoOfFaces];
                    m_NeighborFace = new int[NoOfFaces];
                    ArrayTools.SetAll(m_NeighborFace, int.MinValue);
                }

                // search geometrical neigbor for respective face (if not cached)
                // ==============================================================
                if (m_NeighborFace[iFace] == int.MinValue) {
                    int[,] F2V = Kref.FaceToVertexIndices;
                    Debug.Assert(F2V.GetLength(0) == NoOfFaces);

                    int[] GlbVtxIdx = this.GlobalVerticeInd;

                    //    this.GlobalVerticeInd.Select(idx => )
                    int[] FaceVertice = new int[F2V.GetLength(1)];
                    for (int i = 0; i < FaceVertice.Length; i++)
                        FaceVertice[i] = GlbVtxIdx[F2V[iFace, i]];

                    SubdivisionTreeNode[] Neighbors = GetLevel();
                    int[] NeighFaceVertice = new int[F2V.GetLength(1)];
                    bool bFound = false;
                    Debug.Assert(Neighbors.ContainsExactly(this));

                    for (int iNeigh = 0; iNeigh < Neighbors.Length; iNeigh++) {
                        if (object.ReferenceEquals(Neighbors[iNeigh], this))
                            continue;
                        int[] NeighGlbVtxIdx = Neighbors[iNeigh].GlobalVerticeInd;

                        for(int iNeighFace = 0; iNeighFace < NoOfFaces; iNeighFace++) {
                            for (int i = 0; i < FaceVertice.Length; i++)
                                NeighFaceVertice[i] = NeighGlbVtxIdx[F2V[iNeighFace, i]];

                            if(FaceVertice.SetEquals(NeighFaceVertice)) {
                                m_NeighborIndex[iFace] = iNeigh;
                                m_NeighborFace[iFace] = iNeighFace;
                                bFound = true;
                                break;
                            }

                        }

                        if (bFound)
                            break;
                    }

                    if(!bFound) {
                        // no neighbor for respective face
                        m_NeighborIndex[iFace] = -1;

                        double[] Normal = Kref.FaceNormals.GetRow(iFace);
                        double[] TrfNormal = new double[Normal.Length];
                        Debug.Assert(Normal.Length == Kref.SpatialDimension);
                        this.Trafo2Root.Matrix.gemv(1.0, Normal, 0.0, TrfNormal, true);

                        int OrgFace = -1;
                        for (int iOrgFace = 0; iOrgFace < Kref.NoOfFaces; iOrgFace++) {
                            double alpha = GenericBlas.Angle(TrfNormal, Kref.FaceNormals.GetRow(iOrgFace)).Abs();
                            if(alpha < (1.0*Math.PI/180.0)) {
                                Debug.Assert(OrgFace < 0);
                                OrgFace = iOrgFace;
                            }
                        }
                        Debug.Assert(OrgFace >= 0);

                        m_NeighborFace[iFace] = OrgFace;
                    }

                }

                // return
                // ======
                return (new Tuple<int, int>(m_NeighborIndex[iFace], m_NeighborFace[iFace]));
            }


            /// <summary>
            /// Implementation of <see cref="GetLeaves"/>;
            /// </summary>
            private void GetLevelRecursive(List<SubdivisionTreeNode> outp, int Depth) {
                if (Depth == 0) {
                    outp.Add(this);
                } else {
                    foreach (var c in Children)
                        c.GetLevelRecursive(outp, Depth - 1);
                }
            }

            RefElement m_RefElement;

            /// <summary>
            /// The reference element which created the root object of this subdivision tree
            /// during the invocation of <see cref="RefElement.GetSubdivisionTree"/>;
            /// </summary>
            public RefElement RefElement {
                get {
                    return m_RefElement;
                }
            }



            /// <summary>
            /// <list type="bullet">
            ///   <item>1st index: Vertex index</item>
            ///   <item>2nd index: spatial dimension</item>
            /// </list>
            /// </summary>
            public NodeSet Vertices;

            /// <summary>
            /// <see cref="GlobalVerticeInd"/>
            /// </summary>
            int[] m_GlobalVerticeInd;

            /// <summary>
            /// Indices of the vertices that belong to this node 
            /// in the <see cref="GlobalVertice"/>-array.
            /// </summary>
            public int[] GlobalVerticeInd {
                get {
                    if (m_GlobalVerticeInd == null) {
                        this.GetRoot().BuildExpensiveVertexLists();
                    }
                    return m_GlobalVerticeInd;
                }
            }

            /// <summary>
            /// <see cref="m_GlobalVertice"/>
            /// </summary>
            double[,] m_GlobalVertice;

            NodeSet m_GlobalVerticeNS;

            /// <summary>
            /// Global vertex list of the whole subdivision tree.
            /// <list type="bullet">
            ///   <item>1st index: global Vertex index</item>
            ///   <item>2nd index: spatial dimension</item>
            /// </list>
            /// </summary>
            /// <remarks>
            /// The returned object is reference-equal for all object in the 
            /// subdivision tree.
            /// </remarks>
            public NodeSet GlobalVertice {
                get {
                    if (m_GlobalVertice == null) {
                        this.GetRoot().BuildExpensiveVertexLists();
                        m_GlobalVerticeNS = new NodeSet(this.m_RefElement, m_GlobalVertice);
                    }
                    return m_GlobalVerticeNS;
                }
            }

            /// <summary>
            /// <see cref="Parrent"/>
            /// </summary>
            SubdivisionTreeNode m_Parrent;

            /// <summary>
            /// Parent node of this
            /// </summary>
            public SubdivisionTreeNode Parrent {
                get {
                    return m_Parrent;
                }
            }

            /// <summary>
            /// Root of the tree/ancestor of all <see cref="Parrent"/>s.
            /// </summary>
            public SubdivisionTreeNode TreeRoot {
                get {
                    if (m_Parrent == null)
                        return this;
                    else
                        return m_Parrent.TreeRoot;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public int GetParrentIndex() {
                return Array.IndexOf<SubdivisionTreeNode>(m_Parrent.Children, this);
            }

            /// <summary>
            /// Affine-linear transformation from the local coordinates of this node
            /// to the coordinate system of the parent node.
            /// </summary>
            public AffineTrafo Trafo2Parrent {
                get;
                private set;
            }

            /// <summary>
            /// The inverse of <see cref="Trafo2Parrent"/>
            /// </summary>
            public AffineTrafo TrafoFromParrent {
                get;
                private set;
            }

            /// <summary>
            /// Affine-linear transformation from the local coordinates of this node
            /// to the coordinate system of the root node of the subdivision tree.
            /// </summary>
            public AffineTrafo Trafo2Root {
                get;
                private set;
            }

            /// <summary>
            /// inverse of <see cref="Trafo2Root"/>
            /// </summary>
            public AffineTrafo TrafoFromRoot {
                get;
                private set;
            }

            /// <summary>
            /// initializes <see cref="Trafo2Parrent"/> and <see cref="Trafo2Root"/>;<br/>
            /// Therefore, it is required that <see cref="m_Parrent"/> and <see cref="Vertices"/>
            /// are correctly initialized. <br/>
            /// If <see cref="m_Parrent"/> is null, the transformations are set to be the idenity;
            /// </summary>
            void InitializeTrafo() {
                int D = Vertices.GetLength(1);

                if (m_Parrent == null) {
                    Trafo2Parrent = AffineTrafo.Identity(D);
                    Trafo2Root = AffineTrafo.Identity(D);
                } else {
                    SubdivisionTreeNode root = GetRoot();

                    double[,] this_Basis_plus1 = GetSubMatrix(Vertices, 0, D + 1, 0, D);
                    double[,] root_Basis_plus1 = GetSubMatrix(root.Vertices, 0, D + 1, 0, D);
                    
                    Trafo2Root = AffineTrafo.FromPoints(root_Basis_plus1, this_Basis_plus1);
                    Trafo2Parrent = (m_Parrent.TrafoFromRoot) * Trafo2Root;
                }

                TrafoFromParrent = Trafo2Parrent.Invert();
                TrafoFromRoot = Trafo2Root.Invert();
            }

            public static double[,] GetSubMatrix(NodeSet inp, int RowStart, int NoOfRows, int ColStart, int NoOfCols) {
                double[,] outp = new double[NoOfRows, NoOfCols];
                for (int i = 0; i < NoOfRows; i++) {
                    for (int j = 0; j < NoOfCols; j++) {
                        outp[i, j] = inp[RowStart + i, ColStart + j];
                    }
                }
                return outp;
            }

            /// <summary>
            /// returns the root of this tree of objects
            /// </summary>
            /// <returns></returns>
            public SubdivisionTreeNode GetRoot() {
                if (m_Parrent == null)
                    return this;
                else
                    return m_Parrent.GetRoot();
            }

            /// <summary>
            /// 
            /// </summary>
            internal void SubdivideRecursive(int Levels) {
                if (Levels == 0)
                    return;

                AffineTrafo[] subdivS = m_RefElement.GetSubdivision();
                Children = new SubdivisionTreeNode[subdivS.Length];
                for (int i = 0; i < subdivS.Length; i++) {
                    Children[i] = new SubdivisionTreeNode(this, subdivS[i]);
                    Children[i].SubdivideRecursive(Levels - 1);
                }
            }

            /// <summary>
            /// Can be null, is 
            /// initialized by calling <see cref="SubdivideRecursive"/>
            /// </summary>
            public SubdivisionTreeNode[] Children {
                get;
                private set;
            }

            /// <summary>
            /// Computes a transformation matrix \f$ A \f$, which expresses the orthonormal polynomials in the root of the subdivision tree,
            /// \f$ \underline{\phi}^{\text{root}} \f$, in terms of the orthonormal polynomials in this leaf, \f$ \underline{\phi}^{\text{leaf}} \f$,
            /// i.e.
            /// \f[
            ///   \underline{\phi}^{\text{root}} = \underline{\phi}^{\text{leaf}} A .
            /// \f]
            /// </summary>
            /// <param name="p">Polynomial degree.</param>
            /// <returns>
            /// The transformation matrix, which is upper-triangular.
            /// </returns>
            public MultidimensionalArray GetBasisTrafo(int p) {

                PolynomialList Polys = m_RefElement.GetOrthonormalPolynomials(p);
                var qr = m_RefElement.GetQuadratureRule(p * 2);

                int Np = Polys.Count;
                int D = m_RefElement.SpatialDimension;
                int Nk = qr.NoOfNodes;

                MultidimensionalArray PolyAtNodes = MultidimensionalArray.Create(Nk, Np);
                Polys.Evaluate(qr.Nodes, PolyAtNodes);

                NodeSet NodesAtRoot = new NodeSet(m_RefElement, Nk, qr.Nodes.SpatialDimension);
                Trafo2Root.Transform(qr.Nodes, NodesAtRoot);
                NodesAtRoot.LockForever();

                MultidimensionalArray PolyAtRootNodes = MultidimensionalArray.Create(Nk, Np);
                Polys.Evaluate(NodesAtRoot, PolyAtRootNodes);

                double onbScale = 1.0/Math.Sqrt(TrafoFromRoot.Matrix.Determinant());

                MultidimensionalArray A = MultidimensionalArray.Create(Np, Np);
                A.Multiply(onbScale, PolyAtNodes, PolyAtRootNodes, qr.Weights, 0.0, "mn", "km", "kn", "k");

#if DEBUG
                // Check that the matrix is upper triangular
                double LowerTriNorm = 0.0;
                double UpperTriNorm = 0.0;
                for(int iRow = 0; iRow < Np; iRow++) {
                    for(int j = 0; j < iRow; j++)
                        LowerTriNorm += Math.Abs(A[iRow, j]);
                    for(int j = iRow; j < Np; j++)
                        UpperTriNorm += Math.Abs(A[iRow, j]);
                }
                Debug.Assert(LowerTriNorm < UpperTriNorm * 1.0e-10);
#endif
                for(int iRow = 0; iRow < Np; iRow++) {
                    for(int j = 0; j < iRow; j++)
                        A[iRow, j] = 0.0; // strictly enforce lower tri.
                }

#if DEBUG
                for(int k = 0; k < Nk; k++) { // test for each node
                    for(int n = 0; n < Np; n++) { // test for each polynomial
                        double RootVal = PolyAtRootNodes[k, n];

                        double LeafVal = 0.0;
                        for(int m = 0; m < Np; m++) {
                            LeafVal += PolyAtNodes[k, m] * A[m, n] / onbScale;
                        }
                        Debug.Assert(Math.Abs(LeafVal - RootVal) < 1.0e-8);
                    }
                }
#endif
                return A;
            }


        }

        /*

        /// <summary>
        /// Transformation from reference coordinates to physical coordinates
        /// </summary>
        /// <param name="ReferenceVerticesIn">
        /// Input; vertices in the local coordinate system of an element;
        /// <list type="bullet">
        ///   <item>1st index: vertex index;</item>
        ///   <item>
        ///   2nd index: spatial coordinate index 0 for 1D and 0,1 for 2D and
        ///   0,1,2 for 3D;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="JacobianOut">
        /// Output; the Jacobian of the transformation defined by
        /// <paramref name="TransformationParams"/>
        /// <list type="bullet">
        ///   <item>
        ///     1st index: must be set to <paramref name="jOffset"/>
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="ReferenceVerticesIn"/>;
        ///   </item>
        ///   <item>
        ///     3rd index: spatial coordinate index 0 for 1D and 0,1 for 2D and 0,1,2 for 3D;
        ///   </item>
        /// </list>            
        /// </param>
        /// <param name="jOffset">
        /// 1st index into <paramref name="JacobianOut"/>
        /// </param>
        /// <param name="Type"><see cref="CellType"/></param>
        /// <param name="TransformationParams">
        /// <see cref="Element.TransformationParams"/>
        /// </param>
        public void JacobianOfTransformation(MultidimensionalArray ReferenceVerticesIn, MultidimensionalArray JacobianOut, int jOffset, CellType Type, MultidimensionalArray TransformationParams) {
            int K = ReferenceVerticesIn.GetLength(0);
            int D = this.SpatialDimension;
            if (ReferenceVerticesIn.Dimension != 2)
                throw new ArgumentOutOfRangeException("ReferenceVerticesIn", "Wrong array dimension.");
            if (ReferenceVerticesIn.GetLength(1) != D)
                throw new ArgumentOutOfRangeException("ReferenceVerticesIn", "Wrong spatial dimension.");
            if (JacobianOut.Dimension != 4)
                throw new ArgumentOutOfRangeException("JacobianOut", "Wrong array dimension.");
            if (JacobianOut.GetLength(1) != K)
                throw new ArgumentOutOfRangeException("JacobianOut", "Wrong number of nodes.");
            if (JacobianOut.GetLength(2) != D)
                throw new ArgumentOutOfRangeException("JacobianOut", "Wrong spatial dimension.");
            if (JacobianOut.GetLength(3) != D)
                throw new ArgumentOutOfRangeException("JacobianOut", "Wrong spatial dimension.");

            var Polys = this.GetInterpolationPolynomials(Type);

            var PolyGradientValues = MultidimensionalArray.Create(Polys.Length, K, D);
            for (int i = 0; i < Polys.Length; i++) {
                if (Polys == null)
                    Console.Error.WriteLine("polys null");
                if (Polys[i] == null)
                    Console.Error.WriteLine("polys[i] null");
                if (PolyGradientValues == null)
                    Console.Error.WriteLine("PolyGradientValues null");
                Polys[i].EvaluateGradient(PolyGradientValues.ExtractSubArrayShallow(i, -1, -1), ReferenceVerticesIn);
            }
            if (TransformationParams.GetLength(0) != Polys.Length)
                throw new ArgumentOutOfRangeException("TransformationParams");
            if (TransformationParams.GetLength(1) != D)
                throw new ArgumentOutOfRangeException("TransformationParams");
            var S = JacobianOut.ExtractSubArrayShallow(
                new int[] { jOffset, 0, 0, 0 },
                new int[] { jOffset - 1, K - 1, D - 1, D - 1 });
            S.Multiply(1.0, PolyGradientValues, TransformationParams, 0.0, "krc", "ikc", "ir");

#if DEBUG
            for(int k = 0; k < S.GetLength(0); k++) {
                double JacobiDet = S.ExtractSubArrayShallow(k, -1, -1).Determinat();
                Debug.Assert(JacobiDet > 0);
            }
#endif
        }

        /// <summary>
        /// Jacobian of the transformation from reference coordinates to
        /// physical coordinates;
        /// </summary>
        /// <param name="ReferenceVerticesIn"></param>
        /// <param name="GlobalVerticesOut">
        /// <list type="bullet">
        ///   <item>
        ///     1st index: must be set to <paramref name="jOffset"/>
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="ReferenceVerticesIn"/>;
        ///   </item>
        ///   <item>
        ///     3rd index: row index of the Jacobian matrix
        ///   </item>
        ///   <item>
        ///     4th index: column index of the Jacobian matrix
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="jOffset"></param>
        /// <param name="Type"></param>
        /// <param name="TransformationParams"></param>
        public void TransformLocal2Global(MultidimensionalArray ReferenceVerticesIn, MultidimensionalArray GlobalVerticesOut, int jOffset, CellType Type, MultidimensionalArray TransformationParams) {
            int K = ReferenceVerticesIn.GetLength(0);
            int D = this.SpatialDimension;
            if (ReferenceVerticesIn.Dimension != 2)
                throw new ArgumentOutOfRangeException("ReferenceVerticesIn", "Wrong array dimension.");
            if (ReferenceVerticesIn.GetLength(1) != D)
                throw new ArgumentOutOfRangeException("ReferenceVerticesIn", "Wrong spatial dimension.");
            if (GlobalVerticesOut.Dimension != 3)
                throw new ArgumentOutOfRangeException("GlobalVerticesOut", "Wrong array dimension.");
            if (GlobalVerticesOut.GetLength(1) != K)
                throw new ArgumentOutOfRangeException("GlobalVerticesOut", "Wrong number of nodes.");
            if (GlobalVerticesOut.GetLength(2) != D)
                throw new ArgumentOutOfRangeException("GlobalVerticesOut", "Wrong spatial dimension.");

            var Polys = this.GetInterpolationPolynomials(Type);
            var PolyValues = MultidimensionalArray.Create(Polys.Length, K);
            for (int i = 0; i < Polys.Length; i++) {
                Polys[i].Evaluate(PolyValues.ExtractSubArrayShallow(i, -1), ReferenceVerticesIn);
            }

            if (TransformationParams.GetLength(1) != D)
                throw new ArgumentOutOfRangeException("TransformationParams");
            if (TransformationParams.GetLength(0) != Polys.Length)
                throw new ArgumentOutOfRangeException("TransformationParams");

            var S = GlobalVerticesOut.ExtractSubArrayShallow(new int[] { jOffset, 0, 0 }, new int[] { jOffset - 1, K - 1, D - 1 });
            S.Multiply(1.0, PolyValues, TransformationParams, 0.0, "kd", "ik", "id");
        }

        /// <summary>
        /// Adjugate Jacobian of the transformation from reference coordinates
        /// to physical coordinates;
        /// </summary>
        /// <param name="ReferenceVerticesIn">
        /// <list type="bullet">
        ///   <item>1st index: node index </item>
        ///   <item>2nd index: spatial dimension; </item>
        /// </list>
        /// </param>
        /// <param name="AdjugateJacobianOut">
        /// <list type="bullet">
        ///   <item>
        ///     1st index: must be set to <paramref name="jOffset"/>
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="ReferenceVerticesIn"/>
        ///   </item>
        ///   <item>
        ///     3rd index: row index of the AdjugateJacobian matrix
        ///     </item>
        ///   <item>
        ///     4th index: column index of the Adjugate Jacobian matrix
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="jOffset"></param>
        /// <param name="Type"></param>
        /// <param name="TransformationParams"></param>
        public void AdjugateJacobianOfTransformation(MultidimensionalArray ReferenceVerticesIn, MultidimensionalArray AdjugateJacobianOut, int jOffset, CellType Type, MultidimensionalArray TransformationParams) {
            int D = ReferenceVerticesIn.GetLength(1); // dim
            int L = ReferenceVerticesIn.GetLength(0); // number of nodes

            Debug.Assert(ReferenceVerticesIn.Dimension == 2);


            Debug.Assert(AdjugateJacobianOut.GetLength(2) == D, "wrong spatial dimension of JacobianOut");
            Debug.Assert(AdjugateJacobianOut.GetLength(3) == D, "wrong spatial dimension of JacobianOut");
            Debug.Assert(ReferenceVerticesIn.GetLength(1) == D, "wrong spatial dimension of ReferenceVerticesIn");
            Debug.Assert(ReferenceVerticesIn.GetLength(0) == AdjugateJacobianOut.GetLength(1), "mismatch in number of vertices per cell.");
            Debug.Assert(TransformationParams.GetLength(1) == D, "wrong spatial dimension of TransformationParams");

            if (D > 1)
                JacobianOfTransformation(ReferenceVerticesIn, AdjugateJacobianOut, jOffset, Type, TransformationParams);

            Adjugate(AdjugateJacobianOut, jOffset, D, L);
        }

        private static void Adjugate(MultidimensionalArray AdjugateJacobianOut, int jOffset, int D, int L) {
            if (D == 1) {
                // 1D - adjungate
                for (int l = 0; l < L; l++)
                    AdjugateJacobianOut[jOffset, l, 0, 0] = 1.0;
            } else if (D == 2) {
                // 2D - adjungate

                for (int l = 0; l < L; l++) {
                    double m00 = AdjugateJacobianOut[jOffset, l, 0, 0];
                    double m01 = AdjugateJacobianOut[jOffset, l, 0, 1];
                    double m10 = AdjugateJacobianOut[jOffset, l, 1, 0];
                    double m11 = AdjugateJacobianOut[jOffset, l, 1, 1];


                    AdjugateJacobianOut[jOffset, l, 0, 0] = m11;
                    AdjugateJacobianOut[jOffset, l, 0, 1] = -m01;
                    AdjugateJacobianOut[jOffset, l, 1, 0] = -m10;
                    AdjugateJacobianOut[jOffset, l, 1, 1] = m00;
                }
            } else if (D == 3) {
                // 3D - adjungate

                for (int l = 0; l < L; l++) {
                    double m00 = AdjugateJacobianOut[jOffset, l, 0, 0];
                    double m01 = AdjugateJacobianOut[jOffset, l, 0, 1];
                    double m02 = AdjugateJacobianOut[jOffset, l, 0, 2];
                    double m10 = AdjugateJacobianOut[jOffset, l, 1, 0];
                    double m11 = AdjugateJacobianOut[jOffset, l, 1, 1];
                    double m12 = AdjugateJacobianOut[jOffset, l, 1, 2];
                    double m20 = AdjugateJacobianOut[jOffset, l, 2, 0];
                    double m21 = AdjugateJacobianOut[jOffset, l, 2, 1];
                    double m22 = AdjugateJacobianOut[jOffset, l, 2, 2];


                    AdjugateJacobianOut[jOffset, l, 0, 0] = m11 * m22 - m12 * m21;
                    AdjugateJacobianOut[jOffset, l, 0, 1] = -m01 * m22 + m02 * m21;
                    AdjugateJacobianOut[jOffset, l, 0, 2] = m01 * m12 - m02 * m11;
                    AdjugateJacobianOut[jOffset, l, 1, 0] = -m10 * m22 + m12 * m20;
                    AdjugateJacobianOut[jOffset, l, 1, 1] = m00 * m22 - m02 * m20;
                    AdjugateJacobianOut[jOffset, l, 1, 2] = -m00 * m12 + m02 * m10;
                    AdjugateJacobianOut[jOffset, l, 2, 0] = m10 * m21 - m11 * m20;
                    AdjugateJacobianOut[jOffset, l, 2, 1] = -m00 * m21 + m01 * m20;
                    AdjugateJacobianOut[jOffset, l, 2, 2] = m00 * m11 - m01 * m10;
                }
            } else {
                throw new NotSupportedException("spatial dimension of " + D + " is not supported yet.");
            }
        }

        /// <summary>
        /// Inverse Jacobian of the transformation from reference coordinates
        /// to physical coordinates;
        /// </summary>
        /// <param name="ReferenceVerticesIn">
        /// <list type="bullet">
        ///   <item>1st index: node index </item>
        ///   <item>2nd index: spatial dimension; </item>
        /// </list>
        /// </param>
        /// <param name="InverseJacobianOut">
        /// <list type="bullet">
        ///   <item>1st index: must be set to <paramref name="jOffset"/></item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="ReferenceVerticesIn"/>
        ///   </item>
        ///   <item>3rd index: row index of the inverse Jacobian matrix</item>
        ///   <item>
        ///   4th index: column index of the inverse Jacobian matrix
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="jOffset"></param>
        /// <param name="Type"></param>
        /// <param name="TransformationParams"></param>
        public void InverseJacobianOfTransformation(MultidimensionalArray ReferenceVerticesIn, MultidimensionalArray InverseJacobianOut, int jOffset, CellType Type, MultidimensionalArray TransformationParams) {
            int D = ReferenceVerticesIn.GetLength(1); // dim
            int L = ReferenceVerticesIn.GetLength(0); // number of nodes

            Debug.Assert(ReferenceVerticesIn.Dimension == 2);
            Debug.Assert(InverseJacobianOut.GetLength(2) == D, "wrong spatial dimension of JacobianOut");
            Debug.Assert(InverseJacobianOut.GetLength(3) == D, "wrong spatial dimension of JacobianOut");
            Debug.Assert(ReferenceVerticesIn.GetLength(1) == D, "wrong spatial dimension of ReferenceVerticesIn");
            Debug.Assert(ReferenceVerticesIn.GetLength(0) == InverseJacobianOut.GetLength(1), "mismatch in number of vertices per cell.");

            JacobianOfTransformation(ReferenceVerticesIn, InverseJacobianOut, jOffset, Type, TransformationParams);

            unsafe {

                double* det = stackalloc double[L];

                for (int l = 0; l < L; l++) {
                    det[l] = Determinant(D, InverseJacobianOut, l);
                }

                Adjugate(InverseJacobianOut, jOffset, D, L);

                for (int l = 0; l < L; l++) {
                    InverseJacobianOut.ExtractSubArrayShallow(jOffset, l, -1, -1).Scale(1.0 / det[l]);
                }
            }
        }

        /// <summary>
        /// Determinant of the Jacobian matrix of the transformation from
        /// reference coordinates to physical coordinates;
        /// </summary>
        /// <param name="ReferenceVerticesIn"></param>
        /// <param name="JacobianDetOut">
        /// <list type="bullet">
        ///   <item>
        ///     1st index: must be set to <paramref name="jOffset"/>
        ///   </item>
        ///   <item>
        ///     2nd index: vertex index, corresponds with the 1st index of
        ///     <paramref name="ReferenceVerticesIn"/>
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="jOffset"></param>
        /// <param name="Type"></param>
        /// <param name="TransformationParams"></param>
        virtual public void JacobianDetTransformation(MultidimensionalArray ReferenceVerticesIn, MultidimensionalArray JacobianDetOut, int jOffset, CellType Type, MultidimensionalArray TransformationParams) {
            int D = ReferenceVerticesIn.GetLength(1); // dim
            int L = ReferenceVerticesIn.GetLength(0); // number of nodes

            Debug.Assert(ReferenceVerticesIn.Dimension == 2);
            Debug.Assert(JacobianDetOut.Dimension == 2);
            Debug.Assert(JacobianDetOut.GetLength(0) > jOffset, "insufficient number of elements in JacobianDetOut");
            Debug.Assert(JacobianDetOut.GetLength(1) == L, "wrong number of nodes for JacobianDetOut");
            Debug.Assert(ReferenceVerticesIn.GetLength(1) == D, "wrong spatial dimension of ReferenceVerticesIn");

            var Jacobian = MultidimensionalArray.Create(1, L, D, D);
            this.JacobianOfTransformation(ReferenceVerticesIn, Jacobian, 0, Type, TransformationParams);


            for (int l = 0; l < L; l++) {
                double det;
                det = Determinant(D, Jacobian, l);

                Debug.Assert(det >= 0, String.Format(
                    "Encountered negative determinant in cell {0}. This indicates an invalid mesh description.",
                    l));

                JacobianDetOut[jOffset, l] = det;
            }
        }

        private static double Determinant(int D, MultidimensionalArray Jacobian, int l) {
            double det;
            switch (D) {
                case 1:
                    det = Jacobian[0, l, 0, 0];
                    break;

                case 2:
                    det = (Jacobian[0, l, 0, 0] * Jacobian[0, l, 1, 1] - Jacobian[0, l, 0, 1] * Jacobian[0, l, 1, 0]);
                    break;

                case 3:
                    det = (Jacobian[0, l, 0, 0] * Jacobian[0, l, 1, 1] * Jacobian[0, l, 2, 2] - Jacobian[0, l, 0, 0] * Jacobian[0, l, 1, 2] * Jacobian[0, l, 2, 1]
                                 - Jacobian[0, l, 1, 0] * Jacobian[0, l, 0, 1] * Jacobian[0, l, 2, 2] + Jacobian[0, l, 1, 0] * Jacobian[0, l, 0, 2] * Jacobian[0, l, 2, 1]
                                 + Jacobian[0, l, 2, 0] * Jacobian[0, l, 0, 1] * Jacobian[0, l, 1, 2] - Jacobian[0, l, 2, 0] * Jacobian[0, l, 0, 2] * Jacobian[0, l, 1, 1]);
                    break;

                default:
                    throw new NotSupportedException();
            }
            return det;
        }

        /// <summary>
        /// Debug code: computation of the Jacobian by finite differences, in order to
        /// test the correctness of <see cref="JacobianOfTransformation"/>.
        /// </summary>
        [Conditional("DEBUG")]
        public void FiniteDifferenceJacobian(MultidimensionalArray vtx, MultidimensionalArray JacobianOut, int jOffset, CellType Type, MultidimensionalArray TransformationParams) {
            int D = vtx.GetLength(1); // dim
            int L = vtx.GetLength(0); // number of nodes

            //var Jacobian = MultidimensionalArray.Create(1, L, D, D);
            //this.JacobianOfTransformation(vtx, Jacobian, 0, MinorCellType, TransformationParam);

            var vtx_plus = MultidimensionalArray.Create(vtx.Lengths);
            var vtx_minus = MultidimensionalArray.Create(vtx.Lengths);


            var vtx_glob_plus = MultidimensionalArray.Create(1, L, D);
            var vtx_glob_minus = MultidimensionalArray.Create(1, L, D);

            //var Jacobian_finiteDiff = MultidimensionalArray.Create(Jacobian.Lengths);
            for (int d = 0; d < D; d++) {
                double delta = 1.0e-3;

                vtx_minus.Set(vtx);
                vtx_plus.Set(vtx);
                for (int l = 0; l < L; l++) {
                    vtx_minus[l, d] -= 0.5 * delta;
                    vtx_plus[l, d] += 0.5 * delta;
                }

                vtx_glob_plus.Clear();
                vtx_glob_minus.Clear();
                this.TransformLocal2Global(vtx_plus, vtx_glob_plus, 0, Type, TransformationParams);
                this.TransformLocal2Global(vtx_minus, vtx_glob_minus, 0, Type, TransformationParams);

                var JacobiCol = JacobianOut.ExtractSubArrayShallow(-1, -1, -1, d);
                JacobiCol.Acc(1.0, vtx_glob_plus);
                JacobiCol.Acc(-1.0, vtx_glob_minus);
                JacobiCol.Scale(1.0 / delta);
            }
        }
         */

        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var _obj = obj as RefElement;
            if (_obj == null)
                return false;
            else
                return this.Equals(_obj);
        }

        /// <summary>
        /// base impl.
        /// </summary>
        public override int GetHashCode() {
            return this.GetType().GetHashCode();
        }

        /// <summary>
        /// equality
        /// </summary>
        public bool Equals(RefElement obj) {
            if (object.ReferenceEquals(this, obj))
                return true;

            if (this.GetType() != obj.GetType()) // eigentlich reicht Typ-Gleichheit schon...
                return false;

            // ... aber ein paar Tests haben noch nie geschadet.
            Debug.Assert(this.SpatialDimension == obj.SpatialDimension);
            Debug.Assert(this.NoOfVertices == obj.NoOfVertices);
            Debug.Assert(ArrayTools.ListEquals(this.m_Vertices.Storage, obj.m_Vertices.Storage));

            return true;
        }

        /// <summary>
        /// returns the interpolation nodes for this reference element/the
        /// specific <paramref name="type"/>
        /// </summary>
        abstract protected void GetInterpolationNodes_NonLin(CellType type, out NodeSet InterpolationNodes, out PolynomialList InterpolationPolynomials, out int[] NodeType, out int[] EntityIndex);


        Dictionary<CellType, int[]> m_NodeTypes = new Dictionary<CellType, int[]>();
        Dictionary<CellType, int[]> m_EntityIndices = new Dictionary<CellType, int[]>();
        Dictionary<CellType, NodeSet> m_InterpolationNodes = new Dictionary<CellType, NodeSet>();
        Dictionary<CellType, PolynomialList> m_InterpolationPolynomials = new Dictionary<CellType, PolynomialList>();

        bool ContainsZero(PolynomialList g) {
            if (g == null)
                return false;

            foreach (var gg in g)
                if (gg == null)
                    return false;

            return true;
        }

        private void InitMinorType(CellType type) {
            if (m_NodeTypes.ContainsKey(type)) {
                Debug.Assert(m_EntityIndices.ContainsKey(type));
                Debug.Assert(m_InterpolationNodes.ContainsKey(type));
                Debug.Assert(m_InterpolationPolynomials.ContainsKey(type));
                Debug.Assert(ContainsZero(m_InterpolationPolynomials[type]), "interpolation poly not found");
            } else {
                NodeSet InterpolationNodes;
                PolynomialList InterpolationPolynomials;
                int[] NodeType;
                int[] EntityIndex;

                if (type.IsLinear()) {
                    this.SelectNodalPolynomials(1, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                } else {
                    GetInterpolationNodes_NonLin(type, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                }

                m_NodeTypes.Add(type, NodeType);
                m_EntityIndices.Add(type, EntityIndex);
                m_InterpolationNodes.Add(type, InterpolationNodes);
                m_InterpolationPolynomials.Add(type, InterpolationPolynomials);
                Debug.Assert(ContainsZero(InterpolationPolynomials), "interpolation poly not found, 3");
                Debug.Assert(ContainsZero(m_InterpolationPolynomials[type]), "interpolation poly not found, 2");
            }
        }

        /// <summary>
        /// returns the interpolation nodes for this reference element/the
        /// specific <paramref name="type"/>
        /// </summary>
        public NodeSet GetInterpolationNodes(CellType type) {
            InitMinorType(type);
            return m_InterpolationNodes[type];
        }

        /// <summary>
        /// returns the interpolation nodes types for this reference
        /// element/the specific <paramref name="type"/>
        /// </summary>
        public int[] GetInterpolationNodes_NodeType(CellType type) {
            InitMinorType(type);
            return m_NodeTypes[type];
        }

        /// <summary>
        /// returns the interpolation entity indices for this reference
        /// element/the specific <paramref name="type"/>
        /// </summary>
        public int[] GetInterpolationNodes_EntityIndices(CellType type) {
            InitMinorType(type);
            return m_EntityIndices[type];
        }

        /// <summary>
        /// Returns the interpolation polynomials for this reference
        /// element/the specific <paramref name="type"/>.
        /// </summary>
        public PolynomialList GetInterpolationPolynomials(CellType type) {
            InitMinorType(type);
            return m_InterpolationPolynomials[type];
        }


        Dictionary<CellType,PolynomialList[]> m_InterpolationPolynomials1stDeriv = new Dictionary<CellType, PolynomialList[]>();


        /// <summary>
        /// Returns the first derivatives/the gradient 
        /// of the interpolation polynomials for this reference
        /// element/the specific <paramref name="type"/>.
        /// </summary>
        public PolynomialList[] GetInterpolationPolynomials1stDeriv(CellType type) {
            InitMinorType(type);

            PolynomialList[] R;
            if(!this.m_InterpolationPolynomials1stDeriv.TryGetValue(type, out R)) {

                int D = this.SpatialDimension;

                R = new PolynomialList[D];
                PolynomialList Polys = this.GetInterpolationPolynomials(type);

                int[] DerivExp = new int[D];

                for(int d = 0; d < D; d++) {
                    DerivExp[d] = 1;

                    int N = Polys.Count;
                    Polynomial[] tmp = new Polynomial[N];
                    for(int n = 0; n < N; n++) {
                        tmp[n] = Polys[n].Derive(DerivExp);
                    }

                    R[d] = new PolynomialList(tmp);

                    DerivExp[d] = 0;
                }
                this.m_InterpolationPolynomials1stDeriv.Add(type, R);
            }
            return R;
        }

        Dictionary<CellType, int> m_InterpolationDegree = new Dictionary<CellType, int>();

        /// <summary>
        /// returns the polynomial degree of the local-to-global -- mapping;
        /// </summary>
        public int GetInterpolationDegree(CellType type) {
            InitMinorType(type);
            int maxP;
            if (!m_InterpolationDegree.TryGetValue(type, out maxP)) {
                var P = GetInterpolationPolynomials(type);
                maxP = P.Max(p => p.AbsoluteDegree);
            }
            return maxP;
        }

        private CellType[] m_SupportedCellTypes = null;

        /// <summary>
        /// returns a list of supported minor cell index types
        /// </summary>
        public IEnumerable<CellType> SupportedCellTypes {
            get {
                if (m_SupportedCellTypes == null) {
                    string[] ee = typeof(CellType).GetEnumNames();
                    string thisname = this.GetType().Name.ToLowerInvariant();
                    IEnumerable<string> eee = ee.Where(t => t.ToLowerInvariant().StartsWith(thisname));
                    m_SupportedCellTypes = eee.Select(t =>
                        ((CellType)Enum.Parse(typeof(CellType), t))).ToArray();
                }

                return m_SupportedCellTypes;
            }
        }

        /// <summary>
        /// helper function to switch array elements;
        /// </summary>
        protected void SwitchNode(ref int source, ref int target) {
            int buffer = source;

            source = target;
            target = buffer;
        }

        /// <summary>
        /// List of known meshing tools / file formants
        /// </summary>
        public enum ExchangeFormats {

            /// <summary>
            /// Gambit Neutral file format ('*.neu'-files)
            /// </summary>
            GambitNeutral,

            /// <summary>
            /// Gmsh Mesh files ('*.msh');
            /// </summary>
            Gmsh,

            /// <summary>
            /// CFD General Notation System (CGNS);
            /// </summary>
            CGNS
        }

        /// <summary>
        /// provides a permutation between BoSSS element numbering convention
        /// and other tools
        /// </summary>
        /// <param name="type"></param>
        /// <param name="conv"></param>
        /// <returns>
        /// the mapping from BoSSS to Foreign elements; <br/>
        /// index <em>i</em>: BoSSS vertex index; <br/>
        /// content <em>c</em>[<em>i</em>]: node index in the 'foreign'
        /// specification;<br/>
        /// If the return value is null, the specific element, determined by
        /// the <paramref name="type"/>, is not supported in the foreign
        /// specification;
        /// </returns>
        public virtual int[] GetForeignElementMapping(CellType type, RefElement.ExchangeFormats conv) {
            int[] permutationArray = new int[Vertices.GetLength(0)];
            for (int i = 0; i < permutationArray.Length; i++) {
                permutationArray[i] = i;
            }
            return permutationArray;
        }

        /// <summary>
        /// converts a BoSSS Element specification (given by the object type of
        /// this object <see cref="Object.GetType"/>, the minor type index
        /// <paramref name="type"/>) into a foreign specification.
        /// </summary>
        /// <param name="type">BoSSS cell type specification</param>
        /// <param name="conv">Foreign format convention</param>
        /// <returns>
        /// false, if the element is not supported in the
        /// <paramref name="conv"/>-specification
        /// </returns>
        /// <param name="ForeignName">
        /// If specified in the <paramref name="conv"/>-convention, the element
        /// name.
        /// </param>
        /// <param name="ForeignTypeConstant">
        /// If specified in the <paramref name="conv"/>-convention, some integer
        /// constant which identifies the element type..
        /// </param>
        abstract public void GetForeignElementType(
            CellType type, ExchangeFormats conv, out string ForeignName, out int ForeignTypeConstant);

        /// <summary>
        /// provides a permutation between BoSSS element numbering convention
        /// and other tools
        /// </summary>
        /// <param name="type">
        /// The BoSSS type of the cell
        /// </param>
        /// <param name="ForeignNodes">
        /// Interpolation nodes, in any user defined sequence; however, the
        /// nodes must be the same as the ones that BoSSS is using;
        /// Use an Affine transformation, if necessary (see
        /// <see cref="AffineTrafo"/>).
        /// </param>
        /// <returns>
        /// the mapping from BoSSS to Foreign elements; <br/>
        /// index <em>i</em>: BoSSS node index; <br/>
        /// content <em>c</em>[<em>i</em>]: node index in the 'foreign'
        /// specification;<br/>
        /// If the return value is null, the specific element, determined by
        /// the <paramref name="type"/>, is not supported in the foreign
        /// specification;
        /// </returns>
        public int[] ComputeForeignElementMapping(CellType type, MultidimensionalArray ForeignNodes) {
            if (ForeignNodes.Dimension != 2)
                throw new ArgumentException();
            if (ForeignNodes.GetLength(1) != this.SpatialDimension)
                throw new ArgumentException();

            var BoSSSNodes = this.GetInterpolationNodes(type);

            if (BoSSSNodes.GetLength(0) != ForeignNodes.GetLength(0))
                throw new ArgumentException();

            int K = BoSSSNodes.GetLength(0);
            bool[] Zugeordnet = new bool[K];

            int[] R = new int[K];
            for (int k = 0; k < K; k++) {
                var B_pt = BoSSSNodes.GetRow(k);

                bool found = false;
                for (int l = 0; l < K; l++) {
                    if (Zugeordnet[l])
                        continue;

                    var F_pt = ForeignNodes.GetRow(l);

                    if (GenericBlas.L2Dist(B_pt, F_pt) < 1.0e-9) {
                        found = true;
                        R[k] = l;
                        Zugeordnet[l] = true;
                        break;
                    }
                }

                if (!found)
                    throw new ArgumentException("unable to assign foreign node.");
            }

            return R;
        }
    }
}
