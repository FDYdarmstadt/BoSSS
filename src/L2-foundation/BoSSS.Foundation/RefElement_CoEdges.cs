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
using ilPSP;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using ilPSP.Utils;
using System.Threading;

namespace BoSSS.Foundation.Grid.RefElements {

    public partial class RefElement {

        /// <summary>
        /// helper structure
        /// </summary>
        class CoFaceInfo {
            public int Index;
            public int[] iVtx;
            public int iFace0 = -1;
            public int iFace1 = -1;
        }

        Tuple<int[,], int[,]> InitCoFaces() {
            int D = this.SpatialDimension;


            var CoFaceList = new List<CoFaceInfo>();




            if (D == 3) {

                if (!(this.FaceRefElement.FaceRefElement.NoOfVertices == 2)) {
                    // co-edges are supposed to be lines
                    throw new NotSupportedException("Don't know what to do.");
                }


                var iY = this.FaceToVertexIndices;
                var iZ = this.FaceRefElement.FaceToVertexIndices;

                if (iZ.GetLength(1) != 2)
                    throw new NotSupportedException("Don't know what to do.");
                if (iY.GetLength(1) != iZ.GetLength(0))
                    throw new NotSupportedException("Don't know what to do.");


                for (int iFace = 0; iFace < this.NoOfFaces; iFace++) {
                    for (int iCoFace = 0; iCoFace < this.FaceRefElement.NoOfFaces; iCoFace++) {

                        int iV0 = iY[iFace, iZ[iCoFace, 0]];
                        int iV1 = iY[iFace, iZ[iCoFace, 1]];
                        if (iV0 > iV1) {
                            int tmp = iV0;
                            iV0 = iV1;
                            iV1 = tmp;
                        }


                        CoFaceInfo CoFace;
                        var L = CoFaceList.Where(ce => (ce.iVtx[0] == iV0 && ce.iVtx[1] == iV1));
                        if (L.Count() > 0) {
                            Debug.Assert(L.Count() == 1, "Error in algorithm.");
                            CoFace = L.First();
                            Debug.Assert(CoFace.iFace0 >= 0);
                            Debug.Assert(CoFace.iFace1 < 0);
                            CoFace.iFace1 = iFace;
                        } else {
                            CoFace = new CoFaceInfo();
                            CoFace.Index = CoFaceList.Count();
                            CoFace.iFace0 = iFace;
                            CoFace.iVtx = new int[] { iV0, iV1 };
                            CoFaceList.Add(CoFace);
                        }


                    }
                }


            } else if (D == 2) {
                // 2D: Co-Edges are equivalent to vertices
                // +++++++++++++++++++++++++++++++++++++++

                if (!(this.FaceRefElement.FaceRefElement.NoOfVertices == 1)) {
                    throw new NotSupportedException("don't know what to do.");
                }

                var face2Vtx = this.FaceToVertexIndices;
                if (face2Vtx.GetLength(1) != 2)
                    throw new ApplicationException("Don't know what to do.");
                Debug.Assert(face2Vtx.GetLength(1) == FaceRefElement.NoOfVertices);


                for (int iCoFace = 0; iCoFace < this.NoOfVertices; iCoFace++) { // loop over vertices
                    //                                                             in 2D: vertex == co-face
                    var CoFace = new CoFaceInfo();
                    CoFace.Index = iCoFace;
                    CoFace.iVtx = new int[] { iCoFace };

                    for (int iFace = 0; iFace < this.NoOfFaces; iFace++) {
                        for (int t = 0; t < 2; t++) {
                            if (face2Vtx[iFace, t] == iCoFace) {
                                Debug.Assert(CoFace.iFace1 < 0);
                                if (CoFace.iFace0 < 0) {
                                    CoFace.iFace0 = iFace;
                                } else {
                                    CoFace.iFace1 = iFace;
                                }
                            }
                        }
                    }


                    CoFaceList.Add(CoFace);
                }
            } else if (D == 1) {
                // 1D: No Co-Edges
                // +++++++++++++++

            } else {
                throw new NotSupportedException("unknown spatial dimension");
            }

            int[,] _CoFaceVerticeIndices, _CoFacesToFaceIndices;

            // Assemble final data structures
            // ==============================
            if(D == 0) {
                _CoFaceVerticeIndices = new int[0, 0];
                _CoFacesToFaceIndices = new int[0, 2];
            } else {
                Debug.Assert(FaceRefElement.FaceRefElement.NoOfVertices == D - 1);
                _CoFaceVerticeIndices = new int[CoFaceList.Count, D - 1];
                _CoFacesToFaceIndices = new int[CoFaceList.Count, 2];

                for (int iCoFace = 0; iCoFace < CoFaceList.Count; iCoFace++) {
                    var CF = CoFaceList[iCoFace];
                    Debug.Assert(CF.iVtx.Length == _CoFaceVerticeIndices.GetLength(1));

                    for (int i = 0; i < _CoFaceVerticeIndices.GetLength(1); i++) {
                        _CoFaceVerticeIndices[iCoFace, i] = CF.iVtx[i];
                    }

                    if(CF.iFace0 == CF.iFace1)
                        throw new ApplicationException("error in algorithm");

                    if (CF.iFace0 < CF.iFace1) {
                        _CoFacesToFaceIndices[iCoFace, 0] = CF.iFace0;
                        _CoFacesToFaceIndices[iCoFace, 1] = CF.iFace1;
                    } else {
                        _CoFacesToFaceIndices[iCoFace, 0] = CF.iFace1;
                        _CoFacesToFaceIndices[iCoFace, 1] = CF.iFace0;
                    }
                }
            }


            return new Tuple<int[,], int[,]>(_CoFaceVerticeIndices, _CoFacesToFaceIndices);
        }

        //Lazy<(int[,] CoFaceVerticeIndices, int[,] CoFacesToFaceIndices)> m_CoFace = new Lazy<(int[,] CoFaceVerticeIndices, int[,] CoFacesToFaceIndices)>(InitCoFaces, true);

        Tuple<int[,],int[,]> m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices; // (int[,] CoFaceVerticeIndices, int[,] CoFacesToFaceIndices)

        /// <summary>
        /// Vertex indices of the co-faces;
        /// - 1st index: co-face index <em>iCoFace</em>
        /// - 2nd index: vertex index of the Co-Face
        /// </summary>
        public int[,] CoFaceVerticeIndices {
            get {
                LazyInitializer.EnsureInitialized(ref m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices, InitCoFaces);
                return m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices.Item1.CloneAs();
            }
        }

        /// <summary>
        /// For each co-face, the indices of the faces whose geometric intersection forms the <em>iCoFace</em>-th co-face.
        /// - 1st index: co-face index <em>iCoFace</em> 
        /// - 2nd index: in {0,1}, corresponds to first and second face.
        /// </summary>
        /// <remarks>
        /// Note that each co-face can be described as the geometric intersection of 2 faces.
        /// </remarks>
        public int[,] CoFaceToFaceIndices {
            get {
                LazyInitializer.EnsureInitialized(ref m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices, InitCoFaces);
                return m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices.Item2.CloneAs();
            }
        }

        int[,] m_CoFaceToFaceFaceIndex;

        /// <summary>
        /// - 1st index: face-index
        /// - 2nd index: face-index of the face
        /// </summary>
        public int[,] FaceToCoFaceIndices {
            get {
                LazyInitializer.EnsureInitialized(ref m_FaceToCoFaceIndices, delegate () {

                    int[,] cf2fi = CoFaceToFaceIndices;
                    int[,] cf2ff = CoFaceToFaceFaceIndex;


                    var _FaceToCoFaceIndices = new int[this.NoOfFaces, this.FaceRefElement.NoOfFaces];
                    _FaceToCoFaceIndices.SetAll(-1);
                    for(int iCoFace = 0; iCoFace < this.NoOfCoFaces; iCoFace++) {
                        for(int iInOt = 0; iInOt < 2; iInOt++) {
                            int iFace = cf2fi[iCoFace, iInOt];
                            int iFaceFace = cf2ff[iCoFace, iInOt];

                            _FaceToCoFaceIndices[iFace, iFaceFace] = iCoFace;
                        }
                    }

                    if(_FaceToCoFaceIndices.Reshape(false).Any(i => i < 0))
                        throw new Exception("internal error");

                    return _FaceToCoFaceIndices;
                });
                return m_FaceToCoFaceIndices.CloneAs();
            }
        }

        int[,] m_FaceToCoFaceIndices;   


        /// <summary>
        /// Correlates with <see cref="CoFaceToFaceIndices"/>:
        /// For each co-face, the face index of the face on the face reference element which it corresponds to.
        /// - 1st index: co-face index <em>iCoFace</em> 
        /// - 2nd index: in {0,1}, corresponds to first and second face.
        /// </summary>
        public int[,] CoFaceToFaceFaceIndex {
            get {
                
                double PointSetComparison_OneWay(IMatrix A, IMatrix B) {
                    double acc = 0;
                    for(int i = 0; i < A.NoOfRows; i++) {
                        B.MindistRow(A.GetRow(i), out double Dmin, out _);
                        acc += Dmin * Dmin;
                    }
                    return acc.Sqrt();
                }

                double PointSetComparison_TwoWay(IMatrix A, IMatrix B) {
                    double ABdist = PointSetComparison_OneWay(A, B);
                    double BAdist = PointSetComparison_OneWay(B, A);

                    return Math.Sqrt(ABdist.Pow2() + BAdist.Pow2());
                }


                LazyInitializer.EnsureInitialized(ref m_CoFaceToFaceFaceIndex, delegate () {
                    int[,] _CoFaceToFaceFaceIndex;
                    if(this.SpatialDimension <= 1) {
                        _CoFaceToFaceFaceIndex = new int[0, 2];
                    } else {
                        var cfvtx = CoFaceVerticeIndices;
                        var cvf2f = CoFaceToFaceIndices;
                        int NoOfCoFaces = cfvtx.GetLength(0);
                        _CoFaceToFaceFaceIndex = new int[NoOfCoFaces, 2];

                        var CoFaceRefElement = FaceRefElement.FaceRefElement;
                        Debug.Assert(CoFaceRefElement.NoOfVertices == cfvtx.GetLength(1));


                        for(int iCoFace = 0; iCoFace < NoOfCoFaces; iCoFace++) {

                            var coFaceVeritices_A = MultidimensionalArray.Create(CoFaceRefElement.NoOfVertices, this.SpatialDimension);
                            for(int iVtx = 0; iVtx < CoFaceRefElement.NoOfVertices; iVtx++) {
                                coFaceVeritices_A.SetRowPt(iVtx, this.Vertices.GetRowPt(cfvtx[iCoFace, iVtx]));
                            }


                            for(int iInOt = 0; iInOt < 2; iInOt++) {
                                int iFace = cvf2f[iCoFace, iInOt];

                                int iFaceFaceFound = -1;
                                for(int iFaceFace = 0; iFaceFace < FaceRefElement.NoOfFaces; iFaceFace++) {

                                    var coFaceVeritices_B = this.GetFaceTrafo(iFace).Transform(FaceRefElement.GetFaceVertices(iFaceFace));
                                    var dist = PointSetComparison_TwoWay(coFaceVeritices_A, coFaceVeritices_B);
                                    if(dist < BLAS.MachineEps.Sqrt()) {
                                        if(iFaceFaceFound >= 0)
                                            throw new ApplicationException("error in algorithm (1)");
                                        iFaceFaceFound = iFaceFace;
                                    }
                                }

                                if(iFaceFaceFound < 0)
                                    throw new ApplicationException("error in algorithm (2)");
                                _CoFaceToFaceFaceIndex[iCoFace, iInOt] = iFaceFaceFound;

                            }
                        }

                        // check:
                        for(int iCoFace = 0; iCoFace < NoOfCoFaces; iCoFace++) {
                            int iFace0 = m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices.Item2[iCoFace, 0];
                            int iCoFc0 = _CoFaceToFaceFaceIndex[iCoFace, 0];

                            int iFace1 = m_CoFace__CoFaceVerticeIndices_CoFacesToFaceIndices.Item2[iCoFace, 1];
                            int iCoFc1 = _CoFaceToFaceFaceIndex[iCoFace, 1];

                            var Vtx0 = this.GetFaceTrafo(iFace0).Transform(FaceRefElement.GetFaceTrafo(iCoFc0).Transform(CoFaceRefElement.Vertices));
                            var Vtx1 = this.GetFaceTrafo(iFace0).Transform(FaceRefElement.GetFaceTrafo(iCoFc0).Transform(CoFaceRefElement.Vertices));
                            double err = Vtx0.L2Dist(Vtx1);
                            if(err > BLAS.MachineEps.Sqrt())
                                throw new ApplicationException("error in algorithm (3)");

                            var __coFaceVeritices = MultidimensionalArray.Create(CoFaceRefElement.NoOfVertices, this.SpatialDimension);
                            for(int iVtx = 0; iVtx < CoFaceRefElement.NoOfVertices; iVtx++) {
                                __coFaceVeritices.SetRowPt(iVtx, this.Vertices.GetRowPt(cfvtx[iCoFace, iVtx]));
                            }

                            double err2 = PointSetComparison_TwoWay(__coFaceVeritices, Vtx0);
                            if(err > BLAS.MachineEps.Sqrt())
                                throw new ApplicationException("error in algorithm (4)");

                        }
                    }
                    return _CoFaceToFaceFaceIndex;
                });

                return m_CoFaceToFaceFaceIndex.CloneAs();
            }
        }

        /// <summary>
        /// number of co-faces
        /// </summary>
        public int NoOfCoFaces {
            get {
                return CoFaceToFaceIndices.GetLength(0);
            }
        }

        /// <summary>
        /// Transformation from the local coordinate system of the co-face to the local coordinate system for this reference element.
        /// </summary>
        virtual public AffineTrafo GetCoFaceTrafo(int CoFaceIndex) {
            if(this.SpatialDimension < 2)
                throw new NotSupportedException("Co-Faces exist only for 2D and higher.");

            int iFace = this.CoFaceToFaceIndices[CoFaceIndex, 0];
            var cf = this.GetFaceTrafo(iFace);

            var cfT = this.FaceRefElement.GetFaceTrafo(CoFaceToFaceFaceIndex[CoFaceIndex, 0]);

            var coT = cf * cfT;
            return coT;
        }





        AffineManifold[] m_FacePlanes;

        /// <summary>
        /// For each face, the affine manifold that represents the plane in which the face is located;
        /// </summary>
        /// <param name="iFace">
        /// face index
        /// </param>
        public AffineManifold GetFacePlane(int iFace) {
            LazyInitializer.EnsureInitialized(ref m_FacePlanes, delegate () {
                var _FacePlanes = new AffineManifold[this.NoOfFaces];

                var fCen = this.FaceCenters;
                var fNor = this.FaceNormals;
                for(int iF = 0; iF < this.NoOfFaces; iF++) {
                    _FacePlanes[iF] = new AffineManifold(fNor.GetRow(iF), fCen.GetRow(iF));
                    Debug.Assert(_FacePlanes[iF].Normal.AbsSquare() > 0.0);
                    Debug.Assert(_FacePlanes[iF].PointDistance(fCen.GetRow(iF)).Abs() < 1.0e-6);
                }
                return _FacePlanes;
            });

            Debug.Assert(m_FacePlanes[iFace].Normal.Dim == this.SpatialDimension);
            Debug.Assert(m_FacePlanes[iFace].Normal.AbsSquare() > 0.0);
            return m_FacePlanes[iFace].CloneAs();
        }
    
    }
}
