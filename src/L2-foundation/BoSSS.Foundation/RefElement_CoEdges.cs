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

        void InitCoFaces() {
            int D = this.SpatialDimension;


            var CoEdgesList = new List<CoFaceInfo>();




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


                        CoFaceInfo CoEdge;
                        var L = CoEdgesList.Where(ce => (ce.iVtx[0] == iV0 && ce.iVtx[1] == iV1));
                        if (L.Count() > 0) {
                            Debug.Assert(L.Count() == 1, "Error in algorithm.");
                            CoEdge = L.First();
                            Debug.Assert(CoEdge.iFace0 >= 0);
                            Debug.Assert(CoEdge.iFace1 < 0);
                            CoEdge.iFace1 = iFace;
                        } else {
                            CoEdge = new CoFaceInfo();
                            CoEdge.Index = CoEdgesList.Count();
                            CoEdge.iFace0 = iFace;
                            CoEdge.iVtx = new int[] { iV0, iV1 };
                            CoEdgesList.Add(CoEdge);
                        }


                    }
                }


            } else if (D==2) {
                // 2D: Co-Edges are equivalent to vertices
                // +++++++++++++++++++++++++++++++++++++++


                if (!(this.FaceRefElement.FaceRefElement.NoOfVertices == 1)) {
                    throw new NotSupportedException("don't know what to do.");
                }

                var iY = this.FaceToVertexIndices;
                if (iY.GetLength(1) != 2)
                    throw new ApplicationException("Don't know what to do.");


                for (int iCoFace = 0; iCoFace < this.NoOfVertices; iCoFace++) { // loop over vertices
                    //                                                             in 2D: vertex == co-edge
                    var CoEdge = new CoFaceInfo();
                    CoEdge.Index = iCoFace;
                    CoEdge.iVtx = new int[] { iCoFace };

                    for (int iFace = 0; iFace < this.NoOfFaces; iFace++) {
                        for (int t = 0; t < 2; t++) {
                            if (iY[iFace, 0] == iCoFace) {
                                Debug.Assert(CoEdge.iFace1 < 0);
                                if (CoEdge.iFace0 < 0) {
                                    CoEdge.iFace0 = iFace;
                                } else {
                                    CoEdge.iFace1 = iFace;
                                }
                            }
                        }
                    }
                }
            } else if (D == 1) {
                // 1D: No Co-Edges
                // +++++++++++++++

            } else {
                throw new NotSupportedException("unknown spatial dimension");
            }


            // Assemble final data structures
            // ==============================
            {
                m_CoFaceVerticeIndices = new int[CoEdgesList.Count, D-1];
                m_CoFaces_FaceIndices = new int[CoEdgesList.Count, 2];

                for (int iCoFace = 0; iCoFace < CoEdgesList.Count; iCoFace++) {
                    var CF = CoEdgesList[iCoFace];
                    Debug.Assert(CF.iVtx.Length == m_CoFaceVerticeIndices.GetLength(1));

                    for (int i = 0; i < m_CoFaceVerticeIndices.GetLength(1); i++) {
                        m_CoFaceVerticeIndices[iCoFace, i] = CF.iVtx[i];
                    }

                    if (CF.iFace0 < CF.iFace1) {
                        m_CoFaces_FaceIndices[iCoFace, 0] = CF.iFace0;
                        m_CoFaces_FaceIndices[iCoFace, 1] = CF.iFace1;
                    } else {
                        m_CoFaces_FaceIndices[iCoFace, 0] = CF.iFace1;
                        m_CoFaces_FaceIndices[iCoFace, 1] = CF.iFace0;
                    }
                }
            }
        }


        int[,] m_CoFaceVerticeIndices;

        int[,] m_CoFaces_FaceIndices;

        /// <summary>
        /// Vertex indices of the co-faces; <br/>
        /// 1st index: co-face index <em>iCoFace</em>  <br/>
        /// 2nd index: vertex index of the Co-Face
        /// </summary>
        public int[,] CoFaceVerticeIndices {
            get {
                if (m_CoFaceVerticeIndices == null)
                    InitCoFaces();
                return ((int[,])(m_CoFaceVerticeIndices.Clone()));
            }
        }

        /// <summary>
        /// For each co-face, the indices of the faces whose geometric intersection forms the <em>iCoFace</em>-th co-face.<br/>
        /// 1st index: co-face index <em>iCoFace</em>  <br/>
        /// 2nd index: in {0,1}, corresponds to first and second face.
        /// </summary>
        /// <remarks>
        /// Note that each co-face can be described as the geometric intersection of 2 faces.
        /// </remarks>
        public int[,] CoFace_FaceIndices {
            get {
                if (m_CoFaces_FaceIndices == null)
                    InitCoFaces();
                return ((int[,])(m_CoFaces_FaceIndices.Clone()));
            }
        }

        /// <summary>
        /// number of co-faces
        /// </summary>
        public int NoOfCoFaces {
            get {
                return CoFace_FaceIndices.GetLength(0);
            }
        }


        AffineManifold[] m_FacePlanes;

        /// <summary>
        /// For each face, the affine manifold that represents the plane in which the face is located;<br/>
        /// </summary>
        /// <param name="iFace">
        /// face index
        /// </param>
        public AffineManifold GetFacePlane(int iFace) {
            if (m_FacePlanes == null) {
                m_FacePlanes = new AffineManifold[this.NoOfFaces];

                var fCen = this.FaceCenters;
                var fNor = this.FaceNormals;
                for (int iF = 0; iF < this.NoOfFaces; iF++) {
                    m_FacePlanes[iF] = new AffineManifold(fNor.GetRow(iF), fCen.GetRow(iF));

                    Debug.Assert(m_FacePlanes[iF].PointDistance(fCen.GetRow(iF)).Abs() < 1.0e-6);
                }
            }

            return m_FacePlanes[iFace].CloneAs();
        }
    
    }
}
