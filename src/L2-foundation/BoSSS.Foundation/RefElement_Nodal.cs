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
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP.Utils;
using System.Collections;
using ilPSP;

namespace BoSSS.Foundation.Grid.RefElements {


    public partial class RefElement {

        /// <summary>
        /// Returns a node-set 
        /// </summary>
        /// <param name="px">
        /// Number of subdivisions in each direction.
        /// </param>
        /// <param name="Nodes">
        /// 2D-array, 1st index: node index, 2nd index: spatial direction 
        /// </param>
        /// <param name="Type">
        /// The location of the node within the cell; <br/>
        /// Index: node index;
        /// <list type="bullet">
        /// <item>in 1D: 0 for cells, 1 for edges</item>
        /// <item>in 2D: 0 for cells/volumes, 1 for edges/faces, 2 for corners</item>
        /// <item>in 3D: 0 for cells/volumes, 1 for faces, 2 for edges, 3 for corners</item>
        /// </list>
        /// </param>
        protected abstract void GetNodeSet(int px, out MultidimensionalArray Nodes, out int[] Type);

        /// <summary>
        /// Returns a node-set 
        /// </summary>
        /// <param name="P">
        /// Number of subdivisions in each direction, along one edge.
        /// </param>
        /// <param name="Nodes">
        /// 2D-array, 1st index: node index, 2nd index: spatial direction 
        /// </param>
        /// <param name="Type">
        /// The location of the node within the cell; <br/>
        /// Index: node index;
        /// <list type="bullet">
        /// <item>in 1D: 0 for cells (volume nodes, volume is 1D), 1 for vertices (vertex nodes, faces are 0D, i.e. points or vertices)</item>
        /// <item>in 2D: 0 for cells (volume nodes, volume is 2D), <br/>
        ///              1 for edges (face nodes, faces are 1D, i.e. lines),  <br/>
        ///              2 for corners (in 2D the vertices correlate with the co-faces, i.e. 0D objects)</item>
        /// <item>in 3D: 0 for cells (volume nodes, volume is 3D), <br/>
        ///              1 for faces (face nodes, the faces are 2D - polygons, e.g. rectangles or triangles),<br/>
        ///              2 for edges (co-face nodes: co-faces are the edges of the faces; these are lines) <br/>
        ///              3 for corners (etc.) </item>
        /// </list>
        /// </param>
        /// <param name="EntityIndex">
        /// If the k-th node is 
        /// <list type="bullet">
        /// <item>a type 0 nodes, i.e. a cell node, <paramref name="EntityIndex"/>[k] is always 0.</item>
        /// <item>a face node, <paramref name="EntityIndex"/>[k] is the index of the face within this reference element.</item>
        /// <item>a co-face node, <paramref name="EntityIndex"/>[k] is the index of the co-face within this reference element</item>
        /// <item>a vertex node, <paramref name="EntityIndex"/>[k] is the index of the vertex within this reference element</item>
        /// </list>
        /// </param>
		/// <param name="TypeFilter">
		/// Node types (see <paramref name="Type"/>) which should be filtered OUT, i.e. nodes with types in <paramref name="TypeFilter"/> will be omitted.
		/// </param>
        public virtual void GetNodeSet(int P, out NodeSet Nodes, out int[] Type, out int[] EntityIndex, params int[] TypeFilter) {
            MultidimensionalArray UnsortNodes;
            int[] UnsortType;
            int D = this.SpatialDimension;
            if (P <= 0)
                throw new ArgumentOutOfRangeException("illegal number of nodes");
            if (P > 1) {
                GetNodeSet(P, out UnsortNodes, out UnsortType);
            } else {
                UnsortNodes = MultidimensionalArray.Create(D+1, D);
                UnsortNodes.Set(this.Vertices.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { D, D - 1 }));
                UnsortType = new int[D+1];
                UnsortType.SetAll(D);
            }
            if (UnsortNodes.GetLength(0) != UnsortType.Length)
                throw new ApplicationException();
            if (UnsortNodes.GetLength(1) != D)
                throw new ApplicationException();

			if(TypeFilter == null) {
				TypeFilter = new int[0];
			} else {
				foreach(int tf in TypeFilter)
					if(tf < 0 || tf > D)
						throw new ArgumentOutOfRangeException("illegal type filter");
			}

            Nodes = new NodeSet(this, UnsortNodes.Lengths[0], UnsortNodes.Lengths[1]);
            int NoOfNodes = UnsortType.Length;
            Type = new int[NoOfNodes];
            EntityIndex = new int[NoOfNodes];

            int K = NoOfNodes - 1;
            BitArray Identified = new BitArray(NoOfNodes);

            // identify vertex nodes (1D, 2D and 3D)
            // =====================================
            {
                var V = this.Vertices;

                int Vert_type = D; // 3D: vertex/corner nodes are type 3
                //                 // 2D: vertex/corner nodes are type 2
                //                 // 1D: vertex/corner nodes are type 1 
                for (int iVertex = V.GetLength(0) - 1; iVertex >= 0; iVertex--) {
                    var v = V.GetRow(iVertex);

                    int idxFnd = int.MinValue;
                    double dist = double.MaxValue;
                    for (int i = 0; i < NoOfNodes; i++) {
                        if (Identified[i])
                            continue;
                        if (UnsortType[i] != Vert_type)
                            continue;
                        var vx = UnsortNodes.GetRow(i);

                        double a = GenericBlas.L2Dist(vx, v);
                        if (a < 1.0e-12) {
                            dist = a;
                            idxFnd = i;
                        }
                    }

                    if (idxFnd < 0) {
                        if (P > 1)
                            throw new ApplicationException("Error in implementation.");
                    } else {
                        Identified[idxFnd] = true;
                        Type[K] = UnsortType[idxFnd];
                        EntityIndex[K] = iVertex;
                        Nodes.SetRow(K, v);
                        K--;
                    }
                }
            }

            // identify Co-Face nodes (edges of faces, only 3D)
            // ================================================
            if (D == 3) {

                int Vert_type = 2; 

                for (int iCoFace = this.NoOfCoFaces - 1; iCoFace >= 0; iCoFace--) { // loop over co-faces
                    var plane0 = this.GetFacePlane(this.CoFace_FaceIndices[iCoFace, 0]);
                    var plane1 = this.GetFacePlane(this.CoFace_FaceIndices[iCoFace, 1]);


                    for (int i = 0; i < NoOfNodes; i++) { // loop over all nodes
                        if (Identified[i])
                            continue;
                        if (UnsortType[i] != Vert_type)
                            continue;
                        var vx = UnsortNodes.GetRow(i);

                        double a1 = plane0.PointDistance(vx);
                        double a2 = plane1.PointDistance(vx);

                        if (Math.Abs(a1) < 1.0e-10 && Math.Abs(a2) < 1.0e-10) {
                            Identified[i] = true;
                            Type[K] = UnsortType[i];
                            EntityIndex[K] = iCoFace;
                            Nodes.SetRow(K, vx);
                            K--;
                        }
                    }
                }
            }

            // identify Face-Nodes (on faces, 2D and 3D)
            // =========================================
            if (D == 2 || D == 3) {
                
                int Vert_type = 1;

                int NoOfFaces = this.NoOfFaces;
                for (int iFace = NoOfFaces - 1; iFace >= 0; iFace--) {

                    var offset = this.FaceCenters.GetRow(iFace);
                    var normal = this.FaceNormals.GetRow(iFace);
                    double dora = GenericBlas.InnerProd(offset, normal);

                    for (int i = 0; i < NoOfNodes; i++) {
                        if (Identified[i])
                            continue;
                        if (UnsortType[i] != Vert_type)
                            continue;
                        var vx = UnsortNodes.GetRow(i);

                        double a = GenericBlas.InnerProd(vx, normal) - dora;
                        if (Math.Abs(a) < 1.0e-10) {
                            Identified[i] = true;
                            Type[K] = UnsortType[i];
                            EntityIndex[K] = iFace;
                            Nodes.SetRow(K, vx);
                            K--;
                        }
                    }
                }
            }

            // identify Volume-Nodes (1D, 2D, 3D)
            // ==================================
            {
                int Vert_type = 0;

                for (int i = 0; i < NoOfNodes; i++) {
                    if (Identified[i]) {
                        if (UnsortType[i] == Vert_type)
                            throw new ApplicationException("Error in algorithm.");
                    } else {
                        if (UnsortType[i] != Vert_type)
                            throw new ApplicationException("Error in algorithm.");

                        var vx = UnsortNodes.GetRow(i);
                        Identified[i] = true;
                        Type[K] = UnsortType[i];
                        EntityIndex[K] = 0;
                        Nodes.SetRow(K, vx);
                        K--;
                    }
                }
            }
            Nodes.LockForever();

			// Apply optional filters
			// ======================
			if(TypeFilter.Length > 0) {
				MultidimensionalArray _Nodes = Nodes;
				int[] _Type = Type; 
				int[] _EntityIndex = EntityIndex;

                
				int NewNoOfNodes = Type.Where(ty => !TypeFilter.Contains(ty)).Count();

				Nodes = new NodeSet(this, NewNoOfNodes, D);
				Type = new int[NewNoOfNodes];
				EntityIndex = new int[NewNoOfNodes];

				int kk = 0;
				for(int k = 0; k < NoOfNodes; k++) {
					if(TypeFilter.Contains(_Type[k]))
						continue;

					Type[kk] = _Type[k];
					EntityIndex[kk] = _EntityIndex[kk];
					Nodes.ExtractSubArrayShallow(kk,-1).Set(_Nodes.ExtractSubArrayShallow(k,-1));

                    kk++;
				}


                Debug.Assert(kk == NewNoOfNodes);

                Nodes.LockForever();
			}
        }

        /// <summary>
        /// Creates a node-set and a corresponding set of nodal polynomials;
        /// </summary>
        public void SelectNodalPolynomials(int px, out NodeSet Nodes, out PolynomialList NodalBasis, out int[] Type, out int[] EntityIndex, Func<Polynomial, bool> ModalBasisSelector = null, int[] NodeTypeFilter = null) {
            MultidimensionalArray a, b;
            SelectNodalPolynomials(px, out Nodes, out NodalBasis, out  Type, out EntityIndex, out a, out b, ModalBasisSelector, NodeTypeFilter);
            Debug.Assert(ContainsZero(NodalBasis), "interpolation poly not found, 4");
        }


        /// <summary>
        /// Creates a node-set and a corresponding set of nodal polynomials;
        /// </summary>
        public void SelectNodalPolynomials(int px, out NodeSet Nodes, out PolynomialList _NodalBasis, out int[] Type, out int[] EntityIndex, out MultidimensionalArray Nodal2Modal, out MultidimensionalArray Modal2Nodal, Func<Polynomial, bool> ModalBasisSelector = null, int[] NodeTypeFilter = null) {

            int D = this.SpatialDimension;

            if (ModalBasisSelector== null)
                ModalBasisSelector = delegate(Polynomial p) {

                    Debug.Assert(p.Coeff.Length == p.Exponents.GetLength(0));
                    Debug.Assert(p.Exponents.GetLength(1) == D);


                    for (int l = 0; l < p.Coeff.Length; l++) {
                        for (int d = 0; d < D; d++) {
                            if (p.Exponents[l, d] > (px - 1))
                                return false;
                        }
                    }

                    return true;
                };

            // Find node set
            // =============
            GetNodeSet(px, out Nodes, out Type, out EntityIndex, NodeTypeFilter);
            int NoOfNodes = Type.Length;

#if DEBUG
            double[,] DIST = new double[NoOfNodes,NoOfNodes];
            for (int j1 = 0; j1 < NoOfNodes; j1++) {
                for (int j2 = 0; j2 < NoOfNodes; j2++) {
                    if(j2 == j1)
                        continue;

                    var node_j1 = Nodes.GetRow(j1);
                    var node_j2 = Nodes.GetRow(j2);

                    DIST[j1, j2] = GenericBlas.L2Dist(node_j1, node_j2);
                    if (DIST[j1, j2] <= 1.0e-8)
                        throw new ApplicationException("internal error.");

                }
            }
#endif


            // Find modal basis of approximation space
            // =======================================

            var ortho_polys = this.OrthonormalPolynomials;
            
            List<Polynomial> Basis = new List<Polynomial>();
            int deg;
            if (px > 1) {
                for (deg = 0; deg <= (px - 1)*D; deg++) {
                    var r = ortho_polys.Where(p => ((p.AbsoluteDegree == deg) && ModalBasisSelector(p)));

                    Basis.AddRange(r);
                    if (Basis.Count >= NoOfNodes)
                        break;
                }
            } else {
                var r = ortho_polys.Where(pol => pol.AbsoluteDegree <= 1);
                Basis.AddRange(r);
                deg = 1;
            }


            if (Basis.Count != NoOfNodes)
                throw new ApplicationException("Basis selection failed.");

            
            // check for basis selection
            // ==========================

            // case Triangle, Tetra: the complete P_{px-1}(x,y) should be selected
            // case Quad, Cube: P_{px-1}(x)*P_{px-1)(y) subset of P_{px-1}(x,y) will be selected
            int NoOfPolys = NoOfNodes;

            
            int[] idx_Basis = new int[NoOfPolys];
            for(int i = 0; i < NoOfPolys; i++) {
                idx_Basis[i] = Array.IndexOf(ortho_polys, Basis[i]);
                Debug.Assert(idx_Basis[i] >= 0);
            }

 
            // Construct nodal Polynomials
            // ===========================

            MultidimensionalArray PolyAtNodes = MultidimensionalArray.Create(NoOfPolys, NoOfNodes);
            for( int i = 0; i < NoOfNodes; i++) {
                Basis[i].Evaluate(PolyAtNodes.ExtractSubArrayShallow(-1, i), Nodes);
            }

            //FullMatrix MtxPolyAtNodes = new FullMatrix(NoOfPolys, NoOfNodes);
            //MtxPolyAtNodes.Set(PolyAtNodes);

            Debug.Assert(!PolyAtNodes.ContainsNanOrInf(), "Nodal Poly generation, illegal value");
            var Sol = PolyAtNodes.GetInverse();
            Debug.Assert(!Sol.ContainsNanOrInf(), "Nodal Poly generation, solution, illegal value");

            Modal2Nodal = MultidimensionalArray.Create(ortho_polys.Where(p => p.AbsoluteDegree <= deg).Count(), NoOfPolys);
            Nodal2Modal = MultidimensionalArray.Create(Modal2Nodal.NoOfCols, Modal2Nodal.NoOfRows);
            MultidimensionalArray _Modal2Nodal = MultidimensionalArray.Create(NoOfPolys, NoOfPolys);
            MultidimensionalArray _Nodal2Modal = MultidimensionalArray.Create(NoOfPolys, NoOfPolys);

            var b = new double[NoOfPolys];
            var rhs = new double[NoOfNodes];
            Polynomial[] NodalBasis = new Polynomial[NoOfNodes];
            for (int k = 0; k < NoOfNodes; k++) {
                Array.Clear(rhs, 0, rhs.Length);
                rhs[k] = 1.0;

                Sol.gemv(1.0, rhs, 0.0, b);

                for (int i = 0; i < NoOfPolys; i++) {
                    if (Math.Abs(b[i]) > 1.0e-12) {
                        if (NodalBasis[k] == null)
                            NodalBasis[k] = b[i]*Basis[i];
                        else
                            NodalBasis[k] = NodalBasis[k] + b[i]*Basis[i];
                    } else {
                        //Console.WriteLine("tresh {0}, {1}", Math.Abs(b[i]), b[i]);
                    }

                    _Modal2Nodal[i, k] = b[i];
                }

                //Console.WriteLine("k: {0}, NoOfPolys {1}, isnull {2} ", k, NoOfPolys, NodalBasis[k] == null);
            }
            _NodalBasis = new PolynomialList(NodalBasis);
            Debug.Assert(ContainsZero(_NodalBasis), "interpolation poly not found, 5");

            _Modal2Nodal.InvertTo(_Nodal2Modal);

            for (int k = 0; k < NoOfNodes; k++) {
                for (int i = 0; i < NoOfPolys; i++) {
                    Modal2Nodal[idx_Basis[i], k] = _Modal2Nodal[i, k];
                    Nodal2Modal[i, idx_Basis[k]] = _Nodal2Modal[i, k];
                
                }
            }


        }

    }
}
