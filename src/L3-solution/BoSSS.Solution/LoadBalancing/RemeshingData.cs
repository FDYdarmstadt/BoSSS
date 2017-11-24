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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {

    /// <summary>
    /// During mesh/grid adaptation (refinement or coarsening) <see cref="Application{T}.MpiRedistributeAndMeshAdapt(int, double)"/>, this is used to
    /// - backup/serialize objects on the original grid
    /// - restore/serialize objects on the refined grid
    /// </summary>
    public class RemeshingData : GridUpdateData {

        /// <summary>
        /// Stores the backup data of DG Fields after mesh adaptation.
        /// - dictionary key: string reference to identify DG field
        /// - 1st value index: local cell index, new grid.
        /// - 2nd value index: enumeration of cells; for refined and conserved cells, this contains only one element; for cells which are coarsened, this corresponds to all cells of the coarsening cluster.
        /// - 3rd value index: DG coordinate index within cell.
        /// </summary>
        Dictionary<string, double[][][]> m_newDGFieldData_GridAdapta;

        /// <summary>
        /// How the cells in the old mesh relate to cells in the new mesh.
        /// </summary>
        GridCorrelation m_Old2NewCorr;

        /// <summary>
        /// Apply the resorting (including mesh adaptation).
        /// </summary>
        public void Resort(GridCorrelation r, GridData NewGrid) {
            using(new FuncTrace()) {
                this.GridAdaptation = false;
                m_Old2NewCorr = r;
                m_OldGrid = null;
                m_OldTracker = null;
                m_NewGrid = NewGrid;
                m_newJ = m_NewGrid.iLogicalCells.NoOfLocalUpdatedCells;

                // fix the sequence in which we serialize/de-serialize fields
                string[] FieldNames = this.m_oldDGFieldData.Keys.ToArray();
                int NoFields = FieldNames.Length;

                // format data for serialization
                double[][] oldFieldsData = new double[m_oldJ][]; // 1st: cell; 2nd enum
                {

                    double[][][] oldFields = FieldNames.Select(rstring => this.m_oldDGFieldData[rstring]).ToArray(); // 1st: field; 2nd: cell; 3rd: DG coord

                    for(int j = 0; j < m_oldJ; j++) {
                        int Ntot = 0;
                        for(int iF = 0; iF < NoFields; iF++)
                            Ntot += oldFields[iF][j].Length;

                        // pack the DG coords for each field into one array of the form
                        // { Np_f1, DGcoords_f1, Np_f2, DGcoords_f2, .... };

                        double[] data_j = new double[Ntot + NoFields];
                        int cnt = 0;
                        for(int iF = 0; iF < NoFields; iF++) {
                            double[] data_iF_j = oldFields[iF][j];
                            int Nj = data_iF_j.Length;
                            data_j[cnt] = Nj;
                            cnt++;
                            for(int n = 0; n < Nj; n++) {
                                data_j[cnt] = data_iF_j[n];
                                cnt++;
                            }
                        }
                        oldFieldsData[j] = data_j;
                    }
                    oldFields = null;
                    m_oldDGFieldData = null;
                }

                // redistribute data
                double[][][] newFieldsData = new double[m_newJ][][];
                {
                    //Resorting.ApplyToVector(oldFieldsData, newFieldsData, m_NewGrid.CellPartitioning);
                    r.ApplyToVector(oldFieldsData, newFieldsData, m_NewGrid.CellPartitioning);
                    oldFieldsData = null;
                    Debug.Assert(newFieldsData.Length == m_newJ);
                }

                // re-format re-distributed data
                m_newDGFieldData_GridAdapta = new Dictionary<string, double[][][]>();
                {
                    double[][][][] newFields = new double[FieldNames.Length][][][]; // indices: [field, cell idx, enum over coarsening, DG mode]
                    for(int iF = 0; iF < NoFields; iF++) {
                        newFields[iF] = new double[m_newJ][][];
                    }
                    for(int j = 0; j < m_newJ; j++) {
                        double[][] data_j = newFieldsData[j];
                        int L = data_j.Length; // cell cluster size (equal 1 for refined or conserved cells)

                        for(int l = 0; l < L; l++) {
                            double[] data_jl = data_j[l];

                            int cnt = 0;
                            for(int iF = 0; iF < NoFields; iF++) {
                                int Nj = (int)data_jl[cnt];
                                cnt++;
                                double[] data_iF_jl = new double[Nj];
                                for(int n = 0; n < Nj; n++) {
                                    data_iF_jl[n] = data_jl[cnt];
                                    cnt++;
                                }
                                newFields[iF][j][l] = data_iF_jl;
                            }
                            Debug.Assert(cnt == data_jl.Length);
                        }
                    }

                    for(int iF = 0; iF < NoFields; iF++) {
                        m_newDGFieldData_GridAdapta.Add(FieldNames[iF], newFields[iF]);
                    }
                }
            }

        }



        /// <summary>
        /// Loads the DG coordinates after grid-redistribution.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        public override void RestoreDGField(DGField f, string Reference) {
            using(new FuncTrace()) {
                int newJ = this.m_newJ;

                GridData NewGrid = (GridData)m_NewGrid;

                int pDeg = f.Basis.Degree; //  Refined_TestData.Basis.Degree;

                //int newJ = RefinedGrid.Cells.NoOfLocalUpdatedCells;
                int[][] TargMappingIdx = new int[newJ][];

                m_Old2NewCorr.GetTargetMappingIndex(TargMappingIdx, NewGrid.CellPartitioning);

                double[][][] ReDistDGCoords = m_newDGFieldData_GridAdapta[Reference];
                Debug.Assert(ReDistDGCoords.Length == newJ);

                //double[][][] ReDistDGCoords = new double[newJ][][];
                //Old2NewCorr.ApplyToVector(TestData_DGCoordinates, ReDistDGCoords, RefinedGrid.CellPartitioning);

                for(int j = 0; j < newJ; j++) {
                    if(TargMappingIdx[j] == null) {
                        Debug.Assert(ReDistDGCoords[j].Length == 1);
                        f.Coordinates.SetRow(j, ReDistDGCoords[j][0]);
                    } else {
                        Debug.Assert(ReDistDGCoords[j].Length == TargMappingIdx[j].Length);

                        int iKref = NewGrid.Cells.GetRefElementIndex(j);

                        if(TargMappingIdx[j].Length == 1) {
                            // refinement

                            MultidimensionalArray Trafo = m_Old2NewCorr.GetSubdivBasisTransform(iKref, TargMappingIdx[j][0], pDeg);

                            double[] Coords_j = f.Coordinates.GetRow(j);
                            Trafo.gemv(1.0, ReDistDGCoords[j][0], 1.0, Coords_j, transpose: false);
                            f.Coordinates.SetRow(j, Coords_j);

                            //Console.WriteLine("Refine projection...");

                        } else {
                            // coarsening

                            int L = ReDistDGCoords[j].Length;
                            double[] Coords_j = f.Coordinates.GetRow(j);
                            for(int l = 0; l < L; l++) {
                                var Trafo = m_Old2NewCorr.GetSubdivBasisTransform(iKref, TargMappingIdx[j][l], pDeg);
                                Trafo.gemv(1.0, ReDistDGCoords[j][l], 1.0, Coords_j, transpose: true);
                            }
                            f.Coordinates.SetRow(j, Coords_j);

                            //Console.WriteLine("Coarsen projection...");
                        }

                    }
                }
            }
        }
    }
}