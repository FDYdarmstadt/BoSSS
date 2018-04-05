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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools.FastMarcher;

namespace BoSSS.Solution.LevelSetTools.FastMarching.GlobalMarcher {

    /// <summary>
    /// Wrapperclass for Fastmarcher. 
    ///This class turns the information contained in GridData into a Graph that the fastmarcher needs.
    ///The value that is used by the fastmarcher to find a marching order is the meanvalue of phi of a cell.
    ///Usage :   Setup the general graph properties with the method Initialize(...). 
    ///          Use BuildInitialAcceptedCells(...) to get the initial accepted Nodes in the Graph for FastMarching. 
    /// 
    /// </summary>

    class MarchingCell : IMarchingNode {

        static ILocalSolver fMSolver;
        static BitArray inUseMask;
        static BitArray reinitField;
        static int[] queueIDList;
        static SinglePhaseField phi;
        static GridData gridDat;

        int jCell;
        MarchingCell[] neighbors = null;
        double phi_CellAverage;

        MarchingCell(int JCell) {
            jCell = JCell;
        }

        MarchingCell(int JCell, double Phi_CellAverage) {
            jCell = JCell;
            phi_CellAverage = Phi_CellAverage;
        }

        //General information
        public static void Initialize(ILocalSolver FMSolver, SinglePhaseField Phi, GridData GridDat, CellMask ReinitField) {
            phi = Phi;
            fMSolver = FMSolver;
            phi = Phi;
            gridDat = GridDat;
            inUseMask = new BitArray(GridDat.Cells.NoOfCells);
            reinitField = ReinitField.GetBitMask();
            queueIDList = new int[GridDat.Cells.NoOfCells];
        }

        //Create an array of Cells from the Accepted CellMask. 
        public static MarchingCell[] BuildInitialAcceptedCells(CellMask Accepted) {
            int[] AcceptedCells_Indices = Accepted.ItemEnum.ToArray();
            int AcceptedCells_Total = AcceptedCells_Indices.Length;
            MarchingCell[] AcceptedCells = new MarchingCell[AcceptedCells_Total];

            for (int i = 0; i < AcceptedCells_Total; ++i) {
                int jCell = AcceptedCells_Indices[i];
                double phi_average = phi.GetMeanValue(jCell);
                AcceptedCells[i] = new MarchingCell(jCell, phi_average);
                AcceptedCells[i].Accept();
            }

            return AcceptedCells;
        }

        //Methods for IMarchingNode
        #region IMarchingNode

        public double Value {
            get {
                return phi_CellAverage;
            }
        }

        public void Accept() {
            inUseMask[jCell] = true; //Update inUseMask
        }

        //Solve in this Cell. Then Calculate the MeanValue for Cellwise fast marching.
        public void CalculateValue() {
            fMSolver.LocalSolve(jCell, inUseMask, phi);
            phi_CellAverage = phi.GetMeanValue(jCell);
        }

        //This method returns the neighbors of each Cell. Each call will wrap the neighboring Datagrid-Cells into a Cell. 
        //To make sure that there aren't any duplicates in the graph, the GetHashCode and Equals methods are overwritten.
        public IMarchingNode[] Neighbors {
            get {
                //GetNeighbors via BoSSS Framework
                //Only calculate neighbors once:
                if (neighbors == null) {
                    Tuple<int, int, int>[] bossNeighbors = gridDat.GetCellNeighboursViaEdges(jCell);
                    LinkedList<MarchingCell> neighborsList = new LinkedList<MarchingCell>();
                    for (int i = 0; i < bossNeighbors.Length; ++i) {
                        int jCellNeighbor = bossNeighbors[i].Item1;
                        if (reinitField[jCellNeighbor]) {   //only use cells in reinitField
                            neighborsList.AddLast( new MarchingCell(jCellNeighbor));
                        } 
                    }
                    neighbors = neighborsList.ToArray();
                }
                return neighbors;
            }
        }

        //Each cell with indices jCell has only one ID (this ID is used/changed by the fastmarcher). 
        public int QueueID {
            get {
                return queueIDList[jCell];
            }

            set {
                queueIDList[jCell] = value;
            }
        }

        //Overwrite hashing stuff- each cell is identified by its indice jCell
        public override int GetHashCode() {
            return jCell;
        }

        //Cells are equal iff jCell matches
        public override bool Equals(object obj) {
            //If parameter cannot be cast to Cell return false:
            MarchingCell p = obj as MarchingCell;
            if ((object)p == null) {
                return false;
            }

            //Return true if the Cell indices match:
            return p.jCell == jCell;
        }

        #endregion

    }
}
