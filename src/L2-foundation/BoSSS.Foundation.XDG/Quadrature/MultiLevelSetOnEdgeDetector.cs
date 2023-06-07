using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG.Quadrature {
    class MultiLevelSetOnEdgeDetector {
        GridData grddat;

        LevelSetTracker.LevelSetData[] levelSetData;
        int iLevSet; // Lev set that determines the hierarchy, e.g. if the "1" LevelSet is chosen and coincides with an edge, the cell should be treated as a "normal" CutCell for the "0" LevSet
        int jLevSet;
        int D;

        public MultiLevelSetOnEdgeDetector(LevelSetTracker.LevelSetData[] levelSetData, CombinedID id) {
            this.levelSetData = levelSetData;
            this.iLevSet = 1;
            this.jLevSet = 0;

            grddat = levelSetData[0].GridDat;
            D = grddat.SpatialDimension;
        }

        /// <summary>
        /// Determine cells in which the <see cref="iLevSet"/> lies on a face
        /// </summary>
        /// <param name="j"></param>
        /// <returns></returns>
        public bool IsSpecialCell(int j) {
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData[0].Region.m_LevSetCoincidingFaces;
            if (CoIncFaces == null)
                return false;
            if (CoIncFaces[j] == null)
                return false;

            foreach (var t in CoIncFaces[j]) {
                if (t.iLevSet == iLevSet)
                    return true; // jetzt geht der Spass los!
            }
            return false;
        }

        /// <summary>
        /// Determines if the cell lies in the part of <see cref="iLevSet"/> where <see cref="jLevSet"/> is active
        /// </summary>
        /// <param name="cell"></param>
        /// <returns></returns>
        public bool IsActiveCell(int cell) {
            if (iLevSet != jLevSet) {
                MultidimensionalArray phiAtCenter = MultidimensionalArray.Create(1, 1);
                levelSetData[iLevSet].LevelSet.Evaluate(cell, 1, grddat.Cells.GetRefElement(cell).Center, phiAtCenter);
                if (phiAtCenter[0, 0] < 0) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return true;
            }
        }

        /// <summary>
        /// Determines if the edge lies in the part of <see cref="iLevSet"/> where <see cref="jLevSet"/> is active
        /// </summary>
        /// <param name="j"></param>
        /// <returns></returns>
        public bool IsActiveEdge(int j) {
            int iCell, jCell;
            iCell = grddat.Edges.CellIndices[j, 0];
            jCell = grddat.Edges.CellIndices[j, 1];
            bool iActive = IsActiveCell(iCell);
            bool jActive = false;
            if(jCell >= 0) {
                jActive = IsActiveCell(jCell);
            }
            return iActive | jActive;
        }

        /// <summary>
        /// Determine edges in which one of the neighbors is a <see cref="IsSpecialCell(int)"/>
        /// </summary>
        /// <param name="j"></param>
        /// <returns></returns>
        public bool IsSpecialEdge(int j) {
            int iCell, jCell;
            iCell = grddat.Edges.CellIndices[j, 0];
            jCell = grddat.Edges.CellIndices[j, 1];
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData[0].Region.m_LevSetCoincidingFaces;
            if (CoIncFaces == null)
                return false;

            bool iSpecial = false;
            if (CoIncFaces[iCell] != null) {
                foreach (var t in CoIncFaces[iCell]) {
                    if (t.iLevSet == iLevSet)
                        iSpecial = true; // jetzt geht der Spass los!
                }
            }

            bool jSpecial = false;
            if (jCell >= 0 && jCell < grddat.Cells.NoOfLocalUpdatedCells) {
                if (CoIncFaces[jCell] == null)
                    return iSpecial;                
                foreach (var t in CoIncFaces[jCell]) {
                    if (t.iLevSet == iLevSet)
                        jSpecial = true; // jetzt geht der Spass los!
                }
            } else {
                jSpecial = true; // local boundary edge, act as if other cell is special
            }
            return iSpecial | jSpecial; // if either neighbor cell is a special we want to employ the special treatment
        }

        ///// <summary>
        ///// return the edge index, for the special edge in a certain cell and the corresponding Edge2CellTrafoIndex
        ///// </summary>
        //public (int iEdge, int e2ctrfindex) GetSpecialEdge(int i) {
        //    (int iLevSet, int iFace)[][] CoIncFaces = levelSetData[0].Region.m_LevSetCoincidingFaces;
        //    if (CoIncFaces == null)
        //        return (-1, -1);
        //    if (CoIncFaces[i] == null)
        //        return (-1, -1);

        //    foreach (var t in CoIncFaces[i]) {
        //        if (t.iLevSet == iLevSet) {
        //            int iEdge = grddat.CellToEdge(i, t.iFace);
        //            int trf;
        //            if (grddat.iGeomEdges.CellIndices[iEdge, 0] == i) {
        //                trf = grddat.iGeomEdges.Edge2CellTrafoIndex[iEdge, 0];
        //            } else {
        //                trf = grddat.iGeomEdges.Edge2CellTrafoIndex[iEdge, 1];
        //            }
        //            return (iEdge, trf); // jetzt geht der Spass los!
        //        }
        //    }

        //    return (-1, -1);
        //}

        /// <summary>
        /// return the face index, for the special face in a certain cell
        /// </summary>
        public int GetSpecialFace(int i) {
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData[0].Region.m_LevSetCoincidingFaces;
            if (CoIncFaces == null)
                return -1;
            if (CoIncFaces[i] == null)
                return -1;

            foreach (var t in CoIncFaces[i]) {
                if (t.iLevSet == iLevSet) {
                    return t.iFace; // jetzt geht der Spass los!
                }
            }

            return -1;
        }
    }    
}
