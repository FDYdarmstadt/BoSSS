using System.Collections;
using System.Collections.Generic;
using System;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class Divider<T>
    {
        Mesh<T> mesh;

        InsideCellEnumerator<T> insideCells;

        public Divider(Mesh<T> mesh, int firstCell_NodeIndice)
        {
            insideCells = new InsideCellEnumerator<T>(mesh, firstCell_NodeIndice);
            this.mesh = mesh;
        }

        public Divider(Mesh<T> mesh)
        {
            insideCells = new InsideCellEnumerator<T>(mesh);
            this.mesh = mesh;
        }

        public MeshCell<T> GetFirst(BoundaryLine boundaryLine)
        {
            insideCells.SetFirstCell((Vector)boundaryLine.Start);
            return insideCells.GetFirstCell();
        }

        public void RemoveOutsideCells()
        {
            IdentifyInsideCells();
            List<MeshCell<T>> insideCells = CollectInsideCells();
            mesh.Cells = insideCells;
            List<T> nodes = CollectNodes(insideCells);
            mesh.Nodes = nodes;
        }

        List<MeshCell<T>> CollectInsideCells()
        {
            List<MeshCell<T>> cells = new List<MeshCell<T>>(mesh.Cells.Count);
            foreach(MeshCell<T> cell in mesh.Cells)
            {
                if(cell.type == MeshCellType.Inside)
                {
                    cells.Add(cell);
                }
            }
            return cells;
        }

        static List<T> CollectNodes(List<MeshCell<T>> cells)
        {
            List<T> nodes = new List<T>(cells.Count);
            foreach(MeshCell<T> cell in cells)
            {
                nodes.Add(cell.Node);
            }
            return nodes;
        }

        void IdentifyInsideCells()
        {
            foreach (MeshCell<T> cell in insideCells.Cells())
            {
                cell.type = MeshCellType.Inside;
            }
        }
    }
}
