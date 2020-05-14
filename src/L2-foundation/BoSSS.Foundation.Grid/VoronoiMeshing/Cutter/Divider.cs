using System.Collections;
using System.Collections.Generic;
using System;
using ilPSP;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Cutter
{
    class Divider<T>
    {
        Mesh<T> mesh;

        InsideCellEnumerator<T> insideCells;

        public Divider(Mesh<T> mesh, Boundary<T> boundary)
        {
            insideCells = new InsideCellEnumerator<T>(mesh, boundary);
            this.mesh = mesh;
        }

        public Divider(Domain<T> mesh)
        {
            insideCells = new InsideCellEnumerator<T>(mesh);
            this.mesh = mesh.Mesh;
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
                if(cell.Type == MeshCellType.Inside)
                {
                    cell.ID = cells.Count;
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
            foreach (MeshCell<T> cell in insideCells.EnumerateCellsInConcentricCircles())
            {
                cell.Type = MeshCellType.Inside;
            }
        }
    }
}
