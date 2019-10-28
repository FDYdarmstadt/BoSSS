using ilPSP.Connectors.Matlab;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MatlabPlotter
    {
        /// <summary>
        /// MatlabColor 
        /// </summary>
        public string CellColor;

        /// <summary>
        /// https://de.mathworks.com/help/matlab/ref/scatter.html
        /// </summary>
        public string NodePlotSettings;

        public string Path; 

        public MatlabPlotter()
        {
            CellColor = "[0 0.4470 0.7410]";
            NodePlotSettings = "30, 'k'";
            Path = "./";
        }

        public void Plot<T>(Mesh<T> mesh, string name)
            where T : ILocatable
        {

            using (BatchmodeConnector bmc = new BatchmodeConnector(Path))
            {
                bmc.Cmd("hold on");
                Plot(bmc, mesh.Cells);
                Plot(bmc, mesh.Nodes);
                bmc.Cmd("hold off");
                bmc.Cmd($"saveas(gcf, '{name}.png')");
                bmc.Execute();
            }
        }

        void Plot<T>(BatchmodeConnector bmc, IList<MeshCell<T>> cells)
        {
            string matlabFillArgument = default(string);
            foreach (MeshCell<T> cell in cells)
            {
                try
                {
                    (double[] x, double[] y) cellCoordinates = GetCoordinatesOf(cell);
                    string x = ToMatlabArray(cellCoordinates.x);
                    string y = ToMatlabArray(cellCoordinates.y);
                    matlabFillArgument += x + "," + y + "," + CellColor + ",";
                }
                catch(Exception e)
                {
                    Console.WriteLine($"Can not plot cell {cell.ID}. {e}");
                }
                
            }
            matlabFillArgument = matlabFillArgument.Remove(matlabFillArgument.Length - 1);
            
            //Example: bmc.Cmd("fill([0 0 1],[0 1 0],'r', [0 1 1], [1 0 1], 'g')");
            bmc.Cmd($"fill({matlabFillArgument})");
        }

        static (double[] x, double[] y) GetCoordinatesOf<T>(MeshCell<T> cell)
        {
            double[] x = new double[cell.Vertices.Length];
            double[] y = new double[cell.Vertices.Length]; 
            for(int i = 0;  i < cell.Vertices.Length; ++i)
            {
                x[i] = cell.Vertices[i].Position.x;
                y[i] = cell.Vertices[i].Position.y;
            }
            return (x, y);
        }

        static string ToMatlabArray(IList<double> array)
        {
            string matlabArray = "[";
            for (int i = 0; i < array.Count; ++i)
            {
                matlabArray += array[i].ToString("G", CultureInfo.InvariantCulture) + " ";
            }
            matlabArray = matlabArray.Remove(matlabArray.Length - 1);
            matlabArray += "]";
            return matlabArray;
        }

        void Plot<T>(BatchmodeConnector bmc, IList<T> nodes)
            where T : ILocatable
        {
            try
            {
                (double[] x, double[] y) nodeCoordinates = GetCoordinatesOf(nodes);
                string x = ToMatlabArray(nodeCoordinates.x);
                string y = ToMatlabArray(nodeCoordinates.y);
                bmc.Cmd($"scatter({x}, {y}, {NodePlotSettings})");
            }
            catch (Exception e)
            {
                Console.WriteLine($"Can not plot nodes. {e}");
            }
        }

        static (double[] x, double[] y) GetCoordinatesOf<T>(IList<T> nodes)
            where T : ILocatable
        {
            double[] x = new double[nodes.Count];
            double[] y = new double[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                x[i] = nodes[i].Position.x;
                y[i] = nodes[i].Position.y;
            }
            return (x, y);
        }
    }
}
