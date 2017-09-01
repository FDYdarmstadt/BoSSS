using BoSSS.Foundation;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System.Drawing;
using System.Globalization;
using System.Windows.Forms;

namespace MatrixVisualizer {

    /// <summary>
    /// Provides functionality to display instances of
    /// <see cref="MultidimensionalArray"/>, <see cref="MsrMatrix"/>,
    /// <see cref="NodeSet"/> and <see cref="BlockDiagonalMatrix"/> in
    /// a <see cref="DataGridView"/>s
    /// </summary>
    public static class VisualizationProvider {

        /// <summary>
        /// Interprets the given <paramref name="obj"/> as some kind of matrix
        /// and displays it in a <see cref="DataGridView"/> using 
        /// <see cref="MatrixGridViewForm"/>
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public static Form GetForm(object obj) {
            MatrixGridViewForm form = new MatrixGridViewForm();
            DataGridView grid = (DataGridView)form.Controls["gridView"];

            IMatrix matrix;
            string errorMessage = null;
            if (obj is MsrMatrix) {
                matrix = obj.As<MsrMatrix>().ToFullMatrixOnProc0();
            } else {
                matrix = obj as IMatrix;

                if (matrix is MultidimensionalArray && matrix.As<MultidimensionalArray>().Dimension != 2) {
                    errorMessage = "Given 'MultidimensionalArray' is not two-dimensional. Consider using 'ExtractSubArrayShallow' to extract a two-dimensional sub-array";
                }
            }

            if (matrix == null) {
                errorMessage = "Could not determine matrix type. This should not have happened";
            }

            if (errorMessage == null) {
                grid.ColumnCount = matrix.NoOfCols;
                grid.RowCount = matrix.NoOfRows;

                for (int i = 0; i < matrix.NoOfRows; i++) {
                    grid.Rows[i].HeaderCell.Value = i.ToString();

                    for (int j = 0; j < matrix.NoOfCols; j++) {
                        grid.Columns[j].HeaderText = j.ToString();

                        if (matrix[i, j] == 0.0) {
                            // Make cells with zeroth look different
                            grid.Rows[i].Cells[j].Value = 0.0;
                            grid.Rows[i].Cells[j].Style.BackColor = Color.LightGray;
                        } else {
                            grid.Rows[i].Cells[j].Value = matrix[i, j].ToString("e", NumberFormatInfo.InvariantInfo);
                        }

                    }
                }
            } else {
                grid.ColumnCount = 1;
                grid.RowCount = 1;
                grid.Rows[0].Cells[0].Value = errorMessage;
            }

            return form;
        }

        /// <summary>
        /// Transforms the given <paramref name="origTarget"/> into some
        /// serializable format.
        /// </summary>
        /// <param name="origTarget"></param>
        /// <returns></returns>
        public static object GetTarget(object origTarget) {
            if (origTarget is IMutableMatrixEx) {
                // Hack: Type is not serializable, so use full matrix instead
                origTarget = origTarget.As<IMutableMatrixEx>().ToFullMatrixOnProc0();
            } else if (origTarget is NodeSet) {
                NodeSet nodeSet = origTarget.Cast<NodeSet>();
                MultidimensionalArray temp = MultidimensionalArray.Create(nodeSet.Lengths);
                temp.Acc(1.0, nodeSet);
                origTarget = temp;
            }

            return origTarget;
        }
    }
}
