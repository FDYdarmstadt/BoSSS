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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Solution;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace ProjectionTest {

    public class Program : Application {

        private int noOfCellsX = -1;

        private int noOfNodesPerEdge = -1;

        private int dgDegree = -1;

        private GridData gridData;

        private NodeSet nodeSet;

        private int aspectRatio;

        public static void Main(string[] args) {
            BoSSS.Solution.Application._Main(args, true, () => new Program());
        }

        protected override GridCommons CreateOrLoadGrid() {
            // Options
            string dataPath = @"\\fdyprime\userspace\mueller\Jan\waves\";
            string databasePath = @"e:\bosss_db\GridOfTomorrow\";

            noOfCellsX = 64;
            noOfNodesPerEdge = 4;
            dgDegree = 6;

            //noOfCellsX = 85;
            //noOfNodesPerEdge = 4;
            //dgDegree = 6;

            Directory.SetCurrentDirectory(dataPath);
            DatabaseInfo database = new DatabaseInfo(databasePath);
            string sessionDescription = "Blafasel";

            // Create the grid
            aspectRatio = 4;
            GridCommons grid = Grid2D.Cartesian2DGrid(
                GenericBlas.Linspace(0.0, 1.0, noOfCellsX + 1),
                GenericBlas.Linspace(-1.5, 2.5, aspectRatio * noOfCellsX + 1),
                CellType.Square_Linear,
                periodicX: true,
                periodicY: false);
            grid.EdgeTagNames.Add(1, "subsonicOutlet");
            grid.EdgeTagNames.Add(2, "supersonicInlet");
            grid.DefineEdgeTags(delegate (double[] x) {
                if (x[1] > 0.0) {
                    return 2;
                } else {
                    return 1;
                }
            });
            gridData = new GridData(grid);

            // Read the values
            MultidimensionalArray rho = ReadVariableValues("rho");
            MultidimensionalArray u = ReadVariableValues("vx1");
            MultidimensionalArray v = ReadVariableValues("vx2");
            MultidimensionalArray p = ReadVariableValues("prs");

            // Assemble node set corresponding to the ordering from Matlab
            double[] nodes1D = GenericBlas.Linspace(-1.0, 1.0, noOfNodesPerEdge);
            int noOfNodesPerCell = noOfNodesPerEdge * noOfNodesPerEdge;
            double[,] localNodes = new double[noOfNodesPerCell, gridData.SpatialDimension];
            for (int i = 0; i < noOfNodesPerEdge; i++) {
                for (int j = 0; j < noOfNodesPerEdge; j++) {
                    int localNodeIndex = i * noOfNodesPerEdge + j;
                    localNodes[localNodeIndex, 0] = nodes1D[i];
                    localNodes[localNodeIndex, 1] = nodes1D[j];
                }
            }
            nodeSet = new NodeSet(gridData.Grid.RefElements[0], localNodes);

            // Interpolate
            //SinglePhaseField[] fields = LeastSquaresInterpolation(rho, u, v, p);
            SinglePhaseField[] fields = SpecFEMInterpolation(rho, u, v, p);

            // Save everything
            database.Controller.DBDriver.SaveGrid(grid, database);

            SessionInfo session = database.Controller.DBDriver.CreateNewSession(database);
            session.Description = sessionDescription;
            session.Save();

            database.Controller.DBDriver.SaveTimestep(
                0.0, 0, session, gridData, fields);

            return grid;
        }

        private MultidimensionalArray ReadVariableValues(string variableName) {
            List<double> list = new List<double>();
            using (var stream = new BinaryReader(new FileStream(variableName + ".dbl", FileMode.Open), Encoding.ASCII)) {
                while (stream.PeekChar() != -1) {
                    list.Add(stream.ReadDouble());
                }
            }

            int noOfNodesX = noOfCellsX * (noOfNodesPerEdge - 1) + 1;
            int noOfNodesY = aspectRatio * noOfCellsX * (noOfNodesPerEdge - 1) + 1;
            Debug.Assert(noOfNodesX * noOfNodesY == list.Count);

            return MultidimensionalArray.CreateWrapper(
                list.ToArray(), noOfNodesY, noOfNodesX);
        }

        public SinglePhaseField[] LeastSquaresInterpolation(MultidimensionalArray rhoValues, MultidimensionalArray uValues, MultidimensionalArray vValues, MultidimensionalArray pValues) {
            // Project DG fields
            SinglePhaseField rhoField = new SinglePhaseField(new Basis(gridData, dgDegree), "rho");
            rhoField.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => rho),
                nodeSet);

            SinglePhaseField u0Field = new SinglePhaseField(new Basis(gridData, dgDegree), "u0");
            u0Field.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => u),
                nodeSet);

            SinglePhaseField u1Field = new SinglePhaseField(new Basis(gridData, dgDegree), "u1");
            u1Field.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => v),
                nodeSet);

            SinglePhaseField pField = new SinglePhaseField(new Basis(gridData, dgDegree), "p");
            pField.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => p),
                nodeSet);

            SinglePhaseField m0Field = new SinglePhaseField(new Basis(gridData, dgDegree), "m0");
            m0Field.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => rho * u),
                nodeSet);

            SinglePhaseField m1Field = new SinglePhaseField(new Basis(gridData, dgDegree), "m1");
            m1Field.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => rho * v),
                nodeSet);

            SinglePhaseField rhoEField = new SinglePhaseField(new Basis(gridData, dgDegree), "rhoE");
            rhoEField.ProjectNodal(
                1.0,
                GetScalarFunction(rhoValues, uValues, vValues, pValues, (rho, u, v, p) => p / 0.4 + 0.5 * rho * (u * u + v * v)),
                nodeSet);

            return new SinglePhaseField[] { rhoField, u0Field, u1Field, pField, m0Field, m1Field, rhoEField };
        }

        private ScalarFunctionEx GetScalarFunction(MultidimensionalArray rho, MultidimensionalArray u, MultidimensionalArray v, MultidimensionalArray p, Func<double, double, double, double, double> func) {
            int noOfNodesPerCell = noOfNodesPerEdge * noOfNodesPerEdge;
            return delegate(int j0, int Len, NodeSet nodes, MultidimensionalArray result) {
                Debug.Assert(noOfNodesPerCell == nodes.GetLength(0));

                for (int i = 0; i < Len; i++) {
                    int cell = i + j0;

                    for (int j = 0; j < noOfNodesPerCell; j++) {
                        var indices = GetNodeIndices(cell, j);
                        result[i, j] = func(
                            rho[indices.Item1, indices.Item2],
                            u[indices.Item1, indices.Item2],
                            v[indices.Item1, indices.Item2],
                            p[indices.Item1, indices.Item2]);
                    }
                }
            };
        }

        public SinglePhaseField[] SpecFEMInterpolation(MultidimensionalArray rhoValues, MultidimensionalArray uValues, MultidimensionalArray vValues, MultidimensionalArray pValues) {
            SpecFemBasis basis = new SpecFemBasis(gridData, noOfNodesPerEdge);

            SinglePhaseField rhoField = SpecFEMInterpolation("rho", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => rho);
            SinglePhaseField u0Field = SpecFEMInterpolation("u0", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => u);
            SinglePhaseField u1Field = SpecFEMInterpolation("u1", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => v);
            SinglePhaseField m0Field = SpecFEMInterpolation("m0", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => rho * u);
            SinglePhaseField m1Field = SpecFEMInterpolation("m1", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => rho * v);
            SinglePhaseField pField = SpecFEMInterpolation("p", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => p);
            SinglePhaseField rhoEField = SpecFEMInterpolation("rhoE", basis, rhoValues, uValues, vValues, pValues, (rho, u, v, p) => p / 0.4 + 0.5 * rho * (u * u + v * v));

            return new SinglePhaseField[] { rhoField, u0Field, u1Field, pField, m0Field, m1Field, rhoEField };
        }

        private SinglePhaseField SpecFEMInterpolation(string fieldName, SpecFemBasis basis, MultidimensionalArray rho, MultidimensionalArray u, MultidimensionalArray v, MultidimensionalArray p, Func<double, double, double, double, double> func) {
            SpecFemField specFEMField = new SpecFemField(basis);
            int noOfNodesPerCell = noOfNodesPerEdge * noOfNodesPerEdge;
            int[,] cellNodeToNode = specFEMField.Basis.CellNode_To_Node;

            Debug.Assert(cellNodeToNode.GetLength(1) == noOfNodesPerCell);
            for (int i = 0; i < gridData.Grid.NoOfUpdateCells; i++) {
                for (int j = 0; j < noOfNodesPerCell; j++) {
                    int matlabNodeIndex = GetLocalMatlabNodeIndex(
                        basis.CellNodes[0][j, 0], basis.CellNodes[0][j, 1]);

                    var indices = GetNodeIndices(i, matlabNodeIndex);
                    specFEMField.Coordinates[cellNodeToNode[i, j]] = func(
                        rho[indices.Item1, indices.Item2],
                        u[indices.Item1, indices.Item2],
                        v[indices.Item1, indices.Item2],
                        p[indices.Item1, indices.Item2]);
                }
            }

            SinglePhaseField dgField = new SinglePhaseField(new Basis(gridData, dgDegree), fieldName);
            specFEMField.AccToDGField(1.0, dgField);
            return dgField;
        }

        private int GetLocalMatlabNodeIndex(double x, double y) {
            int index = -1;
            for (int i = 0; i < nodeSet.NoOfNodes; i++) {
                double dist = (nodeSet[i, 0] - x) * (nodeSet[i, 0] - x) +
                    (nodeSet[i, 1] - y) * (nodeSet[i, 1] - y);

                if (dist < 1e-12) {
                    index = i;
                    break;
                }
            }

            if (index < 0) {
                throw new Exception("Could not find node");
            }

            return index;
        }

        private Tuple<int, int> GetNodeIndices(int cell, int cellNode) {
            int xCell = cell / (aspectRatio * noOfCellsX);
            int yCell = cell % (aspectRatio * noOfCellsX);

            int xNode = cellNode / noOfNodesPerEdge;
            int yNode = cellNode % noOfNodesPerEdge;

            int xNodeIndex = xNode + xCell * (noOfNodesPerEdge - 1);
            int yNodeIndex = yNode + yCell * (noOfNodesPerEdge - 1);

            // Beware of transposition
            return new Tuple<int, int>(yNodeIndex, xNodeIndex);
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            this.TerminationKey = true;
            return dt;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            throw new NotImplementedException();
        }
    }
}
