using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Xml;

namespace BoSSS.Application.ExternalBinding.MatlabCutCellQuadInterface {

    [TestFixture]
    public static class MatlabCutCellQuadInterfaceTests {
        /// <summary>
        /// Basic Testing for external language binding.
        /// </summary>
        [Test]
        public static void circle2D() {
            var app = new BoSSS.Application.ExternalBinding.MatlabCutCellQuadInterface.MatlabCutCellQuadInterface();
            app.BoSSSInitialize();

            double[] xnodesCoarse = { 0, 0.5, 1 };
            double[] xnodesNormal = { 0, 0.25, 0.5, 0.75, 1 };
            var xnodes = xnodesCoarse;
            app.SetDomain(2, xnodes, xnodes);

            double[] center = {  0.5, 0.5 };
            double R = 0.2;
            _2D phiCircle = (double x, double y) => -((x - center[0])* (x - center[0]) + (y - center[1])* (y - center[1]) - R*R);

            app.Submit2DLevelSet(phiCircle);
            app.ProjectLevelSetWithGaussAndStokes(3);
            //app.CompileQuadRules(2, -1);
            app.PlotCurrentState(3);
            app.CompileQuadRules(3, 1);
            double tot = 0;
            for (int jCell = 0; jCell < Math.Pow(xnodes.Length,2); jCell++) {
			    var ret = app.GetQuadRules(jCell, 1);
                if(ret is null)
                    continue;

            app.Submit2DLevelSet(phiCircle);
            app.ProjectLevelSet(3);
            app.CompileQuadRules(2, -1);
            app.CompileQuadRules(2, 1);
			double analyticalSol = R*R*Math.PI;
			double area = 0.0;
			for (int cell = 0; cell < xnodes.Length * xnodes.Length; cell++) {
				var Q = app.GetQuadRules(cell, -1);

                if (Q is null)
                    continue;

				//WriteQuadratureAsVtk("cell" + cell + ".vtp", Q);

				for (int i = 0; i < Q.Lengths[0]; i++) {
					area += Q[i, 2];
				}
			}


			Console.WriteLine($"Error: {area-analyticalSol}, calculate area of circle:{area} analytical:{analyticalSol}");
			app.BoSSSFinalize();
            //app.WriteVolQuadRules(2);
        }

		/// <summary>
		/// Basic Testing for external language binding.
		/// </summary>
		[Test]
		public static void ellipse2D() {
			var app = new BoSSS.Application.ExternalBinding.MatlabCutCellQuadInterface.MatlabCutCellQuadInterface();
			app.BoSSSInitialize();

			double[] xnodesCoarse = { 0, 0.5, 1 };
			double[] xnodesNormal = { 0, 0.25, 0.5, 0.75, 1 };
			double[] xnodesPaper = GenericBlas.Linspace(-2.0, 2.0, 10 + 1);
			var xnodes = xnodesPaper;
			app.SetDomain(2, xnodes, xnodes);

			var rnd = new Random(5);
			double randomValue = 0;// rnd.NextDouble() * 0.2 - 0.1;
			Console.WriteLine($"Random value for perturbation: {randomValue}");
			double[] center = { 0.0, 0.0 };
			double a = 1.5, b = 0.75;
			_2D phiEllipse = (double x, double y) =>
				((x - center[0]) / a) * ((x - center[0]) / a)
			  + ((y - center[1]) / b) * ((y - center[1]) / b)
			  - 1.0 + randomValue / 10000;

			app.Submit2DLevelSet(phiEllipse);
			app.ProjectLevelSet(2);
			for (int o = 0; o < 5; o++) {
				app.CompileQuadRules(o, -1);
				app.CompileQuadRules(o, 1);
				double analyticalSol = Math.PI * a * b;  // = 1.125 * π ≈ 3.5343
				double area = 0.0;
				for (int cell = 0; cell < xnodes.Length * xnodes.Length; cell++) {
					var Q = app.GetQuadRules(cell, -1);

					if (Q is null)
						continue;

					//WriteQuadratureAsVtk("cell" + cell + ".vtp", Q);

					for (int i = 0; i < Q.Lengths[0]; i++) {
						area += Q[i, 2];
					}
				}
				Console.WriteLine(
				  $"{o}-order error: {(area - analyticalSol):0.00e+00}, " +
				  $"calculated area of ellipse: {area}, " +
				  $"analytical: {analyticalSol}"
				);


			}

			app.BoSSSFinalize();
		}

		static void WriteQuadratureAsVtk(string filePath, MultidimensionalArray rule) {
			// rule is a 2D Array: rows = points, columns = [x, y, (z,) weight]
			int nPoints = rule.GetLength(0);
			int nCols = rule.GetLength(1);
			int dim = nCols - 1;  // spatial dimension inferred

			var settings = new XmlWriterSettings { Indent = true };
			using (var writer = XmlWriter.Create(filePath, settings)) {
				writer.WriteStartDocument();
				writer.WriteStartElement("VTKFile");
				writer.WriteAttributeString("type", "PolyData");
				writer.WriteAttributeString("version", "0.1");
				writer.WriteAttributeString("byte_order", "LittleEndian");

				writer.WriteStartElement("PolyData");
				writer.WriteStartElement("Piece");
				writer.WriteAttributeString("NumberOfPoints", nPoints.ToString());
				writer.WriteAttributeString("NumberOfVerts", nPoints.ToString());
				writer.WriteAttributeString("NumberOfLines", "0");
				writer.WriteAttributeString("NumberOfStrips", "0");
				writer.WriteAttributeString("NumberOfPolys", "0");

				// Points
				writer.WriteStartElement("Points");
				writer.WriteStartElement("DataArray");
				writer.WriteAttributeString("type", "Float32");
				writer.WriteAttributeString("Name", "Points");
				writer.WriteAttributeString("NumberOfComponents", "3");
				writer.WriteAttributeString("format", "ascii");
				for (int i = 0; i < nPoints; i++) {
					string x = rule[i, 0].ToString(CultureInfo.InvariantCulture);
					string y = rule[i, 1].ToString(CultureInfo.InvariantCulture);
					string z = dim == 3
						? rule[i, 2].ToString(CultureInfo.InvariantCulture)
						: "0.0";
					writer.WriteString($"{x} {y} {z}\n");
				}
				writer.WriteEndElement(); // DataArray
				writer.WriteEndElement(); // Points

				// Verts
				writer.WriteStartElement("Verts");
				// connectivity
				writer.WriteStartElement("DataArray");
				writer.WriteAttributeString("type", "Int32");
				writer.WriteAttributeString("Name", "connectivity");
				writer.WriteAttributeString("format", "ascii");
				for (int i = 0; i < nPoints; i++)
					writer.WriteString($"{i}\n");
				writer.WriteEndElement();

				// offsets
				writer.WriteStartElement("DataArray");
				writer.WriteAttributeString("type", "Int32");
				writer.WriteAttributeString("Name", "offsets");
				writer.WriteAttributeString("format", "ascii");
				for (int i = 1; i <= nPoints; i++)
					writer.WriteString($"{i}\n");
				writer.WriteEndElement();
				writer.WriteEndElement(); // Verts

				// PointData (weights)
				writer.WriteStartElement("PointData");
				writer.WriteAttributeString("Scalars", "w");
				writer.WriteStartElement("DataArray");
				writer.WriteAttributeString("type", "Float32");
				writer.WriteAttributeString("Name", "w");
				writer.WriteAttributeString("NumberOfComponents", "1");
				writer.WriteAttributeString("format", "ascii");
				for (int i = 0; i < nPoints; i++) {
					string w = rule[i, dim].ToString(CultureInfo.InvariantCulture);
					writer.WriteString($"{w}\n");
				}
				writer.WriteEndElement();
				writer.WriteEndElement(); // PointData

				writer.WriteEndElement(); // Piece
				writer.WriteEndElement(); // PolyData
				writer.WriteEndElement(); // VTKFile
				writer.WriteEndDocument();
			}
		}

	}
}
