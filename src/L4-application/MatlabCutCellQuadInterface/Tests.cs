using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

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

			    Console.WriteLine($"jCell={jCell} with a rule of {ret.GetLength(1)} nodes");

			    for (int i = 0; i < ret.GetLength(0); i++) {
                    Console.WriteLine($"- node {i}: x={ret[i,0]}, y={ret[i,1]}, w={ret[i, 2]}");
					tot += ret[i, 2];
				}
			}
            Console.WriteLine($"Total area: {tot} vs analytical area={Math.PI*Math.Pow(R,2)}");
			app.BoSSSFinalize();
            //app.WriteVolQuadRules(2);
        }



    }
}
