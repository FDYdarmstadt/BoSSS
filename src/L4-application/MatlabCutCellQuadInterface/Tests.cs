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
            var xnodes = xnodesNormal;
            app.SetDomain(2, xnodes, xnodes);

            double[] center = {  0.5, 0.5 };
            double R = 0.2;
            _2D phiCircle = (double x, double y) => (x - center[0])* (x - center[0]) + (y - center[1])* (y - center[1]) - R*R;

            app.SetLevelSets(3, phiCircle);
            app.CompileQuadRules(2, -1);
            app.CompileQuadRules(2, 1);

            app.GetQuadRules(3, 1);
            //app.WriteVolQuadRules(2);
        }



    }
}
