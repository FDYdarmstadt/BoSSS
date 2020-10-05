using BoSSS.Solution.Gnuplot;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Formatting tools for <see cref="Plot2Ddata"/>.
    /// </summary>
    static public class Plot2DdataExtensions {

        /// <summary>
        /// Maybe useful for multi-plots, when just one plot should show the legend for all plots:
        /// It loops over all plots, collects all names and formats and adds a dummy graph to a specific
        /// plot if the name/format pair is not already shown in this plot.
        /// </summary>
        /// <param name="multiplots"></param>
        /// <param name="I">destination plot (where the dummys are added) row</param>
        /// <param name="J">destination plot (where the dummys are added) column</param>
        /// <param name="byName">
        /// - if true: go by <see cref="Plot2Ddata.XYvalues.Name"/>.
        /// - if false: go by <see cref="Plot2Ddata.XYvalues.Format"/>.
        /// </param>
        /// <param name="DummyX">x-value of the dummy plot to add, should be outside visible range, <see cref="Plot2Ddata.XrangeMax"/>.</param>
        /// <param name="DummyY">y-value of the dummy plot to add, should be outside visible range, <see cref="Plot2Ddata.YrangeMax"/>.</param>
        public static void AddDummyPlotsForLegend(this Plot2Ddata[,] multiplots, int I, int J, bool byName = true, double DummyX = 1e55, double DummyY = 1e56) {

            // collect all names & formats
            // ---------------------------

            var names = new List<string>();
            var fomts = new List<PlotFormat>();

            for (int i = 0; i < multiplots.GetLength(0); i++) {
                for (int j = 0; j < multiplots.GetLength(1); j++) {
                    if (multiplots[i, j] != null) {
                        foreach (var p in multiplots[i, j].dataGroups) {

                            bool isthere;
                            if (byName) {
                                isthere = names.Contains(p.Name);
                            } else {
                                isthere = fomts.Contains(p.Format);
                            }

                            if (!isthere) {
                                names.Add(p.Name);
                                fomts.Add(p.Format);
                            }

                        }
                    }
                }
            }
            // see what we have to add
            // -----------------------

            for(int iGraph = 0; iGraph < names.Count; iGraph++) {
                //foreach (var p in multiplots[I, J].dataGroups) {
                //
                //}
                string nmn = names[iGraph];
                PlotFormat fmt = fomts[iGraph];

                bool isThere;
                if(byName) {
                    isThere = multiplots[I, J].dataGroups.Where(gr => gr.Name.Equals(nmn)).Count() > 0;
                } else {
                    isThere = multiplots[I, J].dataGroups.Where(gr => gr.Format.Equals(fmt)).Count() > 0;
                }

                if(!isThere) {
                    Plot2Ddata.XYvalues dummy = new Plot2Ddata.XYvalues(nmn) {
                        Format = fmt,
                        Abscissas = new double[] { DummyX },
                        Values = new double[] { DummyY }
                    };

                    ArrayTools.AddToArray(dummy, ref multiplots[I, J].dataGroups);
                }
            }




        }



    }
}
