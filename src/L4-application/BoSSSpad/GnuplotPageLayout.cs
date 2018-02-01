using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Additional layout options for gnuplot output which are not put into the <see cref="Plot2Ddata"/> container.
    /// </summary>
    public class GnuplotPageLayout {

        /// <summary>
        /// In Cairolatex mode, the horizontal size in centimeters.
        /// </summary>
        public double Cairolatex_xSizeCm = 14.0;

        /// <summary>
        /// In Cairolatex mode, the vertical size in centimeters.
        /// </summary>
        public double Cairolatex_ySizeCm = 10.5;




    }
}
