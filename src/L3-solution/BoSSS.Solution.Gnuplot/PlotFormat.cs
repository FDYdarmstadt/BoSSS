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

using BoSSS.Platform;
using Newtonsoft.Json;
using System;
using System.Linq;
using System.Runtime.Serialization;


namespace BoSSS.Solution.Gnuplot {

    /// <summary>
    /// General Style of a plot, i.e. line type, color, etc.
    /// </summary>
    [Serializable]
    [DataContract]
    public class PlotFormat {

        static bool ContainsRemove(ref string f, string token) {
            token = token.ToLowerInvariant();
            bool ret = f.Contains(token);

            if(ret)
                f = f.Replace(token, "");

            return ret;
        }

        /// <summary>
        /// Sets the plot from a Matlab-like string, e.g. '--kx' (dashed, black, x-marks).
        /// </summary>
        public void FromString(string FormatString) {
            //linestyle
            FormatString = FormatString.ToLowerInvariant();

            if (ContainsRemove(ref FormatString, "Black"))
                this.LineColor = LineColors.Black;
            if (ContainsRemove(ref FormatString, "Red"))
                this.LineColor = LineColors.Red;
            if (ContainsRemove(ref FormatString, "Green"))
                this.LineColor = LineColors.Green;
            if (ContainsRemove(ref FormatString, "Blue"))
                this.LineColor = LineColors.Blue;
            if (ContainsRemove(ref FormatString, "Yellow"))
                this.LineColor = LineColors.Yellow;
            if (ContainsRemove(ref FormatString, "Magenta"))
                this.LineColor = LineColors.Magenta;

            if (ContainsRemove(ref FormatString, "k"))
                this.LineColor = LineColors.Black;
            if (ContainsRemove(ref FormatString, "r"))
                this.LineColor = LineColors.Red;
            if (ContainsRemove(ref FormatString, "g"))
                this.LineColor = LineColors.Green;
            if (ContainsRemove(ref FormatString, "b"))
                this.LineColor = LineColors.Blue;
            if (ContainsRemove(ref FormatString, "y"))
                this.LineColor = LineColors.Yellow;
            if (ContainsRemove(ref FormatString, "m"))
                this.LineColor = LineColors.Magenta;



            if (ContainsRemove(ref FormatString, "--"))
                this.DashType = DashTypes.Dashed;
            if (ContainsRemove(ref FormatString, ":"))
                this.DashType = DashTypes.Dotted;
            if (ContainsRemove(ref FormatString, "-."))
                this.DashType = DashTypes.DotDashed;
            if (ContainsRemove(ref FormatString, "-"))
                this.DashType = DashTypes.Solid;

            if (ContainsRemove(ref FormatString, "+"))
                this.PointType = PointTypes.Plus;
            if (ContainsRemove(ref FormatString, "o"))
                this.PointType = PointTypes.OpenCircle;
            if (ContainsRemove(ref FormatString, "*"))
                this.PointType = PointTypes.Asterisk;
            if (ContainsRemove(ref FormatString, "x"))
                this.PointType = PointTypes.Times;
            if (ContainsRemove(ref FormatString, "s"))
                this.PointType = PointTypes.OpenBox;
            if (ContainsRemove(ref FormatString, "^"))
                this.PointType = PointTypes.OpenLowerTriangle;
            if (ContainsRemove(ref FormatString, "^"))
                this.PointType = PointTypes.LowerTriangle;

        }

               
        /// <summary>
        /// Constructor
        /// </summary>
        public PlotFormat(PlotFormat baseLineFormat = null, LineColors? lineColor = null, DashTypes? dashType = null, double? lineWidth = null, Styles? Style = null, PointTypes? pointType = null, double? pointSize = null) {
            if (baseLineFormat != null) {
                this.lineColor = lineColor ?? baseLineFormat.lineColor;
                this.dashType = dashType ?? baseLineFormat.dashType;
                this.lineWidth = lineWidth ?? baseLineFormat.lineWidth;
                this._Style = Style ?? baseLineFormat.Style;
                this.pointType = pointType ?? baseLineFormat.pointType;
                this.pointSize = pointSize ?? baseLineFormat.pointSize;
            } else {
                this.lineColor = lineColor;
                this.dashType = dashType;
                this.lineWidth = lineWidth;
                this._Style = Style;
                this.pointType = pointType;
                this.pointSize = pointSize;
            }
        }

        [DataMember]
        private LineColors? lineColor;

        /// <summary>
        /// Color
        /// </summary>
        [JsonIgnore]
        public LineColors LineColor {
            get {
                return lineColor ?? LineColors.Red;
            }
            set {
                lineColor = value;
            }
        }

        [DataMember]
        private DashTypes? dashType;

        /// <summary>
        /// In a plot which uses lines, the style of those.
        /// </summary>
        [JsonIgnore]
        public DashTypes DashType {
            get {
                return dashType ?? DashTypes.Solid;
            }
            set {
                dashType = value;
            }
        }

        [DataMember]
        private double? lineWidth;

        /// <summary>
        /// line width
        /// </summary>
        [JsonIgnore]
        public double LineWidth {
            get {
                return lineWidth ?? 1.0;
            }
            set {
                lineWidth = value;
            }
        }

        [DataMember]
        private Styles? _Style;

        /// <summary>
        /// Basic plot type, e.g. lines, points, bars, etc.
        /// </summary>
        [JsonIgnore]
        public Styles Style {
            get {
                return _Style ?? Styles.Lines;
            }
            set {
                _Style = value;
            }
        }

        [DataMember]
        private PointTypes? pointType;

        /// <summary>
        /// Point type, only effective for <see cref="Styles.Points"/> and <see cref="Styles.LinesPoints"/>.
        /// </summary>
        [JsonIgnore]
        public PointTypes PointType {
            get {
                return pointType ?? PointTypes.OpenCircle; // dass man zumindest irgendwas sieht, wenn man auf 'LinesPoints' stellt
            }
            set {
                pointType = value;
            }
        }

        [DataMember]
        private double? pointSize;

        /// <summary>
        /// Size, may depend on output driver (screen or LaTeX),
        /// only effective for <see cref="Styles.Points"/> and <see cref="Styles.LinesPoints"/>.
        /// </summary>
        [JsonIgnore]
        public double PointSize {
            get {
                return pointSize ?? 1.0;
            }
            set {
                pointSize = value;
            }
        }

        /// <summary>
        /// Cloning and setting <see cref="LineColor"/>.
        /// </summary>
        public PlotFormat WithLineColor(LineColors color) {
            return new PlotFormat(baseLineFormat: this, lineColor: color);
        }

        /// <summary>
        /// Cloning and setting <see cref="LineColor"/>.
        /// </summary>
        public PlotFormat WithLineColor(int i) {
            int numberOfLineColors = Enum<LineColors>.GetValues().Count();
            return new PlotFormat(
                baseLineFormat: this,
                lineColor: (LineColors)((i + 1) % numberOfLineColors));
        }

        /// <summary>
        /// Cloning and setting <see cref="LineWidth"/>.
        /// </summary>
        public PlotFormat WithLineWidth(double lineWidth) {
            return new PlotFormat(baseLineFormat: this, lineWidth: LineWidth);
        }

        /// <summary>
        /// Cloning and setting <see cref="DashType"/>.
        /// </summary>
        public PlotFormat WithDashType(DashTypes dashtype) {
            return new PlotFormat(baseLineFormat: this, dashType: dashtype);
        }

        /// <summary>
        /// Cloning and setting <see cref="DashType"/>.
        /// </summary>
        public PlotFormat WithDashType(int i) {
            int numberOfDashTypes = Enum<DashTypes>.GetValues().Count();
            return new PlotFormat(
                baseLineFormat: this,
                dashType: (DashTypes)((i + 1) % numberOfDashTypes));
        }

        /// <summary>
        /// Cloning and setting <see cref="PointSize"/>.
        /// </summary>
        public PlotFormat WithPointSize(double pointSize) {
            return new PlotFormat(baseLineFormat: this, pointSize: pointSize);
        }

        /// <summary>
        /// Cloning and setting <see cref="pointType"/>.
        /// </summary>
        public PlotFormat WithPointType(PointTypes pointType) {
            return new PlotFormat(baseLineFormat: this, pointType: pointType);
        }

        /// <summary>
        /// Cloning and setting <see cref="pointType"/>.
        /// </summary>
        public PlotFormat WithPointType(int i) {
            int numberOfPointTypes = Enum<PointTypes>.GetValues().Count();
            return new PlotFormat(
                baseLineFormat: this,
                pointType: (PointTypes)(i % numberOfPointTypes));
        }

        /// <summary>
        /// Cloning and setting <see cref="Style"/>.
        /// </summary>
        public PlotFormat WithStyle(Styles style) {
            return new PlotFormat(baseLineFormat: this, Style: style);
        }
    }
}
