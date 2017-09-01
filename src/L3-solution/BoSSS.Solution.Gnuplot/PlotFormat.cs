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
using System.Linq;

namespace BoSSS.Solution.Gnuplot {

    public class PlotFormat {
               
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

        private readonly LineColors? lineColor;

        public LineColors LineColor {
            get {
                return lineColor ?? LineColors.Red;
            }
        }

        private readonly DashTypes? dashType;

        public DashTypes DashType {
            get {
                return dashType ?? DashTypes.Solid;
            }
        }

        private readonly double? lineWidth;

        public double LineWidth {
            get {
                return lineWidth ?? 1.0;
            }
        }

        private readonly Styles? _Style;

        public Styles Style {
            get {
                return _Style ?? Styles.Lines;
            }
        }

        private readonly PointTypes? pointType;

        public PointTypes PointType {
            get {
                return pointType ?? PointTypes.OpenCircle; // dass man zumindest irgendwas sieht, wenn man auf 'LinesPoints' stellt
            }
        }

        private readonly double? pointSize;

        public double PointSize {
            get {
                return pointSize ?? 1.0;
            }
        }

        public PlotFormat WithLineColor(LineColors color) {
            return new PlotFormat(baseLineFormat: this, lineColor: color);
        }

        public PlotFormat WithLineColor(int i) {
            int numberOfLineColors = Enum<LineColors>.GetValues().Count();
            return new PlotFormat(
                baseLineFormat: this,
                lineColor: (LineColors)((i + 1) % numberOfLineColors));
        }

        public PlotFormat WithLineWidth(double lineWidth) {
            return new PlotFormat(baseLineFormat: this, lineWidth: LineWidth);
        }

        public PlotFormat WithDashType(DashTypes dashtype) {
            return new PlotFormat(baseLineFormat: this, dashType: dashtype);
        }

        public PlotFormat WithDashType(int i) {
            int numberOfDashTypes = Enum<DashTypes>.GetValues().Count();
            return new PlotFormat(
                baseLineFormat: this,
                dashType: (DashTypes)((i + 1) % numberOfDashTypes));
        }

        public PlotFormat WithPointSize(double pointSize) {
            return new PlotFormat(baseLineFormat: this, pointSize: pointSize);
        }

        public PlotFormat WithPointType(PointTypes pointType) {
            return new PlotFormat(baseLineFormat: this, pointType: pointType);
        }

        public PlotFormat WithPointType(int i) {
            int numberOfPointTypes = Enum<PointTypes>.GetValues().Count();
            return new PlotFormat(
                baseLineFormat: this,
                pointType: (PointTypes)(i % numberOfPointTypes));
        }

        public PlotFormat WithStyle(Styles style) {
            return new PlotFormat(baseLineFormat: this, Style: style);
        }
    }
}
