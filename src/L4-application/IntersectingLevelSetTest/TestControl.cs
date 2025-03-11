using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IntersectingLevelSetTest {
    internal class TestControl : AppControl {

        public TestControl() { }

        public TestControl(Func<double, double, double, double> levelSet0, Func<double, double, double, double> levelSet1) {
            LevelSet0 = levelSet0;
            LevelSet1 = levelSet1;
        }

        public Func<double, double, double, double> LevelSet0 { get; set; }

        public Func<double, double, double, double> LevelSet1 { get; set; }

        public Func<double, double, double, double, double> LevelSet3D_0 { get; set; }

        public Func<double, double, double, double, double> LevelSet3D_1 { get; set; }

        public int Resolution = 2;

        public double ErrorThreshold = 1e-6;

        public int Dimension = 2;

        public double[] ErrorList = new double[4];

    }

    //internal class TestControl3D : AppControl {

    //    public TestControl3D() { }

    //    public TestControl3D(Func<double, double, double, double, double> levelSet3D0, Func<double, double, double, double, double> levelSet3D1) {
    //        LevelSet3D0 = levelSet3D0;
    //        LevelSet3D1 = levelSet3D1;
    //    }

    //    public Func<double, double, double, double, double> LevelSet3D0 { get; set; }

    //    public Func<double, double, double, double, double> LevelSet3D1 { get; set; }

    //    public int Resolution = 2;

    //    public double ErrorThreshold = 1e-6;

    //}
}
